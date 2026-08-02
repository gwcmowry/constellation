use std::fs::File;
use std::io::{self, BufRead, BufReader};
use std::path::Path;

#[derive(Clone, Debug, Eq, PartialEq)]
pub struct FastqRecord {
    pub id: String,
    pub seq: Vec<u8>,
    pub qual: Vec<u8>,
}

pub fn read_fastq(path: impl AsRef<Path>) -> io::Result<Vec<FastqRecord>> {
    let file = File::open(path)?;
    let mut lines = BufReader::new(file).lines();
    let mut records = Vec::new();

    loop {
        let Some(id_line) = lines.next() else {
            break;
        };
        let id_line = id_line?;
        if id_line.trim().is_empty() {
            continue;
        }
        let seq = require_line(&mut lines, "sequence")?;
        let plus = require_line(&mut lines, "separator")?;
        let qual = require_line(&mut lines, "quality")?;

        if !id_line.starts_with('@') {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                format!("FASTQ record id does not start with @: {id_line}"),
            ));
        }
        if !plus.starts_with('+') {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                format!(
                    "FASTQ separator does not start with + for {}",
                    &id_line[1..]
                ),
            ));
        }
        if seq.len() != qual.len() {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                format!("sequence/quality length mismatch for {}", &id_line[1..]),
            ));
        }

        records.push(FastqRecord {
            id: id_line[1..].trim().to_string(),
            seq: seq.into_bytes(),
            qual: qual.into_bytes(),
        });
    }

    Ok(records)
}

fn require_line(
    lines: &mut impl Iterator<Item = io::Result<String>>,
    kind: &str,
) -> io::Result<String> {
    match lines.next() {
        Some(line) => line,
        None => Err(io::Error::new(
            io::ErrorKind::UnexpectedEof,
            format!("missing FASTQ {kind} line"),
        )),
    }
}

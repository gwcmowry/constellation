use std::fs::File;
use std::io::{self, BufRead, BufReader};
use std::path::Path;

#[derive(Clone, Debug, Eq, PartialEq)]
pub struct FastaRecord {
    pub header: String,
    pub seq: Vec<u8>,
}

pub fn read_fasta(path: impl AsRef<Path>) -> io::Result<Vec<FastaRecord>> {
    let file = File::open(path)?;
    let reader = BufReader::new(file);
    let mut records = Vec::new();
    let mut header: Option<String> = None;
    let mut seq = Vec::new();

    for line in reader.lines() {
        let line = line?;
        if let Some(rest) = line.strip_prefix('>') {
            if let Some(previous) = header.replace(rest.trim().to_string()) {
                records.push(FastaRecord {
                    header: previous,
                    seq: std::mem::take(&mut seq),
                });
            }
        } else {
            seq.extend(line.trim().as_bytes());
        }
    }

    if let Some(header) = header {
        records.push(FastaRecord { header, seq });
    }

    Ok(records)
}

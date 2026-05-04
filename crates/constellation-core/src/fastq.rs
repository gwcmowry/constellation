use flate2::read::MultiGzDecoder;
use flate2::write::GzEncoder;
use flate2::Compression;
use std::fs;
use std::io::{Read, Write};
use std::path::Path;
use thiserror::Error;

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct FastqRecord {
    pub id: String,
    pub seq: Vec<u8>,
    pub qual: Vec<u8>,
}

#[derive(Debug, Error)]
pub enum FastqError {
    #[error("FASTQ I/O error: {0}")]
    Io(#[from] std::io::Error),
    #[error("malformed FASTQ near record {record}: {message}")]
    Malformed { record: usize, message: String },
}

pub fn read_fastq(path: impl AsRef<Path>) -> Result<Vec<FastqRecord>, FastqError> {
    let path = path.as_ref();
    let path_text = path.display().to_string();
    let text = if path_text.ends_with(".gz") {
        let file = fs::File::open(path)?;
        let mut decoder = MultiGzDecoder::new(file);
        let mut text = String::new();
        decoder.read_to_string(&mut text)?;
        text
    } else {
        fs::read_to_string(path)?
    };
    let mut lines = text.lines();
    let mut out = Vec::new();
    let mut record = 0;

    while let Some(id_line) = lines.next() {
        record += 1;
        let seq = lines.next().ok_or_else(|| FastqError::Malformed {
            record,
            message: "missing sequence line".to_owned(),
        })?;
        let plus = lines.next().ok_or_else(|| FastqError::Malformed {
            record,
            message: "missing plus line".to_owned(),
        })?;
        let qual = lines.next().ok_or_else(|| FastqError::Malformed {
            record,
            message: "missing quality line".to_owned(),
        })?;

        if !id_line.starts_with('@') {
            return Err(FastqError::Malformed {
                record,
                message: "identifier line must start with @".to_owned(),
            });
        }
        if !plus.starts_with('+') {
            return Err(FastqError::Malformed {
                record,
                message: "separator line must start with +".to_owned(),
            });
        }
        if seq.len() != qual.len() {
            return Err(FastqError::Malformed {
                record,
                message: format!(
                    "sequence length {} != quality length {}",
                    seq.len(),
                    qual.len()
                ),
            });
        }

        out.push(FastqRecord {
            id: id_line[1..]
                .split_whitespace()
                .next()
                .unwrap_or("")
                .to_owned(),
            seq: seq.as_bytes().to_ascii_uppercase(),
            qual: qual.as_bytes().to_vec(),
        });
    }

    Ok(out)
}

pub fn write_fastq(path: impl AsRef<Path>, records: &[FastqRecord]) -> Result<(), std::io::Error> {
    let mut out = String::new();
    for record in records {
        out.push('@');
        out.push_str(&record.id);
        out.push('\n');
        out.push_str(&String::from_utf8_lossy(&record.seq));
        out.push_str("\n+\n");
        out.push_str(&String::from_utf8_lossy(&record.qual));
        out.push('\n');
    }
    let path = path.as_ref();
    if path.display().to_string().ends_with(".gz") {
        let file = fs::File::create(path)?;
        let mut encoder = GzEncoder::new(file, Compression::default());
        encoder.write_all(out.as_bytes())?;
        encoder.finish()?;
        Ok(())
    } else {
        fs::write(path, out)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn gzip_roundtrip() {
        let tmp = tempfile::Builder::new()
            .suffix(".fastq.gz")
            .tempfile()
            .unwrap();
        let records = vec![FastqRecord {
            id: "read0".to_owned(),
            seq: b"ACGT".to_vec(),
            qual: b"IIII".to_vec(),
        }];
        write_fastq(tmp.path(), &records).unwrap();
        assert_eq!(read_fastq(tmp.path()).unwrap(), records);
    }
}

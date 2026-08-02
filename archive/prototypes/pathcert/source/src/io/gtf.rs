use std::collections::HashMap;
use std::fs::File;
use std::io::{self, BufRead, BufReader};
use std::path::Path;

#[derive(Clone, Debug, Eq, PartialEq)]
pub struct ExonRecord {
    pub seqname: String,
    pub start: usize,
    pub end: usize,
    pub strand: char,
    pub gene_id: String,
    pub transcript_id: String,
}

pub fn read_gtf_exons(path: impl AsRef<Path>) -> io::Result<Vec<ExonRecord>> {
    let file = File::open(path)?;
    let reader = BufReader::new(file);
    let mut exons = Vec::new();

    for line in reader.lines() {
        let line = line?;
        let trimmed = line.trim();
        if trimmed.is_empty() || trimmed.starts_with('#') {
            continue;
        }
        let fields: Vec<&str> = trimmed.split('\t').collect();
        if fields.len() < 9 || fields[2] != "exon" {
            continue;
        }
        let attrs = parse_attributes(fields[8]);
        let Some(gene_id) = attrs.get("gene_id").cloned() else {
            continue;
        };
        let Some(transcript_id) = attrs.get("transcript_id").cloned() else {
            continue;
        };
        let start = fields[3].parse::<usize>().map_err(invalid_data)?;
        let end = fields[4].parse::<usize>().map_err(invalid_data)?;
        let strand = fields[6].chars().next().unwrap_or('+');
        exons.push(ExonRecord {
            seqname: fields[0].to_string(),
            start,
            end,
            strand,
            gene_id,
            transcript_id,
        });
    }

    Ok(exons)
}

pub fn parse_attributes(attributes: &str) -> HashMap<String, String> {
    let mut parsed = HashMap::new();
    for item in attributes.split(';') {
        let item = item.trim();
        if item.is_empty() {
            continue;
        }
        let mut pieces = item.splitn(2, char::is_whitespace);
        let Some(key) = pieces.next() else {
            continue;
        };
        let value = pieces
            .next()
            .unwrap_or("")
            .trim()
            .trim_matches('"')
            .to_string();
        if !key.is_empty() && !value.is_empty() {
            parsed.insert(key.to_string(), value);
        }
    }
    parsed
}

fn invalid_data(error: impl ToString) -> io::Error {
    io::Error::new(io::ErrorKind::InvalidData, error.to_string())
}

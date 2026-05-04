use rustc_hash::FxHashMap;
use std::fs;
use std::path::Path;
use thiserror::Error;

#[derive(Debug, Error)]
pub enum GtfError {
    #[error("GTF I/O error: {0}")]
    Io(#[from] std::io::Error),
}

pub fn transcript_gene_map(path: impl AsRef<Path>) -> Result<FxHashMap<String, String>, GtfError> {
    let text = fs::read_to_string(path)?;
    let mut out = FxHashMap::default();
    for line in text.lines() {
        if line.trim().is_empty() || line.starts_with('#') {
            continue;
        }
        let fields: Vec<_> = line.split('\t').collect();
        if fields.len() < 9 {
            continue;
        }
        let attrs = parse_attributes(fields[8]);
        let Some(transcript_id) = attrs.get("transcript_id") else {
            continue;
        };
        let Some(gene_id) = attrs.get("gene_id") else {
            continue;
        };
        out.entry(transcript_id.clone())
            .or_insert_with(|| gene_id.clone());
    }
    Ok(out)
}

pub fn parse_attributes(attributes: &str) -> FxHashMap<String, String> {
    let mut out = FxHashMap::default();
    for part in attributes.split(';') {
        let part = part.trim();
        if part.is_empty() {
            continue;
        }
        let mut split = part.splitn(2, char::is_whitespace);
        let Some(key) = split.next() else {
            continue;
        };
        let Some(value) = split.next() else {
            continue;
        };
        out.insert(key.to_owned(), value.trim().trim_matches('"').to_owned());
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn parses_gtf_attributes() {
        let attrs = parse_attributes(r#"gene_id "GENE"; transcript_id "TX1";"#);
        assert_eq!(attrs.get("gene_id").unwrap(), "GENE");
        assert_eq!(attrs.get("transcript_id").unwrap(), "TX1");
    }
}

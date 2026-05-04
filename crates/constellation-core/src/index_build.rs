use crate::dna::{encode_acgt, iter_kmers_2bit};
use crate::index::{KmerEntry, Posting, TranscriptIndex, TranscriptMeta};
use crate::{GeneId, TranscriptId};
use rayon::prelude::*;
use rustc_hash::FxHashMap;
use std::collections::BTreeMap;
use std::path::Path;
use thiserror::Error;

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct FastaRecord {
    pub header: String,
    pub name: String,
    pub gene: String,
    pub seq: Vec<u8>,
}

#[derive(Debug, Error)]
pub enum BuildIndexError {
    #[error("FASTA I/O error: {0}")]
    Io(#[from] std::io::Error),
    #[error("malformed FASTA: {0}")]
    Malformed(String),
    #[error("k must be between 1 and 32, got {0}")]
    InvalidK(u8),
}

pub fn read_fasta(path: impl AsRef<Path>) -> Result<Vec<FastaRecord>, BuildIndexError> {
    read_fasta_needletail(path)
}

fn read_fasta_needletail(path: impl AsRef<Path>) -> Result<Vec<FastaRecord>, BuildIndexError> {
    let mut reader = needletail::parse_fastx_file(path.as_ref())
        .map_err(|err| BuildIndexError::Malformed(err.to_string()))?;
    let mut records = Vec::new();
    while let Some(record) = reader.next() {
        let record = record.map_err(|err| BuildIndexError::Malformed(err.to_string()))?;
        let header = String::from_utf8_lossy(record.id()).into_owned();
        let seq = record
            .seq()
            .iter()
            .map(|base| base.to_ascii_uppercase())
            .collect();
        records.push(make_record(header, seq));
    }

    if records.is_empty() {
        return Err(BuildIndexError::Malformed("no records found".to_owned()));
    }
    Ok(records)
}

pub fn build_transcript_index(
    fasta: impl AsRef<Path>,
    k: u8,
    max_kmer_frequency: u32,
) -> Result<TranscriptIndex, BuildIndexError> {
    build_transcript_index_with_gene_map(fasta, k, max_kmer_frequency, None)
}

pub fn build_transcript_index_with_gene_map(
    fasta: impl AsRef<Path>,
    k: u8,
    max_kmer_frequency: u32,
    transcript_gene_map: Option<&FxHashMap<String, String>>,
) -> Result<TranscriptIndex, BuildIndexError> {
    if !(1..=32).contains(&k) {
        return Err(BuildIndexError::InvalidK(k));
    }

    let records = read_fasta(fasta)?;
    let mut gene_ids: FxHashMap<String, GeneId> = FxHashMap::default();
    let mut genes = Vec::new();
    let mut transcripts = Vec::new();
    let mut transcript_sequences = Vec::new();
    let mut postings_by_kmer: BTreeMap<u64, Vec<Posting>> = BTreeMap::new();

    for (idx, record) in records.iter().enumerate() {
        let gene_name = transcript_gene_map
            .and_then(|map| map.get(&record.name))
            .cloned()
            .unwrap_or_else(|| record.gene.clone());
        let gene_id = match gene_ids.get(&gene_name) {
            Some(&id) => id,
            None => {
                let id = genes.len() as GeneId;
                gene_ids.insert(gene_name.clone(), id);
                genes.push(gene_name);
                id
            }
        };
        let transcript_id = idx as TranscriptId;
        transcripts.push(TranscriptMeta {
            transcript_id,
            gene_id,
            name: record.name.clone(),
            len: record.seq.len() as u32,
        });
        transcript_sequences.push(String::from_utf8_lossy(&record.seq).into_owned());

        let encoded = encode_acgt(&record.seq);
        for kmer in iter_kmers_2bit(&encoded, k) {
            postings_by_kmer
                .entry(kmer.code)
                .or_default()
                .push(Posting {
                    transcript_id,
                    gene_id,
                    pos: kmer.pos,
                    strand: 0,
                });
        }
    }

    let mut kmers = Vec::with_capacity(postings_by_kmer.len());
    let mut postings = Vec::new();
    for (code, mut kmer_postings) in postings_by_kmer {
        kmer_postings
            .par_sort_by_key(|posting| (posting.transcript_id, posting.pos, posting.strand));
        let start = postings.len() as u64;
        let len = kmer_postings.len() as u32;
        postings.extend(kmer_postings);
        kmers.push(KmerEntry {
            kmer_code: code,
            postings_start: start,
            postings_len: len,
            freq_class: if len > max_kmer_frequency { 1 } else { 0 },
        });
    }

    Ok(TranscriptIndex {
        format_version: 1,
        k,
        max_kmer_frequency,
        genes,
        transcripts,
        transcript_sequences,
        kmers,
        postings,
    })
}

fn make_record(header: String, seq: Vec<u8>) -> FastaRecord {
    let name = header
        .split_whitespace()
        .next()
        .unwrap_or(&header)
        .split('|')
        .next()
        .unwrap_or(&header)
        .to_owned();
    let gene = header
        .split(['|', ' ', '\t'])
        .find_map(|part| {
            part.strip_prefix("gene=")
                .or_else(|| part.strip_prefix("gene:"))
        })
        .map(|gene| gene.to_owned())
        .unwrap_or_else(|| name.clone());
    FastaRecord {
        header,
        name,
        gene,
        seq,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::dna::kmer_code_from_ascii;
    use std::io::Write;

    #[test]
    fn tiny_index_contains_expected_kmers() {
        let mut fasta = tempfile::NamedTempFile::new().unwrap();
        writeln!(
            fasta,
            ">tx1|gene=GENE_A\nACGTACGT\n>tx2|gene=GENE_B\nTTTACG"
        )
        .unwrap();
        let index = build_transcript_index(fasta.path(), 3, 256).unwrap();
        let acg = kmer_code_from_ascii(b"ACG").unwrap();
        let postings = index.lookup(acg);
        assert_eq!(index.genes, vec!["GENE_A", "GENE_B"]);
        assert_eq!(postings.len(), 3);
    }

    #[test]
    fn gtf_gene_map_overrides_fasta_gene() {
        let mut fasta = tempfile::NamedTempFile::new().unwrap();
        writeln!(fasta, ">tx1\nACGTACGT").unwrap();
        let mut map = FxHashMap::default();
        map.insert("tx1".to_owned(), "GENE_FROM_GTF".to_owned());
        let index = build_transcript_index_with_gene_map(fasta.path(), 3, 256, Some(&map)).unwrap();
        assert_eq!(index.genes, vec!["GENE_FROM_GTF"]);
    }

    #[test]
    fn parses_ensembl_cdna_gene_colon_header() {
        let record = make_record(
            "ENST1 cdna chromosome:GRCh38:1:1:10:1 gene:ENSG1 gene_symbol:ABC".to_owned(),
            b"ACGT".to_vec(),
        );
        assert_eq!(record.name, "ENST1");
        assert_eq!(record.gene, "ENSG1");
    }
}

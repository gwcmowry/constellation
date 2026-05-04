use crate::index::{Posting, TranscriptIndex};

pub fn lookup_kmer(index: &TranscriptIndex, kmer_code: u64) -> &[Posting] {
    index.lookup(kmer_code)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::dna::kmer_code_from_ascii;
    use crate::index::{KmerEntry, TranscriptIndex, TranscriptMeta};

    #[test]
    fn known_kmer_returns_expected_postings() {
        let code = kmer_code_from_ascii(b"ACG").unwrap();
        let index = TranscriptIndex {
            format_version: 1,
            k: 3,
            max_kmer_frequency: 256,
            genes: vec!["GENE".to_owned()],
            transcripts: vec![TranscriptMeta {
                transcript_id: 0,
                gene_id: 0,
                name: "tx".to_owned(),
                len: 3,
            }],
            transcript_sequences: vec!["ACG".to_owned()],
            kmers: vec![KmerEntry {
                kmer_code: code,
                postings_start: 0,
                postings_len: 1,
                freq_class: 0,
                transcript_df: 1,
                gene_df: 1,
            }],
            postings: vec![crate::index::Posting {
                transcript_id: 0,
                gene_id: 0,
                pos: 0,
                strand: 0,
            }],
        };
        assert_eq!(lookup_kmer(&index, code).len(), 1);
    }
}

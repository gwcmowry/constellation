use crate::candidate::CandidateLocusBucket;
use crate::index::IndexAccess;
use crate::{GeneId, ReadId, TranscriptId};

#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord)]
pub struct ScoredCandidate {
    pub read_id: ReadId,
    pub transcript_id: TranscriptId,
    pub gene_id: GeneId,
    pub pos: u32,
    pub mismatches: u32,
    pub score: u16,
    pub flags: u16,
}

pub trait CandidateScorer: Send + Sync {
    fn score_bucket(
        &self,
        reads: &[(ReadId, Vec<u8>)],
        index: &dyn IndexAccess,
        bucket: &CandidateLocusBucket,
        out: &mut Vec<ScoredCandidate>,
    );
}

pub fn read_seq_by_id(reads: &[(ReadId, Vec<u8>)], read_id: ReadId) -> Option<&[u8]> {
    let idx = usize::try_from(read_id).ok()?;
    let (stored_id, seq) = reads.get(idx)?;
    if *stored_id == read_id {
        Some(seq.as_slice())
    } else {
        None
    }
}

#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub struct ScoreBlock {
    pub read_ids: Vec<ReadId>,
    pub transcript_ids: Vec<TranscriptId>,
    pub gene_ids: Vec<GeneId>,
    pub pos: Vec<u32>,
    pub len: Vec<u16>,
    pub read_offsets: Vec<u32>,
    pub ref_offsets: Vec<u32>,
    pub read_bases: Vec<u8>,
    pub ref_bases: Vec<u8>,
}

impl ScoreBlock {
    pub fn from_bucket(
        reads: &[(ReadId, Vec<u8>)],
        index: &dyn IndexAccess,
        bucket: &CandidateLocusBucket,
    ) -> Self {
        let mut block = Self::default();
        for hit in &bucket.hits {
            let Some(read_seq) = read_seq_by_id(reads, hit.read_id) else {
                continue;
            };
            let Some(transcript_seq) = index.transcript_seq(hit.transcript_id) else {
                continue;
            };
            let start = hit.pos as usize;
            let end = start + read_seq.len();
            if end > transcript_seq.len() {
                continue;
            }

            block.read_ids.push(hit.read_id);
            block.transcript_ids.push(hit.transcript_id);
            block.gene_ids.push(hit.gene_id);
            block.pos.push(hit.pos);
            block.len.push(read_seq.len().min(u16::MAX as usize) as u16);
            block.read_offsets.push(block.read_bases.len() as u32);
            block.ref_offsets.push(block.ref_bases.len() as u32);
            if hit.strand == 0 {
                block.read_bases.extend(read_seq);
            } else {
                block
                    .read_bases
                    .extend(read_seq.iter().rev().map(|&base| complement_ascii(base)));
            }
            block.ref_bases.extend(&transcript_seq[start..end]);
        }
        block
    }

    pub fn is_empty(&self) -> bool {
        self.read_ids.is_empty()
    }

    pub fn pairs(&self) -> impl Iterator<Item = (&[u8], &[u8])> {
        (0..self.read_ids.len()).map(|idx| {
            let len = self.len[idx] as usize;
            let read_start = self.read_offsets[idx] as usize;
            let ref_start = self.ref_offsets[idx] as usize;
            (
                &self.read_bases[read_start..read_start + len],
                &self.ref_bases[ref_start..ref_start + len],
            )
        })
    }
}

pub fn hamming_ascii(a: &[u8], b: &[u8]) -> u32 {
    let shared = a.len().min(b.len());
    let mut mismatches = a.len().abs_diff(b.len()) as u32;
    for idx in 0..shared {
        if !a[idx].eq_ignore_ascii_case(&b[idx]) {
            mismatches += 1;
        }
    }
    mismatches
}

pub fn hamming_revcomp_ascii(read: &[u8], reference: &[u8]) -> u32 {
    let shared = read.len().min(reference.len());
    let mut mismatches = read.len().abs_diff(reference.len()) as u32;
    for idx in 0..shared {
        let read_base = complement_ascii(read[read.len() - 1 - idx]);
        if !read_base.eq_ignore_ascii_case(&reference[idx]) {
            mismatches += 1;
        }
    }
    mismatches
}

fn complement_ascii(base: u8) -> u8 {
    match base.to_ascii_uppercase() {
        b'A' => b'T',
        b'C' => b'G',
        b'G' => b'C',
        b'T' | b'U' => b'A',
        _ => b'N',
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::candidate::{CandidateBucketKey, CandidateHit};
    use crate::index::{TranscriptIndex, TranscriptMeta};

    #[test]
    fn score_block_materializes_parallel_arrays() {
        let index = TranscriptIndex {
            format_version: 1,
            k: 3,
            max_kmer_frequency: 256,
            genes: vec!["GENE".to_owned()],
            transcripts: vec![TranscriptMeta {
                transcript_id: 0,
                gene_id: 0,
                name: "tx".to_owned(),
                len: 8,
            }],
            transcript_sequences: vec!["ACGTACGT".to_owned()],
            kmers: vec![],
            postings: vec![],
        };
        let bucket = CandidateLocusBucket {
            key: CandidateBucketKey {
                transcript_id: 0,
                pos_bin: 0,
                strand: 0,
            },
            hits: vec![CandidateHit {
                read_id: 0,
                transcript_id: 0,
                gene_id: 0,
                pos: 2,
                strand: 0,
                seed_count: 1,
                seed_score: 1,
            }],
        };
        let block = ScoreBlock::from_bucket(&[(0, b"GTAC".to_vec())], &index, &bucket);
        assert_eq!(block.read_ids, vec![0]);
        assert_eq!(block.pairs().next().unwrap(), (&b"GTAC"[..], &b"GTAC"[..]));
    }

    #[test]
    fn hamming_reverse_complement_match_zero() {
        assert_eq!(hamming_revcomp_ascii(b"AAGT", b"ACTT"), 0);
    }
}

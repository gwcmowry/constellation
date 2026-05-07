use crate::candidate::CandidateLocusBucket;
use crate::index::IndexAccess;
use crate::{GeneId, ReadId, TranscriptId};

#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord)]
pub struct ScoredCandidate {
    pub read_id: ReadId,
    pub transcript_id: TranscriptId,
    pub gene_id: GeneId,
    pub pos: u32,
    pub strand: u8,
    pub mismatches: u32,
    pub score: u16,
    pub flags: u16,
}

pub const SCORE_FLAG_FULL_LENGTH: u16 = 1 << 0;
pub const SCORE_FLAG_TRIMMED_POLY_A: u16 = 1 << 1;
pub const SCORE_FLAG_TRIMMED_POLY_T: u16 = 1 << 2;
pub const SCORE_FLAG_TRIMMED_LOW_QUALITY: u16 = 1 << 3;
pub const SCORE_FLAG_RIGHT_SOFTCLIP: u16 = 1 << 4;
pub const SCORE_FLAG_LEFT_SOFTCLIP: u16 = 1 << 5;
pub const SCORE_FLAG_TRIMMED_TSO: u16 = 1 << 6;
pub const SCORE_FLAG_ANTISENSE: u16 = 1 << 7;
pub const SCORE_FLAG_TARGET_EXON: u16 = 1 << 8;
pub const SCORE_FLAG_TARGET_INTRON: u16 = 1 << 9;
pub const SCORE_FLAG_TARGET_GENE_BODY: u16 = 1 << 10;
pub const SCORE_FLAG_RESCUED: u16 = 1 << 11;

pub const TENX_3P_TSO: &[u8] = b"AAGCAGTGGTATCAACGCAGAGTACATGGG";

pub fn tso_prefix_trim_len(read_seq: &[u8], min_match_len: u8, max_mismatches: u8) -> usize {
    if read_seq.is_empty() {
        return 0;
    }
    let max_len = read_seq.len().min(TENX_3P_TSO.len());
    let min_len = usize::from(min_match_len).min(max_len);
    if min_len == 0 {
        return 0;
    }
    let mut best = 0;
    let mut mismatches = 0_u8;
    for len in 1..=max_len {
        let idx = len - 1;
        if !ascii_base_eq(read_seq[idx], TENX_3P_TSO[idx]) {
            mismatches = mismatches.saturating_add(1);
        }
        if len >= min_len && mismatches <= max_mismatches {
            best = len;
        }
    }
    best
}

#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub enum LibraryStrand {
    #[default]
    Unstranded,
    Forward,
    Reverse,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct ScoreConfig {
    pub max_mismatches: u32,
    pub max_right_softclip: u16,
    pub max_left_softclip: u16,
    pub trim_poly_a: bool,
    pub trim_poly_t: bool,
    pub trim_low_quality_tail: bool,
    pub trim_tso: bool,
    pub max_tso_mismatches: u8,
    pub min_tso_match_len: u8,
    pub library_strand: LibraryStrand,
    pub min_scored_len: u16,
    pub min_tail_phred: u8,
}

impl Default for ScoreConfig {
    fn default() -> Self {
        Self {
            max_mismatches: 6,
            max_right_softclip: 0,
            max_left_softclip: 0,
            trim_poly_a: false,
            trim_poly_t: false,
            trim_low_quality_tail: false,
            trim_tso: false,
            max_tso_mismatches: 3,
            min_tso_match_len: 10,
            library_strand: LibraryStrand::Unstranded,
            min_scored_len: 35,
            min_tail_phred: 10,
        }
    }
}

impl ScoreConfig {
    pub fn uses_default_window(self) -> bool {
        self.max_right_softclip == 0
            && self.max_left_softclip == 0
            && !self.trim_poly_a
            && !self.trim_poly_t
            && !self.trim_low_quality_tail
            && !self.trim_tso
    }
}

#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub struct ScoreFailureStats {
    pub candidates_seen: u64,
    pub candidates_out_of_bounds: u64,
    pub candidates_failed_mismatch: u64,
    pub candidates_failed_min_scored_len: u64,
    pub candidates_failed_trimmed_out_of_bounds: u64,
    pub candidates_full_length_mismatch_only: u64,
    pub candidates_passed_full_length: u64,
    pub candidates_passed_trimmed_or_softclipped: u64,
}

pub trait CandidateScorer: Send + Sync {
    fn score_bucket(
        &self,
        reads: &[(ReadId, Vec<u8>)],
        quals: &[(ReadId, Vec<u8>)],
        index: &dyn IndexAccess,
        bucket: &CandidateLocusBucket,
        out: &mut Vec<ScoredCandidate>,
        read_score_stats: Option<&mut [ScoreFailureStats]>,
    ) -> ScoreFailureStats;
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

pub fn qual_by_read_id(quals: &[(ReadId, Vec<u8>)], read_id: ReadId) -> Option<&[u8]> {
    let idx = usize::try_from(read_id).ok()?;
    let (stored_id, qual) = quals.get(idx)?;
    if *stored_id == read_id {
        Some(qual.as_slice())
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
        if !ascii_base_eq(a[idx], b[idx]) {
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
        if !ascii_base_eq(read_base, reference[idx]) {
            mismatches += 1;
        }
    }
    mismatches
}

#[inline(always)]
fn ascii_base_eq(a: u8, b: u8) -> bool {
    a == b || (a | 0x20) == (b | 0x20)
}

#[inline(always)]
fn complement_ascii(base: u8) -> u8 {
    match base {
        b'A' | b'a' => b'T',
        b'C' | b'c' => b'G',
        b'G' | b'g' => b'C',
        b'T' | b't' | b'U' | b'u' => b'A',
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

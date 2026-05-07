use crate::candidate::CandidateLocusBucket;
use crate::index::{IndexAccess, TargetClass};
use crate::score::{
    hamming_ascii, hamming_revcomp_ascii, qual_by_read_id, read_seq_by_id, tso_prefix_trim_len,
    CandidateScorer, LibraryStrand, ScoreConfig, ScoreFailureStats, ScoredCandidate,
    SCORE_FLAG_ANTISENSE, SCORE_FLAG_FULL_LENGTH, SCORE_FLAG_LEFT_SOFTCLIP,
    SCORE_FLAG_RIGHT_SOFTCLIP, SCORE_FLAG_TARGET_EXON, SCORE_FLAG_TARGET_GENE_BODY,
    SCORE_FLAG_TARGET_INTRON, SCORE_FLAG_TRIMMED_LOW_QUALITY, SCORE_FLAG_TRIMMED_POLY_A,
    SCORE_FLAG_TRIMMED_POLY_T, SCORE_FLAG_TRIMMED_TSO,
};
use crate::ReadId;

#[derive(Debug, Clone, Copy)]
pub struct ScalarScorer {
    pub config: ScoreConfig,
    pub target_class_flags: &'static [u16],
}

impl Default for ScalarScorer {
    fn default() -> Self {
        Self {
            config: ScoreConfig::default(),
            target_class_flags: &[],
        }
    }
}

impl CandidateScorer for ScalarScorer {
    fn score_bucket(
        &self,
        reads: &[(ReadId, Vec<u8>)],
        quals: &[(ReadId, Vec<u8>)],
        index: &dyn IndexAccess,
        bucket: &CandidateLocusBucket,
        out: &mut Vec<ScoredCandidate>,
        mut read_score_stats: Option<&mut [ScoreFailureStats]>,
    ) -> ScoreFailureStats {
        let mut stats = ScoreFailureStats::default();
        for hit in &bucket.hits {
            let mut hit_stats = ScoreFailureStats {
                candidates_seen: 1,
                ..ScoreFailureStats::default()
            };
            let Some(read_seq) = read_seq_by_id(reads, hit.read_id) else {
                add_stats(&mut stats, hit_stats);
                add_read_stats(read_score_stats.as_deref_mut(), hit.read_id, hit_stats);
                continue;
            };
            let Some(transcript_seq) = index.transcript_seq(hit.transcript_id) else {
                add_stats(&mut stats, hit_stats);
                add_read_stats(read_score_stats.as_deref_mut(), hit.read_id, hit_stats);
                continue;
            };
            let qual = qual_by_read_id(quals, hit.read_id);
            let outcome = score_hit(self.config, hit, read_seq, qual, transcript_seq);
            let Some(mut candidate) = outcome.candidate else {
                match outcome.failure {
                    ScoreFailureReason::FullLengthOutOfBounds => {
                        hit_stats.candidates_out_of_bounds += 1;
                    }
                    ScoreFailureReason::FullLengthMismatchOnly => {
                        hit_stats.candidates_full_length_mismatch_only += 1;
                    }
                    ScoreFailureReason::MinScoredLen => {
                        hit_stats.candidates_failed_min_scored_len += 1;
                    }
                    ScoreFailureReason::TrimmedOutOfBounds => {
                        hit_stats.candidates_failed_trimmed_out_of_bounds += 1;
                    }
                    ScoreFailureReason::Mismatch => {
                        hit_stats.candidates_failed_mismatch += 1;
                    }
                    ScoreFailureReason::None => {}
                }
                add_stats(&mut stats, hit_stats);
                add_read_stats(read_score_stats.as_deref_mut(), hit.read_id, hit_stats);
                continue;
            };
            if is_antisense_candidate(self.config, hit) {
                candidate.flags |= SCORE_FLAG_ANTISENSE;
            }
            candidate.flags |= self.target_class_flag(index, hit.transcript_id);
            if candidate.flags & SCORE_FLAG_FULL_LENGTH != 0 {
                hit_stats.candidates_passed_full_length += 1;
            } else {
                hit_stats.candidates_passed_trimmed_or_softclipped += 1;
            }
            out.push(candidate);
            add_stats(&mut stats, hit_stats);
            add_read_stats(read_score_stats.as_deref_mut(), hit.read_id, hit_stats);
        }
        stats
    }
}

impl ScalarScorer {
    fn target_class_flag(
        &self,
        index: &dyn IndexAccess,
        transcript_id: crate::TranscriptId,
    ) -> u16 {
        self.target_class_flags
            .get(transcript_id as usize)
            .copied()
            .unwrap_or_else(|| target_class_flag(index.transcript_target_class(transcript_id)))
    }
}

pub fn target_class_flag(target_class: TargetClass) -> u16 {
    match target_class {
        TargetClass::ExonTranscript => SCORE_FLAG_TARGET_EXON,
        TargetClass::Intron => SCORE_FLAG_TARGET_INTRON,
        TargetClass::GeneBody => SCORE_FLAG_TARGET_GENE_BODY,
        TargetClass::Unknown => 0,
    }
}

fn is_antisense_candidate(config: ScoreConfig, hit: &crate::candidate::CandidateHit) -> bool {
    match config.library_strand {
        LibraryStrand::Unstranded => false,
        LibraryStrand::Forward => hit.strand == 1,
        LibraryStrand::Reverse => hit.strand == 0,
    }
}

fn add_stats(total: &mut ScoreFailureStats, next: ScoreFailureStats) {
    total.candidates_seen += next.candidates_seen;
    total.candidates_out_of_bounds += next.candidates_out_of_bounds;
    total.candidates_failed_mismatch += next.candidates_failed_mismatch;
    total.candidates_failed_min_scored_len += next.candidates_failed_min_scored_len;
    total.candidates_failed_trimmed_out_of_bounds += next.candidates_failed_trimmed_out_of_bounds;
    total.candidates_full_length_mismatch_only += next.candidates_full_length_mismatch_only;
    total.candidates_passed_full_length += next.candidates_passed_full_length;
    total.candidates_passed_trimmed_or_softclipped += next.candidates_passed_trimmed_or_softclipped;
}

fn add_read_stats(
    read_score_stats: Option<&mut [ScoreFailureStats]>,
    read_id: ReadId,
    next: ScoreFailureStats,
) {
    let Some(read_score_stats) = read_score_stats else {
        return;
    };
    let Some(slot) = read_score_stats.get_mut(read_id as usize) else {
        return;
    };
    add_stats(slot, next);
}

fn score_hit(
    config: ScoreConfig,
    hit: &crate::candidate::CandidateHit,
    read_seq: &[u8],
    qual: Option<&[u8]>,
    transcript_seq: &[u8],
) -> ScoreOutcome {
    let tso_trim = tso_trim_len(read_seq, config);
    let (read_window, qual_window) = trim_tso_prefix(read_seq, qual, tso_trim);
    let ref_offset = 0;
    let start = hit.pos as usize + ref_offset;
    let end = start + read_window.len();
    let mut full_length_failure = ScoreFailureReason::FullLengthOutOfBounds;
    if end <= transcript_seq.len() {
        let reference = &transcript_seq[start..end];
        let mismatches = if hit.strand == 0 {
            hamming_ascii(read_window, reference)
        } else {
            hamming_revcomp_ascii(read_window, reference)
        };
        if mismatches <= config.max_mismatches {
            return ScoreOutcome::candidate(ScoredCandidate {
                read_id: hit.read_id,
                transcript_id: hit.transcript_id,
                gene_id: hit.gene_id,
                pos: start as u32,
                strand: hit.strand,
                mismatches,
                score: read_window
                    .len()
                    .saturating_sub(mismatches as usize)
                    .saturating_sub(tso_trim)
                    .min(u16::MAX as usize) as u16,
                flags: if tso_trim > 0 {
                    SCORE_FLAG_TRIMMED_TSO
                } else {
                    SCORE_FLAG_FULL_LENGTH
                },
            });
        }
        full_length_failure = ScoreFailureReason::FullLengthMismatchOnly;
    }
    if config.uses_default_window() {
        return ScoreOutcome::failure(full_length_failure);
    }
    score_trimmed_or_softclipped(
        config,
        hit,
        read_window,
        qual_window,
        ref_offset,
        tso_trim,
        transcript_seq,
    )
}

fn score_trimmed_or_softclipped(
    config: ScoreConfig,
    hit: &crate::candidate::CandidateHit,
    read_seq: &[u8],
    qual: Option<&[u8]>,
    ref_offset: usize,
    already_trimmed: usize,
    transcript_seq: &[u8],
) -> ScoreOutcome {
    let (oriented_read, oriented_qual) = orient_read_and_qual(read_seq, qual, hit.strand);
    let (trimmed_end, trim_flags) = trimmed_end(&oriented_read, oriented_qual.as_deref(), config);
    let mut best: Option<ScoredCandidate> = None;
    let mut saw_len_ok = false;
    let mut saw_bounds_ok = false;
    let mut saw_mismatch_failure = false;
    let max_left = usize::from(config.max_left_softclip).min(trimmed_end);
    let max_right = usize::from(config.max_right_softclip).min(trimmed_end);
    for left_clip in 0..=max_left {
        for extra_right_clip in 0..=max_right {
            let Some(scored_end) = trimmed_end.checked_sub(extra_right_clip) else {
                continue;
            };
            if scored_end <= left_clip {
                continue;
            }
            let scored_len = scored_end - left_clip;
            if scored_len < usize::from(config.min_scored_len) {
                continue;
            }
            saw_len_ok = true;
            let ref_start = hit.pos as usize + ref_offset + left_clip;
            let ref_end = ref_start + scored_len;
            if ref_end > transcript_seq.len() {
                continue;
            }
            saw_bounds_ok = true;
            let read_window = &oriented_read[left_clip..scored_end];
            let reference = &transcript_seq[ref_start..ref_end];
            let mismatches = hamming_ascii(read_window, reference);
            if mismatches > config.max_mismatches {
                saw_mismatch_failure = true;
                continue;
            }
            let total_clipped = oriented_read
                .len()
                .saturating_sub(scored_len)
                .saturating_add(already_trimmed);
            let score = scored_len
                .saturating_sub(mismatches as usize)
                .saturating_sub(total_clipped)
                .min(u16::MAX as usize) as u16;
            let mut flags = trim_flags;
            if already_trimmed > 0 {
                flags |= SCORE_FLAG_TRIMMED_TSO;
            }
            if left_clip > 0 {
                flags |= SCORE_FLAG_LEFT_SOFTCLIP;
            }
            if extra_right_clip > 0 || trimmed_end < oriented_read.len() {
                flags |= SCORE_FLAG_RIGHT_SOFTCLIP;
            }
            let candidate = ScoredCandidate {
                read_id: hit.read_id,
                transcript_id: hit.transcript_id,
                gene_id: hit.gene_id,
                pos: hit.pos + ref_offset as u32 + left_clip as u32,
                strand: hit.strand,
                mismatches,
                score,
                flags,
            };
            if best
                .as_ref()
                .is_none_or(|current| candidate.score > current.score)
            {
                best = Some(candidate);
            }
        }
    }
    if let Some(candidate) = best {
        return ScoreOutcome::candidate(candidate);
    }
    if !saw_len_ok {
        ScoreOutcome::failure(ScoreFailureReason::MinScoredLen)
    } else if !saw_bounds_ok {
        ScoreOutcome::failure(ScoreFailureReason::TrimmedOutOfBounds)
    } else if saw_mismatch_failure {
        ScoreOutcome::failure(ScoreFailureReason::Mismatch)
    } else {
        ScoreOutcome::failure(ScoreFailureReason::Mismatch)
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum ScoreFailureReason {
    None,
    FullLengthOutOfBounds,
    FullLengthMismatchOnly,
    MinScoredLen,
    TrimmedOutOfBounds,
    Mismatch,
}

#[derive(Debug, Clone, PartialEq, Eq)]
struct ScoreOutcome {
    candidate: Option<ScoredCandidate>,
    failure: ScoreFailureReason,
}

impl ScoreOutcome {
    fn candidate(candidate: ScoredCandidate) -> Self {
        Self {
            candidate: Some(candidate),
            failure: ScoreFailureReason::None,
        }
    }

    fn failure(failure: ScoreFailureReason) -> Self {
        Self {
            candidate: None,
            failure,
        }
    }
}

fn orient_read_and_qual(
    read_seq: &[u8],
    qual: Option<&[u8]>,
    strand: u8,
) -> (Vec<u8>, Option<Vec<u8>>) {
    if strand == 0 {
        return (read_seq.to_vec(), qual.map(|q| q.to_vec()));
    }
    let read = read_seq
        .iter()
        .rev()
        .map(|&base| complement_ascii(base))
        .collect();
    let qual = qual.map(|q| q.iter().rev().copied().collect());
    (read, qual)
}

fn trim_tso_prefix<'a>(
    read_seq: &'a [u8],
    qual: Option<&'a [u8]>,
    trim_len: usize,
) -> (&'a [u8], Option<&'a [u8]>) {
    if trim_len == 0 {
        return (read_seq, qual);
    }
    (&read_seq[trim_len..], qual.and_then(|q| q.get(trim_len..)))
}

fn tso_trim_len(read_seq: &[u8], config: ScoreConfig) -> usize {
    if !config.trim_tso {
        return 0;
    }
    tso_prefix_trim_len(
        read_seq,
        config.min_tso_match_len,
        config.max_tso_mismatches,
    )
}

fn trimmed_end(read: &[u8], qual: Option<&[u8]>, config: ScoreConfig) -> (usize, u16) {
    let mut end = read.len();
    let mut flags = 0_u16;
    if config.trim_low_quality_tail {
        while end > 0
            && qual
                .and_then(|q| q.get(end - 1))
                .is_some_and(|&q| q.saturating_sub(33) < config.min_tail_phred)
        {
            end -= 1;
            flags |= SCORE_FLAG_TRIMMED_LOW_QUALITY;
        }
    }
    while end > 0 {
        match read[end - 1] {
            b'A' | b'a' if config.trim_poly_a => {
                end -= 1;
                flags |= SCORE_FLAG_TRIMMED_POLY_A;
            }
            b'T' | b't' if config.trim_poly_t => {
                end -= 1;
                flags |= SCORE_FLAG_TRIMMED_POLY_T;
            }
            _ => break,
        }
    }
    (end, flags)
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
    use crate::candidate::{CandidateBucketKey, CandidateHit, CandidateLocusBucket};
    use crate::index::{TranscriptIndex, TranscriptMeta};
    use crate::score::{
        hamming_ascii, SCORE_FLAG_RIGHT_SOFTCLIP, SCORE_FLAG_TRIMMED_LOW_QUALITY,
        SCORE_FLAG_TRIMMED_POLY_A, SCORE_FLAG_TRIMMED_TSO,
    };

    #[test]
    fn hamming_exact_match_zero() {
        assert_eq!(hamming_ascii(b"ACGT", b"ACGT"), 0);
    }

    #[test]
    fn hamming_one_mismatch() {
        assert_eq!(hamming_ascii(b"ACGT", b"ACAT"), 1);
    }

    #[test]
    fn default_scoring_matches_old_behavior() {
        let index = index_with_transcript("ACGTACGT");
        let bucket = bucket_at(0);
        let reads = vec![(0, b"NNNNNNNN".to_vec())];
        let quals = vec![(0, b"IIIIIIII".to_vec())];
        let mut scored = Vec::new();
        ScalarScorer::default().score_bucket(&reads, &quals, &index, &bucket, &mut scored, None);
        assert!(scored.is_empty());
    }

    #[test]
    fn poly_a_tail_can_softclip() {
        let index = index_with_transcript("ACGTACGT");
        let bucket = bucket_at(0);
        let reads = vec![(0, b"ACGTACGTAAAA".to_vec())];
        let quals = vec![(0, b"IIIIIIIIIIII".to_vec())];
        let mut scored = Vec::new();
        let scorer = ScalarScorer {
            config: ScoreConfig {
                trim_poly_a: true,
                min_scored_len: 8,
                ..ScoreConfig::default()
            },
            ..ScalarScorer::default()
        };
        scorer.score_bucket(&reads, &quals, &index, &bucket, &mut scored, None);
        assert_eq!(scored.len(), 1);
        assert_ne!(scored[0].flags & SCORE_FLAG_TRIMMED_POLY_A, 0);
        assert_ne!(scored[0].flags & SCORE_FLAG_RIGHT_SOFTCLIP, 0);
    }

    #[test]
    fn low_quality_tail_can_trim() {
        let index = index_with_transcript("ACGTACGT");
        let bucket = bucket_at(0);
        let reads = vec![(0, b"ACGTACGTGGGG".to_vec())];
        let quals = vec![(0, b"IIIIIIII!!!!".to_vec())];
        let mut scored = Vec::new();
        let scorer = ScalarScorer {
            config: ScoreConfig {
                trim_low_quality_tail: true,
                min_scored_len: 8,
                ..ScoreConfig::default()
            },
            ..ScalarScorer::default()
        };
        scorer.score_bucket(&reads, &quals, &index, &bucket, &mut scored, None);
        assert_eq!(scored.len(), 1);
        assert_ne!(scored[0].flags & SCORE_FLAG_TRIMMED_LOW_QUALITY, 0);
    }

    #[test]
    fn tso_prefix_can_trim() {
        let index = index_with_transcript("ACGTACGT");
        let bucket = bucket_at(0);
        let reads = vec![(0, b"AAGCAGTGGTACGTACGT".to_vec())];
        let quals = vec![(0, b"IIIIIIIIIIIIIIIIII".to_vec())];
        let mut scored = Vec::new();
        let scorer = ScalarScorer {
            config: ScoreConfig {
                trim_tso: true,
                min_scored_len: 8,
                ..ScoreConfig::default()
            },
            ..ScalarScorer::default()
        };
        scorer.score_bucket(&reads, &quals, &index, &bucket, &mut scored, None);
        assert_eq!(scored.len(), 1);
        assert_eq!(scored[0].pos, 0);
        assert_ne!(scored[0].flags & SCORE_FLAG_TRIMMED_TSO, 0);
    }

    #[test]
    fn full_length_match_beats_softclipped_match() {
        let index = index_with_transcript("ACGTACGT");
        let bucket = bucket_at(0);
        let reads = vec![(0, b"ACGTACGT".to_vec())];
        let quals = vec![(0, b"IIIIIIII".to_vec())];
        let mut scored = Vec::new();
        let scorer = ScalarScorer {
            config: ScoreConfig {
                max_right_softclip: 4,
                min_scored_len: 4,
                ..ScoreConfig::default()
            },
            ..ScalarScorer::default()
        };
        scorer.score_bucket(&reads, &quals, &index, &bucket, &mut scored, None);
        assert_eq!(scored.len(), 1);
        assert_eq!(scored[0].score, 8);
        assert_eq!(scored[0].flags & SCORE_FLAG_RIGHT_SOFTCLIP, 0);
    }

    fn index_with_transcript(seq: &str) -> TranscriptIndex {
        TranscriptIndex {
            format_version: 1,
            k: 3,
            max_kmer_frequency: 256,
            genes: vec!["GENE".to_owned()],
            transcripts: vec![TranscriptMeta {
                transcript_id: 0,
                gene_id: 0,
                name: "tx".to_owned(),
                len: seq.len() as u32,
            }],
            transcript_sequences: vec![seq.to_owned()],
            kmers: vec![],
            postings: vec![],
        }
    }

    fn bucket_at(pos: u32) -> CandidateLocusBucket {
        CandidateLocusBucket {
            key: CandidateBucketKey {
                transcript_id: 0,
                pos_bin: 0,
                strand: 0,
            },
            hits: vec![CandidateHit {
                read_id: 0,
                transcript_id: 0,
                gene_id: 0,
                pos,
                strand: 0,
                seed_count: 1,
                seed_score: 1,
            }],
        }
    }
}

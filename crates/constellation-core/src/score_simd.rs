use crate::candidate::CandidateLocusBucket;
use crate::index::IndexAccess;
use crate::score::{CandidateScorer, ScoreBlock, ScoreConfig, ScoreFailureStats, ScoredCandidate};
use crate::score_scalar::ScalarScorer;
use crate::ReadId;

#[derive(Debug, Clone, Copy)]
pub struct PulpScorer {
    pub config: ScoreConfig,
}

impl Default for PulpScorer {
    fn default() -> Self {
        Self {
            config: ScoreConfig::default(),
        }
    }
}

impl CandidateScorer for PulpScorer {
    fn score_bucket(
        &self,
        reads: &[(ReadId, Vec<u8>)],
        quals: &[(ReadId, Vec<u8>)],
        index: &dyn IndexAccess,
        bucket: &CandidateLocusBucket,
        out: &mut Vec<ScoredCandidate>,
        read_score_stats: Option<&mut [ScoreFailureStats]>,
    ) -> ScoreFailureStats {
        if !self.config.uses_default_window() {
            return ScalarScorer {
                config: self.config,
            }
            .score_bucket(reads, quals, index, bucket, out, read_score_stats);
        }
        if self.config.library_strand != crate::score::LibraryStrand::Unstranded {
            return ScalarScorer {
                config: self.config,
            }
            .score_bucket(reads, quals, index, bucket, out, read_score_stats);
        }
        let before = out.len();
        let block = ScoreBlock::from_bucket(reads, index, bucket);
        let mut stats = ScoreFailureStats {
            candidates_seen: bucket.hits.len() as u64,
            candidates_out_of_bounds: bucket.hits.len().saturating_sub(block.read_ids.len()) as u64,
            ..ScoreFailureStats::default()
        };
        let mut read_score_stats = read_score_stats;
        for (idx, (read_seq, ref_seq)) in block.pairs().enumerate() {
            let mismatches = hamming_ascii_pulp(read_seq, ref_seq);
            if mismatches <= self.config.max_mismatches {
                out.push(ScoredCandidate {
                    read_id: block.read_ids[idx],
                    transcript_id: block.transcript_ids[idx],
                    gene_id: block.gene_ids[idx],
                    pos: block.pos[idx],
                    strand: 0,
                    mismatches,
                    score: read_seq.len().saturating_sub(mismatches as usize) as u16,
                    flags: crate::score::SCORE_FLAG_FULL_LENGTH,
                });
            } else {
                stats.candidates_failed_mismatch += 1;
                if let Some(read_score_stats) = read_score_stats.as_deref_mut() {
                    if let Some(slot) = read_score_stats.get_mut(block.read_ids[idx] as usize) {
                        slot.candidates_seen += 1;
                        slot.candidates_failed_mismatch += 1;
                    }
                }
            }
        }
        stats.candidates_passed_full_length = out.len().saturating_sub(before) as u64;
        stats
    }
}

#[cfg(feature = "simd")]
fn hamming_ascii_pulp(a: &[u8], b: &[u8]) -> u32 {
    use pulp::{Arch, Simd, WithSimd};

    struct Hamming<'a> {
        a: &'a [u8],
        b: &'a [u8],
    }

    impl WithSimd for Hamming<'_> {
        type Output = u32;

        #[inline(always)]
        fn with_simd<S: Simd>(self, simd: S) -> Self::Output {
            let shared = self.a.len().min(self.b.len());
            let mut mismatches = self.a.len().abs_diff(self.b.len()) as u32;
            let (a_head, a_tail) = S::as_simd_u8s(&self.a[..shared]);
            let (b_head, b_tail) = S::as_simd_u8s(&self.b[..shared]);

            for (&a_vec, &b_vec) in a_head.iter().zip(b_head.iter()) {
                let eq = simd.equal_u8s(a_vec, b_vec);
                let neq = simd.not_m8s(eq);
                let neq_bytes = simd.transmute_u8s_m8s(neq);
                let bytes = unsafe {
                    std::slice::from_raw_parts(
                        (&neq_bytes as *const S::u8s).cast::<u8>(),
                        S::U8_LANES,
                    )
                };
                mismatches += bytes.iter().filter(|&&byte| byte != 0).count() as u32;
            }

            mismatches
                + a_tail
                    .iter()
                    .zip(b_tail.iter())
                    .filter(|(left, right)| left != right)
                    .count() as u32
        }
    }

    Arch::new().dispatch(Hamming { a, b })
}

#[cfg(not(feature = "simd"))]
fn hamming_ascii_pulp(a: &[u8], b: &[u8]) -> u32 {
    crate::score::hamming_ascii(a, b)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::candidate::{CandidateBucketKey, CandidateHit, CandidateLocusBucket};
    use crate::index::{TranscriptIndex, TranscriptMeta};
    use crate::score::CandidateScorer;
    use crate::score_scalar::ScalarScorer;

    #[test]
    fn simd_matches_scalar() {
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
                pos: 0,
                strand: 0,
                seed_count: 1,
                seed_score: 1,
            }],
        };
        let reads = vec![(0, b"ACGT".to_vec())];
        let mut scalar = Vec::new();
        let mut simd = Vec::new();
        let quals = vec![(0, b"IIII".to_vec())];
        ScalarScorer::default().score_bucket(&reads, &quals, &index, &bucket, &mut scalar, None);
        PulpScorer::default().score_bucket(&reads, &quals, &index, &bucket, &mut simd, None);
        assert_eq!(scalar, simd);
    }
}

use crate::candidate::CandidateLocusBucket;
use crate::index::IndexAccess;
use crate::score::{CandidateScorer, ScoreBlock, ScoredCandidate};
use crate::ReadId;

#[derive(Debug, Clone, Copy)]
pub struct PulpScorer {
    pub max_mismatches: u32,
}

impl Default for PulpScorer {
    fn default() -> Self {
        Self { max_mismatches: 4 }
    }
}

impl CandidateScorer for PulpScorer {
    fn score_bucket(
        &self,
        reads: &[(ReadId, Vec<u8>)],
        index: &dyn IndexAccess,
        bucket: &CandidateLocusBucket,
        out: &mut Vec<ScoredCandidate>,
    ) {
        let block = ScoreBlock::from_bucket(reads, index, bucket);
        for (idx, (read_seq, ref_seq)) in block.pairs().enumerate() {
            let mismatches = hamming_ascii_pulp(read_seq, ref_seq);
            if mismatches <= self.max_mismatches {
                out.push(ScoredCandidate {
                    read_id: block.read_ids[idx],
                    transcript_id: block.transcript_ids[idx],
                    gene_id: block.gene_ids[idx],
                    pos: block.pos[idx],
                    mismatches,
                    score: read_seq.len().saturating_sub(mismatches as usize) as u16,
                    flags: 0,
                });
            }
        }
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
        ScalarScorer::default().score_bucket(&reads, &index, &bucket, &mut scalar);
        PulpScorer::default().score_bucket(&reads, &index, &bucket, &mut simd);
        assert_eq!(scalar, simd);
    }
}

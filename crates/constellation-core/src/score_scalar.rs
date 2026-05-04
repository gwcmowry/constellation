use crate::candidate::CandidateLocusBucket;
use crate::index::IndexAccess;
use crate::score::{
    hamming_ascii, hamming_revcomp_ascii, read_seq_by_id, CandidateScorer, ScoredCandidate,
};
use crate::ReadId;

#[derive(Debug, Clone, Copy)]
pub struct ScalarScorer {
    pub max_mismatches: u32,
}

impl Default for ScalarScorer {
    fn default() -> Self {
        Self { max_mismatches: 4 }
    }
}

impl CandidateScorer for ScalarScorer {
    fn score_bucket(
        &self,
        reads: &[(ReadId, Vec<u8>)],
        index: &dyn IndexAccess,
        bucket: &CandidateLocusBucket,
        out: &mut Vec<ScoredCandidate>,
    ) {
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
            let reference = &transcript_seq[start..end];
            let mismatches = if hit.strand == 0 {
                hamming_ascii(read_seq, reference)
            } else {
                hamming_revcomp_ascii(read_seq, reference)
            };
            if mismatches <= self.max_mismatches {
                out.push(ScoredCandidate {
                    read_id: hit.read_id,
                    transcript_id: hit.transcript_id,
                    gene_id: hit.gene_id,
                    pos: hit.pos,
                    mismatches,
                    score: read_seq.len().saturating_sub(mismatches as usize) as u16,
                    flags: 0,
                });
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use crate::score::hamming_ascii;

    #[test]
    fn hamming_exact_match_zero() {
        assert_eq!(hamming_ascii(b"ACGT", b"ACGT"), 0);
    }

    #[test]
    fn hamming_one_mismatch() {
        assert_eq!(hamming_ascii(b"ACGT", b"ACAT"), 1);
    }
}

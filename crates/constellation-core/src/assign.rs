use crate::score::ScoredCandidate;
use crate::score::{
    SCORE_FLAG_ANTISENSE, SCORE_FLAG_TARGET_EXON, SCORE_FLAG_TARGET_GENE_BODY,
    SCORE_FLAG_TARGET_INTRON,
};
use crate::{GeneId, ReadId, TranscriptId};
use std::collections::BTreeSet;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum AssignmentType {
    UniqueGene,
    AmbiguousGene,
    AmbiguousTranscriptSameGene,
    AntisenseGene,
    Unmapped,
    LowComplexity,
    LowQuality,
}

impl AssignmentType {
    pub fn as_str(self) -> &'static str {
        match self {
            Self::UniqueGene => "unique_gene",
            Self::AmbiguousGene => "ambiguous_gene",
            Self::AmbiguousTranscriptSameGene => "ambiguous_transcript_same_gene",
            Self::AntisenseGene => "antisense_gene",
            Self::Unmapped => "unmapped",
            Self::LowComplexity => "low_complexity",
            Self::LowQuality => "low_quality",
        }
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Assignment {
    pub read_id: ReadId,
    pub cell_barcode: String,
    pub umi: String,
    pub assignment_type: AssignmentType,
    pub gene_id: Option<GeneId>,
    pub transcript_id: Option<TranscriptId>,
    pub candidate_count: u32,
    pub score: u16,
    pub flags: u16,
}

impl Assignment {
    pub fn tsv_header() -> &'static str {
        "read_id\tcell_barcode\tumi\tassignment_type\tgene_id\ttranscript_id\tcandidate_count\tscore\tflags\n"
    }

    pub fn to_tsv_row(&self) -> String {
        format!(
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n",
            self.read_id,
            self.cell_barcode,
            self.umi,
            self.assignment_type.as_str(),
            self.gene_id
                .map_or_else(|| ".".to_owned(), |id| id.to_string()),
            self.transcript_id
                .map_or_else(|| ".".to_owned(), |id| id.to_string()),
            self.candidate_count,
            self.score,
            self.flags
        )
    }
}

pub fn assign_read(
    read_id: ReadId,
    cell_barcode: String,
    umi: String,
    candidates: &[ScoredCandidate],
) -> Assignment {
    if candidates.is_empty() {
        return Assignment {
            read_id,
            cell_barcode,
            umi,
            assignment_type: AssignmentType::Unmapped,
            gene_id: None,
            transcript_id: None,
            candidate_count: 0,
            score: 0,
            flags: 0,
        };
    }

    let best_score = candidates
        .iter()
        .map(|candidate| candidate.score)
        .max()
        .unwrap_or(0);
    let best: Vec<_> = candidates
        .iter()
        .filter(|candidate| candidate.score == best_score)
        .collect();
    let sense_best: Vec<_> = best
        .iter()
        .copied()
        .filter(|candidate| candidate.flags & SCORE_FLAG_ANTISENSE == 0)
        .collect();
    if !sense_best.is_empty() {
        let sense_best = target_class_priority_best(&sense_best);
        return assign_from_best(
            read_id,
            cell_barcode,
            umi,
            candidates.len() as u32,
            best_score,
            &sense_best,
        );
    }

    let genes: BTreeSet<_> = best.iter().map(|candidate| candidate.gene_id).collect();
    let gene_id = (genes.len() == 1)
        .then(|| genes.iter().next().copied())
        .flatten();
    return Assignment {
        read_id,
        cell_barcode,
        umi,
        assignment_type: AssignmentType::AntisenseGene,
        gene_id,
        transcript_id: None,
        candidate_count: candidates.len() as u32,
        score: best_score,
        flags: best_flags(&best),
    };
}

fn target_class_priority_best<'a>(best: &[&'a ScoredCandidate]) -> Vec<&'a ScoredCandidate> {
    let genes: BTreeSet<_> = best.iter().map(|candidate| candidate.gene_id).collect();
    if genes.len() != 1 {
        return best.to_vec();
    }
    let best_rank = best
        .iter()
        .map(|candidate| target_class_rank(candidate.flags))
        .max()
        .unwrap_or(0);
    best.iter()
        .copied()
        .filter(|candidate| target_class_rank(candidate.flags) == best_rank)
        .collect()
}

fn target_class_rank(flags: u16) -> u8 {
    if flags & SCORE_FLAG_TARGET_EXON != 0 {
        3
    } else if flags & SCORE_FLAG_TARGET_INTRON != 0 {
        2
    } else if flags & SCORE_FLAG_TARGET_GENE_BODY != 0 {
        1
    } else {
        0
    }
}

fn assign_from_best(
    read_id: ReadId,
    cell_barcode: String,
    umi: String,
    candidate_count: u32,
    best_score: u16,
    best: &[&ScoredCandidate],
) -> Assignment {
    let genes: BTreeSet<_> = best.iter().map(|candidate| candidate.gene_id).collect();
    let transcripts: BTreeSet<_> = best
        .iter()
        .map(|candidate| candidate.transcript_id)
        .collect();

    let (assignment_type, gene_id, transcript_id) = if genes.len() == 1 && transcripts.len() == 1 {
        (
            AssignmentType::UniqueGene,
            genes.iter().next().copied(),
            transcripts.iter().next().copied(),
        )
    } else if genes.len() == 1 {
        (
            AssignmentType::AmbiguousTranscriptSameGene,
            genes.iter().next().copied(),
            None,
        )
    } else {
        (AssignmentType::AmbiguousGene, None, None)
    };

    Assignment {
        read_id,
        cell_barcode,
        umi,
        assignment_type,
        gene_id,
        transcript_id,
        candidate_count,
        score: best_score,
        flags: best_flags(best),
    }
}

fn best_flags(best: &[&ScoredCandidate]) -> u16 {
    best.iter()
        .fold(0_u16, |flags, candidate| flags | candidate.flags)
}

pub fn flagged_assignment(
    read_id: ReadId,
    cell_barcode: String,
    umi: String,
    assignment_type: AssignmentType,
    flags: u16,
) -> Assignment {
    Assignment {
        read_id,
        cell_barcode,
        umi,
        assignment_type,
        gene_id: None,
        transcript_id: None,
        candidate_count: 0,
        score: 0,
        flags,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn candidate(read_id: ReadId, gene_id: GeneId, transcript_id: TranscriptId) -> ScoredCandidate {
        ScoredCandidate {
            read_id,
            gene_id,
            transcript_id,
            pos: 0,
            strand: 0,
            mismatches: 0,
            score: 10,
            flags: 0,
        }
    }

    fn candidate_with_flags(
        read_id: ReadId,
        gene_id: GeneId,
        transcript_id: TranscriptId,
        flags: u16,
    ) -> ScoredCandidate {
        ScoredCandidate {
            flags,
            ..candidate(read_id, gene_id, transcript_id)
        }
    }

    #[test]
    fn unique_gene_assignment() {
        let assignment = assign_read(0, "CB".to_owned(), "UMI".to_owned(), &[candidate(0, 1, 2)]);
        assert_eq!(assignment.assignment_type, AssignmentType::UniqueGene);
        assert_eq!(assignment.gene_id, Some(1));
    }

    #[test]
    fn ambiguous_gene_assignment() {
        let candidates = [candidate(0, 1, 2), candidate(0, 2, 3)];
        let assignment = assign_read(0, "CB".to_owned(), "UMI".to_owned(), &candidates);
        assert_eq!(assignment.assignment_type, AssignmentType::AmbiguousGene);
    }

    #[test]
    fn same_gene_multi_transcript_assignment_is_gene_countable() {
        let candidates = [candidate(0, 1, 2), candidate(0, 1, 3)];
        let assignment = assign_read(0, "CB".to_owned(), "UMI".to_owned(), &candidates);
        assert_eq!(
            assignment.assignment_type,
            AssignmentType::AmbiguousTranscriptSameGene
        );
        assert_eq!(assignment.gene_id, Some(1));
        assert_eq!(assignment.transcript_id, None);
    }

    #[test]
    fn sense_candidate_beats_equal_score_antisense_candidate() {
        let sense = candidate(0, 1, 2);
        let mut antisense = candidate(0, 2, 3);
        antisense.flags = SCORE_FLAG_ANTISENSE;
        let assignment = assign_read(0, "CB".to_owned(), "UMI".to_owned(), &[antisense, sense]);
        assert_eq!(assignment.assignment_type, AssignmentType::UniqueGene);
        assert_eq!(assignment.gene_id, Some(1));
    }

    #[test]
    fn antisense_only_best_candidates_are_classified() {
        let mut antisense = candidate(0, 1, 2);
        antisense.flags = SCORE_FLAG_ANTISENSE;
        let assignment = assign_read(0, "CB".to_owned(), "UMI".to_owned(), &[antisense]);
        assert_eq!(assignment.assignment_type, AssignmentType::AntisenseGene);
        assert_eq!(assignment.gene_id, Some(1));
    }

    #[test]
    fn exonic_same_gene_tie_beats_gene_body_tie() {
        let exon = candidate_with_flags(0, 1, 2, SCORE_FLAG_TARGET_EXON);
        let gene_body = candidate_with_flags(0, 1, 3, SCORE_FLAG_TARGET_GENE_BODY);
        let assignment = assign_read(0, "CB".to_owned(), "UMI".to_owned(), &[gene_body, exon]);
        assert_eq!(assignment.assignment_type, AssignmentType::UniqueGene);
        assert_eq!(assignment.transcript_id, Some(2));
        assert_ne!(assignment.flags & SCORE_FLAG_TARGET_EXON, 0);
    }

    #[test]
    fn target_priority_does_not_resolve_multi_gene_tie() {
        let exon = candidate_with_flags(0, 1, 2, SCORE_FLAG_TARGET_EXON);
        let gene_body = candidate_with_flags(0, 2, 3, SCORE_FLAG_TARGET_GENE_BODY);
        let assignment = assign_read(0, "CB".to_owned(), "UMI".to_owned(), &[gene_body, exon]);
        assert_eq!(assignment.assignment_type, AssignmentType::AmbiguousGene);
    }
}

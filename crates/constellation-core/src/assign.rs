use crate::score::ScoredCandidate;
use crate::{GeneId, ReadId, TranscriptId};
use std::collections::BTreeSet;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum AssignmentType {
    UniqueGene,
    AmbiguousGene,
    AmbiguousTranscriptSameGene,
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
        candidate_count: candidates.len() as u32,
        score: best_score,
        flags: 0,
    }
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
            mismatches: 0,
            score: 10,
            flags: 0,
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
}

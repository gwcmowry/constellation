use std::collections::{BTreeMap, HashMap};
use std::io;
use std::path::Path;

use crate::dna::{normalize_seq, revcomp};
use crate::io::fasta::{read_fasta, FastaRecord};
use crate::io::gtf::read_gtf_exons;

#[derive(Clone, Debug, Eq, PartialEq)]
pub struct TranscriptPath {
    pub name: String,
    pub sequence: Vec<u8>,
}

#[derive(Clone, Debug, Eq, PartialEq)]
pub struct Component {
    pub id: u32,
    pub gene_id: String,
    pub transcripts: Vec<TranscriptPath>,
}

impl Component {
    pub fn all_paths_mask(&self) -> u64 {
        if self.transcripts.len() >= 64 {
            u64::MAX
        } else {
            (1u64 << self.transcripts.len()) - 1
        }
    }

    pub fn transcript_names_for_mask(&self, mask: u64) -> Vec<String> {
        self.transcripts
            .iter()
            .enumerate()
            .filter(|(idx, _)| mask & (1u64 << idx) != 0)
            .map(|(_, transcript)| transcript.name.clone())
            .collect()
    }
}

pub fn components_from_transcript_fasta(path: impl AsRef<Path>) -> io::Result<Vec<Component>> {
    let records = read_fasta(path)?;
    components_from_fasta_records(records)
}

pub fn components_from_fasta_records(records: Vec<FastaRecord>) -> io::Result<Vec<Component>> {
    let mut by_gene: BTreeMap<String, Vec<TranscriptPath>> = BTreeMap::new();

    for record in records {
        let name = transcript_name(&record.header).to_string();
        let gene_id = parse_gene_id(&record.header).unwrap_or_else(|| name.clone());
        by_gene.entry(gene_id).or_default().push(TranscriptPath {
            name,
            sequence: normalize_seq(&record.seq),
        });
    }

    let mut components = Vec::new();
    for (gene_id, transcripts) in by_gene {
        push_component_chunks(&mut components, gene_id, transcripts);
    }

    Ok(components)
}

pub fn components_from_genome_gtf(
    genome_path: impl AsRef<Path>,
    gtf_path: impl AsRef<Path>,
) -> io::Result<Vec<Component>> {
    let genome_records = read_fasta(genome_path)?;
    let genome: HashMap<String, Vec<u8>> = genome_records
        .into_iter()
        .map(|record| {
            (
                transcript_name(&record.header).to_string(),
                normalize_seq(&record.seq),
            )
        })
        .collect();
    let exons = read_gtf_exons(gtf_path)?;

    let mut grouped: BTreeMap<String, BTreeMap<String, Vec<_>>> = BTreeMap::new();
    for exon in exons {
        grouped
            .entry(exon.gene_id.clone())
            .or_default()
            .entry(exon.transcript_id.clone())
            .or_default()
            .push(exon);
    }

    let mut components = Vec::new();
    for (component_id, (gene_id, tx_exons)) in grouped.into_iter().enumerate() {
        let mut transcripts = Vec::new();
        for (transcript_id, mut exons) in tx_exons {
            exons.sort_by_key(|exon| exon.start);
            let strand = exons.first().map(|exon| exon.strand).unwrap_or('+');
            if strand == '-' {
                exons.reverse();
            }
            let mut sequence = Vec::new();
            for exon in exons {
                let chrom = genome.get(&exon.seqname).ok_or_else(|| {
                    io::Error::new(
                        io::ErrorKind::InvalidData,
                        format!("GTF references missing genome sequence {}", exon.seqname),
                    )
                })?;
                if exon.start == 0 || exon.end > chrom.len() || exon.start > exon.end {
                    return Err(io::Error::new(
                        io::ErrorKind::InvalidData,
                        format!(
                            "invalid exon coordinates {}:{}-{}",
                            exon.seqname, exon.start, exon.end
                        ),
                    ));
                }
                sequence.extend_from_slice(&chrom[exon.start - 1..exon.end]);
            }
            if strand == '-' {
                sequence = revcomp(&sequence);
            }
            transcripts.push(TranscriptPath {
                name: transcript_id,
                sequence,
            });
        }
        let _ = component_id;
        push_component_chunks(&mut components, gene_id, transcripts);
    }

    Ok(components)
}

fn push_component_chunks(
    components: &mut Vec<Component>,
    gene_id: String,
    transcripts: Vec<TranscriptPath>,
) {
    let chunk_count = transcripts.len().div_ceil(64);
    for (chunk_idx, chunk) in transcripts.chunks(64).enumerate() {
        let component_gene_id = if chunk_count > 1 {
            format!("{gene_id}#chunk{}", chunk_idx + 1)
        } else {
            gene_id.clone()
        };
        components.push(Component {
            id: components.len() as u32,
            gene_id: component_gene_id,
            transcripts: chunk.to_vec(),
        });
    }
}

pub fn toy_exon_skipping_component() -> Component {
    let exon1 = b"AAAACCCCGGGG";
    let exon2 = b"TTTTGGGGCCCC";
    let exon3 = b"ACACACACTATA";
    Component {
        id: 0,
        gene_id: "TOY_SKIP".to_string(),
        transcripts: vec![
            TranscriptPath {
                name: "T1_include".to_string(),
                sequence: [exon1.as_slice(), exon2.as_slice(), exon3.as_slice()].concat(),
            },
            TranscriptPath {
                name: "T2_skip".to_string(),
                sequence: [exon1.as_slice(), exon3.as_slice()].concat(),
            },
        ],
    }
}

pub fn parse_gene_id(header: &str) -> Option<String> {
    for token in header.split(|c: char| c.is_whitespace() || c == '|') {
        if let Some(value) = token.strip_prefix("gene=") {
            return Some(value.to_string());
        }
        if let Some(value) = token.strip_prefix("gene_id=") {
            return Some(value.to_string());
        }
        if let Some(value) = token.strip_prefix("gene:") {
            return Some(value.to_string());
        }
    }

    let marker = "gene_id \"";
    let start = header.find(marker)? + marker.len();
    let end = header[start..].find('"')? + start;
    Some(header[start..end].to_string())
}

pub fn transcript_name(header: &str) -> &str {
    header
        .split_whitespace()
        .next()
        .unwrap_or(header)
        .split('|')
        .next()
        .unwrap_or(header)
}

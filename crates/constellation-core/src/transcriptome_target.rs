use crate::dna::Base;
use crate::gtf::parse_attributes;
use flate2::read::MultiGzDecoder;
use rustc_hash::FxHashMap;
use std::cmp::Reverse;
use std::fs::File;
use std::io::{self, BufRead, BufReader, BufWriter, Write};
use std::path::Path;
use thiserror::Error;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum TranscriptomeTargetKind {
    ExonTranscripts,
    GeneBodies,
    IntronsOnly,
    ExonPlusGeneBody,
}

#[derive(Debug, Error)]
pub enum TranscriptomeTargetError {
    #[error("target I/O error: {0}")]
    Io(#[from] io::Error),
    #[error("malformed FASTA: {0}")]
    MalformedFasta(String),
    #[error("malformed GTF line {line}: {message}")]
    MalformedGtf { line: usize, message: String },
    #[error("no exon transcripts were found in the GTF")]
    NoTranscripts,
    #[error("genome contig {0} referenced by the GTF was not found in FASTA")]
    MissingContig(String),
    #[error("GTF interval {seqname}:{start}-{end} is outside FASTA contig length {len}")]
    IntervalOutOfBounds {
        seqname: String,
        start: u32,
        end: u32,
        len: usize,
    },
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct TranscriptomeTargetStats {
    pub num_genome_contigs: usize,
    pub num_transcripts: usize,
    pub num_exons: usize,
    pub total_bases: usize,
}

#[derive(Debug, Clone)]
struct Exon {
    start: u32,
    end: u32,
}

#[derive(Debug, Clone)]
struct TranscriptExons {
    transcript_id: String,
    gene_id: String,
    seqname: String,
    strand: u8,
    exons: Vec<Exon>,
}

#[derive(Debug, Clone)]
struct GeneModel {
    gene_id: String,
    seqname: String,
    strand: u8,
    explicit_start: Option<u32>,
    explicit_end: Option<u32>,
    exon_start: u32,
    exon_end: u32,
    exons: Vec<Exon>,
}

impl GeneModel {
    fn body_interval(&self) -> Option<Exon> {
        match (self.explicit_start, self.explicit_end) {
            (Some(start), Some(end)) => Some(Exon { start, end }),
            _ if self.exon_start <= self.exon_end => Some(Exon {
                start: self.exon_start,
                end: self.exon_end,
            }),
            _ => None,
        }
    }
}

#[derive(Debug, Clone)]
struct GtfModels {
    transcripts: Vec<TranscriptExons>,
    genes: Vec<GeneModel>,
}

pub fn build_transcriptome_target(
    genome_fasta: impl AsRef<Path>,
    gtf: impl AsRef<Path>,
    out_fasta: impl AsRef<Path>,
    kind: TranscriptomeTargetKind,
) -> Result<TranscriptomeTargetStats, TranscriptomeTargetError> {
    let genome = read_genome_fasta(genome_fasta)?;
    let mut models = read_gtf_models(gtf)?;
    if models.transcripts.is_empty() {
        return Err(TranscriptomeTargetError::NoTranscripts);
    }
    models
        .transcripts
        .sort_by(|a, b| a.transcript_id.cmp(&b.transcript_id));
    models.genes.sort_by(|a, b| a.gene_id.cmp(&b.gene_id));

    let out = File::create(out_fasta)?;
    let mut writer = BufWriter::new(out);
    let mut stats = TranscriptomeTargetStats {
        num_genome_contigs: genome.len(),
        num_transcripts: 0,
        num_exons: 0,
        total_bases: 0,
    };

    match kind {
        TranscriptomeTargetKind::ExonTranscripts => {
            write_exon_transcript_targets(
                &genome,
                &mut writer,
                &mut models.transcripts,
                &mut stats,
            )?;
        }
        TranscriptomeTargetKind::GeneBodies => {
            write_gene_body_targets(&genome, &mut writer, &models.genes, &mut stats)?;
        }
        TranscriptomeTargetKind::IntronsOnly => {
            write_intron_targets(&genome, &mut writer, &models.genes, &mut stats)?;
        }
        TranscriptomeTargetKind::ExonPlusGeneBody => {
            write_exon_transcript_targets(
                &genome,
                &mut writer,
                &mut models.transcripts,
                &mut stats,
            )?;
            write_gene_body_targets(&genome, &mut writer, &models.genes, &mut stats)?;
        }
    }
    writer.flush()?;

    Ok(stats)
}

fn write_exon_transcript_targets(
    genome: &FxHashMap<String, Vec<u8>>,
    writer: &mut impl Write,
    transcripts: &mut [TranscriptExons],
    stats: &mut TranscriptomeTargetStats,
) -> Result<(), TranscriptomeTargetError> {
    for transcript in transcripts {
        if transcript.strand == b'-' {
            transcript.exons.sort_by_key(|exon| Reverse(exon.start));
        } else {
            transcript.exons.sort_by_key(|exon| exon.start);
        }
        let contig = genome
            .get(&transcript.seqname)
            .ok_or_else(|| TranscriptomeTargetError::MissingContig(transcript.seqname.clone()))?;

        let mut seq = Vec::new();
        for exon in &transcript.exons {
            let start = exon.start as usize;
            let end = exon.end as usize;
            if start == 0 || end < start || end > contig.len() {
                return Err(TranscriptomeTargetError::IntervalOutOfBounds {
                    seqname: transcript.seqname.clone(),
                    start: exon.start,
                    end: exon.end,
                    len: contig.len(),
                });
            }
            let slice = &contig[(start - 1)..end];
            if transcript.strand == b'-' {
                push_reverse_complement(slice, &mut seq);
            } else {
                seq.extend_from_slice(slice);
            }
        }

        writeln!(
            writer,
            ">{}|gene:{}|target:exon_transcript|contig:{}|strand:{}|start:{}|end:{}",
            transcript.transcript_id,
            transcript.gene_id,
            transcript.seqname,
            transcript.strand as char,
            transcript
                .exons
                .iter()
                .map(|exon| exon.start)
                .min()
                .unwrap_or(0),
            transcript
                .exons
                .iter()
                .map(|exon| exon.end)
                .max()
                .unwrap_or(0)
        )?;
        write_wrapped_fasta(writer, &seq, 80)?;
        stats.num_transcripts += 1;
        stats.num_exons += transcript.exons.len();
        stats.total_bases += seq.len();
    }
    Ok(())
}

fn write_gene_body_targets(
    genome: &FxHashMap<String, Vec<u8>>,
    writer: &mut impl Write,
    genes: &[GeneModel],
    stats: &mut TranscriptomeTargetStats,
) -> Result<(), TranscriptomeTargetError> {
    for gene in genes {
        let Some(body) = gene.body_interval() else {
            continue;
        };
        let contig = genome
            .get(&gene.seqname)
            .ok_or_else(|| TranscriptomeTargetError::MissingContig(gene.seqname.clone()))?;
        let mut seq = Vec::new();
        append_interval_sequence(contig, &gene.seqname, gene.strand, &body, &mut seq)?;
        writeln!(
            writer,
            ">{}|gene:{}|target:gene_body|contig:{}|strand:{}|start:{}|end:{}",
            gene.gene_id, gene.gene_id, gene.seqname, gene.strand as char, body.start, body.end
        )?;
        write_wrapped_fasta(writer, &seq, 80)?;
        stats.num_transcripts += 1;
        stats.num_exons += 1;
        stats.total_bases += seq.len();
    }
    Ok(())
}

fn write_intron_targets(
    genome: &FxHashMap<String, Vec<u8>>,
    writer: &mut impl Write,
    genes: &[GeneModel],
    stats: &mut TranscriptomeTargetStats,
) -> Result<(), TranscriptomeTargetError> {
    for gene in genes {
        let Some(body) = gene.body_interval() else {
            continue;
        };
        let introns = intron_intervals(&body, &gene.exons);
        if introns.is_empty() {
            continue;
        }
        let contig = genome
            .get(&gene.seqname)
            .ok_or_else(|| TranscriptomeTargetError::MissingContig(gene.seqname.clone()))?;
        let mut ordered_introns = introns;
        if gene.strand == b'-' {
            ordered_introns.sort_by_key(|intron| Reverse(intron.start));
        }
        let mut seq = Vec::new();
        for intron in &ordered_introns {
            append_interval_sequence(contig, &gene.seqname, gene.strand, intron, &mut seq)?;
        }
        writeln!(
            writer,
            ">{}|gene:{}|target:intron|contig:{}|strand:{}|start:{}|end:{}",
            gene.gene_id, gene.gene_id, gene.seqname, gene.strand as char, body.start, body.end
        )?;
        write_wrapped_fasta(writer, &seq, 80)?;
        stats.num_transcripts += 1;
        stats.num_exons += ordered_introns.len();
        stats.total_bases += seq.len();
    }
    Ok(())
}

fn read_genome_fasta(
    path: impl AsRef<Path>,
) -> Result<FxHashMap<String, Vec<u8>>, TranscriptomeTargetError> {
    let mut reader = needletail::parse_fastx_file(path.as_ref())
        .map_err(|err| TranscriptomeTargetError::MalformedFasta(err.to_string()))?;
    let mut genome = FxHashMap::default();
    while let Some(record) = reader.next() {
        let record =
            record.map_err(|err| TranscriptomeTargetError::MalformedFasta(err.to_string()))?;
        let id = String::from_utf8_lossy(record.id())
            .split_whitespace()
            .next()
            .unwrap_or_default()
            .to_owned();
        let seq = record
            .seq()
            .iter()
            .map(|base| base.to_ascii_uppercase())
            .collect();
        genome.insert(id, seq);
    }
    if genome.is_empty() {
        return Err(TranscriptomeTargetError::MalformedFasta(
            "no records found".to_owned(),
        ));
    }
    Ok(genome)
}

fn read_gtf_models(path: impl AsRef<Path>) -> Result<GtfModels, TranscriptomeTargetError> {
    let reader = open_text_maybe_gz(path)?;
    let mut transcripts_by_id: FxHashMap<String, TranscriptExons> = FxHashMap::default();
    let mut genes_by_id: FxHashMap<String, GeneModel> = FxHashMap::default();

    for (idx, line) in reader.lines().enumerate() {
        let line_number = idx + 1;
        let line = line?;
        if line.trim().is_empty() || line.starts_with('#') {
            continue;
        }
        let fields: Vec<_> = line.split('\t').collect();
        if fields.len() < 9 {
            return Err(TranscriptomeTargetError::MalformedGtf {
                line: line_number,
                message: "expected at least 9 tab-delimited fields".to_owned(),
            });
        }
        if fields[2] != "exon" && fields[2] != "gene" {
            continue;
        }

        let start =
            fields[3]
                .parse::<u32>()
                .map_err(|_| TranscriptomeTargetError::MalformedGtf {
                    line: line_number,
                    message: format!("invalid start coordinate {}", fields[3]),
                })?;
        let end = fields[4]
            .parse::<u32>()
            .map_err(|_| TranscriptomeTargetError::MalformedGtf {
                line: line_number,
                message: format!("invalid end coordinate {}", fields[4]),
            })?;
        let strand = match fields[6].as_bytes().first().copied() {
            Some(b'+') => b'+',
            Some(b'-') => b'-',
            _ => {
                return Err(TranscriptomeTargetError::MalformedGtf {
                    line: line_number,
                    message: format!("unsupported strand {}", fields[6]),
                });
            }
        };

        let attrs = parse_attributes(fields[8]);
        let Some(gene_id) = attrs.get("gene_id") else {
            continue;
        };
        let gene_name = gene_id.clone();
        let seqname = fields[0].to_owned();

        let gene = genes_by_id
            .entry(gene_name.clone())
            .or_insert_with(|| GeneModel {
                gene_id: gene_name.clone(),
                seqname: seqname.clone(),
                strand,
                explicit_start: None,
                explicit_end: None,
                exon_start: u32::MAX,
                exon_end: 0,
                exons: Vec::new(),
            });
        if gene.seqname != seqname || gene.strand != strand {
            return Err(TranscriptomeTargetError::MalformedGtf {
                line: line_number,
                message: format!("gene {} has inconsistent seqname or strand", gene.gene_id),
            });
        }

        if fields[2] == "gene" {
            gene.explicit_start = Some(start);
            gene.explicit_end = Some(end);
            continue;
        }

        gene.exon_start = gene.exon_start.min(start);
        gene.exon_end = gene.exon_end.max(end);
        gene.exons.push(Exon { start, end });

        let Some(transcript_id) = attrs.get("transcript_id") else {
            continue;
        };
        let transcript_name = versioned_id(transcript_id, attrs.get("transcript_version"));

        let entry = transcripts_by_id
            .entry(transcript_name.clone())
            .or_insert_with(|| TranscriptExons {
                transcript_id: transcript_name,
                gene_id: gene_name,
                seqname: seqname.clone(),
                strand,
                exons: Vec::new(),
            });
        if entry.seqname != seqname || entry.strand != strand {
            return Err(TranscriptomeTargetError::MalformedGtf {
                line: line_number,
                message: format!(
                    "transcript {} has inconsistent seqname or strand",
                    entry.transcript_id
                ),
            });
        }
        entry.exons.push(Exon { start, end });
    }

    Ok(GtfModels {
        transcripts: transcripts_by_id.into_values().collect(),
        genes: genes_by_id.into_values().collect(),
    })
}

fn open_text_maybe_gz(
    path: impl AsRef<Path>,
) -> Result<Box<dyn BufRead>, TranscriptomeTargetError> {
    let path = path.as_ref();
    let file = File::open(path)?;
    if path.extension().is_some_and(|extension| extension == "gz") {
        Ok(Box::new(BufReader::new(MultiGzDecoder::new(file))))
    } else {
        Ok(Box::new(BufReader::new(file)))
    }
}

fn versioned_id(id: &str, version: Option<&String>) -> String {
    match version {
        Some(version) if !id.contains('.') => format!("{id}.{version}"),
        _ => id.to_owned(),
    }
}

fn push_reverse_complement(seq: &[u8], out: &mut Vec<u8>) {
    out.reserve(seq.len());
    for &base in seq.iter().rev() {
        let complemented = Base::from_ascii(base)
            .map(|base| base.complement().to_ascii())
            .unwrap_or(b'N');
        out.push(complemented);
    }
}

fn append_interval_sequence(
    contig: &[u8],
    seqname: &str,
    strand: u8,
    interval: &Exon,
    out: &mut Vec<u8>,
) -> Result<(), TranscriptomeTargetError> {
    let start = interval.start as usize;
    let end = interval.end as usize;
    if start == 0 || end < start || end > contig.len() {
        return Err(TranscriptomeTargetError::IntervalOutOfBounds {
            seqname: seqname.to_owned(),
            start: interval.start,
            end: interval.end,
            len: contig.len(),
        });
    }
    let slice = &contig[(start - 1)..end];
    if strand == b'-' {
        push_reverse_complement(slice, out);
    } else {
        out.extend_from_slice(slice);
    }
    Ok(())
}

fn intron_intervals(body: &Exon, exons: &[Exon]) -> Vec<Exon> {
    let mut merged_exons = exons
        .iter()
        .filter_map(|exon| {
            let start = exon.start.max(body.start);
            let end = exon.end.min(body.end);
            (start <= end).then_some(Exon { start, end })
        })
        .collect::<Vec<_>>();
    merged_exons.sort_by_key(|exon| exon.start);

    let mut merged = Vec::<Exon>::new();
    for exon in merged_exons {
        match merged.last_mut() {
            Some(last) if exon.start <= last.end.saturating_add(1) => {
                last.end = last.end.max(exon.end);
            }
            _ => merged.push(exon),
        }
    }

    let mut introns = Vec::new();
    let mut cursor = body.start;
    for exon in merged {
        if cursor < exon.start {
            introns.push(Exon {
                start: cursor,
                end: exon.start - 1,
            });
        }
        cursor = cursor.max(exon.end.saturating_add(1));
    }
    if cursor <= body.end {
        introns.push(Exon {
            start: cursor,
            end: body.end,
        });
    }
    introns
}

fn write_wrapped_fasta(writer: &mut impl Write, seq: &[u8], width: usize) -> Result<(), io::Error> {
    for chunk in seq.chunks(width) {
        writer.write_all(chunk)?;
        writer.write_all(b"\n")?;
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Write;

    #[test]
    fn builds_exon_transcripts_with_strand_and_versions() {
        let mut genome = tempfile::NamedTempFile::new().unwrap();
        writeln!(genome, ">1\nAACCGGTTAACCGGTT").unwrap();
        let mut gtf = tempfile::NamedTempFile::new().unwrap();
        writeln!(
            gtf,
            "1\ttest\texon\t2\t4\t.\t+\t.\tgene_id \"GENE1\"; transcript_id \"TX1\"; transcript_version \"3\";"
        )
        .unwrap();
        writeln!(
            gtf,
            "1\ttest\texon\t9\t12\t.\t+\t.\tgene_id \"GENE1\"; transcript_id \"TX1\"; transcript_version \"3\";"
        )
        .unwrap();
        writeln!(
            gtf,
            "1\ttest\texon\t1\t4\t.\t-\t.\tgene_id \"GENE2\"; transcript_id \"TX2\";"
        )
        .unwrap();
        writeln!(
            gtf,
            "1\ttest\texon\t9\t10\t.\t-\t.\tgene_id \"GENE2\"; transcript_id \"TX2\";"
        )
        .unwrap();
        let out = tempfile::NamedTempFile::new().unwrap();

        let stats = build_transcriptome_target(
            genome.path(),
            gtf.path(),
            out.path(),
            TranscriptomeTargetKind::ExonTranscripts,
        )
        .unwrap();
        let text = std::fs::read_to_string(out.path()).unwrap();

        assert_eq!(stats.num_transcripts, 2);
        assert_eq!(stats.num_exons, 4);
        assert!(text.contains(
            ">TX1.3|gene:GENE1|target:exon_transcript|contig:1|strand:+|start:2|end:12\nACCAACC"
        ));
        assert!(text.contains(
            ">TX2|gene:GENE2|target:exon_transcript|contig:1|strand:-|start:1|end:10\nTTGGTT"
        ));
    }

    #[test]
    fn builds_gene_body_and_intron_targets() {
        let mut genome = tempfile::NamedTempFile::new().unwrap();
        writeln!(genome, ">1\nAACCGGTTAACCGGTT").unwrap();
        let mut gtf = tempfile::NamedTempFile::new().unwrap();
        writeln!(gtf, "1\ttest\tgene\t2\t15\t.\t+\t.\tgene_id \"GENE1\";").unwrap();
        writeln!(
            gtf,
            "1\ttest\texon\t2\t4\t.\t+\t.\tgene_id \"GENE1\"; transcript_id \"TX1\";"
        )
        .unwrap();
        writeln!(
            gtf,
            "1\ttest\texon\t9\t12\t.\t+\t.\tgene_id \"GENE1\"; transcript_id \"TX1\";"
        )
        .unwrap();

        let gene_body_out = tempfile::NamedTempFile::new().unwrap();
        let gene_body_stats = build_transcriptome_target(
            genome.path(),
            gtf.path(),
            gene_body_out.path(),
            TranscriptomeTargetKind::GeneBodies,
        )
        .unwrap();
        let gene_body_text = std::fs::read_to_string(gene_body_out.path()).unwrap();
        assert_eq!(gene_body_stats.num_transcripts, 1);
        assert!(gene_body_text.contains(
            ">GENE1|gene:GENE1|target:gene_body|contig:1|strand:+|start:2|end:15\nACCGGTTAACCGGT"
        ));

        let intron_out = tempfile::NamedTempFile::new().unwrap();
        let intron_stats = build_transcriptome_target(
            genome.path(),
            gtf.path(),
            intron_out.path(),
            TranscriptomeTargetKind::IntronsOnly,
        )
        .unwrap();
        let intron_text = std::fs::read_to_string(intron_out.path()).unwrap();
        assert_eq!(intron_stats.num_transcripts, 1);
        assert_eq!(intron_stats.num_exons, 2);
        assert!(intron_text
            .contains(">GENE1|gene:GENE1|target:intron|contig:1|strand:+|start:2|end:15\nGGTTGGT"));
    }

    #[test]
    fn exon_plus_gene_body_emits_both_target_classes() {
        let mut genome = tempfile::NamedTempFile::new().unwrap();
        writeln!(genome, ">1\nAACCGGTTAACCGGTT").unwrap();
        let mut gtf = tempfile::NamedTempFile::new().unwrap();
        writeln!(
            gtf,
            "1\ttest\texon\t2\t4\t.\t+\t.\tgene_id \"GENE1\"; transcript_id \"TX1\";"
        )
        .unwrap();
        writeln!(
            gtf,
            "1\ttest\texon\t9\t12\t.\t+\t.\tgene_id \"GENE1\"; transcript_id \"TX1\";"
        )
        .unwrap();
        let out = tempfile::NamedTempFile::new().unwrap();

        let stats = build_transcriptome_target(
            genome.path(),
            gtf.path(),
            out.path(),
            TranscriptomeTargetKind::ExonPlusGeneBody,
        )
        .unwrap();
        let text = std::fs::read_to_string(out.path()).unwrap();

        assert_eq!(stats.num_transcripts, 2);
        assert!(text.contains("|target:exon_transcript|"));
        assert!(text.contains("|target:gene_body|"));
    }
}

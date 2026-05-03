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

pub fn build_transcriptome_target(
    genome_fasta: impl AsRef<Path>,
    gtf: impl AsRef<Path>,
    out_fasta: impl AsRef<Path>,
    kind: TranscriptomeTargetKind,
) -> Result<TranscriptomeTargetStats, TranscriptomeTargetError> {
    match kind {
        TranscriptomeTargetKind::ExonTranscripts => {
            build_exon_transcript_target(genome_fasta, gtf, out_fasta)
        }
    }
}

fn build_exon_transcript_target(
    genome_fasta: impl AsRef<Path>,
    gtf: impl AsRef<Path>,
    out_fasta: impl AsRef<Path>,
) -> Result<TranscriptomeTargetStats, TranscriptomeTargetError> {
    let genome = read_genome_fasta(genome_fasta)?;
    let mut transcripts = read_gtf_exons(gtf)?;
    if transcripts.is_empty() {
        return Err(TranscriptomeTargetError::NoTranscripts);
    }
    transcripts.sort_by(|a, b| a.transcript_id.cmp(&b.transcript_id));

    let out = File::create(out_fasta)?;
    let mut writer = BufWriter::new(out);
    let mut num_exons = 0_usize;
    let mut total_bases = 0_usize;

    for transcript in &mut transcripts {
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
            ">{}|gene:{}",
            transcript.transcript_id, transcript.gene_id
        )?;
        write_wrapped_fasta(&mut writer, &seq, 80)?;
        num_exons += transcript.exons.len();
        total_bases += seq.len();
    }
    writer.flush()?;

    Ok(TranscriptomeTargetStats {
        num_genome_contigs: genome.len(),
        num_transcripts: transcripts.len(),
        num_exons,
        total_bases,
    })
}

fn read_genome_fasta(
    path: impl AsRef<Path>,
) -> Result<FxHashMap<String, Vec<u8>>, TranscriptomeTargetError> {
    let mut reader = needletail::parse_fastx_file(path.as_ref())
        .map_err(|err| TranscriptomeTargetError::MalformedFasta(err.to_string()))?;
    let mut genome = FxHashMap::default();
    while let Some(record) = reader.next() {
        let record = record.map_err(|err| TranscriptomeTargetError::MalformedFasta(err.to_string()))?;
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

fn read_gtf_exons(
    path: impl AsRef<Path>,
) -> Result<Vec<TranscriptExons>, TranscriptomeTargetError> {
    let reader = open_text_maybe_gz(path)?;
    let mut transcripts_by_id: FxHashMap<String, TranscriptExons> = FxHashMap::default();

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
        if fields[2] != "exon" {
            continue;
        }

        let start = fields[3]
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
        let Some(transcript_id) = attrs.get("transcript_id") else {
            continue;
        };
        let Some(gene_id) = attrs.get("gene_id") else {
            continue;
        };
        let transcript_name = versioned_id(transcript_id, attrs.get("transcript_version"));
        let gene_name = gene_id.clone();
        let seqname = fields[0].to_owned();

        let entry =
            transcripts_by_id
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

    Ok(transcripts_by_id.into_values().collect())
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

fn write_wrapped_fasta(
    writer: &mut impl Write,
    seq: &[u8],
    width: usize,
) -> Result<(), io::Error> {
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
        assert!(text.contains(">TX1.3|gene:GENE1\nACCAACC"));
        assert!(text.contains(">TX2|gene:GENE2\nTTGGTT"));
    }
}

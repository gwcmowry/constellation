use crate::dna::{encode_acgt, iter_kmers_2bit};
use crate::index::{
    CompactIndexTempFiles, CompactTranscriptBuildRecord, KmerEntry, Posting, TranscriptIndex,
    TranscriptMeta,
};
use crate::{GeneId, TranscriptId};
use rayon::prelude::*;
use rustc_hash::{FxHashMap, FxHashSet};
use std::cmp::{Ordering, Reverse};
use std::collections::{BTreeMap, BinaryHeap};
use std::fs::{self, File};
use std::io::{self, BufReader, BufWriter, Read, Write};
use std::path::{Path, PathBuf};
use thiserror::Error;

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct FastaRecord {
    pub header: String,
    pub name: String,
    pub gene: String,
    pub seq: Vec<u8>,
}

#[derive(Debug, Error)]
pub enum BuildIndexError {
    #[error("FASTA I/O error: {0}")]
    Io(#[from] std::io::Error),
    #[error("malformed FASTA: {0}")]
    Malformed(String),
    #[error("k must be between 1 and 32, got {0}")]
    InvalidK(u8),
    #[error("index I/O error: {0}")]
    IndexIo(#[from] crate::index::IndexIoError),
    #[error("streaming index cannot represent more than u32::MAX postings")]
    TooManyPostings,
}

pub fn read_fasta(path: impl AsRef<Path>) -> Result<Vec<FastaRecord>, BuildIndexError> {
    read_fasta_needletail(path)
}

fn read_fasta_needletail(path: impl AsRef<Path>) -> Result<Vec<FastaRecord>, BuildIndexError> {
    let mut reader = needletail::parse_fastx_file(path.as_ref())
        .map_err(|err| BuildIndexError::Malformed(err.to_string()))?;
    let mut records = Vec::new();
    while let Some(record) = reader.next() {
        let record = record.map_err(|err| BuildIndexError::Malformed(err.to_string()))?;
        let header = String::from_utf8_lossy(record.id()).into_owned();
        let seq = record
            .seq()
            .iter()
            .map(|base| base.to_ascii_uppercase())
            .collect();
        records.push(make_record(header, seq));
    }

    if records.is_empty() {
        return Err(BuildIndexError::Malformed("no records found".to_owned()));
    }
    Ok(records)
}

pub fn build_transcript_index(
    fasta: impl AsRef<Path>,
    k: u8,
    max_kmer_frequency: u32,
) -> Result<TranscriptIndex, BuildIndexError> {
    build_transcript_index_with_gene_map(fasta, k, max_kmer_frequency, None)
}

pub fn build_transcript_index_with_gene_map(
    fasta: impl AsRef<Path>,
    k: u8,
    max_kmer_frequency: u32,
    transcript_gene_map: Option<&FxHashMap<String, String>>,
) -> Result<TranscriptIndex, BuildIndexError> {
    if !(1..=32).contains(&k) {
        return Err(BuildIndexError::InvalidK(k));
    }

    let records = read_fasta(fasta)?;
    let mut gene_ids: FxHashMap<String, GeneId> = FxHashMap::default();
    let mut genes = Vec::new();
    let mut transcripts = Vec::new();
    let mut transcript_sequences = Vec::new();
    let mut postings_by_kmer: BTreeMap<u64, Vec<Posting>> = BTreeMap::new();

    for (idx, record) in records.iter().enumerate() {
        let gene_name = transcript_gene_map
            .and_then(|map| map.get(&record.name))
            .cloned()
            .unwrap_or_else(|| record.gene.clone());
        let gene_id = match gene_ids.get(&gene_name) {
            Some(&id) => id,
            None => {
                let id = genes.len() as GeneId;
                gene_ids.insert(gene_name.clone(), id);
                genes.push(gene_name);
                id
            }
        };
        let transcript_id = idx as TranscriptId;
        transcripts.push(TranscriptMeta {
            transcript_id,
            gene_id,
            name: record.name.clone(),
            len: record.seq.len() as u32,
        });
        transcript_sequences.push(String::from_utf8_lossy(&record.seq).into_owned());

        let encoded = encode_acgt(&record.seq);
        for kmer in iter_kmers_2bit(&encoded, k) {
            postings_by_kmer
                .entry(kmer.code)
                .or_default()
                .push(Posting {
                    transcript_id,
                    gene_id,
                    pos: kmer.pos,
                    strand: 0,
                });
        }
    }

    let mut kmers = Vec::with_capacity(postings_by_kmer.len());
    let mut postings = Vec::new();
    for (code, mut kmer_postings) in postings_by_kmer {
        kmer_postings
            .par_sort_by_key(|posting| (posting.transcript_id, posting.pos, posting.strand));
        let start = postings.len() as u64;
        let len = kmer_postings.len() as u32;
        let (transcript_df, gene_df) = kmer_document_frequencies(&kmer_postings);
        postings.extend(kmer_postings);
        kmers.push(KmerEntry {
            kmer_code: code,
            postings_start: start,
            postings_len: len,
            freq_class: if len > max_kmer_frequency { 1 } else { 0 },
            transcript_df,
            gene_df,
        });
    }

    Ok(TranscriptIndex {
        format_version: 1,
        k,
        max_kmer_frequency,
        genes,
        transcripts,
        transcript_sequences,
        kmers,
        postings,
    })
}

pub fn build_compact_transcript_index_streaming(
    fasta: impl AsRef<Path>,
    k: u8,
    max_kmer_frequency: u32,
    out_path: impl AsRef<Path>,
) -> Result<(), BuildIndexError> {
    if !(1..=32).contains(&k) {
        return Err(BuildIndexError::InvalidK(k));
    }

    let out_path = out_path.as_ref();
    let temp_dir = streaming_temp_dir(out_path);
    fs::create_dir(&temp_dir)?;
    let result = build_compact_transcript_index_streaming_inner(
        fasta.as_ref(),
        k,
        max_kmer_frequency,
        out_path,
        &temp_dir,
    );
    let cleanup = fs::remove_dir_all(&temp_dir);
    match (result, cleanup) {
        (Ok(()), Ok(())) => Ok(()),
        (Err(err), _) => Err(err),
        (Ok(()), Err(err)) => Err(BuildIndexError::Io(err)),
    }
}

fn build_compact_transcript_index_streaming_inner(
    fasta: &Path,
    k: u8,
    max_kmer_frequency: u32,
    out_path: &Path,
    temp_dir: &Path,
) -> Result<(), BuildIndexError> {
    const RUN_RECORD_LIMIT: usize = 8_000_000;

    let sequences_path = temp_dir.join("sequences.bin");
    let kmers_path = temp_dir.join("kmers.bin");
    let gene_df_path = temp_dir.join("gene_df.bin");
    let postings_path = temp_dir.join("postings.bin");
    let mut sequences = BufWriter::new(File::create(&sequences_path)?);
    let mut reader = needletail::parse_fastx_file(fasta)
        .map_err(|err| BuildIndexError::Malformed(err.to_string()))?;
    let mut gene_ids: FxHashMap<String, GeneId> = FxHashMap::default();
    let mut genes = Vec::new();
    let mut transcripts = Vec::new();
    let mut run_records = Vec::with_capacity(RUN_RECORD_LIMIT);
    let mut run_paths = Vec::new();
    let mut sequences_len = 0_u64;

    while let Some(record) = reader.next() {
        let record = record.map_err(|err| BuildIndexError::Malformed(err.to_string()))?;
        let header = String::from_utf8_lossy(record.id()).into_owned();
        let seq: Vec<u8> = record
            .seq()
            .iter()
            .map(|base| base.to_ascii_uppercase())
            .collect();
        let record = make_record(header, seq);
        let gene_id = match gene_ids.get(&record.gene) {
            Some(&id) => id,
            None => {
                let id = genes.len() as GeneId;
                gene_ids.insert(record.gene.clone(), id);
                genes.push(record.gene.clone());
                id
            }
        };
        let transcript_id = transcripts.len() as TranscriptId;
        let seq_start = sequences_len;
        sequences.write_all(&record.seq)?;
        sequences_len += record.seq.len() as u64;
        transcripts.push(CompactTranscriptBuildRecord {
            transcript_id,
            gene_id,
            name: record.name,
            seq_start,
            len: record.seq.len() as u32,
        });

        let encoded = encode_acgt(&record.seq);
        for kmer in iter_kmers_2bit(&encoded, k) {
            run_records.push(RunPosting {
                kmer_code: kmer.code,
                transcript_id,
                pos_strand: kmer.pos << 1,
            });
            if run_records.len() >= RUN_RECORD_LIMIT {
                flush_run(temp_dir, &mut run_paths, &mut run_records)?;
            }
        }
    }
    sequences.flush()?;
    if transcripts.is_empty() {
        return Err(BuildIndexError::Malformed("no records found".to_owned()));
    }
    flush_run(temp_dir, &mut run_paths, &mut run_records)?;

    let transcript_gene_ids: Vec<_> = transcripts.iter().map(|record| record.gene_id).collect();
    let (num_kmers, num_postings) = merge_runs(
        &run_paths,
        &kmers_path,
        &gene_df_path,
        &postings_path,
        &transcript_gene_ids,
    )?;
    TranscriptIndex::save_compact_from_temp_files(
        out_path,
        k,
        max_kmer_frequency,
        &genes,
        &transcripts,
        num_kmers,
        num_postings,
        sequences_len,
        CompactIndexTempFiles {
            kmers_path: &kmers_path,
            gene_df_path: &gene_df_path,
            postings_path: &postings_path,
            sequences_path: &sequences_path,
        },
    )?;
    Ok(())
}

fn streaming_temp_dir(out_path: &Path) -> PathBuf {
    let file_name = out_path
        .file_name()
        .and_then(|name| name.to_str())
        .unwrap_or("constellation-index");
    let parent = out_path.parent().unwrap_or_else(|| Path::new("."));
    parent.join(format!(".{file_name}.tmp.{}", std::process::id()))
}

fn kmer_document_frequencies(postings: &[Posting]) -> (u32, u32) {
    let mut transcript_df = 0_u32;
    let mut last_transcript = None;
    let mut genes = FxHashSet::default();
    for posting in postings {
        if last_transcript != Some(posting.transcript_id) {
            transcript_df += 1;
            last_transcript = Some(posting.transcript_id);
        }
        genes.insert(posting.gene_id);
    }
    (transcript_df, genes.len().min(u32::MAX as usize) as u32)
}

fn flush_run(
    temp_dir: &Path,
    run_paths: &mut Vec<PathBuf>,
    records: &mut Vec<RunPosting>,
) -> Result<(), BuildIndexError> {
    if records.is_empty() {
        return Ok(());
    }
    records.par_sort_unstable();
    let path = temp_dir.join(format!("run_{:05}.bin", run_paths.len()));
    let mut out = BufWriter::new(File::create(&path)?);
    for record in records.iter().copied() {
        write_run_posting(&mut out, record)?;
    }
    out.flush()?;
    records.clear();
    run_paths.push(path);
    Ok(())
}

fn merge_runs(
    run_paths: &[PathBuf],
    kmers_path: &Path,
    gene_df_path: &Path,
    postings_path: &Path,
    transcript_gene_ids: &[GeneId],
) -> Result<(u64, u64), BuildIndexError> {
    let mut readers = Vec::with_capacity(run_paths.len());
    let mut heap = BinaryHeap::new();
    for path in run_paths {
        let mut reader = BufReader::new(File::open(path)?);
        if let Some(record) = read_run_posting(&mut reader)? {
            heap.push(Reverse(HeapRecord {
                record,
                run_idx: readers.len(),
            }));
        }
        readers.push(reader);
    }

    let mut kmers = BufWriter::new(File::create(kmers_path)?);
    let mut gene_dfs = BufWriter::new(File::create(gene_df_path)?);
    let mut postings = BufWriter::new(File::create(postings_path)?);
    let mut current_kmer = None::<u64>;
    let mut current_start = 0_u64;
    let mut current_len = 0_u32;
    let mut current_last_transcript = None::<TranscriptId>;
    let mut current_genes: FxHashSet<GeneId> = FxHashSet::default();
    let mut num_kmers = 0_u64;
    let mut num_postings = 0_u64;

    while let Some(Reverse(item)) = heap.pop() {
        if current_kmer != Some(item.record.kmer_code) {
            if let Some(kmer_code) = current_kmer {
                write_compact_kmer_entry(
                    &mut kmers,
                    &mut gene_dfs,
                    kmer_code,
                    current_start,
                    current_len,
                    current_genes.len().min(u32::MAX as usize) as u32,
                )?;
                num_kmers += 1;
            }
            current_kmer = Some(item.record.kmer_code);
            current_start = num_postings;
            current_len = 0;
            current_last_transcript = None;
            current_genes.clear();
        }
        if current_last_transcript != Some(item.record.transcript_id) {
            current_last_transcript = Some(item.record.transcript_id);
        }
        if let Some(&gene_id) = transcript_gene_ids.get(item.record.transcript_id as usize) {
            current_genes.insert(gene_id);
        }
        write_compact_posting(&mut postings, item.record)?;
        num_postings += 1;
        current_len = current_len
            .checked_add(1)
            .ok_or(BuildIndexError::TooManyPostings)?;

        if let Some(next) = read_run_posting(&mut readers[item.run_idx])? {
            heap.push(Reverse(HeapRecord {
                record: next,
                run_idx: item.run_idx,
            }));
        }
    }
    if let Some(kmer_code) = current_kmer {
        write_compact_kmer_entry(
            &mut kmers,
            &mut gene_dfs,
            kmer_code,
            current_start,
            current_len,
            current_genes.len().min(u32::MAX as usize) as u32,
        )?;
        num_kmers += 1;
    }
    kmers.flush()?;
    gene_dfs.flush()?;
    postings.flush()?;
    Ok((num_kmers, num_postings))
}

fn make_record(header: String, seq: Vec<u8>) -> FastaRecord {
    let name = header
        .split_whitespace()
        .next()
        .unwrap_or(&header)
        .split('|')
        .next()
        .unwrap_or(&header)
        .to_owned();
    let gene = header
        .split(['|', ' ', '\t'])
        .find_map(|part| {
            part.strip_prefix("gene=")
                .or_else(|| part.strip_prefix("gene:"))
                .or_else(|| part.strip_prefix("gene_id="))
                .or_else(|| part.strip_prefix("gene_id:"))
        })
        .map(clean_gene_token)
        .or_else(|| parse_quoted_gene_id(&header))
        .unwrap_or_else(|| name.clone());
    FastaRecord {
        header,
        name,
        gene,
        seq,
    }
}

fn clean_gene_token(gene: &str) -> String {
    gene.trim()
        .trim_matches('"')
        .trim_matches(';')
        .trim_matches('"')
        .to_owned()
}

fn parse_quoted_gene_id(header: &str) -> Option<String> {
    let rest = header.split_once("gene_id")?.1.trim_start();
    let rest = rest
        .strip_prefix('=')
        .or_else(|| rest.strip_prefix(':'))
        .unwrap_or(rest)
        .trim_start();
    if let Some(quoted) = rest.strip_prefix('"') {
        return quoted.split_once('"').map(|(gene, _)| gene.to_owned());
    }
    rest.split(['|', ' ', '\t', ';'])
        .find(|token| !token.is_empty())
        .map(clean_gene_token)
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
struct RunPosting {
    kmer_code: u64,
    transcript_id: u32,
    pos_strand: u32,
}

impl Ord for RunPosting {
    fn cmp(&self, other: &Self) -> Ordering {
        (self.kmer_code, self.transcript_id, self.pos_strand).cmp(&(
            other.kmer_code,
            other.transcript_id,
            other.pos_strand,
        ))
    }
}

impl PartialOrd for RunPosting {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
struct HeapRecord {
    record: RunPosting,
    run_idx: usize,
}

impl Ord for HeapRecord {
    fn cmp(&self, other: &Self) -> Ordering {
        self.record
            .cmp(&other.record)
            .then_with(|| self.run_idx.cmp(&other.run_idx))
    }
}

impl PartialOrd for HeapRecord {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

fn write_run_posting(out: &mut impl Write, record: RunPosting) -> io::Result<()> {
    out.write_all(&record.kmer_code.to_ne_bytes())?;
    out.write_all(&record.transcript_id.to_ne_bytes())?;
    out.write_all(&record.pos_strand.to_ne_bytes())?;
    Ok(())
}

fn read_run_posting(input: &mut impl Read) -> io::Result<Option<RunPosting>> {
    let mut bytes = [0_u8; 16];
    let mut read = 0;
    while read < bytes.len() {
        let n = input.read(&mut bytes[read..])?;
        if n == 0 {
            if read == 0 {
                return Ok(None);
            }
            return Err(io::Error::new(
                io::ErrorKind::UnexpectedEof,
                "truncated posting run",
            ));
        }
        read += n;
    }
    let mut kmer = [0_u8; 8];
    kmer.copy_from_slice(&bytes[..8]);
    let mut transcript = [0_u8; 4];
    transcript.copy_from_slice(&bytes[8..12]);
    let mut pos = [0_u8; 4];
    pos.copy_from_slice(&bytes[12..16]);
    Ok(Some(RunPosting {
        kmer_code: u64::from_ne_bytes(kmer),
        transcript_id: u32::from_ne_bytes(transcript),
        pos_strand: u32::from_ne_bytes(pos),
    }))
}

fn write_compact_kmer_entry(
    kmers: &mut impl Write,
    gene_dfs: &mut impl Write,
    kmer_code: u64,
    postings_start: u64,
    postings_len: u32,
    gene_df: u32,
) -> Result<(), BuildIndexError> {
    let postings_start: u32 = postings_start
        .try_into()
        .map_err(|_| BuildIndexError::TooManyPostings)?;
    kmers.write_all(&kmer_code.to_ne_bytes())?;
    kmers.write_all(&postings_start.to_ne_bytes())?;
    kmers.write_all(&postings_len.to_ne_bytes())?;
    gene_dfs.write_all(&(gene_df.min(u16::MAX as u32) as u16).to_ne_bytes())?;
    Ok(())
}

fn write_compact_posting(out: &mut impl Write, record: RunPosting) -> io::Result<()> {
    out.write_all(&record.transcript_id.to_ne_bytes())?;
    out.write_all(&record.pos_strand.to_ne_bytes())?;
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::dna::kmer_code_from_ascii;
    use std::io::Write;

    #[test]
    fn tiny_index_contains_expected_kmers() {
        let mut fasta = tempfile::NamedTempFile::new().unwrap();
        writeln!(
            fasta,
            ">tx1|gene=GENE_A\nACGTACGT\n>tx2|gene=GENE_B\nTTTACG"
        )
        .unwrap();
        let index = build_transcript_index(fasta.path(), 3, 256).unwrap();
        let acg = kmer_code_from_ascii(b"ACG").unwrap();
        let postings = index.lookup(acg);
        assert_eq!(index.genes, vec!["GENE_A", "GENE_B"]);
        assert_eq!(postings.len(), 3);
    }

    #[test]
    fn computes_gene_df_lower_than_raw_postings_for_isoforms() {
        let mut fasta = tempfile::NamedTempFile::new().unwrap();
        writeln!(
            fasta,
            ">tx1|gene=GENE_A\nAAAGGG\n>tx2|gene=GENE_A\nAAATTT\n>tx3|gene=GENE_A\nAAAGAA\n>tx4|gene=GENE_B\nCCCAAA"
        )
        .unwrap();
        let index = build_transcript_index(fasta.path(), 3, 256).unwrap();
        let aaa = kmer_code_from_ascii(b"AAA").unwrap();
        let stats = crate::index::IndexAccess::kmer_stats(&index, aaa).unwrap();
        assert_eq!(stats.raw_postings, 4);
        assert_eq!(stats.transcript_df, 4);
        assert_eq!(stats.gene_df, 2);
    }

    #[test]
    fn gtf_gene_map_overrides_fasta_gene() {
        let mut fasta = tempfile::NamedTempFile::new().unwrap();
        writeln!(fasta, ">tx1\nACGTACGT").unwrap();
        let mut map = FxHashMap::default();
        map.insert("tx1".to_owned(), "GENE_FROM_GTF".to_owned());
        let index = build_transcript_index_with_gene_map(fasta.path(), 3, 256, Some(&map)).unwrap();
        assert_eq!(index.genes, vec!["GENE_FROM_GTF"]);
    }

    #[test]
    fn parses_ensembl_cdna_gene_colon_header() {
        let record = make_record(
            "ENST1 cdna chromosome:GRCh38:1:1:10:1 gene:ENSG1 gene_symbol:ABC".to_owned(),
            b"ACGT".to_vec(),
        );
        assert_eq!(record.name, "ENST1");
        assert_eq!(record.gene, "ENSG1");
    }

    #[test]
    fn parses_gene_id_header_variants() {
        for (header, gene) in [
            ("TX1|gene_id=GENE1", "GENE1"),
            ("TX1|gene_id:GENE2", "GENE2"),
            ("TX1 gene_id \"GENE3\";", "GENE3"),
            ("TX1|gene:GENE4|target:gene_body", "GENE4"),
        ] {
            let record = make_record(header.to_owned(), b"ACGT".to_vec());
            assert_eq!(record.gene, gene);
        }
    }
}

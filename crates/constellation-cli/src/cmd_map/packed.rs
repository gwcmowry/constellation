//! Bounded reusable input slabs and fused gene-EC mapping/output workers.
use super::*;
use constellation_core::dna::{acgt_code, encode_and_check_complexity};
use constellation_core::ec_index::{EcConsensusScratch, EcResolvedMapping};
use std::sync::mpsc::Sender;
use std::sync::Mutex;

#[derive(Default)]
struct Slab {
    data: Vec<u8>,
    // Sequence start, sequence end, quality end. R1 qualities are validated but discarded.
    offsets: Vec<(usize, usize, usize)>,
}

impl Slab {
    fn fill(
        &mut self,
        reader: &mut FastqStreamReader,
        batch_size: usize,
        keep_quality: bool,
    ) -> Result<()> {
        self.data.clear();
        self.offsets.clear();
        for _ in 0..batch_size {
            if !reader.read_record()? {
                break;
            }
            let start = self.data.len();
            self.data
                .extend_from_slice(trim_fastq_line(reader.seq_buf.as_bytes()));
            self.data[start..].make_ascii_uppercase();
            let seq_end = self.data.len();
            if keep_quality {
                self.data
                    .extend_from_slice(trim_fastq_line(reader.qual_buf.as_bytes()));
            }
            self.offsets.push((start, seq_end, self.data.len()));
        }
        Ok(())
    }

    fn seq(&self, i: usize) -> &[u8] {
        let (start, end, _) = self.offsets[i];
        &self.data[start..end]
    }

    fn qual(&self, i: usize) -> &[u8] {
        let (_, start, end) = self.offsets[i];
        &self.data[start..end]
    }
}

struct ReaderPipe {
    rx: Option<Receiver<std::result::Result<Slab, String>>>,
    recycle: Option<Sender<Slab>>,
    handle: Option<JoinHandle<()>>,
}

impl ReaderPipe {
    fn new(path: PathBuf, batch_size: usize, keep_quality: bool) -> Self {
        let (tx, rx) = mpsc::sync_channel(2);
        let (recycle, recycled) = mpsc::channel();
        let handle = thread::spawn(move || {
            let run = || -> Result<()> {
                let mut reader = FastqStreamReader::open(&path)?;
                let mut allocated = 0;
                loop {
                    let mut slab = match recycled.try_recv() {
                        Ok(slab) => slab,
                        Err(_) if allocated < 3 => {
                            allocated += 1;
                            Slab::default()
                        }
                        Err(_) => match recycled.recv() {
                            Ok(slab) => slab,
                            Err(_) => return Ok(()),
                        },
                    };
                    slab.fill(&mut reader, batch_size, keep_quality)?;
                    let eof = slab.offsets.is_empty();
                    if tx.send(Ok(slab)).is_err() || eof {
                        return Ok(());
                    }
                }
            };
            if let Err(error) = run() {
                let _ = tx.send(Err(format!("{}: {error:#}", path.display())));
            }
        });
        Self {
            rx: Some(rx),
            recycle: Some(recycle),
            handle: Some(handle),
        }
    }

    fn receive(&self) -> Result<Slab> {
        self.rx
            .as_ref()
            .unwrap()
            .recv()
            .map_err(|_| anyhow!("FASTQ reader stopped"))?
            .map_err(anyhow::Error::msg)
    }
}

impl Drop for ReaderPipe {
    fn drop(&mut self) {
        // Cancel either a blocked producer send or a wait for a recycled slab before joining.
        self.rx.take();
        self.recycle.take();
        if let Some(handle) = self.handle.take() {
            let _ = handle.join();
        }
    }
}

pub(super) struct Batch {
    r1: Slab,
    r2: Slab,
    first_read: u64,
}

pub(super) struct Pipeline {
    r1: ReaderPipe,
    r2: ReaderPipe,
    next_read: u64,
}

impl Pipeline {
    pub(super) fn new(r1: PathBuf, r2: PathBuf, batch_size: usize) -> Self {
        Self {
            r1: ReaderPipe::new(r1, batch_size, false),
            r2: ReaderPipe::new(r2, batch_size, true),
            next_read: 0,
        }
    }

    pub(super) fn next_batch(&mut self) -> Result<Option<Batch>> {
        let r1 = self.r1.receive()?;
        let r2 = self.r2.receive()?;
        if r1.offsets.len() != r2.offsets.len() {
            return Err(anyhow!(
                "R1/R2 record count mismatch near read {}",
                self.next_read
            ));
        }
        if r1.offsets.is_empty() {
            return Ok(None);
        }
        let first_read = self.next_read;
        self.next_read += r1.offsets.len() as u64;
        Ok(Some(Batch { r1, r2, first_read }))
    }

    pub(super) fn recycle(&self, batch: Batch) {
        let _ = self.r1.recycle.as_ref().unwrap().send(batch.r1);
        let _ = self.r2.recycle.as_ref().unwrap().send(batch.r2);
    }
}

struct Worker {
    consensus: EcConsensusScratch,
    scratch: EcScoreScratch,
    encoded: EncodedSeq,
    resolved: EcResolvedMapping,
    output: Vec<u8>,
    counts: GeneEcDrainResult,
}

#[derive(Default)]
struct OutputChunk {
    output: Vec<u8>,
    counts: GeneEcDrainResult,
}

pub(super) struct Workers {
    states: Vec<Mutex<Worker>>,
    chunks: Vec<OutputChunk>,
}

impl Workers {
    pub(super) fn new(num_genes: usize) -> Self {
        let threads = rayon::current_num_threads().max(1);
        Self {
            states: (0..threads)
                .map(|_| {
                    Mutex::new(Worker {
                        consensus: EcConsensusScratch::default(),
                        scratch: EcScoreScratch::new(num_genes),
                        encoded: EncodedSeq::default(),
                        resolved: EcResolvedMapping::default(),
                        output: Vec::new(),
                        counts: GeneEcDrainResult::default(),
                    })
                })
                .collect(),
            chunks: (0..threads * 4).map(|_| OutputChunk::default()).collect(),
        }
    }

    pub(super) fn map(
        &mut self,
        index: &LoadedEcIndex,
        batch: &Batch,
        args: &MapArgs,
        write_output: bool,
    ) -> Result<GeneEcDrainResult> {
        let chunk_size = batch.r2.offsets.len().div_ceil(self.chunks.len());
        self.chunks
            .par_iter_mut()
            .enumerate()
            .try_for_each(|(chunk_id, chunk)| -> Result<()> {
                // There is no nested Rayon work while this state is borrowed.
                // Locks are normally uncontended, and keep caller-thread use safe.
                let state_id = rayon::current_thread_index().unwrap_or(0) % self.states.len();
                let mut worker = self.states[state_id].lock().unwrap();
                worker.counts = GeneEcDrainResult::default();
                std::mem::swap(&mut worker.output, &mut chunk.output);
                worker.output.clear();
                let start = chunk_id * chunk_size;
                let end = (start + chunk_size).min(batch.r2.offsets.len());
                if args.splicing_state {
                    for i in start..end {
                        worker.map_read::<true>(
                            index,
                            batch.first_read + i as u64,
                            batch.r1.seq(i),
                            batch.r2.seq(i),
                            batch.r2.qual(i),
                            args,
                            write_output,
                        )?;
                    }
                } else {
                    for i in start..end {
                        worker.map_read::<false>(
                            index,
                            batch.first_read + i as u64,
                            batch.r1.seq(i),
                            batch.r2.seq(i),
                            batch.r2.qual(i),
                            args,
                            write_output,
                        )?;
                    }
                }
                std::mem::swap(&mut worker.output, &mut chunk.output);
                chunk.counts = std::mem::take(&mut worker.counts);
                Ok(())
            })?;
        let mut total = GeneEcDrainResult::default();
        for chunk in &mut self.chunks {
            let counts = std::mem::take(&mut chunk.counts);
            add_ec_stats(&mut total.ec_stats, counts.ec_stats);
            total.assignment_counts.add_all(counts.assignment_counts);
            total.total_gene_candidates += counts.total_gene_candidates;
            total.total_bases += counts.total_bases;
            total.num_reads += counts.num_reads;
        }
        Ok(total)
    }

    pub(super) fn write(&self, writer: &mut GeneEcOutputWriter) -> Result<()> {
        for chunk in &self.chunks {
            writer.write_all(&chunk.output)?;
        }
        Ok(())
    }
}

impl Worker {
    #[allow(clippy::too_many_arguments)]
    fn map_read<const SPLICING: bool>(
        &mut self,
        index: &LoadedEcIndex,
        read_id: u64,
        r1: &[u8],
        r2: &[u8],
        qual: &[u8],
        args: &MapArgs,
        write_output: bool,
    ) -> Result<()> {
        let trim = if args.trim_tso {
            tso_prefix_trim_len(r2, args.min_tso_match_len, args.max_tso_mismatches)
        } else {
            0
        };
        let seq = &r2[trim..];
        let qual = &qual[trim..];
        let valid_barcode = r1.len() >= 28 && r1[..28].iter().all(|&b| acgt_code(b) < 4);
        let (ty, flags, gene_id, score, candidate_count) =
            if !valid_barcode || mean_phred_quality(qual) < args.min_mean_quality {
                self.resolved.top_gene_ids.clear();
                (
                    AssignmentType::LowQuality,
                    if valid_barcode { 2 } else { 4 },
                    None,
                    0,
                    0,
                )
            } else if encode_and_check_complexity(seq, index.k(), &mut self.encoded) {
                self.resolved.top_gene_ids.clear();
                (AssignmentType::LowComplexity, 1, None, 0, 0)
            } else {
                let stats = if SPLICING {
                    index.map_encoded_splicing_into(
                        &self.encoded,
                        args.sparse_probe_stride,
                        args.max_postings_per_seed as u32,
                        &mut self.consensus,
                        &mut self.resolved,
                    )
                } else if args.gene_ec_consensus {
                    index.map_encoded_consensus_into(
                        &self.encoded,
                        args.sparse_probe_stride,
                        args.max_postings_per_seed as u32,
                        &mut self.consensus,
                        &mut self.resolved,
                    )
                } else {
                    index.map_encoded_sparse_resolved_into(
                        &self.encoded,
                        args.search_reverse_complement,
                        args.sparse_probe_stride,
                        args.max_postings_per_seed as u32,
                        write_output && matches!(args.output_format, MapOutputFormatArg::GeneEcRad),
                        &mut self.scratch,
                        &mut self.resolved,
                    )
                };
                add_ec_stats(&mut self.counts.ec_stats, stats);
                let m = &self.resolved;
                let ty = if m.best_score == 0 {
                    AssignmentType::Unmapped
                } else if m.gene_id.is_some() {
                    AssignmentType::UniqueGene
                } else {
                    AssignmentType::AmbiguousGene
                };
                let flags = if m.best_score == 0 {
                    0
                } else {
                    (u16::from(m.has_exon_support) * SCORE_FLAG_TARGET_EXON)
                        | (u16::from(m.has_intron_support) * SCORE_FLAG_TARGET_INTRON)
                        | (u16::from(m.has_antisense_support) * SCORE_FLAG_ANTISENSE)
                };
                (ty, flags, m.gene_id, m.score, m.gene_candidate_count)
            };
        self.counts.assignment_counts.add_fields(ty, flags);
        self.counts.num_reads += 1;
        self.counts.total_bases += seq.len() as u64;
        self.counts.total_gene_candidates += candidate_count as u64;
        if !write_output {
            return Ok(());
        }
        // Preserve even short/invalid barcode strings and the historical RAD valid mask.
        let mut padded = [b'N'; 28];
        let n = r1.len().min(28);
        padded[..n].copy_from_slice(&r1[..n]);
        let cb = String::from_utf8_lossy(&padded[..16]);
        let umi = String::from_utf8_lossy(&padded[16..]);
        match args.output_format {
            MapOutputFormatArg::Tsv => {
                write!(
                    &mut self.output,
                    "{read_id}\t{cb}\t{umi}\t{}\t",
                    ty.as_str()
                )?;
                if let Some(gene_id) = gene_id {
                    write!(&mut self.output, "{gene_id}")?;
                } else {
                    self.output.push(b'.');
                }
                write!(&mut self.output, "\t.\t{candidate_count}\t{score}\t{flags}")?;
                if SPLICING {
                    self.output.push(b'\t');
                    for (i, &gene) in self.resolved.top_gene_ids.iter().enumerate() {
                        if i > 0 {
                            self.output.push(b';');
                        }
                        write!(
                            &mut self.output,
                            "{gene}:{}",
                            self.consensus.splicing_state(gene).label()
                        )?;
                    }
                }
                self.output.push(b'\n');
            }
            MapOutputFormatArg::GeneEcRad => {
                let (code, mask) = pack_barcode_umi(&cb, &umi);
                let start = self.output.len();
                write_gene_ec_rad_fields(
                    &mut self.output,
                    code,
                    mask,
                    gene_id,
                    candidate_count,
                    score,
                    flags,
                    ty,
                    &self.resolved.top_gene_ids,
                )?;
                if SPLICING {
                    self.output.truncate(start + 28);
                    for &gene in &self.resolved.top_gene_ids {
                        self.output.extend_from_slice(&gene.to_le_bytes());
                        self.output.push(self.consensus.splicing_state(gene) as u8);
                    }
                }
            }
        }
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn slabs_recycle_and_cancel_on_errors() {
        let dir = tempfile::tempdir().unwrap();
        let left = dir.path().join("left.fq");
        let right = dir.path().join("right.fq.gz");
        let records: String = (0..25).map(|i| format!("@{i}\nacgt\n+\nIIII\n")).collect();
        std::fs::write(&left, &records).unwrap();
        // Concatenated gzip members must remain valid input.
        let mut compressed = Vec::new();
        for part in [&records[..records.len() / 5], &records[records.len() / 5..]] {
            let mut encoder =
                flate2::write::GzEncoder::new(Vec::new(), flate2::Compression::fast());
            encoder.write_all(part.as_bytes()).unwrap();
            compressed.extend(encoder.finish().unwrap());
        }
        std::fs::write(&right, compressed).unwrap();
        let mut pipeline = Pipeline::new(left.clone(), right.clone(), 3);
        let mut count = 0;
        while let Some(batch) = pipeline.next_batch().unwrap() {
            assert_eq!(batch.first_read, count);
            for i in 0..batch.r1.offsets.len() {
                assert_eq!(batch.r1.seq(i), b"ACGT");
                assert_eq!(batch.r2.qual(i), b"IIII");
            }
            count += batch.r1.offsets.len() as u64;
            pipeline.recycle(batch);
        }
        assert_eq!(count, 25);
        drop(pipeline);
        // Early cancellation while producers are blocked must not deadlock.
        drop(Pipeline::new(left.clone(), right.clone(), 1));
        let broken_gzip = dir.path().join("broken.fq.gz");
        let compressed = std::fs::read(&right).unwrap();
        for truncated in [true, false] {
            let mut broken = compressed.clone();
            if truncated {
                broken.truncate(broken.len() - 3);
            } else {
                let crc_byte = broken.len() - 8;
                broken[crc_byte] ^= 1;
            }
            std::fs::write(&broken_gzip, broken).unwrap();
            let mut pipeline = Pipeline::new(left.clone(), broken_gzip.clone(), 3);
            loop {
                match pipeline.next_batch() {
                    Ok(Some(batch)) => pipeline.recycle(batch),
                    Err(_) => break,
                    Ok(None) => panic!("corrupt gzip input accepted"),
                }
            }
        }
        let bad = dir.path().join("bad.fq");
        for malformed in [
            "@x\nACGT\n+\n",
            "@x\nACGT\n+\nII\n",
            "x\nACGT\n+\nIIII\n",
            "@x\nACGT\nx\nIIII\n",
        ] {
            std::fs::write(&bad, malformed).unwrap();
            let mut pipeline = Pipeline::new(left.clone(), bad.clone(), 3);
            assert!(pipeline.next_batch().is_err());
        }
        std::fs::write(&bad, "@x\nACGT\n+\nIIII\n").unwrap();
        let mut pipeline = Pipeline::new(left, bad, 3);
        assert!(pipeline.next_batch().is_err());
    }
}

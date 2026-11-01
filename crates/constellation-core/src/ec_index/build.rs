//! Bounded sequence batches, radix partitioning, and parallel EC construction.
//! Partitions contain complete keys, so deduplication/class union is exact.
use super::*;
use rayon::prelude::*;

struct TempParts(PathBuf);
impl Drop for TempParts {
    fn drop(&mut self) {
        let _ = std::fs::remove_dir_all(&self.0);
    }
}
struct Target {
    gene: GeneId,
    class: u8,
    sequence: Vec<u8>,
}
#[derive(Default)]
struct Part {
    kmers: Vec<EcKmerRecord>,
    offsets: Vec<u32>,
    targets: Vec<EcGeneTarget>,
}

impl EcIndex {
    /// Partition by high key bits instead of performing a serial heap merge.
    /// The output order/EC IDs match the streaming builder exactly.
    pub fn build_partitioned(
        fasta: impl AsRef<Path>,
        t2g: Option<impl AsRef<Path>>,
        k: u8,
        sampling: EcSampling,
    ) -> Result<Self, EcIndexError> {
        if !(1..=32).contains(&k) {
            return Err(EcIndexError::InvalidK(k));
        }
        sampling.validate(k)?;
        let started = std::time::Instant::now();
        let annotations = match t2g {
            Some(p) => load_t2g_target_info(p)?,
            None => FxHashMap::default(),
        };
        let mut genes = Vec::new();
        let mut ids = FxHashMap::default();
        let bits = (2 * k as usize).min(8);
        let shift = 2 * k as usize - bits;
        let parts = 1usize << bits;
        let dir = std::env::temp_dir().join(format!(
            "constellation_parts_{}_{}",
            std::process::id(),
            unique_temp_suffix()
        ));
        std::fs::create_dir(&dir)?;
        let guard = TempParts(dir);
        let paths: Vec<_> = (0..parts)
            .map(|i| guard.0.join(format!("{i}.bin")))
            .collect();
        let mut writers: Vec<_> = paths
            .iter()
            .map(|p| File::create(p).map(|f| BufWriter::with_capacity(64 << 10, f)))
            .collect::<Result<_, _>>()?;
        let mut batch = Vec::new();
        let mut bases = 0;
        let mut reader = needletail::parse_fastx_file(fasta.as_ref())
            .map_err(|e| crate::index_build::BuildIndexError::Malformed(e.to_string()))?;
        while let Some(record) = reader.next() {
            let record = record
                .map_err(|e| crate::index_build::BuildIndexError::Malformed(e.to_string()))?;
            let header = String::from_utf8_lossy(record.id());
            let target = annotations
                .get(parse_target_name(&header))
                .cloned()
                .unwrap_or_else(|| TargetInfo {
                    gene: parse_gene_from_header(&header),
                    class_bits: infer_class_bits(&header),
                });
            let gene = *ids.entry(target.gene.clone()).or_insert_with(|| {
                let id = genes.len() as GeneId;
                genes.push(target.gene);
                id
            });
            let sequence = record.seq().to_vec();
            bases += sequence.len();
            batch.push(Target {
                gene,
                class: target.class_bits,
                sequence,
            });
            if bases >= 32 << 20 {
                partition_batch(&batch, k, sampling, shift, &mut writers)?;
                batch.clear();
                bases = 0;
            }
        }
        partition_batch(&batch, k, sampling, shift, &mut writers)?;
        drop(batch);
        for w in &mut writers {
            w.flush()?;
        }
        drop(writers);
        let temp_bytes: u64 = paths
            .iter()
            .map(|p| std::fs::metadata(p).map(|m| m.len()))
            .collect::<Result<Vec<_>, _>>()?
            .iter()
            .sum();
        build_phase("partition_generate", started);
        if std::env::var_os("CONSTELLATION_PROFILE_BUILD").is_some() {
            eprintln!("build_temp_bytes\t{temp_bytes}");
        }
        let started = std::time::Instant::now();
        let partitions: Vec<Part> = paths
            .par_iter()
            .map(|path| build_part(path))
            .collect::<Result<_, _>>()?;
        build_phase("partition_sort_intern", started);
        let started = std::time::Instant::now();
        let mut output = Self {
            k,
            sampling,
            genes,
            kmers: Vec::with_capacity(partitions.iter().map(|p| p.kmers.len()).sum()),
            ec_offsets: vec![0],
            ec_gene_ids: Vec::new(),
            ec_class_bits: Vec::new(),
        };
        let mut ec_ids: FxHashMap<Vec<EcGeneTarget>, u32> = FxHashMap::default();
        for mut part in partitions {
            let mut remap = Vec::with_capacity(part.offsets.len().saturating_sub(1));
            for range in part.offsets.windows(2) {
                let targets = &part.targets[range[0] as usize..range[1] as usize];
                let id = if let Some(&id) = ec_ids.get(targets) {
                    id
                } else {
                    let id = u32::try_from(ec_ids.len())
                        .map_err(|_| EcIndexError::TooLarge("EC count"))?;
                    output.ec_gene_ids.extend(targets.iter().map(|t| t.gene_id));
                    output
                        .ec_class_bits
                        .extend(targets.iter().map(|t| t.class_bits));
                    output.ec_offsets.push(
                        u32::try_from(output.ec_gene_ids.len())
                            .map_err(|_| EcIndexError::TooLarge("ec_gene_ids"))?,
                    );
                    ec_ids.insert(targets.to_vec(), id);
                    id
                };
                remap.push(id);
            }
            for record in &mut part.kmers {
                record.ec_id = remap[record.ec_id as usize];
            }
            output.kmers.append(&mut part.kmers);
        }
        build_phase("partition_combine", started);
        Ok(output)
    }
}

fn partition_batch(
    batch: &[Target],
    k: u8,
    sampling: EcSampling,
    shift: usize,
    writers: &mut [BufWriter<File>],
) -> Result<(), EcIndexError> {
    if batch.is_empty() {
        return Ok(());
    }
    let chunk = batch.len().div_ceil(rayon::current_num_threads()).max(1);
    let buckets: Vec<Vec<Vec<u128>>> = batch
        .par_chunks(chunk)
        .map(|targets| {
            let mut buckets = vec![Vec::new(); writers.len()];
            for target in targets {
                let encoded = encode_acgt(&target.sequence);
                let mut add = |seed: Kmer| {
                    buckets[(seed.code >> shift) as usize].push(
                        ((seed.code as u128) << 64)
                            | ((target.gene as u128) << 8)
                            | target.class as u128,
                    );
                };
                match sampling {
                    EcSampling::AllKmers => iter_kmers_2bit(&encoded, k).for_each(&mut add),
                    EcSampling::ClosedSyncmer { s } => {
                        iter_closed_syncmers_2bit(&encoded, k, s).for_each(&mut add)
                    }
                }
            }
            buckets
        })
        .collect();
    for buckets in buckets {
        for (writer, records) in writers.iter_mut().zip(buckets) {
            writer.write_all(bytemuck::cast_slice(&records))?;
        }
    }
    Ok(())
}

fn build_part(path: &Path) -> Result<Part, EcIndexError> {
    let mut file = File::open(path)?;
    let size = usize::try_from(file.metadata()?.len())
        .map_err(|_| EcIndexError::TooLarge("partition size"))?;
    if size % 16 != 0 {
        return Err(EcIndexError::Truncated);
    }
    let mut records = vec![0u128; size / 16];
    file.read_exact(bytemuck::cast_slice_mut(&mut records))?;
    drop(file);
    std::fs::remove_file(path)?;
    records.sort_unstable();
    let mut out = Part {
        offsets: vec![0],
        ..Part::default()
    };
    let mut ids: FxHashMap<Vec<EcGeneTarget>, u32> = FxHashMap::default();
    let mut targets: Vec<EcGeneTarget> = Vec::new();
    let mut i = 0;
    while i < records.len() {
        let code = (records[i] >> 64) as u64;
        targets.clear();
        while i < records.len() && (records[i] >> 64) as u64 == code {
            let gene = (records[i] >> 8) as u32;
            let class = records[i] as u8;
            if let Some(last) = targets.last_mut().filter(|t| t.gene_id == gene) {
                last.class_bits |= class;
            } else {
                targets.push(EcGeneTarget {
                    gene_id: gene,
                    class_bits: class,
                });
            }
            i += 1;
        }
        let id = if let Some(&id) = ids.get(&targets) {
            id
        } else {
            let id = u32::try_from(ids.len())
                .map_err(|_| EcIndexError::TooLarge("partition EC count"))?;
            out.targets.extend_from_slice(&targets);
            out.offsets.push(
                u32::try_from(out.targets.len())
                    .map_err(|_| EcIndexError::TooLarge("partition EC targets"))?,
            );
            ids.insert(targets.clone(), id);
            id
        };
        out.kmers.push(EcKmerRecord {
            kmer_code: code,
            ec_id: id,
            gene_df: targets.len() as u32,
        });
    }
    Ok(out)
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn partitioned_matches_streaming_for_classes_repeats_and_all_k_sizes() {
        let dir = tempfile::tempdir().unwrap();
        let fa = dir.path().join("ref.fa");
        let t2g = dir.path().join("t2g.tsv");
        let mut sequence = String::new();
        let mut state = 17u64;
        for _ in 0..2000 {
            state = splitmix64(state);
            sequence.push(b"ACGT"[(state & 3) as usize] as char);
        }
        std::fs::write(
            &fa,
            format!(
                ">a\n{sequence}\n>b\n{}NN{}\n>c\n{}\n>d\n{}\n",
                &sequence[..1500].to_lowercase(),
                &sequence[..800],
                &sequence[300..],
                "ACGT".repeat(800)
            ),
        )
        .unwrap();
        std::fs::write(&t2g, "a\tG\tS\nb\tG\tU\nc\tH\tS\nd\tI\t?\n").unwrap();
        for k in [1, 3, 15, 31, 32] {
            for sampling in [
                EcSampling::AllKmers,
                EcSampling::ClosedSyncmer { s: k.min(9) },
            ] {
                let expected =
                    EcIndex::build_from_fasta_with_t2g_and_sampling(&fa, Some(&t2g), k, sampling)
                        .unwrap();
                let actual = EcIndex::build_partitioned(&fa, Some(&t2g), k, sampling).unwrap();
                assert_eq!(actual, expected);
                let path = dir.path().join("direct.idx");
                actual.save_ptrhash(&path).unwrap();
                let loaded = LoadedEcIndex::load(&path).unwrap();
                for read in [
                    sequence.as_bytes(),
                    b"NNACGTACGTNN",
                    b"ACGTACGTACGTACGTACGTACGTACGTACGTACGT",
                ] {
                    assert_eq!(
                        loaded.map_read_sparse(read, false, 1, u32::MAX),
                        LoadedEcIndex::Owned(expected.clone()).map_read_sparse(
                            read,
                            false,
                            1,
                            u32::MAX
                        )
                    );
                }
            }
        }
    }
    #[test]
    fn empty_reference_has_supported_empty_index() {
        let d = tempfile::tempdir().unwrap();
        let fa = d.path().join("ref.fa");
        std::fs::write(&fa, ">empty\nNNN\n").unwrap();
        let index =
            EcIndex::build_partitioned(&fa, None::<&Path>, 31, EcSampling::ClosedSyncmer { s: 20 })
                .unwrap();
        assert!(index.kmers.is_empty());
        index.save_ptrhash(d.path().join("empty.idx")).unwrap();
        assert_eq!(
            LoadedEcIndex::load(d.path().join("empty.idx"))
                .unwrap()
                .num_genes(),
            1
        );
    }
}

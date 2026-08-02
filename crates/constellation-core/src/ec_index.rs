use crate::dna::{encode_acgt, iter_kmers_2bit, reverse_complement, EncodedSeq};
use crate::GeneId;
use bytemuck::{Pod, Zeroable};
use memmap2::{Advice, Mmap};
use rustc_hash::{FxHashMap, FxHashSet};
use std::cmp::Reverse;
use std::collections::{BTreeMap, BinaryHeap};
use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Read, Write};
use std::path::{Path, PathBuf};
use thiserror::Error;

const EC_MAGIC: &[u8; 8] = b"CSTLEC1\0";
const EC_MMAP_MAGIC: [u8; 8] = *b"CSTLEC2\0";
const EC_MMAP_VERSION_V3: u32 = 3;
const EC_MMAP_VERSION_V4: u32 = 4;
const EC_MMAP_VERSION_V5: u32 = 5;
const EC_MMAP_VERSION_V6: u32 = 6;
const EC_MMAP_VERSION: u32 = 7;
const EC_MMAP_VERSION_V8: u32 = 8;
const EC_PREFIX_BITS: u8 = 24;
const EC_V6_BLOCK_SIZE: usize = 32;
const EC_MPH_GROUP_SIZE: usize = 4;
const EC_MPH_LOAD_NUMERATOR: usize = 123;
const EC_MPH_LOAD_DENOMINATOR: usize = 100;
const EC_BUILD_BUFFER_RECORDS: usize = 8_000_000;
const EC_CLASS_UNKNOWN: u8 = 0;
const EC_CLASS_EXON: u8 = 1;
const EC_CLASS_INTRON: u8 = 2;
const EC_CLASS_ANTISENSE: u8 = 4;
const EC_V4_TARGET_WIDTH: usize = 5;
const EC_V5_TARGET_WIDTH: usize = 5;

#[derive(Debug, Clone, Copy, Default, PartialEq, Eq, Zeroable, Pod)]
#[repr(C)]
pub struct EcKmerRecord {
    pub kmer_code: u64,
    pub ec_id: u32,
    pub gene_df: u32,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct EcIndex {
    pub k: u8,
    pub genes: Vec<String>,
    pub kmers: Vec<EcKmerRecord>,
    pub ec_offsets: Vec<u32>,
    pub ec_gene_ids: Vec<GeneId>,
    pub ec_class_bits: Vec<u8>,
}

#[derive(Debug)]
pub enum LoadedEcIndex {
    Owned(EcIndex),
    Mmap(MmapEcIndex),
}

#[derive(Debug)]
pub struct MmapEcIndex {
    mmap: Mmap,
    header: EcMmapHeader,
    header_v4: Option<EcMmapHeaderV4>,
    header_v5: Option<EcMmapHeaderV5>,
    header_v6: Option<EcMmapHeaderV6>,
    header_v7: Option<EcMmapHeaderV7>,
    header_v8: Option<EcMmapHeaderV8>,
    prefix_offsets: Vec<u32>,
    seed_scores: Vec<u16>,
}

#[derive(Debug, Clone, Copy, Default, PartialEq, Eq, Zeroable, Pod)]
#[repr(C)]
struct EcMmapHeaderV1 {
    magic: [u8; 8],
    version: u32,
    k: u32,
    num_genes: u64,
    num_kmers: u64,
    num_ec_offsets: u64,
    num_ec_gene_ids: u64,
    genes_offset: u64,
    kmers_offset: u64,
    ec_offsets_offset: u64,
    ec_gene_ids_offset: u64,
    gene_names_offset: u64,
    gene_names_len: u64,
}

#[derive(Debug, Clone, Copy, Default, PartialEq, Eq, Zeroable, Pod)]
#[repr(C)]
struct EcMmapHeaderV2 {
    magic: [u8; 8],
    version: u32,
    k: u32,
    num_genes: u64,
    num_kmers: u64,
    num_ec_offsets: u64,
    num_ec_gene_ids: u64,
    genes_offset: u64,
    kmers_offset: u64,
    ec_offsets_offset: u64,
    ec_gene_ids_offset: u64,
    gene_names_offset: u64,
    gene_names_len: u64,
    ec_class_bits_offset: u64,
    ec_class_bits_len: u64,
}

#[derive(Debug, Clone, Copy, Default, PartialEq, Eq, Zeroable, Pod)]
#[repr(C)]
struct EcMmapHeader {
    magic: [u8; 8],
    version: u32,
    k: u32,
    prefix_bits: u32,
    _pad0: u32,
    num_genes: u64,
    num_kmers: u64,
    num_ec_offsets: u64,
    num_ec_gene_ids: u64,
    genes_offset: u64,
    kmers_offset: u64,
    ec_offsets_offset: u64,
    ec_gene_ids_offset: u64,
    gene_names_offset: u64,
    gene_names_len: u64,
    ec_class_bits_offset: u64,
    ec_class_bits_len: u64,
    prefix_offsets_offset: u64,
    prefix_offsets_len: u64,
}

#[derive(Debug, Clone, Copy, Default, PartialEq, Eq, Zeroable, Pod)]
#[repr(C)]
struct EcMmapHeaderV4 {
    magic: [u8; 8],
    version: u32,
    k: u32,
    prefix_bits: u32,
    suffix_bits: u32,
    suffix_width: u32,
    gene_df_width: u32,
    _pad0: u32,
    _pad1: u32,
    num_genes: u64,
    num_kmers: u64,
    num_ec_offsets: u64,
    num_ec_targets: u64,
    genes_offset: u64,
    suffixes_offset: u64,
    ec_ids_offset: u64,
    gene_dfs_offset: u64,
    ec_offsets_offset: u64,
    ec_targets_offset: u64,
    gene_names_offset: u64,
    gene_names_len: u64,
    prefix_offsets_offset: u64,
    prefix_offsets_len: u64,
}

#[derive(Debug, Clone, Copy, Default, PartialEq, Eq, Zeroable, Pod)]
#[repr(C)]
struct EcMmapHeaderV5 {
    magic: [u8; 8],
    version: u32,
    k: u32,
    prefix_bits: u32,
    suffix_bits: u32,
    suffix_width: u32,
    _pad1: u32,
    num_genes: u64,
    num_kmers: u64,
    num_ec_offsets: u64,
    num_ec_targets: u64,
    genes_offset: u64,
    suffixes_offset: u64,
    payloads_offset: u64,
    ec_offsets_offset: u64,
    ec_targets_offset: u64,
    gene_names_offset: u64,
    gene_names_len: u64,
    prefix_offsets_offset: u64,
    prefix_offsets_len: u64,
}

#[derive(Debug, Clone, Copy, Default, PartialEq, Eq, Zeroable, Pod)]
#[repr(C)]
struct EcMmapHeaderV6 {
    magic: [u8; 8],
    version: u32,
    k: u32,
    prefix_bits: u32,
    suffix_bits: u32,
    suffix_width: u32,
    block_size: u32,
    num_genes: u64,
    num_kmers: u64,
    num_blocks: u64,
    num_ec_offsets: u64,
    num_ec_targets: u64,
    genes_offset: u64,
    prefix_offsets_offset: u64,
    prefix_offsets_len: u64,
    block_offsets_offset: u64,
    block_max_suffixes_offset: u64,
    suffixes_offset: u64,
    payloads_offset: u64,
    ec_offsets_offset: u64,
    ec_targets_offset: u64,
    gene_names_offset: u64,
    gene_names_len: u64,
}

#[derive(Debug, Clone, Copy, Default, PartialEq, Eq, Zeroable, Pod)]
#[repr(C)]
struct EcMmapHeaderV7 {
    magic: [u8; 8],
    version: u32,
    k: u32,
    prefix_bits: u32,
    suffix_bits: u32,
    suffix_width: u32,
    block_size: u32,
    target_width: u32,
    _pad0: u32,
    num_genes: u64,
    num_kmers: u64,
    num_blocks: u64,
    num_ec_ranges: u64,
    num_ec_targets: u64,
    genes_offset: u64,
    prefix_offsets_offset: u64,
    prefix_offsets_len: u64,
    block_offsets_offset: u64,
    block_max_suffixes_offset: u64,
    suffixes_offset: u64,
    payloads_offset: u64,
    ec_ranges_offset: u64,
    ec_targets_offset: u64,
    gene_names_offset: u64,
    gene_names_len: u64,
}

#[derive(Debug, Clone, Copy, Default, PartialEq, Eq, Zeroable, Pod)]
#[repr(C)]
struct EcMmapHeaderV8 {
    magic: [u8; 8],
    version: u32,
    k: u32,
    prefix_bits: u32,
    suffix_bits: u32,
    suffix_width: u32,
    group_size: u32,
    num_genes: u64,
    num_slots: u64,
    num_groups: u64,
    num_ec_offsets: u64,
    num_ec_targets: u64,
    genes_offset: u64,
    slot_offsets_offset: u64,
    slot_offsets_len: u64,
    group_offsets_offset: u64,
    group_offsets_len: u64,
    group_displacements_offset: u64,
    suffixes_offset: u64,
    payloads_offset: u64,
    ec_offsets_offset: u64,
    ec_targets_offset: u64,
    gene_names_offset: u64,
    gene_names_len: u64,
}

#[derive(Debug, Clone, Copy, Default, PartialEq, Eq, Zeroable, Pod)]
#[repr(C)]
struct EcV5KmerPayload {
    target_start: u32,
    gene_df: u16,
    target_len: u16,
}

#[derive(Debug, Clone, Copy, Default, PartialEq, Eq, Zeroable, Pod)]
#[repr(C)]
struct EcTargetRange {
    target_start: u32,
    target_len: u16,
    _pad0: u16,
}

#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub struct EcLookupHit {
    kmer_code: u64,
    gene_df: u32,
    target_start: u32,
    target_len: u32,
    ec_id: Option<u32>,
}

impl From<EcMmapHeaderV1> for EcMmapHeader {
    fn from(header: EcMmapHeaderV1) -> Self {
        Self {
            magic: header.magic,
            version: header.version,
            k: header.k,
            prefix_bits: 0,
            _pad0: 0,
            num_genes: header.num_genes,
            num_kmers: header.num_kmers,
            num_ec_offsets: header.num_ec_offsets,
            num_ec_gene_ids: header.num_ec_gene_ids,
            genes_offset: header.genes_offset,
            kmers_offset: header.kmers_offset,
            ec_offsets_offset: header.ec_offsets_offset,
            ec_gene_ids_offset: header.ec_gene_ids_offset,
            gene_names_offset: header.gene_names_offset,
            gene_names_len: header.gene_names_len,
            ec_class_bits_offset: 0,
            ec_class_bits_len: 0,
            prefix_offsets_offset: 0,
            prefix_offsets_len: 0,
        }
    }
}

impl From<EcMmapHeaderV2> for EcMmapHeader {
    fn from(header: EcMmapHeaderV2) -> Self {
        Self {
            magic: header.magic,
            version: header.version,
            k: header.k,
            prefix_bits: 0,
            _pad0: 0,
            num_genes: header.num_genes,
            num_kmers: header.num_kmers,
            num_ec_offsets: header.num_ec_offsets,
            num_ec_gene_ids: header.num_ec_gene_ids,
            genes_offset: header.genes_offset,
            kmers_offset: header.kmers_offset,
            ec_offsets_offset: header.ec_offsets_offset,
            ec_gene_ids_offset: header.ec_gene_ids_offset,
            gene_names_offset: header.gene_names_offset,
            gene_names_len: header.gene_names_len,
            ec_class_bits_offset: header.ec_class_bits_offset,
            ec_class_bits_len: header.ec_class_bits_len,
            prefix_offsets_offset: 0,
            prefix_offsets_len: 0,
        }
    }
}

impl From<EcMmapHeaderV4> for EcMmapHeader {
    fn from(header: EcMmapHeaderV4) -> Self {
        Self {
            magic: header.magic,
            version: header.version,
            k: header.k,
            prefix_bits: header.prefix_bits,
            _pad0: 0,
            num_genes: header.num_genes,
            num_kmers: header.num_kmers,
            num_ec_offsets: header.num_ec_offsets,
            num_ec_gene_ids: header.num_ec_targets,
            genes_offset: header.genes_offset,
            kmers_offset: 0,
            ec_offsets_offset: header.ec_offsets_offset,
            ec_gene_ids_offset: header.ec_targets_offset,
            gene_names_offset: header.gene_names_offset,
            gene_names_len: header.gene_names_len,
            ec_class_bits_offset: 0,
            ec_class_bits_len: 0,
            prefix_offsets_offset: header.prefix_offsets_offset,
            prefix_offsets_len: header.prefix_offsets_len,
        }
    }
}

impl From<EcMmapHeaderV5> for EcMmapHeader {
    fn from(header: EcMmapHeaderV5) -> Self {
        Self {
            magic: header.magic,
            version: header.version,
            k: header.k,
            prefix_bits: header.prefix_bits,
            _pad0: 0,
            num_genes: header.num_genes,
            num_kmers: header.num_kmers,
            num_ec_offsets: header.num_ec_offsets,
            num_ec_gene_ids: header.num_ec_targets,
            genes_offset: header.genes_offset,
            kmers_offset: header.suffixes_offset,
            ec_offsets_offset: header.ec_offsets_offset,
            ec_gene_ids_offset: header.ec_targets_offset,
            gene_names_offset: header.gene_names_offset,
            gene_names_len: header.gene_names_len,
            ec_class_bits_offset: 0,
            ec_class_bits_len: 0,
            prefix_offsets_offset: header.prefix_offsets_offset,
            prefix_offsets_len: header.prefix_offsets_len,
        }
    }
}

impl From<EcMmapHeaderV6> for EcMmapHeader {
    fn from(header: EcMmapHeaderV6) -> Self {
        Self {
            magic: header.magic,
            version: header.version,
            k: header.k,
            prefix_bits: header.prefix_bits,
            _pad0: 0,
            num_genes: header.num_genes,
            num_kmers: header.num_kmers,
            num_ec_offsets: header.num_ec_offsets,
            num_ec_gene_ids: header.num_ec_targets,
            genes_offset: header.genes_offset,
            kmers_offset: 0,
            ec_offsets_offset: header.ec_offsets_offset,
            ec_gene_ids_offset: header.ec_targets_offset,
            gene_names_offset: header.gene_names_offset,
            gene_names_len: header.gene_names_len,
            ec_class_bits_offset: 0,
            ec_class_bits_len: 0,
            prefix_offsets_offset: header.prefix_offsets_offset,
            prefix_offsets_len: header.prefix_offsets_len,
        }
    }
}

impl From<EcMmapHeaderV7> for EcMmapHeader {
    fn from(header: EcMmapHeaderV7) -> Self {
        Self {
            magic: header.magic,
            version: header.version,
            k: header.k,
            prefix_bits: header.prefix_bits,
            _pad0: 0,
            num_genes: header.num_genes,
            num_kmers: header.num_kmers,
            num_ec_offsets: header.num_ec_ranges,
            num_ec_gene_ids: header.num_ec_targets,
            genes_offset: header.genes_offset,
            kmers_offset: 0,
            ec_offsets_offset: header.ec_ranges_offset,
            ec_gene_ids_offset: header.ec_targets_offset,
            gene_names_offset: header.gene_names_offset,
            gene_names_len: header.gene_names_len,
            ec_class_bits_offset: 0,
            ec_class_bits_len: 0,
            prefix_offsets_offset: header.prefix_offsets_offset,
            prefix_offsets_len: header.prefix_offsets_len,
        }
    }
}

impl From<EcMmapHeaderV8> for EcMmapHeader {
    fn from(header: EcMmapHeaderV8) -> Self {
        Self {
            magic: header.magic,
            version: header.version,
            k: header.k,
            prefix_bits: header.prefix_bits,
            _pad0: 0,
            num_genes: header.num_genes,
            num_kmers: header.num_slots,
            num_ec_offsets: header.num_ec_offsets,
            num_ec_gene_ids: header.num_ec_targets,
            genes_offset: header.genes_offset,
            kmers_offset: 0,
            ec_offsets_offset: header.ec_offsets_offset,
            ec_gene_ids_offset: header.ec_targets_offset,
            gene_names_offset: header.gene_names_offset,
            gene_names_len: header.gene_names_len,
            ec_class_bits_offset: 0,
            ec_class_bits_len: 0,
            prefix_offsets_offset: header.slot_offsets_offset,
            prefix_offsets_len: header.slot_offsets_len,
        }
    }
}

#[derive(Debug, Clone, Copy, Default, PartialEq, Eq, Zeroable, Pod)]
#[repr(C)]
struct EcGeneMeta {
    name_start: u32,
    name_len: u32,
}

#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub struct EcMappingStats {
    pub seed_lookups: u64,
    pub distinct_seed_lookups: u64,
    pub seeds_found: u64,
    pub seeds_skipped_high_gene_df: u64,
    pub reads_with_gene_support: u64,
    pub ec_hits: u64,
    pub gene_votes: u64,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct EcReadMapping {
    pub gene_scores: Vec<EcGeneScore>,
    pub best_score: u16,
    pub ec_hits: u32,
}

#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub struct EcGeneScore {
    pub gene_id: GeneId,
    pub total_score: u16,
    pub exon_score: u16,
    pub intron_score: u16,
    pub antisense_score: u16,
    pub unknown_score: u16,
}

/// Reusable per-worker gene-score storage for the read-at-a-time mapper.
///
/// Epoch tags make resetting proportional to the genes touched by a read
/// rather than to the number of genes in the reference. The score and epoch
/// arrays are intentionally dense: for a human reference they fit in roughly
/// one MiB and avoid a hash-table allocation for every read.
#[derive(Debug, Default)]
pub struct EcScoreScratch {
    scores: Vec<EcGeneScore>,
    epochs: Vec<u32>,
    touched_gene_ids: Vec<GeneId>,
    epoch: u32,
}

impl EcScoreScratch {
    pub fn new(num_genes: usize) -> Self {
        Self {
            scores: vec![EcGeneScore::default(); num_genes],
            epochs: vec![0; num_genes],
            touched_gene_ids: Vec::with_capacity(32),
            epoch: 0,
        }
    }

    fn begin_read(&mut self, num_genes: usize) {
        if self.scores.len() != num_genes {
            self.scores.resize(num_genes, EcGeneScore::default());
            self.epochs.resize(num_genes, 0);
            self.epochs.fill(0);
            self.epoch = 0;
        }
        self.touched_gene_ids.clear();
        if self.epoch == u32::MAX {
            self.epochs.fill(0);
            self.epoch = 1;
        } else {
            self.epoch += 1;
        }
    }

    #[inline]
    fn add(&mut self, gene_id: GeneId, seed_score: u16, class_bits: u8) {
        let idx = gene_id as usize;
        if self.epochs[idx] != self.epoch {
            self.epochs[idx] = self.epoch;
            self.scores[idx] = EcGeneScore {
                gene_id,
                ..EcGeneScore::default()
            };
            self.touched_gene_ids.push(gene_id);
        }
        add_to_gene_score(&mut self.scores[idx], seed_score, class_bits);
    }

    fn finish(&self, stats: EcMappingStats) -> (EcReadMapping, EcMappingStats) {
        finish_ec_scores(
            self.touched_gene_ids
                .iter()
                .map(|&gene_id| self.scores[gene_id as usize])
                .collect(),
            stats,
        )
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
struct KmerGeneRecord {
    kmer_code: u64,
    gene_id: GeneId,
    class_bits: u8,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
struct EcGeneTarget {
    gene_id: GeneId,
    class_bits: u8,
}

#[derive(Debug, Clone, PartialEq, Eq)]
struct TargetInfo {
    gene: String,
    class_bits: u8,
}

#[derive(Debug, Error)]
pub enum EcIndexError {
    #[error("EC index I/O error: {0}")]
    Io(#[from] std::io::Error),
    #[error("EC index build error: {0}")]
    Build(#[from] crate::index_build::BuildIndexError),
    #[error("invalid EC index magic")]
    InvalidMagic,
    #[error("truncated EC index")]
    Truncated,
    #[error("misaligned EC index section")]
    Misaligned,
    #[error("EC index version {0} is unsupported")]
    UnsupportedVersion(u32),
    #[error("k must be between 1 and 32, got {0}")]
    InvalidK(u8),
    #[error("EC index field is too large: {0}")]
    TooLarge(&'static str),
}

impl EcIndex {
    pub fn build_from_fasta(fasta: impl AsRef<Path>, k: u8) -> Result<Self, EcIndexError> {
        Self::build_from_fasta_with_t2g(fasta, None::<&Path>, k)
    }

    pub fn build_from_fasta_with_t2g(
        fasta: impl AsRef<Path>,
        t2g_map: Option<impl AsRef<Path>>,
        k: u8,
    ) -> Result<Self, EcIndexError> {
        if !(1..=32).contains(&k) {
            return Err(EcIndexError::InvalidK(k));
        }
        let target_info = match t2g_map {
            Some(path) => load_t2g_target_info(path)?,
            None => FxHashMap::default(),
        };
        let mut gene_ids: FxHashMap<String, GeneId> = FxHashMap::default();
        let mut genes = Vec::new();
        let mut buffer = Vec::with_capacity(EC_BUILD_BUFFER_RECORDS);
        let mut runs = Vec::new();
        let run_prefix = std::env::temp_dir().join(format!(
            "constellation_ec_build_{}_{}",
            std::process::id(),
            unique_temp_suffix()
        ));

        let mut reader = needletail::parse_fastx_file(fasta.as_ref())
            .map_err(|err| crate::index_build::BuildIndexError::Malformed(err.to_string()))?;
        while let Some(record) = reader.next() {
            let record = record
                .map_err(|err| crate::index_build::BuildIndexError::Malformed(err.to_string()))?;
            let header = String::from_utf8_lossy(record.id()).into_owned();
            let target_name = parse_target_name(&header);
            let inferred = target_info
                .get(target_name)
                .cloned()
                .unwrap_or_else(|| TargetInfo {
                    gene: parse_gene_from_header(&header),
                    class_bits: infer_class_bits(&header),
                });
            let gene_id = match gene_ids.get(&inferred.gene) {
                Some(&id) => id,
                None => {
                    let id = genes.len() as GeneId;
                    gene_ids.insert(inferred.gene.clone(), id);
                    genes.push(inferred.gene.clone());
                    id
                }
            };
            let seq = record
                .seq()
                .iter()
                .map(|base| base.to_ascii_uppercase())
                .collect::<Vec<_>>();
            let encoded = encode_acgt(&seq);
            let mut seen_in_gene_target = FxHashSet::default();
            for kmer in iter_kmers_2bit(&encoded, k) {
                if seen_in_gene_target.insert(kmer.code) {
                    buffer.push(KmerGeneRecord {
                        kmer_code: kmer.code,
                        gene_id,
                        class_bits: inferred.class_bits,
                    });
                    if buffer.len() >= EC_BUILD_BUFFER_RECORDS {
                        runs.push(flush_kmer_gene_run(&mut buffer, &run_prefix, runs.len())?);
                    }
                }
            }
        }
        if !buffer.is_empty() {
            runs.push(flush_kmer_gene_run(&mut buffer, &run_prefix, runs.len())?);
        }

        let mut ec_ids: BTreeMap<Vec<EcGeneTarget>, u32> = BTreeMap::new();
        let mut ec_offsets = vec![0_u32];
        let mut ec_gene_ids = Vec::new();
        let mut ec_class_bits = Vec::new();
        let mut kmers = Vec::new();
        merge_kmer_gene_runs(&runs, |kmer_code, gene_targets| {
            let ec_id = match ec_ids.get(&gene_targets) {
                Some(&id) => id,
                None => {
                    let id = ec_ids.len() as u32;
                    ec_gene_ids.extend(gene_targets.iter().map(|target| target.gene_id));
                    ec_class_bits.extend(gene_targets.iter().map(|target| target.class_bits));
                    ec_offsets.push(
                        u32::try_from(ec_gene_ids.len())
                            .map_err(|_| EcIndexError::TooLarge("ec_gene_ids"))?,
                    );
                    ec_ids.insert(gene_targets.clone(), id);
                    id
                }
            };
            kmers.push(EcKmerRecord {
                kmer_code,
                ec_id,
                gene_df: gene_targets.len() as u32,
            });
            Ok(())
        })?;
        for run in runs {
            let _ = std::fs::remove_file(run);
        }

        Ok(Self {
            k,
            genes,
            kmers,
            ec_offsets,
            ec_gene_ids,
            ec_class_bits,
        })
    }

    pub fn load(path: impl AsRef<Path>) -> Result<Self, EcIndexError> {
        let mut reader = BufReader::new(File::open(path)?);
        let mut magic = [0_u8; 8];
        reader.read_exact(&mut magic)?;
        if &magic != EC_MAGIC {
            return Err(EcIndexError::InvalidMagic);
        }
        let k = read_u32(&mut reader)? as u8;
        let num_genes = read_u64(&mut reader)? as usize;
        let num_kmers = read_u64(&mut reader)? as usize;
        let num_ec_offsets = read_u64(&mut reader)? as usize;
        let num_ec_gene_ids = read_u64(&mut reader)? as usize;
        let gene_names_len = read_u64(&mut reader)? as usize;

        let mut gene_name_bytes = vec![0_u8; gene_names_len];
        reader.read_exact(&mut gene_name_bytes)?;
        let genes = split_nul_strings(&gene_name_bytes, num_genes)?;

        let mut ec_offsets = Vec::with_capacity(num_ec_offsets);
        for _ in 0..num_ec_offsets {
            ec_offsets.push(read_u32(&mut reader)?);
        }
        let mut ec_gene_ids = Vec::with_capacity(num_ec_gene_ids);
        for _ in 0..num_ec_gene_ids {
            ec_gene_ids.push(read_u32(&mut reader)?);
        }
        let ec_class_bits = vec![EC_CLASS_UNKNOWN; ec_gene_ids.len()];
        let mut kmers = Vec::with_capacity(num_kmers);
        for _ in 0..num_kmers {
            kmers.push(EcKmerRecord {
                kmer_code: read_u64(&mut reader)?,
                ec_id: read_u32(&mut reader)?,
                gene_df: read_u32(&mut reader)?,
            });
        }

        Ok(Self {
            k,
            genes,
            kmers,
            ec_offsets,
            ec_gene_ids,
            ec_class_bits,
        })
    }

    pub fn load_auto(path: impl AsRef<Path>) -> Result<LoadedEcIndex, EcIndexError> {
        LoadedEcIndex::load(path)
    }

    pub fn save(&self, path: impl AsRef<Path>) -> Result<(), EcIndexError> {
        let mut writer = BufWriter::new(File::create(path)?);
        writer.write_all(EC_MAGIC)?;
        write_u32(&mut writer, self.k as u32)?;
        write_u64(&mut writer, self.genes.len() as u64)?;
        write_u64(&mut writer, self.kmers.len() as u64)?;
        write_u64(&mut writer, self.ec_offsets.len() as u64)?;
        write_u64(&mut writer, self.ec_gene_ids.len() as u64)?;
        let gene_names = join_nul_strings(&self.genes);
        write_u64(&mut writer, gene_names.len() as u64)?;
        writer.write_all(&gene_names)?;
        for &offset in &self.ec_offsets {
            write_u32(&mut writer, offset)?;
        }
        for &gene_id in &self.ec_gene_ids {
            write_u32(&mut writer, gene_id)?;
        }
        for record in &self.kmers {
            write_u64(&mut writer, record.kmer_code)?;
            write_u32(&mut writer, record.ec_id)?;
            write_u32(&mut writer, record.gene_df)?;
        }
        Ok(())
    }

    pub fn save_mmap(&self, path: impl AsRef<Path>) -> Result<(), EcIndexError> {
        let prefix_bits = EC_PREFIX_BITS.min(self.k.saturating_mul(2));
        let prefix_offsets = build_prefix_offsets(&self.kmers, self.k, prefix_bits)?;
        let mut gene_names = Vec::new();
        let mut gene_metas = Vec::with_capacity(self.genes.len());
        for gene in &self.genes {
            let name_start = u32::try_from(gene_names.len())
                .map_err(|_| EcIndexError::TooLarge("gene_names"))?;
            let name_len =
                u32::try_from(gene.len()).map_err(|_| EcIndexError::TooLarge("gene name"))?;
            gene_names.extend_from_slice(gene.as_bytes());
            gene_metas.push(EcGeneMeta {
                name_start,
                name_len,
            });
        }

        save_mmap_locality_from_slices(
            path,
            self.k,
            &gene_metas,
            &gene_names,
            &self.kmers,
            &self.ec_offsets,
            &self.ec_gene_ids,
            &self.ec_class_bits,
            &prefix_offsets,
        )
    }

    #[cfg(test)]
    fn save_mmap_v3(&self, path: impl AsRef<Path>) -> Result<(), EcIndexError> {
        let prefix_bits = EC_PREFIX_BITS.min(self.k.saturating_mul(2));
        let prefix_offsets = build_prefix_offsets(&self.kmers, self.k, prefix_bits)?;
        let mut gene_names = Vec::new();
        let mut gene_metas = Vec::with_capacity(self.genes.len());
        for gene in &self.genes {
            let name_start = u32::try_from(gene_names.len())
                .map_err(|_| EcIndexError::TooLarge("gene_names"))?;
            let name_len =
                u32::try_from(gene.len()).map_err(|_| EcIndexError::TooLarge("gene name"))?;
            gene_names.extend_from_slice(gene.as_bytes());
            gene_metas.push(EcGeneMeta {
                name_start,
                name_len,
            });
        }

        let mut cursor = align_up(std::mem::size_of::<EcMmapHeader>() as u64, 8);
        let prefix_offsets_offset = cursor;
        cursor += byte_len::<u32>(prefix_offsets.len());
        cursor = align_up(cursor, 8);
        let genes_offset = cursor;
        cursor += byte_len::<EcGeneMeta>(gene_metas.len());
        cursor = align_up(cursor, 8);
        let kmers_offset = cursor;
        cursor += byte_len::<EcKmerRecord>(self.kmers.len());
        cursor = align_up(cursor, 4);
        let ec_offsets_offset = cursor;
        cursor += byte_len::<u32>(self.ec_offsets.len());
        cursor = align_up(cursor, 4);
        let ec_gene_ids_offset = cursor;
        cursor += byte_len::<u32>(self.ec_gene_ids.len());
        cursor = align_up(cursor, 8);
        let ec_class_bits_offset = cursor;
        cursor += self.ec_class_bits.len() as u64;
        cursor = align_up(cursor, 8);
        let gene_names_offset = cursor;

        let header = EcMmapHeader {
            magic: EC_MMAP_MAGIC,
            version: EC_MMAP_VERSION_V3,
            k: self.k as u32,
            prefix_bits: prefix_bits as u32,
            _pad0: 0,
            num_genes: self.genes.len() as u64,
            num_kmers: self.kmers.len() as u64,
            num_ec_offsets: self.ec_offsets.len() as u64,
            num_ec_gene_ids: self.ec_gene_ids.len() as u64,
            genes_offset,
            kmers_offset,
            ec_offsets_offset,
            ec_gene_ids_offset,
            gene_names_offset,
            gene_names_len: gene_names.len() as u64,
            ec_class_bits_offset,
            ec_class_bits_len: self.ec_class_bits.len() as u64,
            prefix_offsets_offset,
            prefix_offsets_len: byte_len::<u32>(prefix_offsets.len()),
        };

        let mut file = File::create(path)?;
        file.write_all(bytemuck::bytes_of(&header))?;
        write_padding(
            &mut file,
            std::mem::size_of::<EcMmapHeader>() as u64,
            prefix_offsets_offset,
        )?;
        write_pod_slice(&mut file, &prefix_offsets)?;
        write_padding(
            &mut file,
            prefix_offsets_offset + byte_len::<u32>(prefix_offsets.len()),
            genes_offset,
        )?;
        write_pod_slice(&mut file, &gene_metas)?;
        write_padding(
            &mut file,
            genes_offset + byte_len::<EcGeneMeta>(gene_metas.len()),
            kmers_offset,
        )?;
        write_pod_slice(&mut file, &self.kmers)?;
        write_padding(
            &mut file,
            kmers_offset + byte_len::<EcKmerRecord>(self.kmers.len()),
            ec_offsets_offset,
        )?;
        write_pod_slice(&mut file, &self.ec_offsets)?;
        write_padding(
            &mut file,
            ec_offsets_offset + byte_len::<u32>(self.ec_offsets.len()),
            ec_gene_ids_offset,
        )?;
        write_pod_slice(&mut file, &self.ec_gene_ids)?;
        write_padding(
            &mut file,
            ec_gene_ids_offset + byte_len::<u32>(self.ec_gene_ids.len()),
            ec_class_bits_offset,
        )?;
        file.write_all(&self.ec_class_bits)?;
        write_padding(
            &mut file,
            ec_class_bits_offset + self.ec_class_bits.len() as u64,
            gene_names_offset,
        )?;
        file.write_all(&gene_names)?;
        Ok(())
    }

    pub fn gene_name(&self, gene_id: GeneId) -> Option<&str> {
        self.genes.get(gene_id as usize).map(String::as_str)
    }

    pub fn lookup(&self, kmer_code: u64) -> Option<EcKmerRecord> {
        self.kmers
            .binary_search_by_key(&kmer_code, |record| record.kmer_code)
            .ok()
            .map(|idx| self.kmers[idx])
    }

    pub fn ec_genes(&self, ec_id: u32) -> &[GeneId] {
        let start = self.ec_offsets.get(ec_id as usize).copied().unwrap_or(0) as usize;
        let end = self
            .ec_offsets
            .get(ec_id as usize + 1)
            .copied()
            .unwrap_or(start as u32) as usize;
        &self.ec_gene_ids[start..end]
    }

    pub fn ec_class_bits(&self, ec_id: u32) -> &[u8] {
        let start = self.ec_offsets.get(ec_id as usize).copied().unwrap_or(0) as usize;
        let end = self
            .ec_offsets
            .get(ec_id as usize + 1)
            .copied()
            .unwrap_or(start as u32) as usize;
        &self.ec_class_bits[start..end]
    }

    pub fn map_read(
        &self,
        seq: &[u8],
        search_reverse_complement: bool,
    ) -> (EcReadMapping, EcMappingStats) {
        self.map_read_sparse(seq, search_reverse_complement, 1, u32::MAX)
    }

    pub fn map_read_sparse(
        &self,
        seq: &[u8],
        search_reverse_complement: bool,
        stride: u32,
        max_gene_df: u32,
    ) -> (EcReadMapping, EcMappingStats) {
        let mut stats = EcMappingStats::default();
        let mut gene_scores: FxHashMap<GeneId, EcGeneScore> = FxHashMap::default();
        self.accumulate_seq(
            seq,
            stride.max(1),
            max_gene_df.max(1),
            &mut gene_scores,
            &mut stats,
        );
        if search_reverse_complement {
            let rc = reverse_complement(&encode_acgt(seq));
            let rc = crate::dna::decode_to_vec(&rc);
            self.accumulate_seq(
                &rc,
                stride.max(1),
                max_gene_df.max(1),
                &mut gene_scores,
                &mut stats,
            );
        }
        finish_ec_mapping(gene_scores, stats)
    }

    fn accumulate_seq(
        &self,
        seq: &[u8],
        stride: u32,
        max_gene_df: u32,
        gene_scores: &mut FxHashMap<GeneId, EcGeneScore>,
        stats: &mut EcMappingStats,
    ) {
        let encoded = encode_acgt(seq);
        for kmer in iter_kmers_2bit(&encoded, self.k) {
            if kmer.pos % stride != 0 {
                continue;
            }
            stats.seed_lookups += 1;
            let Some(record) = self.lookup(kmer.code) else {
                continue;
            };
            if record.gene_df > max_gene_df {
                stats.seeds_skipped_high_gene_df += 1;
                continue;
            }
            stats.seeds_found += 1;
            stats.ec_hits += 1;
            let seed_score = ec_seed_score(self.genes.len(), record.gene_df);
            let genes = self.ec_genes(record.ec_id);
            let class_bits = self.ec_class_bits(record.ec_id);
            for (idx, &gene_id) in genes.iter().enumerate() {
                let class_bits = class_bits.get(idx).copied().unwrap_or(EC_CLASS_UNKNOWN);
                add_gene_score(gene_scores, gene_id, seed_score, class_bits);
                stats.gene_votes += 1;
            }
        }
    }
}

impl LoadedEcIndex {
    pub fn load(path: impl AsRef<Path>) -> Result<Self, EcIndexError> {
        let path = path.as_ref();
        let file = File::open(path)?;
        let mmap = unsafe { Mmap::map(&file)? };
        if mmap.starts_with(&EC_MMAP_MAGIC) {
            Ok(Self::Mmap(MmapEcIndex::from_mmap(mmap)?))
        } else if mmap.starts_with(EC_MAGIC) {
            Ok(Self::Owned(EcIndex::load(path)?))
        } else {
            Err(EcIndexError::InvalidMagic)
        }
    }

    pub fn save_mmap(&self, path: impl AsRef<Path>) -> Result<(), EcIndexError> {
        match self {
            Self::Owned(index) => index.save_mmap(path),
            Self::Mmap(index) => index.save_mmap(path),
        }
    }

    pub fn preload(&self) -> std::io::Result<()> {
        match self {
            Self::Owned(_) => Ok(()),
            Self::Mmap(index) => index.preload(),
        }
    }

    pub fn k(&self) -> u8 {
        match self {
            Self::Owned(index) => index.k,
            Self::Mmap(index) => index.k(),
        }
    }

    pub fn kmer_prefix(&self, kmer_code: u64) -> u32 {
        match self {
            Self::Owned(index) => ec_kmer_prefix(kmer_code, index.k, EC_PREFIX_BITS),
            Self::Mmap(index) => index.kmer_prefix(kmer_code),
        }
    }

    pub fn num_genes(&self) -> usize {
        match self {
            Self::Owned(index) => index.genes.len(),
            Self::Mmap(index) => index.num_genes(),
        }
    }

    pub fn gene_name(&self, gene_id: GeneId) -> Option<&str> {
        match self {
            Self::Owned(index) => index.gene_name(gene_id),
            Self::Mmap(index) => index.gene_name(gene_id),
        }
    }

    pub fn lookup(&self, kmer_code: u64) -> Option<EcKmerRecord> {
        match self {
            Self::Owned(index) => index.lookup(kmer_code),
            Self::Mmap(index) => index.lookup(kmer_code),
        }
    }

    pub fn lookup_hit(&self, kmer_code: u64) -> Option<EcLookupHit> {
        match self {
            Self::Owned(index) => index.lookup(kmer_code).map(|record| {
                let offsets = &index.ec_offsets;
                let start = offsets.get(record.ec_id as usize).copied().unwrap_or(0);
                let end = offsets
                    .get(record.ec_id as usize + 1)
                    .copied()
                    .unwrap_or(start);
                EcLookupHit {
                    kmer_code,
                    gene_df: record.gene_df,
                    target_start: start,
                    target_len: end.saturating_sub(start),
                    ec_id: Some(record.ec_id),
                }
            }),
            Self::Mmap(index) => index.lookup_hit(kmer_code),
        }
    }

    pub fn seed_score(&self, gene_df: u32) -> u16 {
        match self {
            Self::Owned(index) => ec_seed_score(index.genes.len(), gene_df),
            Self::Mmap(index) => index.seed_score(gene_df),
        }
    }

    pub fn ec_genes(&self, ec_id: u32) -> &[GeneId] {
        match self {
            Self::Owned(index) => index.ec_genes(ec_id),
            Self::Mmap(index) => index.ec_genes(ec_id),
        }
    }

    pub fn ec_class_bits(&self, ec_id: u32) -> &[u8] {
        match self {
            Self::Owned(index) => index.ec_class_bits(ec_id),
            Self::Mmap(index) => index.ec_class_bits(ec_id),
        }
    }

    pub fn for_each_ec_target(&self, ec_id: u32, mut f: impl FnMut(GeneId, u8)) {
        match self {
            Self::Owned(index) => {
                let genes = index.ec_genes(ec_id);
                let class_bits = index.ec_class_bits(ec_id);
                for (idx, &gene_id) in genes.iter().enumerate() {
                    f(
                        gene_id,
                        class_bits.get(idx).copied().unwrap_or(EC_CLASS_UNKNOWN),
                    );
                }
            }
            Self::Mmap(index) => index.for_each_ec_target(ec_id, f),
        }
    }

    pub fn accumulate_lookup_record(
        &self,
        record: EcKmerRecord,
        max_gene_df: u32,
        gene_scores: &mut FxHashMap<GeneId, EcGeneScore>,
        stats: &mut EcMappingStats,
    ) {
        self.visit_lookup_record(
            record,
            max_gene_df,
            stats,
            |gene_id, seed_score, class_bits| {
                add_gene_score(gene_scores, gene_id, seed_score, class_bits);
            },
        );
    }

    #[inline]
    fn visit_lookup_record(
        &self,
        record: EcKmerRecord,
        max_gene_df: u32,
        stats: &mut EcMappingStats,
        mut visit: impl FnMut(GeneId, u16, u8),
    ) {
        if record.gene_df > max_gene_df {
            stats.seeds_skipped_high_gene_df += 1;
            return;
        }
        stats.seeds_found += 1;
        stats.ec_hits += 1;
        let seed_score = self.seed_score(record.gene_df);
        let mut gene_votes = 0_u64;
        self.for_each_ec_target(record.ec_id, |gene_id, class_bits| {
            visit(gene_id, seed_score, class_bits);
            gene_votes += 1;
        });
        stats.gene_votes += gene_votes;
    }

    pub fn accumulate_lookup_hit(
        &self,
        hit: EcLookupHit,
        max_gene_df: u32,
        gene_scores: &mut FxHashMap<GeneId, EcGeneScore>,
        stats: &mut EcMappingStats,
    ) {
        self.visit_lookup_hit(
            hit,
            max_gene_df,
            stats,
            |gene_id, seed_score, class_bits| {
                add_gene_score(gene_scores, gene_id, seed_score, class_bits);
            },
        );
    }

    #[inline]
    fn visit_lookup_hit(
        &self,
        hit: EcLookupHit,
        max_gene_df: u32,
        stats: &mut EcMappingStats,
        mut visit: impl FnMut(GeneId, u16, u8),
    ) {
        if hit.gene_df > max_gene_df {
            stats.seeds_skipped_high_gene_df += 1;
            return;
        }
        stats.seeds_found += 1;
        stats.ec_hits += 1;
        let seed_score = self.seed_score(hit.gene_df);
        let mut gene_votes = 0_u64;
        match self {
            Self::Owned(index) => {
                let start = hit.target_start as usize;
                let end = start + hit.target_len as usize;
                let class_bits = &index.ec_class_bits;
                for (idx, &gene_id) in index.ec_gene_ids[start..end].iter().enumerate() {
                    visit(
                        gene_id,
                        seed_score,
                        class_bits
                            .get(start + idx)
                            .copied()
                            .unwrap_or(EC_CLASS_UNKNOWN),
                    );
                    gene_votes += 1;
                }
            }
            Self::Mmap(index) => {
                if index.header_v5.is_some()
                    || index.header_v6.is_some()
                    || index.header_v7.is_some()
                    || index.header_v8.is_some()
                {
                    index.for_each_ec_target_range(
                        hit.target_start,
                        hit.target_len,
                        |gene_id, class_bits| {
                            visit(gene_id, seed_score, class_bits);
                            gene_votes += 1;
                        },
                    );
                } else if let Some(ec_id) = hit.ec_id {
                    index.for_each_ec_target(ec_id, |gene_id, class_bits| {
                        visit(gene_id, seed_score, class_bits);
                        gene_votes += 1;
                    });
                }
            }
        }
        stats.gene_votes += gene_votes;
    }

    pub fn map_read_sparse(
        &self,
        seq: &[u8],
        search_reverse_complement: bool,
        stride: u32,
        max_gene_df: u32,
    ) -> (EcReadMapping, EcMappingStats) {
        let encoded = encode_acgt(seq);
        self.map_encoded_sparse(&encoded, search_reverse_complement, stride, max_gene_df)
    }

    pub fn map_encoded_sparse(
        &self,
        encoded: &EncodedSeq,
        search_reverse_complement: bool,
        stride: u32,
        max_gene_df: u32,
    ) -> (EcReadMapping, EcMappingStats) {
        let mut stats = EcMappingStats::default();
        let mut gene_scores: FxHashMap<GeneId, EcGeneScore> = FxHashMap::default();
        accumulate_ec_seq(
            self,
            encoded,
            stride.max(1),
            max_gene_df.max(1),
            &mut gene_scores,
            &mut stats,
        );
        if search_reverse_complement {
            accumulate_ec_seq(
                self,
                &reverse_complement(encoded),
                stride.max(1),
                max_gene_df.max(1),
                &mut gene_scores,
                &mut stats,
            );
        }
        finish_ec_mapping(gene_scores, stats)
    }

    /// Maps an encoded read while reusing dense per-worker score storage.
    ///
    /// The lookup and scoring semantics are identical to [`Self::map_encoded_sparse`];
    /// only the temporary accumulator differs.
    pub fn map_encoded_sparse_with_scratch(
        &self,
        encoded: &EncodedSeq,
        search_reverse_complement: bool,
        stride: u32,
        max_gene_df: u32,
        scratch: &mut EcScoreScratch,
    ) -> (EcReadMapping, EcMappingStats) {
        scratch.begin_read(self.num_genes());
        let mut stats = EcMappingStats::default();
        accumulate_ec_seq_dense(
            self,
            encoded,
            stride.max(1),
            max_gene_df.max(1),
            scratch,
            &mut stats,
        );
        if search_reverse_complement {
            accumulate_ec_seq_dense(
                self,
                &reverse_complement(encoded),
                stride.max(1),
                max_gene_df.max(1),
                scratch,
                &mut stats,
            );
        }
        scratch.finish(stats)
    }
}

impl MmapEcIndex {
    fn from_mmap(mmap: Mmap) -> Result<Self, EcIndexError> {
        if mmap.len() < std::mem::size_of::<EcMmapHeaderV1>() {
            return Err(EcIndexError::Truncated);
        }
        let header_v1: EcMmapHeaderV1 =
            *bytemuck::from_bytes(&mmap[..std::mem::size_of::<EcMmapHeaderV1>()]);
        if header_v1.magic != EC_MMAP_MAGIC {
            return Err(EcIndexError::InvalidMagic);
        }
        let mut header_v4 = None;
        let mut header_v5 = None;
        let mut header_v6 = None;
        let mut header_v7 = None;
        let mut header_v8 = None;
        let header = if header_v1.version == 1 {
            EcMmapHeader::from(header_v1)
        } else if header_v1.version == 2 {
            if mmap.len() < std::mem::size_of::<EcMmapHeaderV2>() {
                return Err(EcIndexError::Truncated);
            }
            let header_v2: EcMmapHeaderV2 =
                *bytemuck::from_bytes(&mmap[..std::mem::size_of::<EcMmapHeaderV2>()]);
            EcMmapHeader::from(header_v2)
        } else if header_v1.version == EC_MMAP_VERSION_V3 {
            if mmap.len() < std::mem::size_of::<EcMmapHeader>() {
                return Err(EcIndexError::Truncated);
            }
            *bytemuck::from_bytes(&mmap[..std::mem::size_of::<EcMmapHeader>()])
        } else if header_v1.version == EC_MMAP_VERSION_V4 {
            if mmap.len() < std::mem::size_of::<EcMmapHeaderV4>() {
                return Err(EcIndexError::Truncated);
            }
            let v4: EcMmapHeaderV4 =
                *bytemuck::from_bytes(&mmap[..std::mem::size_of::<EcMmapHeaderV4>()]);
            header_v4 = Some(v4);
            EcMmapHeader::from(v4)
        } else if header_v1.version == EC_MMAP_VERSION_V5 {
            if mmap.len() < std::mem::size_of::<EcMmapHeaderV5>() {
                return Err(EcIndexError::Truncated);
            }
            let v5: EcMmapHeaderV5 =
                *bytemuck::from_bytes(&mmap[..std::mem::size_of::<EcMmapHeaderV5>()]);
            header_v5 = Some(v5);
            EcMmapHeader::from(v5)
        } else if header_v1.version == EC_MMAP_VERSION_V6 {
            if mmap.len() < std::mem::size_of::<EcMmapHeaderV6>() {
                return Err(EcIndexError::Truncated);
            }
            let v6: EcMmapHeaderV6 =
                *bytemuck::from_bytes(&mmap[..std::mem::size_of::<EcMmapHeaderV6>()]);
            header_v6 = Some(v6);
            EcMmapHeader::from(v6)
        } else if header_v1.version == EC_MMAP_VERSION {
            if mmap.len() < std::mem::size_of::<EcMmapHeaderV7>() {
                return Err(EcIndexError::Truncated);
            }
            let v7: EcMmapHeaderV7 =
                *bytemuck::from_bytes(&mmap[..std::mem::size_of::<EcMmapHeaderV7>()]);
            header_v7 = Some(v7);
            EcMmapHeader::from(v7)
        } else if header_v1.version == EC_MMAP_VERSION_V8 {
            if mmap.len() < std::mem::size_of::<EcMmapHeaderV8>() {
                return Err(EcIndexError::Truncated);
            }
            let v8: EcMmapHeaderV8 =
                *bytemuck::from_bytes(&mmap[..std::mem::size_of::<EcMmapHeaderV8>()]);
            header_v8 = Some(v8);
            EcMmapHeader::from(v8)
        } else {
            return Err(EcIndexError::UnsupportedVersion(header_v1.version));
        };
        let mut this = Self {
            mmap,
            header,
            header_v4,
            header_v5,
            header_v6,
            header_v7,
            header_v8,
            prefix_offsets: Vec::new(),
            seed_scores: Vec::new(),
        };
        this.validate()?;
        if this.header.prefix_bits > 0 {
            this.prefix_offsets = this.prefix_offsets_slice().to_vec();
        }
        this.seed_scores = build_seed_score_table(this.num_genes());
        Ok(this)
    }

    fn save_mmap(&self, path: impl AsRef<Path>) -> Result<(), EcIndexError> {
        if self.header_v6.is_some() {
            let mut file = File::create(path)?;
            file.write_all(&self.mmap)?;
            return Ok(());
        }
        if self.header_v7.is_some() {
            let mut file = File::create(path)?;
            file.write_all(&self.mmap)?;
            return Ok(());
        }
        if self.header_v8.is_some() {
            let mut file = File::create(path)?;
            file.write_all(&self.mmap)?;
            return Ok(());
        }
        if self.header_v5.is_some() {
            if ec_target_dedup_enabled() {
                return save_mmap_v7_from_mmap_v5(path, self);
            }
            if ec_mmap_v6_enabled() {
                return save_mmap_v6_from_mmap_v5(path, self);
            }
            let mut file = File::create(path)?;
            file.write_all(&self.mmap)?;
            return Ok(());
        }
        if self.header_v4.is_some() {
            if ec_mmap_v5_enabled() {
                return save_mmap_v5_from_mmap_v4(path, self);
            }
            let mut file = File::create(path)?;
            file.write_all(&self.mmap)?;
            return Ok(());
        }
        let prefix_bits = EC_PREFIX_BITS.min(self.k().saturating_mul(2));
        let prefix_offsets = build_prefix_offsets(self.kmers(), self.k(), prefix_bits)?;
        let class_bits = if self.header.version >= 2 {
            self.class_bits()
        } else {
            &[]
        };

        save_mmap_locality_from_slices(
            path,
            self.k(),
            self.genes(),
            self.gene_names(),
            self.kmers(),
            self.ec_offsets(),
            self.ec_gene_ids(),
            class_bits,
            &prefix_offsets,
        )
    }

    fn validate(&self) -> Result<(), EcIndexError> {
        self.slice::<EcGeneMeta>(self.header.genes_offset, self.header.num_genes)?;
        if self.header_v7.is_some() {
            self.slice::<EcTargetRange>(self.header.ec_offsets_offset, self.header.num_ec_offsets)?;
        } else if self.header_v8.is_some() {
            self.slice::<u32>(self.header.ec_offsets_offset, self.header.num_ec_offsets)?;
        } else {
            self.slice::<u32>(self.header.ec_offsets_offset, self.header.num_ec_offsets)?;
        }
        if let Some(v8) = self.header_v8 {
            let suffix_width = usize::try_from(v8.suffix_width)
                .map_err(|_| EcIndexError::TooLarge("EC v8 suffix width"))?;
            if v8.suffix_bits > 48 || suffix_width > 8 || v8.group_size == 0 {
                return Err(EcIndexError::Truncated);
            }
            self.bytes(v8.slot_offsets_offset, v8.slot_offsets_len)?;
            self.bytes(v8.group_offsets_offset, v8.group_offsets_len)?;
            self.slice::<u32>(v8.group_displacements_offset, v8.num_groups)?;
            self.bytes(v8.suffixes_offset, v8.num_slots * suffix_width as u64)?;
            self.slice::<EcV5KmerPayload>(v8.payloads_offset, v8.num_slots)?;
            self.bytes(
                v8.ec_targets_offset,
                v8.num_ec_targets * EC_V5_TARGET_WIDTH as u64,
            )?;
        } else if let Some(v7) = self.header_v7 {
            let suffix_width = usize::try_from(v7.suffix_width)
                .map_err(|_| EcIndexError::TooLarge("EC v7 suffix width"))?;
            if v7.suffix_bits > 48
                || suffix_width > 8
                || v7.block_size == 0
                || !matches!(v7.target_width, 4 | 5)
            {
                return Err(EcIndexError::Truncated);
            }
            self.slice::<u32>(v7.block_offsets_offset, v7.num_blocks + 1)?;
            self.bytes(
                v7.block_max_suffixes_offset,
                v7.num_blocks * suffix_width as u64,
            )?;
            self.bytes(v7.suffixes_offset, v7.num_kmers * suffix_width as u64)?;
            self.slice::<EcV5KmerPayload>(v7.payloads_offset, v7.num_kmers)?;
            self.slice::<EcTargetRange>(v7.ec_ranges_offset, v7.num_ec_ranges)?;
            self.bytes(
                v7.ec_targets_offset,
                v7.num_ec_targets * v7.target_width as u64,
            )?;
        } else if let Some(v6) = self.header_v6 {
            let suffix_width = usize::try_from(v6.suffix_width)
                .map_err(|_| EcIndexError::TooLarge("EC v6 suffix width"))?;
            if v6.suffix_bits > 48 || suffix_width > 8 || v6.block_size == 0 {
                return Err(EcIndexError::Truncated);
            }
            self.slice::<u32>(v6.block_offsets_offset, v6.num_blocks + 1)?;
            self.bytes(
                v6.block_max_suffixes_offset,
                v6.num_blocks * suffix_width as u64,
            )?;
            self.bytes(v6.suffixes_offset, v6.num_kmers * suffix_width as u64)?;
            self.slice::<EcV5KmerPayload>(v6.payloads_offset, v6.num_kmers)?;
            self.bytes(
                v6.ec_targets_offset,
                v6.num_ec_targets * EC_V5_TARGET_WIDTH as u64,
            )?;
        } else if let Some(v5) = self.header_v5 {
            let suffix_width = usize::try_from(v5.suffix_width)
                .map_err(|_| EcIndexError::TooLarge("EC v5 suffix width"))?;
            if v5.suffix_bits > 48 || suffix_width > 8 {
                return Err(EcIndexError::Truncated);
            }
            self.bytes(v5.suffixes_offset, v5.num_kmers * suffix_width as u64)?;
            self.slice::<EcV5KmerPayload>(v5.payloads_offset, v5.num_kmers)?;
            self.bytes(
                v5.ec_targets_offset,
                v5.num_ec_targets * EC_V5_TARGET_WIDTH as u64,
            )?;
        } else if let Some(v4) = self.header_v4 {
            let suffix_width = usize::try_from(v4.suffix_width)
                .map_err(|_| EcIndexError::TooLarge("EC v4 suffix width"))?;
            let gene_df_width = usize::try_from(v4.gene_df_width)
                .map_err(|_| EcIndexError::TooLarge("EC v4 gene df width"))?;
            if suffix_width > 8 || !matches!(gene_df_width, 2 | 4) {
                return Err(EcIndexError::Truncated);
            }
            self.bytes(v4.suffixes_offset, v4.num_kmers * suffix_width as u64)?;
            self.slice::<u32>(v4.ec_ids_offset, v4.num_kmers)?;
            self.bytes(v4.gene_dfs_offset, v4.num_kmers * gene_df_width as u64)?;
            self.bytes(
                v4.ec_targets_offset,
                v4.num_ec_targets * EC_V4_TARGET_WIDTH as u64,
            )?;
        } else {
            self.slice::<EcKmerRecord>(self.header.kmers_offset, self.header.num_kmers)?;
            self.slice::<u32>(self.header.ec_gene_ids_offset, self.header.num_ec_gene_ids)?;
        }
        if self.header.version >= 2
            && self.header_v4.is_none()
            && self.header_v5.is_none()
            && self.header_v6.is_none()
            && self.header_v7.is_none()
            && self.header_v8.is_none()
        {
            self.bytes(
                self.header.ec_class_bits_offset,
                self.header.ec_class_bits_len,
            )?;
        }
        if self.header.prefix_bits > 0 {
            let expected_len = byte_len::<u32>((1_usize << self.header.prefix_bits as usize) + 1);
            if self.header.prefix_offsets_len != expected_len {
                return Err(EcIndexError::Truncated);
            }
            self.bytes(
                self.header.prefix_offsets_offset,
                self.header.prefix_offsets_len,
            )?;
        }
        self.bytes(self.header.gene_names_offset, self.header.gene_names_len)?;
        Ok(())
    }

    fn preload(&self) -> std::io::Result<()> {
        self.mmap.advise(Advice::WillNeed)?;
        let mut checksum = 0_u8;
        for offset in (0..self.mmap.len()).step_by(4096) {
            checksum ^= std::hint::black_box(self.mmap[offset]);
        }
        if let Some(&last) = self.mmap.last() {
            checksum ^= std::hint::black_box(last);
        }
        std::hint::black_box(checksum);
        Ok(())
    }

    fn k(&self) -> u8 {
        self.header.k as u8
    }

    fn num_genes(&self) -> usize {
        self.header.num_genes as usize
    }

    fn gene_name(&self, gene_id: GeneId) -> Option<&str> {
        let meta = self.genes().get(gene_id as usize)?;
        let start = meta.name_start as usize;
        let end = start.checked_add(meta.name_len as usize)?;
        std::str::from_utf8(self.gene_names().get(start..end)?).ok()
    }

    fn lookup(&self, kmer_code: u64) -> Option<EcKmerRecord> {
        if self.header_v4.is_some() {
            return self.lookup_v4(kmer_code);
        }
        if self.header_v8.is_some() {
            let hit = self.lookup_v8(kmer_code)?;
            return Some(EcKmerRecord {
                kmer_code,
                ec_id: self.v5_ec_id_for_target_start(hit.target_start)?,
                gene_df: hit.gene_df,
            });
        }
        if self.header_v7.is_some() {
            let hit = self.lookup_v7(kmer_code)?;
            return Some(EcKmerRecord {
                kmer_code,
                ec_id: self.v7_ec_id_for_range(hit.target_start, hit.target_len)?,
                gene_df: hit.gene_df,
            });
        }
        if self.header_v6.is_some() {
            let hit = self.lookup_v6(kmer_code)?;
            return Some(EcKmerRecord {
                kmer_code,
                ec_id: self.v5_ec_id_for_target_start(hit.target_start)?,
                gene_df: hit.gene_df,
            });
        }
        if self.header_v5.is_some() {
            let hit = self.lookup_v5(kmer_code)?;
            return Some(EcKmerRecord {
                kmer_code,
                ec_id: self.v5_ec_id_for_target_start(hit.target_start)?,
                gene_df: hit.gene_df,
            });
        }
        let kmers = self.kmer_lookup_range(kmer_code);
        kmers
            .binary_search_by_key(&kmer_code, |record| record.kmer_code)
            .ok()
            .map(|idx| kmers[idx])
    }

    #[inline]
    fn lookup_hit(&self, kmer_code: u64) -> Option<EcLookupHit> {
        if self.header_v8.is_some() {
            return self.lookup_v8(kmer_code);
        }
        if self.header_v7.is_some() {
            return self.lookup_v7(kmer_code);
        }
        if self.header_v6.is_some() {
            return self.lookup_v6(kmer_code);
        }
        if self.header_v5.is_some() {
            return self.lookup_v5(kmer_code);
        }
        if self.header_v4.is_some() {
            let record = self.lookup_v4(kmer_code)?;
            let offsets = self.ec_offsets();
            let start = offsets.get(record.ec_id as usize).copied().unwrap_or(0);
            let end = offsets
                .get(record.ec_id as usize + 1)
                .copied()
                .unwrap_or(start);
            return Some(EcLookupHit {
                kmer_code,
                gene_df: record.gene_df,
                target_start: start,
                target_len: end.saturating_sub(start),
                ec_id: Some(record.ec_id),
            });
        }
        self.lookup(kmer_code).map(|record| {
            let offsets = self.ec_offsets();
            let start = offsets.get(record.ec_id as usize).copied().unwrap_or(0);
            let end = offsets
                .get(record.ec_id as usize + 1)
                .copied()
                .unwrap_or(start);
            EcLookupHit {
                kmer_code,
                gene_df: record.gene_df,
                target_start: start,
                target_len: end.saturating_sub(start),
                ec_id: Some(record.ec_id),
            }
        })
    }

    #[inline]
    fn lookup_v8(&self, kmer_code: u64) -> Option<EcLookupHit> {
        let v8 = self.header_v8?;
        let prefix = ec_kmer_prefix(kmer_code, self.k(), v8.prefix_bits as u8) as usize;
        let slot_offsets = self.prefix_offsets();
        let slot_start = slot_offsets.get(prefix).copied().unwrap_or(0) as usize;
        let slot_end = slot_offsets
            .get(prefix + 1)
            .copied()
            .unwrap_or(slot_start as u32) as usize;
        let slot_len = slot_end.checked_sub(slot_start)?;
        if slot_len == 0 {
            return None;
        }
        let group_offsets = self.v8_group_offsets();
        let group_start = group_offsets.get(prefix).copied().unwrap_or(0) as usize;
        let group_end = group_offsets
            .get(prefix + 1)
            .copied()
            .unwrap_or(group_start as u32) as usize;
        let group_count = group_end.checked_sub(group_start)?;
        if group_count == 0 {
            return None;
        }
        let suffix = ec_kmer_suffix(kmer_code, self.k(), v8.prefix_bits as u8);
        let group = group_start + (mph_hash1(suffix) as usize % group_count);
        let displacement = self.v8_group_displacements().get(group).copied()? as u64;
        let lane = mph_slot(suffix, displacement, slot_len);
        let idx = slot_start + lane;
        let suffix_width = v8.suffix_width as usize;
        let suffixes = self.v8_suffixes();
        if read_packed_le(&suffixes[idx * suffix_width..idx * suffix_width + suffix_width])
            != suffix
        {
            return None;
        }
        let hit = self.v5_payloads()[idx];
        if hit.target_len == 0 {
            return None;
        }
        Some(EcLookupHit {
            kmer_code,
            gene_df: hit.gene_df as u32,
            target_start: hit.target_start,
            target_len: hit.target_len as u32,
            ec_id: None,
        })
    }

    #[inline]
    fn lookup_v7(&self, kmer_code: u64) -> Option<EcLookupHit> {
        let v7 = self.header_v7?;
        let prefix = ec_kmer_prefix(kmer_code, self.k(), v7.prefix_bits as u8) as usize;
        let prefix_offsets = self.prefix_offsets();
        let block_start = prefix_offsets.get(prefix).copied().unwrap_or(0) as usize;
        let block_end = prefix_offsets
            .get(prefix + 1)
            .copied()
            .unwrap_or(block_start as u32) as usize;
        if block_start == block_end {
            return None;
        }

        let suffix = ec_kmer_suffix(kmer_code, self.k(), v7.prefix_bits as u8);
        let suffix_width = v7.suffix_width as usize;
        let block_max_suffixes = self.v7_block_max_suffixes();
        let mut lo = block_start;
        let mut hi = block_end;
        while lo < hi {
            let mid = lo + (hi - lo) / 2;
            let max_suffix = read_packed_le(
                &block_max_suffixes[mid * suffix_width..mid * suffix_width + suffix_width],
            );
            if max_suffix < suffix {
                lo = mid + 1;
            } else {
                hi = mid;
            }
        }
        if lo >= block_end {
            return None;
        }

        let block_offsets = self.v7_block_offsets();
        let record_start = block_offsets.get(lo).copied().unwrap_or(0) as usize;
        let record_end = block_offsets
            .get(lo + 1)
            .copied()
            .unwrap_or(record_start as u32) as usize;
        let lane = self.v7_find_suffix_in_block(record_start, record_end, suffix)?;
        let hit = self.v5_payloads()[record_start + lane];
        Some(EcLookupHit {
            kmer_code,
            gene_df: hit.gene_df as u32,
            target_start: hit.target_start,
            target_len: hit.target_len as u32,
            ec_id: None,
        })
    }

    #[inline]
    fn lookup_v6(&self, kmer_code: u64) -> Option<EcLookupHit> {
        let v6 = self.header_v6?;
        let prefix = ec_kmer_prefix(kmer_code, self.k(), v6.prefix_bits as u8) as usize;
        let prefix_offsets = self.prefix_offsets();
        let block_start = prefix_offsets.get(prefix).copied().unwrap_or(0) as usize;
        let block_end = prefix_offsets
            .get(prefix + 1)
            .copied()
            .unwrap_or(block_start as u32) as usize;
        if block_start == block_end {
            return None;
        }

        let suffix = ec_kmer_suffix(kmer_code, self.k(), v6.prefix_bits as u8);
        let suffix_width = v6.suffix_width as usize;
        let block_max_suffixes = self.v6_block_max_suffixes();
        let mut lo = block_start;
        let mut hi = block_end;
        while lo < hi {
            let mid = lo + (hi - lo) / 2;
            let max_suffix = read_packed_le(
                &block_max_suffixes[mid * suffix_width..mid * suffix_width + suffix_width],
            );
            if max_suffix < suffix {
                lo = mid + 1;
            } else {
                hi = mid;
            }
        }
        if lo >= block_end {
            return None;
        }

        let block_offsets = self.v6_block_offsets();
        let record_start = block_offsets.get(lo).copied().unwrap_or(0) as usize;
        let record_end = block_offsets
            .get(lo + 1)
            .copied()
            .unwrap_or(record_start as u32) as usize;
        let lane = self.v6_find_suffix_in_block(record_start, record_end, suffix)?;
        let hit = self.v5_payloads()[record_start + lane];
        Some(EcLookupHit {
            kmer_code,
            gene_df: hit.gene_df as u32,
            target_start: hit.target_start,
            target_len: hit.target_len as u32,
            ec_id: None,
        })
    }

    #[inline]
    fn lookup_v5(&self, kmer_code: u64) -> Option<EcLookupHit> {
        let v5 = self.header_v5?;
        let prefix = ec_kmer_prefix(kmer_code, self.k(), v5.prefix_bits as u8) as usize;
        let offsets = self.prefix_offsets();
        let start = offsets.get(prefix).copied().unwrap_or(0) as usize;
        let end = offsets.get(prefix + 1).copied().unwrap_or(start as u32) as usize;
        let suffix = ec_kmer_suffix(kmer_code, self.k(), v5.prefix_bits as u8);
        let suffix_width = v5.suffix_width as usize;
        let suffixes = self.v5_suffixes();
        let mut lo = start;
        let mut hi = end;
        while lo < hi {
            let mid = lo + (hi - lo) / 2;
            let mid_suffix =
                read_packed_le(&suffixes[mid * suffix_width..mid * suffix_width + suffix_width]);
            match mid_suffix.cmp(&suffix) {
                std::cmp::Ordering::Less => lo = mid + 1,
                std::cmp::Ordering::Greater => hi = mid,
                std::cmp::Ordering::Equal => {
                    lo = mid;
                    break;
                }
            }
        }
        if lo >= end
            || read_packed_le(&suffixes[lo * suffix_width..lo * suffix_width + suffix_width])
                != suffix
        {
            return None;
        }
        let hit = self.v5_payloads()[lo];
        Some(EcLookupHit {
            kmer_code,
            gene_df: hit.gene_df as u32,
            target_start: hit.target_start,
            target_len: hit.target_len as u32,
            ec_id: None,
        })
    }

    fn lookup_v4(&self, kmer_code: u64) -> Option<EcKmerRecord> {
        let v4 = self.header_v4?;
        let prefix = ec_kmer_prefix(kmer_code, self.k(), v4.prefix_bits as u8) as usize;
        let offsets = self.prefix_offsets();
        let start = offsets.get(prefix).copied().unwrap_or(0) as usize;
        let end = offsets.get(prefix + 1).copied().unwrap_or(start as u32) as usize;
        let suffix = ec_kmer_suffix(kmer_code, self.k(), v4.prefix_bits as u8);
        let suffix_width = v4.suffix_width as usize;
        let suffixes = self.v4_suffixes();
        let mut lo = start;
        let mut hi = end;
        while lo < hi {
            let mid = lo + (hi - lo) / 2;
            let mid_suffix =
                read_packed_le(&suffixes[mid * suffix_width..mid * suffix_width + suffix_width]);
            match mid_suffix.cmp(&suffix) {
                std::cmp::Ordering::Less => lo = mid + 1,
                std::cmp::Ordering::Greater => hi = mid,
                std::cmp::Ordering::Equal => {
                    lo = mid;
                    break;
                }
            }
        }
        if lo >= end
            || read_packed_le(&suffixes[lo * suffix_width..lo * suffix_width + suffix_width])
                != suffix
        {
            return None;
        }
        let idx = lo;
        Some(EcKmerRecord {
            kmer_code,
            ec_id: self.v4_ec_ids()[idx],
            gene_df: self.v4_gene_df(idx),
        })
    }

    fn seed_score(&self, gene_df: u32) -> u16 {
        self.seed_scores
            .get(gene_df as usize)
            .copied()
            .unwrap_or_else(|| ec_seed_score(self.num_genes(), gene_df))
    }

    fn kmer_prefix(&self, kmer_code: u64) -> u32 {
        let prefix_bits = if self.header.prefix_bits > 0 {
            self.header.prefix_bits as u8
        } else {
            EC_PREFIX_BITS.min(self.k().saturating_mul(2))
        };
        ec_kmer_prefix(kmer_code, self.k(), prefix_bits)
    }

    fn kmer_lookup_range(&self, kmer_code: u64) -> &[EcKmerRecord] {
        if self.header.prefix_bits == 0 {
            return self.kmers();
        }
        let prefix = ec_kmer_prefix(kmer_code, self.k(), self.header.prefix_bits as u8) as usize;
        let offsets = self.prefix_offsets();
        let start = offsets.get(prefix).copied().unwrap_or(0) as usize;
        let end = offsets.get(prefix + 1).copied().unwrap_or(start as u32) as usize;
        &self.kmers()[start..end]
    }

    fn ec_genes(&self, ec_id: u32) -> &[GeneId] {
        assert!(
            self.header_v4.is_none()
                && self.header_v5.is_none()
                && self.header_v6.is_none()
                && self.header_v7.is_none()
                && self.header_v8.is_none(),
            "packed EC gene targets are only available through iteration"
        );
        let offsets = self.ec_offsets();
        let start = offsets.get(ec_id as usize).copied().unwrap_or(0) as usize;
        let end = offsets
            .get(ec_id as usize + 1)
            .copied()
            .unwrap_or(start as u32) as usize;
        &self.ec_gene_ids()[start..end]
    }

    fn ec_class_bits(&self, ec_id: u32) -> &[u8] {
        assert!(
            self.header_v4.is_none()
                && self.header_v5.is_none()
                && self.header_v6.is_none()
                && self.header_v7.is_none()
                && self.header_v8.is_none(),
            "packed EC class targets are only available through iteration"
        );
        let offsets = self.ec_offsets();
        let start = offsets.get(ec_id as usize).copied().unwrap_or(0) as usize;
        let end = offsets
            .get(ec_id as usize + 1)
            .copied()
            .unwrap_or(start as u32) as usize;
        if self.header.version >= 2 {
            &self.class_bits()[start..end]
        } else {
            &[]
        }
    }

    fn for_each_ec_target(&self, ec_id: u32, mut f: impl FnMut(GeneId, u8)) {
        if self.header_v4.is_some() {
            self.for_each_v4_ec_target(ec_id, f);
            return;
        }
        if self.header_v7.is_some() {
            let Some(range) = self.v7_ec_ranges().get(ec_id as usize).copied() else {
                return;
            };
            self.for_each_ec_target_range(range.target_start, range.target_len as u32, f);
            return;
        }
        assert!(
            self.header_v5.is_none()
                && self.header_v6.is_none()
                && self.header_v7.is_none()
                && self.header_v8.is_none(),
            "packed EC targets are addressed by range, not ec_id"
        );
        let genes = self.ec_genes(ec_id);
        let class_bits = self.ec_class_bits(ec_id);
        for (idx, &gene_id) in genes.iter().enumerate() {
            f(
                gene_id,
                class_bits.get(idx).copied().unwrap_or(EC_CLASS_UNKNOWN),
            );
        }
    }

    fn genes(&self) -> &[EcGeneMeta] {
        self.slice(self.header.genes_offset, self.header.num_genes)
            .expect("validated EC gene metadata")
    }

    fn kmers(&self) -> &[EcKmerRecord] {
        assert!(
            self.header_v4.is_none()
                && self.header_v5.is_none()
                && self.header_v6.is_none()
                && self.header_v7.is_none()
                && self.header_v8.is_none(),
            "packed EC kmer records are stored separately"
        );
        self.slice(self.header.kmers_offset, self.header.num_kmers)
            .expect("validated EC kmers")
    }

    fn ec_offsets(&self) -> &[u32] {
        assert!(
            self.header_v7.is_none(),
            "deduped EC target ranges are stored separately"
        );
        self.slice(self.header.ec_offsets_offset, self.header.num_ec_offsets)
            .expect("validated EC offsets")
    }

    fn v5_ec_id_for_target_start(&self, target_start: u32) -> Option<u32> {
        let offsets = self.ec_offsets();
        offsets
            .binary_search(&target_start)
            .ok()
            .and_then(|idx| u32::try_from(idx).ok())
    }

    fn ec_gene_ids(&self) -> &[GeneId] {
        assert!(
            self.header_v4.is_none()
                && self.header_v5.is_none()
                && self.header_v6.is_none()
                && self.header_v7.is_none()
                && self.header_v8.is_none(),
            "packed EC gene targets are stored as packed target records"
        );
        self.slice(self.header.ec_gene_ids_offset, self.header.num_ec_gene_ids)
            .expect("validated EC gene ids")
    }

    fn gene_names(&self) -> &[u8] {
        self.bytes(self.header.gene_names_offset, self.header.gene_names_len)
            .expect("validated EC gene names")
    }

    fn class_bits(&self) -> &[u8] {
        assert!(
            self.header_v4.is_none()
                && self.header_v5.is_none()
                && self.header_v6.is_none()
                && self.header_v7.is_none()
                && self.header_v8.is_none(),
            "packed EC class bits are stored as packed target records"
        );
        self.bytes(
            self.header.ec_class_bits_offset,
            self.header.ec_class_bits_len,
        )
        .expect("validated EC class bits")
    }

    fn prefix_offsets(&self) -> &[u32] {
        if !self.prefix_offsets.is_empty() {
            return &self.prefix_offsets;
        }
        self.prefix_offsets_slice()
    }

    fn prefix_offsets_slice(&self) -> &[u32] {
        self.slice(
            self.header.prefix_offsets_offset,
            self.header.prefix_offsets_len / std::mem::size_of::<u32>() as u64,
        )
        .expect("validated EC prefix offsets")
    }

    fn slice<T: Pod>(&self, offset: u64, count: u64) -> Result<&[T], EcIndexError> {
        if offset % std::mem::align_of::<T>() as u64 != 0 {
            return Err(EcIndexError::Misaligned);
        }
        let len = byte_len::<T>(count as usize);
        let bytes = self.bytes(offset, len)?;
        bytemuck::try_cast_slice(bytes).map_err(|_| EcIndexError::Misaligned)
    }

    fn bytes(&self, offset: u64, len: u64) -> Result<&[u8], EcIndexError> {
        let start = offset as usize;
        let end = start
            .checked_add(len as usize)
            .ok_or(EcIndexError::Truncated)?;
        self.mmap.get(start..end).ok_or(EcIndexError::Truncated)
    }

    fn v5_suffixes(&self) -> &[u8] {
        let v5 = self.header_v5.expect("v5 header");
        self.bytes(v5.suffixes_offset, v5.num_kmers * v5.suffix_width as u64)
            .expect("validated EC v5 suffixes")
    }

    fn v5_payloads(&self) -> &[EcV5KmerPayload] {
        if let Some(v8) = self.header_v8 {
            return self
                .slice(v8.payloads_offset, v8.num_slots)
                .expect("validated EC v8 payloads");
        }
        if let Some(v7) = self.header_v7 {
            return self
                .slice(v7.payloads_offset, v7.num_kmers)
                .expect("validated EC v7 payloads");
        }
        if let Some(v6) = self.header_v6 {
            return self
                .slice(v6.payloads_offset, v6.num_kmers)
                .expect("validated EC v6 payloads");
        }
        let v5 = self.header_v5.expect("v5/v6 header");
        self.slice(v5.payloads_offset, v5.num_kmers)
            .expect("validated EC v5/v6 payloads")
    }

    fn v8_group_offsets(&self) -> &[u32] {
        let v8 = self.header_v8.expect("v8 header");
        self.slice(
            v8.group_offsets_offset,
            v8.group_offsets_len / std::mem::size_of::<u32>() as u64,
        )
        .expect("validated EC v8 group offsets")
    }

    fn v8_group_displacements(&self) -> &[u32] {
        let v8 = self.header_v8.expect("v8 header");
        self.slice(v8.group_displacements_offset, v8.num_groups)
            .expect("validated EC v8 group displacements")
    }

    fn v8_suffixes(&self) -> &[u8] {
        let v8 = self.header_v8.expect("v8 header");
        self.bytes(v8.suffixes_offset, v8.num_slots * v8.suffix_width as u64)
            .expect("validated EC v8 suffixes")
    }

    fn v7_block_offsets(&self) -> &[u32] {
        let v7 = self.header_v7.expect("v7 header");
        self.slice(v7.block_offsets_offset, v7.num_blocks + 1)
            .expect("validated EC v7 block offsets")
    }

    fn v7_block_max_suffixes(&self) -> &[u8] {
        let v7 = self.header_v7.expect("v7 header");
        self.bytes(
            v7.block_max_suffixes_offset,
            v7.num_blocks * v7.suffix_width as u64,
        )
        .expect("validated EC v7 block max suffixes")
    }

    fn v7_suffixes(&self) -> &[u8] {
        let v7 = self.header_v7.expect("v7 header");
        self.bytes(v7.suffixes_offset, v7.num_kmers * v7.suffix_width as u64)
            .expect("validated EC v7 suffixes")
    }

    fn v7_ec_ranges(&self) -> &[EcTargetRange] {
        let v7 = self.header_v7.expect("v7 header");
        self.slice(v7.ec_ranges_offset, v7.num_ec_ranges)
            .expect("validated EC v7 target ranges")
    }

    fn v7_ec_id_for_range(&self, target_start: u32, target_len: u32) -> Option<u32> {
        let target_len = u16::try_from(target_len).ok()?;
        self.v7_ec_ranges()
            .iter()
            .position(|range| range.target_start == target_start && range.target_len == target_len)
            .and_then(|idx| u32::try_from(idx).ok())
    }

    fn v6_block_offsets(&self) -> &[u32] {
        let v6 = self.header_v6.expect("v6 header");
        self.slice(v6.block_offsets_offset, v6.num_blocks + 1)
            .expect("validated EC v6 block offsets")
    }

    fn v6_block_max_suffixes(&self) -> &[u8] {
        let v6 = self.header_v6.expect("v6 header");
        self.bytes(
            v6.block_max_suffixes_offset,
            v6.num_blocks * v6.suffix_width as u64,
        )
        .expect("validated EC v6 block max suffixes")
    }

    fn v6_suffixes(&self) -> &[u8] {
        let v6 = self.header_v6.expect("v6 header");
        self.bytes(v6.suffixes_offset, v6.num_kmers * v6.suffix_width as u64)
            .expect("validated EC v6 suffixes")
    }

    #[inline]
    fn v6_find_suffix_in_block(
        &self,
        record_start: usize,
        record_end: usize,
        suffix: u64,
    ) -> Option<usize> {
        let v6 = self.header_v6.expect("v6 header");
        let suffix_width = v6.suffix_width as usize;
        let len = record_end.checked_sub(record_start)?;
        if len == 0 {
            return None;
        }
        if suffix_width == 0 {
            return Some(0);
        }
        let suffixes = self.v6_suffixes();
        let block_base = record_start.checked_mul(suffix_width)?;
        let mut lane = 0;
        while lane + 8 <= len {
            let mut equal_high_bits = 0x8080_8080_8080_8080_u64;
            for byte_idx in 0..suffix_width {
                let plane_start = block_base + byte_idx * len + lane;
                let plane = suffixes.get(plane_start..plane_start + 8)?;
                let word = u64::from_le_bytes(plane.try_into().ok()?);
                let query = u64::from_le_bytes([((suffix >> (byte_idx * 8)) & 0xff) as u8; 8]);
                equal_high_bits &= zero_byte_high_bits(word ^ query);
            }
            while equal_high_bits != 0 {
                let candidate = equal_high_bits.trailing_zeros() as usize / 8;
                if read_v6_planar_suffix(suffixes, block_base, len, lane + candidate, suffix_width)
                    == suffix
                {
                    return Some(lane + candidate);
                }
                equal_high_bits &= !(0x80_u64 << (candidate * 8));
            }
            lane += 8;
        }
        while lane < len {
            if read_v6_planar_suffix(suffixes, block_base, len, lane, suffix_width) == suffix {
                return Some(lane);
            }
            lane += 1;
        }
        None
    }

    #[inline]
    fn v7_find_suffix_in_block(
        &self,
        record_start: usize,
        record_end: usize,
        suffix: u64,
    ) -> Option<usize> {
        let v7 = self.header_v7.expect("v7 header");
        let suffix_width = v7.suffix_width as usize;
        let len = record_end.checked_sub(record_start)?;
        if len == 0 {
            return None;
        }
        if suffix_width == 0 {
            return Some(0);
        }
        let suffixes = self.v7_suffixes();
        let block_base = record_start.checked_mul(suffix_width)?;
        let mut lane = 0;
        while lane + 8 <= len {
            let mut equal_high_bits = 0x8080_8080_8080_8080_u64;
            for byte_idx in 0..suffix_width {
                let plane_start = block_base + byte_idx * len + lane;
                let plane = suffixes.get(plane_start..plane_start + 8)?;
                let word = u64::from_le_bytes(plane.try_into().ok()?);
                let query = u64::from_le_bytes([((suffix >> (byte_idx * 8)) & 0xff) as u8; 8]);
                equal_high_bits &= zero_byte_high_bits(word ^ query);
            }
            while equal_high_bits != 0 {
                let candidate = equal_high_bits.trailing_zeros() as usize / 8;
                if read_v6_planar_suffix(suffixes, block_base, len, lane + candidate, suffix_width)
                    == suffix
                {
                    return Some(lane + candidate);
                }
                equal_high_bits &= !(0x80_u64 << (candidate * 8));
            }
            lane += 8;
        }
        while lane < len {
            if read_v6_planar_suffix(suffixes, block_base, len, lane, suffix_width) == suffix {
                return Some(lane);
            }
            lane += 1;
        }
        None
    }

    fn v4_suffixes(&self) -> &[u8] {
        let v4 = self.header_v4.expect("v4 header");
        self.bytes(v4.suffixes_offset, v4.num_kmers * v4.suffix_width as u64)
            .expect("validated EC v4 suffixes")
    }

    fn v4_ec_ids(&self) -> &[u32] {
        let v4 = self.header_v4.expect("v4 header");
        self.slice(v4.ec_ids_offset, v4.num_kmers)
            .expect("validated EC v4 ec ids")
    }

    fn v4_gene_df(&self, idx: usize) -> u32 {
        let v4 = self.header_v4.expect("v4 header");
        let width = v4.gene_df_width as usize;
        let start = idx * width;
        let dfs = self
            .bytes(v4.gene_dfs_offset, v4.num_kmers * width as u64)
            .expect("validated EC v4 gene dfs");
        match width {
            2 => u16::from_le_bytes([dfs[start], dfs[start + 1]]) as u32,
            4 => u32::from_le_bytes([dfs[start], dfs[start + 1], dfs[start + 2], dfs[start + 3]]),
            _ => unreachable!("validated EC v4 gene df width"),
        }
    }

    fn for_each_v4_ec_target(&self, ec_id: u32, mut f: impl FnMut(GeneId, u8)) {
        let v4 = self.header_v4.expect("v4 header");
        let offsets = self.ec_offsets();
        let start = offsets.get(ec_id as usize).copied().unwrap_or(0) as usize;
        let end = offsets
            .get(ec_id as usize + 1)
            .copied()
            .unwrap_or(start as u32) as usize;
        let targets = self
            .bytes(
                v4.ec_targets_offset,
                v4.num_ec_targets * EC_V4_TARGET_WIDTH as u64,
            )
            .expect("validated EC v4 targets");
        for idx in start..end {
            let offset = idx * EC_V4_TARGET_WIDTH;
            let gene_id = u32::from_le_bytes([
                targets[offset],
                targets[offset + 1],
                targets[offset + 2],
                targets[offset + 3],
            ]);
            f(gene_id, targets[offset + 4]);
        }
    }

    #[inline]
    fn for_each_ec_target_range(
        &self,
        target_start: u32,
        target_len: u32,
        mut f: impl FnMut(GeneId, u8),
    ) {
        let (ec_targets_offset, num_ec_targets, target_width) = if let Some(v8) = self.header_v8 {
            (v8.ec_targets_offset, v8.num_ec_targets, EC_V5_TARGET_WIDTH)
        } else if let Some(v7) = self.header_v7 {
            (
                v7.ec_targets_offset,
                v7.num_ec_targets,
                v7.target_width as usize,
            )
        } else if let Some(v6) = self.header_v6 {
            (v6.ec_targets_offset, v6.num_ec_targets, EC_V5_TARGET_WIDTH)
        } else {
            let v5 = self.header_v5.expect("v5/v6 header");
            (v5.ec_targets_offset, v5.num_ec_targets, EC_V5_TARGET_WIDTH)
        };
        let targets = self
            .bytes(ec_targets_offset, num_ec_targets * target_width as u64)
            .expect("validated EC v5/v6 targets");
        let start = target_start as usize;
        let end = start + target_len as usize;
        for idx in start..end {
            let offset = idx * target_width;
            let gene_id = if target_width == 4 {
                u32::from_le_bytes([targets[offset], targets[offset + 1], targets[offset + 2], 0])
            } else {
                u32::from_le_bytes([
                    targets[offset],
                    targets[offset + 1],
                    targets[offset + 2],
                    targets[offset + 3],
                ])
            };
            f(gene_id, targets[offset + target_width - 1]);
        }
    }
}

fn accumulate_ec_seq(
    index: &LoadedEcIndex,
    encoded: &EncodedSeq,
    stride: u32,
    max_gene_df: u32,
    gene_scores: &mut FxHashMap<GeneId, EcGeneScore>,
    stats: &mut EcMappingStats,
) {
    let use_lookup_hit = matches!(
        index,
        LoadedEcIndex::Mmap(mmap)
            if mmap.header_v5.is_some()
                || mmap.header_v6.is_some()
                || mmap.header_v7.is_some()
                || mmap.header_v8.is_some()
    );
    for kmer in iter_kmers_2bit(encoded, index.k()) {
        if kmer.pos % stride != 0 {
            continue;
        }
        stats.seed_lookups += 1;
        if use_lookup_hit {
            let Some(hit) = index.lookup_hit(kmer.code) else {
                continue;
            };
            index.accumulate_lookup_hit(hit, max_gene_df, gene_scores, stats);
        } else {
            let Some(record) = index.lookup(kmer.code) else {
                continue;
            };
            index.accumulate_lookup_record(record, max_gene_df, gene_scores, stats);
        }
    }
}

fn accumulate_ec_seq_dense(
    index: &LoadedEcIndex,
    encoded: &EncodedSeq,
    stride: u32,
    max_gene_df: u32,
    scratch: &mut EcScoreScratch,
    stats: &mut EcMappingStats,
) {
    let use_lookup_hit = matches!(
        index,
        LoadedEcIndex::Mmap(mmap)
            if mmap.header_v5.is_some()
                || mmap.header_v6.is_some()
                || mmap.header_v7.is_some()
                || mmap.header_v8.is_some()
    );
    for kmer in iter_kmers_2bit(encoded, index.k()) {
        if kmer.pos % stride != 0 {
            continue;
        }
        stats.seed_lookups += 1;
        if use_lookup_hit {
            let Some(hit) = index.lookup_hit(kmer.code) else {
                continue;
            };
            index.visit_lookup_hit(
                hit,
                max_gene_df,
                stats,
                |gene_id, seed_score, class_bits| {
                    scratch.add(gene_id, seed_score, class_bits);
                },
            );
        } else {
            let Some(record) = index.lookup(kmer.code) else {
                continue;
            };
            index.visit_lookup_record(
                record,
                max_gene_df,
                stats,
                |gene_id, seed_score, class_bits| {
                    scratch.add(gene_id, seed_score, class_bits);
                },
            );
        }
    }
}

pub fn finish_ec_mapping(
    gene_scores: FxHashMap<GeneId, EcGeneScore>,
    stats: EcMappingStats,
) -> (EcReadMapping, EcMappingStats) {
    finish_ec_scores(gene_scores.into_values().collect(), stats)
}

fn finish_ec_scores(
    mut gene_scores: Vec<EcGeneScore>,
    mut stats: EcMappingStats,
) -> (EcReadMapping, EcMappingStats) {
    gene_scores.sort_unstable_by_key(|score| {
        (
            std::cmp::Reverse(score.total_score),
            std::cmp::Reverse(score.exon_score),
            std::cmp::Reverse(score.intron_score),
            score.gene_id,
        )
    });
    let best_score = gene_scores.first().map_or(0, |score| score.total_score);
    if best_score > 0 {
        stats.reads_with_gene_support = 1;
    }
    (
        EcReadMapping {
            gene_scores,
            best_score,
            ec_hits: stats.ec_hits.min(u32::MAX as u64) as u32,
        },
        stats,
    )
}

#[inline]
fn add_gene_score(
    gene_scores: &mut FxHashMap<GeneId, EcGeneScore>,
    gene_id: GeneId,
    seed_score: u16,
    class_bits: u8,
) {
    let score = gene_scores.entry(gene_id).or_insert(EcGeneScore {
        gene_id,
        ..EcGeneScore::default()
    });
    add_to_gene_score(score, seed_score, class_bits);
}

#[inline]
fn add_to_gene_score(score: &mut EcGeneScore, seed_score: u16, class_bits: u8) {
    score.total_score = score.total_score.saturating_add(seed_score);
    if class_bits & EC_CLASS_EXON != 0 {
        score.exon_score = score.exon_score.saturating_add(seed_score);
    }
    if class_bits & EC_CLASS_INTRON != 0 {
        score.intron_score = score.intron_score.saturating_add(seed_score);
    }
    if class_bits & EC_CLASS_ANTISENSE != 0 {
        score.antisense_score = score.antisense_score.saturating_add(seed_score);
    }
    if class_bits == EC_CLASS_UNKNOWN {
        score.unknown_score = score.unknown_score.saturating_add(seed_score);
    }
}

fn flush_kmer_gene_run(
    buffer: &mut Vec<KmerGeneRecord>,
    prefix: &Path,
    run_id: usize,
) -> Result<PathBuf, EcIndexError> {
    buffer.sort_unstable();
    buffer.dedup();
    let path = prefix.with_extension(format!("run{run_id}.bin"));
    let mut writer = BufWriter::new(File::create(&path)?);
    for record in buffer.iter().copied() {
        write_u64(&mut writer, record.kmer_code)?;
        write_u32(&mut writer, record.gene_id)?;
        writer.write_all(&[record.class_bits])?;
    }
    buffer.clear();
    Ok(path)
}

fn merge_kmer_gene_runs(
    runs: &[PathBuf],
    mut emit: impl FnMut(u64, Vec<EcGeneTarget>) -> Result<(), EcIndexError>,
) -> Result<(), EcIndexError> {
    let mut readers = Vec::with_capacity(runs.len());
    let mut heap = BinaryHeap::new();
    for (idx, path) in runs.iter().enumerate() {
        let mut reader = BufReader::new(File::open(path)?);
        if let Some(record) = read_kmer_gene_record(&mut reader)? {
            heap.push(Reverse((record, idx)));
        }
        readers.push(reader);
    }

    while let Some(Reverse((first, first_idx))) = heap.pop() {
        let kmer_code = first.kmer_code;
        let mut gene_class_bits: BTreeMap<GeneId, u8> = BTreeMap::new();
        *gene_class_bits.entry(first.gene_id).or_default() |= first.class_bits;
        if let Some(record) = read_kmer_gene_record(&mut readers[first_idx])? {
            heap.push(Reverse((record, first_idx)));
        }

        while heap
            .peek()
            .is_some_and(|Reverse((record, _))| record.kmer_code == kmer_code)
        {
            let Reverse((record, idx)) = heap.pop().expect("heap had kmer");
            *gene_class_bits.entry(record.gene_id).or_default() |= record.class_bits;
            if let Some(next) = read_kmer_gene_record(&mut readers[idx])? {
                heap.push(Reverse((next, idx)));
            }
        }
        let targets = gene_class_bits
            .into_iter()
            .map(|(gene_id, class_bits)| EcGeneTarget {
                gene_id,
                class_bits,
            })
            .collect();
        emit(kmer_code, targets)?;
    }

    Ok(())
}

fn read_kmer_gene_record(reader: &mut impl Read) -> Result<Option<KmerGeneRecord>, EcIndexError> {
    let mut kmer_bytes = [0_u8; 8];
    match reader.read_exact(&mut kmer_bytes) {
        Ok(()) => {}
        Err(err) if err.kind() == std::io::ErrorKind::UnexpectedEof => return Ok(None),
        Err(err) => return Err(err.into()),
    }
    let gene_id = read_u32(reader)?;
    let mut class = [0_u8; 1];
    reader.read_exact(&mut class)?;
    Ok(Some(KmerGeneRecord {
        kmer_code: u64::from_le_bytes(kmer_bytes),
        gene_id,
        class_bits: class[0],
    }))
}

fn load_t2g_target_info(
    path: impl AsRef<Path>,
) -> Result<FxHashMap<String, TargetInfo>, EcIndexError> {
    let reader = BufReader::new(File::open(path)?);
    let mut out = FxHashMap::default();
    for line in reader.lines() {
        let line = line?;
        if line.trim().is_empty() || line.starts_with('#') {
            continue;
        }
        let fields = line.split('\t').collect::<Vec<_>>();
        if fields.len() < 2 {
            continue;
        }
        let class_bits = fields
            .get(2)
            .map(|value| class_bits_from_token(value))
            .unwrap_or(EC_CLASS_UNKNOWN);
        out.insert(
            fields[0].to_owned(),
            TargetInfo {
                gene: fields[1].to_owned(),
                class_bits,
            },
        );
    }
    Ok(out)
}

fn parse_target_name(header: &str) -> &str {
    header.split_whitespace().next().unwrap_or(header)
}

fn infer_class_bits(header: &str) -> u8 {
    let lower = header.to_ascii_lowercase();
    if lower.contains("antisense") || lower.ends_with("-a") {
        EC_CLASS_ANTISENSE
    } else if lower.contains("intron") || lower.contains("|u") || lower.ends_with("-u") {
        EC_CLASS_INTRON
    } else if lower.contains("exon") || lower.contains("transcript") || lower.ends_with("-s") {
        EC_CLASS_EXON
    } else {
        EC_CLASS_UNKNOWN
    }
}

fn class_bits_from_token(token: &str) -> u8 {
    match token.trim().to_ascii_uppercase().as_str() {
        "S" | "EXON" | "SPLICED" | "TRANSCRIPT" => EC_CLASS_EXON,
        "U" | "INTRON" | "UNSPLICED" => EC_CLASS_INTRON,
        "A" | "ANTISENSE" => EC_CLASS_ANTISENSE,
        _ => EC_CLASS_UNKNOWN,
    }
}

fn parse_gene_from_header(header: &str) -> String {
    let name = header
        .split_whitespace()
        .next()
        .unwrap_or(header)
        .to_owned();
    header
        .split(['|', ' ', '\t'])
        .find_map(|part| {
            part.strip_prefix("gene=")
                .or_else(|| part.strip_prefix("gene:"))
                .or_else(|| part.strip_prefix("gene_id="))
                .or_else(|| part.strip_prefix("gene_id:"))
        })
        .map(clean_gene_token)
        .or_else(|| parse_quoted_gene_id(header))
        .unwrap_or(name)
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

fn unique_temp_suffix() -> u128 {
    std::time::SystemTime::now()
        .duration_since(std::time::UNIX_EPOCH)
        .map(|duration| duration.as_nanos())
        .unwrap_or(0)
}

fn ec_seed_score(num_genes: usize, gene_df: u32) -> u16 {
    if gene_df == 0 {
        return 0;
    }
    let idf = ((num_genes as f64 + 1.0) / (gene_df as f64 + 1.0)).ln();
    (1.0 + idf * 64.0).round().clamp(1.0, u16::MAX as f64) as u16
}

#[inline]
fn build_seed_score_table(num_genes: usize) -> Vec<u16> {
    (0..=num_genes)
        .map(|gene_df| ec_seed_score(num_genes, gene_df as u32))
        .collect()
}

fn build_prefix_offsets(
    kmers: &[EcKmerRecord],
    k: u8,
    prefix_bits: u8,
) -> Result<Vec<u32>, EcIndexError> {
    let buckets = 1_usize
        .checked_shl(prefix_bits as u32)
        .ok_or(EcIndexError::TooLarge("EC prefix directory"))?;
    let mut offsets = vec![0_u32; buckets + 1];
    let mut current_prefix = 0_usize;
    for (idx, record) in kmers.iter().enumerate() {
        let prefix = ec_kmer_prefix(record.kmer_code, k, prefix_bits) as usize;
        while current_prefix <= prefix {
            offsets[current_prefix] =
                u32::try_from(idx).map_err(|_| EcIndexError::TooLarge("EC kmer offsets"))?;
            current_prefix += 1;
        }
    }
    let end = u32::try_from(kmers.len()).map_err(|_| EcIndexError::TooLarge("EC kmer offsets"))?;
    while current_prefix <= buckets {
        offsets[current_prefix] = end;
        current_prefix += 1;
    }
    Ok(offsets)
}

fn save_mmap_locality_from_slices(
    path: impl AsRef<Path>,
    k: u8,
    gene_metas: &[EcGeneMeta],
    gene_names: &[u8],
    kmers: &[EcKmerRecord],
    ec_offsets: &[u32],
    ec_gene_ids: &[GeneId],
    ec_class_bits: &[u8],
    prefix_offsets: &[u32],
) -> Result<(), EcIndexError> {
    if ec_mph_enabled() && can_save_mmap_v5_from_slices(k, kmers, ec_offsets) {
        save_mmap_v8_from_slices(
            path,
            k,
            gene_metas,
            gene_names,
            kmers,
            ec_offsets,
            ec_gene_ids,
            ec_class_bits,
            prefix_offsets,
        )
    } else if ec_target_dedup_enabled() && can_save_mmap_v5_from_slices(k, kmers, ec_offsets) {
        if ec_target_dedup_is_smaller(ec_offsets, ec_gene_ids, ec_class_bits)? {
            save_mmap_v7_from_slices(
                path,
                k,
                gene_metas,
                gene_names,
                kmers,
                ec_offsets,
                ec_gene_ids,
                ec_class_bits,
                prefix_offsets,
            )
        } else {
            save_mmap_v6_from_slices(
                path,
                k,
                gene_metas,
                gene_names,
                kmers,
                ec_offsets,
                ec_gene_ids,
                ec_class_bits,
                prefix_offsets,
            )
        }
    } else if ec_mmap_v6_enabled() && can_save_mmap_v5_from_slices(k, kmers, ec_offsets) {
        save_mmap_v6_from_slices(
            path,
            k,
            gene_metas,
            gene_names,
            kmers,
            ec_offsets,
            ec_gene_ids,
            ec_class_bits,
            prefix_offsets,
        )
    } else if ec_mmap_v5_enabled() && can_save_mmap_v5_from_slices(k, kmers, ec_offsets) {
        save_mmap_v5_from_slices(
            path,
            k,
            gene_metas,
            gene_names,
            kmers,
            ec_offsets,
            ec_gene_ids,
            ec_class_bits,
            prefix_offsets,
        )
    } else {
        save_mmap_v4_from_slices(
            path,
            k,
            gene_metas,
            gene_names,
            kmers,
            ec_offsets,
            ec_gene_ids,
            ec_class_bits,
            prefix_offsets,
        )
    }
}

fn ec_mph_enabled() -> bool {
    matches!(
        std::env::var("CONSTELLATION_EC_MPH").as_deref(),
        Ok("1") | Ok("true") | Ok("TRUE") | Ok("yes") | Ok("YES")
    )
}

fn ec_target_dedup_enabled() -> bool {
    matches!(
        std::env::var("CONSTELLATION_EC_TARGET_DEDUP").as_deref(),
        Ok("1") | Ok("true") | Ok("TRUE") | Ok("yes") | Ok("YES")
    )
}

fn ec_mmap_v6_enabled() -> bool {
    matches!(
        std::env::var("CONSTELLATION_EC_MMAP_V6").as_deref(),
        Ok("1") | Ok("true") | Ok("TRUE") | Ok("yes") | Ok("YES")
    )
}

fn ec_mmap_v5_enabled() -> bool {
    matches!(
        std::env::var("CONSTELLATION_EC_MMAP_V5").as_deref(),
        Ok("1") | Ok("true") | Ok("TRUE") | Ok("yes") | Ok("YES")
    )
}

fn can_save_mmap_v5_from_slices(k: u8, kmers: &[EcKmerRecord], ec_offsets: &[u32]) -> bool {
    let prefix_bits = EC_PREFIX_BITS.min(k.saturating_mul(2));
    let suffix_bits = k.saturating_mul(2).saturating_sub(prefix_bits);
    suffix_bits <= 48
        && kmers.iter().all(|record| {
            record.gene_df <= u16::MAX as u32
                && ec_offsets
                    .get(record.ec_id as usize)
                    .zip(ec_offsets.get(record.ec_id as usize + 1))
                    .map(|(&start, &end)| end >= start && end - start <= u16::MAX as u32)
                    .unwrap_or(false)
        })
}

struct DedupedEcTargets {
    ec_ranges: Vec<EcTargetRange>,
    gene_ids: Vec<GeneId>,
    class_bits: Vec<u8>,
}

struct MphEcLayout {
    slot_offsets: Vec<u32>,
    group_offsets: Vec<u32>,
    group_displacements: Vec<u32>,
    suffixes: Vec<u64>,
    payloads: Vec<EcV5KmerPayload>,
}

fn save_mmap_v8_from_slices(
    path: impl AsRef<Path>,
    k: u8,
    gene_metas: &[EcGeneMeta],
    gene_names: &[u8],
    kmers: &[EcKmerRecord],
    ec_offsets: &[u32],
    ec_gene_ids: &[GeneId],
    ec_class_bits: &[u8],
    kmer_prefix_offsets: &[u32],
) -> Result<(), EcIndexError> {
    let prefix_bits = EC_PREFIX_BITS.min(k.saturating_mul(2));
    let suffix_bits = k.saturating_mul(2).saturating_sub(prefix_bits);
    let suffix_width = suffix_bits.div_ceil(8) as usize;
    let layout = build_mph_ec_layout(k, prefix_bits, kmers, ec_offsets, kmer_prefix_offsets)?;
    let num_slots = layout.suffixes.len() as u64;
    let num_groups = layout.group_displacements.len() as u64;
    let num_targets = ec_gene_ids.len() as u64;

    let mut cursor = align_up(std::mem::size_of::<EcMmapHeaderV8>() as u64, 8);
    let slot_offsets_offset = cursor;
    cursor += byte_len::<u32>(layout.slot_offsets.len());
    cursor = align_up(cursor, 4);
    let group_offsets_offset = cursor;
    cursor += byte_len::<u32>(layout.group_offsets.len());
    cursor = align_up(cursor, 4);
    let group_displacements_offset = cursor;
    cursor += byte_len::<u32>(layout.group_displacements.len());
    cursor = align_up(cursor, 8);
    let genes_offset = cursor;
    cursor += byte_len::<EcGeneMeta>(gene_metas.len());
    cursor = align_up(cursor, 8);
    let suffixes_offset = cursor;
    cursor += checked_mul_u64(num_slots, suffix_width as u64, "EC v8 suffixes")?;
    cursor = align_up(cursor, 8);
    let payloads_offset = cursor;
    cursor += byte_len::<EcV5KmerPayload>(layout.payloads.len());
    cursor = align_up(cursor, 4);
    let ec_offsets_offset = cursor;
    cursor += byte_len::<u32>(ec_offsets.len());
    cursor = align_up(cursor, 4);
    let ec_targets_offset = cursor;
    cursor += checked_mul_u64(num_targets, EC_V5_TARGET_WIDTH as u64, "EC v8 targets")?;
    cursor = align_up(cursor, 8);
    let gene_names_offset = cursor;

    let header = EcMmapHeaderV8 {
        magic: EC_MMAP_MAGIC,
        version: EC_MMAP_VERSION_V8,
        k: k as u32,
        prefix_bits: prefix_bits as u32,
        suffix_bits: suffix_bits as u32,
        suffix_width: suffix_width as u32,
        group_size: EC_MPH_GROUP_SIZE as u32,
        num_genes: gene_metas.len() as u64,
        num_slots,
        num_groups,
        num_ec_offsets: ec_offsets.len() as u64,
        num_ec_targets: num_targets,
        genes_offset,
        slot_offsets_offset,
        slot_offsets_len: byte_len::<u32>(layout.slot_offsets.len()),
        group_offsets_offset,
        group_offsets_len: byte_len::<u32>(layout.group_offsets.len()),
        group_displacements_offset,
        suffixes_offset,
        payloads_offset,
        ec_offsets_offset,
        ec_targets_offset,
        gene_names_offset,
        gene_names_len: gene_names.len() as u64,
    };

    let mut file = File::create(path)?;
    file.write_all(bytemuck::bytes_of(&header))?;
    write_padding(
        &mut file,
        std::mem::size_of::<EcMmapHeaderV8>() as u64,
        slot_offsets_offset,
    )?;
    write_pod_slice(&mut file, &layout.slot_offsets)?;
    write_padding(
        &mut file,
        slot_offsets_offset + byte_len::<u32>(layout.slot_offsets.len()),
        group_offsets_offset,
    )?;
    write_pod_slice(&mut file, &layout.group_offsets)?;
    write_padding(
        &mut file,
        group_offsets_offset + byte_len::<u32>(layout.group_offsets.len()),
        group_displacements_offset,
    )?;
    write_pod_slice(&mut file, &layout.group_displacements)?;
    write_padding(
        &mut file,
        group_displacements_offset + byte_len::<u32>(layout.group_displacements.len()),
        genes_offset,
    )?;
    write_pod_slice(&mut file, gene_metas)?;
    write_padding(
        &mut file,
        genes_offset + byte_len::<EcGeneMeta>(gene_metas.len()),
        suffixes_offset,
    )?;
    write_packed_suffix_values(&mut file, &layout.suffixes, suffix_width)?;
    write_padding(
        &mut file,
        suffixes_offset + checked_mul_u64(num_slots, suffix_width as u64, "EC v8 suffixes")?,
        payloads_offset,
    )?;
    write_pod_slice(&mut file, &layout.payloads)?;
    write_padding(
        &mut file,
        payloads_offset + byte_len::<EcV5KmerPayload>(layout.payloads.len()),
        ec_offsets_offset,
    )?;
    write_pod_slice(&mut file, ec_offsets)?;
    write_padding(
        &mut file,
        ec_offsets_offset + byte_len::<u32>(ec_offsets.len()),
        ec_targets_offset,
    )?;
    write_v4_targets(&mut file, ec_gene_ids, ec_class_bits)?;
    write_padding(
        &mut file,
        ec_targets_offset
            + checked_mul_u64(num_targets, EC_V5_TARGET_WIDTH as u64, "EC v8 targets")?,
        gene_names_offset,
    )?;
    file.write_all(gene_names)?;
    Ok(())
}

fn save_mmap_v7_from_slices(
    path: impl AsRef<Path>,
    k: u8,
    gene_metas: &[EcGeneMeta],
    gene_names: &[u8],
    kmers: &[EcKmerRecord],
    ec_offsets: &[u32],
    ec_gene_ids: &[GeneId],
    ec_class_bits: &[u8],
    kmer_prefix_offsets: &[u32],
) -> Result<(), EcIndexError> {
    let prefix_bits = EC_PREFIX_BITS.min(k.saturating_mul(2));
    let suffix_bits = k.saturating_mul(2).saturating_sub(prefix_bits);
    let suffix_width = suffix_bits.div_ceil(8) as usize;
    let block_size = EC_V6_BLOCK_SIZE;
    let (prefix_offsets, block_offsets, block_max_suffixes) = build_v6_blocks(
        k,
        prefix_bits,
        suffix_width,
        block_size,
        kmers,
        kmer_prefix_offsets,
    )?;
    let deduped_targets = dedup_ec_targets(ec_offsets, ec_gene_ids, ec_class_bits)?;
    let num_kmers = kmers.len() as u64;
    let num_blocks = block_max_suffixes.len() as u64;
    let num_targets = deduped_targets.gene_ids.len() as u64;
    let target_width = ec_target_width(&deduped_targets.gene_ids);

    let mut cursor = align_up(std::mem::size_of::<EcMmapHeaderV7>() as u64, 8);
    let prefix_offsets_offset = cursor;
    cursor += byte_len::<u32>(prefix_offsets.len());
    cursor = align_up(cursor, 8);
    let genes_offset = cursor;
    cursor += byte_len::<EcGeneMeta>(gene_metas.len());
    cursor = align_up(cursor, 4);
    let block_offsets_offset = cursor;
    cursor += byte_len::<u32>(block_offsets.len());
    cursor = align_up(cursor, 8);
    let block_max_suffixes_offset = cursor;
    cursor += checked_mul_u64(num_blocks, suffix_width as u64, "EC v7 block max suffixes")?;
    cursor = align_up(cursor, 8);
    let suffixes_offset = cursor;
    cursor += checked_mul_u64(num_kmers, suffix_width as u64, "EC v7 suffixes")?;
    cursor = align_up(cursor, 8);
    let payloads_offset = cursor;
    cursor += byte_len::<EcV5KmerPayload>(kmers.len());
    cursor = align_up(cursor, 8);
    let ec_ranges_offset = cursor;
    cursor += byte_len::<EcTargetRange>(deduped_targets.ec_ranges.len());
    cursor = align_up(cursor, 4);
    let ec_targets_offset = cursor;
    cursor += checked_mul_u64(num_targets, target_width as u64, "EC v7 targets")?;
    cursor = align_up(cursor, 8);
    let gene_names_offset = cursor;

    let header = EcMmapHeaderV7 {
        magic: EC_MMAP_MAGIC,
        version: EC_MMAP_VERSION,
        k: k as u32,
        prefix_bits: prefix_bits as u32,
        suffix_bits: suffix_bits as u32,
        suffix_width: suffix_width as u32,
        block_size: block_size as u32,
        target_width: target_width as u32,
        _pad0: 0,
        num_genes: gene_metas.len() as u64,
        num_kmers,
        num_blocks,
        num_ec_ranges: deduped_targets.ec_ranges.len() as u64,
        num_ec_targets: num_targets,
        genes_offset,
        prefix_offsets_offset,
        prefix_offsets_len: byte_len::<u32>(prefix_offsets.len()),
        block_offsets_offset,
        block_max_suffixes_offset,
        suffixes_offset,
        payloads_offset,
        ec_ranges_offset,
        ec_targets_offset,
        gene_names_offset,
        gene_names_len: gene_names.len() as u64,
    };

    let mut file = File::create(path)?;
    file.write_all(bytemuck::bytes_of(&header))?;
    write_padding(
        &mut file,
        std::mem::size_of::<EcMmapHeaderV7>() as u64,
        prefix_offsets_offset,
    )?;
    write_pod_slice(&mut file, &prefix_offsets)?;
    write_padding(
        &mut file,
        prefix_offsets_offset + byte_len::<u32>(prefix_offsets.len()),
        genes_offset,
    )?;
    write_pod_slice(&mut file, gene_metas)?;
    write_padding(
        &mut file,
        genes_offset + byte_len::<EcGeneMeta>(gene_metas.len()),
        block_offsets_offset,
    )?;
    write_pod_slice(&mut file, &block_offsets)?;
    write_padding(
        &mut file,
        block_offsets_offset + byte_len::<u32>(block_offsets.len()),
        block_max_suffixes_offset,
    )?;
    write_v6_block_max_suffixes(&mut file, &block_max_suffixes, suffix_width)?;
    write_padding(
        &mut file,
        block_max_suffixes_offset
            + checked_mul_u64(num_blocks, suffix_width as u64, "EC v7 block max suffixes")?,
        suffixes_offset,
    )?;
    write_v6_suffix_blocks(
        &mut file,
        kmers,
        k,
        prefix_bits,
        suffix_width,
        &block_offsets,
    )?;
    write_padding(
        &mut file,
        suffixes_offset + checked_mul_u64(num_kmers, suffix_width as u64, "EC v7 suffixes")?,
        payloads_offset,
    )?;
    write_v5_payloads_from_ranges(&mut file, kmers, &deduped_targets.ec_ranges)?;
    write_padding(
        &mut file,
        payloads_offset + byte_len::<EcV5KmerPayload>(kmers.len()),
        ec_ranges_offset,
    )?;
    write_pod_slice(&mut file, &deduped_targets.ec_ranges)?;
    write_padding(
        &mut file,
        ec_ranges_offset + byte_len::<EcTargetRange>(deduped_targets.ec_ranges.len()),
        ec_targets_offset,
    )?;
    write_v7_targets(
        &mut file,
        &deduped_targets.gene_ids,
        &deduped_targets.class_bits,
        target_width,
    )?;
    write_padding(
        &mut file,
        ec_targets_offset + checked_mul_u64(num_targets, target_width as u64, "EC v7 targets")?,
        gene_names_offset,
    )?;
    file.write_all(gene_names)?;
    Ok(())
}

fn save_mmap_v6_from_slices(
    path: impl AsRef<Path>,
    k: u8,
    gene_metas: &[EcGeneMeta],
    gene_names: &[u8],
    kmers: &[EcKmerRecord],
    ec_offsets: &[u32],
    ec_gene_ids: &[GeneId],
    ec_class_bits: &[u8],
    kmer_prefix_offsets: &[u32],
) -> Result<(), EcIndexError> {
    let prefix_bits = EC_PREFIX_BITS.min(k.saturating_mul(2));
    let suffix_bits = k.saturating_mul(2).saturating_sub(prefix_bits);
    let suffix_width = suffix_bits.div_ceil(8) as usize;
    let block_size = EC_V6_BLOCK_SIZE;
    let (prefix_offsets, block_offsets, block_max_suffixes) = build_v6_blocks(
        k,
        prefix_bits,
        suffix_width,
        block_size,
        kmers,
        kmer_prefix_offsets,
    )?;
    let num_kmers = kmers.len() as u64;
    let num_blocks = block_max_suffixes.len() as u64;
    let num_targets = ec_gene_ids.len() as u64;

    let mut cursor = align_up(std::mem::size_of::<EcMmapHeaderV6>() as u64, 8);
    let prefix_offsets_offset = cursor;
    cursor += byte_len::<u32>(prefix_offsets.len());
    cursor = align_up(cursor, 8);
    let genes_offset = cursor;
    cursor += byte_len::<EcGeneMeta>(gene_metas.len());
    cursor = align_up(cursor, 4);
    let block_offsets_offset = cursor;
    cursor += byte_len::<u32>(block_offsets.len());
    cursor = align_up(cursor, 8);
    let block_max_suffixes_offset = cursor;
    cursor += checked_mul_u64(num_blocks, suffix_width as u64, "EC v6 block max suffixes")?;
    cursor = align_up(cursor, 8);
    let suffixes_offset = cursor;
    cursor += checked_mul_u64(num_kmers, suffix_width as u64, "EC v6 suffixes")?;
    cursor = align_up(cursor, 8);
    let payloads_offset = cursor;
    cursor += byte_len::<EcV5KmerPayload>(kmers.len());
    cursor = align_up(cursor, 4);
    let ec_offsets_offset = cursor;
    cursor += byte_len::<u32>(ec_offsets.len());
    cursor = align_up(cursor, 4);
    let ec_targets_offset = cursor;
    cursor += checked_mul_u64(num_targets, EC_V5_TARGET_WIDTH as u64, "EC v6 targets")?;
    cursor = align_up(cursor, 8);
    let gene_names_offset = cursor;

    let header = EcMmapHeaderV6 {
        magic: EC_MMAP_MAGIC,
        version: EC_MMAP_VERSION_V6,
        k: k as u32,
        prefix_bits: prefix_bits as u32,
        suffix_bits: suffix_bits as u32,
        suffix_width: suffix_width as u32,
        block_size: block_size as u32,
        num_genes: gene_metas.len() as u64,
        num_kmers,
        num_blocks,
        num_ec_offsets: ec_offsets.len() as u64,
        num_ec_targets: num_targets,
        genes_offset,
        prefix_offsets_offset,
        prefix_offsets_len: byte_len::<u32>(prefix_offsets.len()),
        block_offsets_offset,
        block_max_suffixes_offset,
        suffixes_offset,
        payloads_offset,
        ec_offsets_offset,
        ec_targets_offset,
        gene_names_offset,
        gene_names_len: gene_names.len() as u64,
    };

    let mut file = File::create(path)?;
    file.write_all(bytemuck::bytes_of(&header))?;
    write_padding(
        &mut file,
        std::mem::size_of::<EcMmapHeaderV6>() as u64,
        prefix_offsets_offset,
    )?;
    write_pod_slice(&mut file, &prefix_offsets)?;
    write_padding(
        &mut file,
        prefix_offsets_offset + byte_len::<u32>(prefix_offsets.len()),
        genes_offset,
    )?;
    write_pod_slice(&mut file, gene_metas)?;
    write_padding(
        &mut file,
        genes_offset + byte_len::<EcGeneMeta>(gene_metas.len()),
        block_offsets_offset,
    )?;
    write_pod_slice(&mut file, &block_offsets)?;
    write_padding(
        &mut file,
        block_offsets_offset + byte_len::<u32>(block_offsets.len()),
        block_max_suffixes_offset,
    )?;
    write_v6_block_max_suffixes(&mut file, &block_max_suffixes, suffix_width)?;
    write_padding(
        &mut file,
        block_max_suffixes_offset
            + checked_mul_u64(num_blocks, suffix_width as u64, "EC v6 block max suffixes")?,
        suffixes_offset,
    )?;
    write_v6_suffix_blocks(
        &mut file,
        kmers,
        k,
        prefix_bits,
        suffix_width,
        &block_offsets,
    )?;
    write_padding(
        &mut file,
        suffixes_offset + checked_mul_u64(num_kmers, suffix_width as u64, "EC v6 suffixes")?,
        payloads_offset,
    )?;
    write_v5_payloads(&mut file, kmers, ec_offsets)?;
    write_padding(
        &mut file,
        payloads_offset + byte_len::<EcV5KmerPayload>(kmers.len()),
        ec_offsets_offset,
    )?;
    write_pod_slice(&mut file, ec_offsets)?;
    write_padding(
        &mut file,
        ec_offsets_offset + byte_len::<u32>(ec_offsets.len()),
        ec_targets_offset,
    )?;
    write_v4_targets(&mut file, ec_gene_ids, ec_class_bits)?;
    write_padding(
        &mut file,
        ec_targets_offset
            + checked_mul_u64(num_targets, EC_V5_TARGET_WIDTH as u64, "EC v6 targets")?,
        gene_names_offset,
    )?;
    file.write_all(gene_names)?;
    Ok(())
}

fn save_mmap_v5_from_slices(
    path: impl AsRef<Path>,
    k: u8,
    gene_metas: &[EcGeneMeta],
    gene_names: &[u8],
    kmers: &[EcKmerRecord],
    ec_offsets: &[u32],
    ec_gene_ids: &[GeneId],
    ec_class_bits: &[u8],
    prefix_offsets: &[u32],
) -> Result<(), EcIndexError> {
    let prefix_bits = EC_PREFIX_BITS.min(k.saturating_mul(2));
    let suffix_bits = k.saturating_mul(2).saturating_sub(prefix_bits);
    let suffix_width = suffix_bits.div_ceil(8) as usize;
    let num_kmers = kmers.len() as u64;
    let num_targets = ec_gene_ids.len() as u64;

    let mut cursor = align_up(std::mem::size_of::<EcMmapHeaderV5>() as u64, 8);
    let prefix_offsets_offset = cursor;
    cursor += byte_len::<u32>(prefix_offsets.len());
    cursor = align_up(cursor, 8);
    let genes_offset = cursor;
    cursor += byte_len::<EcGeneMeta>(gene_metas.len());
    cursor = align_up(cursor, 8);
    let suffixes_offset = cursor;
    cursor += checked_mul_u64(num_kmers, suffix_width as u64, "EC v5 suffixes")?;
    cursor = align_up(cursor, 8);
    let payloads_offset = cursor;
    cursor += byte_len::<EcV5KmerPayload>(kmers.len());
    cursor = align_up(cursor, 4);
    let ec_offsets_offset = cursor;
    cursor += byte_len::<u32>(ec_offsets.len());
    cursor = align_up(cursor, 4);
    let ec_targets_offset = cursor;
    cursor += checked_mul_u64(num_targets, EC_V5_TARGET_WIDTH as u64, "EC v5 targets")?;
    cursor = align_up(cursor, 8);
    let gene_names_offset = cursor;

    let header = EcMmapHeaderV5 {
        magic: EC_MMAP_MAGIC,
        version: EC_MMAP_VERSION_V5,
        k: k as u32,
        prefix_bits: prefix_bits as u32,
        suffix_bits: suffix_bits as u32,
        suffix_width: suffix_width as u32,
        _pad1: 0,
        num_genes: gene_metas.len() as u64,
        num_kmers,
        num_ec_offsets: ec_offsets.len() as u64,
        num_ec_targets: num_targets,
        genes_offset,
        suffixes_offset,
        payloads_offset,
        ec_offsets_offset,
        ec_targets_offset,
        gene_names_offset,
        gene_names_len: gene_names.len() as u64,
        prefix_offsets_offset,
        prefix_offsets_len: byte_len::<u32>(prefix_offsets.len()),
    };

    let mut file = File::create(path)?;
    file.write_all(bytemuck::bytes_of(&header))?;
    write_padding(
        &mut file,
        std::mem::size_of::<EcMmapHeaderV5>() as u64,
        prefix_offsets_offset,
    )?;
    write_pod_slice(&mut file, prefix_offsets)?;
    write_padding(
        &mut file,
        prefix_offsets_offset + byte_len::<u32>(prefix_offsets.len()),
        genes_offset,
    )?;
    write_pod_slice(&mut file, gene_metas)?;
    write_padding(
        &mut file,
        genes_offset + byte_len::<EcGeneMeta>(gene_metas.len()),
        suffixes_offset,
    )?;
    write_v4_suffixes(&mut file, kmers, k, prefix_bits, suffix_width)?;
    write_padding(
        &mut file,
        suffixes_offset + checked_mul_u64(num_kmers, suffix_width as u64, "EC v5 suffixes")?,
        payloads_offset,
    )?;
    write_v5_payloads(&mut file, kmers, ec_offsets)?;
    write_padding(
        &mut file,
        payloads_offset + byte_len::<EcV5KmerPayload>(kmers.len()),
        ec_offsets_offset,
    )?;
    write_pod_slice(&mut file, ec_offsets)?;
    write_padding(
        &mut file,
        ec_offsets_offset + byte_len::<u32>(ec_offsets.len()),
        ec_targets_offset,
    )?;
    write_v4_targets(&mut file, ec_gene_ids, ec_class_bits)?;
    write_padding(
        &mut file,
        ec_targets_offset
            + checked_mul_u64(num_targets, EC_V5_TARGET_WIDTH as u64, "EC v5 targets")?,
        gene_names_offset,
    )?;
    file.write_all(gene_names)?;
    Ok(())
}

fn write_v5_payloads(
    writer: &mut impl Write,
    kmers: &[EcKmerRecord],
    ec_offsets: &[u32],
) -> Result<(), EcIndexError> {
    let mut buffer = Vec::with_capacity(1 << 20);
    for record in kmers {
        let start = ec_offsets
            .get(record.ec_id as usize)
            .copied()
            .ok_or(EcIndexError::Truncated)?;
        let end = ec_offsets
            .get(record.ec_id as usize + 1)
            .copied()
            .ok_or(EcIndexError::Truncated)?;
        let payload = EcV5KmerPayload {
            target_start: start,
            gene_df: u16::try_from(record.gene_df)
                .map_err(|_| EcIndexError::TooLarge("EC v5 gene df"))?,
            target_len: u16::try_from(end.saturating_sub(start))
                .map_err(|_| EcIndexError::TooLarge("EC v5 target length"))?,
        };
        buffer.extend_from_slice(bytemuck::bytes_of(&payload));
        if buffer.len() >= 1 << 20 {
            writer.write_all(&buffer)?;
            buffer.clear();
        }
    }
    if !buffer.is_empty() {
        writer.write_all(&buffer)?;
    }
    Ok(())
}

fn write_v5_payloads_from_ranges(
    writer: &mut impl Write,
    kmers: &[EcKmerRecord],
    ec_ranges: &[EcTargetRange],
) -> Result<(), EcIndexError> {
    let mut buffer = Vec::with_capacity(1 << 20);
    for record in kmers {
        let range = ec_ranges
            .get(record.ec_id as usize)
            .ok_or(EcIndexError::Truncated)?;
        let payload = EcV5KmerPayload {
            target_start: range.target_start,
            gene_df: u16::try_from(record.gene_df)
                .map_err(|_| EcIndexError::TooLarge("EC v7 gene df"))?,
            target_len: range.target_len,
        };
        buffer.extend_from_slice(bytemuck::bytes_of(&payload));
        if buffer.len() >= 1 << 20 {
            writer.write_all(&buffer)?;
            buffer.clear();
        }
    }
    if !buffer.is_empty() {
        writer.write_all(&buffer)?;
    }
    Ok(())
}

fn build_mph_ec_layout(
    k: u8,
    prefix_bits: u8,
    kmers: &[EcKmerRecord],
    ec_offsets: &[u32],
    kmer_prefix_offsets: &[u32],
) -> Result<MphEcLayout, EcIndexError> {
    let mut slot_offsets = Vec::with_capacity(kmer_prefix_offsets.len());
    let mut group_offsets = Vec::with_capacity(kmer_prefix_offsets.len());
    let mut total_slots = 0_usize;
    let mut total_groups = 0_usize;
    for bucket in kmer_prefix_offsets.windows(2) {
        slot_offsets
            .push(u32::try_from(total_slots).map_err(|_| EcIndexError::TooLarge("EC v8 slots"))?);
        group_offsets
            .push(u32::try_from(total_groups).map_err(|_| EcIndexError::TooLarge("EC v8 groups"))?);
        let len = (bucket[1] as usize).saturating_sub(bucket[0] as usize);
        if len > 0 {
            total_slots += mph_slot_count(len);
            total_groups += len.div_ceil(EC_MPH_GROUP_SIZE);
        }
    }
    slot_offsets
        .push(u32::try_from(total_slots).map_err(|_| EcIndexError::TooLarge("EC v8 slots"))?);
    group_offsets
        .push(u32::try_from(total_groups).map_err(|_| EcIndexError::TooLarge("EC v8 groups"))?);

    let mut group_displacements = vec![0_u32; total_groups];
    let mut suffixes = vec![0_u64; total_slots];
    let mut payloads = vec![EcV5KmerPayload::default(); total_slots];
    for prefix in 0..kmer_prefix_offsets.len().saturating_sub(1) {
        let record_start = kmer_prefix_offsets[prefix] as usize;
        let record_end = kmer_prefix_offsets[prefix + 1] as usize;
        let len = record_end.saturating_sub(record_start);
        if len == 0 {
            continue;
        }
        let slot_start = slot_offsets[prefix] as usize;
        let slot_len = slot_offsets[prefix + 1] as usize - slot_start;
        let group_start = group_offsets[prefix] as usize;
        let group_count = group_offsets[prefix + 1] as usize - group_start;
        build_mph_bucket(
            k,
            prefix_bits,
            &kmers[record_start..record_end],
            ec_offsets,
            &mut group_displacements[group_start..group_start + group_count],
            &mut suffixes[slot_start..slot_start + slot_len],
            &mut payloads[slot_start..slot_start + slot_len],
        )?;
    }

    Ok(MphEcLayout {
        slot_offsets,
        group_offsets,
        group_displacements,
        suffixes,
        payloads,
    })
}

fn build_mph_bucket(
    k: u8,
    prefix_bits: u8,
    records: &[EcKmerRecord],
    ec_offsets: &[u32],
    displacements: &mut [u32],
    suffix_slots: &mut [u64],
    payload_slots: &mut [EcV5KmerPayload],
) -> Result<(), EcIndexError> {
    let slot_len = suffix_slots.len();
    let group_count = displacements.len();
    let mut groups = vec![Vec::<usize>::new(); group_count];
    let suffixes = records
        .iter()
        .map(|record| ec_kmer_suffix(record.kmer_code, k, prefix_bits))
        .collect::<Vec<_>>();
    for (idx, &suffix) in suffixes.iter().enumerate() {
        groups[mph_hash1(suffix) as usize % group_count].push(idx);
    }
    let mut order = (0..group_count).collect::<Vec<_>>();
    order.sort_unstable_by_key(|&group| std::cmp::Reverse(groups[group].len()));
    let mut occupied = vec![false; slot_len];

    for group in order {
        if groups[group].is_empty() {
            continue;
        }
        let mut displacement = 0_u32;
        let found = loop {
            if displacement == u32::MAX {
                break false;
            }
            let mut slots = Vec::with_capacity(groups[group].len());
            let mut ok = true;
            for &idx in &groups[group] {
                let slot = mph_slot(suffixes[idx], displacement as u64, slot_len);
                if occupied[slot] || slots.contains(&slot) {
                    ok = false;
                    break;
                }
                slots.push(slot);
            }
            if ok {
                for &idx in &groups[group] {
                    let slot = mph_slot(suffixes[idx], displacement as u64, slot_len);
                    occupied[slot] = true;
                    suffix_slots[slot] = suffixes[idx];
                    payload_slots[slot] = ec_payload_for_record(records[idx], ec_offsets)?;
                }
                break true;
            }
            displacement = displacement.wrapping_add(1);
        };
        if !found {
            return Err(EcIndexError::TooLarge("EC v8 displacement search"));
        }
        displacements[group] = displacement;
    }
    Ok(())
}

fn ec_payload_for_record(
    record: EcKmerRecord,
    ec_offsets: &[u32],
) -> Result<EcV5KmerPayload, EcIndexError> {
    let start = ec_offsets
        .get(record.ec_id as usize)
        .copied()
        .ok_or(EcIndexError::Truncated)?;
    let end = ec_offsets
        .get(record.ec_id as usize + 1)
        .copied()
        .ok_or(EcIndexError::Truncated)?;
    Ok(EcV5KmerPayload {
        target_start: start,
        gene_df: u16::try_from(record.gene_df)
            .map_err(|_| EcIndexError::TooLarge("EC v8 gene df"))?,
        target_len: u16::try_from(end.saturating_sub(start))
            .map_err(|_| EcIndexError::TooLarge("EC v8 target length"))?,
    })
}

fn mph_slot_count(len: usize) -> usize {
    let loaded = len
        .saturating_mul(EC_MPH_LOAD_NUMERATOR)
        .div_ceil(EC_MPH_LOAD_DENOMINATOR);
    loaded.max(len)
}

#[inline(always)]
fn mph_slot(suffix: u64, displacement: u64, slot_len: usize) -> usize {
    splitmix64(suffix ^ displacement.wrapping_mul(0x9e37_79b9_7f4a_7c15)) as usize % slot_len
}

fn write_v7_payloads_from_v5_payloads(
    writer: &mut impl Write,
    payloads: &[EcV5KmerPayload],
    ec_offsets: &[u32],
    ec_ranges: &[EcTargetRange],
) -> Result<(), EcIndexError> {
    let mut buffer = Vec::with_capacity(1 << 20);
    for payload in payloads {
        let ec_id = ec_offsets
            .binary_search(&payload.target_start)
            .map_err(|_| EcIndexError::Truncated)?;
        let range = ec_ranges.get(ec_id).ok_or(EcIndexError::Truncated)?;
        let rewritten = EcV5KmerPayload {
            target_start: range.target_start,
            gene_df: payload.gene_df,
            target_len: range.target_len,
        };
        buffer.extend_from_slice(bytemuck::bytes_of(&rewritten));
        if buffer.len() >= 1 << 20 {
            writer.write_all(&buffer)?;
            buffer.clear();
        }
    }
    if !buffer.is_empty() {
        writer.write_all(&buffer)?;
    }
    Ok(())
}

fn unpack_v5_targets(
    index: &MmapEcIndex,
    ec_targets_offset: u64,
    num_ec_targets: u64,
) -> Result<(Vec<GeneId>, Vec<u8>), EcIndexError> {
    let targets = index.bytes(
        ec_targets_offset,
        num_ec_targets * EC_V5_TARGET_WIDTH as u64,
    )?;
    let mut gene_ids = Vec::with_capacity(num_ec_targets as usize);
    let mut class_bits = Vec::with_capacity(num_ec_targets as usize);
    for target in targets.chunks_exact(EC_V5_TARGET_WIDTH) {
        gene_ids.push(u32::from_le_bytes([
            target[0], target[1], target[2], target[3],
        ]));
        class_bits.push(target[4]);
    }
    Ok((gene_ids, class_bits))
}

fn dedup_ec_targets(
    ec_offsets: &[u32],
    ec_gene_ids: &[GeneId],
    ec_class_bits: &[u8],
) -> Result<DedupedEcTargets, EcIndexError> {
    let forward = dedup_ec_targets_with_orientation(ec_offsets, ec_gene_ids, ec_class_bits, false)?;
    let reversed = dedup_ec_targets_with_orientation(ec_offsets, ec_gene_ids, ec_class_bits, true)?;
    if reversed.gene_ids.len() < forward.gene_ids.len() {
        Ok(reversed)
    } else {
        Ok(forward)
    }
}

fn ec_target_dedup_is_smaller(
    ec_offsets: &[u32],
    ec_gene_ids: &[GeneId],
    ec_class_bits: &[u8],
) -> Result<bool, EcIndexError> {
    let deduped = dedup_ec_targets(ec_offsets, ec_gene_ids, ec_class_bits)?;
    Ok(deduped_ec_target_storage_is_smaller(
        ec_offsets.len(),
        ec_gene_ids.len(),
        &deduped,
    )?)
}

fn deduped_ec_target_storage_is_smaller(
    legacy_ec_offsets_len: usize,
    legacy_ec_targets_len: usize,
    deduped: &DedupedEcTargets,
) -> Result<bool, EcIndexError> {
    let legacy_bytes = byte_len::<u32>(legacy_ec_offsets_len)
        + checked_mul_u64(
            legacy_ec_targets_len as u64,
            EC_V5_TARGET_WIDTH as u64,
            "EC target bytes",
        )?;
    let target_width = ec_target_width(&deduped.gene_ids);
    let deduped_bytes = byte_len::<EcTargetRange>(deduped.ec_ranges.len())
        + checked_mul_u64(
            deduped.gene_ids.len() as u64,
            target_width as u64,
            "deduped EC target bytes",
        )?;
    Ok(deduped_bytes < legacy_bytes)
}

fn dedup_ec_targets_with_orientation(
    ec_offsets: &[u32],
    ec_gene_ids: &[GeneId],
    ec_class_bits: &[u8],
    reverse_targets: bool,
) -> Result<DedupedEcTargets, EcIndexError> {
    if ec_offsets.is_empty() {
        return Ok(DedupedEcTargets {
            ec_ranges: Vec::new(),
            gene_ids: Vec::new(),
            class_bits: Vec::new(),
        });
    }
    let ec_count = ec_offsets.len() - 1;
    let mut ec_order = (0..ec_count).collect::<Vec<_>>();
    ec_order.sort_unstable_by_key(|&ec_id| {
        let start = ec_offsets[ec_id] as usize;
        let end = ec_offsets[ec_id + 1] as usize;
        std::cmp::Reverse(end.saturating_sub(start))
    });

    let mut ec_ranges = vec![EcTargetRange::default(); ec_count];
    let mut gene_ids = Vec::new();
    let mut class_bits = Vec::new();
    let mut suffixes: FxHashMap<(u64, u32), Vec<u32>> = FxHashMap::default();
    for ec_id in ec_order {
        let list = ec_target_list(
            ec_offsets,
            ec_gene_ids,
            ec_class_bits,
            ec_id,
            reverse_targets,
        )?;
        let target_len =
            u16::try_from(list.len()).map_err(|_| EcIndexError::TooLarge("EC v7 target length"))?;
        if list.is_empty() {
            ec_ranges[ec_id] = EcTargetRange {
                target_start: 0,
                target_len,
                _pad0: 0,
            };
            continue;
        }
        let key = (hash_ec_target_suffix(&list), list.len() as u32);
        if let Some(start) =
            find_deduped_target_suffix(&suffixes, key, &gene_ids, &class_bits, &list)
        {
            ec_ranges[ec_id] = EcTargetRange {
                target_start: start,
                target_len,
                _pad0: 0,
            };
            continue;
        }

        let target_start =
            u32::try_from(gene_ids.len()).map_err(|_| EcIndexError::TooLarge("EC v7 targets"))?;
        for target in &list {
            gene_ids.push(target.gene_id);
            class_bits.push(target.class_bits);
        }
        register_ec_target_suffixes(&mut suffixes, &list, target_start)?;
        ec_ranges[ec_id] = EcTargetRange {
            target_start,
            target_len,
            _pad0: 0,
        };
    }
    Ok(DedupedEcTargets {
        ec_ranges,
        gene_ids,
        class_bits,
    })
}

fn ec_target_list(
    ec_offsets: &[u32],
    ec_gene_ids: &[GeneId],
    ec_class_bits: &[u8],
    ec_id: usize,
    reverse_targets: bool,
) -> Result<Vec<EcGeneTarget>, EcIndexError> {
    let start = *ec_offsets.get(ec_id).ok_or(EcIndexError::Truncated)? as usize;
    let end = *ec_offsets.get(ec_id + 1).ok_or(EcIndexError::Truncated)? as usize;
    if end < start || end > ec_gene_ids.len() {
        return Err(EcIndexError::Truncated);
    }
    let mut targets = (start..end)
        .map(|idx| EcGeneTarget {
            gene_id: ec_gene_ids[idx],
            class_bits: ec_class_bits.get(idx).copied().unwrap_or(EC_CLASS_UNKNOWN),
        })
        .collect::<Vec<_>>();
    if reverse_targets {
        targets.reverse();
    }
    Ok(targets)
}

fn find_deduped_target_suffix(
    suffixes: &FxHashMap<(u64, u32), Vec<u32>>,
    key: (u64, u32),
    gene_ids: &[GeneId],
    class_bits: &[u8],
    list: &[EcGeneTarget],
) -> Option<u32> {
    let candidates = suffixes.get(&key)?;
    candidates.iter().copied().find(|&start| {
        let start = start as usize;
        let end = start + list.len();
        end <= gene_ids.len()
            && list.iter().enumerate().all(|(idx, target)| {
                gene_ids[start + idx] == target.gene_id
                    && class_bits
                        .get(start + idx)
                        .copied()
                        .unwrap_or(EC_CLASS_UNKNOWN)
                        == target.class_bits
            })
    })
}

fn register_ec_target_suffixes(
    suffixes: &mut FxHashMap<(u64, u32), Vec<u32>>,
    list: &[EcGeneTarget],
    target_start: u32,
) -> Result<(), EcIndexError> {
    let mut hash = EC_TARGET_HASH_OFFSET;
    for (idx, target) in list.iter().enumerate().rev() {
        hash = ec_target_hash_step(hash, *target);
        let suffix_len = u32::try_from(list.len() - idx)
            .map_err(|_| EcIndexError::TooLarge("EC v7 target suffix"))?;
        let suffix_start = target_start
            .checked_add(idx as u32)
            .ok_or(EcIndexError::TooLarge("EC v7 target suffix"))?;
        suffixes
            .entry((hash, suffix_len))
            .or_default()
            .push(suffix_start);
    }
    Ok(())
}

const EC_TARGET_HASH_OFFSET: u64 = 0xcbf2_9ce4_8422_2325;
const EC_TARGET_HASH_PRIME: u64 = 0x0000_0100_0000_01b3;

fn hash_ec_target_suffix(list: &[EcGeneTarget]) -> u64 {
    list.iter()
        .rev()
        .fold(EC_TARGET_HASH_OFFSET, |hash, &target| {
            ec_target_hash_step(hash, target)
        })
}

fn ec_target_hash_step(hash: u64, target: EcGeneTarget) -> u64 {
    let value = ((target.gene_id as u64) << 8) | target.class_bits as u64;
    (hash ^ value).wrapping_mul(EC_TARGET_HASH_PRIME)
}

fn build_v6_blocks(
    k: u8,
    prefix_bits: u8,
    suffix_width: usize,
    block_size: usize,
    kmers: &[EcKmerRecord],
    kmer_prefix_offsets: &[u32],
) -> Result<(Vec<u32>, Vec<u32>, Vec<u64>), EcIndexError> {
    debug_assert!(block_size > 0);
    let mut prefix_offsets = vec![0_u32; kmer_prefix_offsets.len()];
    let mut block_offsets = Vec::new();
    let mut block_max_suffixes = Vec::new();
    let bucket_count = kmer_prefix_offsets.len().saturating_sub(1);
    for prefix in 0..bucket_count {
        prefix_offsets[prefix] = u32::try_from(block_max_suffixes.len())
            .map_err(|_| EcIndexError::TooLarge("EC v6 prefix block offsets"))?;
        let start = kmer_prefix_offsets[prefix] as usize;
        let end = kmer_prefix_offsets[prefix + 1] as usize;
        let mut record_start = start;
        while record_start < end {
            let record_end = (record_start + block_size).min(end);
            block_offsets.push(
                u32::try_from(record_start)
                    .map_err(|_| EcIndexError::TooLarge("EC v6 block offsets"))?,
            );
            let max_suffix = ec_kmer_suffix(kmers[record_end - 1].kmer_code, k, prefix_bits);
            if suffix_width < 8 && max_suffix >= (1_u64 << (suffix_width * 8)) {
                return Err(EcIndexError::TooLarge("EC v6 suffix"));
            }
            block_max_suffixes.push(max_suffix);
            record_start = record_end;
        }
    }
    prefix_offsets[bucket_count] = u32::try_from(block_max_suffixes.len())
        .map_err(|_| EcIndexError::TooLarge("EC v6 prefix block offsets"))?;
    block_offsets.push(
        u32::try_from(kmers.len()).map_err(|_| EcIndexError::TooLarge("EC v6 block offsets"))?,
    );
    Ok((prefix_offsets, block_offsets, block_max_suffixes))
}

fn build_v6_blocks_from_prefix_offsets(
    kmer_prefix_offsets: &[u32],
    block_size: usize,
    num_kmers: usize,
    mut max_suffix_at: impl FnMut(usize) -> u64,
) -> Result<(Vec<u32>, Vec<u32>, Vec<u64>), EcIndexError> {
    debug_assert!(block_size > 0);
    let mut prefix_offsets = vec![0_u32; kmer_prefix_offsets.len()];
    let mut block_offsets = Vec::new();
    let mut block_max_suffixes = Vec::new();
    let bucket_count = kmer_prefix_offsets.len().saturating_sub(1);
    for prefix in 0..bucket_count {
        prefix_offsets[prefix] = u32::try_from(block_max_suffixes.len())
            .map_err(|_| EcIndexError::TooLarge("EC v6 prefix block offsets"))?;
        let start = kmer_prefix_offsets[prefix] as usize;
        let end = kmer_prefix_offsets[prefix + 1] as usize;
        let mut record_start = start;
        while record_start < end {
            let record_end = (record_start + block_size).min(end);
            block_offsets.push(
                u32::try_from(record_start)
                    .map_err(|_| EcIndexError::TooLarge("EC v6 block offsets"))?,
            );
            block_max_suffixes.push(max_suffix_at(record_end - 1));
            record_start = record_end;
        }
    }
    prefix_offsets[bucket_count] = u32::try_from(block_max_suffixes.len())
        .map_err(|_| EcIndexError::TooLarge("EC v6 prefix block offsets"))?;
    block_offsets
        .push(u32::try_from(num_kmers).map_err(|_| EcIndexError::TooLarge("EC v6 block offsets"))?);
    Ok((prefix_offsets, block_offsets, block_max_suffixes))
}

fn write_v6_block_max_suffixes(
    writer: &mut impl Write,
    block_max_suffixes: &[u64],
    suffix_width: usize,
) -> Result<(), EcIndexError> {
    if suffix_width == 0 {
        return Ok(());
    }
    let mut buffer = Vec::with_capacity(1 << 20);
    for &suffix in block_max_suffixes {
        push_packed_le(&mut buffer, suffix, suffix_width);
        if buffer.len() >= 1 << 20 {
            writer.write_all(&buffer)?;
            buffer.clear();
        }
    }
    if !buffer.is_empty() {
        writer.write_all(&buffer)?;
    }
    Ok(())
}

fn write_packed_suffix_values(
    writer: &mut impl Write,
    suffixes: &[u64],
    suffix_width: usize,
) -> Result<(), EcIndexError> {
    if suffix_width == 0 {
        return Ok(());
    }
    let mut buffer = Vec::with_capacity(1 << 20);
    for &suffix in suffixes {
        push_packed_le(&mut buffer, suffix, suffix_width);
        if buffer.len() >= 1 << 20 {
            writer.write_all(&buffer)?;
            buffer.clear();
        }
    }
    if !buffer.is_empty() {
        writer.write_all(&buffer)?;
    }
    Ok(())
}

fn write_v6_suffix_blocks_from_row_major(
    writer: &mut impl Write,
    suffixes: &[u8],
    suffix_width: usize,
    block_offsets: &[u32],
) -> Result<(), EcIndexError> {
    if suffix_width == 0 {
        return Ok(());
    }
    let mut buffer = Vec::with_capacity(1 << 20);
    for block in block_offsets.windows(2) {
        let start = block[0] as usize;
        let end = block[1] as usize;
        for byte_idx in 0..suffix_width {
            for idx in start..end {
                buffer.push(suffixes[idx * suffix_width + byte_idx]);
                if buffer.len() >= 1 << 20 {
                    writer.write_all(&buffer)?;
                    buffer.clear();
                }
            }
        }
    }
    if !buffer.is_empty() {
        writer.write_all(&buffer)?;
    }
    Ok(())
}

fn write_v6_suffix_blocks(
    writer: &mut impl Write,
    kmers: &[EcKmerRecord],
    k: u8,
    prefix_bits: u8,
    suffix_width: usize,
    block_offsets: &[u32],
) -> Result<(), EcIndexError> {
    if suffix_width == 0 {
        return Ok(());
    }
    let mut buffer = Vec::with_capacity(1 << 20);
    for block in block_offsets.windows(2) {
        let start = block[0] as usize;
        let end = block[1] as usize;
        let len = end.saturating_sub(start);
        for byte_idx in 0..suffix_width {
            for record in &kmers[start..end] {
                let suffix = ec_kmer_suffix(record.kmer_code, k, prefix_bits);
                buffer.push(((suffix >> (byte_idx * 8)) & 0xff) as u8);
                if buffer.len() >= 1 << 20 {
                    writer.write_all(&buffer)?;
                    buffer.clear();
                }
            }
            debug_assert_eq!(len, end - start);
        }
    }
    if !buffer.is_empty() {
        writer.write_all(&buffer)?;
    }
    Ok(())
}

fn save_mmap_v6_from_mmap_v5(
    path: impl AsRef<Path>,
    index: &MmapEcIndex,
) -> Result<(), EcIndexError> {
    let v5 = index.header_v5.expect("v5 header");
    let prefix_offsets = index.prefix_offsets();
    let gene_metas = index.genes();
    let gene_names = index.gene_names();
    let ec_offsets = index.ec_offsets();
    let suffix_width = v5.suffix_width as usize;
    let block_size = EC_V6_BLOCK_SIZE;
    let (v6_prefix_offsets, block_offsets, block_max_suffixes) =
        build_v6_blocks_from_prefix_offsets(
            prefix_offsets,
            block_size,
            v5.num_kmers as usize,
            |idx| {
                let suffixes = index.v5_suffixes();
                read_packed_le(&suffixes[idx * suffix_width..idx * suffix_width + suffix_width])
            },
        )?;
    let num_blocks = block_max_suffixes.len() as u64;

    let mut cursor = align_up(std::mem::size_of::<EcMmapHeaderV6>() as u64, 8);
    let prefix_offsets_offset = cursor;
    cursor += byte_len::<u32>(v6_prefix_offsets.len());
    cursor = align_up(cursor, 8);
    let genes_offset = cursor;
    cursor += byte_len::<EcGeneMeta>(gene_metas.len());
    cursor = align_up(cursor, 4);
    let block_offsets_offset = cursor;
    cursor += byte_len::<u32>(block_offsets.len());
    cursor = align_up(cursor, 8);
    let block_max_suffixes_offset = cursor;
    cursor += checked_mul_u64(num_blocks, suffix_width as u64, "EC v6 block max suffixes")?;
    cursor = align_up(cursor, 8);
    let suffixes_offset = cursor;
    cursor += checked_mul_u64(v5.num_kmers, suffix_width as u64, "EC v6 suffixes")?;
    cursor = align_up(cursor, 8);
    let payloads_offset = cursor;
    cursor += byte_len::<EcV5KmerPayload>(v5.num_kmers as usize);
    cursor = align_up(cursor, 4);
    let ec_offsets_offset = cursor;
    cursor += byte_len::<u32>(ec_offsets.len());
    cursor = align_up(cursor, 4);
    let ec_targets_offset = cursor;
    cursor += checked_mul_u64(
        v5.num_ec_targets,
        EC_V5_TARGET_WIDTH as u64,
        "EC v6 targets",
    )?;
    cursor = align_up(cursor, 8);
    let gene_names_offset = cursor;

    let header = EcMmapHeaderV6 {
        magic: EC_MMAP_MAGIC,
        version: EC_MMAP_VERSION_V6,
        k: index.k() as u32,
        prefix_bits: v5.prefix_bits,
        suffix_bits: v5.suffix_bits,
        suffix_width: v5.suffix_width,
        block_size: block_size as u32,
        num_genes: gene_metas.len() as u64,
        num_kmers: v5.num_kmers,
        num_blocks,
        num_ec_offsets: ec_offsets.len() as u64,
        num_ec_targets: v5.num_ec_targets,
        genes_offset,
        prefix_offsets_offset,
        prefix_offsets_len: byte_len::<u32>(v6_prefix_offsets.len()),
        block_offsets_offset,
        block_max_suffixes_offset,
        suffixes_offset,
        payloads_offset,
        ec_offsets_offset,
        ec_targets_offset,
        gene_names_offset,
        gene_names_len: gene_names.len() as u64,
    };

    let mut file = File::create(path)?;
    file.write_all(bytemuck::bytes_of(&header))?;
    write_padding(
        &mut file,
        std::mem::size_of::<EcMmapHeaderV6>() as u64,
        prefix_offsets_offset,
    )?;
    write_pod_slice(&mut file, &v6_prefix_offsets)?;
    write_padding(
        &mut file,
        prefix_offsets_offset + byte_len::<u32>(v6_prefix_offsets.len()),
        genes_offset,
    )?;
    write_pod_slice(&mut file, gene_metas)?;
    write_padding(
        &mut file,
        genes_offset + byte_len::<EcGeneMeta>(gene_metas.len()),
        block_offsets_offset,
    )?;
    write_pod_slice(&mut file, &block_offsets)?;
    write_padding(
        &mut file,
        block_offsets_offset + byte_len::<u32>(block_offsets.len()),
        block_max_suffixes_offset,
    )?;
    write_v6_block_max_suffixes(&mut file, &block_max_suffixes, suffix_width)?;
    write_padding(
        &mut file,
        block_max_suffixes_offset
            + checked_mul_u64(num_blocks, suffix_width as u64, "EC v6 block max suffixes")?,
        suffixes_offset,
    )?;
    write_v6_suffix_blocks_from_row_major(
        &mut file,
        index.v5_suffixes(),
        suffix_width,
        &block_offsets,
    )?;
    write_padding(
        &mut file,
        suffixes_offset + checked_mul_u64(v5.num_kmers, suffix_width as u64, "EC v6 suffixes")?,
        payloads_offset,
    )?;
    write_pod_slice(&mut file, index.v5_payloads())?;
    write_padding(
        &mut file,
        payloads_offset + byte_len::<EcV5KmerPayload>(v5.num_kmers as usize),
        ec_offsets_offset,
    )?;
    write_pod_slice(&mut file, ec_offsets)?;
    write_padding(
        &mut file,
        ec_offsets_offset + byte_len::<u32>(ec_offsets.len()),
        ec_targets_offset,
    )?;
    let targets = index.bytes(
        v5.ec_targets_offset,
        v5.num_ec_targets * EC_V5_TARGET_WIDTH as u64,
    )?;
    file.write_all(targets)?;
    write_padding(
        &mut file,
        ec_targets_offset
            + checked_mul_u64(
                v5.num_ec_targets,
                EC_V5_TARGET_WIDTH as u64,
                "EC v6 targets",
            )?,
        gene_names_offset,
    )?;
    file.write_all(gene_names)?;
    Ok(())
}

fn save_mmap_v7_from_mmap_v5(
    path: impl AsRef<Path>,
    index: &MmapEcIndex,
) -> Result<(), EcIndexError> {
    let v5 = index.header_v5.expect("v5 header");
    let prefix_offsets = index.prefix_offsets();
    let gene_metas = index.genes();
    let gene_names = index.gene_names();
    let ec_offsets = index.ec_offsets();
    let suffix_width = v5.suffix_width as usize;
    let block_size = EC_V6_BLOCK_SIZE;
    let (v7_prefix_offsets, block_offsets, block_max_suffixes) =
        build_v6_blocks_from_prefix_offsets(
            prefix_offsets,
            block_size,
            v5.num_kmers as usize,
            |idx| {
                let suffixes = index.v5_suffixes();
                read_packed_le(&suffixes[idx * suffix_width..idx * suffix_width + suffix_width])
            },
        )?;
    let (ec_gene_ids, ec_class_bits) =
        unpack_v5_targets(index, v5.ec_targets_offset, v5.num_ec_targets)?;
    let deduped_targets = dedup_ec_targets(ec_offsets, &ec_gene_ids, &ec_class_bits)?;
    if !deduped_ec_target_storage_is_smaller(
        ec_offsets.len(),
        v5.num_ec_targets as usize,
        &deduped_targets,
    )? {
        return save_mmap_v6_from_mmap_v5(path, index);
    }
    let num_blocks = block_max_suffixes.len() as u64;
    let num_targets = deduped_targets.gene_ids.len() as u64;
    let target_width = ec_target_width(&deduped_targets.gene_ids);

    let mut cursor = align_up(std::mem::size_of::<EcMmapHeaderV7>() as u64, 8);
    let prefix_offsets_offset = cursor;
    cursor += byte_len::<u32>(v7_prefix_offsets.len());
    cursor = align_up(cursor, 8);
    let genes_offset = cursor;
    cursor += byte_len::<EcGeneMeta>(gene_metas.len());
    cursor = align_up(cursor, 4);
    let block_offsets_offset = cursor;
    cursor += byte_len::<u32>(block_offsets.len());
    cursor = align_up(cursor, 8);
    let block_max_suffixes_offset = cursor;
    cursor += checked_mul_u64(num_blocks, suffix_width as u64, "EC v7 block max suffixes")?;
    cursor = align_up(cursor, 8);
    let suffixes_offset = cursor;
    cursor += checked_mul_u64(v5.num_kmers, suffix_width as u64, "EC v7 suffixes")?;
    cursor = align_up(cursor, 8);
    let payloads_offset = cursor;
    cursor += byte_len::<EcV5KmerPayload>(v5.num_kmers as usize);
    cursor = align_up(cursor, 8);
    let ec_ranges_offset = cursor;
    cursor += byte_len::<EcTargetRange>(deduped_targets.ec_ranges.len());
    cursor = align_up(cursor, 4);
    let ec_targets_offset = cursor;
    cursor += checked_mul_u64(num_targets, target_width as u64, "EC v7 targets")?;
    cursor = align_up(cursor, 8);
    let gene_names_offset = cursor;

    let header = EcMmapHeaderV7 {
        magic: EC_MMAP_MAGIC,
        version: EC_MMAP_VERSION,
        k: index.k() as u32,
        prefix_bits: v5.prefix_bits,
        suffix_bits: v5.suffix_bits,
        suffix_width: v5.suffix_width,
        block_size: block_size as u32,
        target_width: target_width as u32,
        _pad0: 0,
        num_genes: gene_metas.len() as u64,
        num_kmers: v5.num_kmers,
        num_blocks,
        num_ec_ranges: deduped_targets.ec_ranges.len() as u64,
        num_ec_targets: num_targets,
        genes_offset,
        prefix_offsets_offset,
        prefix_offsets_len: byte_len::<u32>(v7_prefix_offsets.len()),
        block_offsets_offset,
        block_max_suffixes_offset,
        suffixes_offset,
        payloads_offset,
        ec_ranges_offset,
        ec_targets_offset,
        gene_names_offset,
        gene_names_len: gene_names.len() as u64,
    };

    let mut file = File::create(path)?;
    file.write_all(bytemuck::bytes_of(&header))?;
    write_padding(
        &mut file,
        std::mem::size_of::<EcMmapHeaderV7>() as u64,
        prefix_offsets_offset,
    )?;
    write_pod_slice(&mut file, &v7_prefix_offsets)?;
    write_padding(
        &mut file,
        prefix_offsets_offset + byte_len::<u32>(v7_prefix_offsets.len()),
        genes_offset,
    )?;
    write_pod_slice(&mut file, gene_metas)?;
    write_padding(
        &mut file,
        genes_offset + byte_len::<EcGeneMeta>(gene_metas.len()),
        block_offsets_offset,
    )?;
    write_pod_slice(&mut file, &block_offsets)?;
    write_padding(
        &mut file,
        block_offsets_offset + byte_len::<u32>(block_offsets.len()),
        block_max_suffixes_offset,
    )?;
    write_v6_block_max_suffixes(&mut file, &block_max_suffixes, suffix_width)?;
    write_padding(
        &mut file,
        block_max_suffixes_offset
            + checked_mul_u64(num_blocks, suffix_width as u64, "EC v7 block max suffixes")?,
        suffixes_offset,
    )?;
    write_v6_suffix_blocks_from_row_major(
        &mut file,
        index.v5_suffixes(),
        suffix_width,
        &block_offsets,
    )?;
    write_padding(
        &mut file,
        suffixes_offset + checked_mul_u64(v5.num_kmers, suffix_width as u64, "EC v7 suffixes")?,
        payloads_offset,
    )?;
    write_v7_payloads_from_v5_payloads(
        &mut file,
        index.v5_payloads(),
        ec_offsets,
        &deduped_targets.ec_ranges,
    )?;
    write_padding(
        &mut file,
        payloads_offset + byte_len::<EcV5KmerPayload>(v5.num_kmers as usize),
        ec_ranges_offset,
    )?;
    write_pod_slice(&mut file, &deduped_targets.ec_ranges)?;
    write_padding(
        &mut file,
        ec_ranges_offset + byte_len::<EcTargetRange>(deduped_targets.ec_ranges.len()),
        ec_targets_offset,
    )?;
    write_v7_targets(
        &mut file,
        &deduped_targets.gene_ids,
        &deduped_targets.class_bits,
        target_width,
    )?;
    write_padding(
        &mut file,
        ec_targets_offset + checked_mul_u64(num_targets, target_width as u64, "EC v7 targets")?,
        gene_names_offset,
    )?;
    file.write_all(gene_names)?;
    Ok(())
}

fn save_mmap_v5_from_mmap_v4(
    path: impl AsRef<Path>,
    index: &MmapEcIndex,
) -> Result<(), EcIndexError> {
    if !can_save_mmap_v5_from_mmap_v4(index) {
        let mut file = File::create(path)?;
        file.write_all(&index.mmap)?;
        return Ok(());
    }

    let v4 = index.header_v4.expect("v4 header");
    let prefix_offsets = index.prefix_offsets();
    let gene_metas = index.genes();
    let gene_names = index.gene_names();
    let ec_offsets = index.ec_offsets();
    let num_kmers = v4.num_kmers;
    let num_targets = v4.num_ec_targets;
    let suffix_width = v4.suffix_width as usize;

    let mut cursor = align_up(std::mem::size_of::<EcMmapHeaderV5>() as u64, 8);
    let prefix_offsets_offset = cursor;
    cursor += byte_len::<u32>(prefix_offsets.len());
    cursor = align_up(cursor, 8);
    let genes_offset = cursor;
    cursor += byte_len::<EcGeneMeta>(gene_metas.len());
    cursor = align_up(cursor, 8);
    let suffixes_offset = cursor;
    cursor += checked_mul_u64(num_kmers, suffix_width as u64, "EC v5 suffixes")?;
    cursor = align_up(cursor, 8);
    let payloads_offset = cursor;
    cursor += checked_mul_u64(
        num_kmers,
        std::mem::size_of::<EcV5KmerPayload>() as u64,
        "EC v5 payloads",
    )?;
    cursor = align_up(cursor, 4);
    let ec_offsets_offset = cursor;
    cursor += byte_len::<u32>(ec_offsets.len());
    cursor = align_up(cursor, 4);
    let ec_targets_offset = cursor;
    cursor += checked_mul_u64(num_targets, EC_V5_TARGET_WIDTH as u64, "EC v5 targets")?;
    cursor = align_up(cursor, 8);
    let gene_names_offset = cursor;

    let header = EcMmapHeaderV5 {
        magic: EC_MMAP_MAGIC,
        version: EC_MMAP_VERSION_V5,
        k: index.k() as u32,
        prefix_bits: v4.prefix_bits,
        suffix_bits: v4.suffix_bits,
        suffix_width: v4.suffix_width,
        _pad1: 0,
        num_genes: gene_metas.len() as u64,
        num_kmers,
        num_ec_offsets: ec_offsets.len() as u64,
        num_ec_targets: num_targets,
        genes_offset,
        suffixes_offset,
        payloads_offset,
        ec_offsets_offset,
        ec_targets_offset,
        gene_names_offset,
        gene_names_len: gene_names.len() as u64,
        prefix_offsets_offset,
        prefix_offsets_len: byte_len::<u32>(prefix_offsets.len()),
    };

    let mut file = File::create(path)?;
    file.write_all(bytemuck::bytes_of(&header))?;
    write_padding(
        &mut file,
        std::mem::size_of::<EcMmapHeaderV5>() as u64,
        prefix_offsets_offset,
    )?;
    write_pod_slice(&mut file, prefix_offsets)?;
    write_padding(
        &mut file,
        prefix_offsets_offset + byte_len::<u32>(prefix_offsets.len()),
        genes_offset,
    )?;
    write_pod_slice(&mut file, gene_metas)?;
    write_padding(
        &mut file,
        genes_offset + byte_len::<EcGeneMeta>(gene_metas.len()),
        suffixes_offset,
    )?;
    let suffixes = index.bytes(v4.suffixes_offset, v4.num_kmers * v4.suffix_width as u64)?;
    file.write_all(suffixes)?;
    write_padding(
        &mut file,
        suffixes_offset + checked_mul_u64(num_kmers, suffix_width as u64, "EC v5 suffixes")?,
        payloads_offset,
    )?;
    write_v5_payloads_from_mmap_v4(&mut file, index)?;
    write_padding(
        &mut file,
        payloads_offset
            + checked_mul_u64(
                num_kmers,
                std::mem::size_of::<EcV5KmerPayload>() as u64,
                "EC v5 payloads",
            )?,
        ec_offsets_offset,
    )?;
    write_pod_slice(&mut file, ec_offsets)?;
    write_padding(
        &mut file,
        ec_offsets_offset + byte_len::<u32>(ec_offsets.len()),
        ec_targets_offset,
    )?;
    let targets = index.bytes(
        v4.ec_targets_offset,
        v4.num_ec_targets * EC_V4_TARGET_WIDTH as u64,
    )?;
    file.write_all(targets)?;
    write_padding(
        &mut file,
        ec_targets_offset
            + checked_mul_u64(num_targets, EC_V5_TARGET_WIDTH as u64, "EC v5 targets")?,
        gene_names_offset,
    )?;
    file.write_all(gene_names)?;
    Ok(())
}

fn can_save_mmap_v5_from_mmap_v4(index: &MmapEcIndex) -> bool {
    let Some(v4) = index.header_v4 else {
        return false;
    };
    if v4.suffix_bits > 48 || v4.gene_df_width > 2 {
        return false;
    }
    let ec_offsets = index.ec_offsets();
    for &ec_id in index.v4_ec_ids() {
        let Some((&start, &end)) = ec_offsets
            .get(ec_id as usize)
            .zip(ec_offsets.get(ec_id as usize + 1))
        else {
            return false;
        };
        if end < start || end - start > u16::MAX as u32 {
            return false;
        }
    }
    true
}

fn write_v5_payloads_from_mmap_v4(
    writer: &mut impl Write,
    index: &MmapEcIndex,
) -> Result<(), EcIndexError> {
    let ec_ids = index.v4_ec_ids();
    let ec_offsets = index.ec_offsets();
    let mut buffer = Vec::with_capacity(1 << 20);
    for (idx, &ec_id) in ec_ids.iter().enumerate() {
        let start = ec_offsets
            .get(ec_id as usize)
            .copied()
            .ok_or(EcIndexError::Truncated)?;
        let end = ec_offsets
            .get(ec_id as usize + 1)
            .copied()
            .ok_or(EcIndexError::Truncated)?;
        let payload = EcV5KmerPayload {
            target_start: start,
            gene_df: u16::try_from(index.v4_gene_df(idx))
                .map_err(|_| EcIndexError::TooLarge("EC v5 gene df"))?,
            target_len: u16::try_from(end.saturating_sub(start))
                .map_err(|_| EcIndexError::TooLarge("EC v5 target length"))?,
        };
        buffer.extend_from_slice(bytemuck::bytes_of(&payload));
        if buffer.len() >= 1 << 20 {
            writer.write_all(&buffer)?;
            buffer.clear();
        }
    }
    if !buffer.is_empty() {
        writer.write_all(&buffer)?;
    }
    Ok(())
}

fn save_mmap_v4_from_slices(
    path: impl AsRef<Path>,
    k: u8,
    gene_metas: &[EcGeneMeta],
    gene_names: &[u8],
    kmers: &[EcKmerRecord],
    ec_offsets: &[u32],
    ec_gene_ids: &[GeneId],
    ec_class_bits: &[u8],
    prefix_offsets: &[u32],
) -> Result<(), EcIndexError> {
    let prefix_bits = EC_PREFIX_BITS.min(k.saturating_mul(2));
    let suffix_bits = k.saturating_mul(2).saturating_sub(prefix_bits);
    let suffix_width = suffix_bits.div_ceil(8) as usize;
    let gene_df_width = if kmers.iter().all(|record| record.gene_df <= u16::MAX as u32) {
        2_usize
    } else {
        4_usize
    };
    let num_kmers = kmers.len() as u64;
    let num_targets = ec_gene_ids.len() as u64;

    let mut cursor = align_up(std::mem::size_of::<EcMmapHeaderV4>() as u64, 8);
    let prefix_offsets_offset = cursor;
    cursor += byte_len::<u32>(prefix_offsets.len());
    cursor = align_up(cursor, 8);
    let genes_offset = cursor;
    cursor += byte_len::<EcGeneMeta>(gene_metas.len());
    cursor = align_up(cursor, 8);
    let suffixes_offset = cursor;
    cursor += checked_mul_u64(num_kmers, suffix_width as u64, "EC v4 suffixes")?;
    cursor = align_up(cursor, 4);
    let ec_ids_offset = cursor;
    cursor += byte_len::<u32>(kmers.len());
    cursor = align_up(cursor, 4);
    let gene_dfs_offset = cursor;
    cursor += checked_mul_u64(num_kmers, gene_df_width as u64, "EC v4 gene dfs")?;
    cursor = align_up(cursor, 4);
    let ec_offsets_offset = cursor;
    cursor += byte_len::<u32>(ec_offsets.len());
    cursor = align_up(cursor, 4);
    let ec_targets_offset = cursor;
    cursor += checked_mul_u64(num_targets, EC_V4_TARGET_WIDTH as u64, "EC v4 targets")?;
    cursor = align_up(cursor, 8);
    let gene_names_offset = cursor;

    let header = EcMmapHeaderV4 {
        magic: EC_MMAP_MAGIC,
        version: EC_MMAP_VERSION_V4,
        k: k as u32,
        prefix_bits: prefix_bits as u32,
        suffix_bits: suffix_bits as u32,
        suffix_width: suffix_width as u32,
        gene_df_width: gene_df_width as u32,
        _pad0: 0,
        _pad1: 0,
        num_genes: gene_metas.len() as u64,
        num_kmers,
        num_ec_offsets: ec_offsets.len() as u64,
        num_ec_targets: num_targets,
        genes_offset,
        suffixes_offset,
        ec_ids_offset,
        gene_dfs_offset,
        ec_offsets_offset,
        ec_targets_offset,
        gene_names_offset,
        gene_names_len: gene_names.len() as u64,
        prefix_offsets_offset,
        prefix_offsets_len: byte_len::<u32>(prefix_offsets.len()),
    };

    let mut file = File::create(path)?;
    file.write_all(bytemuck::bytes_of(&header))?;
    write_padding(
        &mut file,
        std::mem::size_of::<EcMmapHeaderV4>() as u64,
        prefix_offsets_offset,
    )?;
    write_pod_slice(&mut file, prefix_offsets)?;
    write_padding(
        &mut file,
        prefix_offsets_offset + byte_len::<u32>(prefix_offsets.len()),
        genes_offset,
    )?;
    write_pod_slice(&mut file, gene_metas)?;
    write_padding(
        &mut file,
        genes_offset + byte_len::<EcGeneMeta>(gene_metas.len()),
        suffixes_offset,
    )?;
    write_v4_suffixes(&mut file, kmers, k, prefix_bits, suffix_width)?;
    write_padding(
        &mut file,
        suffixes_offset + checked_mul_u64(num_kmers, suffix_width as u64, "EC v4 suffixes")?,
        ec_ids_offset,
    )?;
    write_v4_ec_ids(&mut file, kmers)?;
    write_padding(
        &mut file,
        ec_ids_offset + byte_len::<u32>(kmers.len()),
        gene_dfs_offset,
    )?;
    write_v4_gene_dfs(&mut file, kmers, gene_df_width)?;
    write_padding(
        &mut file,
        gene_dfs_offset + checked_mul_u64(num_kmers, gene_df_width as u64, "EC v4 gene dfs")?,
        ec_offsets_offset,
    )?;
    write_pod_slice(&mut file, ec_offsets)?;
    write_padding(
        &mut file,
        ec_offsets_offset + byte_len::<u32>(ec_offsets.len()),
        ec_targets_offset,
    )?;
    write_v4_targets(&mut file, ec_gene_ids, ec_class_bits)?;
    write_padding(
        &mut file,
        ec_targets_offset
            + checked_mul_u64(num_targets, EC_V4_TARGET_WIDTH as u64, "EC v4 targets")?,
        gene_names_offset,
    )?;
    file.write_all(gene_names)?;
    Ok(())
}

fn write_v4_suffixes(
    writer: &mut impl Write,
    kmers: &[EcKmerRecord],
    k: u8,
    prefix_bits: u8,
    suffix_width: usize,
) -> Result<(), EcIndexError> {
    if suffix_width == 0 {
        return Ok(());
    }
    let mut buffer = Vec::with_capacity(1 << 20);
    for record in kmers {
        let suffix = ec_kmer_suffix(record.kmer_code, k, prefix_bits);
        push_packed_le(&mut buffer, suffix, suffix_width);
        if buffer.len() >= 1 << 20 {
            writer.write_all(&buffer)?;
            buffer.clear();
        }
    }
    if !buffer.is_empty() {
        writer.write_all(&buffer)?;
    }
    Ok(())
}

fn write_v4_ec_ids(writer: &mut impl Write, kmers: &[EcKmerRecord]) -> Result<(), EcIndexError> {
    let mut buffer = Vec::with_capacity(1 << 20);
    for record in kmers {
        buffer.extend_from_slice(&record.ec_id.to_le_bytes());
        if buffer.len() >= 1 << 20 {
            writer.write_all(&buffer)?;
            buffer.clear();
        }
    }
    if !buffer.is_empty() {
        writer.write_all(&buffer)?;
    }
    Ok(())
}

fn write_v4_gene_dfs(
    writer: &mut impl Write,
    kmers: &[EcKmerRecord],
    gene_df_width: usize,
) -> Result<(), EcIndexError> {
    let mut buffer = Vec::with_capacity(1 << 20);
    for record in kmers {
        match gene_df_width {
            2 => buffer.extend_from_slice(&(record.gene_df as u16).to_le_bytes()),
            4 => buffer.extend_from_slice(&record.gene_df.to_le_bytes()),
            _ => unreachable!("validated EC v4 gene df width"),
        }
        if buffer.len() >= 1 << 20 {
            writer.write_all(&buffer)?;
            buffer.clear();
        }
    }
    if !buffer.is_empty() {
        writer.write_all(&buffer)?;
    }
    Ok(())
}

fn write_v4_targets(
    writer: &mut impl Write,
    ec_gene_ids: &[GeneId],
    ec_class_bits: &[u8],
) -> Result<(), EcIndexError> {
    let mut buffer = Vec::with_capacity(1 << 20);
    for (idx, &gene_id) in ec_gene_ids.iter().enumerate() {
        buffer.extend_from_slice(&gene_id.to_le_bytes());
        buffer.push(ec_class_bits.get(idx).copied().unwrap_or(EC_CLASS_UNKNOWN));
        if buffer.len() >= 1 << 20 {
            writer.write_all(&buffer)?;
            buffer.clear();
        }
    }
    if !buffer.is_empty() {
        writer.write_all(&buffer)?;
    }
    Ok(())
}

fn write_v7_targets(
    writer: &mut impl Write,
    ec_gene_ids: &[GeneId],
    ec_class_bits: &[u8],
    target_width: usize,
) -> Result<(), EcIndexError> {
    if target_width == EC_V5_TARGET_WIDTH {
        return write_v4_targets(writer, ec_gene_ids, ec_class_bits);
    }
    debug_assert_eq!(target_width, 4);
    let mut buffer = Vec::with_capacity(1 << 20);
    for (idx, &gene_id) in ec_gene_ids.iter().enumerate() {
        if gene_id > 0x00ff_ffff {
            return Err(EcIndexError::TooLarge("EC v7 packed gene id"));
        }
        let bytes = gene_id.to_le_bytes();
        buffer.extend_from_slice(&bytes[..3]);
        buffer.push(ec_class_bits.get(idx).copied().unwrap_or(EC_CLASS_UNKNOWN));
        if buffer.len() >= 1 << 20 {
            writer.write_all(&buffer)?;
            buffer.clear();
        }
    }
    if !buffer.is_empty() {
        writer.write_all(&buffer)?;
    }
    Ok(())
}

fn ec_target_width(ec_gene_ids: &[GeneId]) -> usize {
    if ec_gene_ids.iter().all(|&gene_id| gene_id <= 0x00ff_ffff) {
        4
    } else {
        EC_V5_TARGET_WIDTH
    }
}

#[inline(always)]
fn mph_hash1(value: u64) -> u64 {
    splitmix64(value ^ 0x9e37_79b9_7f4a_7c15)
}

#[inline(always)]
fn splitmix64(mut value: u64) -> u64 {
    value = value.wrapping_add(0x9e37_79b9_7f4a_7c15);
    value = (value ^ (value >> 30)).wrapping_mul(0xbf58_476d_1ce4_e5b9);
    value = (value ^ (value >> 27)).wrapping_mul(0x94d0_49bb_1331_11eb);
    value ^ (value >> 31)
}

#[inline(always)]
fn ec_kmer_prefix(kmer_code: u64, k: u8, prefix_bits: u8) -> u32 {
    if prefix_bits == 0 {
        return 0;
    }
    let total_bits = k.saturating_mul(2);
    let shift = total_bits.saturating_sub(prefix_bits) as u32;
    let mask = if prefix_bits >= 32 {
        u32::MAX as u64
    } else {
        (1_u64 << prefix_bits) - 1
    };
    ((kmer_code >> shift) & mask) as u32
}

#[inline(always)]
fn ec_kmer_suffix(kmer_code: u64, k: u8, prefix_bits: u8) -> u64 {
    let total_bits = k.saturating_mul(2);
    let suffix_bits = total_bits.saturating_sub(prefix_bits);
    if suffix_bits == 0 {
        0
    } else if suffix_bits >= 64 {
        kmer_code
    } else {
        kmer_code & ((1_u64 << suffix_bits) - 1)
    }
}

#[inline(always)]
fn push_packed_le(out: &mut Vec<u8>, mut value: u64, width: usize) {
    for _ in 0..width {
        out.push((value & 0xff) as u8);
        value >>= 8;
    }
}

#[inline(always)]
fn read_packed_le(bytes: &[u8]) -> u64 {
    bytes.iter().enumerate().fold(0_u64, |value, (idx, &byte)| {
        value | ((byte as u64) << (idx * 8))
    })
}

#[inline(always)]
fn read_v6_planar_suffix(
    suffixes: &[u8],
    block_base: usize,
    block_len: usize,
    lane: usize,
    suffix_width: usize,
) -> u64 {
    let mut value = 0_u64;
    for byte_idx in 0..suffix_width {
        value |= (suffixes[block_base + byte_idx * block_len + lane] as u64) << (byte_idx * 8);
    }
    value
}

#[inline(always)]
fn zero_byte_high_bits(value: u64) -> u64 {
    !value & value.wrapping_sub(0x0101_0101_0101_0101) & 0x8080_8080_8080_8080
}

fn join_nul_strings(values: &[String]) -> Vec<u8> {
    let mut out = Vec::new();
    for value in values {
        out.extend(value.as_bytes());
        out.push(0);
    }
    out
}

fn split_nul_strings(bytes: &[u8], expected: usize) -> Result<Vec<String>, EcIndexError> {
    let values = bytes
        .split(|&byte| byte == 0)
        .filter(|part| !part.is_empty())
        .map(|part| String::from_utf8_lossy(part).into_owned())
        .collect::<Vec<_>>();
    if values.len() != expected {
        return Err(EcIndexError::Truncated);
    }
    Ok(values)
}

fn read_u32(reader: &mut impl Read) -> Result<u32, EcIndexError> {
    let mut bytes = [0_u8; 4];
    reader.read_exact(&mut bytes)?;
    Ok(u32::from_le_bytes(bytes))
}

fn read_u64(reader: &mut impl Read) -> Result<u64, EcIndexError> {
    let mut bytes = [0_u8; 8];
    reader.read_exact(&mut bytes)?;
    Ok(u64::from_le_bytes(bytes))
}

fn write_u32(writer: &mut impl Write, value: u32) -> Result<(), EcIndexError> {
    writer.write_all(&value.to_le_bytes())?;
    Ok(())
}

fn write_u64(writer: &mut impl Write, value: u64) -> Result<(), EcIndexError> {
    writer.write_all(&value.to_le_bytes())?;
    Ok(())
}

fn align_up(value: u64, alignment: u64) -> u64 {
    debug_assert!(alignment.is_power_of_two());
    (value + alignment - 1) & !(alignment - 1)
}

fn byte_len<T>(count: usize) -> u64 {
    (count * std::mem::size_of::<T>()) as u64
}

fn checked_mul_u64(a: u64, b: u64, field: &'static str) -> Result<u64, EcIndexError> {
    a.checked_mul(b).ok_or(EcIndexError::TooLarge(field))
}

fn write_padding(writer: &mut impl Write, current: u64, target: u64) -> Result<(), EcIndexError> {
    if current > target {
        return Err(EcIndexError::TooLarge("section offset"));
    }
    const ZEROES: [u8; 64] = [0; 64];
    let mut remaining = target - current;
    while remaining > 0 {
        let chunk = remaining.min(ZEROES.len() as u64) as usize;
        writer.write_all(&ZEROES[..chunk])?;
        remaining -= chunk as u64;
    }
    Ok(())
}

fn write_pod_slice<T: Pod>(writer: &mut impl Write, values: &[T]) -> Result<(), EcIndexError> {
    writer.write_all(bytemuck::cast_slice(values))?;
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Write;

    #[test]
    fn ec_index_maps_unique_gene_support() {
        let mut fasta = tempfile::NamedTempFile::new().unwrap();
        writeln!(
            fasta,
            ">tx1|gene:GENE_A|target:exon_transcript\nACGTACGT\n>tx2|gene:GENE_B|target:exon_transcript\nTTTTACGT"
        )
        .unwrap();
        let index = EcIndex::build_from_fasta(fasta.path(), 4).unwrap();
        let (mapping, stats) = index.map_read(b"ACGTACGT", false);
        assert!(stats.seed_lookups > 0);
        assert_eq!(
            index.gene_name(mapping.gene_scores[0].gene_id),
            Some("GENE_A")
        );
    }

    #[test]
    fn ec_index_round_trips_binary() {
        let mut fasta = tempfile::NamedTempFile::new().unwrap();
        writeln!(fasta, ">tx1|gene:GENE_A\nACGTACGT").unwrap();
        let index = EcIndex::build_from_fasta(fasta.path(), 3).unwrap();
        let path = tempfile::NamedTempFile::new().unwrap();
        index.save(path.path()).unwrap();
        let loaded = EcIndex::load(path.path()).unwrap();
        assert_eq!(loaded.k, 3);
        assert_eq!(loaded.genes, index.genes);
        assert_eq!(loaded.kmers, index.kmers);
    }

    #[test]
    fn ec_index_round_trips_mmap() {
        let mut fasta = tempfile::NamedTempFile::new().unwrap();
        writeln!(
            fasta,
            ">tx1|gene:GENE_A\nACGTACGT\n>tx2|gene:GENE_B\nTTTTACGT"
        )
        .unwrap();
        let index = EcIndex::build_from_fasta(fasta.path(), 3).unwrap();
        let path = tempfile::NamedTempFile::new().unwrap();
        index.save_mmap(path.path()).unwrap();
        let loaded = LoadedEcIndex::load(path.path()).unwrap();
        assert!(matches!(
            &loaded,
            LoadedEcIndex::Mmap(index) if index.header.version == EC_MMAP_VERSION_V4
        ));
        assert_eq!(loaded.k(), 3);
        assert_eq!(loaded.num_genes(), index.genes.len());
        assert_eq!(loaded.gene_name(0), Some("GENE_A"));
        assert_eq!(
            loaded.lookup(index.kmers[0].kmer_code),
            Some(index.kmers[0])
        );

        let (owned_mapping, owned_stats) = index.map_read_sparse(b"ACGTACGT", false, 1, u32::MAX);
        let (mmap_mapping, mmap_stats) = loaded.map_read_sparse(b"ACGTACGT", false, 1, u32::MAX);
        assert_eq!(mmap_stats, owned_stats);
        assert_eq!(mmap_mapping, owned_mapping);
    }

    #[test]
    fn ec_index_round_trips_v5_mmap() {
        let mut fasta = tempfile::NamedTempFile::new().unwrap();
        writeln!(
            fasta,
            ">tx1|gene:GENE_A\nACGTACGT\n>tx2|gene:GENE_B\nTTTTACGT"
        )
        .unwrap();
        let index = EcIndex::build_from_fasta(fasta.path(), 3).unwrap();
        let prefix_bits = EC_PREFIX_BITS.min(index.k.saturating_mul(2));
        let prefix_offsets = build_prefix_offsets(&index.kmers, index.k, prefix_bits).unwrap();
        let (gene_metas, gene_names) = test_gene_metadata(&index);
        let path = tempfile::NamedTempFile::new().unwrap();
        save_mmap_v5_from_slices(
            path.path(),
            index.k,
            &gene_metas,
            &gene_names,
            &index.kmers,
            &index.ec_offsets,
            &index.ec_gene_ids,
            &index.ec_class_bits,
            &prefix_offsets,
        )
        .unwrap();
        let loaded = LoadedEcIndex::load(path.path()).unwrap();
        assert!(matches!(
            &loaded,
            LoadedEcIndex::Mmap(index) if index.header.version == EC_MMAP_VERSION_V5
        ));
        assert_eq!(
            loaded.lookup(index.kmers[0].kmer_code),
            Some(index.kmers[0])
        );
        let (owned_mapping, owned_stats) = index.map_read_sparse(b"ACGTACGT", false, 1, u32::MAX);
        let (mmap_mapping, mmap_stats) = loaded.map_read_sparse(b"ACGTACGT", false, 1, u32::MAX);
        assert_eq!(mmap_stats, owned_stats);
        assert_eq!(mmap_mapping, owned_mapping);
    }

    #[test]
    fn ec_index_round_trips_v6_mmap() {
        let mut fasta = tempfile::NamedTempFile::new().unwrap();
        writeln!(
            fasta,
            ">tx1|gene:GENE_A\nACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT\n>tx2|gene:GENE_B\nACGTACGTACGTACGTACGTACGTTTTTACGTACGTACGT\n>tx3|gene:GENE_C\nTTTTACGTACGTACGTACGTACGTACGTACGTACGTACGT"
        )
        .unwrap();
        let index = EcIndex::build_from_fasta(fasta.path(), 10).unwrap();
        let prefix_bits = EC_PREFIX_BITS.min(index.k.saturating_mul(2));
        let prefix_offsets = build_prefix_offsets(&index.kmers, index.k, prefix_bits).unwrap();
        let (gene_metas, gene_names) = test_gene_metadata(&index);
        let path = tempfile::NamedTempFile::new().unwrap();
        save_mmap_v6_from_slices(
            path.path(),
            index.k,
            &gene_metas,
            &gene_names,
            &index.kmers,
            &index.ec_offsets,
            &index.ec_gene_ids,
            &index.ec_class_bits,
            &prefix_offsets,
        )
        .unwrap();
        let loaded = LoadedEcIndex::load(path.path()).unwrap();
        assert!(matches!(
            &loaded,
            LoadedEcIndex::Mmap(index) if index.header.version == EC_MMAP_VERSION_V6
        ));
        for record in &index.kmers {
            assert_eq!(loaded.lookup(record.kmer_code), Some(*record));
            assert!(loaded.lookup_hit(record.kmer_code).is_some());
        }
        assert_eq!(loaded.lookup(0), index.lookup(0));
        let seq = b"ACGTACGTACGTACGTACGT";
        let (owned_mapping, owned_stats) = index.map_read_sparse(seq, false, 1, u32::MAX);
        let (mmap_mapping, mmap_stats) = loaded.map_read_sparse(seq, false, 1, u32::MAX);
        assert_eq!(mmap_stats, owned_stats);
        assert_eq!(mmap_mapping, owned_mapping);
    }

    #[test]
    fn ec_index_v6_finds_suffixes_across_multiple_blocks() {
        let k = 16_u8;
        let prefix_bits = EC_PREFIX_BITS;
        let suffix_bits = k * 2 - prefix_bits;
        let prefix = 0x00ab_cdu64;
        let kmers = (0..96_u64)
            .map(|suffix| EcKmerRecord {
                kmer_code: (prefix << suffix_bits) | suffix,
                ec_id: 0,
                gene_df: 1,
            })
            .collect::<Vec<_>>();
        let index = EcIndex {
            k,
            genes: vec!["GENE_A".to_owned()],
            kmers,
            ec_offsets: vec![0, 1],
            ec_gene_ids: vec![0],
            ec_class_bits: vec![EC_CLASS_EXON],
        };
        let prefix_offsets = build_prefix_offsets(&index.kmers, index.k, prefix_bits).unwrap();
        let (gene_metas, gene_names) = test_gene_metadata(&index);
        let path = tempfile::NamedTempFile::new().unwrap();
        save_mmap_v6_from_slices(
            path.path(),
            index.k,
            &gene_metas,
            &gene_names,
            &index.kmers,
            &index.ec_offsets,
            &index.ec_gene_ids,
            &index.ec_class_bits,
            &prefix_offsets,
        )
        .unwrap();
        let loaded = LoadedEcIndex::load(path.path()).unwrap();
        assert!(matches!(
            &loaded,
            LoadedEcIndex::Mmap(index) if index.header.version == EC_MMAP_VERSION_V6
                && index.header_v6.is_some()
                && index.header_v6.unwrap().num_blocks == 3
        ));
        for record in &index.kmers {
            assert_eq!(loaded.lookup(record.kmer_code), Some(*record));
        }
        assert_eq!(loaded.lookup((prefix << suffix_bits) | 127), None);
    }

    #[test]
    fn ec_index_v6_finds_wide_suffixes_across_multiple_blocks() {
        let k = 31_u8;
        let prefix_bits = EC_PREFIX_BITS;
        let suffix_bits = k * 2 - prefix_bits;
        let prefix = 0x0034_56u64;
        let kmers = (0..96_u64)
            .map(|idx| EcKmerRecord {
                kmer_code: (prefix << suffix_bits) | (0x1_0000_0000 + idx * 257),
                ec_id: 0,
                gene_df: 1,
            })
            .collect::<Vec<_>>();
        let index = EcIndex {
            k,
            genes: vec!["GENE_A".to_owned()],
            kmers,
            ec_offsets: vec![0, 1],
            ec_gene_ids: vec![0],
            ec_class_bits: vec![EC_CLASS_EXON],
        };
        let prefix_offsets = build_prefix_offsets(&index.kmers, index.k, prefix_bits).unwrap();
        let (gene_metas, gene_names) = test_gene_metadata(&index);
        let path = tempfile::NamedTempFile::new().unwrap();
        save_mmap_v6_from_slices(
            path.path(),
            index.k,
            &gene_metas,
            &gene_names,
            &index.kmers,
            &index.ec_offsets,
            &index.ec_gene_ids,
            &index.ec_class_bits,
            &prefix_offsets,
        )
        .unwrap();
        let loaded = LoadedEcIndex::load(path.path()).unwrap();
        for record in &index.kmers {
            assert_eq!(loaded.lookup(record.kmer_code), Some(*record));
        }
        assert_eq!(loaded.lookup((prefix << suffix_bits) | 0x1_0000_ffff), None);
    }

    #[test]
    fn ec_index_converts_v5_to_v6_mmap() {
        let k = 31_u8;
        let prefix_bits = EC_PREFIX_BITS;
        let suffix_bits = k * 2 - prefix_bits;
        let prefix = 0x0034_56u64;
        let kmers = (0..96_u64)
            .map(|idx| EcKmerRecord {
                kmer_code: (prefix << suffix_bits) | (0x1_0000_0000 + idx * 257),
                ec_id: 0,
                gene_df: 1,
            })
            .collect::<Vec<_>>();
        let index = EcIndex {
            k,
            genes: vec!["GENE_A".to_owned()],
            kmers,
            ec_offsets: vec![0, 1],
            ec_gene_ids: vec![0],
            ec_class_bits: vec![EC_CLASS_EXON],
        };
        let prefix_bits = EC_PREFIX_BITS.min(index.k.saturating_mul(2));
        let prefix_offsets = build_prefix_offsets(&index.kmers, index.k, prefix_bits).unwrap();
        let (gene_metas, gene_names) = test_gene_metadata(&index);
        let v5_path = tempfile::NamedTempFile::new().unwrap();
        save_mmap_v5_from_slices(
            v5_path.path(),
            index.k,
            &gene_metas,
            &gene_names,
            &index.kmers,
            &index.ec_offsets,
            &index.ec_gene_ids,
            &index.ec_class_bits,
            &prefix_offsets,
        )
        .unwrap();
        let v5 = match LoadedEcIndex::load(v5_path.path()).unwrap() {
            LoadedEcIndex::Mmap(index) => index,
            LoadedEcIndex::Owned(_) => panic!("expected mmap index"),
        };
        let v6_path = tempfile::NamedTempFile::new().unwrap();
        save_mmap_v6_from_mmap_v5(v6_path.path(), &v5).unwrap();
        let loaded = LoadedEcIndex::load(v6_path.path()).unwrap();
        assert!(matches!(
            &loaded,
            LoadedEcIndex::Mmap(index) if index.header.version == EC_MMAP_VERSION_V6
        ));
        for record in &index.kmers {
            assert_eq!(loaded.lookup(record.kmer_code), Some(*record));
        }
        assert_eq!(loaded.lookup((prefix << suffix_bits) | 0x1_0000_ffff), None);
    }

    #[test]
    fn ec_index_v7_dedups_target_suffixes_without_changing_scores() {
        let k = 16_u8;
        let prefix_bits = EC_PREFIX_BITS;
        let suffix_bits = k * 2 - prefix_bits;
        let prefix = 0x0034_56u64;
        let mut ec_offsets = Vec::new();
        let mut ec_gene_ids = Vec::new();
        let mut ec_class_bits = Vec::new();
        ec_offsets.push(0);
        for start_gene in 0..16_u32 {
            for gene_id in start_gene..20 {
                ec_gene_ids.push(gene_id);
                ec_class_bits.push(if gene_id % 2 == 0 {
                    EC_CLASS_EXON
                } else {
                    EC_CLASS_INTRON
                });
            }
            ec_offsets.push(ec_gene_ids.len() as u32);
        }
        let kmers = (0..16_u64)
            .map(|idx| {
                let start = ec_offsets[idx as usize] as usize;
                let end = ec_offsets[idx as usize + 1] as usize;
                EcKmerRecord {
                    kmer_code: (prefix << suffix_bits) | idx,
                    ec_id: idx as u32,
                    gene_df: (end - start) as u32,
                }
            })
            .collect::<Vec<_>>();
        let index = EcIndex {
            k,
            genes: (0..20).map(|idx| format!("GENE_{idx}")).collect(),
            kmers,
            ec_offsets,
            ec_gene_ids,
            ec_class_bits,
        };
        let prefix_offsets = build_prefix_offsets(&index.kmers, index.k, prefix_bits).unwrap();
        let (gene_metas, gene_names) = test_gene_metadata(&index);
        let v6_path = tempfile::NamedTempFile::new().unwrap();
        let v7_path = tempfile::NamedTempFile::new().unwrap();
        save_mmap_v6_from_slices(
            v6_path.path(),
            index.k,
            &gene_metas,
            &gene_names,
            &index.kmers,
            &index.ec_offsets,
            &index.ec_gene_ids,
            &index.ec_class_bits,
            &prefix_offsets,
        )
        .unwrap();
        save_mmap_v7_from_slices(
            v7_path.path(),
            index.k,
            &gene_metas,
            &gene_names,
            &index.kmers,
            &index.ec_offsets,
            &index.ec_gene_ids,
            &index.ec_class_bits,
            &prefix_offsets,
        )
        .unwrap();
        assert!(
            std::fs::metadata(v7_path.path()).unwrap().len()
                < std::fs::metadata(v6_path.path()).unwrap().len()
        );
        let loaded = LoadedEcIndex::load(v7_path.path()).unwrap();
        let LoadedEcIndex::Mmap(mmap) = &loaded else {
            panic!("expected mmap index");
        };
        assert_eq!(mmap.header.version, EC_MMAP_VERSION);
        assert!(mmap.header_v7.is_some());
        assert!(mmap.header_v7.unwrap().num_ec_targets < index.ec_gene_ids.len() as u64);

        let owned = LoadedEcIndex::Owned(index.clone());
        for record in &index.kmers {
            assert_eq!(loaded.lookup(record.kmer_code), Some(*record));
            let hit = loaded.lookup_hit(record.kmer_code).expect("lookup hit");
            let mut record_scores = FxHashMap::default();
            let mut hit_scores = FxHashMap::default();
            let mut record_stats = EcMappingStats::default();
            let mut hit_stats = EcMappingStats::default();
            owned.accumulate_lookup_record(
                *record,
                u32::MAX,
                &mut record_scores,
                &mut record_stats,
            );
            loaded.accumulate_lookup_hit(hit, u32::MAX, &mut hit_scores, &mut hit_stats);
            assert_eq!(hit_stats, record_stats);
            assert_eq!(hit_scores, record_scores);
        }
    }

    #[test]
    fn ec_index_v8_mph_round_trips_lookup_and_mapping() {
        let mut fasta = tempfile::NamedTempFile::new().unwrap();
        writeln!(
            fasta,
            ">tx1|gene:GENE_A\nACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT\n>tx2|gene:GENE_B\nACGTACGTACGTACGTACGTACGTTTTTACGTACGTACGT\n>tx3|gene:GENE_C\nTTTTACGTACGTACGTACGTACGTACGTACGTACGTACGT"
        )
        .unwrap();
        let index = EcIndex::build_from_fasta(fasta.path(), 10).unwrap();
        let prefix_bits = EC_PREFIX_BITS.min(index.k.saturating_mul(2));
        let prefix_offsets = build_prefix_offsets(&index.kmers, index.k, prefix_bits).unwrap();
        let (gene_metas, gene_names) = test_gene_metadata(&index);
        let path = tempfile::NamedTempFile::new().unwrap();
        save_mmap_v8_from_slices(
            path.path(),
            index.k,
            &gene_metas,
            &gene_names,
            &index.kmers,
            &index.ec_offsets,
            &index.ec_gene_ids,
            &index.ec_class_bits,
            &prefix_offsets,
        )
        .unwrap();
        let loaded = LoadedEcIndex::load(path.path()).unwrap();
        let LoadedEcIndex::Mmap(mmap) = &loaded else {
            panic!("expected mmap index");
        };
        assert_eq!(mmap.header.version, EC_MMAP_VERSION_V8);
        assert!(mmap.header_v8.is_some());
        for record in &index.kmers {
            assert_eq!(loaded.lookup(record.kmer_code), Some(*record));
            assert!(loaded.lookup_hit(record.kmer_code).is_some());
        }
        assert_eq!(loaded.lookup(0), index.lookup(0));
        let seq = b"ACGTACGTACGTACGTACGT";
        let (owned_mapping, owned_stats) = index.map_read_sparse(seq, false, 1, u32::MAX);
        let (mmap_mapping, mmap_stats) = loaded.map_read_sparse(seq, false, 1, u32::MAX);
        assert_eq!(mmap_stats, owned_stats);
        assert_eq!(mmap_mapping, owned_mapping);
    }

    #[test]
    fn ec_lookup_hit_accumulates_like_lookup_record() {
        let mut fasta = tempfile::NamedTempFile::new().unwrap();
        writeln!(
            fasta,
            ">tx1|gene:GENE_A\nACGTACGT\n>tx2|gene:GENE_B\nTTTTACGT"
        )
        .unwrap();
        let index = EcIndex::build_from_fasta(fasta.path(), 3).unwrap();
        let prefix_bits = EC_PREFIX_BITS.min(index.k.saturating_mul(2));
        let prefix_offsets = build_prefix_offsets(&index.kmers, index.k, prefix_bits).unwrap();
        let (gene_metas, gene_names) = test_gene_metadata(&index);
        let path = tempfile::NamedTempFile::new().unwrap();
        save_mmap_v4_from_slices(
            path.path(),
            index.k,
            &gene_metas,
            &gene_names,
            &index.kmers,
            &index.ec_offsets,
            &index.ec_gene_ids,
            &index.ec_class_bits,
            &prefix_offsets,
        )
        .unwrap();
        let loaded = LoadedEcIndex::load(path.path()).unwrap();
        let encoded = encode_acgt(b"ACGTACGTACGT");
        let mut record_scores = FxHashMap::default();
        let mut hit_scores = FxHashMap::default();
        let mut record_stats = EcMappingStats::default();
        let mut hit_stats = EcMappingStats::default();
        for kmer in iter_kmers_2bit(&encoded, loaded.k()) {
            let Some(record) = loaded.lookup(kmer.code) else {
                continue;
            };
            let hit = loaded.lookup_hit(kmer.code).expect("lookup hit");
            loaded.accumulate_lookup_record(
                record,
                u32::MAX,
                &mut record_scores,
                &mut record_stats,
            );
            loaded.accumulate_lookup_hit(hit, u32::MAX, &mut hit_scores, &mut hit_stats);
        }
        assert_eq!(hit_stats, record_stats);
        assert_eq!(hit_scores, record_scores);
    }

    #[test]
    fn ec_index_encoded_mapping_matches_ascii_mapping() {
        let mut fasta = tempfile::NamedTempFile::new().unwrap();
        writeln!(
            fasta,
            ">tx1|gene:GENE_A\nACGTACGT\n>tx2|gene:GENE_B\nTTTTACGT"
        )
        .unwrap();
        let index = EcIndex::build_from_fasta(fasta.path(), 3).unwrap();
        let path = tempfile::NamedTempFile::new().unwrap();
        index.save_mmap(path.path()).unwrap();
        let loaded = LoadedEcIndex::load(path.path()).unwrap();
        let seq = b"ACGTACGT";
        let encoded = encode_acgt(seq);
        assert_eq!(
            loaded.map_encoded_sparse(&encoded, false, 1, u32::MAX),
            loaded.map_read_sparse(seq, false, 1, u32::MAX)
        );
        assert_eq!(
            loaded.map_encoded_sparse(&encoded, true, 1, u32::MAX),
            loaded.map_read_sparse(seq, true, 1, u32::MAX)
        );

        let mut scratch = EcScoreScratch::new(loaded.num_genes());
        for (seq, reverse_complement, stride, max_gene_df) in [
            (b"ACGTACGT".as_slice(), false, 1, u32::MAX),
            (b"TTTTACGTACGT".as_slice(), true, 1, u32::MAX),
            (b"ACGTNNNNACGT".as_slice(), false, 2, 1),
            (b"GGGGGGGG".as_slice(), true, 3, u32::MAX),
        ] {
            let encoded = encode_acgt(seq);
            assert_eq!(
                loaded.map_encoded_sparse_with_scratch(
                    &encoded,
                    reverse_complement,
                    stride,
                    max_gene_df,
                    &mut scratch,
                ),
                loaded.map_encoded_sparse(&encoded, reverse_complement, stride, max_gene_df),
                "sequence={:?}",
                String::from_utf8_lossy(seq),
            );
        }
    }

    #[test]
    fn ec_index_loads_v3_mmap() {
        let mut fasta = tempfile::NamedTempFile::new().unwrap();
        writeln!(
            fasta,
            ">tx1|gene:GENE_A\nACGTACGT\n>tx2|gene:GENE_B\nTTTTACGT"
        )
        .unwrap();
        let index = EcIndex::build_from_fasta(fasta.path(), 3).unwrap();
        let path = tempfile::NamedTempFile::new().unwrap();
        index.save_mmap_v3(path.path()).unwrap();
        let loaded = LoadedEcIndex::load(path.path()).unwrap();
        assert!(matches!(
            &loaded,
            LoadedEcIndex::Mmap(index) if index.header.version == EC_MMAP_VERSION_V3
        ));
        assert_eq!(
            loaded.lookup(index.kmers[0].kmer_code),
            Some(index.kmers[0])
        );
        let (owned_mapping, owned_stats) = index.map_read_sparse(b"ACGTACGT", false, 1, u32::MAX);
        let (mmap_mapping, mmap_stats) = loaded.map_read_sparse(b"ACGTACGT", false, 1, u32::MAX);
        assert_eq!(mmap_stats, owned_stats);
        assert_eq!(mmap_mapping, owned_mapping);
    }

    #[test]
    fn ec_index_loads_v4_mmap() {
        let mut fasta = tempfile::NamedTempFile::new().unwrap();
        writeln!(
            fasta,
            ">tx1|gene:GENE_A\nACGTACGT\n>tx2|gene:GENE_B\nTTTTACGT"
        )
        .unwrap();
        let index = EcIndex::build_from_fasta(fasta.path(), 3).unwrap();
        let prefix_bits = EC_PREFIX_BITS.min(index.k.saturating_mul(2));
        let prefix_offsets = build_prefix_offsets(&index.kmers, index.k, prefix_bits).unwrap();
        let (gene_metas, gene_names) = test_gene_metadata(&index);
        let path = tempfile::NamedTempFile::new().unwrap();
        save_mmap_v4_from_slices(
            path.path(),
            index.k,
            &gene_metas,
            &gene_names,
            &index.kmers,
            &index.ec_offsets,
            &index.ec_gene_ids,
            &index.ec_class_bits,
            &prefix_offsets,
        )
        .unwrap();
        let loaded = LoadedEcIndex::load(path.path()).unwrap();
        assert!(matches!(
            &loaded,
            LoadedEcIndex::Mmap(index) if index.header.version == EC_MMAP_VERSION_V4
        ));
        assert_eq!(
            loaded.lookup(index.kmers[0].kmer_code),
            Some(index.kmers[0])
        );
        let (owned_mapping, owned_stats) = index.map_read_sparse(b"ACGTACGT", false, 1, u32::MAX);
        let (mmap_mapping, mmap_stats) = loaded.map_read_sparse(b"ACGTACGT", false, 1, u32::MAX);
        assert_eq!(mmap_stats, owned_stats);
        assert_eq!(mmap_mapping, owned_mapping);
    }

    fn test_gene_metadata(index: &EcIndex) -> (Vec<EcGeneMeta>, Vec<u8>) {
        let mut gene_names = Vec::new();
        let mut gene_metas = Vec::with_capacity(index.genes.len());
        for gene in &index.genes {
            let name_start = gene_names.len() as u32;
            let name_len = gene.len() as u32;
            gene_names.extend_from_slice(gene.as_bytes());
            gene_metas.push(EcGeneMeta {
                name_start,
                name_len,
            });
        }
        (gene_metas, gene_names)
    }
}

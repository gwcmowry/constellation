use crate::dna::{encode_acgt, iter_kmers_2bit, reverse_complement};
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
const EC_MMAP_VERSION: u32 = 3;
const EC_PREFIX_BITS: u8 = 24;
const EC_BUILD_BUFFER_RECORDS: usize = 8_000_000;
const EC_CLASS_UNKNOWN: u8 = 0;
const EC_CLASS_EXON: u8 = 1;
const EC_CLASS_INTRON: u8 = 2;
const EC_CLASS_ANTISENSE: u8 = 4;

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

#[derive(Debug, Clone, Copy, Default, PartialEq, Eq, Zeroable, Pod)]
#[repr(C)]
struct EcGeneMeta {
    name_start: u32,
    name_len: u32,
}

#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub struct EcMappingStats {
    pub seed_lookups: u64,
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
            version: EC_MMAP_VERSION,
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

    fn seed_score(&self, gene_df: u32) -> u16 {
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

    pub fn map_read_sparse(
        &self,
        seq: &[u8],
        search_reverse_complement: bool,
        stride: u32,
        max_gene_df: u32,
    ) -> (EcReadMapping, EcMappingStats) {
        let mut stats = EcMappingStats::default();
        let mut gene_scores: FxHashMap<GeneId, EcGeneScore> = FxHashMap::default();
        accumulate_ec_seq(
            self,
            seq,
            stride.max(1),
            max_gene_df.max(1),
            &mut gene_scores,
            &mut stats,
        );
        if search_reverse_complement {
            let rc = reverse_complement(&encode_acgt(seq));
            let rc = crate::dna::decode_to_vec(&rc);
            accumulate_ec_seq(
                self,
                &rc,
                stride.max(1),
                max_gene_df.max(1),
                &mut gene_scores,
                &mut stats,
            );
        }
        finish_ec_mapping(gene_scores, stats)
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
        let header = if header_v1.version == 1 {
            EcMmapHeader::from(header_v1)
        } else if header_v1.version == 2 {
            if mmap.len() < std::mem::size_of::<EcMmapHeaderV2>() {
                return Err(EcIndexError::Truncated);
            }
            let header_v2: EcMmapHeaderV2 =
                *bytemuck::from_bytes(&mmap[..std::mem::size_of::<EcMmapHeaderV2>()]);
            EcMmapHeader::from(header_v2)
        } else if header_v1.version == EC_MMAP_VERSION {
            if mmap.len() < std::mem::size_of::<EcMmapHeader>() {
                return Err(EcIndexError::Truncated);
            }
            *bytemuck::from_bytes(&mmap[..std::mem::size_of::<EcMmapHeader>()])
        } else {
            return Err(EcIndexError::UnsupportedVersion(header_v1.version));
        };
        let mut this = Self {
            mmap,
            header,
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
        let prefix_bits = EC_PREFIX_BITS.min(self.k().saturating_mul(2));
        let prefix_offsets = build_prefix_offsets(self.kmers(), self.k(), prefix_bits)?;
        let class_bits = if self.header.version >= 2 {
            self.class_bits()
        } else {
            &[]
        };

        let mut cursor = align_up(std::mem::size_of::<EcMmapHeader>() as u64, 8);
        let prefix_offsets_offset = cursor;
        cursor += byte_len::<u32>(prefix_offsets.len());
        cursor = align_up(cursor, 8);
        let genes_offset = cursor;
        cursor += byte_len::<EcGeneMeta>(self.genes().len());
        cursor = align_up(cursor, 8);
        let kmers_offset = cursor;
        cursor += byte_len::<EcKmerRecord>(self.kmers().len());
        cursor = align_up(cursor, 4);
        let ec_offsets_offset = cursor;
        cursor += byte_len::<u32>(self.ec_offsets().len());
        cursor = align_up(cursor, 4);
        let ec_gene_ids_offset = cursor;
        cursor += byte_len::<u32>(self.ec_gene_ids().len());
        cursor = align_up(cursor, 8);
        let ec_class_bits_offset = cursor;
        cursor += class_bits.len() as u64;
        cursor = align_up(cursor, 8);
        let gene_names_offset = cursor;

        let header = EcMmapHeader {
            magic: EC_MMAP_MAGIC,
            version: EC_MMAP_VERSION,
            k: self.k() as u32,
            prefix_bits: prefix_bits as u32,
            _pad0: 0,
            num_genes: self.header.num_genes,
            num_kmers: self.header.num_kmers,
            num_ec_offsets: self.header.num_ec_offsets,
            num_ec_gene_ids: self.header.num_ec_gene_ids,
            genes_offset,
            kmers_offset,
            ec_offsets_offset,
            ec_gene_ids_offset,
            gene_names_offset,
            gene_names_len: self.header.gene_names_len,
            ec_class_bits_offset,
            ec_class_bits_len: class_bits.len() as u64,
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
        write_pod_slice(&mut file, self.genes())?;
        write_padding(
            &mut file,
            genes_offset + byte_len::<EcGeneMeta>(self.genes().len()),
            kmers_offset,
        )?;
        write_pod_slice(&mut file, self.kmers())?;
        write_padding(
            &mut file,
            kmers_offset + byte_len::<EcKmerRecord>(self.kmers().len()),
            ec_offsets_offset,
        )?;
        write_pod_slice(&mut file, self.ec_offsets())?;
        write_padding(
            &mut file,
            ec_offsets_offset + byte_len::<u32>(self.ec_offsets().len()),
            ec_gene_ids_offset,
        )?;
        write_pod_slice(&mut file, self.ec_gene_ids())?;
        write_padding(
            &mut file,
            ec_gene_ids_offset + byte_len::<u32>(self.ec_gene_ids().len()),
            ec_class_bits_offset,
        )?;
        file.write_all(class_bits)?;
        write_padding(
            &mut file,
            ec_class_bits_offset + class_bits.len() as u64,
            gene_names_offset,
        )?;
        file.write_all(self.gene_names())?;
        Ok(())
    }

    fn validate(&self) -> Result<(), EcIndexError> {
        self.slice::<EcGeneMeta>(self.header.genes_offset, self.header.num_genes)?;
        self.slice::<EcKmerRecord>(self.header.kmers_offset, self.header.num_kmers)?;
        self.slice::<u32>(self.header.ec_offsets_offset, self.header.num_ec_offsets)?;
        self.slice::<u32>(self.header.ec_gene_ids_offset, self.header.num_ec_gene_ids)?;
        if self.header.version >= 2 {
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
        let kmers = self.kmer_lookup_range(kmer_code);
        kmers
            .binary_search_by_key(&kmer_code, |record| record.kmer_code)
            .ok()
            .map(|idx| kmers[idx])
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
        let offsets = self.ec_offsets();
        let start = offsets.get(ec_id as usize).copied().unwrap_or(0) as usize;
        let end = offsets
            .get(ec_id as usize + 1)
            .copied()
            .unwrap_or(start as u32) as usize;
        &self.ec_gene_ids()[start..end]
    }

    fn ec_class_bits(&self, ec_id: u32) -> &[u8] {
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

    fn genes(&self) -> &[EcGeneMeta] {
        self.slice(self.header.genes_offset, self.header.num_genes)
            .expect("validated EC gene metadata")
    }

    fn kmers(&self) -> &[EcKmerRecord] {
        self.slice(self.header.kmers_offset, self.header.num_kmers)
            .expect("validated EC kmers")
    }

    fn ec_offsets(&self) -> &[u32] {
        self.slice(self.header.ec_offsets_offset, self.header.num_ec_offsets)
            .expect("validated EC offsets")
    }

    fn ec_gene_ids(&self) -> &[GeneId] {
        self.slice(self.header.ec_gene_ids_offset, self.header.num_ec_gene_ids)
            .expect("validated EC gene ids")
    }

    fn gene_names(&self) -> &[u8] {
        self.bytes(self.header.gene_names_offset, self.header.gene_names_len)
            .expect("validated EC gene names")
    }

    fn class_bits(&self) -> &[u8] {
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
}

fn accumulate_ec_seq(
    index: &LoadedEcIndex,
    seq: &[u8],
    stride: u32,
    max_gene_df: u32,
    gene_scores: &mut FxHashMap<GeneId, EcGeneScore>,
    stats: &mut EcMappingStats,
) {
    let encoded = encode_acgt(seq);
    for kmer in iter_kmers_2bit(&encoded, index.k()) {
        if kmer.pos % stride != 0 {
            continue;
        }
        stats.seed_lookups += 1;
        let Some(record) = index.lookup(kmer.code) else {
            continue;
        };
        if record.gene_df > max_gene_df {
            stats.seeds_skipped_high_gene_df += 1;
            continue;
        }
        stats.seeds_found += 1;
        stats.ec_hits += 1;
        let seed_score = index.seed_score(record.gene_df);
        let genes = index.ec_genes(record.ec_id);
        let class_bits = index.ec_class_bits(record.ec_id);
        for (idx, &gene_id) in genes.iter().enumerate() {
            let class_bits = class_bits.get(idx).copied().unwrap_or(EC_CLASS_UNKNOWN);
            add_gene_score(gene_scores, gene_id, seed_score, class_bits);
            stats.gene_votes += 1;
        }
    }
}

fn finish_ec_mapping(
    gene_scores: FxHashMap<GeneId, EcGeneScore>,
    mut stats: EcMappingStats,
) -> (EcReadMapping, EcMappingStats) {
    let mut gene_scores: Vec<_> = gene_scores.into_values().collect();
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
}

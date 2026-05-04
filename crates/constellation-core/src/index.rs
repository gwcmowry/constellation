use crate::{GeneId, TranscriptId};
use bytemuck::{Pod, Zeroable};
use memmap2::Mmap;
use serde::{Deserialize, Serialize};
use std::fs::File;
use std::io::Write;
use std::path::Path;
use thiserror::Error;

const LEGACY_MAGIC: &[u8; 8] = b"CSTLIDX1";
const COMPACT_MAGIC: [u8; 8] = *b"CSTLCB2\0";
const COMPACT_VERSION: u32 = 2;

#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub struct TranscriptMeta {
    pub transcript_id: TranscriptId,
    pub gene_id: GeneId,
    pub name: String,
    pub len: u32,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub struct Posting {
    pub transcript_id: TranscriptId,
    pub gene_id: GeneId,
    pub pos: u32,
    pub strand: u8,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct SeedPosting {
    pub transcript_id: TranscriptId,
    pub pos: u32,
    pub strand: u8,
}

#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub struct KmerEntry {
    pub kmer_code: u64,
    pub postings_start: u64,
    pub postings_len: u32,
    pub freq_class: u8,
}

#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub struct TranscriptIndex {
    pub format_version: u32,
    pub k: u8,
    pub max_kmer_frequency: u32,
    pub genes: Vec<String>,
    pub transcripts: Vec<TranscriptMeta>,
    pub transcript_sequences: Vec<String>,
    pub kmers: Vec<KmerEntry>,
    pub postings: Vec<Posting>,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct IndexStats {
    pub num_transcripts: usize,
    pub num_genes: usize,
    pub num_distinct_kmers: usize,
    pub num_postings: usize,
    pub k: u8,
    pub max_postings_per_kmer: u32,
    pub high_frequency_kmers: usize,
}

#[derive(Debug, Error)]
pub enum IndexIoError {
    #[error("index I/O error: {0}")]
    Io(#[from] std::io::Error),
    #[error("index JSON error: {0}")]
    Json(#[from] serde_json::Error),
    #[error("invalid index magic")]
    InvalidMagic,
    #[error("truncated index file")]
    Truncated,
    #[error("misaligned compact index section")]
    Misaligned,
    #[error("compact index version {0} is unsupported")]
    UnsupportedVersion(u32),
    #[error("index is too large for compact field: {0}")]
    TooLarge(&'static str),
}

#[derive(Debug, Clone, Copy, Zeroable, Pod)]
#[repr(C)]
struct CompactHeader {
    magic: [u8; 8],
    version: u32,
    k: u32,
    max_kmer_frequency: u32,
    _pad0: u32,
    num_genes: u64,
    num_transcripts: u64,
    num_kmers: u64,
    num_postings: u64,
    genes_offset: u64,
    transcripts_offset: u64,
    kmers_offset: u64,
    postings_offset: u64,
    gene_names_offset: u64,
    gene_names_len: u64,
    transcript_names_offset: u64,
    transcript_names_len: u64,
    sequences_offset: u64,
    sequences_len: u64,
}

#[derive(Debug, Clone, Copy, Zeroable, Pod)]
#[repr(C)]
struct CompactGeneMeta {
    name_start: u32,
    name_len: u32,
}

#[derive(Debug, Clone, Copy, Zeroable, Pod)]
#[repr(C)]
struct CompactTranscriptMeta {
    transcript_id: u32,
    gene_id: u32,
    name_start: u32,
    name_len: u32,
    seq_start: u64,
    len: u32,
    _pad0: u32,
}

#[derive(Debug, Clone, Copy, Zeroable, Pod)]
#[repr(C)]
struct CompactKmerEntry {
    kmer_code: u64,
    postings_start: u32,
    postings_len: u32,
}

#[derive(Debug, Clone, Copy, Zeroable, Pod)]
#[repr(C)]
struct CompactPosting {
    transcript_id: u32,
    pos_strand: u32,
}

impl CompactPosting {
    fn from_posting(value: Posting) -> Result<Self, IndexIoError> {
        if value.pos > (u32::MAX >> 1) {
            return Err(IndexIoError::TooLarge("posting position"));
        }
        Ok(Self {
            transcript_id: value.transcript_id,
            pos_strand: (value.pos << 1) | u32::from(value.strand & 1),
        })
    }

    fn to_posting(self, transcripts: &[CompactTranscriptMeta]) -> Posting {
        let transcript_id = self.transcript_id;
        let gene_id = transcripts
            .get(transcript_id as usize)
            .map_or(0, |meta| meta.gene_id);
        Posting {
            transcript_id,
            gene_id,
            pos: self.pos_strand >> 1,
            strand: (self.pos_strand & 1) as u8,
        }
    }
}

#[derive(Debug)]
pub struct CompactTranscriptIndex {
    mmap: Mmap,
    header: CompactHeader,
}

#[derive(Debug)]
pub enum LoadedIndex {
    Owned(TranscriptIndex),
    Compact(CompactTranscriptIndex),
}

pub struct PostingIter<'a>(PostingIterInner<'a>);

enum PostingIterInner<'a> {
    Owned(std::iter::Copied<std::slice::Iter<'a, Posting>>),
    Compact {
        postings: std::iter::Copied<std::slice::Iter<'a, CompactPosting>>,
        transcripts: &'a [CompactTranscriptMeta],
    },
    Empty,
}

pub struct SeedPostingIter<'a>(SeedPostingIterInner<'a>);

enum SeedPostingIterInner<'a> {
    Owned(std::iter::Copied<std::slice::Iter<'a, Posting>>),
    Compact(std::iter::Copied<std::slice::Iter<'a, CompactPosting>>),
    Empty,
}

impl Iterator for PostingIter<'_> {
    type Item = Posting;

    fn next(&mut self) -> Option<Self::Item> {
        match &mut self.0 {
            PostingIterInner::Owned(iter) => iter.next(),
            PostingIterInner::Compact {
                postings,
                transcripts,
            } => postings
                .next()
                .map(|posting| posting.to_posting(transcripts)),
            PostingIterInner::Empty => None,
        }
    }
}

impl Iterator for SeedPostingIter<'_> {
    type Item = SeedPosting;

    fn next(&mut self) -> Option<Self::Item> {
        match &mut self.0 {
            SeedPostingIterInner::Owned(iter) => iter.next().map(|posting| SeedPosting {
                transcript_id: posting.transcript_id,
                pos: posting.pos,
                strand: posting.strand,
            }),
            SeedPostingIterInner::Compact(iter) => iter.next().map(|posting| SeedPosting {
                transcript_id: posting.transcript_id,
                pos: posting.pos_strand >> 1,
                strand: (posting.pos_strand & 1) as u8,
            }),
            SeedPostingIterInner::Empty => None,
        }
    }
}

pub trait IndexAccess: Send + Sync {
    fn k(&self) -> u8;
    fn max_kmer_frequency(&self) -> u32;
    fn num_genes(&self) -> usize;
    fn gene_name(&self, gene_id: GeneId) -> Option<&str>;
    fn transcript_gene_id(&self, transcript_id: TranscriptId) -> Option<GeneId>;
    fn transcript_seq(&self, transcript_id: TranscriptId) -> Option<&[u8]>;
    fn posting_count(&self, kmer_code: u64) -> usize;
    fn seed_postings(&self, kmer_code: u64) -> SeedPostingIter<'_>;
    fn postings(&self, kmer_code: u64) -> PostingIter<'_>;
    fn stats(&self) -> IndexStats;
}

impl TranscriptIndex {
    pub fn save_json(&self, path: impl AsRef<Path>) -> Result<(), IndexIoError> {
        let bytes = serde_json::to_vec_pretty(self)?;
        std::fs::write(path, bytes)?;
        Ok(())
    }

    pub fn save_auto(&self, path: impl AsRef<Path>) -> Result<(), IndexIoError> {
        let path = path.as_ref();
        if path
            .extension()
            .is_some_and(|extension| extension == "json")
        {
            return self.save_json(path);
        }
        self.save_compact(path)
    }

    pub fn save_compact(&self, path: impl AsRef<Path>) -> Result<(), IndexIoError> {
        let mut gene_names = Vec::new();
        let mut genes = Vec::with_capacity(self.genes.len());
        for gene in &self.genes {
            let start = gene_names.len() as u32;
            gene_names.extend_from_slice(gene.as_bytes());
            genes.push(CompactGeneMeta {
                name_start: start,
                name_len: gene.len() as u32,
            });
        }

        let mut transcript_names = Vec::new();
        let mut sequences = Vec::new();
        let mut transcripts = Vec::with_capacity(self.transcripts.len());
        for meta in &self.transcripts {
            let name_start = transcript_names.len() as u32;
            transcript_names.extend_from_slice(meta.name.as_bytes());
            let seq_start = sequences.len() as u64;
            let seq = self
                .transcript_sequences
                .get(meta.transcript_id as usize)
                .map_or(&b""[..], |seq| seq.as_bytes());
            sequences.extend(seq.iter().map(|base| base.to_ascii_uppercase()));
            transcripts.push(CompactTranscriptMeta {
                transcript_id: meta.transcript_id,
                gene_id: meta.gene_id,
                name_start,
                name_len: meta.name.len() as u32,
                seq_start,
                len: seq.len() as u32,
                _pad0: 0,
            });
        }

        let mut kmers = Vec::with_capacity(self.kmers.len());
        for entry in &self.kmers {
            let postings_start: u32 = entry
                .postings_start
                .try_into()
                .map_err(|_| IndexIoError::TooLarge("postings_start"))?;
            kmers.push(CompactKmerEntry {
                kmer_code: entry.kmer_code,
                postings_start,
                postings_len: entry.postings_len,
            });
        }
        let mut postings = Vec::with_capacity(self.postings.len());
        for posting in &self.postings {
            postings.push(CompactPosting::from_posting(*posting)?);
        }

        let mut cursor = align_up(std::mem::size_of::<CompactHeader>() as u64, 8);
        let genes_offset = cursor;
        cursor += byte_len::<CompactGeneMeta>(genes.len());
        cursor = align_up(cursor, 8);
        let transcripts_offset = cursor;
        cursor += byte_len::<CompactTranscriptMeta>(transcripts.len());
        cursor = align_up(cursor, 8);
        let kmers_offset = cursor;
        cursor += byte_len::<CompactKmerEntry>(kmers.len());
        cursor = align_up(cursor, 8);
        let postings_offset = cursor;
        cursor += byte_len::<CompactPosting>(postings.len());
        cursor = align_up(cursor, 8);
        let gene_names_offset = cursor;
        cursor += gene_names.len() as u64;
        cursor = align_up(cursor, 8);
        let transcript_names_offset = cursor;
        cursor += transcript_names.len() as u64;
        cursor = align_up(cursor, 8);
        let sequences_offset = cursor;

        let header = CompactHeader {
            magic: COMPACT_MAGIC,
            version: COMPACT_VERSION,
            k: self.k as u32,
            max_kmer_frequency: self.max_kmer_frequency,
            _pad0: 0,
            num_genes: genes.len() as u64,
            num_transcripts: transcripts.len() as u64,
            num_kmers: kmers.len() as u64,
            num_postings: postings.len() as u64,
            genes_offset,
            transcripts_offset,
            kmers_offset,
            postings_offset,
            gene_names_offset,
            gene_names_len: gene_names.len() as u64,
            transcript_names_offset,
            transcript_names_len: transcript_names.len() as u64,
            sequences_offset,
            sequences_len: sequences.len() as u64,
        };

        let mut file = File::create(path)?;
        file.write_all(bytemuck::bytes_of(&header))?;
        write_padding(
            &mut file,
            std::mem::size_of::<CompactHeader>() as u64,
            genes_offset,
        )?;
        write_pod_slice(&mut file, &genes)?;
        write_padding(
            &mut file,
            genes_offset + byte_len::<CompactGeneMeta>(genes.len()),
            transcripts_offset,
        )?;
        write_pod_slice(&mut file, &transcripts)?;
        write_padding(
            &mut file,
            transcripts_offset + byte_len::<CompactTranscriptMeta>(transcripts.len()),
            kmers_offset,
        )?;
        write_pod_slice(&mut file, &kmers)?;
        write_padding(
            &mut file,
            kmers_offset + byte_len::<CompactKmerEntry>(kmers.len()),
            postings_offset,
        )?;
        write_pod_slice(&mut file, &postings)?;
        write_padding(
            &mut file,
            postings_offset + byte_len::<CompactPosting>(postings.len()),
            gene_names_offset,
        )?;
        file.write_all(&gene_names)?;
        write_padding(
            &mut file,
            gene_names_offset + gene_names.len() as u64,
            transcript_names_offset,
        )?;
        file.write_all(&transcript_names)?;
        write_padding(
            &mut file,
            transcript_names_offset + transcript_names.len() as u64,
            sequences_offset,
        )?;
        file.write_all(&sequences)?;
        Ok(())
    }

    pub fn load_json(path: impl AsRef<Path>) -> Result<Self, IndexIoError> {
        let bytes = std::fs::read(path)?;
        Ok(serde_json::from_slice(&bytes)?)
    }

    pub fn load_auto(path: impl AsRef<Path>) -> Result<LoadedIndex, IndexIoError> {
        LoadedIndex::load(path)
    }

    pub fn load_owned_auto(path: impl AsRef<Path>) -> Result<Self, IndexIoError> {
        let path = path.as_ref();
        let file = File::open(path)?;
        let mmap = unsafe { Mmap::map(&file)? };
        if mmap.starts_with(LEGACY_MAGIC) {
            return Self::from_legacy_mmap(&mmap);
        }
        if mmap.starts_with(&COMPACT_MAGIC) {
            return CompactTranscriptIndex::from_mmap(mmap)?.to_owned_index();
        }
        Ok(serde_json::from_slice(&mmap)?)
    }

    fn from_legacy_mmap(mmap: &[u8]) -> Result<Self, IndexIoError> {
        if mmap.len() < LEGACY_MAGIC.len() + 8 {
            return Err(IndexIoError::Truncated);
        }
        if &mmap[..LEGACY_MAGIC.len()] != LEGACY_MAGIC {
            return Err(IndexIoError::InvalidMagic);
        }
        let mut len_bytes = [0_u8; 8];
        len_bytes.copy_from_slice(&mmap[LEGACY_MAGIC.len()..LEGACY_MAGIC.len() + 8]);
        let payload_len = u64::from_le_bytes(len_bytes) as usize;
        let start = LEGACY_MAGIC.len() + 8;
        let end = start
            .checked_add(payload_len)
            .ok_or(IndexIoError::Truncated)?;
        if end > mmap.len() {
            return Err(IndexIoError::Truncated);
        }
        Ok(serde_json::from_slice(&mmap[start..end])?)
    }

    pub fn lookup(&self, kmer_code: u64) -> &[Posting] {
        match self
            .kmers
            .binary_search_by_key(&kmer_code, |entry| entry.kmer_code)
        {
            Ok(idx) => {
                let entry = &self.kmers[idx];
                let start = entry.postings_start as usize;
                let end = start + entry.postings_len as usize;
                &self.postings[start..end]
            }
            Err(_) => &[],
        }
    }
}

impl LoadedIndex {
    pub fn load(path: impl AsRef<Path>) -> Result<Self, IndexIoError> {
        let path = path.as_ref();
        let file = File::open(path)?;
        let mmap = unsafe { Mmap::map(&file)? };
        if mmap.starts_with(&COMPACT_MAGIC) {
            Ok(Self::Compact(CompactTranscriptIndex::from_mmap(mmap)?))
        } else if mmap.starts_with(LEGACY_MAGIC) {
            Ok(Self::Owned(TranscriptIndex::from_legacy_mmap(&mmap)?))
        } else {
            Ok(Self::Owned(serde_json::from_slice(&mmap)?))
        }
    }
}

impl CompactTranscriptIndex {
    fn from_mmap(mmap: Mmap) -> Result<Self, IndexIoError> {
        if mmap.len() < std::mem::size_of::<CompactHeader>() {
            return Err(IndexIoError::Truncated);
        }
        let header: CompactHeader =
            *bytemuck::from_bytes(&mmap[..std::mem::size_of::<CompactHeader>()]);
        if header.magic != COMPACT_MAGIC {
            return Err(IndexIoError::InvalidMagic);
        }
        if header.version != COMPACT_VERSION {
            return Err(IndexIoError::UnsupportedVersion(header.version));
        }
        let this = Self { mmap, header };
        this.validate()?;
        Ok(this)
    }

    fn validate(&self) -> Result<(), IndexIoError> {
        self.slice::<CompactGeneMeta>(self.header.genes_offset, self.header.num_genes)?;
        self.slice::<CompactTranscriptMeta>(
            self.header.transcripts_offset,
            self.header.num_transcripts,
        )?;
        self.slice::<CompactKmerEntry>(self.header.kmers_offset, self.header.num_kmers)?;
        self.slice::<CompactPosting>(self.header.postings_offset, self.header.num_postings)?;
        self.bytes(self.header.gene_names_offset, self.header.gene_names_len)?;
        self.bytes(
            self.header.transcript_names_offset,
            self.header.transcript_names_len,
        )?;
        self.bytes(self.header.sequences_offset, self.header.sequences_len)?;
        Ok(())
    }

    fn slice<T: Pod>(&self, offset: u64, count: u64) -> Result<&[T], IndexIoError> {
        if offset % std::mem::align_of::<T>() as u64 != 0 {
            return Err(IndexIoError::Misaligned);
        }
        let len = byte_len::<T>(count as usize);
        let bytes = self.bytes(offset, len)?;
        bytemuck::try_cast_slice(bytes).map_err(|_| IndexIoError::Misaligned)
    }

    fn bytes(&self, offset: u64, len: u64) -> Result<&[u8], IndexIoError> {
        let start = offset as usize;
        let end = start
            .checked_add(len as usize)
            .ok_or(IndexIoError::Truncated)?;
        self.mmap.get(start..end).ok_or(IndexIoError::Truncated)
    }

    fn genes(&self) -> &[CompactGeneMeta] {
        self.slice(self.header.genes_offset, self.header.num_genes)
            .expect("validated compact genes")
    }

    fn transcripts(&self) -> &[CompactTranscriptMeta] {
        self.slice(self.header.transcripts_offset, self.header.num_transcripts)
            .expect("validated compact transcripts")
    }

    fn kmers(&self) -> &[CompactKmerEntry] {
        self.slice(self.header.kmers_offset, self.header.num_kmers)
            .expect("validated compact kmers")
    }

    fn postings_slice(&self) -> &[CompactPosting] {
        self.slice(self.header.postings_offset, self.header.num_postings)
            .expect("validated compact postings")
    }

    fn gene_names(&self) -> &[u8] {
        self.bytes(self.header.gene_names_offset, self.header.gene_names_len)
            .expect("validated compact gene names")
    }

    fn transcript_names(&self) -> &[u8] {
        self.bytes(
            self.header.transcript_names_offset,
            self.header.transcript_names_len,
        )
        .expect("validated compact transcript names")
    }

    fn sequences(&self) -> &[u8] {
        self.bytes(self.header.sequences_offset, self.header.sequences_len)
            .expect("validated compact sequences")
    }

    fn entry(&self, kmer_code: u64) -> Option<CompactKmerEntry> {
        self.kmers()
            .binary_search_by_key(&kmer_code, |entry| entry.kmer_code)
            .ok()
            .map(|idx| self.kmers()[idx])
    }

    fn to_owned_index(&self) -> Result<TranscriptIndex, IndexIoError> {
        let genes = (0..self.num_genes())
            .map(|idx| self.gene_name(idx as GeneId).unwrap_or("").to_owned())
            .collect();
        let transcripts = self
            .transcripts()
            .iter()
            .map(|meta| TranscriptMeta {
                transcript_id: meta.transcript_id,
                gene_id: meta.gene_id,
                name: std::str::from_utf8(
                    &self.transcript_names()
                        [meta.name_start as usize..(meta.name_start + meta.name_len) as usize],
                )
                .unwrap_or("")
                .to_owned(),
                len: meta.len,
            })
            .collect();
        let transcript_sequences = self
            .transcripts()
            .iter()
            .map(|meta| {
                String::from_utf8_lossy(
                    &self.sequences()
                        [meta.seq_start as usize..(meta.seq_start + meta.len as u64) as usize],
                )
                .into_owned()
            })
            .collect();
        let kmers = self
            .kmers()
            .iter()
            .map(|entry| KmerEntry {
                kmer_code: entry.kmer_code,
                postings_start: entry.postings_start as u64,
                postings_len: entry.postings_len,
                freq_class: if entry.postings_len > self.max_kmer_frequency() {
                    1
                } else {
                    0
                },
            })
            .collect();
        let compact_transcripts = self.transcripts();
        let postings = self
            .postings_slice()
            .iter()
            .copied()
            .map(|posting| posting.to_posting(compact_transcripts))
            .collect();
        Ok(TranscriptIndex {
            format_version: self.header.version,
            k: self.k(),
            max_kmer_frequency: self.max_kmer_frequency(),
            genes,
            transcripts,
            transcript_sequences,
            kmers,
            postings,
        })
    }
}

impl IndexAccess for TranscriptIndex {
    fn k(&self) -> u8 {
        self.k
    }

    fn max_kmer_frequency(&self) -> u32 {
        self.max_kmer_frequency
    }

    fn num_genes(&self) -> usize {
        self.genes.len()
    }

    fn gene_name(&self, gene_id: GeneId) -> Option<&str> {
        self.genes.get(gene_id as usize).map(String::as_str)
    }

    fn transcript_gene_id(&self, transcript_id: TranscriptId) -> Option<GeneId> {
        self.transcripts
            .get(transcript_id as usize)
            .map(|meta| meta.gene_id)
    }

    fn transcript_seq(&self, transcript_id: TranscriptId) -> Option<&[u8]> {
        self.transcript_sequences
            .get(transcript_id as usize)
            .map(|seq| seq.as_bytes())
    }

    fn posting_count(&self, kmer_code: u64) -> usize {
        self.kmers
            .binary_search_by_key(&kmer_code, |entry| entry.kmer_code)
            .ok()
            .map_or(0, |idx| self.kmers[idx].postings_len as usize)
    }

    fn seed_postings(&self, kmer_code: u64) -> SeedPostingIter<'_> {
        match self
            .kmers
            .binary_search_by_key(&kmer_code, |entry| entry.kmer_code)
        {
            Ok(idx) => {
                let entry = &self.kmers[idx];
                let start = entry.postings_start as usize;
                let end = start + entry.postings_len as usize;
                SeedPostingIter(SeedPostingIterInner::Owned(
                    self.postings[start..end].iter().copied(),
                ))
            }
            Err(_) => SeedPostingIter(SeedPostingIterInner::Empty),
        }
    }

    fn postings(&self, kmer_code: u64) -> PostingIter<'_> {
        match self
            .kmers
            .binary_search_by_key(&kmer_code, |entry| entry.kmer_code)
        {
            Ok(idx) => {
                let entry = &self.kmers[idx];
                let start = entry.postings_start as usize;
                let end = start + entry.postings_len as usize;
                PostingIter(PostingIterInner::Owned(
                    self.postings[start..end].iter().copied(),
                ))
            }
            Err(_) => PostingIter(PostingIterInner::Empty),
        }
    }

    fn stats(&self) -> IndexStats {
        IndexStats {
            num_transcripts: self.transcripts.len(),
            num_genes: self.genes.len(),
            num_distinct_kmers: self.kmers.len(),
            num_postings: self.postings.len(),
            k: self.k,
            max_postings_per_kmer: self
                .kmers
                .iter()
                .map(|entry| entry.postings_len)
                .max()
                .unwrap_or(0),
            high_frequency_kmers: self
                .kmers
                .iter()
                .filter(|entry| entry.postings_len > self.max_kmer_frequency)
                .count(),
        }
    }
}

impl IndexAccess for CompactTranscriptIndex {
    fn k(&self) -> u8 {
        self.header.k as u8
    }

    fn max_kmer_frequency(&self) -> u32 {
        self.header.max_kmer_frequency
    }

    fn num_genes(&self) -> usize {
        self.header.num_genes as usize
    }

    fn gene_name(&self, gene_id: GeneId) -> Option<&str> {
        let meta = self.genes().get(gene_id as usize)?;
        std::str::from_utf8(
            &self.gene_names()
                [meta.name_start as usize..(meta.name_start + meta.name_len) as usize],
        )
        .ok()
    }

    fn transcript_gene_id(&self, transcript_id: TranscriptId) -> Option<GeneId> {
        self.transcripts()
            .get(transcript_id as usize)
            .map(|meta| meta.gene_id)
    }

    fn transcript_seq(&self, transcript_id: TranscriptId) -> Option<&[u8]> {
        let meta = self.transcripts().get(transcript_id as usize)?;
        let start = meta.seq_start as usize;
        let end = start + meta.len as usize;
        self.sequences().get(start..end)
    }

    fn posting_count(&self, kmer_code: u64) -> usize {
        self.entry(kmer_code)
            .map_or(0, |entry| entry.postings_len as usize)
    }

    fn seed_postings(&self, kmer_code: u64) -> SeedPostingIter<'_> {
        match self.entry(kmer_code) {
            Some(entry) => {
                let start = entry.postings_start as usize;
                let end = start + entry.postings_len as usize;
                SeedPostingIter(SeedPostingIterInner::Compact(
                    self.postings_slice()[start..end].iter().copied(),
                ))
            }
            None => SeedPostingIter(SeedPostingIterInner::Empty),
        }
    }

    fn postings(&self, kmer_code: u64) -> PostingIter<'_> {
        match self.entry(kmer_code) {
            Some(entry) => {
                let start = entry.postings_start as usize;
                let end = start + entry.postings_len as usize;
                PostingIter(PostingIterInner::Compact {
                    postings: self.postings_slice()[start..end].iter().copied(),
                    transcripts: self.transcripts(),
                })
            }
            None => PostingIter(PostingIterInner::Empty),
        }
    }

    fn stats(&self) -> IndexStats {
        IndexStats {
            num_transcripts: self.header.num_transcripts as usize,
            num_genes: self.header.num_genes as usize,
            num_distinct_kmers: self.header.num_kmers as usize,
            num_postings: self.header.num_postings as usize,
            k: self.k(),
            max_postings_per_kmer: self
                .kmers()
                .iter()
                .map(|entry| entry.postings_len)
                .max()
                .unwrap_or(0),
            high_frequency_kmers: self
                .kmers()
                .iter()
                .filter(|entry| entry.postings_len > self.max_kmer_frequency())
                .count(),
        }
    }
}

impl IndexAccess for LoadedIndex {
    fn k(&self) -> u8 {
        match self {
            Self::Owned(index) => index.k(),
            Self::Compact(index) => index.k(),
        }
    }

    fn max_kmer_frequency(&self) -> u32 {
        match self {
            Self::Owned(index) => index.max_kmer_frequency(),
            Self::Compact(index) => index.max_kmer_frequency(),
        }
    }

    fn num_genes(&self) -> usize {
        match self {
            Self::Owned(index) => index.num_genes(),
            Self::Compact(index) => index.num_genes(),
        }
    }

    fn gene_name(&self, gene_id: GeneId) -> Option<&str> {
        match self {
            Self::Owned(index) => index.gene_name(gene_id),
            Self::Compact(index) => index.gene_name(gene_id),
        }
    }

    fn transcript_gene_id(&self, transcript_id: TranscriptId) -> Option<GeneId> {
        match self {
            Self::Owned(index) => index.transcript_gene_id(transcript_id),
            Self::Compact(index) => index.transcript_gene_id(transcript_id),
        }
    }

    fn transcript_seq(&self, transcript_id: TranscriptId) -> Option<&[u8]> {
        match self {
            Self::Owned(index) => index.transcript_seq(transcript_id),
            Self::Compact(index) => index.transcript_seq(transcript_id),
        }
    }

    fn posting_count(&self, kmer_code: u64) -> usize {
        match self {
            Self::Owned(index) => index.posting_count(kmer_code),
            Self::Compact(index) => index.posting_count(kmer_code),
        }
    }

    fn seed_postings(&self, kmer_code: u64) -> SeedPostingIter<'_> {
        match self {
            Self::Owned(index) => index.seed_postings(kmer_code),
            Self::Compact(index) => index.seed_postings(kmer_code),
        }
    }

    fn postings(&self, kmer_code: u64) -> PostingIter<'_> {
        match self {
            Self::Owned(index) => index.postings(kmer_code),
            Self::Compact(index) => index.postings(kmer_code),
        }
    }

    fn stats(&self) -> IndexStats {
        match self {
            Self::Owned(index) => IndexAccess::stats(index),
            Self::Compact(index) => IndexAccess::stats(index),
        }
    }
}

fn write_pod_slice<T: Pod>(file: &mut File, values: &[T]) -> Result<(), IndexIoError> {
    file.write_all(bytemuck::cast_slice(values))?;
    Ok(())
}

fn write_padding(file: &mut File, from: u64, to: u64) -> Result<(), IndexIoError> {
    if to > from {
        file.write_all(&vec![0_u8; (to - from) as usize])?;
    }
    Ok(())
}

fn align_up(value: u64, align: u64) -> u64 {
    (value + align - 1) / align * align
}

fn byte_len<T>(count: usize) -> u64 {
    (count * std::mem::size_of::<T>()) as u64
}

#[cfg(test)]
mod tests {
    use super::*;

    fn tiny_index() -> TranscriptIndex {
        TranscriptIndex {
            format_version: 1,
            k: 3,
            max_kmer_frequency: 256,
            genes: vec!["GENE".to_owned()],
            transcripts: vec![TranscriptMeta {
                transcript_id: 0,
                gene_id: 0,
                name: "tx".to_owned(),
                len: 4,
            }],
            transcript_sequences: vec!["ACGT".to_owned()],
            kmers: vec![KmerEntry {
                kmer_code: 6,
                postings_start: 0,
                postings_len: 1,
                freq_class: 0,
            }],
            postings: vec![Posting {
                transcript_id: 0,
                gene_id: 0,
                pos: 1,
                strand: 0,
            }],
        }
    }

    #[test]
    fn compact_mmap_roundtrip() {
        let index = tiny_index();
        let tmp = tempfile::NamedTempFile::new().unwrap();
        index.save_auto(tmp.path()).unwrap();
        let loaded = LoadedIndex::load(tmp.path()).unwrap();
        assert_eq!(loaded.k(), 3);
        assert_eq!(loaded.gene_name(0), Some("GENE"));
        assert_eq!(loaded.transcript_seq(0), Some(&b"ACGT"[..]));
        assert_eq!(loaded.posting_count(6), 1);
        assert_eq!(loaded.postings(6).collect::<Vec<_>>(), index.postings);
    }

    #[test]
    fn legacy_binary_mmap_roundtrip_to_owned() {
        let index = tiny_index();
        let tmp = tempfile::NamedTempFile::new().unwrap();
        let payload = serde_json::to_vec(&index).unwrap();
        let mut file = File::create(tmp.path()).unwrap();
        file.write_all(LEGACY_MAGIC).unwrap();
        file.write_all(&(payload.len() as u64).to_le_bytes())
            .unwrap();
        file.write_all(&payload).unwrap();
        assert_eq!(TranscriptIndex::load_owned_auto(tmp.path()).unwrap(), index);
    }
}

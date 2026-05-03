# CacheBatch RNA-seq Mapper: Codex-Ready Project Brief

**Working title:** `cachebatch`  
**Date:** 2026-05-03  
**Primary language:** Rust  
**Initial target:** 10x-style short-read single-cell RNA-seq gene-expression preprocessing and quantification  
**MVP goal:** demonstrate that cache-aware read batching + candidate-locus regrouping + SIMD filtering can reduce random index probes and candidate scoring work compared with a read-at-a-time baseline.

---

## 0. Prompt for Codex Agent

You are building a Rust prototype called `cachebatch`: a cache-aware, batch-oriented RNA-seq/scRNA-seq mapper/quantifier. The goal is not to clone STAR, Cell Ranger, kallisto, or Salmon. The goal is to prototype a new dataflow architecture:

```text
FASTQ reads
  -> compact 2-bit read records
  -> minimizer/syncmer/LSH sketch records
  -> sketch-sorted read buckets
  -> seed lookup
  -> candidate-locus buckets
  -> scalar + SIMD filters
  -> gene/transcript compatibility assignments
  -> benchmark reports
```

Start by creating a Cargo workspace with a core library and CLI. Implement a correctness-first scalar path, then add a `pulp` SIMD path for the regular scoring/filtering kernels. Keep all SIMD behind a trait so scalar and SIMD outputs can be compared byte-for-byte. The first milestone should compile, run unit tests, build a tiny transcript index from FASTA, map synthetic reads against it, and report benchmark metrics.

Do **not** begin with a full spliced genomic aligner. Do **not** attempt BAM/CRAM output in the first pass. Do **not** overfit to STAR semantics initially. Build a small but measurable prototype that tests the central hypothesis: batching reads and regrouping by candidate locus improves cache locality and throughput.

---

## 1. Scientific / Engineering Hypothesis

Current RNA-seq alignment and pseudoalignment tools often process reads in a read-centric way:

```text
read_i -> seed lookup -> candidate locations -> score/align -> output
read_j -> seed lookup -> candidate locations -> score/align -> output
...
```

This causes repeated random access into large index structures and candidate lists. For human-scale transcriptome/genome indexes, the global index is far larger than L1/L2 cache. The project hypothesis is:

> The biggest early win may come from cache-aware scheduling: group reads by compact sequence sketches, then regroup by candidate locus so the same reference/index chunks are reused many times before eviction. SIMD should then be used inside these regular, cache-resident scoring blocks.

The intended breakthrough is not merely “add SIMD to dynamic programming.” The intended breakthrough is changing the dataflow so SIMD and cache locality become natural.

---

## 2. Initial Scope

### 2.1 MVP Input

Support these inputs first:

```text
cachebatch index \
  --transcripts transcripts.fa \
  --k 21 \
  --out transcript.cbidx

cachebatch map \
  --index transcript.cbidx \
  --r1 sample_R1.fastq.gz \
  --r2 sample_R2.fastq.gz \
  --chemistry tenx-3p-v3 \
  --out assignments.tsv
```

For MVP, `transcripts.fa` can have headers that encode gene IDs, for example:

```text
>ENST000001|gene=GENE_A
ACGT...
```

GTF parsing can be added later.

### 2.2 MVP Output

Start with assignment output, not a production count matrix:

```tsv
read_id	cell_barcode	umi	assignment_type	gene_id	transcript_id	candidate_count	score	flags
```

Assignment types:

```text
unique_gene
ambiguous_gene
ambiguous_transcript_same_gene
unmapped
low_complexity
low_quality
```

Later add:

```text
cachebatch count -> sparse matrix / Matrix Market / BUS-like output
```

### 2.3 Initial Non-goals

Do not implement these in the first milestone:

- full splice-aware genomic alignment;
- novel splice-junction discovery;
- BAM/CRAM output;
- allele-specific assignment;
- full Cell Ranger compatibility;
- UMI error correction beyond exact deduplication;
- production-grade GTF/GFF handling;
- AVX-512 hand assembly.

---

## 3. Recommended Rust Stack

Use Rust as the main implementation language.

### 3.1 Core Dependencies

Suggested dependencies:

```toml
[dependencies]
anyhow = "1"
thiserror = "2"
clap = { version = "4", features = ["derive"] }
rayon = "1"
memmap2 = "0.9"
pulp = "0.22"
needletail = "0.7"
rustc-hash = "2"
smallvec = "1"
bytemuck = { version = "1", features = ["derive"] }
serde = { version = "1", features = ["derive"] }
serde_json = "1"

[dev-dependencies]
criterion = "0.8"
proptest = "1"
tempfile = "3"
```

Optional later:

```toml
bio = "3"          # useful reference algorithms / FM-index / alignment utilities
rten-simd = "0.24" # alternative SIMD abstraction, broader cross-architecture story
niffler = "3"      # compressed file abstraction, if needed
```

Rationale:

- `pulp` provides a safe abstraction over SIMD and runtime dispatch to vectorized implementations.
- `rayon` provides data-parallel iterators, `par_sort`, and custom thread-pool support.
- `memmap2` enables immutable file-backed index slices.
- `needletail` is a fast FASTA/FASTQ parser.
- `criterion` is the initial microbenchmark harness.
- `bio` is useful as a reference / utilities library but should not dictate the project architecture.

Relevant current docs:

- `pulp`: https://docs.rs/pulp/latest/pulp/
- `rayon`: https://docs.rs/rayon/latest/rayon/
- `memmap2`: https://docs.rs/memmap2/latest/memmap2/
- `needletail`: https://docs.rs/needletail/latest/needletail/
- `rust-bio`: https://docs.rs/bio/latest/bio/
- `criterion`: https://docs.rs/criterion/latest/criterion/
- `rten-simd`: https://docs.rs/rten-simd/latest/rten_simd/

---

## 4. Repository Layout

Create this workspace:

```text
cachebatch/
  Cargo.toml
  README.md
  crates/
    cachebatch-core/
      Cargo.toml
      src/
        lib.rs
        dna.rs
        fastq.rs
        chemistry.rs
        sketch.rs
        minimizer.rs
        index.rs
        index_build.rs
        lookup.rs
        batch.rs
        candidate.rs
        score.rs
        score_scalar.rs
        score_simd.rs
        assign.rs
        metrics.rs
        simulate.rs
    cachebatch-cli/
      Cargo.toml
      src/
        main.rs
        cmd_index.rs
        cmd_map.rs
        cmd_simulate.rs
        cmd_bench_report.rs
  benches/
    encode.rs
    sketch.rs
    lookup.rs
    bucket.rs
    score.rs
    end_to_end_synthetic.rs
  tests/
    tiny_transcriptome.fa
    tiny_reads_R1.fastq
    tiny_reads_R2.fastq
    expected_assignments.tsv
  scripts/
    perf_map.sh
    run_public_benchmarks.sh
    compare_assignments.py
  docs/
    architecture.md
    benchmark_plan.md
```

Root `Cargo.toml`:

```toml
[workspace]
members = [
    "crates/cachebatch-core",
    "crates/cachebatch-cli",
]
resolver = "2"

[profile.release]
debug = true
lto = "thin"
codegen-units = 1

[profile.bench]
debug = true
```

---

## 5. Architecture

### 5.1 Dataflow Overview

```text
             +------------------+
FASTQ R1 --->| CB/UMI parser    |----+
             +------------------+    |
                                     v
             +------------------+  CompactRead
FASTQ R2 --->| 2-bit encoder    |----+
             +------------------+
                                     |
                                     v
             +------------------+  SketchRecord
             | sketch/minimizer |----+
             +------------------+    |
                                     v
             +------------------+  SketchBucket
             | radix partition  |----+
             +------------------+    |
                                     v
             +------------------+  CandidateHit
             | seed lookup      |----+
             +------------------+    |
                                     v
             +------------------+  CandidateLocusBucket
             | locus regrouping |----+
             +------------------+    |
                                     v
             +------------------+  ScoredCandidate
             | scalar/SIMD      |----+
             | filters          |
             +------------------+
                                     |
                                     v
             +------------------+  Assignment
             | compatibility    |-----> TSV / counts later
             | assignment       |
             +------------------+
```

### 5.2 Key Principle

Read similarity is only a means to an end. The real objective is cache coherence.

Do not spend too much time finding exact nearest-neighbor reads. Instead:

1. create cheap read sketch keys;
2. radix-sort or bucket by those keys;
3. perform a few seed lookups;
4. regroup work by candidate reference/transcript locus;
5. run SIMD filters over contiguous SoA buffers.

---

## 6. Core Data Structures

Use compact integer IDs and contiguous arrays. Prefer structure-of-arrays for hot SIMD paths.

### 6.1 DNA Encoding

Use 2-bit encoding:

```text
A = 00
C = 01
G = 10
T = 11
N = invalid/masked bit
```

Implement:

```rust
pub struct EncodedSeq {
    pub words: Vec<u64>,       // 32 bases per u64
    pub len: u32,
    pub n_mask_words: Vec<u64> // 1 bit/base mask, or empty if no Ns
}
```

Functions:

```rust
encode_acgt(seq: &[u8]) -> EncodedSeq
base_at(encoded: &EncodedSeq, pos: usize) -> Option<Base>
reverse_complement(encoded: &EncodedSeq) -> EncodedSeq
iter_kmers_2bit(encoded: &EncodedSeq, k: u8) -> impl Iterator<Item = Kmer>
hamming_2bit(a: &EncodedSeq, b: &EncodedSeq) -> u32
```

Tests:

- encode/decode round trip;
- reverse-complement twice returns original;
- k-mer iteration agrees with string implementation;
- Ns are skipped or masked correctly.

### 6.2 Read Records

```rust
pub type ReadId = u64;
pub type GeneId = u32;
pub type TranscriptId = u32;

pub struct CompactRead {
    pub read_id: ReadId,
    pub seq_offset: u64,       // offset into packed sequence arena
    pub len: u16,
    pub qual_offset: u64,      // optional; can be 0 if qualities not stored
    pub cb: u64,               // packed 16 bp barcode if known
    pub umi: u64,              // packed 10-12 bp UMI if known
    pub flags: u16,
}
```

For 10x 3′ v3:

```text
R1: 28 bp = 16 bp cell barcode + 12 bp UMI
R2: transcript/cDNA read, commonly ~90–100 bp depending on sequencing setup
```

Make `chemistry.rs` responsible for parsing R1.

### 6.3 Sketch Records

```rust
pub struct SketchRecord {
    pub read_id: ReadId,
    pub key_primary: u64,
    pub key_secondary: u64,
    pub rare_anchor: u64,
    pub gc_bin: u8,
    pub len: u16,
    pub flags: u16,
}
```

Potential keys:

```text
key_primary   = hash(minimizer_0, minimizer_1, len_bin, low_complexity_flag)
key_secondary = simhash/sketch hash over selected k-mers
rare_anchor   = lowest-frequency selected seed according to transcript index
```

Do not implement a dense k-mer vector in the MVP. Use sketches and radix-sortable keys.

### 6.4 Transcript Index

MVP index format:

```rust
pub struct TranscriptMeta {
    pub transcript_id: TranscriptId,
    pub gene_id: GeneId,
    pub name_offset: u32,
    pub seq_offset: u64,
    pub len: u32,
}

pub struct Posting {
    pub transcript_id: TranscriptId,
    pub gene_id: GeneId,
    pub pos: u32,
    pub strand: u8,
}

pub struct KmerEntry {
    pub kmer_hash_or_code: u64,
    pub postings_start: u64,
    pub postings_len: u32,
    pub freq_class: u8,
}
```

Build an index where:

- transcript sequences are concatenated and 2-bit encoded;
- observed k-mers point to contiguous posting-list slices;
- k-mer entries are sorted by k-mer code/hash;
- posting lists are sorted by `(transcript_id, pos)`;
- high-frequency k-mers are tagged or excluded from sketch keys.

For the MVP, a simple `HashMap<u64, Vec<Posting>>` builder is acceptable, but write the serialized index as sorted arrays so lookup can be binary search or later a compact hash table.

### 6.5 Candidate Hits and Buckets

```rust
pub struct CandidateHit {
    pub read_id: ReadId,
    pub transcript_id: TranscriptId,
    pub gene_id: GeneId,
    pub pos: u32,
    pub strand: u8,
    pub seed_count: u16,
    pub seed_score: u16,
}

pub struct CandidateBucketKey {
    pub transcript_id: TranscriptId,
    pub pos_bin: u32,  // e.g. pos / 64 or pos / 128
    pub strand: u8,
}
```

Pipeline:

```text
SketchBucket -> seed lookups -> CandidateHit[] -> sort by CandidateBucketKey -> CandidateLocusBucket[]
```

This candidate-locus regrouping is the main cache optimization.

---

## 7. Batching Strategy

### 7.1 Do Not Use Dense Vector Similarity First

A dense vector over all possible k-mers is too large and sparse. A 100 bp read has only ~70–80 useful k-mers for k around 21–31, while the space of possible k-mers is enormous.

Start with cheap sketching:

```text
selected rare k-mers
minimizers
syncmers, optional later
SimHash / MinHash, optional later
quality-filtered anchors
low-complexity mask
poly-A/T flag
GC bin
```

### 7.2 Initial Bucket Formation

Algorithm:

```text
for each read:
    encode R2
    extract valid k-mers
    drop k-mers containing N or low-quality bases
    drop/penalize high-frequency k-mers according to index
    choose rare anchor(s)
    compute minimizer/simhash keys
    emit SketchRecord

sort SketchRecord[] by (key_primary, key_secondary)
split into SketchBuckets of target size, e.g. 256–4096 reads
```

Use `rayon` parallel sorting where easy.

### 7.3 Candidate-Locus Regrouping

Within each sketch bucket:

```text
for read in bucket:
    probe a limited number of informative seeds
    produce candidate hits

sort candidate hits by (transcript_id, pos_bin, strand)
for each candidate-locus bucket:
    load reference/transcript window once
    score many reads against that window
```

This matters more than the first read-similarity sort. The CPU cache benefits when many operations touch the same transcript/reference region.

---

## 8. Scoring and SIMD

### 8.1 Scoring Trait

Define a trait so scalar and SIMD implementations are interchangeable:

```rust
pub trait CandidateScorer: Send + Sync {
    fn score_bucket(
        &self,
        reads: &ReadStore,
        index: &TranscriptIndex,
        bucket: &CandidateLocusBucket,
        out: &mut Vec<ScoredCandidate>,
    );
}
```

Implement:

```rust
pub struct ScalarScorer;
pub struct PulpScorer;
```

The scalar scorer is the correctness oracle. The SIMD scorer must match scalar output for deterministic tests.

### 8.2 SIMD Kernels to Implement First

Start with regular filters, not full DP:

1. **Packed 2-bit Hamming filter**  
   Compare read window and transcript window. Allow a configurable mismatch threshold.

2. **Seed-overlap score**  
   Count how many selected seeds support the same candidate transcript/locus.

3. **Quality-aware mask filter**  
   Ignore or downweight bases below quality threshold.

4. **Small banded edit-distance filter**  
   Later. Start scalar, then vectorize if profiling says it matters.

### 8.3 `pulp` Usage Pattern

Dispatch once per worker or per large batch, not once per read.

```rust
use pulp::{Arch, WithSimd, Simd};

pub struct HammingKernel<'a> {
    // references to SoA buffers
}

impl<'a> WithSimd for HammingKernel<'a> {
    type Output = ();

    #[inline(always)]
    fn with_simd<S: Simd>(self, simd: S) -> Self::Output {
        // vectorized loop here
        // fall back to scalar tail for remainder
    }
}

pub fn run_hamming_kernel(kernel: HammingKernel<'_>) {
    let arch = Arch::new();
    arch.dispatch(kernel);
}
```

If a kernel needs specialized instructions not exposed by `pulp`, keep that in a tiny `unsafe` `std::arch` module and preserve the scalar fallback.

### 8.4 Data Layout for SIMD

Avoid array-of-structs in hot loops.

Bad:

```rust
Vec<ReadCandidatePair>
```

Better:

```rust
pub struct ScoreBlock {
    pub read_ids: Vec<ReadId>,
    pub read_word0: Vec<u64>,
    pub read_word1: Vec<u64>,
    pub ref_word0: Vec<u64>,
    pub ref_word1: Vec<u64>,
    pub qual_mask0: Vec<u64>,
    pub pos: Vec<u32>,
}
```

A candidate-locus bucket should be transformed into a SoA `ScoreBlock` before SIMD filtering.

---

## 9. Assignment Logic

For each read, collect scored candidates and emit one assignment.

Initial rules:

```text
if no candidate passes threshold:
    unmapped
else if exactly one gene among passing candidates:
    unique_gene
else if multiple transcripts but same gene:
    ambiguous_transcript_same_gene
else:
    ambiguous_gene
```

Later add posterior scoring:

```text
P(gene | seeds, candidate scores, read quality, expression prior)
```

For now, keep priors off by default to avoid bias. Implement the assignment module in a way that can later support probabilistic early stopping.

---

## 10. Cache / Memory Strategy

### 10.1 Principle

The whole index cannot fit in L1/L2. The goal is to make the active working set for one candidate-locus bucket fit in L1/L2.

A good hot block should contain:

```text
- one transcript/reference window or a small set of adjacent windows;
- 16–4096 reads depending on stage;
- candidate metadata arrays;
- packed sequence words;
- score output buffers.
```

### 10.2 Implementation Rules

- Use contiguous arrays for postings and candidate hits.
- Avoid pointer-heavy structures in hot paths.
- Use `u32` IDs where possible.
- Sort by keys to improve locality.
- Limit allocations inside per-bucket loops.
- Use reusable thread-local scratch buffers.
- Exclude or downweight very frequent k-mers.
- Measure cache misses rather than guessing.

### 10.3 Metrics to Track

At minimum:

```text
reads/sec
bases/sec
candidate hits/read
candidate loci/read
seed lookups/read
scored candidates/read
unique_gene rate
ambiguous_gene rate
unmapped rate
peak RSS
wall-clock time
CPU time
```

With Linux `perf`:

```text
cycles
instructions
cache-references
cache-misses
L1-dcache-loads
L1-dcache-load-misses
LLC-loads
LLC-load-misses
branch-misses
```

Example script:

```bash
#!/usr/bin/env bash
set -euo pipefail
perf stat \
  -e cycles,instructions,cache-references,cache-misses,branches,branch-misses,L1-dcache-loads,L1-dcache-load-misses,LLC-loads,LLC-load-misses \
  "$@"
```

---

## 11. CLI Design

### 11.1 `index`

```bash
cachebatch index \
  --transcripts transcripts.fa \
  --k 21 \
  --max-kmer-frequency 256 \
  --out transcript.cbidx
```

Responsibilities:

- parse transcript FASTA;
- derive transcript and gene IDs;
- encode sequences;
- build k-mer postings;
- tag high-frequency k-mers;
- serialize index metadata + arrays.

### 11.2 `simulate`

```bash
cachebatch simulate \
  --transcripts transcripts.fa \
  --num-reads 100000 \
  --read-len 91 \
  --error-rate 0.005 \
  --chemistry tenx-3p-v3 \
  --out-prefix sim
```

Responsibilities:

- generate R1/R2 FASTQ;
- generate ground-truth TSV;
- support exact reads, substitution errors, Ns, and multimapping stress tests.

### 11.3 `map`

```bash
cachebatch map \
  --index transcript.cbidx \
  --r1 sim_R1.fastq.gz \
  --r2 sim_R2.fastq.gz \
  --chemistry tenx-3p-v3 \
  --batch-size 2048 \
  --score-mode scalar \
  --out assignments.tsv
```

Also support:

```bash
--score-mode pulp
--emit-metrics metrics.json
--max-seeds-per-read 8
--max-postings-per-seed 256
--candidate-bin-size 64
--threads N
```

### 11.4 `bench-report`

```bash
cachebatch bench-report \
  --assignments assignments.tsv \
  --truth truth.tsv \
  --metrics metrics.json
```

Responsibilities:

- compare assignments to truth;
- report mapping accuracy;
- report ambiguity rates;
- summarize speed and cache metrics if available.

---

## 12. Implementation Milestones

### Milestone 0 — Scaffold and Tiny Tests

Deliverables:

- Cargo workspace compiles.
- CLI has `index`, `simulate`, `map` subcommands stubbed.
- `dna.rs` supports 2-bit encoding, decoding, reverse complement, k-mer iteration.
- Unit/property tests for DNA encoding pass.

Acceptance criteria:

```bash
cargo test
cargo run -p cachebatch-cli -- --help
```

### Milestone 1 — Transcript Index Builder

Deliverables:

- parse FASTA;
- derive `gene_id` from header when available;
- build `k`-mer postings;
- ignore k-mers containing N;
- serialize/deserialze index.

Acceptance criteria:

```bash
cachebatch index --transcripts tests/tiny_transcriptome.fa --k 15 --out /tmp/tiny.cbidx
cachebatch inspect-index --index /tmp/tiny.cbidx
```

`inspect-index` should report:

```text
num_transcripts
num_genes
num_distinct_kmers
num_postings
k
max/posting frequency histogram
```

### Milestone 2 — Read Parser, Chemistry Parser, Simulator

Deliverables:

- read paired R1/R2 FASTQ;
- parse 10x 3′ v3 barcode/UMI;
- simulate reads with ground truth;
- exact no-error simulated reads can be mapped by simple seed lookup.

Acceptance criteria:

```bash
cachebatch simulate --transcripts tests/tiny_transcriptome.fa --num-reads 1000 --read-len 91 --out-prefix /tmp/sim
```

### Milestone 3 — Sketching and Read Bucketing

Deliverables:

- extract selected seeds per read;
- compute primary/secondary sketch keys;
- radix/sort into sketch buckets;
- output metrics on bucket sizes and seed frequency.

Acceptance criteria:

- deterministic output for fixed seed;
- no excessive allocations inside per-read loop;
- benchmark reads/sec for sketching.

### Milestone 4 — Candidate Generation and Locus Regrouping

Deliverables:

- seed lookup returns candidate hits;
- candidate hits sorted by `(transcript_id, pos_bin, strand)`;
- candidate-locus buckets created;
- metrics collected:

```text
candidate hits/read
candidate buckets/read
mean bucket size
median bucket size
max bucket size
postings skipped due to frequency
```

Acceptance criteria:

- candidate generator recall ≥ 99.5% on synthetic no-error reads;
- candidate generator recall ≥ 98% on synthetic reads with 0.5% substitution error, assuming enough seeds.

### Milestone 5 — Scalar Scoring and Assignment

Deliverables:

- scalar Hamming filter against candidate windows;
- assignment rules implemented;
- assignment TSV output;
- synthetic benchmark report.

Acceptance criteria:

- no-error synthetic unique reads: ≥ 99% correct gene assignment;
- reads from identical/shared transcripts are marked ambiguous rather than falsely unique;
- output deterministic across runs.

### Milestone 6 — `pulp` SIMD Scoring

Deliverables:

- `PulpScorer` implements the same trait as `ScalarScorer`;
- differential tests compare scalar vs SIMD outputs;
- Criterion benchmark compares scalar vs SIMD for scoring blocks.

Acceptance criteria:

```bash
cargo test simd_matches_scalar
cargo bench --bench score
```

Report:

```text
scalar ns/read-candidate
pulp ns/read-candidate
speedup
```

### Milestone 7 — Cache-Aware Benchmarking

Deliverables:

- read-at-a-time baseline mode;
- sketch-bucket mode;
- sketch-bucket + candidate-locus regroup mode;
- perf wrapper script;
- benchmark report comparing cache misses and throughput.

Acceptance criteria:

The benchmark report must show:

```text
mode
reads/sec
candidate hits/read
scored candidates/read
L1 miss rate if available
LLC miss rate if available
peak RSS
```

### Milestone 8 — Public Dataset Smoke Test

Deliverables:

- run on a public 10x PBMC subset;
- compare against at least one baseline result from STARsolo, kallisto|bustools, or Cell Ranger when available;
- report gene-count concordance.

Suggested public datasets:

- 10x Genomics 10k PBMCs v3 chemistry dataset.
- 10x Genomics 5k Human PBMCs, 3′ v3.1 dataset.

---

## 13. Testing Plan

### 13.1 Unit Tests

Required tests:

```text
dna::encode_decode_roundtrip
dna::reverse_complement_involution
dna::kmer_iterator_matches_string_reference
chemistry::parse_tenx_3p_v3_r1
index::tiny_index_contains_expected_kmers
lookup::known_kmer_returns_expected_postings
sketch::deterministic_keys_for_fixed_read
candidate::candidate_locus_sorting_is_stable
score_scalar::hamming_exact_match_zero
score_scalar::hamming_one_mismatch
assign::unique_gene_assignment
assign::ambiguous_gene_assignment
```

### 13.2 Property Tests

Use `proptest`:

- random A/C/G/T strings encode/decode correctly;
- reverse-complement twice equals original;
- scalar hamming equals naive string hamming;
- random transcript/read simulation maps exact reads to a candidate containing the true origin;
- sorted candidate buckets preserve all candidate hits.

### 13.3 Differential Tests

Critical:

```text
PulpScorer == ScalarScorer
```

For random score blocks:

```rust
assert_eq!(scalar_scores, pulp_scores);
```

Run these tests under different CPU feature availability when possible. At minimum, ensure the scalar fallback is always available.

### 13.4 Golden Tests

Create a tiny transcriptome:

```text
GENE_A transcript A1: unique region + shared region
GENE_A transcript A2: different unique region + shared region
GENE_B transcript B1: unique region
GENE_C transcript C1: intentionally similar/paralogous region
```

Create reads that should produce:

```text
unique_gene
ambiguous_transcript_same_gene
ambiguous_gene
unmapped
low_complexity
```

Compare exact assignment TSV to `tests/expected_assignments.tsv`.

### 13.5 Regression Tests

Add a regression test for every bug found. Bugs in this domain are often silent correctness bugs, not crashes.

---

## 14. Benchmark Plan

### 14.1 Microbenchmarks with Criterion

Benchmarks:

1. `encode.rs`
   - FASTQ sequence to 2-bit encoding throughput.
   - Report bases/sec.

2. `sketch.rs`
   - k-mer extraction + minimizer/sketch computation.
   - Report reads/sec and ns/read.

3. `lookup.rs`
   - seed lookup random order vs sketch-bucketed order.
   - Report lookups/sec and cache behavior via external `perf`.

4. `bucket.rs`
   - sort/radix partition sketch records.
   - sort candidate hits into candidate-locus buckets.

5. `score.rs`
   - scalar vs `pulp` Hamming filter.
   - batch sizes: 16, 64, 256, 1024, 4096.
   - read lengths: 50, 91, 100, 150.

6. `end_to_end_synthetic.rs`
   - full pipeline on simulated reads.

### 14.2 End-to-End Modes to Compare

Implement three modes:

```text
mode A: read-at-a-time baseline
mode B: sketch-bucketed reads, no locus regroup
mode C: sketch-bucketed + candidate-locus regroup
mode D: mode C + pulp scoring
```

Required report:

```text
mode
num_reads
wall_seconds
reads_per_second
bases_per_second
peak_rss_mb
seed_lookups_per_read
candidate_hits_per_read
scored_candidates_per_read
unique_gene_rate
ambiguous_gene_rate
unmapped_rate
correct_gene_assignment_rate on synthetic data
L1_miss_rate if perf available
LLC_miss_rate if perf available
```

### 14.3 Synthetic Accuracy Benchmarks

Generate synthetic datasets:

```text
S1: exact reads, no errors
S2: substitution error rate 0.1%
S3: substitution error rate 0.5%
S4: substitution error rate 1.0%
S5: reads from paralogous/shared transcript regions
S6: high-expression skew: 10 genes dominate 80% of reads
S7: low-complexity/poly-A stress set
```

Acceptance targets for MVP:

```text
S1 exact unique reads: >= 99% correct gene assignment
S2: >= 98% correct gene assignment for uniquely identifiable reads
S3: >= 95% correct gene assignment for uniquely identifiable reads
Ambiguous truth reads: marked ambiguous, not falsely unique, >= 95%
Candidate generator true-origin recall S1: >= 99.5%
Candidate generator true-origin recall S3: >= 98%
```

### 14.4 Public Dataset Benchmarks

Use public 10x PBMC datasets as smoke tests after synthetic correctness is stable.

Suggested datasets:

- 10x Genomics 10k PBMCs from a Healthy Donor, v3 chemistry.
- 10x Genomics 5k Human PBMCs, 3′ v3.1.

Compare against:

```text
STARsolo
kallisto|bustools / kb-python
Cell Ranger output, if available from dataset page
```

Metrics:

```text
wall-clock time
peak RSS
reads/sec
gene count matrix Pearson correlation against baseline
Spearman correlation against baseline
per-cell total UMI correlation
fraction of reads assigned to genes
fraction ambiguous/unmapped
```

Initial public-data target:

```text
Gene-level count correlation with baseline >= 0.95 on common confidently assigned genes.
```

Do not claim superiority from public-data benchmarks until the assignment semantics are well matched.

### 14.5 Cache Benchmarks

Run with Linux `perf`:

```bash
scripts/perf_map.sh cachebatch map \
  --index transcript.cbidx \
  --r1 sim_R1.fastq.gz \
  --r2 sim_R2.fastq.gz \
  --chemistry tenx-3p-v3 \
  --mode candidate-locus \
  --score-mode pulp \
  --emit-metrics metrics.json \
  --out assignments.tsv
```

Compare:

```text
read-at-a-time vs bucketed
random read order vs sketch-sorted order
candidate unsorted vs candidate-locus sorted
scalar vs pulp
```

Key claim to test:

```text
candidate-locus regrouping should reduce LLC misses per scored candidate and improve reads/sec.
```

---

## 15. Correctness Hazards

Handle or explicitly flag:

- high-frequency k-mers causing giant posting lists;
- low-complexity/poly-A sequences;
- Ns and low-quality bases;
- shared exons/transcript regions;
- paralogous genes and pseudogene-like sequences;
- short reads with insufficient informative seeds;
- reads crossing splice junctions if using transcript sequences;
- barcode/UMI parse errors;
- ambiguous gene assignment;
- output determinism under parallel execution.

For now, prefer conservative ambiguity over false unique assignment.

---

## 16. Research Extensions After MVP

Do not implement until the core benchmark exists.

### 16.1 Probability-Aware Early Stopping

Maintain a posterior over candidate genes/transcripts:

```text
P(gene | seed evidence, score evidence, read quality, optional expression prior)
```

Stop seed lookup once posterior mass is concentrated enough, but keep a conservative fallback for ambiguous/rare cases.

### 16.2 IDF-Weighted Seeds

Weight seeds by inverse document frequency:

```text
w(kmer) = log(N / df(kmer))
```

Ignore or downweight seeds with high transcriptome/genome frequency.

### 16.3 Exact or Near-Exact Read Deduplication

For scRNA-seq, many R2 sequences may repeat. Consider mapping unique R2 sequences once and propagating assignments back to CB/UMI records. Preserve enough metadata to avoid incorrect UMI/barcode behavior.

### 16.4 Genomic / Splice-Aware Index

Extend from transcriptome FASTA to:

```text
exons + introns + splice junction windows
```

This can support pre-mRNA/nucleus datasets and better mimic STARsolo/Cell Ranger exonic/intronic logic.

### 16.5 Better SIMD / Architecture-Specific Kernels

If `pulp` is insufficient for a hot kernel:

- add `std::arch` AVX2/AVX-512 implementation;
- keep the unsafe code small;
- preserve scalar and `pulp` fallbacks;
- benchmark before and after.

---

## 17. First Concrete Implementation Task

Codex should begin with this exact sequence:

1. Create the Cargo workspace and crates.
2. Implement `dna.rs` with 2-bit encoding, decoding, reverse complement, and k-mer iteration.
3. Add unit/property tests for `dna.rs`.
4. Implement a tiny FASTA transcript parser using `needletail`.
5. Implement an in-memory transcript k-mer index builder.
6. Add `cachebatch index` CLI that writes a JSON or simple binary prototype index. JSON is acceptable for the first compile, but binary arrays should follow soon.
7. Implement `simulate` for synthetic R1/R2 FASTQ and ground truth.
8. Implement a naive read-at-a-time mapper.
9. Implement sketch records and sketch sorting.
10. Implement candidate-locus regrouping.
11. Implement scalar scoring.
12. Implement `pulp` scoring and scalar-vs-SIMD differential tests.
13. Add Criterion benchmarks.
14. Add `scripts/perf_map.sh`.

Do not skip tests. The most important early success is a stable correctness harness plus benchmark modes that can test the cache hypothesis.

---

## 18. Expected Early Result

The first prototype does not need to beat STARsolo or kallisto. It needs to answer:

```text
Does sketch batching + candidate-locus regrouping reduce cache misses and scored candidates per read compared with a read-at-a-time implementation?
```

A successful early result would be:

```text
- same assignment results as scalar baseline on synthetic data;
- lower LLC misses per scored candidate in candidate-locus mode;
- measurable scoring speedup from pulp on regular candidate buckets;
- clear benchmark report showing where time is spent.
```

Only after that should the project expand toward production-grade single-cell quantification.

---

## 19. README Summary Paragraph

Use this in the initial README:

> `cachebatch` is an experimental Rust RNA-seq/scRNA-seq mapper focused on cache-aware batching. Instead of processing each read independently, it sketches reads into similarity buckets, performs seed lookup, regroups candidate hits by transcript/reference locus, and scores many read-candidate pairs against cache-resident sequence windows. The project uses scalar reference implementations for correctness and `pulp`-based SIMD kernels for regular scoring filters. The first target is transcriptome-level 10x-style gene-expression assignment with synthetic and public PBMC benchmarks.


# Constellation

`constellation` is an experimental Rust RNA-seq/scRNA-seq mapper focused on cache-aware batching. Instead of processing each read independently, it sketches reads into similarity buckets, performs seed lookup, regroups candidate hits by transcript/reference locus, and scores many read-candidate pairs against cache-resident sequence windows. The project uses scalar reference implementations for correctness and `pulp`-based SIMD kernels for regular scoring filters. The first target is transcriptome-level 10x-style gene-expression assignment with synthetic and public PBMC benchmarks.

This repository was initialized from the original `cachebatch` project brief, with the project and CLI renamed to `constellation`.

## Current MVP

- Cargo workspace with `constellation-core` and `constellation-cli`.
- 2-bit DNA encoding, decoding, reverse complement, k-mer iteration, and tests.
- Transcript FASTA parser and in-memory k-mer index builder.
- Compact mmap-backed `.cbidx` serialization, with JSON loading retained for prototype compatibility.
- Streaming compact-index builder for large FASTA targets that avoid holding all postings in memory.
- `constellation index`, `inspect-index`, `simulate`, `map`, and `bench-report` commands.
- Deterministic synthetic paired FASTQ generation.
- Read-at-a-time, sketch-bucket, and candidate-locus mapping modes for cache-locality comparisons.
- Gzip and plain FASTQ input support.
- `needletail` FASTA parsing, `rayon` sorting, and mmap-backed index loading.
- Optional `--gtf` transcript-to-gene mapping during indexing.
- GTF/genome-derived target generation for exon transcripts, gene bodies, introns-only, and exon-plus-gene-body targets.
- Low-quality and low-complexity read classification.
- Rare-anchor seed selection that prefers lower-frequency k-mers before candidate lookup.
- Gene-level k-mer document frequency in compact v3 indexes, stored as a compact side array for `gene-idf` seed planning.
- Seed-batched retrieval mode that groups query seed occurrences by k-mer and loads each selected posting list once per batch.
- Hot candidate-search tables use 16-byte k-mer entries and 8-byte postings so reference slices can fit more cleanly into L1/L2 cache.
- Scalar scoring plus a `pulp` feature-gated SIMD Hamming path with scalar-vs-SIMD differential tests.
- Configurable mismatch, poly-A/poly-T trimming, low-quality tail trimming, and softclip fallback scoring.
- Assignment TSV output, JSON metrics, and truth-aware benchmark reporting for simulated reads.
- UMI-deduplicated `count` output as Matrix Market plus barcode/features TSVs.
- Split assignment metrics for unique gene, same-gene multi-transcript, multi-gene ambiguous, and gene-countable rates.
- Per-read unmapped diagnostics and aggregate score-failure decomposition.

## Current Performance Snapshot

Latest local benchmark: 100k read pairs from 10x PBMC 10k v3, lane L001 subset, mapped with release build against an Ensembl 93 exon-plus-gene-body compact index.

```text
config                         max_mismatches=6, right softclip=16, polyA/polyT/lowQ trim
wall time                      1.32 s
throughput                     75.7k reads/s
unique gene rate               26.34%
same-gene multi-transcript     41.72%
multi-gene ambiguous           10.62%
gene-countable rate            68.07%
unmapped rate                  20.51%
```

Mismatch sweep on the same target showed stable speed and mostly stable gene identities:

```text
max mismatches   gene-countable   unmapped   throughput
4                66.54%           22.44%     73.6k reads/s
6                68.07%           20.51%     75.7k reads/s
8                69.29%           19.05%     76.4k reads/s
10               70.22%           17.90%     77.4k reads/s
```

For `max_mismatches=6`, 99.84% of reads already countable at `max_mismatches=4` kept the same gene call. Aggregate gene counts had Pearson correlation around 0.942 against the public 10x filtered matrix, using the full 10x matrix as a coarse reference.

The current full exon-plus-gene-body index on disk was built before compact v3 gene-DF support, so `--seed-planner gene-idf` falls back to raw posting frequency on that file. Rebuild the large index to get real gene-level IDF seed planning.

## Example

```bash
cargo run -p constellation-cli -- index \
  --transcripts tests/tiny_transcriptome.fa \
  --k 15 \
  --out /tmp/tiny.cbidx

cargo run -p constellation-cli -- inspect-index --index /tmp/tiny.cbidx

cargo run -p constellation-cli -- simulate \
  --transcripts tests/tiny_transcriptome.fa \
  --num-reads 1000 \
  --read-len 40 \
  --out-prefix /tmp/sim

cargo run -p constellation-cli -- map \
  --index /tmp/tiny.cbidx \
  --r1 /tmp/sim_R1.fastq \
  --r2 /tmp/sim_R2.fastq \
  --chemistry tenx-3p-v3 \
  --mode candidate-locus \
  --seed-planner raw-frequency \
  --retrieval-mode per-read \
  --out /tmp/assignments.tsv \
  --emit-metrics /tmp/metrics.json

cargo run -p constellation-cli -- bench-report \
  --assignments /tmp/assignments.tsv \
  --truth /tmp/sim_truth.tsv \
  --index /tmp/tiny.cbidx \
  --metrics /tmp/metrics.json

cargo run -p constellation-cli -- count \
  --assignments /tmp/assignments.tsv \
  --index /tmp/tiny.cbidx \
  --out-prefix /tmp/counts
```

To compile and run the `pulp` SIMD scoring path through the CLI:

```bash
cargo run -p constellation-cli --features simd -- map \
  --index /tmp/tiny.cbidx \
  --r1 /tmp/sim_R1.fastq \
  --r2 /tmp/sim_R2.fastq \
  --mode candidate-locus \
  --score-mode pulp \
  --out /tmp/assignments.pulp.tsv
```

Mapping modes:

```text
read-at-a-time   baseline mode; scores each read's candidate hits independently
sketch-bucket    sketch-sorted read order without candidate-locus regrouping
candidate-locus  sketch-sorted reads plus candidate hit regrouping by transcript/locus
```

Seed planners:

```text
raw-frequency  prefer selected seeds with fewer raw postings
gene-idf       prefer low gene-document-frequency seeds and weight candidate votes by gene-level IDF
```

Retrieval modes:

```text
per-read      load selected posting lists independently for each read
seed-batched  sort selected query seeds by k-mer, load each posting list once, then reduce candidate votes
```

Target kinds:

```bash
cargo run -p constellation-cli -- build-transcriptome-target \
  --genome /path/to/genome.fa \
  --gtf /path/to/annotation.gtf \
  --target-kind exon-plus-gene-body \
  --out /tmp/target.fa
```

```text
exon-transcripts
gene-bodies
introns-only
exon-plus-gene-body
```

Simulation scenarios:

```text
exact
low-complexity
high-expression-skew
```

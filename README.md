# Constellation

`constellation` is an experimental Rust RNA-seq/scRNA-seq mapper focused on cache-aware batching. Instead of processing each read independently, it sketches reads into similarity buckets, performs seed lookup, regroups candidate hits by transcript/reference locus, and scores many read-candidate pairs against cache-resident sequence windows. The project uses scalar reference implementations for correctness and `pulp`-based SIMD kernels for regular scoring filters. The first target is transcriptome-level 10x-style gene-expression assignment with synthetic and public PBMC benchmarks.

This repository was initialized from the original `cachebatch` project brief, with the project and CLI renamed to `constellation`.

## Current MVP

- Cargo workspace with `constellation-core` and `constellation-cli`.
- 2-bit DNA encoding, decoding, reverse complement, k-mer iteration, and tests.
- Transcript FASTA parser and in-memory k-mer index builder.
- Compact mmap-backed `.cbidx` serialization, with JSON loading retained for prototype compatibility.
- `constellation index`, `inspect-index`, `simulate`, `map`, and `bench-report` commands.
- Deterministic synthetic paired FASTQ generation.
- Read-at-a-time, sketch-bucket, and candidate-locus mapping modes for cache-locality comparisons.
- Gzip and plain FASTQ input support.
- `needletail` FASTA parsing, `rayon` sorting, and mmap-backed index loading.
- Optional `--gtf` transcript-to-gene mapping during indexing.
- Low-quality and low-complexity read classification.
- Rare-anchor seed selection that prefers lower-frequency k-mers before candidate lookup.
- Hot candidate-search tables use 16-byte k-mer entries and 8-byte postings so reference slices can fit more cleanly into L1/L2 cache.
- Scalar scoring plus a `pulp` feature-gated SIMD Hamming path with scalar-vs-SIMD differential tests.
- Assignment TSV output, JSON metrics, and truth-aware benchmark reporting for simulated reads.
- UMI-deduplicated `count` output as Matrix Market plus barcode/features TSVs.
- GTF/genome-derived exon transcript FASTA generation with `build-transcriptome-target`.

## Current Performance Snapshot

Latest local benchmark: 100k read pairs from 10x PBMC 10k v3, lane L001 subset, mapped with release build against the Ensembl 93 GTF-derived exon transcript compact index.

```text
wall time             0.571 s
throughput            175k reads/s
unique gene rate      13.65%
ambiguous gene rate   35.59%
unmapped rate         49.96%
low complexity rate    0.73%
low quality rate       0.07%
```

Function-family timing from that run:

```text
candidate generation   240 ms   generate_candidate_hits_with_quality_stats_parallel
scoring                113 ms   CandidateScorer::score_bucket plus candidate sorting/grouping
FASTQ loading           95 ms   read_fastq_r1/read_fastq_r2
bucket building         32 ms   make_candidate_locus_buckets
preprocess/sketching    21 ms   prepare_reads_parallel plus sketch sorting
assignment              14 ms   assignments_by_read
output writing          22 ms   assignment TSV and metrics output
index load            <0.1 ms   mmap-backed compact index
```

The previous Ensembl cDNA compact index on the same 100k subset had an unmapped rate of about 53.0%, so the GTF-derived target is a modest improvement but does not fully explain the remaining high unmapped fraction.

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

Simulation scenarios:

```text
exact
low-complexity
high-expression-skew
```

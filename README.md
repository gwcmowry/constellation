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
- Sparse-probe candidate search with fallback for faster average-case lookup.
- Streamed seed-batched mapping chunks so candidate generation, scoring, and assignment run on bounded working sets.
- Gene-level WAND-lite candidate pruning with configurable score ratio and max loci per retained gene.
- Reverse-complement seed lookup generated on the fly from encoded query k-mers, without building a larger reverse-complement index.
- Hot candidate-search tables use 16-byte k-mer entries and 8-byte postings so reference slices can fit more cleanly into L1/L2 cache.
- Scalar scoring plus a `pulp` feature-gated SIMD Hamming path with scalar-vs-SIMD differential tests.
- Configurable mismatch, poly-A/poly-T trimming, TSO trimming, low-quality tail trimming, and softclip fallback scoring.
- Assignment TSV output, JSON metrics, and truth-aware benchmark reporting for simulated reads.
- UMI-deduplicated `count` output as Matrix Market plus barcode/features TSVs.
- Split assignment metrics for unique gene, same-gene multi-transcript, multi-gene ambiguous, antisense-gene diagnostics, and gene-countable rates.
- Per-read unmapped diagnostics and aggregate score-failure decomposition.

## Current Performance Snapshot

Latest local benchmark: 100k read pairs from 10x PBMC 10k v3, lane L001 subset, mapped with release build against the Ensembl 93 exon-plus-gene-body compact v3 index:

```text
benchdata/reference/ensembl93/human_ensembl93_exon_plus_gene_body.k21.genedf.cbidx
```

The current fast path is seed-batched, sparse-probe candidate search, streamed chunks, gene-level WAND-lite pruning, reverse-complement search, TSO/polyA/polyT/low-quality trimming, and softclip scoring. Use the internal metrics JSON wall-clock fields for timing; process-level `/usr/bin/time` includes mmap/load/teardown effects that are not comparable with the mapper hot path.

```text
config                         max_mismatches=6, right softclip=16, RC + TSO/polyA/polyT/lowQ trim
internal wall time             0.669 s
mapping hot-path time          0.512 s
throughput                     149.5k reads/s
unique gene rate               42.77%
same-gene multi-transcript     33.53%
multi-gene ambiguous           10.07%
gene-countable rate            76.30%
unmapped rate                  12.68%
low-complexity rate            0.88%
```

Reverse-complement search is a large recall gain but costs throughput. On the same 100k subset:

```text
mode             gene-countable   unmapped   internal throughput
forward only     68.03%           20.41%     172.4k reads/s
RC enabled       73.98%           15.52%     165.0k reads/s
RC + TSO trim    76.30%           12.68%     149.5k reads/s
```

Cell Ranger comparison is currently only against `filtered_feature_bc_matrix`, which is a cell/gene UMI matrix rather than read-level mapping truth. On the 100k PBMC subset, about 91.5% of Constellation UMIs land in Cell Ranger filtered barcodes, all Cell Ranger genes are present in the Ensembl 93 reference, and about 5.6% of Constellation UMIs are assigned to genes outside the Cell Ranger feature set. Log gene-total correlation against the filtered matrix is about 0.51 Pearson and 0.56 Spearman on this tiny read subset, so use it as a rough sanity check rather than a mapping-rate target.

## Handoff Status

Implemented recall and speed work from the Codex mapping brief:

```text
done   gene-countable same-gene transcript ambiguity
done   split assignment metrics and score-failure diagnostics
done   exon-plus-gene-body target generation for intron/pre-mRNA support
done   compact v3 gene-DF index support and gene-idf seed planning
done   sparse-probe candidate search with full fallback
done   streamed seed-batched candidate generation/scoring/assignment
done   WAND-lite and gene-level WAND-lite pruning
done   on-the-fly reverse-complement lookup
done   TSO trimming and antisense diagnostic assignment class
```

Open architecture questions for the next handoff:

```text
1. Cell Ranger comparability:
   add a feature whitelist/evaluation mode based on Cell Ranger features.tsv so assignments to genes outside the 10x feature set can be separated from real recall loss.

2. Target class priority:
   track exon, intron, and gene-body hits separately and prefer exonic confident hits over intronic/gene-body hits when both explain the read.

3. Barcode and UMI correction:
   implement 10x whitelist barcode correction and molecule-level UMI correction before comparing counts to Cell Ranger output.

4. Rescue alignment:
   add a bounded fallback for reads with weak exact-k seed support, such as shorter/spaced seeds plus local alignment, instead of only loosening the full fast path.

5. Reference architecture:
   evaluate a Cell Ranger-like feature reference or splici-style reference so high genome mapping rates do not come from reads that are not countable transcriptome evidence.
```

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

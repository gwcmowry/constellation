# PathCert

PathCert is a research prototype for mapping RNA-seq reads to path-equivalence
certificates over per-locus transcript components.

This crate is intentionally separate from the parent Constellation workspace.
The first implementation proves the MVP data model:

- transcript components grouped by gene or locus
- compact `u64` pathsets for components with up to 64 transcript paths
- sorted seed lookup from k-mer to component-local pathset
- component-bucketed read resolution by pathset intersection
- stable certificate IDs and EC count aggregation

## Build

```bash
cargo test
cargo run -- index --transcripts ../tests/tiny_transcriptome.fa --k 7 --out /tmp/tiny.pcidx
cargo run -- map --index /tmp/tiny.pcidx --reads ../tests/tiny_reads_R2.fastq --out /tmp/tiny.ec.tsv --summary /tmp/tiny.summary.tsv
```

## Current Scope

The mapper is single-end and exact-seed based. GTF/FASTA component construction
is present for small annotations, but the first-class tested path is the toy
component and transcript FASTA path. Paired-end constraints, slow validation,
SIMD kernels, and memory-mapped binary index layouts are explicit follow-up
work rather than hidden assumptions in this MVP.

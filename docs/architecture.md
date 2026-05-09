# Architecture

The current architecture is the gene-EC mapper. Older positional compact-index modes remain only as historical code paths and should not be used for current performance comparisons.

## Mapper Flow

```text
paired FASTQ streams
  -> bounded R1/R2 decode batches
  -> 10x barcode/UMI parse and R2 quality checks
  -> sparse k-mer probes
  -> prefix24 mmap EC-index lookup
  -> per-gene score accumulation
  -> gene assignment
  -> Constellation RAD-like binary or TSV output
  -> optional streaming zstd compression
```

The important runtime property is that FASTQ decoding is pipelined ahead of mapping, so the mapper mostly waits on EC lookup and assignment work rather than gzip input.

## Index

The active index format is the mmap EC index built by `constellation index-ec`. The current best human reference used in benchmarks is:

```text
/tmp/human_ensembl93_splici_r91.k31.t2g.prefix24.mmap.ecidx
```

The prefix24 table reduces lookup work but the index still has a high RSS footprint because postings are explicit and are accessed through large mmap-backed arrays. Current perf profiles show `LoadedEcIndex::lookup` dominates dTLB misses, so future index work should prioritize locality and translation pressure.

## Output

The preferred mapper output is `--output-format gene-ec-rad`. It is a Constellation-specific compact per-read format, not alevin-fry-compatible RAD. Streaming zstd is available with `--output-compression zstd`; it is cheap enough for normal benchmark runs but does not remove the need for molecule-level aggregation.

## Current Constraints

- Memory is the main gap versus simpleaf/alevin-fry: about `25.3G` RSS for Constellation versus about `3.93G` RSS for simpleaf on the 100M PBMC benchmark.
- Per-read output is still larger than ideal even after zstd.
- The output semantics are not yet identical to alevin-fry because Constellation emits per-read gene/EC assignments rather than a full molecule-resolution RAD/counting pipeline.

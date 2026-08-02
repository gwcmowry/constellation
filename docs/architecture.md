# Architecture

The current architecture is the gene-EC mapper. Older positional compact-index modes remain only as historical code paths and should not be used for current performance comparisons.

## Mapper Flow

```text
paired FASTQ streams
  -> bounded R1/R2 decode batches
  -> 10x barcode/UMI parse and R2 quality checks
  -> sparse k-mer probes
  -> prefix24 mmap EC-index lookup
  -> per-worker dense epoch-tagged gene score accumulation
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

The prefix24 table reduces lookup work but the index still has a high RSS footprint because postings are explicit and are accessed through large mmap-backed arrays. The controlled 10M path uses about 16.8GiB max RSS for a 17GB index. Historical perf profiles show `LoadedEcIndex::lookup` dominates dTLB misses, so future index work must reduce the random-page working set rather than merely reorder the same probes. The failed larger-layout and batching families are recorded in `docs/failed_performance_experiments.md`.

The direct mapper reuses one dense score/epoch array per Rayon worker. Reset cost is proportional to genes touched by the read, and score updates avoid a per-read hash table. This was exact but only a modest speedup, confirming that accumulator writes are secondary to index access.

## Output

The preferred mapper output is `--output-format gene-ec-rad`. It is a Constellation-specific compact per-read format, not alevin-fry-compatible RAD. Streaming zstd is available with `--output-compression zstd`; it is cheap enough for normal benchmark runs but does not remove the need for molecule-level aggregation.

## Current Constraints

- Memory is the main gap versus simpleaf/alevin-fry: the current controlled 10M mapper uses about 16.8GiB, while the older 100M comparison reached about 25.3G versus about 3.93G for simpleaf.
- Per-read output is still larger than ideal even after zstd.
- The output semantics are not yet identical to alevin-fry because Constellation emits per-read gene/EC assignments rather than a full molecule-resolution RAD/counting pipeline.
- The current `count` path does not directly consume the preferred gene-EC RAD-like output, so FASTQ-to-matrix is not yet one coherent frontier pipeline.

# Architecture

The current architecture is the gene-EC mapper. Older positional compact-index modes remain only as historical code paths and should not be used for current performance comparisons.

## Mapper Flow

```text
paired FASTQ streams
  -> bounded R1/R2 decode batches
  -> 10x barcode/UMI parse and R2 quality checks
  -> sequence-intrinsic closed-syncmer probes
  -> prefix27 mmap EC-index lookup
  -> per-worker dense epoch-tagged gene score accumulation
  -> allocation-free direct gene resolution
  -> Constellation RAD-like binary or TSV output
  -> optional streaming zstd compression
```

The important runtime property is that FASTQ decoding is pipelined ahead of mapping, so the mapper mostly waits on EC lookup and assignment work rather than gzip input.

## Index

The active index format is the mmap EC index built by `constellation index-ec --syncmer-s 20`. The current best human reference used in benchmarks is:

```text
/tmp/human_ensembl93_splici_r91.k31.syncmer-s20.prefix27.mmap.ecidx
```

The hashed closed-syncmer rule is sequence- and strand-intrinsic: reference and
query select the same exact k-mer independently of absolute transcript/read
position. A rolling canonical s-mer code and fixed monotonic queue make selection
linear and allocation-free. Sampling metadata is stored in the v4 reserved header
fields; incompatible lexicographic prototype indexes are rejected rather than
silently queried with different semantics.

The 27-bit directory averages roughly two sampled keys per bucket on the current
human index. The controlled 10M path uses a 3.538GB index and at most 3.964GiB
mapping RSS, versus 17.892GB and 16.823GiB for dense prefix24. Three-run median
throughput is 1.068M reads/s versus 1.049M for dense. Prefix24 and prefix26 sampled
variants are superseded; the rejected densities and rescue cascade are recorded
in `docs/failed_performance_experiments.md`.

The direct mapper reuses one dense score/epoch array per Rayon worker. Reset cost
is proportional to genes touched by the read, and score updates avoid a per-read
hash table. Direct resolution scans touched scores in place; only output modes
that actually need tied top-gene IDs allocate that vector.

## Output

The preferred mapper output is `--output-format gene-ec-rad`. It is a Constellation-specific compact per-read format, not alevin-fry-compatible RAD. Streaming zstd is available with `--output-compression zstd`; it is cheap enough for normal benchmark runs but does not remove the need for molecule-level aggregation.

## Current Constraints

- Mapper memory is no longer the defining gap on the controlled workload, but a matched end-to-end comparison against simpleaf/alevin-fry still needs to be rerun with the sampled index.
- Per-read output is still larger than ideal even after zstd.
- The output semantics are not yet identical to alevin-fry because Constellation emits per-read gene/EC assignments rather than a full molecule-resolution RAD/counting pipeline.
- The current `count` path does not directly consume the preferred gene-EC RAD-like output, so FASTQ-to-matrix is not yet one coherent frontier pipeline.

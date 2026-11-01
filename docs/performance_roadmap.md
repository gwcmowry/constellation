# Performance and frontier roadmap

## Scaling and usability update 2026-11-01

The [100M/500M comparison](../experiments/2026-11-01-large-scale/README.md)
now measures complete native mapping and counting, including USA output, and
links the independently measured index builds. At 500M, the conservative
workflow takes 322.25s versus 742.65s for piscem/alevin-fry and 1533.78s for
kallisto/bustools (one run per configuration on twelve physical cores).

The product workflow is now `index-ec` once per reference, then `quant` once per
sample, with `--threads`, a bundled installation example, automatic metrics and
intermediate cleanup. `quant` orchestrates the existing stages; it does not fuse
them. The default conservative and optional exact-UMI USA policies remain explicit.

Next performance priorities are serial evidence ingestion, retained counting
memory, and avoiding the CSTRAD2 round trip. Distribution could be simplified
further with tested prebuilt executables and versioned reference/index bundles;
the supported installation path currently builds from source. Dated sections
below preserve the earlier reasoning and should be read in that context.

## Build and counting update 2026-10-31

The partitioned builder and direct v9 writer remove the serial external merge
and intermediate v4 conversion from the recommended build command. Slot-range
partitioning also replaces the serial random-write fill. The state counter now
keeps packed, cell-partitioned evidence and resolves cells independently, avoiding
a global per-molecule hash table and repeated sparse tree insertion.

See the [build-to-matrix experiment](../experiments/2026-10-31-build-count/README.md)
for fresh matched timings, memory, full index verification and native USA parity.
The optional USA resolver recovers additional assignments through explicit
spliced preference, with a measured source-origin precision cost. It does not
solve missing discriminating seeds. The next accuracy priority remains targeted
reference/path verification of ambiguous and weak reads.

## Comparison and quality update 2026-10-30

The matched forward-strand comparison now covers piscem 0.24.0, kallisto 0.52.0,
two PBMC lanes, and fresh 91bp/150bp genome/GTF simulations. The optional
`--gene-ec-consensus` rule retains all genes supported by every distinct matched
seed, requires two seeds, and avoids annotation-based tie resolution. Its median
10M-pair time is 4.405s versus piscem's 13.958s. On the 91bp holdout, unique
precision improves from 99.064% to 99.945%, but correct unique recall falls from
90.155% to 85.181% (piscem: 87.344%). At 150bp, consensus recall is 88.977%,
slightly above both competitors, with 99.984% precision.

Prioritize targeted reference/path verification to recover weak short reads
without restoring forced calls, then bring the bounded prefetch strategy into
consensus lookup. A decoy-aware reference and native molecule-resolution
comparison remain necessary before claiming broader biological equivalence.
The [full report](../experiments/2026-10-30-consensus/README.md) corrects the
orientation mismatch in the initial October 27 competitor comparison and records
remaining scope limits. This is not yet a universal SOTA result.

## Historical update 2026-10-27

The packed pipeline and exact v9 prefetch path are retained. On the matched
10M PBMC/12-core workload, no-output throughput is 3.641M reads/s and RAD+zstd
throughput is 3.357M reads/s, improvements of 2.90x and 3.00x over the original
PtrHash binary. Complete output bytes and all assignment metrics are unchanged;
compressed-output max RSS is 4.265GiB with the same 4.428GB index. Reusable slabs,
one persistent state per worker, fused preparation/output, bounded prefetching
and zlib-rs decoding replaced the dominant allocation and latency costs.

The sequence-aware unitig index and shared hot-key table remain proposed work,
not implemented optimizations. Exact EC-run aggregation was tested and removed
because it did not improve end-to-end time. The remaining large product gap is
still molecule/count-matrix integration; this mapper speedup does not resolve
output equivalence with alevin-fry counting. See
`experiments/2026-10-27-packed-pipeline/README.md` for measurements and gates.

## Historical stocktake (2026-08-02)

Constellation is a fast gene-EC mapping prototype, but it is not yet a frontier
single-cell tool worth releasing as a complete product.

On the controlled 10M PBMC workload at 12 threads, the retained hashed closed-
syncmer `k=31,s=20` prefix27 path has a three-run median of 1.068M reads/s without
output, 78.50747% gene-countable, 9.49510% ambiguous, and 11.26884% unmapped.
Its index is 3.538GB and mapping max RSS is 3.964GiB. The matched dense prefix24
control has a 1.049M reads/s median, 78.53896% countable, a 17.892GB index, and
16.823GiB max RSS. Sampling therefore delivered 5.06x smaller index, 4.24x lower
RSS, and 1.82% higher median throughput, with -0.03149 percentage-point countable
and +0.22683-point ambiguity deltas.

The sampled rule mixes canonical rolling s-mers and selects closed syncmers with
a fixed-size monotonic queue. The rule is self-described in the mmap header and
does not compose with absolute read-position stride. A 27-bit prefix directory
was the best measured balance: prefix24 was slower, prefix26 nearly tied dense,
and prefix27 matched/exceeded dense speed while remaining compact.

The exact allocation-free resolver also scans epoch-tagged scratch scores in
place instead of allocating and sorting every candidate vector. It reduced the
mapping stage by roughly 4% in the paired single runs, while differential tests
and dense aggregate rates remained identical.

There are also product blockers beyond mapping speed:

- the preferred gene-EC RAD-like output is Constellation-specific and the current
  `count` path does not consume it directly;
- molecule-level barcode correction, UMI resolution, and count-matrix semantics
  are not yet output-equivalent to alevin-fry/Cell Ranger;
- the sampled human index build is about 2:06 and 4.90GiB peak RSS, but build
  reproducibility and temporary-space reporting still need a packaged script;
- the prior count comparison was only about 0.94 correlation on common genes;
- unprivileged hardware-counter profiling is blocked by `perf_event_paranoid=4`;
  `sudo perf stat` works and its raw dense/p24 comparison is archived, while
  cargo-flamegraph/perf-record sampling remains unavailable to the agent.

## Status of high-leverage bets

### 1. Sequence-intrinsic sampled index

Completed. The result exceeded the size and memory targets but not the original
1.5x throughput target. It is retained because p27 is a strict practical Pareto
improvement over dense on size, RSS, build cost, and measured median speed, with
only a 0.03149-point countable change. Lexicographic ordering, sparse hashed
`s=11`, dense hashed `s=21`, two-index rescue, and p24/p26 directories are closed;
see the failure ledger and raw archive before revisiting them.

### 2. Shared immutable LFU hot-key table

The 1M-read seed trace contained 5,996,875 occurrences and 3,384,891 unique keys
(1.772× reuse). The top 1,000/10,000/100,000/1,000,000 keys represented
5.20%/17.38%/35.12%/60.23% of occurrences, while 83.27% of unique keys were
singletons.

That explains why the rejected direct-mapped cache lost: it admitted cold keys,
duplicated state per worker, and polluted cache. A genuinely different test is
an immutable, shared table containing only frequency-selected hot keys learned
from an initial sample or a reusable chemistry/reference profile. It needs one
copy, no writes in the map phase, and a cheap membership filter. Include the
sampling/prepass cost in end-to-end time; remove it unless total runtime wins.

### 3. Molecule-first output and counting

The standalone state counter and `quant` now provide FASTQ-to-matrix processing.
The remaining proposal here is fused, bounded-memory aggregation:

Stop treating a variable-length record per read as the final product. Carry the
already packed 16-base barcode and 12-base UMI through mapping, aggregate
`(barcode, UMI, gene/EC)` records in bounded partitions, apply explicit barcode
and UMI correction, and emit a count-ready molecule stream or matrix. Make
`count` consume the preferred mapper output directly.

This can reduce allocations, output size, and downstream I/O while closing the
largest release-readiness gap. Benchmark complete FASTQ-to-matrix wall time and
matrix concordance, not mapper-only time.

### 4. Allocation-free direct resolution

Completed and retained. The direct mapper scans touched dense scores without an
owned all-candidate vector; RAD top-gene IDs are collected only when output is
actually enabled. The old sorted and new direct results are differentially tested.

## Release bar

Do not call the tool frontier or prepare a public release until one benchmark
script demonstrates, on public data and at the 12-thread cap:

1. FASTQ-to-count-matrix wall time competitive with the best current pipelines.
2. Count-matrix and molecule-level concordance under explicitly matched reference
   and filtering semantics.
3. Preserve the sampled path's substantial memory/index-size improvement and
   validate it on the full end-to-end workload.
4. Warm and cold results, index build cost, output size, and three-run variance.
5. Byte-level or semantic regression tests for every accepted optimization.

See `docs/failed_performance_experiments.md` before starting any experiment.

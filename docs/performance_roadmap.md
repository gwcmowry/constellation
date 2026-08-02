# Performance and frontier roadmap

## Current stocktake (2026-08-02)

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

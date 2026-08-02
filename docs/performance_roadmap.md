# Performance and frontier roadmap

## Current stocktake (2026-08-02)

Constellation is a fast gene-EC mapping prototype, but it is not yet a frontier
single-cell tool worth releasing as a complete product.

On the controlled 10M PBMC workload at 12 threads, the accepted direct mapper is
about 1.044–1.052M reads/s without output, 78.54% gene-countable, 9.27% ambiguous,
and 11.46% unmapped. Max RSS is about 17.64 million KiB (16.8 GiB). The 17GB v4
index remains the defining constraint. A cold-cache run took 28.5 seconds and
incurred about 3.0 million major faults, versus roughly 9.5 seconds warm.

The accepted 2026-08-02 change reuses dense epoch-tagged gene score arrays per
Rayon worker and avoids full sketch hashing plus a duplicate quality scan on the
direct path. It is byte-identical to the previous gene-EC RAD output. The warm
no-output improvement is modest (roughly 1–2%); one zstd-output pair improved
from 0.948M to 0.987M reads/s, but that larger delta needs more repetitions.

Regularizing only score writes did not move the ceiling much. The historical
profile still gives the best explanation: index lookup accounted for 34.2% of
cycles and 55.6% of dTLB-load misses. A per-worker lookup cache and huge-page
hint both regressed, so the next meaningful improvement must reduce the number
of random pages required, not add another layer in front of the same 17GB table.

There are also product blockers beyond mapping speed:

- the preferred gene-EC RAD-like output is Constellation-specific and the current
  `count` path does not consume it directly;
- molecule-level barcode correction, UMI resolution, and count-matrix semantics
  are not yet output-equivalent to alevin-fry/Cell Ranger;
- the human index build historically took about 13 minutes, 26.7GB RSS, and
  24.5GB of temporary runs;
- the prior count comparison was only about 0.94 correlation on common genes;
- reproducible hardware-counter profiling is blocked on this host by
  `perf_event_paranoid=4` (cargo-flamegraph and `perf stat` cannot collect data).

## New high-leverage bets

### 1. Sequence-intrinsic sampled index

Prototype a closed-syncmer or strobemer-selected EC index instead of storing
every reference k-mer. This is materially different from the failed read-position
stride experiments: the same sequence-intrinsic rule selects keys in the
reference and reads, so the index itself becomes smaller and its working set can
approach memory/cache capacity.

Start with two or three selector densities and build them as separate, simple
v4-style indexes. Do not invent another payload encoding in the first experiment.
Measure index size, build time, warm/cold RSS, page faults, throughput, and exact
assignment deltas. A first pass can use the sampled index; unresolved reads may
use a separately measured rescue strategy only if it does not repeat the failed
absolute-stride rescue policy.

Acceptance target: at least 2× smaller index and 1.5× warm throughput while
retaining countable/ambiguous/unmapped rates within a predeclared tolerance. A
smaller result is unlikely to change the product position enough.

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

The dense score accumulator still creates and sorts an owned `Vec<EcGeneScore>`
for each supported read. Add a specialized direct-mapper result that scans the
touched dense scores to compute the best total/class tie resolution without
sorting all candidates. When RAD top-gene IDs are required, collect and sort only
the tied top genes. Reuse output buffers where ownership permits.

This is exact and low risk, but it is a secondary lever; reject it if it cannot
produce a stable mapper-only gain above noise.

## Release bar

Do not call the tool frontier or prepare a public release until one benchmark
script demonstrates, on public data and at the 12-thread cap:

1. FASTQ-to-count-matrix wall time competitive with the best current pipelines.
2. Count-matrix and molecule-level concordance under explicitly matched reference
   and filtering semantics.
3. A substantial memory/index-size improvement over the current 17GB/16.8GiB path.
4. Warm and cold results, index build cost, output size, and three-run variance.
5. Byte-level or semantic regression tests for every accepted optimization.

See `docs/failed_performance_experiments.md` before starting any experiment.

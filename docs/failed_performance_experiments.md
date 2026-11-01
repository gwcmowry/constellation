# Failed and superseded performance paths

The 2026-10-30 quality experiment first required three distinct matched seeds
spanning at least k/2 start positions, with gene-set intersection. On development
positives this reached 99.977% precision but only 80.068% correct unique recall.
It was not retained for the 91bp workload; the two-seed rule recovered recall
while maintaining approximately 99.95% precision and was tested on fresh holdouts.
Raw results: `experiments/2026-10-30-consensus/raw/strict-development-*`.
Revisit the stricter rule only for a changed read-length/error model or an explicit
preference for further abstention.

The 2026-10-27 packed pipeline retained bounded within-read v9 prefetching,
which differs from the global v4 grouping below. A further exact aggregation
of repeated EC payloads did not improve runtime (2.888 s versus 2.852 s for
its immediate no-output control, single measured runs after warmup), and was
removed. Raw measurements and the rejected patch are in
`experiments/2026-10-27-packed-pipeline/`. Revisit aggregation only with a changed
cost model, such as avoiding index probes through verified sequence traversal.

This is the authoritative negative-results ledger. Treat every row as **do not
retry** unless the underlying premise changes and the new experiment states that
change explicitly. Raw results from the 2026-08-02 audits are in
`experiments/archive/2026-08-02-memory-locality-audit/raw/` and
`experiments/archive/2026-08-02-syncmer-index/raw/`; older large outputs remain
in the git-ignored `benchdata/` corpus.

All CPU rates below used at most 12 cores unless marked otherwise. The comparison
workload is 10 million PBMC 10k v3 read pairs, a warm v4 prefix24 EC index, batch
size 65,536, and no assignment output. Assignment rates are fractions of reads.

## Rejected in the 2026-08-02 audit

| Path | Result | Why it is closed |
|---|---:|---|
| Read-at-a-time v4 baseline | 1.035M reads/s; 78.54% countable; 9.27% ambiguous; 11.46% unmapped; 17.58GB max RSS | Control, not a rejected path. |
| Seed-batched traversal | 0.609M reads/s | 41% slower with identical assignments. Probe materialization, grouping, and per-read finalization cost more than any locality gained. |
| Prefix-tiled traversal | 0.695M reads/s | 33% slower with identical assignments. Partitioning and scattered per-read accumulation erased lookup locality gains. |
| v5 packed suffix/payload layout | 1.057–1.077M reads/s; about 20.7GB RSS; 21.0GB index | Only 2–4% faster in these runs while adding roughly 3.1GB RSS and a larger index. It is not a frontier tradeoff. |
| v6 blocked suffix lookup | 0.755–0.997M reads/s; about 21.2GB RSS | Slower and larger than v4. Blocking did not overcome metadata and search overhead. |
| v7 target-range deduplication | 0.756–0.989M reads/s; about 20.9GB RSS | Slower and still substantially larger in memory than v4. |
| v8 minimal-perfect-hash family | Build stopped; no valid result | The build was intentionally stopped rather than spend more time in the already unsuccessful larger-layout family. An empty timing file is retained. Do not cite it as a benchmark. |
| Probe stride 16 | 1.200M reads/s; 76.31% countable; 12.95% unmapped | Faster only by discarding useful evidence: countable rate fell 2.23 percentage points. |
| Probe stride 24 | 1.276M reads/s; 74.94% countable; 14.15% unmapped | Faster only by a larger quality regression. |
| Stride 16 plus rescue at 12 | 1.098M reads/s; 77.34% countable; 10.92% ambiguous | Not Pareto-optimal: it remained less countable and more ambiguous than stride 12. |
| Stride 24 plus rescue at 12 | 1.173M reads/s; 76.38% countable; 11.43% ambiguous | Restored baseline unmapped rate by converting too many reads to ambiguity and remained less countable. |
| Per-worker direct-mapped EC lookup cache, including misses | 1.003M reads/s | About 4% slower than the accepted dense-score version. Cache metadata and LLC pollution outweighed avoided index probes. Code was removed. |
| `MADV_HUGEPAGE` on the file-backed index mmap | 1.015M reads/s; unchanged RSS | About 3% slower than the accepted version on a host configured for `madvise` THP. Code was removed. |

The first unbounded/concurrent traversal files (`constellation_audit_v4_*` without
`12t`/`serial`) are retained but invalid for comparison: some overlapped and the
read-at-a-time run predates the 12-thread cap. Use the `12t_*_serial` files.

## Sampled-index paths explored on 2026-08-02

The retained result is hashed canonical closed-syncmer `k=31,s=20` with a
prefix27 v4 directory: 1.068M reads/s median, 3.538GB index, 3.964GiB max RSS,
78.50747% countable, 9.49510% ambiguous, and 11.26884% unmapped. Do not replace it
with one of the following variants without a new premise.

| Path | Result | Why it is closed |
|---|---:|---|
| Lexicographic closed syncmer `s=11` | 1.157M reads/s; 2.127GB index; 76.91349% countable | Faster/smaller but loses 1.625 points countable. Lexicographic s-mer order is composition-biased and its on-disk marker is now explicitly rejected. |
| Lexicographic `s=15` | 1.020M reads/s; 2.578GB; 77.80861% countable | Slower than dense while losing 0.730 points countable. |
| Lexicographic `s=21` | 0.844M reads/s; 3.858GB; 79.20886% countable | Better aggregate quality but 19% slower than the dense control. |
| Hashed canonical `s=11` | 1.208M reads/s; 1.789GB; 76.02365% countable | Fastest sampled endpoint, but loses 2.515 points countable. |
| Hashed canonical `s=21` | 0.912M reads/s; 3.340GB; 78.74331% countable | Clears baseline countable rate but is about 13% slower. |
| `s=11` primary plus `s=21` rescue for unmapped reads | 0.996M reads/s; 77.91391% countable; 10.68219% ambiguous; 5.744GiB RSS | Rescue converted many reads to ambiguity, remained below baseline countable, and required two indexes. Code was removed. |
| Hashed `s=20`, prefix24 | 0.995M reads/s median after direct resolution; 3.069GB | Same assignments as retained p27 but about 6.8% slower; short suffix searches still did avoidable comparisons. |
| Hashed `s=20`, prefix26 | 1.042M reads/s median; 3.270GB; 3.460GiB RSS | Nearly matches dense but is dominated by prefix27's higher throughput at an acceptable 268MB additional index cost. |

An allocation-free direct resolver was retained because it is exact and reduced
mapping-stage time by about 4%; it is not a failed sampling variant. The empty
`syncmer_hashed_s20.perf.data` file is also not a profile: `perf record` was
blocked by host policy. The two `*_perf_stat.txt` files were collected manually
with sudo and are valid counter summaries for dense versus pre-resolver p24.

## Older paths that are also closed

The older 100k-read corpus under `benchdata/pbmc_10k_v3/` contains repeated
variants of these ideas. The archived plans describe their motivation, but they
are not an active backlog:

- sketch-sorted read buckets and candidate-locus regrouping;
- hot/cold or two-tier positional indexes and manual locality sorting;
- seed-batched posting retrieval, posting prefetch, and WAND/early-stop variants;
- small-vector/hash accumulator substitutions and parallel preparation variants;
- exon-first/global rescue policies that trade unmapped reads for ambiguity;
- SIMD/pulp candidate scoring and soft-clip fast paths in the legacy positional
  mapper;
- reverse-complement, TSO/poly-A trimming, and permissive scoring as a speed path.

Those experiments either lost outright, improved only a legacy path that is no
longer the recommended architecture, or failed the combined speed/quality/memory
bar. The current gene-EC mapper performs no candidate-window scoring, so reviving
legacy SIMD scoring cannot address its measured hot path.

Compression is also closed as a primary speed target: zstd level 3 was under 1%
of historical cycle samples and current output writing is a small fraction of
wall time. More threads are not the next experiment either; all project work is
capped at 12 threads, and 6/8/12-core scaling is already healthy.

## Revisit gate

A closed path may be revisited only when all of the following are written down
before implementation:

1. The premise that is new (different index density, workload, hardware feature,
   or algorithm—not a renamed version of the same layout).
2. A correctness oracle, preferably byte-identical output or a predeclared
   quality tolerance.
3. A 12-thread baseline and an acceptance threshold covering wall time, max RSS,
   index size/build cost, countable rate, ambiguity, and unmapped rate.
4. A plan to remove the code if it fails.

# Hashed closed-syncmer index experiment (2026-08-02)

## Outcome

Retain hashed canonical closed syncmers at `k=31,s=20` with the v4 prefix27
directory and allocation-free direct gene resolution.

All CPU work was serial across experiments, pinned with `taskset -c 0-11`, and
used at most 12 Rayon/build threads. The workload was 10 million PBMC 10k v3
read pairs, batch size 65,536, warm mmap indexes, and no assignment output.

| Path | Index | Median reads/s | Countable | Ambiguous | Unmapped | Max map RSS |
|---|---:|---:|---:|---:|---:|---:|
| Dense v4 prefix24 | 17.892GB | 1.048994M | 78.53896% | 9.26827% | 11.46418% | 16.823GiB |
| Hashed `s=20`, prefix24 | 3.069GB | 0.995293M | 78.50747% | 9.49510% | 11.26884% | 3.082GiB |
| Hashed `s=20`, prefix26 | 3.270GB | 1.041952M | 78.50747% | 9.49510% | 11.26884% | 3.460GiB |
| **Hashed `s=20`, prefix27** | **3.538GB** | **1.068055M** | **78.50747%** | **9.49510%** | **11.26884%** | **3.964GiB** |

Relative to dense, retained p27 is 5.06x smaller, uses 4.24x less map RSS, and
is 1.82% faster by median throughput. Countable changes by -0.03149 percentage
points; ambiguity changes by +0.22683 points and unmapped by -0.19534 points.
Its human-reference build took 2:06.07 and 4.90GiB peak RSS.

## Current runtime breakdown with real output

One retained-p27 run wrote 10M uncompressed gene-EC RAD records (369,461,060
bytes) in 9.9035s internal wall time, or 1.0097M reads/s. This is a single
diagnostic run, not a replacement for the three-run no-output comparison.

| Instrumented stage | Seconds | Wall share |
|---|---:|---:|
| Read preparation and mapping | 6.9287 | 69.96% |
| Residual/overlapped pipeline | 2.4396 | 24.63% |
| Index load | 0.2251 | 2.27% |
| Assignment/top-gene record construction | 0.2008 | 2.03% |
| Output write | 0.0740 | 0.75% |
| Time blocked waiting for FASTQ batches | 0.0353 | 0.36% |

The residual is wall time not owned by an explicit timer. It includes background
FASTQ decompression, synchronization, and orchestration; `fastq_load_seconds`
measures only foreground blocking because decompression is pipelined. External
`/usr/bin/time` reported 77.30s user CPU, 0.74s system CPU, 782% mean CPU, and
3.970GiB max RSS.

Compared with the three-run no-output p27 median, real output added 0.541s wall
and 0.124s to the mapping stage (top-gene collection); actual file writing was
only 0.074s. Molecule-first integration may improve downstream count runtime and
intermediate size, but it is not expected to be a large mapper-only speed lever.

## Mapping CPU profile

A user-space cycle profile of the retained p27 mapper was captured over the
same 10M-pair, no-output workload. The run lasted 9.335s and collected 38K
samples (approximately 398.8 billion cycles) with zero lost samples. Percentages
below are self cycles across the whole process, so they include the pipelined
FASTQ reader as well as mapper workers.

| Function or subsystem | Self cycles | Finding |
|---|---:|---|
| `MmapEcIndex::lookup_v4` | 27.88% | Memory latency dominates. In its local annotation, the first random prefix-directory load accounts for 13.71%, and the first packed-suffix load for 52.14%. |
| `ClosedSyncmerIter::next` | 15.98% | The iterator recomputes masks on every yielded-syncmer call and uses modulo 33 for its ring. The mixer itself has 3.81% inclusive process cost. |
| `is_low_complexity_read` | 6.89% | This is a separate full scan of R2 before mapping. |
| `encode_acgt` | 6.76% | This is another full scan plus encoded-buffer construction before syncmer traversal. |
| Dense EC target accumulation | 5.87% | Work after successful index lookups; retained direct gene resolution already avoids an intermediate allocation. |
| `EcScoreScratch::new` | 5.21% | Rayon `map_init` constructs and zero-fills dense score/epoch arrays repeatedly for each batch/job rather than retaining exactly one scratch per worker. Inclusive cost is 6.43%. |
| `EcScoreScratch::add` | 3.81% | Per-hit dense score updates. |
| R1 barcode/UMI parsing | 3.26% | Produces two packed values and also allocates uppercase `String`s for every read. |
| gzip inflate | about 3.6% inclusive | Not a primary mapper bottleneck. |

Allocator cost is split among libc symbols (`_int_malloc` 2.33% self,
`_int_free` 1.08%, `malloc_consolidate` 1.07%, `cfree` 1.04%, plus smaller
entries). Some of that should disappear with persistent score scratch and
allocation-free barcode/UMI handling.

### New optimization order

1. Retain one `EcScoreScratch` per Rayon worker across batches. This is exact,
   localized, and directly targets at least 5.21% self cycles plus allocator
   overhead. Accept only if three-run wall throughput improves by at least 3%
   with byte-identical assignment output and identical rates.
2. Strength-reduce the closed-syncmer iterator without changing selection:
   validate and precompute masks/reverse-complement shift at construction, and
   use a 32-slot ring with `& 31` because `k-s+1 <= 32`. Accept only with
   exhaustive differential tests and at least 2% three-run wall improvement.
3. Fuse ASCII classification, low-complexity statistics, encoding, and closed
   syncmer traversal into one R2 pass. The separate low-complexity and encoding
   functions alone consume 13.65% self cycles, making this the largest
   non-lookup opportunity, but it is a larger correctness-sensitive change.
4. Only after those changes, test a bounded per-read two-phase lookup: collect
   the small set of selected k-mers, issue software prefetches for their prefix
   directory and suffix locations, then resolve/accumulate them. The annotation
   shows that comparison logic is not the problem; any new lookup design must
   reduce cache/TLB stalls and beat the p27 baseline, not merely shorten the
   binary-search loop. This is distinct from the previously failed
   cross-read/locality batching paths.
5. Pack barcode and UMI once without per-read `String` creation when assignments
   are skipped or when the output path can consume packed values directly.

These percentages are not additive speedup promises: removing one bottleneck
changes the denominator and can expose another. A reasonable experiment target
for items 1-3 together is a 12-20% mapper wall-time reduction; anything larger
requires changing the 27.88% random-lookup cost.

## Selector and format

The retained selector mixes canonical rolling s-mer encodings and selects a
k-mer when the minimum order value occurs at either closed-syncmer endpoint. It
is sequence- and strand-intrinsic, single-pass, fixed-queue, and allocation-free.
The v4 reserved fields store a `SYNH` marker and `s`; loaders reject the earlier
lexicographic `SYNC` prototype instead of silently changing query semantics.

Sampled indexes ignore absolute `--sparse-probe-stride`; applying both policies
would compound sampling and was not accepted. Only read-at-a-time traversal is
supported for a sampled index.

## Rejected and superseded variants

| Variant | Throughput | Countable | Disposition |
|---|---:|---:|---|
| Lexicographic `s=11` | 1.156757M/s | 76.91349% | Quality loss; marker rejected. |
| Lexicographic `s=15` | 1.019976M/s | 77.80861% | Slower and lower quality. |
| Lexicographic `s=21` | 0.843972M/s | 79.20886% | Quality up, speed down sharply. |
| Hashed `s=11` | 1.207670M/s | 76.02365% | Fast/small endpoint with unacceptable quality loss. |
| Hashed `s=21` | 0.912297M/s | 78.74331% | Quality clears baseline; speed does not. |
| Lex `s=11` then lex `s=21` unmapped rescue | 0.996078M/s | 77.91391% | Two-index RSS and excess ambiguity; code removed. |
| Hashed `s=20`, prefix24 | 0.995293M/s median | 78.50747% | Superseded by p27. |
| Hashed `s=20`, prefix26 | 1.041952M/s median | 78.50747% | Superseded by p27. |

Do not repeat these paths without satisfying the revisit gate in
`docs/failed_performance_experiments.md`.

## Reproduction

Build the retained index:

```bash
CARGO_BUILD_JOBS=12 RAYON_NUM_THREADS=12 taskset -c 0-11 \
  target/release/constellation index-ec \
  --transcripts benchdata/reference/ensembl93/simpleaf_splici_r91_k31/ref/roers_ref.fa \
  --t2g-map benchdata/reference/ensembl93/simpleaf_splici_r91_k31/ref/t2g_3col.tsv \
  --k 31 --syncmer-s 20 --format mmap \
  --out /tmp/constellation_syncmer_hashed_s20_p27.ecidx
```

Run one no-output measurement:

```bash
RAYON_NUM_THREADS=12 taskset -c 0-11 \
  target/release/constellation map \
  --index /tmp/constellation_syncmer_hashed_s20_p27.ecidx \
  --r1 /tmp/constellation_audit_pbmc10m_R1.fastq.gz \
  --r2 /tmp/constellation_audit_pbmc10m_R2.fastq.gz \
  --batch-size 65536 --skip-assignments \
  --emit-metrics /tmp/constellation_syncmer_p27.metrics.json \
  --out /tmp/constellation_syncmer_p27.cstrad
```

Capture the cycle profile (requires host perf permission):

```bash
sudo perf record -e cycles:u -F 499 -g --call-graph dwarf,16384 \
  -o /tmp/constellation_p27_mapping.perf.data -- \
  env RAYON_NUM_THREADS=12 taskset -c 0-11 \
  target/release/constellation map \
  --index /tmp/constellation_syncmer_hashed_s20_p27.ecidx \
  --r1 /tmp/constellation_audit_pbmc10m_R1.fastq.gz \
  --r2 /tmp/constellation_audit_pbmc10m_R2.fastq.gz \
  --batch-size 65536 --skip-assignments \
  --out /tmp/constellation_p27_mapping.cstrad
```

## Raw artifact labels

- `resolved_dense_v4_r{1,2,3}.*` and
  `resolved_hashed_s20_p{26,27}_r{1,2,3}.*` are the final matched sets.
- `resolved_hashed_s20_r{1,2,3}.*` is the prefix24 set after direct resolution.
- `syncmer_hashed_s{11,20,21}_*` contains hashed-density build/warmup trials.
- `syncmer_s{11,15,21}_*` contains the superseded lexicographic trials.
- `syncmer_s11_rescue_s21_*` is the rejected two-index cascade.
- `dense_v4_perf_stat.txt` and `syncmer_hashed_s20_perf_stat.txt` are valid
  three-repeat sudo `perf stat -d -d -d` reports for dense and pre-resolver p24.
- `syncmer_hashed_s20.perf.data` is zero bytes because `perf record` was blocked
  by `perf_event_paranoid=4`; it is retained only as a clearly labeled failure.
- `current_p27_rad_breakdown.*` is the single uncompressed-output stage-timing
  run described above; its 369MB output remains in `/tmp` and is not archived.
- `current_p27_mapping.self.txt` and `current_p27_mapping.inclusive.txt` are the
  compact symbol reports for the successful user-space cycle profile.
- `current_p27_mapping.lookup_v4.annotate.txt` and
  `current_p27_mapping.syncmer_next.annotate.txt` are the two hot-function
  disassembly annotations used for the conclusions above.
- `raw/perf.data.current_p27_mapping` is the successful 606MB raw capture. It is
  deliberately ignored by Git but retained locally; the original also remains
  at `/tmp/constellation_p27_mapping.perf.data`.

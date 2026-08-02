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

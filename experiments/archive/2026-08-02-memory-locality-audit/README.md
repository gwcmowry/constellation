# 2026-08-02 memory-locality audit

This archive preserves the raw JSON metrics and `/usr/bin/time -v` reports from
the August 2026 stocktake. No experiment logs were deleted. Large FASTQs and
indexes remain outside git in `/tmp`; the durable 152GiB historical corpus remains
under the git-ignored `benchdata/` directory.

## Controls

- Host: Ryzen 9 7900X, 64GB DDR5.
- Workload: 10,000,000 PBMC 10k v3 paired reads.
- Index: human splici v4 prefix24 mmap EC index, approximately 17GB.
- Batch size: 65,536.
- CPU cap: `taskset -c 0-11`, `RAYON_NUM_THREADS=12`; build jobs also capped at 12.
- Headline comparisons: warm page cache, no assignment output, identical mapping
  configuration unless a row names the changed parameter.

## Valid results

| Experiment | Reads/s | Countable | Ambiguous | Unmapped | Status |
|---|---:|---:|---:|---:|---|
| v4 read-at-a-time baseline | 1.035M | 78.54% | 9.27% | 11.46% | Reference (`12t_v4_read_serial`) |
| v4 seed-batched | 0.609M | 78.54% | 9.27% | 11.46% | Rejected |
| v4 prefix-tiled | 0.695M | 78.54% | 9.27% | 11.46% | Rejected |
| v5 packed layout, two warm runs | 1.057–1.077M | 78.54% | 9.27% | 11.46% | Rejected: roughly 20.7GB RSS and larger index |
| v6 blocked layout, best warm run | 0.997M | 78.54% | 9.27% | 11.46% | Rejected |
| v7 deduplicated layout, best warm run | 0.989M | 78.54% | 9.27% | 11.46% | Rejected |
| stride 16 | 1.200M | 76.31% | 10.02% | 12.95% | Rejected quality tradeoff |
| stride 24 | 1.276M | 74.94% | 10.19% | 14.15% | Rejected quality tradeoff |
| stride 16 + rescue 12 | 1.098M | 77.34% | 10.92% | 11.01% | Rejected |
| stride 24 + rescue 12 | 1.173M | 76.38% | 11.43% | 11.46% | Rejected |
| dense epoch scores + cheap preprocessing, three warm runs | 1.044–1.052M | 78.54% | 9.27% | 11.46% | Accepted, modest 1–2% gain |
| per-worker direct-mapped lookup cache | 1.003M | 78.54% | 9.27% | 11.46% | Rejected; code removed |
| file mmap `MADV_HUGEPAGE` | 1.015M | 78.54% | 9.27% | 11.46% | Rejected; code removed |

The accepted dense-score path used 17,642,676–17,644,464 KiB max RSS in its two
recorded warm repetitions, essentially unchanged from the 17,578,156 KiB v4
baseline. Its zstd-output run reached 0.987M reads/s versus 0.948M in the prior
single run. The larger output-path delta requires more repetitions before it is
treated as stable.

The before/after compressed outputs were byte-identical:

```text
7829b876ea967ee7ebc270446e0016acf7b243723eb800810099e0760c92a5e5  before
7829b876ea967ee7ebc270446e0016acf7b243723eb800810099e0760c92a5e5  after
```

## Invalid, incomplete, and diagnostic files

- `constellation_audit_v4_seed_batched.*` and
  `constellation_audit_v4_prefix_tiled.*` were launched concurrently. They are
  retained but invalid for performance comparison.
- `constellation_audit_v4_read_at_a_time.*` predates the enforced 12-thread cap
  and is not comparable to the controlled results.
- `constellation_audit_12t_v4_r1.*` is a valid cold-cache diagnostic: 28.47s,
  0.351M reads/s, and about 3.0 million major faults. It is not a warm baseline.
- `constellation_audit_build_v8.time.txt` is empty because the build was stopped.
  It is an incomplete attempt, not a performance result.
- `constellation_audit_failed_perf.data` is empty. Both cargo-flamegraph and perf
  hardware-counter collection were blocked by host `perf_event_paranoid=4`, even
  when attempted with elevation. No profile can be reconstructed from this file.
- `after_lookup_cache_r1.*` and `after_hugepage_hint_r1.*` are retained negative
  results even though their code was removed.

## Seed reuse trace

A separate 1M-read query trace counted 5,996,875 selected seed occurrences and
3,384,891 unique keys, a 1.77166 reuse factor. Occurrence share was 5.196% for
the top 1,000 keys, 17.380% for the top 10,000, 35.122% for the top 100,000, and
60.231% for the top 1,000,000. Of unique keys, 83.27% appeared only once.

This trace motivated the rejected direct-mapped cache and now supports only the
more selective shared immutable LFU experiment described in
`docs/performance_roadmap.md`.

## File naming

- `*.metrics.json`: Constellation internal metrics.
- `*.time.txt`: `/usr/bin/time -v` process resource report and exact command.
- `constellation_audit_convert_v*.time.txt`: index conversion resource report.
- `after_dense_scratch_*`: accepted implementation measurements.
- `after_lookup_cache_*`, `after_hugepage_hint_*`: rejected new locality tests.

The rejection rationale and revisit rules live at
`docs/failed_performance_experiments.md`.

# PtrHash exact EC-index experiment (2026-08-02)

## Premise

Test a global, non-minimal PtrHash static perfect hash over the retained hashed
closed-syncmer key set (`k=31,s=20`). Each slot stores the full 62-bit canonical
k-mer for exact negative-query rejection and the direct EC target payload.

This is not the stopped v8 prototype. V8 retained a 24-bit prefix directory,
built an independent displacement table inside every prefix bucket, used 1.23x
slots, and never completed a valid production build. PtrHash uses one global,
modern fixed-pilot MPHF designed for billion-scale integer and k-mer keys, has a
streaming-prefetch query path, and requires roughly 2.4 bits per indexed key for
the hash function itself.

## Predeclared gates

- All builds and runs use at most 12 threads and CPUs 0-11.
- Exhaustively compare every indexed positive key with the source index and test
  a large deterministic set of negative keys; no probabilistic membership test
  is accepted.
- Mapping assignment rates and diagnostics must match the retained p27 index
  exactly; compare emitted assignments byte-for-byte when practical.
- Re-run a warm retained-p27 baseline and the candidate three times each on the
  same 10M PBMC workload, batch size 65,536, with assignments disabled.
- Accept only for at least 10% higher median end-to-end throughput, index size no
  larger than 5GB, max mapping RSS no larger than 6GiB, and a practical build
  time and memory footprint on the documented workstation.
- If it fails, remove the production code and dependency, archive the raw logs
  and this result, and add the path to `docs/failed_performance_experiments.md`.

## Result: retained

The candidate passed every predeclared gate.

| Metric | Retained p27 | PtrHash v9 | Change |
|---|---:|---:|---:|
| Median end-to-end wall time, 10M reads | 9.147 s | 7.902 s | -13.6% |
| Median mapping time | 6.564 s | 5.428 s | -17.3% |
| Median throughput | 1.0933M reads/s | 1.2655M reads/s | +15.8% |
| Median max mapping RSS | 3.96 GiB | 4.36 GiB | +10.1% |
| Index size | 3,538,492,341 B | 4,427,647,629 B | +25.1% |

The audited build took 61.50 seconds and peaked at 10,316,360 KiB RSS. It
exhaustively re-queried every retained positive key and payload, then checked
one million deterministic keys confirmed absent by the source v4 index. The
10M mapping runs had no differences in any non-timing metric. Full assignment
output was also byte-identical (`cmp`), with both RAD files exactly 369,461,060
bytes. All 76 workspace tests passed.

The three-run `perf stat -d -d -d` profile measured 1.28 IPC, 5.14% L1D load
misses, and 16.02% dTLB load misses. The archived older-format profiles had
35.85--37.37% dTLB load misses, which is consistent with the intended reduction
in directory/TLB pressure, but is contextual rather than a controlled p27
counter comparison. A `perf record` attempt produced no samples because the
host kept `kernel.perf_event_paranoid=4` even for the tool's elevated run; the
failure and rerun command are recorded in `perf-record-failed.txt`.

## Interpretation and next work

This result validates the global fixed-pilot/static-perfect-hash direction and
is materially different from the stopped v8 per-prefix displacement table. The
tradeoff is a 25% larger index and about 0.4 GiB more mapping RSS. The remaining
16% dTLB miss rate is still substantial. A plausible next experiment is to use
PtrHash's streaming/prefetch lookup API over small per-worker seed batches, but
only if it preserves read-local scoring and is tested as a separate gated
experiment; the archived global seed-batching attempts must not be repeated.

Raw build, benchmark, metrics, assignment timing, and counter files are stored
beside this README.

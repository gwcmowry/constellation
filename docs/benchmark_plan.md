# Benchmark Plan

Use the gene-EC mapper for current performance work. The older positional compact-index modes are legacy and should not be used for headline numbers. Read `docs/failed_performance_experiments.md` before proposing a variant. All CPU experiments and builds are capped at 12 threads/cores and performance runs must execute serially, not concurrently.

## Baseline Run

```bash
CARGO_BUILD_JOBS=12 RAYON_NUM_THREADS=12 taskset -c 0-11 \
  cargo run --release -p constellation-cli -- index-ec \
  --transcripts /path/to/reference.fa \
  --t2g-map /path/to/t2g.tsv \
  --k 31 \
  --syncmer-s 20 \
  --format mmap \
  --out /tmp/reference.k31.syncmer-s20.mmap.ecidx

RAYON_NUM_THREADS=12 taskset -c 0-11 \
  target/release/constellation map \
  --index /tmp/reference.k31.syncmer-s20.mmap.ecidx \
  --r1 /path/to/R1.fastq.gz \
  --r2 /path/to/R2.fastq.gz \
  --batch-size 65536 \
  --output-format gene-ec-rad \
  --output-compression zstd \
  --zstd-level 3 \
  --emit-metrics /tmp/constellation.metrics.json \
  --out /tmp/constellation.cstrad.zst
```

## Metrics

Report these for every benchmark:

```text
wall time
internal wall time
reads/sec
max RSS
mapping_seconds
assignment_seconds
write_seconds
unique_gene_rate
ambiguous_gene_rate
unmapped_rate
output size
index size and build time
warm/cold cache state
major page faults
```

Run at least one warmup followed by three serial measured repetitions. Preserve
each raw metrics JSON and `/usr/bin/time -v` report in a dated directory under
`experiments/archive/`. For exact changes, compare compressed output hashes or a
stronger semantic oracle before accepting speed results.

For profiling runs, also capture:

```bash
perf stat -d -d -d -- target/release/constellation map ...
perf record -F 999 -g -- target/release/constellation map ...
perf report --stdio --no-children --sort comm,dso,symbol
```

The current bottleneck to watch is the total cost of sampled EC lookup and gene
vote accumulation. The smaller working set lowers page faults and dTLB miss rate,
but sparse densities that require too many probes can still increase absolute
dTLB misses and instruction count. Zstd compression is currently below 1% of
cycle samples and should not be treated as the primary optimization target.

The current host has `perf_event_paranoid=4`, so cargo-flamegraph and perf hardware
counters are blocked until host policy changes. Preserve failed/empty profiler
artifacts with an explicit label; never present them as profiles.

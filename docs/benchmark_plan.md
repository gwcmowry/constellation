# Benchmark Plan

Use the gene-EC mapper for current performance work. The older positional compact-index modes are legacy and should not be used for headline numbers.

## Baseline Run

```bash
cargo run --release -p constellation-cli -- index-ec \
  --transcripts /path/to/reference.fa \
  --t2g-map /path/to/t2g.tsv \
  --k 31 \
  --format mmap \
  --out /tmp/reference.prefix24.mmap.ecidx

cargo run --release -p constellation-cli -- map \
  --index /tmp/reference.prefix24.mmap.ecidx \
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
```

For profiling runs, also capture:

```bash
perf stat -d -d -d -- target/release/constellation map ...
perf record -F 999 -g -- target/release/constellation map ...
perf report --stdio --no-children --sort comm,dso,symbol
```

The current bottleneck to watch is `LoadedEcIndex::lookup`, especially dTLB-load-miss concentration. Zstd compression is currently below 1% of cycle samples and should not be treated as the primary optimization target.

# Benchmark Plan

Initial benchmark modes:

```text
read-at-a-time baseline
sketch-bucketed reads
sketch-bucketed + candidate-locus regrouping
candidate-locus regrouping + SIMD scoring
```

Metrics to track:

```text
mode
score_mode
reads/sec
bases/sec
seed lookups/read
seed candidates considered/read
candidate hits/read
candidate buckets/read
postings skipped due to frequency
scored candidates/read
mean/median/max candidate bucket size
unique_gene rate
ambiguous_gene rate
low_complexity rate
low_quality rate
unmapped rate
L1 miss rate if perf stat input is supplied
LLC miss rate if perf stat input is supplied
```

Current smoke benchmark flow:

```bash
cargo run -p constellation-cli -- index \
  --transcripts tests/tiny_transcriptome.fa \
  --k 15 \
  --out /tmp/tiny.cbidx

cargo run -p constellation-cli -- simulate \
  --transcripts tests/tiny_transcriptome.fa \
  --num-reads 1000 \
  --read-len 40 \
  --out-prefix /tmp/sim

cargo run -p constellation-cli -- map \
  --index /tmp/tiny.cbidx \
  --r1 /tmp/sim_R1.fastq \
  --r2 /tmp/sim_R2.fastq \
  --mode candidate-locus \
  --score-mode scalar \
  --emit-metrics /tmp/metrics.json \
  --out /tmp/assignments.tsv

cargo run -p constellation-cli -- bench-report \
  --assignments /tmp/assignments.tsv \
  --truth /tmp/sim_truth.tsv \
  --index /tmp/tiny.cbidx \
  --metrics /tmp/metrics.json

cargo run -p constellation-cli -- count \
  --assignments /tmp/assignments.tsv \
  --index /tmp/tiny.cbidx \
  --out-prefix /tmp/counts
```

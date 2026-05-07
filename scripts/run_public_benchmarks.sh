#!/usr/bin/env bash
set -euo pipefail

if [[ $# -lt 5 ]]; then
  cat >&2 <<'USAGE'
usage:
  scripts/run_public_benchmarks.sh TRANSCRIPTS_FA R1_FASTQ R2_FASTQ OUT_DIR K

This smoke wrapper expects the public dataset files to already be downloaded.
It builds a transcript index, maps the supplied 10x-style FASTQs, and emits
assignments plus metrics for later comparison against external baselines.
USAGE
  exit 2
fi

transcripts=$1
r1=$2
r2=$3
out_dir=$4
k=$5

mkdir -p "$out_dir"

cargo run -p constellation-cli -- index \
  --transcripts "$transcripts" \
  --k "$k" \
  --out "$out_dir/transcripts.cbidx"

cargo run -p constellation-cli -- map \
  --index "$out_dir/transcripts.cbidx" \
  --r1 "$r1" \
  --r2 "$r2" \
  --chemistry tenx-3p-v3 \
  --mode candidate-locus \
  --score-mode scalar \
  --emit-metrics "$out_dir/metrics.json" \
  --out "$out_dir/assignments.tsv"

cargo run -p constellation-cli -- count \
  --assignments "$out_dir/assignments.tsv" \
  --index "$out_dir/transcripts.cbidx" \
  --out-prefix "$out_dir/counts"

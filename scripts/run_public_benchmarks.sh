#!/usr/bin/env bash
set -euo pipefail

if [[ $# -lt 5 ]]; then
  cat >&2 <<'USAGE'
usage:
  scripts/run_public_benchmarks.sh TRANSCRIPTS_FA R1_FASTQ R2_FASTQ OUT_DIR K [T2G_MAP]

This smoke wrapper expects the public dataset files to already be downloaded.
It builds the current gene-EC index, maps the supplied 10x-style FASTQs, and
emits Constellation RAD-like output plus metrics for external comparison.
USAGE
  exit 2
fi

transcripts=$1
r1=$2
r2=$3
out_dir=$4
k=$5
t2g=${6:-}

mkdir -p "$out_dir"

index_args=()
if [[ -n "$t2g" ]]; then
  index_args+=(--t2g-map "$t2g")
fi

cargo run --release -p constellation-cli -- index-ec \
  --transcripts "$transcripts" \
  "${index_args[@]}" \
  --k "$k" \
  --format mmap \
  --out "$out_dir/reference.ecidx"

cargo run --release -p constellation-cli -- map \
  --index "$out_dir/reference.ecidx" \
  --r1 "$r1" \
  --r2 "$r2" \
  --chemistry tenx-3p-v3 \
  --batch-size 65536 \
  --output-format gene-ec-rad \
  --output-compression zstd \
  --zstd-level 3 \
  --emit-metrics "$out_dir/metrics.json" \
  --out "$out_dir/assignments.cstrad.zst"

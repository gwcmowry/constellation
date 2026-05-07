#!/usr/bin/env bash
set -euo pipefail

repo_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
export ALEVIN_FRY_HOME="${ALEVIN_FRY_HOME:-$repo_dir/.af-home}"

simpleaf="${SIMPLEAF:-$repo_dir/.conda-af/bin/simpleaf}"
genome="${GENOME:-$repo_dir/benchdata/reference/ensembl93/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz}"
gtf="${GTF:-$repo_dir/benchdata/reference/ensembl93/Homo_sapiens.GRCh38.93.gtf}"
r1="${R1:-$repo_dir/benchdata/pbmc_10k_v3/pbmc_10k_v3_100k_L001_R1.fastq.gz}"
r2="${R2:-$repo_dir/benchdata/pbmc_10k_v3/pbmc_10k_v3_100k_L001_R2.fastq.gz}"
threads="${THREADS:-16}"
index_dir="${INDEX_DIR:-$repo_dir/benchdata/reference/ensembl93/simpleaf_splici_r91_k31}"
out_dir="${OUT_DIR:-$repo_dir/benchdata/pbmc_10k_v3/simpleaf_splici_r91_k31_100k_cr_like}"
work_dir="${WORK_DIR:-/tmp/simpleaf_splici_r91_work}"

if [[ ! -x "$simpleaf" ]]; then
  echo "simpleaf not found at $simpleaf" >&2
  exit 1
fi

mkdir -p "$ALEVIN_FRY_HOME"

if [[ ! -e "$ALEVIN_FRY_HOME/config.json" ]]; then
  "$simpleaf" set-paths \
    --piscem "$repo_dir/.conda-af/bin/piscem" \
    --alevin-fry "$repo_dir/.conda-af/bin/alevin-fry"
fi

if [[ ! -d "$index_dir/index" ]]; then
  "$simpleaf" index \
    --fasta "$genome" \
    --gtf "$gtf" \
    --rlen 91 \
    --threads "$threads" \
    --output "$index_dir" \
    --work-dir "$work_dir"
fi

"$simpleaf" quant \
  --index "$index_dir/index/piscem_idx" \
  --reads1 "$r1" \
  --reads2 "$r2" \
  --chemistry 10xv3 \
  --expected-ori fw \
  --threads "$threads" \
  --expect-cells 10000 \
  --resolution cr-like \
  --t2g-map "$index_dir/ref/t2g_3col.tsv" \
  --output "$out_dir"

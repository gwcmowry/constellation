#!/usr/bin/env bash
set -euo pipefail

if [[ $# -lt 5 ]]; then
  cat >&2 <<'USAGE'
usage:
  scripts/benchmark_gene_ec_perf.sh EC_INDEX R1_FASTQ R2_FASTQ OUT_DIR LABEL [CONSTELLATION_BIN]

Runs the current gene-EC mapper with internal metrics and perf stat counters.
The output directory receives:
  LABEL.metrics.json
  LABEL.assignments.tsv
  LABEL.perf.txt
USAGE
  exit 2
fi

index=$1
r1=$2
r2=$3
out_dir=$4
label=$5
bin=${6:-target/release/constellation}

mkdir -p "$out_dir"

metrics="$out_dir/$label.metrics.json"
assignments="$out_dir/$label.assignments.tsv"
perf_out="$out_dir/$label.perf.txt"
skip_args=()
if [[ "${SKIP_ASSIGNMENTS:-0}" == "1" ]]; then
  skip_args+=(--skip-assignments)
fi

perf stat -d -d -d -o "$perf_out" \
  "$bin" map \
    --mode gene-ec \
    --index "$index" \
    --r1 "$r1" \
    --r2 "$r2" \
    --chemistry tenx-3p-v3 \
    --search-reverse-complement \
    --trim-tso \
    --min-mean-quality 10 \
    "${skip_args[@]}" \
    --emit-metrics "$metrics" \
    --out "$assignments"

python3 - "$metrics" "$perf_out" <<'PY'
import json
import re
import sys
from pathlib import Path

metrics = json.loads(Path(sys.argv[1]).read_text())
perf = Path(sys.argv[2]).read_text()

values = {}
for line in perf.splitlines():
    stripped = line.strip()
    match = re.match(r"([0-9,]+)\s+([A-Za-z0-9_.:-]+)", stripped)
    if match:
        values[match.group(2)] = int(match.group(1).replace(",", ""))

reads = max(metrics["num_reads"], 1)
cycles = values.get("cycles", 0)
instructions = values.get("instructions", 0)
dtlb = values.get("dTLB-load-misses", 0)
page_faults = values.get("page-faults", 0)

print(f"num_reads\t{metrics['num_reads']}")
print(f"wall_seconds\t{metrics['wall_seconds']:.6f}")
print(f"index_load_seconds\t{metrics['index_load_seconds']:.6f}")
print(f"fastq_load_seconds\t{metrics['fastq_load_seconds']:.6f}")
print(f"preprocess_seconds\t{metrics['preprocess_seconds']:.6f}")
print(f"mapping_seconds\t{metrics['mapping_seconds']:.6f}")
print(f"write_seconds\t{metrics['write_seconds']:.6f}")
print(f"reads_per_second_wall\t{metrics['num_reads'] / metrics['wall_seconds']:.2f}")
print(f"reads_per_second_mapping\t{metrics['num_reads'] / metrics['mapping_seconds']:.2f}")
if cycles:
    print(f"cycles_per_read\t{cycles / reads:.2f}")
if instructions and cycles:
    print(f"ipc\t{instructions / cycles:.4f}")
if dtlb:
    print(f"dtlb_load_misses_per_read\t{dtlb / reads:.2f}")
if page_faults:
    print(f"page_faults_per_read\t{page_faults / reads:.4f}")
print(f"gene_countable_rate\t{metrics['gene_countable_rate']:.6f}")
print(f"ambiguous_gene_rate\t{metrics['ambiguous_gene_rate']:.6f}")
print(f"unmapped_rate\t{metrics['unmapped_rate']:.6f}")
PY

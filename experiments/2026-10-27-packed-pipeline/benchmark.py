"""Serial full-pipeline benchmark, with exact aggregate-metric checks."""
import fcntl
import json
import os
from pathlib import Path
import subprocess
import sys

root = Path(__file__).resolve().parents[2]
os.chdir(root)
data = Path("benchdata/throughput-2026-10-27")
raw = Path(__file__).parent / "raw"
label, binary, mode = sys.argv[1:4]
repeats = int(sys.argv[4]) if len(sys.argv) > 4 else 3
lock = (data / "benchmark.lock").open("w")
fcntl.flock(lock, fcntl.LOCK_EX)
reference = json.loads((root / "experiments/2026-10-27-read-throughput/raw/baseline-1.metrics.json").read_text())


def semantic(m):
    return {k: v for k, v in m.items() if not (k.endswith("seconds") or k.endswith("per_second"))}


for repeat in ["warmup", *range(1, repeats + 1)]:
    prefix = raw / f"{label}-{mode}-{repeat}"
    output = data / f"{label}.{mode}"
    cmd = ["/usr/bin/time", "-v", "-o", f"{prefix}.time", "taskset", "-c", "0-11",
           binary, "map", "--index", str(data / "reference.ptrhash.ecidx"),
           "--r1", str(data / "pbmc10m_R1.fastq.gz"), "--r2", str(data / "pbmc10m_R2.fastq.gz"),
           "--batch-size", "65536", "--emit-metrics", f"{prefix}.metrics.json",
           "--out", str(output), "--output-format", "tsv" if mode == "tsv" else "gene-ec-rad"]
    if mode == "none":
        cmd += ["--skip-assignments"]
    if mode == "zstd":
        cmd += ["--output-compression", "zstd", "--zstd-level", "3"]
    with open(f"{prefix}.log", "w") as log:
        subprocess.run(cmd, env=dict(os.environ, RAYON_NUM_THREADS="12"), stdout=log, stderr=log, check=True)
    metrics = json.loads(Path(f"{prefix}.metrics.json").read_text())
    actual, expected = semantic(metrics), semantic(reference)
    assert actual == expected, {k: (actual[k], expected[k]) for k in actual if actual[k] != expected[k]}
    print(json.dumps({"run": f"{label}-{mode}-{repeat}", **{k: metrics[k] for k in (
        "wall_seconds", "reads_per_second", "mapping_seconds", "fastq_load_seconds", "write_seconds")}}), flush=True)

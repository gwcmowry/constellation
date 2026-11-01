"""Run serial, warm-cache pipeline ablations; preserve timings and metrics."""
import json
import fcntl
import os
from pathlib import Path
import subprocess
import sys

root = Path(__file__).resolve().parents[2]
os.chdir(root)
data = Path("benchdata/throughput-2026-10-27")
raw = Path(__file__).parent / "raw"
lock = (data / "benchmark.lock").open("w")
fcntl.flock(lock, fcntl.LOCK_EX)
label, binary, allocator, compression = sys.argv[1:5]
env = dict(os.environ, RAYON_NUM_THREADS="12")
if allocator == "tcmalloc":
    env["LD_PRELOAD"] = "/lib/x86_64-linux-gnu/libtcmalloc_minimal.so.4"
for repeat in ("warmup", "1", "2", "3"):
    prefix = raw / f"{label}-{repeat}"
    suffix = ".fastq.gz" if compression == "gzip" else ".fastq"
    cmd = [
        "/usr/bin/time", "-v", "-o", f"{prefix}.time", "taskset", "-c", "0-11",
        binary, "map", "--index", str(data / "reference.ptrhash.ecidx"),
        "--r1", str(data / f"pbmc10m_R1{suffix}"),
        "--r2", str(data / f"pbmc10m_R2{suffix}"),
        "--batch-size", "65536", "--skip-assignments", "--emit-metrics",
        f"{prefix}.metrics.json", "--out", str(data / "unused.rad"),
    ]
    with open(f"{prefix}.log", "w") as log:
        subprocess.run(cmd, env=env, stdout=log, stderr=log, check=True)
    metrics = json.loads(Path(f"{prefix}.metrics.json").read_text())
    print(json.dumps({"run": f"{label}-{repeat}", **{
        k: metrics[k] for k in ["wall_seconds", "mapping_seconds", "reads_per_second"]
    }, "function_seconds": metrics["function_seconds"]}), flush=True)

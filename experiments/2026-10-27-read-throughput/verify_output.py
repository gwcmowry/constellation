"""Compare complete 10M-read RAD output for baseline and allocator diagnostic."""
import fcntl
import hashlib
import json
import os
from pathlib import Path
import subprocess

root = Path(__file__).resolve().parents[2]
os.chdir(root)
data = Path("benchdata/throughput-2026-10-27")
raw = Path(__file__).parent / "raw"
lock = (data / "benchmark.lock").open("w")
fcntl.flock(lock, fcntl.LOCK_EX)
outputs = []
for label, binary, allocator in [
    ("baseline", str(data / "constellation.baseline"), False),
    ("diagnostic-tcmalloc", str(data / "constellation.diagnostic"), True),
]:
    env = dict(os.environ, RAYON_NUM_THREADS="12")
    if allocator:
        env["LD_PRELOAD"] = "/lib/x86_64-linux-gnu/libtcmalloc_minimal.so.4"
    out = data / f"verify-{label}.rad"
    cmd = [
        "/usr/bin/time", "-v", "-o", str(raw / f"verify-{label}.time"),
        "taskset", "-c", "0-11", binary, "map", "--index",
        str(data / "reference.ptrhash.ecidx"), "--r1", str(data / "pbmc10m_R1.fastq.gz"),
        "--r2", str(data / "pbmc10m_R2.fastq.gz"), "--batch-size", "65536",
        "--output-format", "gene-ec-rad", "--out", str(out),
    ]
    with (raw / f"verify-{label}.log").open("w") as log:
        subprocess.run(cmd, env=env, stdout=log, stderr=log, check=True)
    outputs.append(out)
    print(f"{label} output complete", flush=True)
digest = hashlib.sha256()
with outputs[0].open("rb") as left, outputs[1].open("rb") as right:
    while chunk := left.read(8 << 20):
        assert chunk == right.read(len(chunk)), "RAD output differs"
        digest.update(chunk)
    assert not right.read(1), "Output lengths differ"
result = {"reads": 10_000_000, "byte_identical": True,
          "bytes": outputs[0].stat().st_size, "sha256": digest.hexdigest()}
(raw / "output-verification.json").write_text(json.dumps(result, indent=2) + "\n")
print(json.dumps(result), flush=True)

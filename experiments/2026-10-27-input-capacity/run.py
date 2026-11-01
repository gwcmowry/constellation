"""Run the release-only input diagnostic using the actual production reader."""
import json
import os
from pathlib import Path
import subprocess
import fcntl

root = Path(__file__).resolve().parents[2]
raw = Path(__file__).parent / "raw"
data = root / "benchdata/throughput-2026-10-27"
lock = (data / "benchmark.lock").open("w")
fcntl.flock(lock, fcntl.LOCK_EX)
artifacts = [json.loads(line) for line in (raw / "build.jsonl").read_text().splitlines()]
executable = next(row["executable"] for row in reversed(artifacts)
                  if row.get("reason") == "compiler-artifact" and row.get("executable")
                  and row.get("profile", {}).get("test"))
env = dict(os.environ, RAYON_NUM_THREADS="12",
           CST_INPUT_R1=str(data / "pbmc10m_R1.fastq.gz"),
           CST_INPUT_R2=str(data / "pbmc10m_R2.fastq.gz"),
           CST_INPUT_JSON=str(raw / "measurements.json"))
cmd = ["/usr/bin/time", "-v", "-o", str(raw / "diagnostic.time"),
       "taskset", "-c", "0-11", executable,
       "input_capacity_diagnostic", "--ignored", "--nocapture", "--test-threads=1"]
(raw / "command.json").write_text(json.dumps(cmd, indent=2) + "\n")
with (raw / "diagnostic.log").open("w") as log:
    subprocess.run(cmd, env=env, cwd=root, stdout=log, stderr=log, check=True)
print("Input diagnostic complete")

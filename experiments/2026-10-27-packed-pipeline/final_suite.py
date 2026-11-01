"""Run the retained baseline and candidate serially, then compare complete output."""
import hashlib
import json
from pathlib import Path
import subprocess
import sys

root = Path(__file__).resolve().parents[2]
data = root / "benchdata/throughput-2026-10-27"
experiment = Path(__file__).parent
for label, mode, repeats in [
    ("baseline", "none", 3), ("optimized", "none", 3),
    ("baseline", "zstd", 3), ("optimized", "zstd", 3),
    ("optimized", "rad", 1), ("baseline", "tsv", 1), ("optimized", "tsv", 1),
]:
    binary = data / f"constellation.{label}"
    subprocess.run([sys.executable, str(experiment / "benchmark.py"), label, str(binary), mode, str(repeats)], check=True, cwd=root)


def compare(left, right):
    digest = hashlib.sha256()
    size = 0
    while block := left.read(1 << 20):
        # A pipe's read may return fewer bytes; fill exactly one comparison block.
        other = bytearray()
        while len(other) < len(block):
            piece = right.read(len(block) - len(other))
            assert piece, "candidate output truncated"
            other.extend(piece)
        assert block == other, f"output mismatch near byte {size}"
        size += len(block)
        digest.update(block)
    assert not right.read(1), "candidate output longer than reference"
    return {"byte_identical": True, "bytes": size, "sha256": digest.hexdigest()}


checks = {}
with (data / "verify-baseline.rad").open("rb") as a, (data / "optimized.rad").open("rb") as b:
    checks["rad"] = compare(a, b)
for label in ["baseline", "optimized"]:
    with subprocess.Popen(["zstd", "-q", "-d", "-c", str(data / f"{label}.zstd")], stdout=subprocess.PIPE) as decompressor:
        with (data / "verify-baseline.rad").open("rb") as expected:
            checks[f"{label}_zstd_decoded"] = compare(expected, decompressor.stdout)
        assert decompressor.wait() == 0
with (data / "baseline.tsv").open("rb") as a, (data / "optimized.tsv").open("rb") as b:
    checks["tsv"] = compare(a, b)
(experiment / "raw/output-equality.json").write_text(json.dumps(checks, indent=2) + "\n")
print(json.dumps(checks), flush=True)

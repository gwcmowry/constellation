"""Summarize final serial measurements, including process CPU time and RSS."""
import json
from pathlib import Path
import re
from statistics import median

root = Path(__file__).resolve().parents[2]
folder = Path(__file__).parent / "raw"
summary = {}
for label, mode, n in [(a, b, 3) for a in ["baseline", "optimized"] for b in ["none", "zstd"]] + [
    ("optimized", "rad", 1), ("baseline", "tsv", 1), ("optimized", "tsv", 1)
]:
    metrics = [json.loads((folder / f"{label}-{mode}-{i}.metrics.json").read_text()) for i in range(1, n + 1)]
    times = [(folder / f"{label}-{mode}-{i}.time").read_text() for i in range(1, n + 1)]
    def values(pattern):
        return [float(re.search(pattern, text)[1]) for text in times]
    summary[f"{label}-{mode}"] = {
        "repetitions": n,
        "median_wall_seconds": median(m["wall_seconds"] for m in metrics),
        "wall_range": [min(m["wall_seconds"] for m in metrics), max(m["wall_seconds"] for m in metrics)],
        "median_reads_per_second": median(m["reads_per_second"] for m in metrics),
        "median_max_rss_kib": median(values(r"Maximum resident set size \(kbytes\): (\d+)")),
        "median_cpu_seconds": median(a + b for a, b in zip(values(r"User time \(seconds\): ([\d.]+)"), values(r"System time \(seconds\): ([\d.]+)"))),
        "major_faults": values(r"Major \(requiring I/O\) page faults: (\d+)"),
        "median_mapping_seconds": median(m["mapping_seconds"] for m in metrics),
        "median_input_wait_seconds": median(m["fastq_load_seconds"] for m in metrics),
        "median_write_seconds": median(m["write_seconds"] for m in metrics),
    }
    if mode != "none":
        summary[f"{label}-{mode}"]["output_bytes"] = (root / "benchdata/throughput-2026-10-27" / f"{label}.{mode}").stat().st_size
for mode in ["none", "zstd", "tsv"]:
    summary[f"optimized-{mode}"]["speedup"] = summary[f"baseline-{mode}"]["median_wall_seconds"] / summary[f"optimized-{mode}"]["median_wall_seconds"]
summary["output_equality"] = json.loads((folder / "output-equality.json").read_text())
(folder / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
for label, values in summary.items():
    if label != "output_equality":
        print(label, json.dumps(values))

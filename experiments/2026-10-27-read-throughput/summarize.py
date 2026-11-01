"""Summarize the controlled runs and check all non-timing metrics for equality."""
import json
from pathlib import Path
import re
from statistics import median

raw = Path(__file__).parent / "raw"
baseline = json.loads((raw / "baseline-1.metrics.json").read_text())


def semantic(metrics):
    return {k: v for k, v in metrics.items() if not (
        k.endswith("seconds") or k.endswith("per_second") or k == "function_seconds"
    )}


summary = {}
for label in ["baseline", "instrumented", "tcmalloc", "uncompressed"]:
    rows = [json.loads((raw / f"{label}-{i}.metrics.json").read_text()) for i in (1, 2, 3)]
    for row in rows:
        assert semantic(row) == semantic(baseline), f"Assignment/statistics mismatch in {label}"
    rss = []
    major_faults = []
    for i in (1, 2, 3):
        report = (raw / f"{label}-{i}.time").read_text()
        rss.append(int(re.search(r"Maximum resident set size \(kbytes\): (\d+)", report)[1]))
        major_faults.append(int(re.search(r"Major \(requiring I/O\) page faults: (\d+)", report)[1]))
    summary[label] = {
        "wall_seconds": [r["wall_seconds"] for r in rows],
        "median_wall_seconds": median(r["wall_seconds"] for r in rows),
        "median_reads_per_second": median(r["reads_per_second"] for r in rows),
        "median_max_rss_kib": median(rss),
        "major_faults": major_faults,
        "median_function_seconds": {
            name: median(r["function_seconds"][name] for r in rows)
            for name in rows[0]["function_seconds"]
        },
        "semantic_metrics_equal": True,
    }
historical = Path(__file__).parent.parent / "2026-08-02-ptrhash-index/ptrhash-r1.metrics.json"
assert semantic(json.loads(historical.read_text())) == semantic(baseline)
summary["historical_semantic_metrics_equal"] = True
(raw / "pipeline-summary.json").write_text(json.dumps(summary, indent=2) + "\n")
print(json.dumps(summary, indent=2))

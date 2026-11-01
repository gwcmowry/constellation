"""Summarize input-only measurements and compare with mapper throughput."""
import json
from pathlib import Path
from statistics import median

experiment = Path(__file__).parent
inputs = json.loads((experiment.parent / "2026-10-27-read-throughput/raw/inputs.json").read_text())
rows = json.loads((experiment / "raw/measurements.json").read_text())
result = {}
for row in rows:
    label = row["label"]
    mates = ["R1", "R2"] if label.startswith("paired_") else [label[:2]]
    decoded_bytes = sum(inputs[m]["raw_bytes"] for m in mates)
    if "decode_only" in label:
        assert row["bytes"] == decoded_bytes
    else:
        assert row["records"] == 10_000_000
    seconds = median(row["seconds"])
    result[label] = {
        "median_seconds": seconds,
        "seconds_range": [min(row["seconds"]), max(row["seconds"])],
        "million_records_per_second": 10 / seconds,
        "record_unit": "read pairs" if len(mates) == 2 else "mate records",
        "record_count_assumed_from_fixture": "decode_only" in label,
        "decompressed_input_GB_per_second": decoded_bytes / seconds / 1e9,
    }
pipeline_rate = result["paired_pipeline_gzip"]["million_records_per_second"]
mapper = json.loads((experiment.parent / "2026-10-27-packed-pipeline/raw/summary.json").read_text())
comparison = {
    mode: {"mapper_million_pairs_per_second": mapper[f"optimized-{mode}"]["median_reads_per_second"] / 1e6,
           "isolated_reader_headroom_factor": pipeline_rate / (mapper[f"optimized-{mode}"]["median_reads_per_second"] / 1e6)}
    for mode in ["none", "zstd"]
}
(experiment / "raw/summary.json").write_text(json.dumps({"measurements": result, "comparison": comparison}, indent=2) + "\n")
for label, row in result.items():
    unit = "pairs" if row["record_unit"] == "read pairs" else "records"
    print(f'{label:28s} {row["median_seconds"]:.4f} s  {row["million_records_per_second"]:.3f} M {unit}/s  {row["decompressed_input_GB_per_second"]:.3f} GB/s')
print(json.dumps(comparison, indent=2))

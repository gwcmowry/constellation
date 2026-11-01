"""Recreate the first 10M L001 pairs, with raw and gzip inputs for ablations."""
import gzip
import hashlib
import json
from pathlib import Path

root = Path(__file__).resolve().parents[2]
dest = root / "benchdata/throughput-2026-10-27"
source = root / "benchdata/pbmc_10k_v3/pbmc_10k_v3_fastqs"
manifest = {}
for mate in (1, 2):
    src = source / f"pbmc_10k_v3_S1_L001_R{mate}_001.fastq.gz"
    raw = dest / f"pbmc10m_R{mate}.fastq"
    digest = hashlib.sha256()
    with gzip.open(src, "rb") as inp, raw.open("wb") as out:
        for batch in range(100):
            block = bytearray()
            for _ in range(100_000):
                lines = [inp.readline() for _ in range(4)]
                if not all(lines):
                    raise ValueError("Input has fewer than 10M complete records")
                block.extend(b"".join(lines))
            out.write(block)
            digest.update(block)
    compressed = raw.with_suffix(".fastq.gz")
    with raw.open("rb") as inp, gzip.open(compressed, "wb", compresslevel=1) as out:
        while block := inp.read(8 << 20):
            out.write(block)
    manifest[f"R{mate}"] = {
        "source": str(src.relative_to(root)),
        "reads": 10_000_000,
        "raw_bytes": raw.stat().st_size,
        "raw_sha256": digest.hexdigest(),
        "gzip_level": 1,
        "gzip_bytes": compressed.stat().st_size,
    }
    print(f"R{mate} complete", flush=True)
(Path(__file__).parent / "raw/inputs.json").write_text(json.dumps(manifest, indent=2) + "\n")

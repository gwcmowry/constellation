#!/usr/bin/env python3
import argparse
import gzip
from pathlib import Path


def opener(path, mode):
    return gzip.open(path, mode) if str(path).endswith(".gz") else open(path, mode)


def copy_records(src, dst, records):
    with opener(src, "rt") as inp, gzip.open(dst, "wt") as out:
        for _ in range(records):
            block = [inp.readline() for _ in range(4)]
            if not block[0]:
                break
            if any(line == "" for line in block):
                raise ValueError(f"truncated FASTQ record in {src}")
            out.writelines(block)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--r1", required=True)
    parser.add_argument("--r2", required=True)
    parser.add_argument("--records", type=int, required=True)
    parser.add_argument("--out-prefix", required=True)
    args = parser.parse_args()

    prefix = Path(args.out_prefix)
    copy_records(args.r1, f"{prefix}_R1.fastq.gz", args.records)
    copy_records(args.r2, f"{prefix}_R2.fastq.gz", args.records)


if __name__ == "__main__":
    main()

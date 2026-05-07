#!/usr/bin/env python3
import argparse
import gzip


def opener(path):
    return gzip.open(path, "rt") if str(path).endswith(".gz") else open(path, "rt")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", required=True)
    parser.add_argument("--records", type=int, required=True)
    parser.add_argument("--out", required=True)
    args = parser.parse_args()

    written = 0
    keep = False
    with opener(args.input) as inp, open(args.out, "wt") as out:
        for line in inp:
            if line.startswith(">"):
                if written >= args.records:
                    break
                written += 1
                keep = True
            if keep:
                out.write(line)


if __name__ == "__main__":
    main()

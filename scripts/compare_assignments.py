#!/usr/bin/env python3
import argparse
import csv


def read_assignments(path):
    rows = {}
    with open(path, newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            rows[row["read_id"]] = row
    return rows


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--left", required=True)
    parser.add_argument("--right", required=True)
    args = parser.parse_args()

    left = read_assignments(args.left)
    right = read_assignments(args.right)
    common = sorted(set(left) & set(right), key=int)
    same_type = 0
    same_gene = 0
    for read_id in common:
        same_type += left[read_id]["assignment_type"] == right[read_id]["assignment_type"]
        same_gene += left[read_id]["gene_id"] == right[read_id]["gene_id"]

    n = len(common)
    print(f"common_reads\t{n}")
    print(f"assignment_type_agreement\t{same_type / n if n else 0:.6f}")
    print(f"gene_id_agreement\t{same_gene / n if n else 0:.6f}")


if __name__ == "__main__":
    main()

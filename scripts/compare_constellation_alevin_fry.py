#!/usr/bin/env python3
from __future__ import annotations

import argparse
import math
import re
from collections import defaultdict
from pathlib import Path


def read_constellation_features(path: Path) -> list[str]:
    genes = []
    with path.open() as handle:
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            genes.append(fields[1])
    return genes


def read_alevin_features(path: Path) -> list[str]:
    return [base_gene(line.strip()) for line in path.open() if line.strip()]


def base_gene(gene: str) -> str:
    return re.sub(r"-[ASU]$", "", gene)


def constellation_gene_totals(mtx: Path, features: list[str]) -> dict[str, float]:
    totals: dict[str, float] = defaultdict(float)
    with mtx.open() as handle:
        shape_seen = False
        for line in handle:
            if line.startswith("%"):
                continue
            if not shape_seen:
                shape_seen = True
                continue
            row, _col, value = line.split()[:3]
            totals[features[int(row) - 1]] += float(value)
    return totals


def alevin_gene_totals(mtx: Path, features: list[str]) -> dict[str, float]:
    totals: dict[str, float] = defaultdict(float)
    with mtx.open() as handle:
        shape_seen = False
        for line in handle:
            if line.startswith("%"):
                continue
            if not shape_seen:
                shape_seen = True
                continue
            _row, col, value = line.split()[:3]
            totals[features[int(col) - 1]] += float(value)
    return totals


def pearson(xs: list[float], ys: list[float]) -> float:
    if len(xs) < 2:
        return float("nan")
    mx = sum(xs) / len(xs)
    my = sum(ys) / len(ys)
    vx = sum((x - mx) ** 2 for x in xs)
    vy = sum((y - my) ** 2 for y in ys)
    if vx == 0 or vy == 0:
        return float("nan")
    return sum((x - mx) * (y - my) for x, y in zip(xs, ys)) / math.sqrt(vx * vy)


def ranks(values: list[float]) -> list[float]:
    order = sorted(range(len(values)), key=values.__getitem__)
    out = [0.0] * len(values)
    i = 0
    while i < len(order):
        j = i + 1
        while j < len(order) and values[order[j]] == values[order[i]]:
            j += 1
        rank = (i + j - 1) / 2.0
        for k in range(i, j):
            out[order[k]] = rank
        i = j
    return out


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--constellation-mtx", required=True, type=Path)
    parser.add_argument("--constellation-features", required=True, type=Path)
    parser.add_argument("--alevin-mtx", required=True, type=Path)
    parser.add_argument("--alevin-features", required=True, type=Path)
    args = parser.parse_args()

    constellation_features = read_constellation_features(args.constellation_features)
    alevin_features = read_alevin_features(args.alevin_features)
    constellation = constellation_gene_totals(args.constellation_mtx, constellation_features)
    alevin = alevin_gene_totals(args.alevin_mtx, alevin_features)

    common = sorted(set(constellation) & set(alevin))
    xs = [math.log1p(constellation[gene]) for gene in common]
    ys = [math.log1p(alevin[gene]) for gene in common]

    print(f"constellation_total_umis\t{int(sum(constellation.values()))}")
    print(f"alevin_fry_total_umis\t{int(sum(alevin.values()))}")
    print(f"constellation_nonzero_genes\t{len(constellation)}")
    print(f"alevin_fry_nonzero_genes\t{len(alevin)}")
    print(f"common_nonzero_genes\t{len(common)}")
    print(f"pearson_log_common\t{pearson(xs, ys):.6f}")
    print(f"spearman_log_common\t{pearson(ranks(xs), ranks(ys)):.6f}")
    print("top_constellation\t" + repr(sorted(constellation.items(), key=lambda x: x[1], reverse=True)[:10]))
    print("top_alevin_fry\t" + repr(sorted(alevin.items(), key=lambda x: x[1], reverse=True)[:10]))


if __name__ == "__main__":
    main()

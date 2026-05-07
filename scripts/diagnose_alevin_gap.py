#!/usr/bin/env python3
from __future__ import annotations

import argparse
import collections
import re
from pathlib import Path


COUNTABLE = {"unique_gene", "ambiguous_transcript_same_gene"}


def base_gene(gene: str) -> str:
    return re.sub(r"-[ASU]$", "", gene)


def read_constellation_features(path: Path | None) -> list[str]:
    if path is None:
        return []
    genes = []
    with path.open() as handle:
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            genes.append(fields[1])
    return genes


def read_assignments(path: Path, features: list[str]) -> dict[int, dict[str, str]]:
    assignments = {}
    with path.open() as handle:
        header = next(handle).rstrip("\n").split("\t")
        cols = {name: idx for idx, name in enumerate(header)}
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            read_id = int(fields[cols["read_id"]])
            gene_id = fields[cols["gene_id"]]
            gene = "."
            if gene_id != "." and features:
                idx = int(gene_id)
                if idx < len(features):
                    gene = features[idx]
            assignments[read_id] = {
                "cell_barcode": fields[cols["cell_barcode"]],
                "umi": fields[cols["umi"]],
                "assignment_type": fields[cols["assignment_type"]],
                "gene": gene,
                "candidate_count": fields[cols["candidate_count"]],
                "score": fields[cols["score"]],
                "flags": fields[cols["flags"]],
            }
    return assignments


def read_unmapped_diagnostics(path: Path | None) -> dict[int, str]:
    if path is None:
        return {}
    reasons = {}
    with path.open() as handle:
        header = next(handle).rstrip("\n").split("\t")
        cols = {name: idx for idx, name in enumerate(header)}
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            reasons[int(fields[cols["read_id"]])] = fields[cols["reason"]]
    return reasons


def read_t2g(path: Path) -> dict[str, tuple[str, str]]:
    t2g = {}
    with path.open() as handle:
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            if len(fields) >= 2:
                target_class = fields[2] if len(fields) > 2 else "."
                t2g[fields[0]] = (base_gene(fields[1]), target_class)
    return t2g


def parse_rad_view(path: Path, t2g: dict[str, tuple[str, str]]) -> dict[int, dict[str, object]]:
    reads: dict[int, dict[str, object]] = {}
    with path.open() as handle:
        for line in handle:
            if not line.startswith("ID:"):
                continue
            fields = {}
            target = None
            for token in line.rstrip("\n").split("\t"):
                if ":" in token:
                    key, value = token.split(":", 1)
                    fields[key] = value
                else:
                    target = token
            if target is None:
                continue
            read_id = int(fields["ID"])
            gene, target_class = t2g.get(target, (base_gene(target), "."))
            entry = reads.setdefault(
                read_id,
                {
                    "targets": set(),
                    "genes": set(),
                    "classes": set(),
                    "dirs": set(),
                    "nh": int(fields.get("NH", "0") or 0),
                    "cb": fields.get("CB", ""),
                    "umi": fields.get("UMI", ""),
                },
            )
            entry["targets"].add(target)
            entry["genes"].add(gene)
            entry["classes"].add(target_class)
            entry["dirs"].add(fields.get("DIR", "."))
            entry["nh"] = max(int(entry["nh"]), int(fields.get("NH", "0") or 0))
    return reads


def class_label(classes: set[str]) -> str:
    classes = {c for c in classes if c and c != "."}
    if not classes:
        return "unknown"
    return "+".join(sorted(classes))


def print_counter(name: str, counter: collections.Counter, limit: int | None = None) -> None:
    print(f"\n{name}")
    items = counter.most_common(limit)
    for key, value in items:
        print(f"{key}\t{value}")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--assignments", required=True, type=Path)
    parser.add_argument("--constellation-features", type=Path)
    parser.add_argument("--unmapped-diagnostics", type=Path)
    parser.add_argument("--alevin-view", required=True, type=Path)
    parser.add_argument("--t2g", required=True, type=Path)
    args = parser.parse_args()

    features = read_constellation_features(args.constellation_features)
    assignments = read_assignments(args.assignments, features)
    reasons = read_unmapped_diagnostics(args.unmapped_diagnostics)
    alevin = parse_rad_view(args.alevin_view, read_t2g(args.t2g))

    all_read_ids = set(assignments)
    constellation_countable = {
        read_id
        for read_id, row in assignments.items()
        if row["assignment_type"] in COUNTABLE
    }
    alevin_mapped = set(alevin)

    print(f"total_constellation_reads_seen\t{len(all_read_ids)}")
    print(f"constellation_countable\t{len(constellation_countable)}")
    print(f"alevin_mapped\t{len(alevin_mapped)}")
    print("note\tRAD_IDs_are_mapped_record_ordinals_not_original_FASTQ_read_ids")

    print_counter(
        "constellation_assignment_by_read",
        collections.Counter(row["assignment_type"] for row in assignments.values()),
    )
    print_counter(
        "constellation_noncountable_diagnostic_reason_by_read",
        collections.Counter(reasons.values()),
    )
    print_counter(
        "alevin_gene_cardinality_by_mapped_record",
        collections.Counter(len(row["genes"]) for row in alevin.values()),
        limit=20,
    )
    print_counter(
        "alevin_target_class_by_mapped_record",
        collections.Counter(class_label(row["classes"]) for row in alevin.values()),
    )
    print_counter(
        "alevin_nh_by_mapped_record",
        collections.Counter(row["nh"] for row in alevin.values()),
        limit=20,
    )

    constellation_molecules: dict[tuple[str, str], collections.Counter] = collections.defaultdict(
        collections.Counter
    )
    constellation_molecule_genes: dict[tuple[str, str], collections.Counter] = (
        collections.defaultdict(collections.Counter)
    )
    for row in assignments.values():
        key = (row.get("cell_barcode", ""), row.get("umi", ""))
        # Older dictionaries from read_assignments do not expose CB/UMI.
        if not key[0]:
            continue
        constellation_molecules[key][row["assignment_type"]] += 1
        if row["assignment_type"] in COUNTABLE and row["gene"] != ".":
            constellation_molecule_genes[key][row["gene"]] += 1

    alevin_molecules: dict[tuple[str, str], collections.Counter] = collections.defaultdict(
        collections.Counter
    )
    alevin_molecule_classes: dict[tuple[str, str], collections.Counter] = (
        collections.defaultdict(collections.Counter)
    )
    for row in alevin.values():
        key = (str(row["cb"]), str(row["umi"]))
        for gene in row["genes"]:
            alevin_molecules[key][gene] += 1
        alevin_molecule_classes[key][class_label(row["classes"])] += 1

    if constellation_molecules:
        c_keys = set(constellation_molecules)
        c_countable_keys = set(constellation_molecule_genes)
        a_keys = set(alevin_molecules)
        print(f"\nconstellation_molecules_any_assignment\t{len(c_keys)}")
        print(f"constellation_molecules_countable\t{len(c_countable_keys)}")
        print(f"alevin_molecules_mapped\t{len(a_keys)}")
        print(f"molecules_both_countable_and_alevin_mapped\t{len(c_countable_keys & a_keys)}")
        print(f"molecules_alevin_mapped_not_constellation_countable\t{len(a_keys - c_countable_keys)}")
        print(f"molecules_constellation_countable_not_alevin_mapped\t{len(c_countable_keys - a_keys)}")

        print_counter(
            "constellation_assignment_for_alevin_mapped_not_countable_molecules",
            collections.Counter(
                "+".join(sorted(constellation_molecules.get(key, {})))
                if key in constellation_molecules
                else "absent_from_constellation"
                for key in a_keys - c_countable_keys
            ),
            limit=30,
        )
        print_counter(
            "alevin_target_class_for_alevin_mapped_not_countable_molecules",
            collections.Counter(
                alevin_molecule_classes[key].most_common(1)[0][0] for key in a_keys - c_countable_keys
            ),
        )
        top_gap_genes = collections.Counter()
        for key in a_keys - c_countable_keys:
            if len(alevin_molecules[key]) == 1:
                top_gap_genes.update(alevin_molecules[key])
        print_counter(
            "top_single_genes_for_alevin_mapped_not_countable_molecules",
            top_gap_genes,
            limit=20,
        )


if __name__ == "__main__":
    main()

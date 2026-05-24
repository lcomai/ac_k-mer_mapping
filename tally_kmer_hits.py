#!/usr/bin/env python3
"""Tally mapped genome positions per k-mer from a mapping TSV."""

from __future__ import annotations

import argparse
import csv
import sys
from collections import Counter
from pathlib import Path


KMER_COLUMN_CANDIDATES = ("kmer", "k-mer")


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Count the number of mapped genome positions per k-mer from a mapping TSV."
        )
    )
    parser.add_argument("mapping_tsv", type=Path, help="Mapping TSV from map_kmers_ac.py.")
    parser.add_argument("output_tsv", type=Path, help="Output count TSV.")
    parser.add_argument(
        "--kmers",
        type=Path,
        default=None,
        help="Optional original k-mer list. Includes valid unmapped k-mers with count 0.",
    )
    parser.add_argument(
        "--kmer-column",
        default=None,
        help="Name of the k-mer column. Default: auto-detect kmer or k-mer.",
    )
    return parser.parse_args(argv)


def read_kmer_list(kmers_file: Path) -> list[str]:
    kmers: list[str] = []
    seen: set[str] = set()
    with kmers_file.open() as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if not line or line.startswith("#"):
                continue
            kmer = line.split()[0].upper()
            if kmer not in seen:
                seen.add(kmer)
                kmers.append(kmer)
    return kmers


def detect_kmer_column(fieldnames: list[str] | None, requested: str | None) -> str:
    if not fieldnames:
        raise ValueError("Mapping TSV is missing a header row.")
    if requested:
        if requested not in fieldnames:
            raise ValueError(
                f"Requested k-mer column {requested!r} was not found. "
                f"Available columns: {', '.join(fieldnames)}"
            )
        return requested
    for candidate in KMER_COLUMN_CANDIDATES:
        if candidate in fieldnames:
            return candidate
    raise ValueError(
        "Could not find a k-mer column. Expected one of "
        f"{', '.join(KMER_COLUMN_CANDIDATES)}; available columns: {', '.join(fieldnames)}"
    )


def tally_mapping(mapping_tsv: Path, kmer_column: str | None = None) -> Counter:
    counts: Counter = Counter()
    with mapping_tsv.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        column = detect_kmer_column(reader.fieldnames, kmer_column)
        for row in reader:
            kmer = row.get(column, "").strip().upper()
            if kmer:
                counts[kmer] += 1
    return counts


def write_counts(output_tsv: Path, counts: Counter, ordered_kmers: list[str] | None) -> None:
    output_tsv.parent.mkdir(parents=True, exist_ok=True)
    if ordered_kmers is None:
        rows = sorted(counts.items(), key=lambda item: (-item[1], item[0]))
    else:
        rows = [(kmer, counts.get(kmer, 0)) for kmer in ordered_kmers]

    with output_tsv.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["kmer", "mapped_positions"])
        writer.writerows(rows)


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    try:
        counts = tally_mapping(args.mapping_tsv, args.kmer_column)
        ordered_kmers = read_kmer_list(args.kmers) if args.kmers else None
        write_counts(args.output_tsv, counts, ordered_kmers)
    except (FileNotFoundError, ValueError) as error:
        print(f"Error: {error}", file=sys.stderr)
        return 1

    total_kmers = len(ordered_kmers) if ordered_kmers is not None else len(counts)
    print(
        f"Wrote counts for {total_kmers} k-mers and "
        f"{sum(counts.values())} mapped positions to {args.output_tsv}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

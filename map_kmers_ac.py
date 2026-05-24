#!/usr/bin/env python3
"""Map exact k-mer instances in FASTA files with an Aho-Corasick automaton."""

from __future__ import annotations

import argparse
import csv
import sys
from collections import Counter, defaultdict
from pathlib import Path
from typing import Iterable, Iterator


FASTA_SUFFIXES = (".fa", ".fasta", ".fna", ".fas", ".fsa")
DNA_BASES = set("ACGT")
DNA_BASES_WITH_N = set("ACGTN")
COMPLEMENT = str.maketrans("ACGTN", "TGCAN")


def revcomp(seq: str) -> str:
    """Return the reverse complement of an uppercase DNA sequence."""
    return seq.translate(COMPLEMENT)[::-1]


def clean_tsv_field(value: str) -> str:
    """Keep FASTA headers usable as single TSV fields."""
    return value.replace("\t", " ").strip()


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Map exact k-mer occurrences in one FASTA file or all FASTA files in "
            "a directory. Positions are 1-based inclusive coordinates."
        )
    )
    parser.add_argument("kmers", type=Path, help="Text file with one k-mer per line.")
    parser.add_argument(
        "genome",
        type=Path,
        help="Reference FASTA file, or a directory containing FASTA files.",
    )
    parser.add_argument("output", type=Path, help="Output TSV path.")
    parser.add_argument(
        "--recursive",
        action="store_true",
        help="Search FASTA files recursively when genome is a directory.",
    )
    parser.add_argument(
        "--fasta-glob",
        default=None,
        help=(
            "When genome is a directory, scan FASTA files matching this glob, "
            "for example '*Col-PEK1.5*.fasta'."
        ),
    )
    parser.add_argument(
        "--all-fastas",
        action="store_true",
        help=(
            "When genome is a directory, allow scanning every discovered FASTA. "
            "Without this flag, directories with multiple FASTAs are rejected."
        ),
    )
    parser.add_argument(
        "--name-mode",
        choices=("id", "full"),
        default="id",
        help=(
            "How FASTA record names are written. 'id' uses the first header token; "
            "'full' preserves the full header line. Default: id."
        ),
    )
    parser.add_argument(
        "--allow-n",
        action="store_true",
        help="Allow N in input k-mers as a literal N, not as a wildcard.",
    )
    parser.add_argument(
        "--forward-only",
        action="store_true",
        help="Search only the input k-mer sequence, not its reverse complement.",
    )
    parser.add_argument(
        "--unmatched-output",
        type=Path,
        default=None,
        help="Path for unmatched k-mers. Default: <output>.no_match",
    )
    parser.add_argument(
        "--invalid-output",
        type=Path,
        default=None,
        help="Path for invalid or skipped input lines. Default: <output>.invalid.tsv",
    )
    parser.add_argument(
        "--summary-output",
        type=Path,
        default=None,
        help="Path for run summary TSV. Default: <output>.summary.tsv",
    )
    return parser.parse_args(argv)


def discover_fasta_paths(
    genome: Path,
    recursive: bool = False,
    fasta_glob: str | None = None,
    all_fastas: bool = False,
) -> list[Path]:
    if genome.is_file():
        return [genome]
    if not genome.is_dir():
        raise FileNotFoundError(f"Genome path does not exist: {genome}")

    if fasta_glob:
        iterator = genome.rglob(fasta_glob) if recursive else genome.glob(fasta_glob)
        fasta_paths = [path for path in iterator if path.is_file()]
    else:
        iterator = genome.rglob("*") if recursive else genome.iterdir()
        fasta_paths = [
            path
            for path in iterator
            if path.is_file() and path.suffix.lower() in FASTA_SUFFIXES
        ]
    if not fasta_paths:
        if fasta_glob:
            raise FileNotFoundError(
                f"No FASTA files matching {fasta_glob!r} found in {genome}"
            )
        suffixes = ", ".join(FASTA_SUFFIXES)
        raise FileNotFoundError(f"No FASTA files with suffixes {suffixes} found in {genome}")
    if len(fasta_paths) > 1 and not (all_fastas or fasta_glob):
        files = "\n  ".join(str(path) for path in sorted(fasta_paths))
        raise ValueError(
            "Genome directory contains multiple FASTA files. Use a specific FASTA "
            "path, pass --fasta-glob to select one set, or pass --all-fastas to "
            f"scan all of them:\n  {files}"
        )
    return sorted(fasta_paths)


def read_kmers(
    kmers_file: Path, allow_n: bool = False
) -> tuple[list[str], list[tuple[int, str, str]], Counter]:
    valid_bases = DNA_BASES_WITH_N if allow_n else DNA_BASES
    kmers: list[str] = []
    invalid: list[tuple[int, str, str]] = []
    counts: Counter = Counter()

    with kmers_file.open() as handle:
        for line_number, raw_line in enumerate(handle, start=1):
            line = raw_line.strip()
            if not line or line.startswith("#"):
                continue
            kmer = line.split()[0].upper()
            if not set(kmer) <= valid_bases:
                invalid.append((line_number, line, "non_acgt_base"))
                continue
            if len(kmer) == 0:
                invalid.append((line_number, line, "empty_kmer"))
                continue
            counts[kmer] += 1
            if counts[kmer] == 1:
                kmers.append(kmer)

    if not kmers:
        raise ValueError(f"No valid k-mers found in {kmers_file}")
    return kmers, invalid, counts


def build_automaton(
    kmers: Iterable[str], forward_only: bool = False
):
    try:
        import ahocorasick
    except ImportError as error:
        raise RuntimeError(
            "Missing dependency 'pyahocorasick'. Install with "
            "'python3 -m pip install -r requirements.txt'."
        ) from error

    pattern_hits: dict[str, list[tuple[str, str]]] = defaultdict(list)

    for kmer in kmers:
        pattern_hits[kmer].append((kmer, "+"))
        if not forward_only:
            rc = revcomp(kmer)
            if rc != kmer:
                pattern_hits[rc].append((kmer, "-"))

    automaton = ahocorasick.Automaton()
    for pattern, hits in pattern_hits.items():
        automaton.add_word(pattern, (pattern, hits))
    automaton.make_automaton()
    return automaton


def fasta_records(path: Path, name_mode: str = "id") -> Iterator[tuple[str, str]]:
    header: str | None = None
    seq_chunks: list[str] = []

    with path.open() as handle:
        for raw_line in handle:
            line = raw_line.rstrip("\n\r")
            if line.startswith(">"):
                if header is not None:
                    yield format_record_name(header, name_mode), "".join(seq_chunks)
                header = line[1:].strip()
                seq_chunks = []
            elif header is not None:
                seq_chunks.append(line.strip())

    if header is not None:
        yield format_record_name(header, name_mode), "".join(seq_chunks)


def format_record_name(header: str, name_mode: str) -> str:
    if name_mode == "full":
        return clean_tsv_field(header)
    return clean_tsv_field(header.split()[0])


def default_sidecar(output: Path, suffix: str) -> Path:
    return Path(str(output) + suffix)


def write_invalid(path: Path, invalid: list[tuple[int, str, str]]) -> None:
    with path.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["line", "value", "reason"])
        writer.writerows(invalid)


def write_unmatched(path: Path, kmers: Iterable[str]) -> None:
    with path.open("w") as handle:
        for kmer in sorted(kmers):
            handle.write(f"{kmer}\n")


def write_summary(path: Path, summary: dict[str, int | str]) -> None:
    with path.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["metric", "value"])
        for key, value in summary.items():
            writer.writerow([key, value])


def map_kmers(args: argparse.Namespace) -> dict[str, int | str]:
    fasta_paths = discover_fasta_paths(
        args.genome, args.recursive, args.fasta_glob, args.all_fastas
    )
    kmers, invalid, kmer_counts = read_kmers(args.kmers, args.allow_n)
    automaton = build_automaton(kmers, args.forward_only)

    args.output.parent.mkdir(parents=True, exist_ok=True)
    matched: set[str] = set()
    records_scanned = 0
    bases_scanned = 0
    hits_written = 0

    with args.output.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["kmer", "seq_id", "start", "end", "strand", "source_fasta"])

        for fasta_path in fasta_paths:
            for seq_id, sequence in fasta_records(fasta_path, args.name_mode):
                records_scanned += 1
                seq = sequence.upper()
                bases_scanned += len(seq)
                for end_zero, (_, hits) in automaton.iter(seq):
                    for original_kmer, strand in hits:
                        start = end_zero - len(original_kmer) + 2
                        end = end_zero + 1
                        writer.writerow(
                            [
                                original_kmer,
                                seq_id,
                                start,
                                end,
                                strand,
                                str(fasta_path),
                            ]
                        )
                        matched.add(original_kmer)
                        hits_written += 1

    unmatched = set(kmers) - matched
    unmatched_output = args.unmatched_output or default_sidecar(args.output, ".no_match")
    invalid_output = args.invalid_output or default_sidecar(args.output, ".invalid.tsv")
    summary_output = args.summary_output or default_sidecar(args.output, ".summary.tsv")

    write_unmatched(unmatched_output, unmatched)
    write_invalid(invalid_output, invalid)

    summary = {
        "input_kmer_lines": sum(kmer_counts.values()) + len(invalid),
        "unique_valid_kmers": len(kmers),
        "duplicate_valid_kmer_lines": sum(count - 1 for count in kmer_counts.values()),
        "invalid_or_skipped_lines": len(invalid),
        "fasta_files_scanned": len(fasta_paths),
        "records_scanned": records_scanned,
        "bases_scanned": bases_scanned,
        "hits_written": hits_written,
        "matched_kmers": len(matched),
        "unmatched_kmers": len(unmatched),
        "output_tsv": str(args.output),
        "unmatched_output": str(unmatched_output),
        "invalid_output": str(invalid_output),
    }
    write_summary(summary_output, summary)
    summary["summary_output"] = str(summary_output)
    return summary


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    try:
        summary = map_kmers(args)
    except (FileNotFoundError, RuntimeError, ValueError) as error:
        print(f"Error: {error}", file=sys.stderr)
        return 1
    print(
        "Mapped {matched_kmers}/{unique_valid_kmers} unique k-mers; "
        "wrote {hits_written} hits to {output_tsv}".format(**summary)
    )
    print(f"Summary: {summary['summary_output']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

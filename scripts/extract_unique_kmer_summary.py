#!/usr/bin/env python3
"""Extract BBNorm/BBNorm-rs unique-kmer estimates from stderr logs."""

from __future__ import annotations

import csv
import re
import sys
from pathlib import Path


FIELDS = [
    "tool",
    "log",
    "reads_in",
    "bases_in",
    "reads_kept",
    "reads_tossed",
    "estimated_unique_kmers",
    "low_depth_max",
    "low_depth_kmers",
    "high_depth_min",
    "high_depth_kmers",
    "hist_raw_kmers",
    "hist_unique_kmers",
    "rhist_reads",
    "rhist_bases",
    "sketch_tables",
    "sketch_total_cells",
    "sketch_memory_bytes",
    "sketch_layouts",
    "input_cardinality_k",
    "input_cardinality_buckets",
    "input_cardinality_unique_kmers",
    "output_cardinality_k",
    "output_cardinality_buckets",
    "output_cardinality_unique_kmers",
    "countup_spill_initial_runs",
    "countup_spill_merge_runs",
    "countup_spill_final_runs",
    "countup_spill_bytes_written",
    "countup_spill_peak_live_bytes",
    "countup_spill_final_live_bytes",
]


def parse_int(text: str) -> str:
    return text.replace(",", "")


def parse_scaled_number(text: str, unit: str = "") -> int:
    multipliers = {
        "": 1,
        "K": 1_000,
        "M": 1_000_000,
        "G": 1_000_000_000,
        "T": 1_000_000_000_000,
    }
    value = float(text.replace(",", ""))
    return int(round(value * multipliers.get(unit.upper(), 1)))


def parse_memory_bytes(text: str, unit: str) -> int:
    multipliers = {
        "B": 1,
        "KB": 1024,
        "MB": 1024**2,
        "GB": 1024**3,
        "TB": 1024**4,
    }
    value = float(text.replace(",", ""))
    return int(round(value * multipliers.get(unit.upper(), 1)))


def add_sketch_layout(row: dict[str, str], fields: dict[str, str]) -> None:
    cells = fields.get("cells", "")
    memory_bytes = fields.get("memory_bytes", "")
    if row["sketch_tables"]:
        row["sketch_tables"] = str(int(row["sketch_tables"]) + 1)
    else:
        row["sketch_tables"] = "1"
    if cells.isdigit():
        prior = int(row["sketch_total_cells"] or "0")
        row["sketch_total_cells"] = str(prior + int(cells))
    if memory_bytes.isdigit():
        prior = int(row["sketch_memory_bytes"] or "0")
        row["sketch_memory_bytes"] = str(prior + int(memory_bytes))
    layout = ";".join(f"{key}={value}" for key, value in fields.items() if value != "")
    row["sketch_layouts"] = (
        f"{row['sketch_layouts']}|{layout}" if row["sketch_layouts"] else layout
    )


def parse_rust_layout(text: str) -> dict[str, str]:
    fields: dict[str, str] = {}
    for part in text.split(";"):
        if "=" not in part:
            continue
        key, value = part.split("=", 1)
        fields[key.strip()] = value.strip()
    return fields


def add_hist_totals(row: dict[str, str], hist_path: Path | None) -> None:
    if hist_path is None:
        return
    raw_kmers = 0
    unique_kmers = 0
    try:
        lines = hist_path.read_text(encoding="utf-8", errors="replace").splitlines()
    except FileNotFoundError:
        return
    for line in lines:
        if not line or line.startswith("#"):
            continue
        fields = line.split("\t")
        if len(fields) < 3:
            continue
        try:
            raw_kmers += int(fields[1])
            unique_kmers += int(fields[2])
        except ValueError:
            continue
    row["hist_raw_kmers"] = str(raw_kmers)
    row["hist_unique_kmers"] = str(unique_kmers)


def add_rhist_totals(row: dict[str, str], rhist_path: Path | None) -> None:
    if rhist_path is None:
        return
    reads = 0
    bases = 0
    try:
        lines = rhist_path.read_text(encoding="utf-8", errors="replace").splitlines()
    except FileNotFoundError:
        return
    for line in lines:
        if not line or line.startswith("#"):
            continue
        fields = line.split("\t")
        if len(fields) < 3:
            continue
        try:
            reads += int(fields[1])
            bases += int(fields[2])
        except ValueError:
            continue
    row["rhist_reads"] = str(reads)
    row["rhist_bases"] = str(bases)


def parse_log(
    tool: str,
    log_path: Path,
    hist_path: Path | None = None,
    rhist_path: Path | None = None,
) -> dict[str, str]:
    row = {field: "" for field in FIELDS}
    row["tool"] = tool
    row["log"] = str(log_path)
    try:
        lines = log_path.read_text(encoding="utf-8", errors="replace").splitlines()
    except FileNotFoundError:
        return row

    rust_summary = re.compile(
        r"Processed\s+([\d,]+)\s+reads\s+\(([\d,]+)\s+bases\);\s+"
        r"kept\s+([\d,]+)\s+reads,\s+tossed\s+([\d,]+)\s+reads;\s+"
        r"input unique kmers:\s+([\d,]+)"
        r"(?:;\s+input unique kmers depth 1-([\d,]+):\s+([\d,]+);\s+"
        r"input unique kmers depth\s+([\d,]+)\+:\s+([\d,]+))?"
    )
    java_low = re.compile(r"Estimated kmers of depth 1-([\d,]+):\s*([\d,]+)")
    java_high = re.compile(r"Estimated kmers of depth\s+([\d,]+)\+\s*:?\s*([\d,]+)")
    java_total = re.compile(r"Estimated unique kmers:\s*([\d,]+)")
    java_reads = re.compile(r"Total reads in:\s*([\d,]+)")
    java_bases = re.compile(r"Total bases in:\s*([\d,]+)")
    rust_layout = re.compile(r"Sketch layout:\s*(.*)")
    rust_countup_spill = re.compile(r"Count-up spill:\s*(.*)")
    rust_cardinality = re.compile(r"Cardinality estimate:\s*(.*)")
    java_bits = re.compile(r"bits per cell:\s*([\d,]+)")
    java_prefilter_bits = re.compile(r"prefilter bits:\s*([\d,]+)")
    java_made_table = re.compile(
        r"Made\s+(prefilter|hash table):\s*.*?hashes\s*=\s*([\d,]+)"
        r".*?mem\s*=\s*([\d.]+)\s*([KMGT]?B)"
        r".*?cells\s*=\s*([\d.]+)\s*([KMGT]?)",
        re.IGNORECASE,
    )
    main_bits = ""
    prefilter_bits = ""

    for line in lines:
        if match := rust_summary.search(line):
            row["reads_in"] = parse_int(match.group(1))
            row["bases_in"] = parse_int(match.group(2))
            row["reads_kept"] = parse_int(match.group(3))
            row["reads_tossed"] = parse_int(match.group(4))
            row["estimated_unique_kmers"] = parse_int(match.group(5))
            if match.group(6):
                row["low_depth_max"] = parse_int(match.group(6))
                row["low_depth_kmers"] = parse_int(match.group(7))
                row["high_depth_min"] = parse_int(match.group(8))
                row["high_depth_kmers"] = parse_int(match.group(9))
            continue
        if match := java_low.search(line):
            row["low_depth_max"] = parse_int(match.group(1))
            row["low_depth_kmers"] = parse_int(match.group(2))
            continue
        if match := java_high.search(line):
            row["high_depth_min"] = parse_int(match.group(1))
            row["high_depth_kmers"] = parse_int(match.group(2))
            continue
        if match := java_total.search(line):
            row["estimated_unique_kmers"] = parse_int(match.group(1))
            continue
        if match := java_reads.search(line):
            row["reads_in"] = parse_int(match.group(1))
            continue
        if match := java_bases.search(line):
            row["bases_in"] = parse_int(match.group(1))
            continue
        if match := rust_layout.search(line):
            add_sketch_layout(row, parse_rust_layout(match.group(1)))
            continue
        if match := rust_countup_spill.search(line):
            fields = parse_rust_layout(match.group(1))
            row["countup_spill_initial_runs"] = fields.get("initial_runs", "")
            row["countup_spill_merge_runs"] = fields.get("merge_runs", "")
            row["countup_spill_final_runs"] = fields.get("final_runs", "")
            row["countup_spill_bytes_written"] = fields.get("bytes_written", "")
            row["countup_spill_peak_live_bytes"] = fields.get("peak_live_bytes", "")
            row["countup_spill_final_live_bytes"] = fields.get("final_live_bytes", "")
            continue
        if match := rust_cardinality.search(line):
            fields = parse_rust_layout(match.group(1))
            scope = fields.get("scope", "")
            if scope in {"input", "output"}:
                prefix = f"{scope}_cardinality"
                row[f"{prefix}_k"] = fields.get("k", "")
                row[f"{prefix}_buckets"] = fields.get("buckets", "")
                row[f"{prefix}_unique_kmers"] = fields.get("unique_kmers", "")
            continue
        if match := java_bits.search(line):
            main_bits = parse_int(match.group(1))
            continue
        if match := java_prefilter_bits.search(line):
            prefilter_bits = parse_int(match.group(1))
            continue
        if match := java_made_table.search(line):
            label = match.group(1).lower()
            is_prefilter = label == "prefilter"
            add_sketch_layout(
                row,
                {
                    "table": "input_prefilter" if is_prefilter else "input_main",
                    "kind": "java",
                    "cells": str(parse_scaled_number(match.group(5), match.group(6))),
                    "hashes": parse_int(match.group(2)),
                    "bits": prefilter_bits if is_prefilter else main_bits,
                    "memory_bytes": str(parse_memory_bytes(match.group(3), match.group(4))),
                },
            )
            continue

    add_hist_totals(row, hist_path)
    add_rhist_totals(row, rhist_path)
    return row


def parse_arg(arg: str) -> tuple[str, Path, Path | None, Path | None]:
    if "=" in arg:
        tool, paths_text = arg.split("=", 1)
        paths = paths_text.split(",")
        return (
            tool,
            Path(paths[0]),
            Path(paths[1]) if len(paths) > 1 and paths[1] else None,
            Path(paths[2]) if len(paths) > 2 and paths[2] else None,
        )
    paths = arg.split(",")
    path = Path(paths[0])
    return (
        path.stem,
        path,
        Path(paths[1]) if len(paths) > 1 and paths[1] else None,
        Path(paths[2]) if len(paths) > 2 and paths[2] else None,
    )


def main() -> int:
    if len(sys.argv) < 2:
        print(
            "Usage: scripts/extract_unique_kmer_summary.py tool=stderr.log[,hist.tsv[,rhist.tsv]] [...]",
            file=sys.stderr,
        )
        return 2

    writer = csv.DictWriter(
        sys.stdout, fieldnames=FIELDS, delimiter="\t", lineterminator="\n"
    )
    writer.writeheader()
    for arg in sys.argv[1:]:
        tool, log_path, hist_path, rhist_path = parse_arg(arg)
        writer.writerow(parse_log(tool, log_path, hist_path, rhist_path))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

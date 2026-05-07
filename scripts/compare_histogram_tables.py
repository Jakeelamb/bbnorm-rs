#!/usr/bin/env python3
"""Summarize numeric deltas between Java and Rust histogram TSV files."""

from __future__ import annotations

import csv
import sys
from pathlib import Path


FIELDS = [
    "label",
    "status",
    "left_path",
    "right_path",
    "rows_left",
    "rows_right",
    "first_diff_key",
    "col2_left",
    "col2_right",
    "col2_delta",
    "col2_abs_delta_sum",
    "col2_abs_delta_ppm",
    "col2_max_abs_delta",
    "col3_left",
    "col3_right",
    "col3_delta",
    "col3_abs_delta_sum",
    "col3_abs_delta_ppm",
    "col3_max_abs_delta",
]


def read_rows(path: Path) -> dict[str, tuple[int, int]]:
    rows: dict[str, tuple[int, int]] = {}
    for line in path.read_text(encoding="utf-8", errors="replace").splitlines():
        if not line or line.startswith("#"):
            continue
        fields = line.split("\t")
        if len(fields) < 3:
            continue
        try:
            rows[fields[0]] = (int(fields[1]), int(fields[2]))
        except ValueError:
            continue
    return rows


def delta_ppm(abs_delta: int, left_total: int, right_total: int) -> str:
    denominator = left_total if left_total > 0 else right_total
    if denominator <= 0:
        return "0" if abs_delta == 0 else ""
    return str((abs_delta * 1_000_000 + denominator - 1) // denominator)


def compare(label: str, left_path: Path, right_path: Path) -> dict[str, str]:
    row = {field: "" for field in FIELDS}
    row["label"] = label
    row["left_path"] = str(left_path)
    row["right_path"] = str(right_path)
    if not left_path.exists() or not right_path.exists():
        row["status"] = "missing"
        return row

    left = read_rows(left_path)
    right = read_rows(right_path)
    keys = sorted(set(left) | set(right), key=lambda value: int(value) if value.isdigit() else value)
    row["rows_left"] = str(len(left))
    row["rows_right"] = str(len(right))

    left2 = right2 = left3 = right3 = 0
    abs_delta2 = abs_delta3 = 0
    max_delta2 = max_delta3 = 0
    first_diff_key = ""
    for key in keys:
        l2, l3 = left.get(key, (0, 0))
        r2, r3 = right.get(key, (0, 0))
        if not first_diff_key and (l2, l3) != (r2, r3):
            first_diff_key = key
        left2 += l2
        right2 += r2
        left3 += l3
        right3 += r3
        d2 = abs(r2 - l2)
        d3 = abs(r3 - l3)
        abs_delta2 += d2
        abs_delta3 += d3
        max_delta2 = max(max_delta2, d2)
        max_delta3 = max(max_delta3, d3)

    row["status"] = "identical" if not first_diff_key else "different"
    row["first_diff_key"] = first_diff_key
    row["col2_left"] = str(left2)
    row["col2_right"] = str(right2)
    row["col2_delta"] = str(right2 - left2)
    row["col2_abs_delta_sum"] = str(abs_delta2)
    row["col2_abs_delta_ppm"] = delta_ppm(abs_delta2, left2, right2)
    row["col2_max_abs_delta"] = str(max_delta2)
    row["col3_left"] = str(left3)
    row["col3_right"] = str(right3)
    row["col3_delta"] = str(right3 - left3)
    row["col3_abs_delta_sum"] = str(abs_delta3)
    row["col3_abs_delta_ppm"] = delta_ppm(abs_delta3, left3, right3)
    row["col3_max_abs_delta"] = str(max_delta3)
    return row


def parse_arg(arg: str) -> tuple[str, Path, Path]:
    if "=" in arg:
        label, paths_text = arg.split("=", 1)
    else:
        label, paths_text = "hist", arg
    paths = paths_text.split(",", 1)
    if len(paths) != 2:
        raise ValueError(f"expected label=left.tsv,right.tsv, got: {arg}")
    return label, Path(paths[0]), Path(paths[1])


def main() -> int:
    if len(sys.argv) < 2:
        print(
            "Usage: scripts/compare_histogram_tables.py label=left.tsv,right.tsv [...]",
            file=sys.stderr,
        )
        return 2
    writer = csv.DictWriter(
        sys.stdout, fieldnames=FIELDS, delimiter="\t", lineterminator="\n"
    )
    writer.writeheader()
    try:
        for arg in sys.argv[1:]:
            writer.writerow(compare(*parse_arg(arg)))
    except ValueError as err:
        print(f"error: {err}", file=sys.stderr)
        return 2
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

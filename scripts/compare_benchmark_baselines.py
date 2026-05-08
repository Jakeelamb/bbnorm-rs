#!/usr/bin/env python3
"""Compare aggregate benchmark summaries from two baseline artifact dirs."""

from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path


DEFAULT_METRICS = [
    "elapsed_seconds",
    "max_rss_kb",
    "stage_input_counting",
    "stage_normalize",
    "stage_input_hist",
    "stage_table_creation",
    "stage_table_read",
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("baseline", type=Path)
    parser.add_argument("candidate", type=Path)
    parser.add_argument(
        "--metrics",
        default=",".join(DEFAULT_METRICS),
        help="Comma-separated metric names to compare",
    )
    return parser.parse_args()


def read_summary(path: Path) -> dict[tuple[str, str], dict[str, str]]:
    summary = path / "summary.tsv"
    if not summary.exists():
        raise SystemExit(f"missing summary.tsv: {summary}")
    with summary.open("r", encoding="utf-8", newline="") as handle:
        return {
            (row["variant"], row["metric"]): row
            for row in csv.DictReader(handle, delimiter="\t")
        }


def parse_float(value: str) -> float | None:
    if value == "":
        return None
    try:
        return float(value)
    except ValueError:
        return None


def format_float(value: float | None) -> str:
    if value is None:
        return ""
    return f"{value:.6f}"


def main() -> int:
    args = parse_args()
    metrics = [part.strip() for part in args.metrics.split(",") if part.strip()]
    baseline = read_summary(args.baseline)
    candidate = read_summary(args.candidate)
    keys = sorted(set(baseline) | set(candidate))

    writer = csv.DictWriter(
        sys.stdout,
        fieldnames=[
            "variant",
            "metric",
            "baseline_median",
            "candidate_median",
            "delta",
            "ratio_candidate_baseline",
            "baseline_n",
            "candidate_n",
        ],
        delimiter="\t",
        lineterminator="\n",
    )
    writer.writeheader()
    for variant, metric in keys:
        if metric not in metrics:
            continue
        base_row = baseline.get((variant, metric), {})
        cand_row = candidate.get((variant, metric), {})
        base = parse_float(base_row.get("median", ""))
        cand = parse_float(cand_row.get("median", ""))
        delta = cand - base if base is not None and cand is not None else None
        ratio = cand / base if base not in {None, 0.0} and cand is not None else None
        writer.writerow(
            {
                "variant": variant,
                "metric": metric,
                "baseline_median": format_float(base),
                "candidate_median": format_float(cand),
                "delta": format_float(delta),
                "ratio_candidate_baseline": format_float(ratio),
                "baseline_n": base_row.get("n", ""),
                "candidate_n": cand_row.get("n", ""),
            }
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

#!/usr/bin/env python3
"""Emit a three-column Java CPU / Rust CPU / Rust GPU comparison table.

The GPU column is intentionally labeled as a counting probe until GPU-assisted
counts are integrated into the full bbnorm-rs normalization pipeline.
"""

from __future__ import annotations

import argparse
import csv
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_CPU_SUMMARY = (
    ROOT / "tmp/parity_acceptance_publish_ready_20260508/human_bounded_core_500k/summary.tsv"
)
DEFAULT_GPU_REPORT = ROOT / "tmp/cuda_full_read_dataset_20260508_173427/report.tsv"
DEFAULT_MODE = "default"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cpu-summary", type=Path, default=DEFAULT_CPU_SUMMARY)
    parser.add_argument("--gpu-report", type=Path, default=DEFAULT_GPU_REPORT)
    parser.add_argument("--mode", default=DEFAULT_MODE)
    parser.add_argument(
        "--out",
        type=Path,
        default=ROOT / "tmp/cpu_gpu_three_column_comparison.tsv",
    )
    return parser.parse_args()


def read_cpu_row(path: Path, mode: str) -> dict[str, str]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            if row.get("mode") == mode:
                return row
    raise SystemExit(f"mode {mode!r} not found in {path}")


def read_report(path: Path) -> dict[str, str]:
    values = {}
    with path.open("r", encoding="utf-8") as handle:
        for line in handle:
            line = line.strip()
            if not line or "\t" not in line:
                continue
            key, value = line.split("\t", 1)
            values[key] = value
    return values


def fmt_seconds(value: str) -> str:
    if value == "":
        return ""
    return f"{float(value):.3f} s"


def fmt_ms(value: str) -> str:
    if value == "":
        return ""
    return f"{float(value):.3f} ms"


def fmt_gib_from_kb(value: str) -> str:
    if value == "":
        return ""
    return f"{int(value) / 1024 / 1024:.2f} GiB"


def fmt_int(value: str) -> str:
    if value == "":
        return ""
    return f"{int(value):,}"


def write_comparison(cpu: dict[str, str], gpu: dict[str, str], out: Path) -> None:
    read_pairs = fmt_int(gpu["reads"])
    total_reads = fmt_int(str(int(gpu["reads"]) * 2))
    rows = [
        ("scope", "full tool", "full tool", "counting probe only"),
        ("read_pairs", read_pairs, read_pairs, read_pairs),
        ("reads_processed", total_reads, total_reads, total_reads),
        (
            "kmers_processed",
            fmt_int(gpu["kmers"]),
            fmt_int(gpu["kmers"]),
            fmt_int(gpu["kmers"]),
        ),
        (
            "end_to_end_wall",
            fmt_seconds(cpu["java_seconds"]),
            fmt_seconds(cpu["rust_seconds"]),
            fmt_seconds(gpu["pipeline_wall_seconds"]),
        ),
        (
            "counting_or_table_phase",
            fmt_seconds(cpu["java_table_creation_s"]),
            fmt_seconds(cpu["rust_input_counting_s"]),
            fmt_ms(gpu["gpu_total_timed_ms"]),
        ),
        (
            "peak_rss",
            fmt_gib_from_kb(cpu["java_rss_kb"]),
            fmt_gib_from_kb(cpu["rust_rss_kb"]),
            "not equivalent",
        ),
        (
            "unique_kmers",
            fmt_int(cpu["java_unique_kmers"]),
            fmt_int(cpu["rust_unique_kmers"]),
            fmt_int(gpu["unique"]),
        ),
        (
            "validation",
            cpu["drift_gate"],
            cpu["drift_gate"],
            f"checksum {gpu['gpu_checksum']}; {gpu['match']}",
        ),
        ("integrated_normalization", "yes", "yes", "no"),
    ]

    out.parent.mkdir(parents=True, exist_ok=True)
    with out.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["metric", "java_bbnorm_cpu", "rust_bbnorm_cpu", "rust_bbnorm_gpu"])
        writer.writerows(rows)


def main() -> int:
    args = parse_args()
    cpu = read_cpu_row(args.cpu_summary, args.mode)
    gpu = read_report(args.gpu_report)
    write_comparison(cpu, gpu, args.out)
    print(args.out.read_text(encoding="utf-8"), end="")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

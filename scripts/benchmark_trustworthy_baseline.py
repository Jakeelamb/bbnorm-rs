#!/usr/bin/env python3
"""Run repeated Java/Rust BBNorm baseline benchmarks with metadata.

This harness is for publishable benchmark evidence, not quick smoke testing.
It records host/tool/input metadata, runs Java plus deterministic and
nondeterministic Rust variants, and writes raw per-run TSVs plus aggregate
median/p10/p90 summaries.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
import platform
import shlex
import socket
import subprocess
import sys
import time
from datetime import datetime, timezone
from pathlib import Path
from statistics import mean


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_R1 = ROOT / "tmp/human_benchmark_8threads/human_GRCh38_500k_R1.fq.gz"
DEFAULT_R2 = ROOT / "tmp/human_benchmark_8threads/human_GRCh38_500k_R2.fq.gz"
MEASURE_SCRIPT = ROOT / "scripts/measure_command.py"

RUST_STAGE_PREFIX = "Stage timing: name="
JAVA_STAGE_PATTERNS = {
    "table_creation": "Table creation time:",
    "table_read": "Table read time:",
    "total": "Total time:",
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--r1", type=Path, default=DEFAULT_R1)
    parser.add_argument("--r2", type=Path, default=DEFAULT_R2)
    parser.add_argument(
        "--outdir",
        type=Path,
        default=ROOT
        / "tmp"
        / f"trustworthy_baseline_{datetime.now(timezone.utc).strftime('%Y%m%d_%H%M%S')}",
    )
    parser.add_argument("--repeats", type=int, default=3)
    parser.add_argument("--reads", default="100000")
    parser.add_argument("--table-reads", default=None)
    parser.add_argument("--threads", default="8")
    parser.add_argument("--zipthreads", default="1")
    parser.add_argument("--java-xmx", default="4g")
    parser.add_argument("--mem", default="4279m")
    parser.add_argument("--k", default="31")
    parser.add_argument("--min", default="5")
    parser.add_argument("--target", default="40")
    parser.add_argument("--max", default="80")
    parser.add_argument("--bits", default="")
    parser.add_argument(
        "--autocountmin-reads",
        default="1",
        help="Rust autocountminreads trigger; use empty/off/none/0 to omit",
    )
    parser.add_argument("--timeout", default="8m")
    parser.add_argument("--sample-interval", type=float, default=0.05)
    parser.add_argument("--extra-args", default="")
    parser.add_argument("--java-extra-args", default="")
    parser.add_argument("--rust-extra-args", default="")
    parser.add_argument(
        "--variants",
        default="java,rust_deterministic,rust_nondeterministic",
        help="Comma-separated subset of java,rust_deterministic,rust_nondeterministic",
    )
    parser.add_argument("--write-outputs", action="store_true")
    parser.add_argument("--skip-build", action="store_true")
    parser.add_argument("--skip-input-sha256", action="store_true")
    parser.add_argument(
        "--report-only",
        action="store_true",
        help="Regenerate report.md from an existing outdir without running benchmarks",
    )
    parser.add_argument("--classpath", default=str(ROOT / "vendor/BBTools-master/current"))
    args = parser.parse_args()

    if args.repeats < 1:
        parser.error("--repeats must be at least 1")
    if args.report_only:
        return args
    if not args.r1.exists():
        parser.error(f"missing R1 input: {args.r1}")
    if args.r2 and not args.r2.exists():
        parser.error(f"missing R2 input: {args.r2}")
    selected = [part.strip() for part in args.variants.split(",") if part.strip()]
    allowed = {"java", "rust_deterministic", "rust_nondeterministic"}
    unknown = sorted(set(selected) - allowed)
    if unknown:
        parser.error(f"unknown variants: {', '.join(unknown)}")
    args.variants = selected
    args.table_reads = args.table_reads or args.reads
    return args


def run_text(cmd: list[str]) -> str:
    try:
        completed = subprocess.run(
            cmd,
            cwd=ROOT,
            text=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            check=False,
        )
    except FileNotFoundError:
        return ""
    return completed.stdout.strip()


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def file_metadata(path: Path, include_sha256: bool) -> dict[str, object]:
    stat = path.stat()
    data: dict[str, object] = {
        "path": str(path),
        "bytes": stat.st_size,
        "mtime": stat.st_mtime,
    }
    if include_sha256:
        data["sha256"] = sha256_file(path)
    return data


def split_args(value: str) -> list[str]:
    return shlex.split(value) if value.strip() else []


def limit_arg(name: str, value: str) -> list[str]:
    lowered = value.strip().lower()
    if lowered in {"", "all", "none", "0"}:
        return []
    return [f"{name}={value}"]


def sequence_output_args(run_dir: Path, variant: str, paired: bool, write_outputs: bool) -> list[str]:
    if not write_outputs:
        args = ["out=null", "outt=null"]
        if paired:
            args.extend(["out2=null", "outt2=null"])
        return args
    args = [
        f"out={run_dir / f'{variant}.keep1.fq.gz'}",
        f"outt={run_dir / f'{variant}.toss1.fq.gz'}",
    ]
    if paired:
        args.extend(
            [
                f"out2={run_dir / f'{variant}.keep2.fq.gz'}",
                f"outt2={run_dir / f'{variant}.toss2.fq.gz'}",
            ]
        )
    return args


def common_bbnorm_args(args: argparse.Namespace) -> list[str]:
    common = [
        f"in={args.r1}",
        "passes=1",
        f"target={args.target}",
        f"max={args.max}",
        f"min={args.min}",
        f"k={args.k}",
        f"threads={args.threads}",
        f"zipthreads={args.zipthreads}",
        "overwrite=t",
    ]
    if args.r2:
        common.append(f"in2={args.r2}")
    if args.bits:
        common.append(f"bits={args.bits}")
    common.extend(limit_arg("reads", args.reads))
    common.extend(limit_arg("tablereads", args.table_reads))
    common.extend(split_args(args.extra_args))
    return common


def variant_command(args: argparse.Namespace, variant: str, run_dir: Path) -> list[str]:
    common = common_bbnorm_args(args)
    outputs = [
        f"hist={run_dir / f'{variant}.hist.tsv'}",
        f"rhist={run_dir / f'{variant}.rhist.tsv'}",
    ]
    outputs.extend(sequence_output_args(run_dir, variant, bool(args.r2), args.write_outputs))
    if variant == "java":
        return [
            "java",
            f"-Xmx{args.java_xmx}",
            "-cp",
            args.classpath,
            "jgi.KmerNormalize",
            *common,
            *split_args(args.java_extra_args),
            *outputs,
        ]
    rust_args = [
        "target/release/bbnorm-rs",
        *common,
        f"mem={args.mem}",
        "autocountmin=t",
        *split_args(args.rust_extra_args),
    ]
    if args.autocountmin_reads.strip().lower() not in {"", "off", "none", "0"}:
        rust_args.append(f"autocountminreads={args.autocountmin_reads}")
    if variant == "rust_nondeterministic":
        rust_args.append("deterministic=f")
    return [*rust_args, *outputs]


def run_measured(
    cmd: list[str],
    run_dir: Path,
    variant: str,
    timeout: str,
    sample_interval: float,
) -> dict[str, str]:
    metrics = run_dir / f"{variant}.metrics.tsv"
    stdout = run_dir / f"{variant}.stdout.log"
    stderr = run_dir / f"{variant}.stderr.log"
    measure_cmd = [
        sys.executable,
        str(MEASURE_SCRIPT),
        "--metrics",
        str(metrics),
        "--stdout",
        str(stdout),
        "--stderr",
        str(stderr),
        "--sample-interval",
        str(sample_interval),
    ]
    if timeout:
        measure_cmd.extend(["--timeout", timeout])
    completed = subprocess.run([*measure_cmd, "--", *cmd], cwd=ROOT, check=False)
    if not metrics.exists():
        return {
            "elapsed_seconds": "",
            "max_rss_kb": "",
            "status": str(completed.returncode),
            "stdout": str(stdout),
            "stderr": str(stderr),
        }
    fields = metrics.read_text(encoding="utf-8").strip().split("\t")
    while len(fields) < 3:
        fields.append("")
    return {
        "elapsed_seconds": fields[0],
        "max_rss_kb": fields[1],
        "status": fields[2],
        "stdout": str(stdout),
        "stderr": str(stderr),
    }


def parse_stage_timings(stderr_path: Path, variant: str) -> dict[str, float]:
    if not stderr_path.exists():
        return {}
    text = stderr_path.read_text(encoding="utf-8", errors="replace")
    stages: dict[str, float] = {}
    if variant.startswith("rust"):
        for line in text.splitlines():
            if RUST_STAGE_PREFIX not in line:
                continue
            after = line.split(RUST_STAGE_PREFIX, 1)[1]
            if "; seconds=" not in after:
                continue
            name, seconds = after.split("; seconds=", 1)
            try:
                stages[name.strip()] = float(seconds.strip())
            except ValueError:
                continue
    else:
        for stage, marker in JAVA_STAGE_PATTERNS.items():
            for line in text.splitlines():
                if marker not in line:
                    continue
                tail = line.split(marker, 1)[1].strip().split()
                if not tail:
                    continue
                try:
                    stages[stage] = float(tail[0])
                except ValueError:
                    pass
                break
    for line in text.splitlines():
        if line.startswith("Elapsed (wall clock) time (seconds):"):
            try:
                stages["wall_clock"] = float(line.rsplit(":", 1)[1].strip())
            except ValueError:
                pass
    return stages


def read_hist(path: Path) -> dict[str, tuple[int, int]]:
    rows: dict[str, tuple[int, int]] = {}
    if not path.exists():
        return rows
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


def compare_hist(label: str, java_path: Path, rust_path: Path) -> dict[str, str]:
    row = {
        "label": label,
        "status": "missing",
        "col2_abs_delta_sum": "",
        "col2_abs_delta_ppm": "",
        "col3_abs_delta_sum": "",
        "col3_abs_delta_ppm": "",
        "first_diff_key": "",
    }
    if not java_path.exists() or not rust_path.exists():
        return row
    java = read_hist(java_path)
    rust = read_hist(rust_path)
    keys = sorted(set(java) | set(rust), key=lambda value: int(value) if value.isdigit() else value)
    left2 = right2 = left3 = right3 = 0
    abs_delta2 = abs_delta3 = 0
    first = ""
    for key in keys:
        l2, l3 = java.get(key, (0, 0))
        r2, r3 = rust.get(key, (0, 0))
        if not first and (l2, l3) != (r2, r3):
            first = key
        left2 += l2
        right2 += r2
        left3 += l3
        right3 += r3
        abs_delta2 += abs(r2 - l2)
        abs_delta3 += abs(r3 - l3)
    row.update(
        {
            "status": "identical" if not first else "different",
            "col2_abs_delta_sum": str(abs_delta2),
            "col2_abs_delta_ppm": delta_ppm(abs_delta2, left2, right2),
            "col3_abs_delta_sum": str(abs_delta3),
            "col3_abs_delta_ppm": delta_ppm(abs_delta3, left3, right3),
            "first_diff_key": first,
        }
    )
    return row


def quantile(values: list[float], q: float) -> float:
    ordered = sorted(values)
    if not ordered:
        return float("nan")
    if len(ordered) == 1:
        return ordered[0]
    position = (len(ordered) - 1) * q
    lower = int(position)
    upper = min(lower + 1, len(ordered) - 1)
    fraction = position - lower
    return ordered[lower] * (1.0 - fraction) + ordered[upper] * fraction


def fmt_float(value: float) -> str:
    if value != value:
        return ""
    return f"{value:.6f}"


def write_tsv(path: Path, rows: list[dict[str, object]], fieldnames: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=fieldnames,
            delimiter="\t",
            lineterminator="\n",
            extrasaction="ignore",
        )
        writer.writeheader()
        writer.writerows(rows)


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def collect_metadata(args: argparse.Namespace) -> dict[str, object]:
    include_sha = not args.skip_input_sha256
    return {
        "timestamp_utc": datetime.now(timezone.utc).isoformat(),
        "host": socket.gethostname(),
        "platform": platform.platform(),
        "python": sys.version,
        "cwd": str(ROOT),
        "git_commit": run_text(["git", "rev-parse", "HEAD"]),
        "git_branch": run_text(["git", "branch", "--show-current"]),
        "git_status_short": run_text(["git", "status", "--short", "--branch"]),
        "rustc": run_text(["rustc", "--version"]),
        "cargo": run_text(["cargo", "--version"]),
        "java": run_text(["java", "-version"]),
        "inputs": {
            "r1": file_metadata(args.r1, include_sha),
            "r2": file_metadata(args.r2, include_sha) if args.r2 else None,
        },
        "settings": {
            "repeats": args.repeats,
            "reads": args.reads,
            "table_reads": args.table_reads,
            "threads": args.threads,
            "zipthreads": args.zipthreads,
            "java_xmx": args.java_xmx,
            "mem": args.mem,
            "k": args.k,
            "min": args.min,
            "target": args.target,
            "max": args.max,
            "bits": args.bits,
            "autocountmin_reads": args.autocountmin_reads,
            "timeout": args.timeout,
            "extra_args": args.extra_args,
            "java_extra_args": args.java_extra_args,
            "rust_extra_args": args.rust_extra_args,
            "variants": args.variants,
            "write_outputs": args.write_outputs,
        },
    }


def aggregate(raw_rows: list[dict[str, object]]) -> list[dict[str, str]]:
    numeric_metrics = [
        "elapsed_seconds",
        "max_rss_kb",
        "stage_input_counting",
        "stage_input_main_counting",
        "stage_input_exact_counting",
        "stage_normalize",
        "stage_input_hist",
        "stage_table_creation",
        "stage_table_read",
        "stage_total",
        "stage_wall_clock",
    ]
    rows: list[dict[str, str]] = []
    variants = sorted({str(row["variant"]) for row in raw_rows})
    for variant in variants:
        subset = [row for row in raw_rows if row["variant"] == variant and str(row["status"]) == "0"]
        for metric in numeric_metrics:
            values: list[float] = []
            for row in subset:
                value = row.get(metric, "")
                if value in {"", None}:
                    continue
                try:
                    values.append(float(value))
                except (TypeError, ValueError):
                    continue
            if not values:
                continue
            rows.append(
                {
                    "variant": variant,
                    "metric": metric,
                    "n": str(len(values)),
                    "min": fmt_float(min(values)),
                    "p10": fmt_float(quantile(values, 0.10)),
                    "median": fmt_float(quantile(values, 0.50)),
                    "p90": fmt_float(quantile(values, 0.90)),
                    "max": fmt_float(max(values)),
                    "mean": fmt_float(mean(values)),
                }
            )
    return rows


def summary_lookup(summary_rows: list[dict[str, str]]) -> dict[tuple[str, str], dict[str, str]]:
    return {(row["variant"], row["metric"]): row for row in summary_rows}


def median_value(
    lookup: dict[tuple[str, str], dict[str, str]], variant: str, metric: str
) -> float | None:
    value = lookup.get((variant, metric), {}).get("median", "")
    if value == "":
        return None
    try:
        return float(value)
    except ValueError:
        return None


def ratio_text(numerator: float | None, denominator: float | None) -> str:
    if numerator is None or denominator is None or denominator == 0.0:
        return ""
    return f"{numerator / denominator:.3f}"


def comparison_medians(comparison_rows: list[dict[str, object]]) -> list[dict[str, str]]:
    grouped: dict[tuple[str, str], dict[str, list[float]]] = {}
    for row in comparison_rows:
        key = (str(row.get("rust_variant", "")), str(row.get("label", "")))
        bucket = grouped.setdefault(
            key,
            {
                "col2_abs_delta_ppm": [],
                "col3_abs_delta_ppm": [],
                "col2_abs_delta_sum": [],
                "col3_abs_delta_sum": [],
            },
        )
        for field, values in bucket.items():
            try:
                values.append(float(str(row.get(field, ""))))
            except ValueError:
                pass
    result: list[dict[str, str]] = []
    for (variant, label), values_by_field in sorted(grouped.items()):
        result.append(
            {
                "rust_variant": variant,
                "label": label,
                "median_col2_abs_delta_sum": fmt_float(
                    quantile(values_by_field["col2_abs_delta_sum"], 0.5)
                ),
                "median_col2_abs_delta_ppm": fmt_float(
                    quantile(values_by_field["col2_abs_delta_ppm"], 0.5)
                ),
                "median_col3_abs_delta_sum": fmt_float(
                    quantile(values_by_field["col3_abs_delta_sum"], 0.5)
                ),
                "median_col3_abs_delta_ppm": fmt_float(
                    quantile(values_by_field["col3_abs_delta_ppm"], 0.5)
                ),
            }
        )
    return result


def write_report(
    path: Path,
    metadata: dict[str, object],
    summary_rows: list[dict[str, str]],
    comparison_rows: list[dict[str, object]],
) -> None:
    lookup = summary_lookup(summary_rows)
    variants = sorted({row["variant"] for row in summary_rows})
    java_elapsed = median_value(lookup, "java", "elapsed_seconds")
    java_rss = median_value(lookup, "java", "max_rss_kb")

    lines = [
        "# Benchmark Report",
        "",
        f"- Git commit: `{metadata.get('git_commit', '')}`",
        f"- Git branch: `{metadata.get('git_branch', '')}`",
        f"- Host: `{metadata.get('host', '')}`",
        f"- Timestamp UTC: `{metadata.get('timestamp_utc', '')}`",
        "",
        "## Median Summary",
        "",
        "| Variant | Wall s | RSS MiB | Input Counting s | Normalize s | Wall / Java | RSS / Java |",
        "| --- | ---: | ---: | ---: | ---: | ---: | ---: |",
    ]
    for variant in variants:
        elapsed = median_value(lookup, variant, "elapsed_seconds")
        rss = median_value(lookup, variant, "max_rss_kb")
        input_counting = median_value(lookup, variant, "stage_input_counting")
        normalize = median_value(lookup, variant, "stage_normalize")
        if variant == "java":
            input_counting = median_value(lookup, variant, "stage_table_creation")
            normalize = median_value(lookup, variant, "stage_table_read")
        lines.append(
            "| "
            + " | ".join(
                [
                    variant,
                    fmt_float(elapsed) if elapsed is not None else "",
                    fmt_float(rss / 1024.0) if rss is not None else "",
                    fmt_float(input_counting) if input_counting is not None else "",
                    fmt_float(normalize) if normalize is not None else "",
                    ratio_text(elapsed, java_elapsed),
                    ratio_text(rss, java_rss),
                ]
            )
            + " |"
        )

    comparison_summary = comparison_medians(comparison_rows)
    if comparison_summary:
        lines.extend(
            [
                "",
                "## Median Java/Rust Drift",
                "",
                "| Rust Variant | Table | Col2 Abs Delta | Col2 PPM | Col3 Abs Delta | Col3 PPM |",
                "| --- | --- | ---: | ---: | ---: | ---: |",
            ]
        )
        for row in comparison_summary:
            lines.append(
                "| "
                + " | ".join(
                    [
                        row["rust_variant"],
                        row["label"],
                        row["median_col2_abs_delta_sum"],
                        row["median_col2_abs_delta_ppm"],
                        row["median_col3_abs_delta_sum"],
                        row["median_col3_abs_delta_ppm"],
                    ]
                )
                + " |"
            )

    lines.extend(
        [
            "",
            "## Files",
            "",
            "- `metadata.json`: environment, git, input, and setting metadata.",
            "- `commands.tsv`: exact commands for each repeat and variant.",
            "- `raw_runs.tsv`: raw wall-clock, RSS, and selected stage timings.",
            "- `stage_timings.tsv`: long-form stage timing table.",
            "- `comparisons.tsv`: per-repeat Java/Rust histogram drift.",
            "- `summary.tsv`: aggregate min/p10/median/p90/max/mean metrics.",
            "",
        ]
    )
    path.write_text("\n".join(lines), encoding="utf-8")


def main() -> int:
    args = parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)

    if args.report_only:
        metadata_path = args.outdir / "metadata.json"
        summary_path = args.outdir / "summary.tsv"
        comparisons_path = args.outdir / "comparisons.tsv"
        if not metadata_path.exists() or not summary_path.exists():
            raise SystemExit(
                f"--report-only requires existing metadata.json and summary.tsv in {args.outdir}"
            )
        metadata = json.loads(metadata_path.read_text(encoding="utf-8"))
        summary_rows = read_tsv(summary_path)
        comparison_rows: list[dict[str, object]] = (
            list(read_tsv(comparisons_path)) if comparisons_path.exists() else []
        )
        write_report(args.outdir / "report.md", metadata, summary_rows, comparison_rows)
        print(f"Report written to {args.outdir / 'report.md'}")
        return 0

    if not args.skip_build and any(variant.startswith("rust") for variant in args.variants):
        subprocess.run(["cargo", "build", "--release", "--quiet"], cwd=ROOT, check=True)

    metadata = collect_metadata(args)
    (args.outdir / "metadata.json").write_text(
        json.dumps(metadata, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )

    command_rows: list[dict[str, object]] = []
    raw_rows: list[dict[str, object]] = []
    stage_rows: list[dict[str, object]] = []
    comparison_rows: list[dict[str, object]] = []

    for repeat in range(1, args.repeats + 1):
        repeat_java_dir = args.outdir / "runs" / f"{repeat:02d}_java"
        for variant in args.variants:
            run_dir = args.outdir / "runs" / f"{repeat:02d}_{variant}"
            run_dir.mkdir(parents=True, exist_ok=True)
            cmd = variant_command(args, variant, run_dir)
            command_rows.append(
                {
                    "repeat": repeat,
                    "variant": variant,
                    "command": shlex.join(cmd),
                    "run_dir": str(run_dir),
                }
            )
            print(f"repeat {repeat}/{args.repeats}: running {variant}", flush=True)
            started = time.perf_counter()
            result = run_measured(cmd, run_dir, variant, args.timeout, args.sample_interval)
            run_elapsed = time.perf_counter() - started
            stages = parse_stage_timings(Path(result["stderr"]), variant)
            raw_row: dict[str, object] = {
                "repeat": repeat,
                "variant": variant,
                "elapsed_seconds": result["elapsed_seconds"],
                "max_rss_kb": result["max_rss_kb"],
                "status": result["status"],
                "run_elapsed_seconds": f"{run_elapsed:.6f}",
                "run_dir": str(run_dir),
            }
            for stage, seconds in stages.items():
                raw_row[f"stage_{stage}"] = f"{seconds:.6f}"
                stage_rows.append(
                    {
                        "repeat": repeat,
                        "variant": variant,
                        "stage": stage,
                        "seconds": f"{seconds:.6f}",
                    }
                )
            raw_rows.append(raw_row)

        if "java" in args.variants:
            repeat_java_dir = args.outdir / "runs" / f"{repeat:02d}_java"
            for rust_variant in [v for v in args.variants if v.startswith("rust")]:
                rust_dir = args.outdir / "runs" / f"{repeat:02d}_{rust_variant}"
                for label in ["hist", "rhist"]:
                    row = compare_hist(
                        label,
                        repeat_java_dir / f"java.{label}.tsv",
                        rust_dir / f"{rust_variant}.{label}.tsv",
                    )
                    row.update({"repeat": repeat, "rust_variant": rust_variant})
                    comparison_rows.append(row)

    raw_fields = [
        "repeat",
        "variant",
        "elapsed_seconds",
        "max_rss_kb",
        "status",
        "run_elapsed_seconds",
        "stage_input_counting",
        "stage_input_main_counting",
        "stage_input_exact_counting",
        "stage_normalize",
        "stage_input_hist",
        "stage_table_creation",
        "stage_table_read",
        "stage_total",
        "stage_wall_clock",
        "run_dir",
    ]
    write_tsv(args.outdir / "raw_runs.tsv", raw_rows, raw_fields)
    write_tsv(args.outdir / "stage_timings.tsv", stage_rows, ["repeat", "variant", "stage", "seconds"])
    write_tsv(args.outdir / "commands.tsv", command_rows, ["repeat", "variant", "command", "run_dir"])
    write_tsv(
        args.outdir / "comparisons.tsv",
        comparison_rows,
        [
            "repeat",
            "rust_variant",
            "label",
            "status",
            "col2_abs_delta_sum",
            "col2_abs_delta_ppm",
            "col3_abs_delta_sum",
            "col3_abs_delta_ppm",
            "first_diff_key",
        ],
    )
    summary_rows = aggregate(raw_rows)
    write_tsv(
        args.outdir / "summary.tsv",
        summary_rows,
        ["variant", "metric", "n", "min", "p10", "median", "p90", "max", "mean"],
    )
    write_report(args.outdir / "report.md", metadata, summary_rows, comparison_rows)

    print(f"\nBenchmark baseline written to {args.outdir}")
    print("Key aggregate rows:")
    for row in summary_rows:
        if row["metric"] in {"elapsed_seconds", "stage_input_counting", "stage_normalize", "max_rss_kb"}:
            print(
                f"{row['variant']}\t{row['metric']}\tn={row['n']}"
                f"\tmedian={row['median']}\tp10={row['p10']}\tp90={row['p90']}"
            )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

#!/usr/bin/env python3
"""Run the BBNorm Java/Rust parity acceptance matrix.

The matrix is intentionally a thin coordinator over the existing benchmark
harnesses. It keeps the policy layer separate from the low-level artifact
collection so every row can be classified as exact, bounded-drift, or accepted
Rust-over-Java divergence.
"""

from __future__ import annotations

import csv
import os
import subprocess
import sys
from datetime import datetime
from pathlib import Path
from typing import Iterable


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_MANIFEST = ROOT / "docs" / "parity_acceptance_matrix.tsv"
DEFAULT_OUT = ROOT / "tmp" / f"parity_acceptance_matrix_{datetime.now():%Y%m%d_%H%M%S}"

SUMMARY_FIELDS = [
    "row_id",
    "dataset",
    "mode",
    "expectation",
    "pass",
    "verdict",
    "reason",
    "java_status",
    "rust_status",
    "hist_cmp",
    "rhist_cmp",
    "sequence_cmp",
    "drift_gate",
    "drift_classification",
    "hist_raw_delta_ppm",
    "hist_unique_delta_ppm",
    "rhist_reads_delta_ppm",
    "rhist_bases_delta_ppm",
    "java_seconds",
    "rust_seconds",
    "java_rss_kb",
    "rust_rss_kb",
    "rust_input_counting_s",
    "rust_input_main_counting_s",
    "rust_input_prefilter_counting_s",
    "rust_normalize_s",
    "artifact_dir",
    "notes",
]


def read_manifest(path: Path) -> list[dict[str, str]]:
    lines = [
        line
        for line in path.read_text(encoding="utf-8").splitlines()
        if line.strip() and not line.lstrip().startswith("#")
    ]
    if not lines:
        return []
    return list(csv.DictReader(lines, delimiter="\t"))


def wanted_rows(rows: Iterable[dict[str, str]]) -> list[dict[str, str]]:
    requested = {
        part.strip()
        for part in os.environ.get("ROW_CASES", "").replace(",", " ").split()
        if part.strip()
    }
    include_disabled = os.environ.get("MATRIX_INCLUDE_DISABLED", "0") == "1"
    selected = []
    for row in rows:
        row_id = row.get("row_id", "")
        if requested and row_id not in requested:
            continue
        if not include_disabled and row.get("enabled", "1") != "1":
            continue
        selected.append(row)
    return selected


def row_env(row: dict[str, str]) -> dict[str, str]:
    env = os.environ.copy()
    mapping = {
        "modes": "MODE_CASES",
        "drift_gate_profile": "DRIFT_GATE_PROFILE",
        "require_identical_comparisons": "REQUIRE_IDENTICAL_COMPARISONS",
        "expected_failure_modes": "EXPECTED_FAILURE_MODES",
        "skip_expected_failure_java": "SKIP_EXPECTED_FAILURE_JAVA",
        "reads": "READS",
        "table_reads": "TABLE_READS",
        "threads": "THREADS",
        "zipthreads": "ZIPTHREADS",
        "java_xmx": "JAVA_XMX",
        "mem": "MEM",
        "rust_mem_auto_from_java": "RUST_MEM_AUTO_FROM_JAVA",
        "rust_mem_auto_max_bytes": "RUST_MEM_AUTO_MAX_BYTES",
        "write_outputs": "WRITE_OUTPUTS",
        "java_max_rss_kb": "JAVA_MAX_RSS_KB",
        "rust_max_rss_kb": "RUST_MAX_RSS_KB",
        "timeout": "TIMEOUT",
    }
    for source, target in mapping.items():
        value = row.get(source, "")
        if value != "":
            env[target] = value
    env["ALLOW_MODE_FAILURES"] = "1"
    return env


def read_tsv(path: Path) -> list[dict[str, str]]:
    if not path.exists():
        return []
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def read_key_value_tsv(path: Path) -> dict[str, str]:
    if not path.exists():
        return {}
    values: dict[str, str] = {}
    for line in path.read_text(encoding="utf-8", errors="replace").splitlines():
        if not line:
            continue
        key, _, value = line.partition("\t")
        values[key] = value
    return values


def sequence_cmp(mode_artifact: Path, require_sequence: bool) -> str:
    if not require_sequence:
        return ""
    comparisons = read_key_value_tsv(mode_artifact / "comparison.tsv")
    suffixes = ["keep1.fq.gz", "toss1.fq.gz", "keep2.fq.gz", "toss2.fq.gz"]
    present = [
        (suffix, comparisons.get(f"{suffix}_cmp", ""))
        for suffix in suffixes
        if f"{suffix}_cmp" in comparisons
    ]
    if not present:
        return "missing"
    if all(value == "identical" for _, value in present):
        return "identical"
    return ",".join(f"{suffix}={value}" for suffix, value in present)


def classify(
    row: dict[str, str], mode_row: dict[str, str], mode_artifact: Path
) -> tuple[bool, str, str, str]:
    expectation = row.get("expectation", "exact")
    require_sequence = row.get("require_sequence_outputs", "0") == "1"
    seq_cmp = sequence_cmp(mode_artifact, require_sequence)
    java_status = mode_row.get("java_status", "")
    rust_status = mode_row.get("rust_status", "")
    hist_cmp = mode_row.get("hist_cmp", "")
    rhist_cmp = mode_row.get("rhist_cmp", "")
    drift_gate = mode_row.get("drift_gate", "")
    drift_classification = mode_row.get("drift_classification", "")
    status = mode_row.get("status", "")

    if rust_status != "0":
        return False, "fail_rust", f"rust_status={rust_status}", seq_cmp

    if expectation == "exact":
        reasons = []
        if status != "0":
            reasons.append(f"mode_status={status}")
        if java_status != "0":
            reasons.append(f"java_status={java_status}")
        if hist_cmp != "identical":
            reasons.append(f"hist_cmp={hist_cmp}")
        if rhist_cmp != "identical":
            reasons.append(f"rhist_cmp={rhist_cmp}")
        if require_sequence and seq_cmp != "identical":
            reasons.append(f"sequence_cmp={seq_cmp}")
        if reasons:
            return False, "fail_exact", ",".join(reasons), seq_cmp
        return True, "exact_match", "", seq_cmp

    if expectation == "bounded":
        if java_status != "0":
            return False, "fail_bounded", f"java_status={java_status}", seq_cmp
        if drift_gate == "fail":
            return (
                False,
                "fail_bounded",
                mode_row.get("drift_gate_reason", "drift_gate=fail"),
                seq_cmp,
            )
        if hist_cmp == "identical" and rhist_cmp == "identical":
            return True, "exact_match", "", seq_cmp
        return True, "bounded_drift", drift_classification or drift_gate, seq_cmp

    if expectation == "rust_better_java_crashes":
        if java_status != "0" or drift_classification == "skipped_java_failed":
            return True, "rust_better_java_crashes", "", seq_cmp
        return False, "fail_expected_java_failure", "java_completed", seq_cmp

    if expectation == "unsupported":
        return True, "unsupported", "declared unsupported", seq_cmp

    return False, "fail_manifest", f"unknown expectation={expectation}", seq_cmp


def missing_dataset_row(row: dict[str, str], reason: str) -> dict[str, str]:
    return {
        "row_id": row.get("row_id", ""),
        "dataset": row.get("dataset", ""),
        "mode": "",
        "expectation": row.get("expectation", ""),
        "pass": "0" if os.environ.get("STRICT_DATASETS", "0") == "1" else "1",
        "verdict": "skipped_missing_dataset",
        "reason": reason,
        "notes": row.get("notes", ""),
    }


def run_row(row: dict[str, str], outdir: Path) -> list[dict[str, str]]:
    row_id = row["row_id"]
    r1 = ROOT / row["r1"]
    r2_text = row.get("r2", "")
    r2 = ROOT / r2_text if r2_text else None
    if not r1.exists():
        return [missing_dataset_row(row, f"missing r1: {row['r1']}")]
    if r2 is not None and not r2.exists():
        return [missing_dataset_row(row, f"missing r2: {r2_text}")]

    row_out = outdir / row_id
    row_out.parent.mkdir(parents=True, exist_ok=True)
    cmd = [
        str(ROOT / "scripts" / "benchmark_java_rust_modes.sh"),
        str(r1),
        str(r2 or ""),
        str(row_out),
    ]
    print(f"Running matrix row {row_id}: modes={row.get('modes', '')}", flush=True)
    completed = subprocess.run(cmd, cwd=ROOT, env=row_env(row), check=False)
    mode_rows = read_tsv(row_out / "summary.tsv")
    if not mode_rows:
        return [
            {
                "row_id": row_id,
                "dataset": row.get("dataset", ""),
                "mode": "",
                "expectation": row.get("expectation", ""),
                "pass": "0",
                "verdict": "fail_harness",
                "reason": f"harness_exit={completed.returncode}; missing summary.tsv",
                "artifact_dir": str(row_out),
                "notes": row.get("notes", ""),
            }
        ]

    results = []
    for mode_row in mode_rows:
        mode = mode_row.get("mode", "")
        mode_artifact = Path(mode_row.get("artifact_dir") or row_out / mode)
        passed, verdict, reason, seq_cmp = classify(row, mode_row, mode_artifact)
        results.append(
            {
                "row_id": row_id,
                "dataset": row.get("dataset", ""),
                "mode": mode,
                "expectation": row.get("expectation", ""),
                "pass": "1" if passed else "0",
                "verdict": verdict,
                "reason": reason,
                "java_status": mode_row.get("java_status", ""),
                "rust_status": mode_row.get("rust_status", ""),
                "hist_cmp": mode_row.get("hist_cmp", ""),
                "rhist_cmp": mode_row.get("rhist_cmp", ""),
                "sequence_cmp": seq_cmp,
                "drift_gate": mode_row.get("drift_gate", ""),
                "drift_classification": mode_row.get("drift_classification", ""),
                "hist_raw_delta_ppm": mode_row.get("hist_raw_delta_ppm", ""),
                "hist_unique_delta_ppm": mode_row.get("hist_unique_delta_ppm", ""),
                "rhist_reads_delta_ppm": mode_row.get("rhist_reads_delta_ppm", ""),
                "rhist_bases_delta_ppm": mode_row.get("rhist_bases_delta_ppm", ""),
                "java_seconds": mode_row.get("java_seconds", ""),
                "rust_seconds": mode_row.get("rust_seconds", ""),
                "java_rss_kb": mode_row.get("java_rss_kb", ""),
                "rust_rss_kb": mode_row.get("rust_rss_kb", ""),
                "rust_input_counting_s": mode_row.get("rust_input_counting_s", ""),
                "rust_input_main_counting_s": mode_row.get("rust_input_main_counting_s", ""),
                "rust_input_prefilter_counting_s": mode_row.get(
                    "rust_input_prefilter_counting_s", ""
                ),
                "rust_normalize_s": mode_row.get("rust_normalize_s", ""),
                "artifact_dir": str(mode_artifact),
                "notes": row.get("notes", ""),
            }
        )
    return results


def main(argv: list[str]) -> int:
    manifest = DEFAULT_MANIFEST
    outdir = DEFAULT_OUT
    if len(argv) > 1:
        first = Path(argv[1])
        first_resolved = first if first.is_absolute() else ROOT / first
        if first_resolved.is_file():
            manifest = first
            if len(argv) > 2:
                outdir = Path(argv[2])
        else:
            outdir = first
    if len(argv) > 3:
        print(
            "Usage: scripts/parity_acceptance_matrix.py [manifest.tsv] [outdir]",
            file=sys.stderr,
        )
        return 2
    if not manifest.is_absolute():
        manifest = ROOT / manifest
    if not outdir.is_absolute():
        outdir = ROOT / outdir
    rows = wanted_rows(read_manifest(manifest))
    if not rows:
        print(f"No matrix rows selected from {manifest}", file=sys.stderr)
        return 2
    outdir.mkdir(parents=True, exist_ok=True)
    all_results: list[dict[str, str]] = []
    for row in rows:
        all_results.extend(run_row(row, outdir))

    summary_path = outdir / "acceptance_summary.tsv"
    with summary_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=SUMMARY_FIELDS, delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        for result in all_results:
            writer.writerow({field: result.get(field, "") for field in SUMMARY_FIELDS})

    print(summary_path)
    with summary_path.open(encoding="utf-8") as handle:
        sys.stdout.write(handle.read())

    failures = [row for row in all_results if row.get("pass") != "1"]
    if failures and os.environ.get("ALLOW_MATRIX_FAILURES", "0") != "1":
        print(f"{len(failures)} matrix verdict(s) failed", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv))

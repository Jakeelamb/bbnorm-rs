#!/usr/bin/env python3
"""Profile one exact trustworthy-baseline command shape.

This intentionally reuses benchmark_trustworthy_baseline.variant_command so
profiling cannot silently drift away from the benchmark harness knobs.
"""

from __future__ import annotations

import argparse
import shutil
import shlex
import subprocess
import sys
from pathlib import Path

import benchmark_trustworthy_baseline as baseline


ROOT = Path(__file__).resolve().parents[1]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--variant", default="rust_deterministic")
    parser.add_argument("--outdir", type=Path, required=True)
    parser.add_argument("--tool", choices=["perf", "flamegraph"], default="perf")
    parser.add_argument("--frequency", default="499")
    parser.add_argument("--flamegraph-output", default="")
    parser.add_argument("--skip-build", action="store_true")
    parser.add_argument(
        "baseline_args",
        nargs=argparse.REMAINDER,
        help="Arguments forwarded to benchmark_trustworthy_baseline.py after --",
    )
    args = parser.parse_args()
    if args.baseline_args and args.baseline_args[0] == "--":
        args.baseline_args = args.baseline_args[1:]
    if args.variant == "java" and args.tool == "flamegraph":
        parser.error("--tool flamegraph only supports Rust variants")
    allowed = {"java", "rust_deterministic", "rust_nondeterministic"}
    if args.variant not in allowed:
        parser.error(f"--variant must be one of: {', '.join(sorted(allowed))}")
    return args


def baseline_namespace(forwarded: list[str], outdir: Path, variant: str) -> argparse.Namespace:
    old_argv = sys.argv
    try:
        sys.argv = [
            "benchmark_trustworthy_baseline.py",
            "--outdir",
            str(outdir),
            "--repeats",
            "1",
            "--variants",
            variant,
            *forwarded,
        ]
        return baseline.parse_args()
    finally:
        sys.argv = old_argv


def write_command(path: Path, command: list[str]) -> None:
    path.write_text(shlex.join(command) + "\n", encoding="utf-8")


def timeout_prefix(timeout: str) -> list[str]:
    lowered = timeout.strip().lower()
    if lowered in {"", "off", "none", "0"}:
        return []
    return ["timeout", timeout]


def main() -> int:
    args = parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)

    profile_run_dir = args.outdir / "profile_run"
    profile_run_dir.mkdir(parents=True, exist_ok=True)
    baseline_args = baseline_namespace(args.baseline_args, args.outdir, args.variant)

    if not args.skip_build and args.variant.startswith("rust"):
        subprocess.run(["cargo", "build", "--release", "--quiet"], cwd=ROOT, check=True)

    command = baseline.variant_command(baseline_args, args.variant, profile_run_dir)
    write_command(args.outdir / "profile_command.txt", command)

    if args.tool == "perf":
        perf_data = args.outdir / f"{args.variant}.perf.data"
        profile_cmd = [
            *timeout_prefix(baseline_args.timeout),
            "perf",
            "record",
            "-F",
            args.frequency,
            "-g",
            "--output",
            str(perf_data),
            "--",
            *command,
        ]
        write_command(args.outdir / "profile_invocation.txt", profile_cmd)
        return subprocess.run(profile_cmd, cwd=ROOT, check=False).returncode

    flamegraph_output = (
        Path(args.flamegraph_output)
        if args.flamegraph_output
        else args.outdir / f"{args.variant}.flamegraph.svg"
    )
    profile_cmd = [
        *timeout_prefix(baseline_args.timeout),
        "env",
        "CARGO_PROFILE_RELEASE_DEBUG=true",
        "cargo",
        "flamegraph",
        "--freq",
        args.frequency,
        "--output",
        str(flamegraph_output),
        "--bin",
        "bbnorm-rs",
        "--",
        *command[1:],
    ]
    write_command(args.outdir / "profile_invocation.txt", profile_cmd)
    completed = subprocess.run(profile_cmd, cwd=ROOT, check=False)
    flamegraph_perf_data = ROOT / "perf.data"
    if flamegraph_perf_data.exists():
        shutil.move(
            str(flamegraph_perf_data),
            str(args.outdir / f"{args.variant}.flamegraph.perf.data"),
        )
    return completed.returncode


if __name__ == "__main__":
    raise SystemExit(main())

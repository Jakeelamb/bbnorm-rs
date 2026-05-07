#!/usr/bin/env python3
"""Extract coarse Java/Rust BBNorm stage timings from stderr logs."""
from __future__ import annotations

import csv
import re
import sys
from pathlib import Path

RUST_STAGE_RE = re.compile(r"Stage timing:\s*name=([^;]+);\s*seconds=([0-9.]+)")
JAVA_PATTERNS = [
    ("table_creation", re.compile(r"Table creation time:\s*([0-9.]+)\s+seconds")),
    ("table_read", re.compile(r"Table read time:\s*([0-9.]+)\s+seconds")),
    ("total", re.compile(r"Total time:\s*([0-9.]+)\s+seconds")),
]
MEASURE_WALL_RE = re.compile(r"Elapsed \(wall clock\) time \(seconds\):\s*([0-9.]+)")


def parse_spec(spec: str) -> tuple[str, Path]:
    if "=" not in spec:
        raise SystemExit(f"Expected tool=stderr.log, got {spec!r}")
    tool, path = spec.split("=", 1)
    if not tool:
        raise SystemExit(f"Missing tool name in {spec!r}")
    return tool, Path(path)


def emit_rows(tool: str, path: Path) -> list[dict[str, str]]:
    if not path.exists():
        return []
    text = path.read_text(encoding="utf-8", errors="replace")
    rows: list[dict[str, str]] = []
    for match in RUST_STAGE_RE.finditer(text):
        rows.append({"tool": tool, "stage": match.group(1), "seconds": match.group(2)})
    if rows:
        return rows
    for stage, pattern in JAVA_PATTERNS:
        if match := pattern.search(text):
            rows.append({"tool": tool, "stage": stage, "seconds": match.group(1)})
    if match := MEASURE_WALL_RE.search(text):
        rows.append({"tool": tool, "stage": "wall_clock", "seconds": match.group(1)})
    return rows


def main(argv: list[str]) -> int:
    if len(argv) < 2:
        raise SystemExit("Usage: extract_stage_timings.py tool=stderr.log [tool=stderr.log ...]")
    writer = csv.DictWriter(
        sys.stdout,
        fieldnames=["tool", "stage", "seconds"],
        delimiter="\t",
        lineterminator="\n",
    )
    writer.writeheader()
    for spec in argv[1:]:
        tool, path = parse_spec(spec)
        writer.writerows(emit_rows(tool, path))
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv))

#!/usr/bin/env python3
"""Run a command while recording elapsed seconds and peak RSS.

This intentionally avoids requiring GNU /usr/bin/time. It samples the process
tree through /proc so benchmark scripts still get useful RSS numbers on systems
where package install or GNU time is unavailable.
"""

from __future__ import annotations

import argparse
import os
import resource
import signal
import subprocess
import sys
import time
from pathlib import Path


def parse_timeout(value: str | None) -> float | None:
    if not value:
        return None
    raw = value.strip().lower()
    multiplier = 1.0
    if raw.endswith("ms"):
        raw = raw[:-2]
        multiplier = 0.001
    elif raw.endswith("s"):
        raw = raw[:-1]
    elif raw.endswith("m"):
        raw = raw[:-1]
        multiplier = 60.0
    elif raw.endswith("h"):
        raw = raw[:-1]
        multiplier = 3600.0
    elif raw.endswith("d"):
        raw = raw[:-1]
        multiplier = 86400.0
    return float(raw) * multiplier


def children_of(pid: int) -> list[int]:
    task_dir = Path(f"/proc/{pid}/task")
    children: list[int] = []
    try:
        tasks = list(task_dir.iterdir())
    except OSError:
        return children
    for task in tasks:
        try:
            data = (task / "children").read_text(encoding="utf-8").strip()
        except OSError:
            continue
        if data:
            children.extend(int(item) for item in data.split())
    return children


def process_tree(root_pid: int) -> set[int]:
    seen: set[int] = set()
    stack = [root_pid]
    while stack:
        pid = stack.pop()
        if pid in seen:
            continue
        seen.add(pid)
        stack.extend(child for child in children_of(pid) if child not in seen)
    return seen


def rss_kb(pid: int) -> int:
    try:
        with open(f"/proc/{pid}/statm", "r", encoding="utf-8") as handle:
            fields = handle.read().split()
    except OSError:
        return 0
    if len(fields) < 2:
        return 0
    try:
        resident_pages = int(fields[1])
    except ValueError:
        return 0
    return resident_pages * (os.sysconf("SC_PAGE_SIZE") // 1024)


def tree_rss_kb(root_pid: int) -> int:
    return sum(rss_kb(pid) for pid in process_tree(root_pid))


def terminate_group(proc: subprocess.Popen[object]) -> None:
    try:
        os.killpg(proc.pid, signal.SIGTERM)
    except OSError:
        return
    try:
        proc.wait(timeout=5)
        return
    except subprocess.TimeoutExpired:
        pass
    try:
        os.killpg(proc.pid, signal.SIGKILL)
    except OSError:
        pass
    proc.wait()


def open_output(path: str | None):
    if path is None:
        return None
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    return open(path, "wb")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--metrics", required=True, help="TSV file for elapsed, RSS, and status")
    parser.add_argument("--stdout", help="File for command stdout")
    parser.add_argument("--stderr", help="File for command stderr and monitor footer")
    parser.add_argument("--timeout", help="Kill the command after this duration, e.g. 30m or 4h")
    parser.add_argument(
        "--max-rss-kb",
        type=int,
        default=0,
        help="Kill the process tree if sampled RSS exceeds this many kilobytes",
    )
    parser.add_argument("--sample-interval", type=float, default=0.05)
    parser.add_argument("cmd", nargs=argparse.REMAINDER)
    args = parser.parse_args()
    if args.max_rss_kb < 0:
        parser.error("--max-rss-kb must be non-negative")

    cmd = args.cmd
    if cmd and cmd[0] == "--":
        cmd = cmd[1:]
    if not cmd:
        parser.error("missing command after --")

    timeout_seconds = parse_timeout(args.timeout)
    metrics_path = Path(args.metrics)
    metrics_path.parent.mkdir(parents=True, exist_ok=True)

    stdout_handle = open_output(args.stdout)
    stderr_handle = open_output(args.stderr)
    start = time.perf_counter()
    max_rss = 0
    timed_out = False
    rss_exceeded = False

    try:
        proc = subprocess.Popen(
            cmd,
            stdout=stdout_handle if stdout_handle is not None else None,
            stderr=stderr_handle if stderr_handle is not None else None,
            start_new_session=True,
        )
        while True:
            current_rss = tree_rss_kb(proc.pid)
            max_rss = max(max_rss, current_rss)
            status = proc.poll()
            elapsed = time.perf_counter() - start
            if status is not None:
                break
            if args.max_rss_kb and current_rss > args.max_rss_kb:
                rss_exceeded = True
                terminate_group(proc)
                status = 125
                break
            if timeout_seconds is not None and elapsed >= timeout_seconds:
                timed_out = True
                terminate_group(proc)
                status = 124
                break
            time.sleep(max(args.sample_interval, 0.001))

        max_rss = max(max_rss, tree_rss_kb(proc.pid), resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss)
        elapsed = time.perf_counter() - start
        status = 125 if rss_exceeded else 124 if timed_out else int(status if status is not None else proc.returncode)
    finally:
        if stdout_handle is not None:
            stdout_handle.close()
        if stderr_handle is not None:
            stderr_handle.close()

    with metrics_path.open("w", encoding="utf-8") as handle:
        handle.write(f"{elapsed:.6f}\t{max_rss}\t{status}\n")

    if args.stderr:
        with open(args.stderr, "a", encoding="utf-8") as handle:
            handle.write(f"Elapsed (wall clock) time (seconds): {elapsed:.6f}\n")
            handle.write(f"Maximum resident set size (kbytes): {max_rss}\n")
            handle.write(f"Exit status: {status}\n")
            handle.write(f"Timed out: {'true' if timed_out else 'false'}\n")
            handle.write(f"RSS guard limit (kbytes): {args.max_rss_kb or 'unlimited'}\n")
            handle.write(f"RSS guard exceeded: {'true' if rss_exceeded else 'false'}\n")

    return status


if __name__ == "__main__":
    sys.exit(main())

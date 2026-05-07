#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT"

DATA1="vendor/BBTools-master/resources/sample1.fq.gz"
DATA2="vendor/BBTools-master/resources/sample2.fq.gz"
OUT="${1:-tmp/thread_scaling_benchmark}"
THREAD_CASES="${THREAD_CASES:-1 2 auto}"
PYTHON="${PYTHON:-python}"
mkdir -p "$OUT"
rm -f "$OUT"/*

COMMON=(
  "in=$DATA1"
  "in2=$DATA2"
  "passes=1"
  "keepall=t"
  "lowbindepth=1"
  "highbindepth=1"
  "k=31"
  "minq=0"
  "minprob=0"
  "min=0"
  "minkmers=1"
  "target=999999999"
  "max=999999999"
  "overwrite=t"
  "bits=32"
)

measure() {
  local metrics_file="$1"
  shift
  "$PYTHON" - "$metrics_file" "$@" <<'PY'
import resource
import subprocess
import sys
import time

metrics_file = sys.argv[1]
cmd = sys.argv[2:]
start = time.perf_counter()
proc = subprocess.run(cmd)
elapsed = time.perf_counter() - start
usage = resource.getrusage(resource.RUSAGE_CHILDREN)
with open(metrics_file, "w", encoding="utf-8") as handle:
    handle.write(f"{elapsed:.6f}\t{usage.ru_maxrss}\t{proc.returncode}\n")
sys.exit(proc.returncode)
PY
}

write_result() {
  local label="$1"
  local metrics_file="$2"
  local elapsed rss status
  read -r elapsed rss status < "$metrics_file"
  printf '%s\t%s\t%s\t%s\n' "$label" "$elapsed" "$rss" "$status" >> "$OUT/results.tsv"
}

label_for_threads() {
  local label="threads_$1"
  printf '%s\n' "${label//[^A-Za-z0-9_]/_}"
}

run_case() {
  local threads="$1"
  local label
  label="$(label_for_threads "$threads")"
  printf 'Running Rust bbnorm-rs with threads=%s...\n' "$threads"
  measure "$OUT/${label}.metrics.tsv" \
    target/release/bbnorm-rs \
    "${COMMON[@]}" \
    "threads=$threads" \
    "out=$OUT/${label}.keep1.fq" "out2=$OUT/${label}.keep2.fq" \
    "outlow=$OUT/${label}.low1.fq" "outlow2=$OUT/${label}.low2.fq" \
    "outmid=$OUT/${label}.mid1.fq" "outmid2=$OUT/${label}.mid2.fq" \
    "outhigh=$OUT/${label}.high1.fq" "outhigh2=$OUT/${label}.high2.fq" \
    "hist=$OUT/${label}.hist.tsv" "rhist=$OUT/${label}.rhist.tsv" \
    >"$OUT/${label}.stdout.log" 2>"$OUT/${label}.stderr.log"
  write_result "$label" "$OUT/${label}.metrics.tsv"
}

compare_to_baseline() {
  local baseline="$1"
  local label="$2"
  for suffix in keep1.fq keep2.fq low1.fq low2.fq mid1.fq mid2.fq high1.fq high2.fq hist.tsv rhist.tsv; do
    cmp "$OUT/$baseline.$suffix" "$OUT/$label.$suffix"
  done
}

cargo build --release --quiet

printf 'thread_case\telapsed_seconds\tmax_rss_kb\tstatus\n' > "$OUT/results.tsv"
for threads in $THREAD_CASES; do
  run_case "$threads"
done

baseline_label=""
for threads in $THREAD_CASES; do
  label="$(label_for_threads "$threads")"
  if [[ -z "$baseline_label" ]]; then
    baseline_label="$label"
  fi
  compare_to_baseline "$baseline_label" "$label"
done

"$PYTHON" - "$OUT/results.tsv" <<'PY'
import csv
import sys

rows = list(csv.DictReader(open(sys.argv[1], encoding="utf-8"), delimiter="\t"))
print("\nThread scaling summary:")
print("thread_case\telapsed_seconds\tmax_rss_kb")
for row in rows:
    print(f"{row['thread_case']}\t{float(row['elapsed_seconds']):.6f}\t{row['max_rss_kb']}")
fastest = min(rows, key=lambda row: float(row["elapsed_seconds"]))
lowest_rss = min(rows, key=lambda row: int(row["max_rss_kb"]))
print(
    f"\nFastest: {fastest['thread_case']} "
    f"({float(fastest['elapsed_seconds']):.6f}s)"
)
print(
    f"Lowest RSS: {lowest_rss['thread_case']} "
    f"({int(lowest_rss['max_rss_kb'])} KB)"
)
PY

printf '\nThread scaling benchmark passed. Per-run metrics: %s\n' "$OUT/results.tsv"
printf 'Outputs and logs: %s\n' "$OUT"

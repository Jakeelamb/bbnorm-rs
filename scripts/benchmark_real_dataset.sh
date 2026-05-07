#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT"

DATA1="vendor/BBTools-master/resources/sample1.fq.gz"
DATA2="vendor/BBTools-master/resources/sample2.fq.gz"
CP="vendor/BBTools-master/current"
OUT="${1:-tmp/real_benchmark}"
REPEATS="${REPEATS:-3}"
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
  "threads=1"
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
  local tool="$1"
  local iter="$2"
  local metrics_file="$3"
  local elapsed rss status
  read -r elapsed rss status < "$metrics_file"
  printf '%s\t%s\t%s\t%s\t%s\n' "$tool" "$iter" "$elapsed" "$rss" "$status" >> "$OUT/results.tsv"
}

compare_iteration_outputs() {
  local iter="$1"
  for suffix in keep1.fq keep2.fq low1.fq low2.fq mid1.fq mid2.fq high1.fq high2.fq hist.tsv rhist.tsv; do
    cmp "$OUT/java.${iter}.${suffix}" "$OUT/rust.${iter}.${suffix}"
  done
}

cargo build --release --quiet

printf 'tool\titeration\telapsed_seconds\tmax_rss_kb\tstatus\n' > "$OUT/results.tsv"

for iter in $(seq 1 "$REPEATS"); do
  printf 'Benchmark iteration %s/%s: Java BBNorm...\n' "$iter" "$REPEATS"
  measure "$OUT/java.${iter}.metrics.tsv" \
    java -Xmx1g -cp "$CP" jgi.KmerNormalize \
    "${COMMON[@]}" \
    "out=$OUT/java.${iter}.keep1.fq" "out2=$OUT/java.${iter}.keep2.fq" \
    "outlow=$OUT/java.${iter}.low1.fq" "outlow2=$OUT/java.${iter}.low2.fq" \
    "outmid=$OUT/java.${iter}.mid1.fq" "outmid2=$OUT/java.${iter}.mid2.fq" \
    "outhigh=$OUT/java.${iter}.high1.fq" "outhigh2=$OUT/java.${iter}.high2.fq" \
    "hist=$OUT/java.${iter}.hist.tsv" "rhist=$OUT/java.${iter}.rhist.tsv" \
    >"$OUT/java.${iter}.stdout.log" 2>"$OUT/java.${iter}.stderr.log"
  write_result java "$iter" "$OUT/java.${iter}.metrics.tsv"

  printf 'Benchmark iteration %s/%s: Rust bbnorm-rs...\n' "$iter" "$REPEATS"
  measure "$OUT/rust.${iter}.metrics.tsv" \
    target/release/bbnorm-rs \
    "${COMMON[@]}" \
    "out=$OUT/rust.${iter}.keep1.fq" "out2=$OUT/rust.${iter}.keep2.fq" \
    "outlow=$OUT/rust.${iter}.low1.fq" "outlow2=$OUT/rust.${iter}.low2.fq" \
    "outmid=$OUT/rust.${iter}.mid1.fq" "outmid2=$OUT/rust.${iter}.mid2.fq" \
    "outhigh=$OUT/rust.${iter}.high1.fq" "outhigh2=$OUT/rust.${iter}.high2.fq" \
    "hist=$OUT/rust.${iter}.hist.tsv" "rhist=$OUT/rust.${iter}.rhist.tsv" \
    >"$OUT/rust.${iter}.stdout.log" 2>"$OUT/rust.${iter}.stderr.log"
  write_result rust "$iter" "$OUT/rust.${iter}.metrics.tsv"

  compare_iteration_outputs "$iter"
done

"$PYTHON" - "$OUT/results.tsv" <<'PY'
from collections import defaultdict
import csv
import sys

path = sys.argv[1]
rows = list(csv.DictReader(open(path, encoding="utf-8"), delimiter="\t"))
by_tool = defaultdict(list)
for row in rows:
    by_tool[row["tool"]].append(row)
print("\nBenchmark summary:")
print("tool\trepeats\tavg_seconds\tbest_seconds\tmax_rss_kb")
for tool in sorted(by_tool):
    data = by_tool[tool]
    seconds = [float(row["elapsed_seconds"]) for row in data]
    rss = [int(row["max_rss_kb"]) for row in data]
    print(f"{tool}\t{len(data)}\t{sum(seconds)/len(seconds):.6f}\t{min(seconds):.6f}\t{max(rss)}")
PY

printf '\nBenchmark parity passed. Per-run metrics: %s\n' "$OUT/results.tsv"
printf 'Outputs and logs: %s\n' "$OUT"

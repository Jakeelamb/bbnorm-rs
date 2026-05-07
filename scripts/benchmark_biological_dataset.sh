#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT"

BIO_ROOT="${BIO_ROOT:-/home/jake/Projects/biological data}"
DATA1="${DATA1:-$BIO_ROOT/reads/short_reads/scer_s288c_pe_srr23631023/SRR23631023_1.fastq.gz}"
DATA2="${DATA2:-$BIO_ROOT/reads/short_reads/scer_s288c_pe_srr23631023/SRR23631023_2.fastq.gz}"
OUT="${1:-tmp/biological_dataset_benchmark}"
THREAD_CASES="${THREAD_CASES:-1 2 auto}"
READS="${READS:-250000}"
TABLE_READS="${TABLE_READS:-$READS}"
KEEP_OUTPUTS="${KEEP_OUTPUTS:-1}"
MAX_RSS_KB="${MAX_RSS_KB:-}"
TIMEOUT="${TIMEOUT:-}"
PYTHON="${PYTHON:-python3}"
MEASURE_SCRIPT="${MEASURE_SCRIPT:-scripts/measure_command.py}"
EXTRA_ARGS="${EXTRA_ARGS:-}"
mkdir -p "$OUT"
rm -f "$OUT"/*

cleanup_outputs() {
  if [[ "$KEEP_OUTPUTS" == "0" || "$KEEP_OUTPUTS" == "false" || "$KEEP_OUTPUTS" == "f" ]]; then
    rm -f "$OUT"/*.fq "$OUT"/*.fastq
  fi
}
trap cleanup_outputs EXIT

if [[ ! -f "$DATA1" || ! -f "$DATA2" ]]; then
  printf 'Missing paired dataset files:\n  DATA1=%s\n  DATA2=%s\n' "$DATA1" "$DATA2" >&2
  exit 2
fi

case "${MAX_RSS_KB,,}" in
  ""|"0"|"none"|"off"|"unlimited")
    MAX_RSS_KB=
    ;;
esac
if [[ -n "$MAX_RSS_KB" && ! "$MAX_RSS_KB" =~ ^[0-9]+$ ]]; then
  printf 'MAX_RSS_KB must be an integer number of KB, or 0/off/none/unlimited; got %s\n' "$MAX_RSS_KB" >&2
  exit 2
fi

COMMON=(
  "in=$DATA1"
  "in2=$DATA2"
  "passes=1"
  "reads=$READS"
  "tablereads=$TABLE_READS"
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

EXTRA_ARGS_ARRAY=()
if [[ -n "$EXTRA_ARGS" ]]; then
  read -r -a EXTRA_ARGS_ARRAY <<< "$EXTRA_ARGS"
fi

measure() {
  local metrics_file="$1"
  local stdout_file="$2"
  local stderr_file="$3"
  shift 3
  local measure_args=(
    "$PYTHON"
    "$MEASURE_SCRIPT"
    "--metrics" "$metrics_file"
    "--stdout" "$stdout_file"
    "--stderr" "$stderr_file"
  )
  if [[ -n "$TIMEOUT" ]]; then
    measure_args+=("--timeout" "$TIMEOUT")
  fi
  if [[ -n "$MAX_RSS_KB" ]]; then
    measure_args+=("--max-rss-kb" "$MAX_RSS_KB")
  fi
  measure_args+=("--" "$@")
  "${measure_args[@]}"
}

write_result() {
  local label="$1"
  local metrics_file="$2"
  local elapsed rss status
  read -r elapsed rss status < "$metrics_file"
  printf '%s\t%s\t%s\t%s\n' "$label" "$elapsed" "$rss" "$status" >> "$OUT/results.tsv"
}

enforce_max_rss() {
  local label="$1"
  local metrics_file="$2"
  local elapsed rss status
  if [[ -z "$MAX_RSS_KB" ]]; then
    return 0
  fi
  read -r elapsed rss status < "$metrics_file"
  if (( rss > MAX_RSS_KB )); then
    printf 'Memory guard failed for %s: max RSS %s KB exceeded MAX_RSS_KB=%s KB\n' "$label" "$rss" "$MAX_RSS_KB" >&2
    exit 3
  fi
}

label_for_threads() {
  local label="threads_$1"
  printf '%s\n' "${label//[^A-Za-z0-9_]/_}"
}

run_case() {
  local threads="$1"
  local label
  label="$(label_for_threads "$threads")"
  printf 'Running Rust bbnorm-rs on %s/%s reads with tablereads=%s and threads=%s...\n' "$(basename "$DATA1")" "$READS" "$TABLE_READS" "$threads"
  set +e
  measure "$OUT/${label}.metrics.tsv" "$OUT/${label}.stdout.log" "$OUT/${label}.stderr.log" \
    target/release/bbnorm-rs \
    "${COMMON[@]}" \
    "${EXTRA_ARGS_ARRAY[@]}" \
    "threads=$threads" \
    "out=$OUT/${label}.keep1.fq" "out2=$OUT/${label}.keep2.fq" \
    "outlow=$OUT/${label}.low1.fq" "outlow2=$OUT/${label}.low2.fq" \
    "outmid=$OUT/${label}.mid1.fq" "outmid2=$OUT/${label}.mid2.fq" \
    "outhigh=$OUT/${label}.high1.fq" "outhigh2=$OUT/${label}.high2.fq" \
    "hist=$OUT/${label}.hist.tsv" "rhist=$OUT/${label}.rhist.tsv"
  local status=$?
  set -e
  write_result "$label" "$OUT/${label}.metrics.tsv"
  enforce_max_rss "$label" "$OUT/${label}.metrics.tsv"
  if (( status != 0 )); then
    printf 'Benchmark case %s failed with status %s; see %s and %s\n' "$label" "$status" "$OUT/${label}.stdout.log" "$OUT/${label}.stderr.log" >&2
    exit "$status"
  fi
}

compare_to_baseline() {
  local baseline="$1"
  local label="$2"
  for suffix in keep1.fq keep2.fq low1.fq low2.fq mid1.fq mid2.fq high1.fq high2.fq hist.tsv rhist.tsv; do
    cmp "$OUT/$baseline.$suffix" "$OUT/$label.$suffix"
  done
}

cargo build --release --quiet

printf 'dataset1\t%s\n' "$DATA1" > "$OUT/dataset.tsv"
printf 'dataset2\t%s\n' "$DATA2" >> "$OUT/dataset.tsv"
printf 'reads\t%s\n' "$READS" >> "$OUT/dataset.tsv"
printf 'tablereads\t%s\n' "$TABLE_READS" >> "$OUT/dataset.tsv"
printf 'keep_outputs\t%s\n' "$KEEP_OUTPUTS" >> "$OUT/dataset.tsv"
printf 'max_rss_kb_guard\t%s\n' "${MAX_RSS_KB:-none}" >> "$OUT/dataset.tsv"
printf 'timeout\t%s\n' "${TIMEOUT:-none}" >> "$OUT/dataset.tsv"
printf 'measure_script\t%s\n' "$MEASURE_SCRIPT" >> "$OUT/dataset.tsv"
printf 'extra_args\t%s\n' "${EXTRA_ARGS:-none}" >> "$OUT/dataset.tsv"
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
print("\nBiological dataset thread scaling summary:")
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

printf '\nBiological dataset benchmark passed. Per-run metrics: %s\n' "$OUT/results.tsv"
printf 'Dataset metadata: %s\n' "$OUT/dataset.tsv"
printf 'Outputs and logs: %s\n' "$OUT"

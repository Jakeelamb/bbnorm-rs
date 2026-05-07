#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT"

BIO_ROOT="${BIO_ROOT:-/home/jake/Projects/biological data}"
OUT="${1:-tmp/biological_matrix_benchmark}"
MATRIX_READS="${MATRIX_READS:-1000}"
MATRIX_THREAD_CASES="${MATRIX_THREAD_CASES:-1 2 auto}"
MATRIX_CASES="${MATRIX_CASES:-all}"
MATRIX_EXTRA_ARGS="${MATRIX_EXTRA_ARGS:-}"
KEEP_OUTPUTS="${KEEP_OUTPUTS:-0}"
MAX_RSS_KB="${MAX_RSS_KB:-}"
PYTHON="${PYTHON:-python}"

CASES=(
  "scer_s288c_srr23631023|$BIO_ROOT/reads/short_reads/scer_s288c_pe_srr23631023/SRR23631023_1.fastq.gz|$BIO_ROOT/reads/short_reads/scer_s288c_pe_srr23631023/SRR23631023_2.fastq.gz"
  "scer_s288c_err915337|$BIO_ROOT/reads/short_reads/scer_s288c_pe_err915337/ERR915337_1.fastq.gz|$BIO_ROOT/reads/short_reads/scer_s288c_pe_err915337/ERR915337_2.fastq.gz"
  "spombe_972_srr17530188|$BIO_ROOT/reads/short_reads/spombe_972_pe_srr17530188/SRR17530188_1.fastq.gz|$BIO_ROOT/reads/short_reads/spombe_972_pe_srr17530188/SRR17530188_2.fastq.gz"
  "ecoli_mg1655_drr023054|$BIO_ROOT/reads/short_reads/ecoli_mg1655_pe_drr023054/DRR023054_1.fastq.gz|$BIO_ROOT/reads/short_reads/ecoli_mg1655_pe_drr023054/DRR023054_2.fastq.gz"
  "ecoli_mg1655_drr217208|$BIO_ROOT/reads/short_reads/ecoli_mg1655_pe_drr217208/DRR217208_1.fastq.gz|$BIO_ROOT/reads/short_reads/ecoli_mg1655_pe_drr217208/DRR217208_2.fastq.gz"
  "ecoli_mg1655_srr13921545|$BIO_ROOT/reads/short_reads/ecoli_mg1655_pe_srr13921545/SRR13921545_1.fastq.gz|$BIO_ROOT/reads/short_reads/ecoli_mg1655_pe_srr13921545/SRR13921545_2.fastq.gz"
  "bacillus_subtilis_drr066522|$BIO_ROOT/reads/short_reads/bacillus_subtilis_168_pe_drr066522/DRR066522_1.fastq.gz|$BIO_ROOT/reads/short_reads/bacillus_subtilis_168_pe_drr066522/DRR066522_2.fastq.gz"
)

list_cases() {
  local entry name data1 data2
  printf 'case\tdataset1\tdataset2\n'
  for entry in "${CASES[@]}"; do
    IFS='|' read -r name data1 data2 <<< "$entry"
    printf '%s\t%s\t%s\n' "$name" "$data1" "$data2"
  done
}

if [[ "$MATRIX_CASES" == "list" ]]; then
  list_cases
  exit 0
fi

rm -rf "$OUT"
mkdir -p "$OUT"

case_selected() {
  local name="$1"
  local wanted
  for wanted in $MATRIX_CASES; do
    if [[ "$wanted" == "all" || "$wanted" == "$name" ]]; then
      return 0
    fi
  done
  return 1
}

printf 'case\tstatus\tfastest_thread\tfastest_seconds\tmax_rss_kb\toutput\n' > "$OUT/summary.tsv"

matched_cases=0
for entry in "${CASES[@]}"; do
  IFS='|' read -r name data1 data2 <<< "$entry"
  case_out="$OUT/$name"

  if ! case_selected "$name"; then
    continue
  fi
  matched_cases=$((matched_cases + 1))

  if [[ ! -f "$data1" || ! -f "$data2" ]]; then
    printf '%s\tskipped\tNA\tNA\tNA\t%s\n' "$name" "$case_out" >> "$OUT/summary.tsv"
    printf 'Skipping %s; missing paired dataset files.\n' "$name" >&2
    continue
  fi

  printf 'Running biological matrix case %s with reads=%s threads=%s extra_args=%s...\n' "$name" "$MATRIX_READS" "$MATRIX_THREAD_CASES" "${MATRIX_EXTRA_ARGS:-none}"
  READS="$MATRIX_READS" \
    TABLE_READS="$MATRIX_READS" \
    THREAD_CASES="$MATRIX_THREAD_CASES" \
    EXTRA_ARGS="$MATRIX_EXTRA_ARGS" \
    KEEP_OUTPUTS="$KEEP_OUTPUTS" \
    MAX_RSS_KB="$MAX_RSS_KB" \
    DATA1="$data1" \
    DATA2="$data2" \
    PYTHON="$PYTHON" \
    scripts/benchmark_biological_dataset.sh "$case_out"

  "$PYTHON" - "$case_out/results.tsv" "$name" "$case_out" >> "$OUT/summary.tsv" <<'PY'
import csv
import sys

path, name, output = sys.argv[1:4]
rows = list(csv.DictReader(open(path, encoding="utf-8"), delimiter="\t"))
fastest = min(rows, key=lambda row: float(row["elapsed_seconds"]))
max_rss = max(int(row["max_rss_kb"]) for row in rows)
print(
    f"{name}\tpassed\t{fastest['thread_case']}\t"
    f"{float(fastest['elapsed_seconds']):.6f}\t{max_rss}\t{output}"
)
PY
done

if [[ "$matched_cases" -eq 0 ]]; then
  printf 'No matrix cases matched MATRIX_CASES=%s\n' "$MATRIX_CASES" >&2
  printf 'Available cases:\n' >&2
  for entry in "${CASES[@]}"; do
    IFS='|' read -r name _ _ <<< "$entry"
    printf '  %s\n' "$name" >&2
  done
  exit 2
fi

printf '\nBiological matrix benchmark summary:\n'
column -t -s $'\t' "$OUT/summary.tsv" || cat "$OUT/summary.tsv"
printf '\nBiological matrix benchmark complete. Summary: %s\n' "$OUT/summary.tsv"

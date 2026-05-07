#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT"

OUT="${1:-tmp/biological_multipass_truncated_parity}"
BIO_ROOT="${BIO_ROOT:-/home/jake/Projects/biological data}"
DATA1="${DATA1:-$BIO_ROOT/reads/short_reads/scer_s288c_pe_srr23631023/SRR23631023_1.fastq.gz}"
DATA2="${DATA2:-$BIO_ROOT/reads/short_reads/scer_s288c_pe_srr23631023/SRR23631023_2.fastq.gz}"
READS="${READS:-1000}"
TABLE_READS="${TABLE_READS:-$READS}"

rm -rf "$OUT"
mkdir -p "$OUT"
printf 'passes\tstatus\toutput\n' > "$OUT/summary.tsv"

run_case() {
  local passes="$1"
  local case_out="$OUT/passes_$passes"

  printf '\nRunning truncated biological Java parity for passes=%s...\n' "$passes"
  BIO_ROOT="$BIO_ROOT" \
  DATA1="$DATA1" \
  DATA2="$DATA2" \
  READS="$READS" \
  TABLE_READS="$TABLE_READS" \
  EXTRA_ARGS="passes=$passes" \
  scripts/parity_biological_dataset_truncated.sh "$case_out"
  printf '%s\tpassed\t%s\n' "$passes" "$case_out" >> "$OUT/summary.tsv"
}

run_case 2
run_case 3
run_case 4

printf '\nTruncated biological multipass parity summary:\n'
column -t -s $'\t' "$OUT/summary.tsv" || cat "$OUT/summary.tsv"
printf '\nTruncated biological multipass parity passed. Outputs and logs: %s\n' "$OUT"

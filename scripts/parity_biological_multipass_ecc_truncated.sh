#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT"

OUT="${1:-tmp/biological_multipass_ecc_truncated_parity}"
BIO_ROOT="${BIO_ROOT:-/home/jake/Projects/biological data}"
DATA1="${DATA1:-$BIO_ROOT/reads/short_reads/scer_s288c_pe_srr23631023/SRR23631023_1.fastq.gz}"
DATA2="${DATA2:-$BIO_ROOT/reads/short_reads/scer_s288c_pe_srr23631023/SRR23631023_2.fastq.gz}"
READS="${READS:-1000}"
TABLE_READS="${TABLE_READS:-$READS}"

rm -rf "$OUT"
mkdir -p "$OUT"
printf 'mode\tstatus\toutput\n' > "$OUT/summary.tsv"

run_case() {
  local mode_name="$1"
  local extra_args="$2"
  local case_out="$OUT/$mode_name"

  printf '\nRunning truncated biological multipass ECC Java parity for %s...\n' "$mode_name"
  BIO_ROOT="$BIO_ROOT" \
  DATA1="$DATA1" \
  DATA2="$DATA2" \
  READS="$READS" \
  TABLE_READS="$TABLE_READS" \
  EXTRA_ARGS="$extra_args" \
  scripts/parity_biological_dataset_truncated.sh "$case_out"
  printf '%s\tpassed\t%s\n' "$mode_name" "$case_out" >> "$OUT/summary.tsv"
}

run_case passes2_ecc "passes=2 ecc=t keepall=t ecco=f target=999999999 max=999999999"
run_case passes2_ecc_markuncorrectable "passes=2 ecc=t keepall=t ecco=f target=999999999 max=999999999 markuncorrectableerrors=t"

printf '\nTruncated biological multipass ECC parity summary:\n'
column -t -s $'\t' "$OUT/summary.tsv" || cat "$OUT/summary.tsv"
printf '\nTruncated biological multipass ECC parity passed. Outputs and logs: %s\n' "$OUT"

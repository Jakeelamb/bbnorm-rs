#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT"

OUT="${1:-tmp/working_modes_benchmark_smoke}"
BIO_ROOT="${BIO_ROOT:-/home/jake/Projects/biological data}"
DATA1="${DATA1:-$BIO_ROOT/reads/short_reads/scer_s288c_pe_srr23631023/SRR23631023_1.fastq.gz}"
DATA2="${DATA2:-$BIO_ROOT/reads/short_reads/scer_s288c_pe_srr23631023/SRR23631023_2.fastq.gz}"
THREAD_CASES="${THREAD_CASES:-1 2 auto}"
READS="${READS:-1000}"
TABLE_READS="${TABLE_READS:-$READS}"
KEEP_OUTPUTS="${KEEP_OUTPUTS:-0}"
MAX_RSS_KB="${MAX_RSS_KB:-1000000}"

rm -rf "$OUT"
mkdir -p "$OUT"

printf 'mode\tstatus\toutput\n' > "$OUT/summary.tsv"

run_mode() {
  local mode_name="$1"
  local extra_args="$2"
  local mode_out="$OUT/$mode_name"

  printf '\nRunning working-mode benchmark smoke for %s...\n' "$mode_name"
  DATA1="$DATA1" \
    DATA2="$DATA2" \
    THREAD_CASES="$THREAD_CASES" \
    READS="$READS" \
    TABLE_READS="$TABLE_READS" \
    KEEP_OUTPUTS="$KEEP_OUTPUTS" \
    MAX_RSS_KB="$MAX_RSS_KB" \
    EXTRA_ARGS="$extra_args" \
    scripts/benchmark_biological_dataset.sh "$mode_out"
  printf '%s\tpassed\t%s\n' "$mode_name" "$mode_out" >> "$OUT/summary.tsv"
}

run_mode default ""
run_mode countup "countup=t"
run_mode countup_ecc "countup=t ecc=t markuncorrectableerrors=t"
run_mode long_kmer "k=40"
run_mode long_kmer_fixspikes "k=40 fixspikes=t"
run_mode multipass "passes=2"
run_mode multipass_ecc "passes=2 ecc=t markuncorrectableerrors=t"

printf '\nWorking-mode benchmark smoke summary:\n'
column -t -s $'\t' "$OUT/summary.tsv" || cat "$OUT/summary.tsv"
printf '\nWorking-mode benchmark smoke passed. Outputs and logs: %s\n' "$OUT"

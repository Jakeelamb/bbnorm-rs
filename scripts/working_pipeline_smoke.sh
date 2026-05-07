#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT"

OUT="${1:-tmp/working_pipeline_smoke}"
BIO_ROOT="${BIO_ROOT:-/home/jake/Projects/biological data}"
ALT_BIO_DATA1="${ALT_BIO_DATA1:-$BIO_ROOT/reads/short_reads/ecoli_mg1655_pe_drr023054/DRR023054_1.fastq.gz}"
ALT_BIO_DATA2="${ALT_BIO_DATA2:-$BIO_ROOT/reads/short_reads/ecoli_mg1655_pe_drr023054/DRR023054_2.fastq.gz}"
ALT_BIO_ERROR_DATA1="${ALT_BIO_ERROR_DATA1:-$BIO_ROOT/reads/short_reads/spombe_972_pe_srr17530188/SRR17530188_1.fastq.gz}"
ALT_BIO_ERROR_DATA2="${ALT_BIO_ERROR_DATA2:-$BIO_ROOT/reads/short_reads/spombe_972_pe_srr17530188/SRR17530188_2.fastq.gz}"
THREAD_CASES="${THREAD_CASES:-1 2 auto}"
MATRIX_READS="${MATRIX_READS:-1000}"
MATRIX_CASES="${MATRIX_CASES:-all}"
MATRIX_THREAD_CASES="${MATRIX_THREAD_CASES:-$THREAD_CASES}"
MATRIX_EXTRA_ARGS="${MATRIX_EXTRA_ARGS:-}"
MAX_RSS_KB="${MAX_RSS_KB:-}"
COUNTUP_STRESS_READS="${COUNTUP_STRESS_READS:-1000}"
COUNTUP_STRESS_THREADS="${COUNTUP_STRESS_THREADS:-$THREAD_CASES}"
SIDE_STRESS_READS="${SIDE_STRESS_READS:-1000}"
SIDE_STRESS_THREADS="${SIDE_STRESS_THREADS:-$THREAD_CASES}"
SIDE_OUTPUT_STATS_STRESS_READS="${SIDE_OUTPUT_STATS_STRESS_READS:-1000}"
SIDE_OUTPUT_STATS_STRESS_THREADS="${SIDE_OUTPUT_STATS_STRESS_THREADS:-$THREAD_CASES}"
LONG_KMER_STRESS_READS="${LONG_KMER_STRESS_READS:-1000}"
LONG_KMER_STRESS_THREADS="${LONG_KMER_STRESS_THREADS:-$THREAD_CASES}"
OVERLAP_STRESS_THREADS="${OVERLAP_STRESS_THREADS:-$THREAD_CASES}"
WORKING_MODES_BENCHMARK_READS="${WORKING_MODES_BENCHMARK_READS:-1000}"
WORKING_MODES_BENCHMARK_THREADS="${WORKING_MODES_BENCHMARK_THREADS:-$THREAD_CASES}"
BIOLOGICAL_JAVA_PARITY_READS="${BIOLOGICAL_JAVA_PARITY_READS:-1000}"
BIOLOGICAL_MULTIPASS_JAVA_PARITY_READS="${BIOLOGICAL_MULTIPASS_JAVA_PARITY_READS:-1000}"
BIOLOGICAL_MULTIPASS_ECC_JAVA_PARITY_READS="${BIOLOGICAL_MULTIPASS_ECC_JAVA_PARITY_READS:-1000}"
BIOLOGICAL_COUNTUP_JAVA_GUARD_READS="${BIOLOGICAL_COUNTUP_JAVA_GUARD_READS:-1000}"
BIOLOGICAL_ALT_JAVA_GUARD_READS="${BIOLOGICAL_ALT_JAVA_GUARD_READS:-1000}"
BIOLOGICAL_ALT_ERROR_JAVA_GUARD_READS="${BIOLOGICAL_ALT_ERROR_JAVA_GUARD_READS:-1000}"
BIOLOGICAL_JAVA_WORKING_MODES_READS="${BIOLOGICAL_JAVA_WORKING_MODES_READS:-1000}"

rm -rf "$OUT"
mkdir -p "$OUT"
printf 'stage\tstatus\toutput\n' > "$OUT/summary.tsv"
{
  printf 'key\tvalue\n'
  printf 'bio_root\t%s\n' "$BIO_ROOT"
  printf 'alt_bio_data1\t%s\n' "$ALT_BIO_DATA1"
  printf 'alt_bio_data2\t%s\n' "$ALT_BIO_DATA2"
  printf 'alt_bio_error_data1\t%s\n' "$ALT_BIO_ERROR_DATA1"
  printf 'alt_bio_error_data2\t%s\n' "$ALT_BIO_ERROR_DATA2"
  printf 'thread_cases\t%s\n' "$THREAD_CASES"
  printf 'matrix_reads\t%s\n' "$MATRIX_READS"
  printf 'matrix_cases\t%s\n' "$MATRIX_CASES"
  printf 'matrix_thread_cases\t%s\n' "$MATRIX_THREAD_CASES"
  printf 'matrix_extra_args\t%s\n' "${MATRIX_EXTRA_ARGS:-none}"
  printf 'countup_stress_reads\t%s\n' "$COUNTUP_STRESS_READS"
  printf 'countup_stress_threads\t%s\n' "$COUNTUP_STRESS_THREADS"
  printf 'side_stress_reads\t%s\n' "$SIDE_STRESS_READS"
  printf 'side_stress_threads\t%s\n' "$SIDE_STRESS_THREADS"
  printf 'side_output_stats_stress_reads\t%s\n' "$SIDE_OUTPUT_STATS_STRESS_READS"
  printf 'side_output_stats_stress_threads\t%s\n' "$SIDE_OUTPUT_STATS_STRESS_THREADS"
  printf 'long_kmer_stress_reads\t%s\n' "$LONG_KMER_STRESS_READS"
  printf 'long_kmer_stress_threads\t%s\n' "$LONG_KMER_STRESS_THREADS"
  printf 'overlap_stress_threads\t%s\n' "$OVERLAP_STRESS_THREADS"
  printf 'working_modes_benchmark_reads\t%s\n' "$WORKING_MODES_BENCHMARK_READS"
  printf 'working_modes_benchmark_threads\t%s\n' "$WORKING_MODES_BENCHMARK_THREADS"
  printf 'biological_java_parity_reads\t%s\n' "$BIOLOGICAL_JAVA_PARITY_READS"
  printf 'biological_multipass_java_parity_reads\t%s\n' "$BIOLOGICAL_MULTIPASS_JAVA_PARITY_READS"
  printf 'biological_multipass_ecc_java_parity_reads\t%s\n' "$BIOLOGICAL_MULTIPASS_ECC_JAVA_PARITY_READS"
  printf 'biological_countup_java_guard_reads\t%s\n' "$BIOLOGICAL_COUNTUP_JAVA_GUARD_READS"
  printf 'biological_alt_java_guard_reads\t%s\n' "$BIOLOGICAL_ALT_JAVA_GUARD_READS"
  printf 'biological_alt_error_java_guard_reads\t%s\n' "$BIOLOGICAL_ALT_ERROR_JAVA_GUARD_READS"
  printf 'biological_java_working_modes_reads\t%s\n' "$BIOLOGICAL_JAVA_WORKING_MODES_READS"
  printf 'max_rss_kb\t%s\n' "${MAX_RSS_KB:-none}"
} > "$OUT/config.tsv"

printf 'Running component smoke...\n'
scripts/component_smoke.sh "$OUT/component"
printf 'component_smoke\tpassed\t%s\n' "$OUT/component" >> "$OUT/summary.tsv"

printf '\nRunning fallback behavior smoke...\n'
scripts/fallback_smoke.sh "$OUT/fallback"
printf 'fallback_smoke\tpassed\t%s\n' "$OUT/fallback" >> "$OUT/summary.tsv"

printf '\nRunning bundled phiX thread-scaling smoke...\n'
THREAD_CASES="$THREAD_CASES" scripts/benchmark_thread_scaling.sh "$OUT/thread_scaling"
printf 'thread_scaling\tpassed\t%s\n' "$OUT/thread_scaling" >> "$OUT/summary.tsv"

default_bio_1="$BIO_ROOT/reads/short_reads/scer_s288c_pe_srr23631023/SRR23631023_1.fastq.gz"
default_bio_2="$BIO_ROOT/reads/short_reads/scer_s288c_pe_srr23631023/SRR23631023_2.fastq.gz"
if [[ -f "$default_bio_1" && -f "$default_bio_2" ]]; then
  printf '\nRunning count-up biological stress smoke...\n'
  BIO_ROOT="$BIO_ROOT" \
    READS="$COUNTUP_STRESS_READS" \
    TABLE_READS="$COUNTUP_STRESS_READS" \
    THREAD_CASES="$COUNTUP_STRESS_THREADS" \
    scripts/parity_countup_biological_stress.sh "$OUT/countup_biological_stress"
  printf 'countup_biological_stress\tpassed\t%s\n' "$OUT/countup_biological_stress" >> "$OUT/summary.tsv"

  printf '\nRunning side-routing biological stress smoke...\n'
  BIO_ROOT="$BIO_ROOT" \
    READS="$SIDE_STRESS_READS" \
    TABLE_READS="$SIDE_STRESS_READS" \
    THREAD_CASES="$SIDE_STRESS_THREADS" \
    scripts/parity_side_routing_biological_stress.sh "$OUT/side_routing_biological_stress"
  printf 'side_routing_biological_stress\tpassed\t%s\n' "$OUT/side_routing_biological_stress" >> "$OUT/summary.tsv"

  printf '\nRunning side-output stats biological stress smoke...\n'
  BIO_ROOT="$BIO_ROOT" \
    READS="$SIDE_OUTPUT_STATS_STRESS_READS" \
    TABLE_READS="$SIDE_OUTPUT_STATS_STRESS_READS" \
    THREAD_CASES="$SIDE_OUTPUT_STATS_STRESS_THREADS" \
    scripts/parity_side_output_stats_biological_stress.sh "$OUT/side_output_stats_biological_stress"
  printf 'side_output_stats_biological_stress\tpassed\t%s\n' "$OUT/side_output_stats_biological_stress" >> "$OUT/summary.tsv"

  printf '\nRunning long-kmer biological stress smoke...\n'
  BIO_ROOT="$BIO_ROOT" \
    READS="$LONG_KMER_STRESS_READS" \
    TABLE_READS="$LONG_KMER_STRESS_READS" \
    THREAD_CASES="$LONG_KMER_STRESS_THREADS" \
    scripts/parity_long_kmer_biological_stress.sh "$OUT/long_kmer_biological_stress"
  printf 'long_kmer_biological_stress\tpassed\t%s\n' "$OUT/long_kmer_biological_stress" >> "$OUT/summary.tsv"

  printf '\nRunning overlap-ECC biological stress smoke...\n'
  BIO_ROOT="$BIO_ROOT" \
    THREAD_CASES="$OVERLAP_STRESS_THREADS" \
    scripts/parity_overlap_ecc_biological_stress.sh "$OUT/overlap_ecc_biological_stress"
  printf 'overlap_ecc_biological_stress\tpassed\t%s\n' "$OUT/overlap_ecc_biological_stress" >> "$OUT/summary.tsv"

  printf '\nRunning truncated biological Java parity smoke...\n'
  BIO_ROOT="$BIO_ROOT" \
    READS="$BIOLOGICAL_JAVA_PARITY_READS" \
    TABLE_READS="$BIOLOGICAL_JAVA_PARITY_READS" \
  scripts/parity_biological_dataset_truncated.sh "$OUT/biological_java_parity"
  printf 'biological_java_parity\tpassed\t%s\n' "$OUT/biological_java_parity" >> "$OUT/summary.tsv"

  printf '\nRunning truncated biological multipass Java parity smoke...\n'
  BIO_ROOT="$BIO_ROOT" \
    READS="$BIOLOGICAL_MULTIPASS_JAVA_PARITY_READS" \
    TABLE_READS="$BIOLOGICAL_MULTIPASS_JAVA_PARITY_READS" \
  scripts/parity_biological_multipass_truncated.sh "$OUT/biological_multipass_java_parity"
  printf 'biological_multipass_java_parity\tpassed\t%s\n' "$OUT/biological_multipass_java_parity" >> "$OUT/summary.tsv"

  printf '\nRunning truncated biological multipass ECC Java parity smoke...\n'
  BIO_ROOT="$BIO_ROOT" \
    READS="$BIOLOGICAL_MULTIPASS_ECC_JAVA_PARITY_READS" \
    TABLE_READS="$BIOLOGICAL_MULTIPASS_ECC_JAVA_PARITY_READS" \
    scripts/parity_biological_multipass_ecc_truncated.sh "$OUT/biological_multipass_ecc_java_parity"
  printf 'biological_multipass_ecc_java_parity\tpassed\t%s\n' "$OUT/biological_multipass_ecc_java_parity" >> "$OUT/summary.tsv"

  printf '\nRunning truncated biological countup Java guard...\n'
  BIO_ROOT="$BIO_ROOT" \
    READS="$BIOLOGICAL_COUNTUP_JAVA_GUARD_READS" \
    TABLE_READS="$BIOLOGICAL_COUNTUP_JAVA_GUARD_READS" \
    scripts/parity_biological_countup_java_guard.sh "$OUT/biological_countup_java_guard"
  printf 'biological_countup_java_guard\tpassed\t%s\n' "$OUT/biological_countup_java_guard" >> "$OUT/summary.tsv"

  if [[ -f "$ALT_BIO_DATA1" && -f "$ALT_BIO_DATA2" ]]; then
    printf '\nRunning alternate-dataset biological Java guard...\n'
    DATA1="$ALT_BIO_DATA1" \
      DATA2="$ALT_BIO_DATA2" \
      READS="$BIOLOGICAL_ALT_JAVA_GUARD_READS" \
      TABLE_READS="$BIOLOGICAL_ALT_JAVA_GUARD_READS" \
      scripts/parity_biological_pairstreamer_java_guard.sh "$OUT/biological_alt_java_guard"
    printf 'biological_alt_java_guard\tpassed\t%s\n' "$OUT/biological_alt_java_guard" >> "$OUT/summary.tsv"
  else
    printf '\nSkipping alternate-dataset biological Java guard; ALT_BIO_DATA1=%s ALT_BIO_DATA2=%s\n' "$ALT_BIO_DATA1" "$ALT_BIO_DATA2"
    printf 'biological_alt_java_guard\tskipped\tALT_BIO_DATA1=%s ALT_BIO_DATA2=%s\n' "$ALT_BIO_DATA1" "$ALT_BIO_DATA2" >> "$OUT/summary.tsv"
  fi

  if [[ -f "$ALT_BIO_ERROR_DATA1" && -f "$ALT_BIO_ERROR_DATA2" ]]; then
    printf '\nRunning alternate-dataset biological Java error-state guard...\n'
    DATA1="$ALT_BIO_ERROR_DATA1" \
      DATA2="$ALT_BIO_ERROR_DATA2" \
      READS="$BIOLOGICAL_ALT_ERROR_JAVA_GUARD_READS" \
      TABLE_READS="$BIOLOGICAL_ALT_ERROR_JAVA_GUARD_READS" \
      scripts/parity_biological_errorstate_java_guard.sh "$OUT/biological_alt_error_java_guard"
    printf 'biological_alt_error_java_guard\tpassed\t%s\n' "$OUT/biological_alt_error_java_guard" >> "$OUT/summary.tsv"
  else
    printf '\nSkipping alternate-dataset biological Java error-state guard; ALT_BIO_ERROR_DATA1=%s ALT_BIO_ERROR_DATA2=%s\n' "$ALT_BIO_ERROR_DATA1" "$ALT_BIO_ERROR_DATA2"
    printf 'biological_alt_error_java_guard\tskipped\tALT_BIO_ERROR_DATA1=%s ALT_BIO_ERROR_DATA2=%s\n' "$ALT_BIO_ERROR_DATA1" "$ALT_BIO_ERROR_DATA2" >> "$OUT/summary.tsv"
  fi

  printf '\nRunning truncated biological working-mode Java parity probe...\n'
  BIO_ROOT="$BIO_ROOT" \
    READS="$BIOLOGICAL_JAVA_WORKING_MODES_READS" \
    TABLE_READS="$BIOLOGICAL_JAVA_WORKING_MODES_READS" \
    scripts/parity_biological_working_modes_truncated.sh "$OUT/biological_java_working_modes"
  printf 'biological_java_working_modes_probe\tpassed\t%s\n' "$OUT/biological_java_working_modes" >> "$OUT/summary.tsv"

  printf '\nRunning working-mode benchmark smoke...\n'
  BIO_ROOT="$BIO_ROOT" \
    THREAD_CASES="$WORKING_MODES_BENCHMARK_THREADS" \
    READS="$WORKING_MODES_BENCHMARK_READS" \
    TABLE_READS="$WORKING_MODES_BENCHMARK_READS" \
    MAX_RSS_KB="${MAX_RSS_KB:-1000000}" \
    scripts/benchmark_working_modes_smoke.sh "$OUT/working_modes_benchmark"
  printf 'working_modes_benchmark\tpassed\t%s\n' "$OUT/working_modes_benchmark" >> "$OUT/summary.tsv"

  printf '\nRunning biological RSS guard smoke...\n'
  BIO_ROOT="$BIO_ROOT" scripts/benchmark_biological_guard_smoke.sh "$OUT/biological_guard"
  printf 'biological_guard\tpassed\t%s\n' "$OUT/biological_guard" >> "$OUT/summary.tsv"

  printf '\nRunning biological matrix smoke...\n'
  BIO_ROOT="$BIO_ROOT" \
    MATRIX_READS="$MATRIX_READS" \
    MATRIX_CASES="$MATRIX_CASES" \
    MATRIX_THREAD_CASES="$MATRIX_THREAD_CASES" \
    MATRIX_EXTRA_ARGS="$MATRIX_EXTRA_ARGS" \
    MAX_RSS_KB="$MAX_RSS_KB" \
    scripts/benchmark_biological_matrix.sh "$OUT/biological_matrix"
  printf 'biological_matrix\tpassed\t%s\n' "$OUT/biological_matrix" >> "$OUT/summary.tsv"
else
  printf '\nSkipping biological smoke checks; default paired dataset files are missing under BIO_ROOT=%s\n' "$BIO_ROOT"
  printf 'countup_biological_stress\tskipped\tBIO_ROOT=%s\n' "$BIO_ROOT" >> "$OUT/summary.tsv"
  printf 'side_routing_biological_stress\tskipped\tBIO_ROOT=%s\n' "$BIO_ROOT" >> "$OUT/summary.tsv"
  printf 'side_output_stats_biological_stress\tskipped\tBIO_ROOT=%s\n' "$BIO_ROOT" >> "$OUT/summary.tsv"
  printf 'long_kmer_biological_stress\tskipped\tBIO_ROOT=%s\n' "$BIO_ROOT" >> "$OUT/summary.tsv"
  printf 'overlap_ecc_biological_stress\tskipped\tBIO_ROOT=%s\n' "$BIO_ROOT" >> "$OUT/summary.tsv"
  printf 'biological_java_parity\tskipped\tBIO_ROOT=%s\n' "$BIO_ROOT" >> "$OUT/summary.tsv"
  printf 'biological_multipass_java_parity\tskipped\tBIO_ROOT=%s\n' "$BIO_ROOT" >> "$OUT/summary.tsv"
  printf 'biological_multipass_ecc_java_parity\tskipped\tBIO_ROOT=%s\n' "$BIO_ROOT" >> "$OUT/summary.tsv"
  printf 'biological_countup_java_guard\tskipped\tBIO_ROOT=%s\n' "$BIO_ROOT" >> "$OUT/summary.tsv"
  printf 'biological_alt_java_guard\tskipped\tBIO_ROOT=%s\n' "$BIO_ROOT" >> "$OUT/summary.tsv"
  printf 'biological_alt_error_java_guard\tskipped\tBIO_ROOT=%s\n' "$BIO_ROOT" >> "$OUT/summary.tsv"
  printf 'biological_java_working_modes_probe\tskipped\tBIO_ROOT=%s\n' "$BIO_ROOT" >> "$OUT/summary.tsv"
  printf 'working_modes_benchmark\tskipped\tBIO_ROOT=%s\n' "$BIO_ROOT" >> "$OUT/summary.tsv"
  printf 'biological_guard\tskipped\tBIO_ROOT=%s\n' "$BIO_ROOT" >> "$OUT/summary.tsv"
  printf 'biological_matrix\tskipped\tBIO_ROOT=%s\n' "$BIO_ROOT" >> "$OUT/summary.tsv"
fi

printf '\nWorking pipeline smoke summary:\n'
column -t -s $'\t' "$OUT/summary.tsv" || cat "$OUT/summary.tsv"
printf '\nWorking pipeline smoke config:\n'
column -t -s $'\t' "$OUT/config.tsv" || cat "$OUT/config.tsv"
printf '\nWorking pipeline smoke passed. Outputs and logs: %s\n' "$OUT"

#!/usr/bin/env bash
set -euo pipefail

if [[ $# -lt 1 ]]; then
  cat >&2 <<'USAGE'
Usage: scripts/benchmark_giant_safe.sh <reads_R1.fq[.gz]> [reads_R2.fq[.gz]] [outdir]

Runs bbnorm-rs in a giant-dataset-safe mode: release binary, bounded count-min,
8 threads by default, null sequence outputs by default, and peak-RSS/time logs.
Override with environment variables such as THREADS=16, MEM=4g, TARGET=40,
MAX=80, K=31, READS=1000000, TABLE_READS=1000000, WRITE_OUTPUTS=1,
TIMEOUT=4h, EXTRA_ARGS='prefilterfraction=0.35'. When READS is set and
TABLE_READS is not, table building is capped to the same read limit so smoke
runs stay genuinely bounded. If pigz/unpigz are installed in ~/.local/bin, this
script adds that directory to PATH. Optional count-up spill guards are available
as MAX_COUNTUP_SPILL_BYTES_WRITTEN, MAX_COUNTUP_SPILL_PEAK_LIVE_BYTES,
MAX_COUNTUP_SPILL_INITIAL_RUNS, MAX_COUNTUP_SPILL_MERGE_RUNS,
MAX_COUNTUP_SPILL_FINAL_RUNS, and MAX_COUNTUP_SPILL_FINAL_LIVE_BYTES.
USAGE
  exit 2
fi

R1=$1
R2=${2:-}
OUT=${3:-tmp/giant_safe_benchmark_$(date +%Y%m%d_%H%M%S)}
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT"
export PATH="$HOME/.local/bin:$PATH"
THREADS=${THREADS:-8}
MEM=${MEM:-2g}
TARGET=${TARGET:-40}
MAX=${MAX:-80}
K=${K:-31}
MIN=${MIN:-5}
BITS=${BITS:-32}
READS=${READS:-}
TABLE_READS=${TABLE_READS:-${READS:-}}
TIMEOUT=${TIMEOUT:-}
WRITE_OUTPUTS=${WRITE_OUTPUTS:-0}
ZIPTHREADS=${ZIPTHREADS:-$THREADS}
EXTRA_ARGS=${EXTRA_ARGS:-}
MAX_RSS_KB=${MAX_RSS_KB:-}
PYTHON=${PYTHON:-python3}
MEASURE_SCRIPT=${MEASURE_SCRIPT:-scripts/measure_command.py}
UNIQUE_SUMMARY_SCRIPT=${UNIQUE_SUMMARY_SCRIPT:-scripts/extract_unique_kmer_summary.py}
MAX_COUNTUP_SPILL_INITIAL_RUNS=${MAX_COUNTUP_SPILL_INITIAL_RUNS:-}
MAX_COUNTUP_SPILL_MERGE_RUNS=${MAX_COUNTUP_SPILL_MERGE_RUNS:-}
MAX_COUNTUP_SPILL_FINAL_RUNS=${MAX_COUNTUP_SPILL_FINAL_RUNS:-}
MAX_COUNTUP_SPILL_BYTES_WRITTEN=${MAX_COUNTUP_SPILL_BYTES_WRITTEN:-}
MAX_COUNTUP_SPILL_PEAK_LIVE_BYTES=${MAX_COUNTUP_SPILL_PEAK_LIVE_BYTES:-}
MAX_COUNTUP_SPILL_FINAL_LIVE_BYTES=${MAX_COUNTUP_SPILL_FINAL_LIVE_BYTES:-}

mkdir -p "$OUT"

validate_optional_integer() {
  local name="$1"
  local value="${!name:-}"
  if [[ -n "$value" && ! "$value" =~ ^[0-9]+$ ]]; then
    printf '%s must be an integer count when set, got: %s\n' "$name" "$value" >&2
    exit 2
  fi
}

case "${MAX_RSS_KB,,}" in
  ""|"0"|"none"|"off"|"unlimited")
    MAX_RSS_KB=
    ;;
esac
if [[ -n "$MAX_RSS_KB" && ! "$MAX_RSS_KB" =~ ^[0-9]+$ ]]; then
  echo "MAX_RSS_KB must be an integer number of kilobytes, or 0/off/none/unlimited, got: $MAX_RSS_KB" >&2
  exit 2
fi
validate_optional_integer MAX_COUNTUP_SPILL_INITIAL_RUNS
validate_optional_integer MAX_COUNTUP_SPILL_MERGE_RUNS
validate_optional_integer MAX_COUNTUP_SPILL_FINAL_RUNS
validate_optional_integer MAX_COUNTUP_SPILL_BYTES_WRITTEN
validate_optional_integer MAX_COUNTUP_SPILL_PEAK_LIVE_BYTES
validate_optional_integer MAX_COUNTUP_SPILL_FINAL_LIVE_BYTES

if [[ ! -r "$R1" ]]; then
  echo "Missing readable R1 input: $R1" >&2
  exit 2
fi
if [[ -n "$R2" && ! -r "$R2" ]]; then
  echo "Missing readable R2 input: $R2" >&2
  exit 2
fi

cargo build --release

ARGS=(
  "in=$R1"
  "passes=1"
  "target=$TARGET"
  "max=$MAX"
  "min=$MIN"
  "k=$K"
  "threads=$THREADS"
  "zipthreads=$ZIPTHREADS"
  "mem=$MEM"
  "bits=$BITS"
  "autocountmin=t"
  "overwrite=t"
  "hist=$OUT/input.hist.tsv"
  "rhist=$OUT/input.rhist.tsv"
)

if [[ -n "$R2" ]]; then
  ARGS+=("in2=$R2")
fi
if [[ -n "$READS" ]]; then
  ARGS+=("reads=$READS")
fi
if [[ -n "$TABLE_READS" ]]; then
  ARGS+=("tablereads=$TABLE_READS")
fi
if [[ -n "$EXTRA_ARGS" ]]; then
  read -r -a EXTRA_ARG_ARRAY <<< "$EXTRA_ARGS"
  ARGS+=("${EXTRA_ARG_ARRAY[@]}")
fi
if [[ -n "$MAX_COUNTUP_SPILL_INITIAL_RUNS" ]]; then
  ARGS+=("maxcountupspillinitialruns=$MAX_COUNTUP_SPILL_INITIAL_RUNS")
fi
if [[ -n "$MAX_COUNTUP_SPILL_MERGE_RUNS" ]]; then
  ARGS+=("maxcountupspillmergeruns=$MAX_COUNTUP_SPILL_MERGE_RUNS")
fi
if [[ -n "$MAX_COUNTUP_SPILL_FINAL_RUNS" ]]; then
  ARGS+=("maxcountupspillfinalruns=$MAX_COUNTUP_SPILL_FINAL_RUNS")
fi
if [[ -n "$MAX_COUNTUP_SPILL_PEAK_LIVE_BYTES" ]]; then
  ARGS+=("maxcountupspillbytes=$MAX_COUNTUP_SPILL_PEAK_LIVE_BYTES")
fi
if [[ -n "$MAX_COUNTUP_SPILL_FINAL_LIVE_BYTES" ]]; then
  ARGS+=("maxcountupspillfinallivebytes=$MAX_COUNTUP_SPILL_FINAL_LIVE_BYTES")
fi
if [[ -n "$MAX_COUNTUP_SPILL_BYTES_WRITTEN" ]]; then
  ARGS+=("maxcountupspillwritebytes=$MAX_COUNTUP_SPILL_BYTES_WRITTEN")
fi

if [[ "$WRITE_OUTPUTS" == "1" ]]; then
  ARGS+=("out=$OUT/keep1.fq.gz" "outt=$OUT/toss1.fq.gz")
  if [[ -n "$R2" ]]; then
    ARGS+=("out2=$OUT/keep2.fq.gz" "outt2=$OUT/toss2.fq.gz")
  fi
else
  ARGS+=("out=null" "outt=null")
  if [[ -n "$R2" ]]; then
    ARGS+=("out2=null" "outt2=null")
  fi
fi

{
  printf 'timestamp\t%s\n' "$(date -Is)"
  printf 'host\t%s\n' "$(hostname)"
  printf 'threads\t%s\n' "$THREADS"
  printf 'zipthreads\t%s\n' "$ZIPTHREADS"
  printf 'mem\t%s\n' "$MEM"
  printf 'target\t%s\n' "$TARGET"
  printf 'max\t%s\n' "$MAX"
  printf 'k\t%s\n' "$K"
  printf 'bits\t%s\n' "$BITS"
  printf 'reads_limit\t%s\n' "${READS:-all}"
  printf 'table_reads_limit\t%s\n' "${TABLE_READS:-all}"
  printf 'write_outputs\t%s\n' "$WRITE_OUTPUTS"
  printf 'extra_args\t%s\n' "${EXTRA_ARGS:-}"
  printf 'max_rss_kb_guard\t%s\n' "${MAX_RSS_KB:-}"
  printf 'pigz\t%s\n' "$(command -v pigz || true)"
  printf 'unpigz\t%s\n' "$(command -v unpigz || true)"
  printf 'measure_script\t%s\n' "$MEASURE_SCRIPT"
  printf 'unique_summary_script\t%s\n' "$UNIQUE_SUMMARY_SCRIPT"
  printf 'max_countup_spill_initial_runs\t%s\n' "$MAX_COUNTUP_SPILL_INITIAL_RUNS"
  printf 'max_countup_spill_merge_runs\t%s\n' "$MAX_COUNTUP_SPILL_MERGE_RUNS"
  printf 'max_countup_spill_final_runs\t%s\n' "$MAX_COUNTUP_SPILL_FINAL_RUNS"
  printf 'max_countup_spill_bytes_written\t%s\n' "$MAX_COUNTUP_SPILL_BYTES_WRITTEN"
  printf 'max_countup_spill_peak_live_bytes\t%s\n' "$MAX_COUNTUP_SPILL_PEAK_LIVE_BYTES"
  printf 'max_countup_spill_final_live_bytes\t%s\n' "$MAX_COUNTUP_SPILL_FINAL_LIVE_BYTES"
  printf 'r1\t%s\n' "$R1"
  printf 'r2\t%s\n' "${R2:-}"
  printf 'r1_bytes\t%s\n' "$(stat -c %s "$R1")"
  if [[ -n "$R2" ]]; then printf 'r2_bytes\t%s\n' "$(stat -c %s "$R2")"; fi
  printf 'free_before\t%s\n' "$(free -h | tr '\n' ' ' | sed 's/[[:space:]]\+/ /g')"
  printf 'df_before\t%s\n' "$(df -h . | tr '\n' ' ' | sed 's/[[:space:]]\+/ /g')"
  printf 'command\t'
  printf '%q ' target/release/bbnorm-rs "${ARGS[@]}"
  printf '\n'
} > "$OUT/environment.tsv"

CMD=(target/release/bbnorm-rs "${ARGS[@]}")

set +e
if [[ -x "$MEASURE_SCRIPT" || -f "$MEASURE_SCRIPT" ]]; then
  MEASURE_ARGS=(
    "$PYTHON"
    "$MEASURE_SCRIPT"
    "--metrics" "$OUT/metrics.tsv"
    "--stdout" "$OUT/stdout.log"
    "--stderr" "$OUT/stderr.time.log"
  )
  if [[ -n "$TIMEOUT" ]]; then
    MEASURE_ARGS+=("--timeout" "$TIMEOUT")
  fi
  if [[ -n "$MAX_RSS_KB" ]]; then
    MEASURE_ARGS+=("--max-rss-kb" "$MAX_RSS_KB")
  fi
  MEASURE_ARGS+=("--" "${CMD[@]}")
  "${MEASURE_ARGS[@]}"
  STATUS=$?
elif [[ -x /usr/bin/time ]]; then
  /usr/bin/time -v "${CMD[@]}" > "$OUT/stdout.log" 2> "$OUT/stderr.time.log"
  STATUS=$?
else
  START_SECONDS=$(date +%s)
  "${CMD[@]}" > "$OUT/stdout.log" 2> "$OUT/stderr.time.log"
  STATUS=$?
  END_SECONDS=$(date +%s)
  {
    printf 'Elapsed (wall clock) time (seconds): %s\n' "$((END_SECONDS - START_SECONDS))"
    printf 'Maximum resident set size (kbytes): unavailable\n'
    printf 'Exit status: %s\n' "$STATUS"
  } >> "$OUT/stderr.time.log"
fi
set -e

{
  printf 'free_after\t%s\n' "$(free -h | tr '\n' ' ' | sed 's/[[:space:]]\+/ /g')"
  printf 'df_after\t%s\n' "$(df -h . | tr '\n' ' ' | sed 's/[[:space:]]\+/ /g')"
} >> "$OUT/environment.tsv"

awk '
  /Elapsed \(wall clock\) time/ {print "elapsed\t" $0}
  /Maximum resident set size/ {print "max_rss_kb\t" $6}
  /File system outputs/ {print "fs_outputs\t" $4}
  /Exit status/ {print "exit_status\t" $3}
  /RSS guard exceeded/ {print "rss_guard_exceeded\t" $4}
' "$OUT/stderr.time.log" > "$OUT/time_summary.tsv"

if [[ -s "$OUT/metrics.tsv" ]]; then
  {
    awk 'BEGIN {print "metric\tvalue"} {print "elapsed_seconds\t" $1; print "max_rss_kb\t" $2; print "exit_status\t" $3}' \
      "$OUT/metrics.tsv"
    awk '
      /Timed out/ {print "timed_out\t" $3}
      /RSS guard limit/ {print "rss_guard_limit_kb\t" $5}
      /RSS guard exceeded/ {print "rss_guard_exceeded\t" $4}
    ' "$OUT/stderr.time.log"
  } > "$OUT/time_summary.tsv"
fi

"$PYTHON" "$UNIQUE_SUMMARY_SCRIPT" \
  "rust=$OUT/stderr.time.log,$OUT/input.hist.tsv,$OUT/input.rhist.tsv" \
  > "$OUT/unique_kmers.tsv"

unique_summary_field() {
  local column="$1"
  awk -v column="$column" '
    BEGIN { FS="\t" }
    NR==1 {
      for (i=1; i<=NF; i++) { idx[$i]=i }
      next
    }
    $1=="rust" && (column in idx) { print $idx[column]; exit }
  ' "$OUT/unique_kmers.tsv"
}

countup_spill_guard_enabled=0
for limit in \
  "$MAX_COUNTUP_SPILL_INITIAL_RUNS" \
  "$MAX_COUNTUP_SPILL_MERGE_RUNS" \
  "$MAX_COUNTUP_SPILL_FINAL_RUNS" \
  "$MAX_COUNTUP_SPILL_BYTES_WRITTEN" \
  "$MAX_COUNTUP_SPILL_PEAK_LIVE_BYTES" \
  "$MAX_COUNTUP_SPILL_FINAL_LIVE_BYTES"
do
  if [[ -n "$limit" ]]; then
    countup_spill_guard_enabled=1
  fi
done

spill_exceeds_limit() {
  local value="$1"
  local limit="$2"
  [[ -n "$limit" ]] || return 1
  [[ "$value" =~ ^[0-9]+$ ]] || value=0
  (( value > limit ))
}

if [[ "$countup_spill_guard_enabled" == "1" ]]; then
  observed_spill_initial_runs="$(unique_summary_field countup_spill_initial_runs)"
  observed_spill_merge_runs="$(unique_summary_field countup_spill_merge_runs)"
  observed_spill_final_runs="$(unique_summary_field countup_spill_final_runs)"
  observed_spill_bytes_written="$(unique_summary_field countup_spill_bytes_written)"
  observed_spill_peak_live_bytes="$(unique_summary_field countup_spill_peak_live_bytes)"
  observed_spill_final_live_bytes="$(unique_summary_field countup_spill_final_live_bytes)"
  spill_guard_status=ok
  spill_guard_reason=
  if spill_exceeds_limit "${observed_spill_initial_runs:-}" "$MAX_COUNTUP_SPILL_INITIAL_RUNS"; then
    spill_guard_status=exceeded
    spill_guard_reason="${spill_guard_reason:+$spill_guard_reason,}countup_spill_initial_runs>${MAX_COUNTUP_SPILL_INITIAL_RUNS}"
  fi
  if spill_exceeds_limit "${observed_spill_merge_runs:-}" "$MAX_COUNTUP_SPILL_MERGE_RUNS"; then
    spill_guard_status=exceeded
    spill_guard_reason="${spill_guard_reason:+$spill_guard_reason,}countup_spill_merge_runs>${MAX_COUNTUP_SPILL_MERGE_RUNS}"
  fi
  if spill_exceeds_limit "${observed_spill_final_runs:-}" "$MAX_COUNTUP_SPILL_FINAL_RUNS"; then
    spill_guard_status=exceeded
    spill_guard_reason="${spill_guard_reason:+$spill_guard_reason,}countup_spill_final_runs>${MAX_COUNTUP_SPILL_FINAL_RUNS}"
  fi
  if spill_exceeds_limit "${observed_spill_bytes_written:-}" "$MAX_COUNTUP_SPILL_BYTES_WRITTEN"; then
    spill_guard_status=exceeded
    spill_guard_reason="${spill_guard_reason:+$spill_guard_reason,}countup_spill_bytes_written>${MAX_COUNTUP_SPILL_BYTES_WRITTEN}"
  fi
  if spill_exceeds_limit "${observed_spill_peak_live_bytes:-}" "$MAX_COUNTUP_SPILL_PEAK_LIVE_BYTES"; then
    spill_guard_status=exceeded
    spill_guard_reason="${spill_guard_reason:+$spill_guard_reason,}countup_spill_peak_live_bytes>${MAX_COUNTUP_SPILL_PEAK_LIVE_BYTES}"
  fi
  if spill_exceeds_limit "${observed_spill_final_live_bytes:-}" "$MAX_COUNTUP_SPILL_FINAL_LIVE_BYTES"; then
    spill_guard_status=exceeded
    spill_guard_reason="${spill_guard_reason:+$spill_guard_reason,}countup_spill_final_live_bytes>${MAX_COUNTUP_SPILL_FINAL_LIVE_BYTES}"
  fi
  {
    printf 'metric\tobserved\tlimit\n'
    printf 'countup_spill_initial_runs\t%s\t%s\n' "${observed_spill_initial_runs:-0}" "$MAX_COUNTUP_SPILL_INITIAL_RUNS"
    printf 'countup_spill_merge_runs\t%s\t%s\n' "${observed_spill_merge_runs:-0}" "$MAX_COUNTUP_SPILL_MERGE_RUNS"
    printf 'countup_spill_final_runs\t%s\t%s\n' "${observed_spill_final_runs:-0}" "$MAX_COUNTUP_SPILL_FINAL_RUNS"
    printf 'countup_spill_bytes_written\t%s\t%s\n' "${observed_spill_bytes_written:-0}" "$MAX_COUNTUP_SPILL_BYTES_WRITTEN"
    printf 'countup_spill_peak_live_bytes\t%s\t%s\n' "${observed_spill_peak_live_bytes:-0}" "$MAX_COUNTUP_SPILL_PEAK_LIVE_BYTES"
    printf 'countup_spill_final_live_bytes\t%s\t%s\n' "${observed_spill_final_live_bytes:-0}" "$MAX_COUNTUP_SPILL_FINAL_LIVE_BYTES"
    printf 'status\t%s\t\n' "$spill_guard_status"
    printf 'reason\t%s\t\n' "$spill_guard_reason"
  } > "$OUT/countup_spill_guard.tsv"
  if [[ "$spill_guard_status" == "exceeded" ]]; then
    echo "Giant-safe benchmark exceeded count-up spill guard ($spill_guard_reason): $OUT" >&2
    STATUS=125
  fi
fi

if [[ -n "$MAX_RSS_KB" ]]; then
  OBSERVED_RSS_KB="$(awk -F '\t' '$1=="max_rss_kb" {print $2; exit}' "$OUT/time_summary.tsv")"
  {
    printf 'metric\tvalue\n'
    printf 'max_rss_kb\t%s\n' "${OBSERVED_RSS_KB:-unknown}"
    printf 'max_rss_kb_guard\t%s\n' "$MAX_RSS_KB"
  } > "$OUT/rss_guard.tsv"
  if [[ "$OBSERVED_RSS_KB" =~ ^[0-9]+$ ]] && (( OBSERVED_RSS_KB > MAX_RSS_KB )); then
    printf 'status\texceeded\n' >> "$OUT/rss_guard.tsv"
    echo "Giant-safe benchmark exceeded MAX_RSS_KB=$MAX_RSS_KB with max RSS ${OBSERVED_RSS_KB} KB: $OUT" >&2
    STATUS=125
  else
    printf 'status\tok\n' >> "$OUT/rss_guard.tsv"
  fi
fi

if [[ "$STATUS" -ne 0 ]]; then
  echo "Giant-safe benchmark failed: $OUT (exit $STATUS)" >&2
  exit "$STATUS"
fi

echo "Giant-safe benchmark complete: $OUT"

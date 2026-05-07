#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT"
export PATH="$HOME/.local/bin:$PATH"

if [[ $# -lt 1 ]]; then
  cat >&2 <<'USAGE'
Usage: scripts/benchmark_java_rust_human.sh <reads_R1.fq[.gz]> [reads_R2.fq[.gz]] [outdir]

Runs vendored Java BBNorm and bbnorm-rs with matched arguments while capturing
elapsed seconds, peak RSS, logs, histograms, and output summaries. Defaults are
deliberately desktop-safe: READS=10000, TABLE_READS=$READS, WRITE_OUTPUTS=0,
THREADS=8, ZIPTHREADS=$THREADS, JAVA_XMX=24g, TIMEOUT=2h.

Set READS=all TABLE_READS=all only on a machine with enough memory and swap.
Set WRITE_OUTPUTS=1 to emit normalized FASTQ files for direct sequence-output
comparison; null outputs are used by default to avoid filling disk. Use
EXTRA_ARGS='prefilter=t k=40' for common Java/Rust mode flags, and
JAVA_EXTRA_ARGS / RUST_EXTRA_ARGS for tool-specific flags.
Optional count-up spill caps are available as MAX_COUNTUP_SPILL_INITIAL_RUNS,
MAX_COUNTUP_SPILL_MERGE_RUNS, MAX_COUNTUP_SPILL_FINAL_RUNS,
MAX_COUNTUP_SPILL_BYTES_WRITTEN, MAX_COUNTUP_SPILL_PEAK_LIVE_BYTES, and
MAX_COUNTUP_SPILL_FINAL_LIVE_BYTES.
Set RUST_MEM_AUTO_FROM_JAVA=1 to run Java first, parse its actual KCountArray
bytes, and append a Rust mem=<recommended> override. Use
RUST_MEM_AUTO_MAX_BYTES to cap that auto override on desktop-safe probes.
If Rust receives explicit sketch sizing such as sketchmemory, cells, or
matrixbits, auto-memory is skipped because those controls override mem sizing.
USAGE
  exit 2
fi

R1=$1
R2=${2:-}
OUT=${3:-tmp/java_rust_human_benchmark_$(date +%Y%m%d_%H%M%S)}
THREADS=${THREADS:-8}
ZIPTHREADS=${ZIPTHREADS:-$THREADS}
READS=${READS:-10000}
TABLE_READS=${TABLE_READS:-$READS}
WRITE_OUTPUTS=${WRITE_OUTPUTS:-0}
TIMEOUT=${TIMEOUT:-2h}
JAVA_XMX=${JAVA_XMX:-24g}
K=${K:-31}
MIN=${MIN:-5}
TARGET=${TARGET:-40}
MAX=${MAX:-80}
BITS=${BITS:-32}
MEM=${MEM:-2g}
EXTRA_ARGS=${EXTRA_ARGS:-}
JAVA_EXTRA_ARGS=${JAVA_EXTRA_ARGS:-}
RUST_EXTRA_ARGS=${RUST_EXTRA_ARGS:-}
ALLOW_JAVA_FAILURE=${ALLOW_JAVA_FAILURE:-0}
SKIP_JAVA=${SKIP_JAVA:-0}
RUST_MEM_AUTO_FROM_JAVA=${RUST_MEM_AUTO_FROM_JAVA:-0}
RUST_MEM_AUTO_MAX_BYTES=${RUST_MEM_AUTO_MAX_BYTES:-}
MAX_RSS_KB=${MAX_RSS_KB:-}
JAVA_MAX_RSS_KB=${JAVA_MAX_RSS_KB:-$MAX_RSS_KB}
RUST_MAX_RSS_KB=${RUST_MAX_RSS_KB:-$MAX_RSS_KB}
MAX_COUNTUP_SPILL_INITIAL_RUNS=${MAX_COUNTUP_SPILL_INITIAL_RUNS:-}
MAX_COUNTUP_SPILL_MERGE_RUNS=${MAX_COUNTUP_SPILL_MERGE_RUNS:-}
MAX_COUNTUP_SPILL_FINAL_RUNS=${MAX_COUNTUP_SPILL_FINAL_RUNS:-}
MAX_COUNTUP_SPILL_BYTES_WRITTEN=${MAX_COUNTUP_SPILL_BYTES_WRITTEN:-}
MAX_COUNTUP_SPILL_PEAK_LIVE_BYTES=${MAX_COUNTUP_SPILL_PEAK_LIVE_BYTES:-}
MAX_COUNTUP_SPILL_FINAL_LIVE_BYTES=${MAX_COUNTUP_SPILL_FINAL_LIVE_BYTES:-}
PYTHON=${PYTHON:-python3}
MEASURE_SCRIPT=${MEASURE_SCRIPT:-scripts/measure_command.py}
UNIQUE_SUMMARY_SCRIPT=${UNIQUE_SUMMARY_SCRIPT:-scripts/extract_unique_kmer_summary.py}
HISTOGRAM_DIFF_SCRIPT=${HISTOGRAM_DIFF_SCRIPT:-scripts/compare_histogram_tables.py}
STAGE_TIMINGS_SCRIPT=${STAGE_TIMINGS_SCRIPT:-scripts/extract_stage_timings.py}
CP="${CP:-vendor/BBTools-master/current}"

mkdir -p "$OUT"

if [[ ! -r "$R1" ]]; then
  echo "Missing readable R1 input: $R1" >&2
  exit 2
fi
if [[ -n "$R2" && ! -r "$R2" ]]; then
  echo "Missing readable R2 input: $R2" >&2
  exit 2
fi
if [[ ! -d "$CP" ]]; then
  echo "Missing vendored BBTools classpath directory: $CP" >&2
  exit 2
fi

normalize_limit_value() {
  case "${1,,}" in
    ""|"all"|"none"|"0") echo "" ;;
    *) echo "$1" ;;
  esac
}

normalize_rss_guard_value() {
  case "${1,,}" in
    ""|"0"|"none"|"off"|"unlimited") echo "" ;;
    *) echo "$1" ;;
  esac
}

validate_optional_integer() {
  local name="$1"
  local value="${!name:-}"
  if [[ -n "$value" && ! "$value" =~ ^[0-9]+$ ]]; then
    printf '%s must be an integer count when set, got: %s\n' "$name" "$value" >&2
    exit 2
  fi
}

READS_ARG="$(normalize_limit_value "$READS")"
TABLE_READS_ARG="$(normalize_limit_value "$TABLE_READS")"
MAX_RSS_KB="$(normalize_rss_guard_value "$MAX_RSS_KB")"
JAVA_MAX_RSS_KB="$(normalize_rss_guard_value "$JAVA_MAX_RSS_KB")"
RUST_MAX_RSS_KB="$(normalize_rss_guard_value "$RUST_MAX_RSS_KB")"
if [[ "$SKIP_JAVA" != "0" && "$SKIP_JAVA" != "1" ]]; then
  echo "SKIP_JAVA must be 0 or 1, got: $SKIP_JAVA" >&2
  exit 2
fi
if [[ "$RUST_MEM_AUTO_FROM_JAVA" != "0" && "$RUST_MEM_AUTO_FROM_JAVA" != "1" ]]; then
  echo "RUST_MEM_AUTO_FROM_JAVA must be 0 or 1, got: $RUST_MEM_AUTO_FROM_JAVA" >&2
  exit 2
fi
if [[ -n "$RUST_MEM_AUTO_MAX_BYTES" && ! "$RUST_MEM_AUTO_MAX_BYTES" =~ ^[0-9]+$ ]]; then
  echo "RUST_MEM_AUTO_MAX_BYTES must be an integer byte count, got: $RUST_MEM_AUTO_MAX_BYTES" >&2
  exit 2
fi
validate_optional_integer MAX_COUNTUP_SPILL_BYTES_WRITTEN
validate_optional_integer MAX_COUNTUP_SPILL_PEAK_LIVE_BYTES
validate_optional_integer MAX_COUNTUP_SPILL_FINAL_LIVE_BYTES
validate_optional_integer MAX_COUNTUP_SPILL_INITIAL_RUNS
validate_optional_integer MAX_COUNTUP_SPILL_MERGE_RUNS
validate_optional_integer MAX_COUNTUP_SPILL_FINAL_RUNS
COMMON_EXTRA_ARG_ARRAY=()
JAVA_EXTRA_ARG_ARRAY=()
RUST_EXTRA_ARG_ARRAY=()
if [[ -n "$EXTRA_ARGS" ]]; then
  read -r -a COMMON_EXTRA_ARG_ARRAY <<< "$EXTRA_ARGS"
fi
if [[ -n "$JAVA_EXTRA_ARGS" ]]; then
  read -r -a JAVA_EXTRA_ARG_ARRAY <<< "$JAVA_EXTRA_ARGS"
fi
if [[ -n "$RUST_EXTRA_ARGS" ]]; then
  read -r -a RUST_EXTRA_ARG_ARRAY <<< "$RUST_EXTRA_ARGS"
fi

cargo build --release --quiet

COMMON=(
  "in=$R1"
  "passes=1"
  "target=$TARGET"
  "max=$MAX"
  "min=$MIN"
  "k=$K"
  "threads=$THREADS"
  "zipthreads=$ZIPTHREADS"
  "bits=$BITS"
  "overwrite=t"
)
if [[ -n "$R2" ]]; then
  COMMON+=("in2=$R2")
fi
if [[ -n "$READS_ARG" ]]; then
  COMMON+=("reads=$READS_ARG")
fi
if [[ -n "$TABLE_READS_ARG" ]]; then
  COMMON+=("tablereads=$TABLE_READS_ARG")
fi

JAVA_COMMON=("${COMMON[@]}")
RUST_COMMON=("${COMMON[@]}" "mem=$MEM" "autocountmin=t")
if [[ -n "$MAX_COUNTUP_SPILL_INITIAL_RUNS" ]]; then
  RUST_COMMON+=("maxcountupspillinitialruns=$MAX_COUNTUP_SPILL_INITIAL_RUNS")
fi
if [[ -n "$MAX_COUNTUP_SPILL_MERGE_RUNS" ]]; then
  RUST_COMMON+=("maxcountupspillmergeruns=$MAX_COUNTUP_SPILL_MERGE_RUNS")
fi
if [[ -n "$MAX_COUNTUP_SPILL_FINAL_RUNS" ]]; then
  RUST_COMMON+=("maxcountupspillfinalruns=$MAX_COUNTUP_SPILL_FINAL_RUNS")
fi
if [[ -n "$MAX_COUNTUP_SPILL_PEAK_LIVE_BYTES" ]]; then
  RUST_COMMON+=("maxcountupspillbytes=$MAX_COUNTUP_SPILL_PEAK_LIVE_BYTES")
fi
if [[ -n "$MAX_COUNTUP_SPILL_FINAL_LIVE_BYTES" ]]; then
  RUST_COMMON+=("maxcountupspillfinallivebytes=$MAX_COUNTUP_SPILL_FINAL_LIVE_BYTES")
fi
if [[ -n "$MAX_COUNTUP_SPILL_BYTES_WRITTEN" ]]; then
  RUST_COMMON+=("maxcountupspillwritebytes=$MAX_COUNTUP_SPILL_BYTES_WRITTEN")
fi

sequence_outputs_for_tool() {
  local tool="$1"
  if [[ "$WRITE_OUTPUTS" == "1" ]]; then
    printf 'out=%s/%s.keep1.fq.gz\n' "$OUT" "$tool"
    printf 'outt=%s/%s.toss1.fq.gz\n' "$OUT" "$tool"
    if [[ -n "$R2" ]]; then
      printf 'out2=%s/%s.keep2.fq.gz\n' "$OUT" "$tool"
      printf 'outt2=%s/%s.toss2.fq.gz\n' "$OUT" "$tool"
    fi
  else
    printf 'out=null\n'
    printf 'outt=null\n'
    if [[ -n "$R2" ]]; then
      printf 'out2=null\n'
      printf 'outt2=null\n'
    fi
  fi
}

run_measured() {
  local tool="$1"
  shift
  local metrics="$OUT/$tool.metrics.tsv"
  local stdout="$OUT/$tool.stdout.log"
  local stderr="$OUT/$tool.stderr.log"
  local cmd=("$@")
  local cap_var
  local cap
  cap_var="$(printf '%s_MAX_RSS_KB' "${tool^^}")"
  cap="${!cap_var:-}"
  local measure_args=(
    "$PYTHON"
    "$MEASURE_SCRIPT"
    "--metrics" "$metrics"
    "--stdout" "$stdout"
    "--stderr" "$stderr"
  )
  if [[ -n "$TIMEOUT" ]]; then
    measure_args+=("--timeout" "$TIMEOUT")
  fi
  if [[ -n "$cap" ]]; then
    if ! [[ "$cap" =~ ^[0-9]+$ ]]; then
      printf '%s must be an integer number of kilobytes, got: %s\n' "$cap_var" "$cap" >&2
      return 2
    fi
    measure_args+=("--max-rss-kb" "$cap")
  fi
  printf 'Running %s...\n' "$tool"
  set +e
  "${measure_args[@]}" -- "${cmd[@]}"
  local status=$?
  set -e
  if [[ -n "$cap" ]]; then
    local observed
    observed="$(awk 'NR==1 {print $2}' "$metrics" 2>/dev/null || true)"
    {
      printf 'tool\tmax_rss_kb\tmax_rss_kb_guard\tstatus\n'
      if [[ "$observed" =~ ^[0-9]+$ ]] && (( observed > cap )); then
        printf '%s\t%s\t%s\texceeded\n' "$tool" "$observed" "$cap"
      else
        printf '%s\t%s\t%s\tok\n' "$tool" "${observed:-unknown}" "$cap"
      fi
    } > "$OUT/$tool.rss_guard.tsv"
    if [[ "$observed" =~ ^[0-9]+$ ]] && (( observed > cap )); then
      printf '%s exceeded %s=%s with max RSS %s KB. Logs: %s %s\n' "$tool" "$cap_var" "$cap" "$observed" "$stdout" "$stderr" >&2
      return 125
    fi
  fi
  if [[ "$status" -ne 0 ]]; then
    printf '%s failed with exit %s. Logs: %s %s\n' "$tool" "$status" "$stdout" "$stderr" >&2
    return "$status"
  fi
}

compare_sequence_file() {
  local java_file="$1"
  local rust_file="$2"
  [[ -f "$java_file" && -f "$rust_file" ]] || return 1
  "$PYTHON" - "$java_file" "$rust_file" <<'PY'
import gzip
import sys

def open_maybe_gzip(path):
    if path.endswith(".gz"):
        return gzip.open(path, "rb")
    return open(path, "rb")

left_path, right_path = sys.argv[1], sys.argv[2]
with open_maybe_gzip(left_path) as left, open_maybe_gzip(right_path) as right:
    while True:
        left_chunk = left.read(1024 * 1024)
        right_chunk = right.read(1024 * 1024)
        if left_chunk != right_chunk:
            sys.exit(1)
        if not left_chunk:
            sys.exit(0)
PY
}

comparison_status() {
  local java_file="$1"
  local rust_file="$2"
  if [[ "$java_status" -ne 0 ]]; then
    echo 'skipped_java_failed'
  elif [[ ! -f "$java_file" || ! -f "$rust_file" ]]; then
    echo 'missing'
  elif cmp -s "$java_file" "$rust_file"; then
    echo 'identical'
  else
    echo 'different'
  fi
}

unique_summary_field() {
  local file="$1"
  local tool="$2"
  local column="$3"
  [[ -s "$file" ]] || return 0
  awk -v tool="$tool" -v column="$column" '
    BEGIN { FS="\t" }
    NR==1 {
      for (i=1; i<=NF; i++) { idx[$i]=i }
      next
    }
    $1==tool && (column in idx) { print $idx[column]; exit }
  ' "$file"
}

rust_args_have_explicit_sketch_sizing() {
  local arg key lower
  for arg in "${COMMON_EXTRA_ARG_ARRAY[@]}" "${RUST_EXTRA_ARG_ARRAY[@]}"; do
    [[ "$arg" == *=* ]] || continue
    key="${arg%%=*}"
    key="${key#--}"
    key="${key#-}"
    lower="${key,,}"
    case "$lower" in
      cells|matrixbits|sketchmemory|sketchmem|countminmemory|countminmem|cmem|\
      prefiltercells|prefiltercell|precells|precell|filtercells|filtercell|\
      prefiltermemory|filtermemory|prefiltermem|filtermem|filtermemoryoverride)
        return 0
        ;;
    esac
  done
  return 1
}

suggest_rust_mem_bytes_for_table() {
  local table_bytes="$1"
  [[ "$table_bytes" =~ ^[0-9]+$ && "$table_bytes" -gt 0 ]] || return 0
  local hist_bytes
  hist_bytes="$(rust_hist_reserved_bytes)"
  awk -v target="$table_bytes" -v hist="$hist_bytes" '
    BEGIN {
      needed = target + hist;
      by_fraction = needed / 0.45;
      by_headroom = needed / 0.73 + 96000000;
      best = by_fraction < by_headroom ? by_fraction : by_headroom;
      printf "%d\n", int((best + 999999) / 1000000) * 1000000;
    }
  '
}

rust_hist_reserved_bytes() {
  local hist_len="${HIST_LEN:-1048577}"
  if [[ ! "$hist_len" =~ ^[0-9]+$ ]]; then
    hist_len=1048577
  elif (( hist_len < 1 )); then
    hist_len=1048577
  fi
  # Rust histograms now use compact Rayon reducers and the writer streams
  # zero-bin folding, so recommendations reserve dense whole-histogram buffers
  # rather than one histlen buffer per worker.
  awk -v hist_len="$hist_len" 'BEGIN { printf "%d\n", hist_len * 8 * 3 }'
}

format_bytes_mib() {
  local bytes="$1"
  [[ "$bytes" =~ ^[0-9]+$ ]] || return 0
  awk -v bytes="$bytes" 'BEGIN { printf "%dm\n", int((bytes + 999999) / 1000000) }'
}

mapfile -t JAVA_OUTPUTS < <(sequence_outputs_for_tool java)
mapfile -t RUST_OUTPUTS < <(sequence_outputs_for_tool rust)

JAVA_CMD=(
  java "-Xmx$JAVA_XMX" -cp "$CP" jgi.KmerNormalize
  "${JAVA_COMMON[@]}"
  "${COMMON_EXTRA_ARG_ARRAY[@]}"
  "${JAVA_EXTRA_ARG_ARRAY[@]}"
  "hist=$OUT/java.hist.tsv"
  "rhist=$OUT/java.rhist.tsv"
  "${JAVA_OUTPUTS[@]}"
)
RUST_CMD=(
  target/release/bbnorm-rs
  "${RUST_COMMON[@]}"
  "${COMMON_EXTRA_ARG_ARRAY[@]}"
  "${RUST_EXTRA_ARG_ARRAY[@]}"
  "hist=$OUT/rust.hist.tsv"
  "rhist=$OUT/rust.rhist.tsv"
  "${RUST_OUTPUTS[@]}"
)

{
  printf 'timestamp\t%s\n' "$(date -Is)"
  printf 'host\t%s\n' "$(hostname)"
  printf 'r1\t%s\n' "$R1"
  printf 'r2\t%s\n' "${R2:-}"
  printf 'r1_bytes\t%s\n' "$(stat -c %s "$R1")"
  if [[ -n "$R2" ]]; then printf 'r2_bytes\t%s\n' "$(stat -c %s "$R2")"; fi
  printf 'threads\t%s\n' "$THREADS"
  printf 'zipthreads\t%s\n' "$ZIPTHREADS"
  printf 'reads\t%s\n' "$READS"
  printf 'tablereads\t%s\n' "$TABLE_READS"
  printf 'write_outputs\t%s\n' "$WRITE_OUTPUTS"
  printf 'extra_args\t%s\n' "${EXTRA_ARGS:-}"
  printf 'java_extra_args\t%s\n' "${JAVA_EXTRA_ARGS:-}"
  printf 'rust_extra_args\t%s\n' "${RUST_EXTRA_ARGS:-}"
  printf 'allow_java_failure\t%s\n' "$ALLOW_JAVA_FAILURE"
  printf 'skip_java\t%s\n' "$SKIP_JAVA"
  printf 'rust_mem_auto_from_java\t%s\n' "$RUST_MEM_AUTO_FROM_JAVA"
  printf 'rust_mem_auto_max_bytes\t%s\n' "$RUST_MEM_AUTO_MAX_BYTES"
  printf 'java_xmx\t%s\n' "$JAVA_XMX"
  printf 'rust_mem\t%s\n' "$MEM"
  printf 'max_rss_kb_guard\t%s\n' "${MAX_RSS_KB:-}"
  printf 'java_max_rss_kb_guard\t%s\n' "${JAVA_MAX_RSS_KB:-}"
  printf 'rust_max_rss_kb_guard\t%s\n' "${RUST_MAX_RSS_KB:-}"
  printf 'max_countup_spill_initial_runs\t%s\n' "$MAX_COUNTUP_SPILL_INITIAL_RUNS"
  printf 'max_countup_spill_merge_runs\t%s\n' "$MAX_COUNTUP_SPILL_MERGE_RUNS"
  printf 'max_countup_spill_final_runs\t%s\n' "$MAX_COUNTUP_SPILL_FINAL_RUNS"
  printf 'max_countup_spill_bytes_written\t%s\n' "$MAX_COUNTUP_SPILL_BYTES_WRITTEN"
  printf 'max_countup_spill_peak_live_bytes\t%s\n' "$MAX_COUNTUP_SPILL_PEAK_LIVE_BYTES"
  printf 'max_countup_spill_final_live_bytes\t%s\n' "$MAX_COUNTUP_SPILL_FINAL_LIVE_BYTES"
  printf 'pigz\t%s\n' "$(command -v pigz || true)"
  printf 'unpigz\t%s\n' "$(command -v unpigz || true)"
  printf 'measure_script\t%s\n' "$MEASURE_SCRIPT"
  printf 'unique_summary_script\t%s\n' "$UNIQUE_SUMMARY_SCRIPT"
  printf 'histogram_diff_script\t%s\n' "$HISTOGRAM_DIFF_SCRIPT"
  printf 'stage_timings_script\t%s\n' "$STAGE_TIMINGS_SCRIPT"
  printf 'free_before\t%s\n' "$(free -h | tr '\n' ' ' | sed 's/[[:space:]]\+/ /g')"
  printf 'df_before\t%s\n' "$(df -h . | tr '\n' ' ' | sed 's/[[:space:]]\+/ /g')"
  printf 'java_command\t'
  printf '%q ' "${JAVA_CMD[@]}"
  printf '\n'
  printf 'rust_command\t'
  printf '%q ' "${RUST_CMD[@]}"
  printf '\n'
} > "$OUT/environment.tsv"

if [[ "$SKIP_JAVA" == "1" ]]; then
  printf 'Skipping Java leg because SKIP_JAVA=1; running Rust only.\n' >&2
  printf '0\t0\t126\n' > "$OUT/java.metrics.tsv"
  : > "$OUT/java.stdout.log"
  printf 'Java leg skipped by SKIP_JAVA=1.\n' > "$OUT/java.stderr.log"
  java_status=126
elif run_measured java "${JAVA_CMD[@]}"; then
  java_status=0
else
  java_status=$?
fi
if [[ "$java_status" -ne 0 ]]; then
  if [[ "$ALLOW_JAVA_FAILURE" == "1" ]]; then
    printf 'Java failed with exit %s, but ALLOW_JAVA_FAILURE=1; continuing to Rust.\n' \
      "$java_status" >&2
  else
    exit "$java_status"
  fi
fi

rust_mem_auto_status=off
rust_mem_auto_java_sketch_bytes=
rust_mem_auto_recommended_bytes=
rust_mem_auto_recommended=
if [[ "$RUST_MEM_AUTO_FROM_JAVA" == "1" ]]; then
  if rust_args_have_explicit_sketch_sizing; then
    rust_mem_auto_status=skipped_explicit_sketch
  elif [[ "$java_status" -ne 0 ]]; then
    rust_mem_auto_status=skipped_java_failed
  elif [[ "$SKIP_JAVA" == "1" ]]; then
    rust_mem_auto_status=skipped_java
  else
    "$PYTHON" "$UNIQUE_SUMMARY_SCRIPT" \
      "java=$OUT/java.stderr.log,$OUT/java.hist.tsv,$OUT/java.rhist.tsv" \
      > "$OUT/java.unique_kmers.pre_rust.tsv"
    rust_mem_auto_java_sketch_bytes="$(
      unique_summary_field "$OUT/java.unique_kmers.pre_rust.tsv" java sketch_memory_bytes
    )"
    rust_mem_auto_recommended_bytes="$(
      suggest_rust_mem_bytes_for_table "$rust_mem_auto_java_sketch_bytes"
    )"
    rust_mem_auto_recommended="$(format_bytes_mib "$rust_mem_auto_recommended_bytes")"
    if [[ -z "$rust_mem_auto_recommended_bytes" || -z "$rust_mem_auto_recommended" ]]; then
      rust_mem_auto_status=skipped_no_java_sketch
    elif [[ -n "$RUST_MEM_AUTO_MAX_BYTES" ]] \
      && (( rust_mem_auto_recommended_bytes > RUST_MEM_AUTO_MAX_BYTES )); then
      rust_mem_auto_status=skipped_above_cap
    else
      RUST_CMD+=("mem=$rust_mem_auto_recommended")
      rust_mem_auto_status=applied
    fi
  fi
fi
{
  printf 'rust_mem_auto_status\t%s\n' "$rust_mem_auto_status"
  printf 'rust_mem_auto_java_sketch_bytes\t%s\n' "$rust_mem_auto_java_sketch_bytes"
  printf 'rust_mem_auto_recommended_bytes\t%s\n' "$rust_mem_auto_recommended_bytes"
  printf 'rust_mem_auto_recommended\t%s\n' "$rust_mem_auto_recommended"
  printf 'rust_command_final\t'
  printf '%q ' "${RUST_CMD[@]}"
  printf '\n'
} >> "$OUT/environment.tsv"

run_measured rust "${RUST_CMD[@]}"

{
  printf 'tool\telapsed_seconds\tmax_rss_kb\tstatus\n'
  awk -v tool=java '{print tool "\t" $1 "\t" $2 "\t" $3}' "$OUT/java.metrics.tsv"
  awk -v tool=rust '{print tool "\t" $1 "\t" $2 "\t" $3}' "$OUT/rust.metrics.tsv"
} > "$OUT/results.tsv"

{
  printf 'hist_cmp\t%s\n' "$(comparison_status "$OUT/java.hist.tsv" "$OUT/rust.hist.tsv")"
  printf 'rhist_cmp\t%s\n' "$(comparison_status "$OUT/java.rhist.tsv" "$OUT/rust.rhist.tsv")"
} > "$OUT/comparison.tsv"

"$PYTHON" "$UNIQUE_SUMMARY_SCRIPT" \
  "java=$OUT/java.stderr.log,$OUT/java.hist.tsv,$OUT/java.rhist.tsv" \
  "rust=$OUT/rust.stderr.log,$OUT/rust.hist.tsv,$OUT/rust.rhist.tsv" \
  > "$OUT/unique_kmers.tsv"
"$PYTHON" "$HISTOGRAM_DIFF_SCRIPT" \
  "hist=$OUT/java.hist.tsv,$OUT/rust.hist.tsv" \
  "rhist=$OUT/java.rhist.tsv,$OUT/rust.rhist.tsv" \
  > "$OUT/histogram_diffs.tsv"
"$PYTHON" "$STAGE_TIMINGS_SCRIPT" \
  "java=$OUT/java.stderr.log" \
  "rust=$OUT/rust.stderr.log" \
  > "$OUT/stage_timings.tsv"

if [[ "$WRITE_OUTPUTS" == "1" ]]; then
  for suffix in keep1.fq.gz toss1.fq.gz keep2.fq.gz toss2.fq.gz; do
    java_file="$OUT/java.$suffix"
    rust_file="$OUT/rust.$suffix"
    if [[ -f "$java_file" || -f "$rust_file" ]]; then
      if [[ "$java_status" -ne 0 ]]; then
        printf '%s_cmp\tskipped_java_failed\n' "$suffix" >> "$OUT/comparison.tsv"
      elif compare_sequence_file "$java_file" "$rust_file"; then
        printf '%s_cmp\tidentical\n' "$suffix" >> "$OUT/comparison.tsv"
      else
        printf '%s_cmp\tdifferent\n' "$suffix" >> "$OUT/comparison.tsv"
      fi
    fi
  done
fi

{
  printf 'free_after\t%s\n' "$(free -h | tr '\n' ' ' | sed 's/[[:space:]]\+/ /g')"
  printf 'df_after\t%s\n' "$(df -h . | tr '\n' ' ' | sed 's/[[:space:]]\+/ /g')"
} >> "$OUT/environment.tsv"

"$PYTHON" - "$OUT/results.tsv" "$OUT/comparison.tsv" <<'PY'
import csv
import sys

rows = list(csv.DictReader(open(sys.argv[1], encoding="utf-8"), delimiter="\t"))
print("\nJava/Rust benchmark summary:")
print("tool\telapsed_seconds\tmax_rss_kb\tstatus")
for row in rows:
    print(f"{row['tool']}\t{float(row['elapsed_seconds']):.6f}\t{row['max_rss_kb']}\t{row['status']}")
print("\nOutput comparison:")
for line in open(sys.argv[2], encoding="utf-8"):
    print(line.rstrip())
PY

"$PYTHON" - "$OUT/unique_kmers.tsv" <<'PY'
import csv
import sys

rows = list(csv.DictReader(open(sys.argv[1], encoding="utf-8"), delimiter="\t"))
if rows:
    print("\nUnique-kmer estimate summary:")
    print("tool\testimated_unique_kmers\thist_unique_kmers\tlow_depth_kmers\thigh_depth_kmers\tsketch_tables\tsketch_total_cells\tsketch_memory_bytes")
    for row in rows:
        print(
            "\t".join(
                [
                    row["tool"],
                    row["estimated_unique_kmers"],
                    row["hist_unique_kmers"],
                    row["low_depth_kmers"],
                    row["high_depth_kmers"],
                    row["sketch_tables"],
                    row["sketch_total_cells"],
                    row["sketch_memory_bytes"],
                ]
            )
        )
PY

"$PYTHON" - "$OUT/histogram_diffs.tsv" <<'PY'
import csv
import sys

rows = list(csv.DictReader(open(sys.argv[1], encoding="utf-8"), delimiter="\t"))
if rows:
    print("\nHistogram delta summary:")
    print("label\tstatus\tcol2_abs_delta_sum\tcol2_abs_delta_ppm\tcol3_abs_delta_sum\tcol3_abs_delta_ppm\tfirst_diff_key")
    for row in rows:
        print(
            "\t".join(
                [
                    row["label"],
                    row["status"],
                    row["col2_abs_delta_sum"],
                    row["col2_abs_delta_ppm"],
                    row["col3_abs_delta_sum"],
                    row["col3_abs_delta_ppm"],
                    row["first_diff_key"],
                ]
            )
        )
PY

"$PYTHON" - "$OUT/stage_timings.tsv" <<'PY'
import csv
import sys

rows = list(csv.DictReader(open(sys.argv[1], encoding="utf-8"), delimiter="\t"))
if rows:
    print("\nStage timing summary:")
    print("tool\tstage\tseconds")
    for row in rows:
        print("\t".join([row["tool"], row["stage"], row["seconds"]]))
PY

printf '\nJava/Rust benchmark complete: %s\n' "$OUT"

#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT"

if [[ $# -lt 1 ]]; then
  cat >&2 <<'USAGE'
Usage: scripts/benchmark_java_rust_modes.sh <reads_R1.fq[.gz]> [reads_R2.fq[.gz]] [outdir]

Runs a desktop-safe matrix of Java BBNorm vs bbnorm-rs modes by delegating to
scripts/benchmark_java_rust_human.sh. Outputs one subdirectory per mode plus a
summary.tsv with status, elapsed seconds, peak RSS, and hist/rhist comparison.

Defaults are intentionally conservative. Override MODE_CASES to choose modes,
e.g. MODE_CASES='default prefilter k40 k40_fixspikes passes2 countup_tossbadreads'. For a custom
mode name, provide MODE_ARGS_<name>, for example:
  MODE_CASES='default mymode' MODE_ARGS_mymode='prefilter=t hashes=8'
Set MODE_PROFILE=bounded_core for the current stable approximate-sketch matrix,
MODE_PROFILE=countup_expected for count-up modes where vendored Java is expected
to fail but Rust should complete, or MODE_PROFILE=production_probe for both.

By default, any failed mode or non-identical hist/rhist comparison makes this
wrapper exit nonzero after writing summary.tsv. Set ALLOW_MODE_FAILURES=1 for
exploratory runs against known-unstable Java modes, or
REQUIRE_IDENTICAL_COMPARISONS=0 for approximate-sketch comparison sweeps.
Set EXPECTED_FAILURE_MODES='countup countup_prefilter ...' for known Java-failing modes; those
cases run the Rust leg with ALLOW_JAVA_FAILURE=1 and only fail the wrapper if
Rust fails or if Java unexpectedly succeeds with non-identical hist/rhist.
Set SKIP_EXPECTED_FAILURE_JAVA=1 to skip Java entirely for those known-broken
modes; count-up mode profiles enable that by default to keep giant probes safe.
For approximate-sketch sweeps with REQUIRE_IDENTICAL_COMPARISONS=0, optional
guards such as MAX_HIST_ABS_RAW_DELTA, MAX_HIST_ABS_UNIQUE_DELTA,
MAX_RHIST_ABS_READS_DELTA, and MAX_RHIST_ABS_BASES_DELTA fail the wrapper when
histogram drift exceeds the configured absolute-count cap.
Relative drift guards are also available as integer parts-per-million caps:
MAX_HIST_RAW_DELTA_PPM, MAX_HIST_UNIQUE_DELTA_PPM,
MAX_RHIST_READS_DELTA_PPM, and MAX_RHIST_BASES_DELTA_PPM.
Sketch-geometry guards are available as minimum Rust-vs-Java parts-per-million
ratios: MIN_RUST_JAVA_SKETCH_CELL_PPM and MIN_RUST_JAVA_SKETCH_MEMORY_PPM.
Runs are also classified as underprovisioned when the smaller Rust-vs-Java
cell/memory ratio is below SKETCH_UNDERPROVISIONED_PPM, default 500000.
Rust count-up spill guards are available as integer caps:
MAX_COUNTUP_SPILL_INITIAL_RUNS, MAX_COUNTUP_SPILL_MERGE_RUNS,
MAX_COUNTUP_SPILL_FINAL_RUNS, MAX_COUNTUP_SPILL_BYTES_WRITTEN,
MAX_COUNTUP_SPILL_PEAK_LIVE_BYTES, and MAX_COUNTUP_SPILL_FINAL_LIVE_BYTES.
Set DRIFT_GATE_PROFILE=bounded for a broad approximate-sketch matrix guard, or
DRIFT_GATE_PROFILE=strict10k for the tighter 10k default-mode regression guard.
USAGE
  exit 2
fi

R1=$1
R2=${2:-}
OUT=${3:-tmp/java_rust_modes_$(date +%Y%m%d_%H%M%S)}
MODE_CASES_WAS_SET=0
if [[ "${MODE_CASES+x}" == "x" ]]; then
  MODE_CASES_WAS_SET=1
fi
EXPECTED_FAILURE_MODES_WAS_SET=0
if [[ "${EXPECTED_FAILURE_MODES+x}" == "x" ]]; then
  EXPECTED_FAILURE_MODES_WAS_SET=1
fi
DRIFT_GATE_PROFILE_WAS_SET=0
if [[ "${DRIFT_GATE_PROFILE+x}" == "x" ]]; then
  DRIFT_GATE_PROFILE_WAS_SET=1
fi
SKIP_EXPECTED_FAILURE_JAVA_WAS_SET=0
if [[ "${SKIP_EXPECTED_FAILURE_JAVA+x}" == "x" ]]; then
  SKIP_EXPECTED_FAILURE_JAVA_WAS_SET=1
fi
MODE_PROFILE=${MODE_PROFILE:-none}
MODE_CASES=${MODE_CASES:-"default prefilter k40 k40_fixspikes passes2"}
HARNESS=${HARNESS:-scripts/benchmark_java_rust_human.sh}
ALLOW_MODE_FAILURES=${ALLOW_MODE_FAILURES:-0}
REQUIRE_IDENTICAL_COMPARISONS_WAS_SET=0
if [[ "${REQUIRE_IDENTICAL_COMPARISONS+x}" == "x" ]]; then
  REQUIRE_IDENTICAL_COMPARISONS_WAS_SET=1
fi
REQUIRE_IDENTICAL_COMPARISONS=${REQUIRE_IDENTICAL_COMPARISONS:-1}
EXPECTED_FAILURE_MODES=${EXPECTED_FAILURE_MODES:-}
DRIFT_GATE_PROFILE=${DRIFT_GATE_PROFILE:-none}
SKIP_EXPECTED_FAILURE_JAVA=${SKIP_EXPECTED_FAILURE_JAVA:-0}
MAX_HIST_ABS_RAW_DELTA=${MAX_HIST_ABS_RAW_DELTA:-}
MAX_HIST_ABS_UNIQUE_DELTA=${MAX_HIST_ABS_UNIQUE_DELTA:-}
MAX_RHIST_ABS_READS_DELTA=${MAX_RHIST_ABS_READS_DELTA:-}
MAX_RHIST_ABS_BASES_DELTA=${MAX_RHIST_ABS_BASES_DELTA:-}
MAX_HIST_RAW_DELTA_PPM=${MAX_HIST_RAW_DELTA_PPM:-}
MAX_HIST_UNIQUE_DELTA_PPM=${MAX_HIST_UNIQUE_DELTA_PPM:-}
MAX_RHIST_READS_DELTA_PPM=${MAX_RHIST_READS_DELTA_PPM:-}
MAX_RHIST_BASES_DELTA_PPM=${MAX_RHIST_BASES_DELTA_PPM:-}
MIN_RUST_JAVA_SKETCH_CELL_PPM=${MIN_RUST_JAVA_SKETCH_CELL_PPM:-}
MIN_RUST_JAVA_SKETCH_MEMORY_PPM=${MIN_RUST_JAVA_SKETCH_MEMORY_PPM:-}
SKETCH_UNDERPROVISIONED_PPM=${SKETCH_UNDERPROVISIONED_PPM:-500000}
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

set_default_if_empty() {
  local name="$1"
  local value="$2"
  if [[ -z "${!name:-}" ]]; then
    printf -v "$name" '%s' "$value"
  fi
}

append_unique_words() {
  local existing="$1"
  shift
  local word
  local output="$existing"
  for word in "$@"; do
    [[ -n "$word" ]] || continue
    case " $output " in
      *" $word "*) ;;
      *) output="${output:+$output }$word" ;;
    esac
  done
  printf '%s\n' "$output"
}

apply_mode_profile() {
  local countup_modes="countup countup_prefilter countup_tossbadreads countup_prefilter_tossbadreads"
  case "${MODE_PROFILE,,}" in
    ""|"none"|"off")
      ;;
    "bounded_core")
      if [[ "$MODE_CASES_WAS_SET" == "0" ]]; then
        MODE_CASES="default prefilter k40_fixspikes passes2"
      fi
      if [[ "$DRIFT_GATE_PROFILE_WAS_SET" == "0" ]]; then
        DRIFT_GATE_PROFILE=bounded
      fi
      ;;
    "countup_expected")
      if [[ "$MODE_CASES_WAS_SET" == "0" ]]; then
        MODE_CASES="$countup_modes"
      fi
      if [[ "$EXPECTED_FAILURE_MODES_WAS_SET" == "0" ]]; then
        EXPECTED_FAILURE_MODES="$(append_unique_words "$EXPECTED_FAILURE_MODES" $countup_modes)"
      fi
      if [[ "$SKIP_EXPECTED_FAILURE_JAVA_WAS_SET" == "0" ]]; then
        SKIP_EXPECTED_FAILURE_JAVA=1
      fi
      ;;
    "production_probe")
      if [[ "$MODE_CASES_WAS_SET" == "0" ]]; then
        MODE_CASES="default prefilter k40_fixspikes passes2 $countup_modes"
      fi
      if [[ "$EXPECTED_FAILURE_MODES_WAS_SET" == "0" ]]; then
        EXPECTED_FAILURE_MODES="$(append_unique_words "$EXPECTED_FAILURE_MODES" $countup_modes)"
      fi
      if [[ "$DRIFT_GATE_PROFILE_WAS_SET" == "0" ]]; then
        DRIFT_GATE_PROFILE=bounded
      fi
      if [[ "$SKIP_EXPECTED_FAILURE_JAVA_WAS_SET" == "0" ]]; then
        SKIP_EXPECTED_FAILURE_JAVA=1
      fi
      ;;
    *)
      printf 'Unknown MODE_PROFILE=%s; expected none, bounded_core, countup_expected, or production_probe.\n' \
        "$MODE_PROFILE" >&2
      exit 2
      ;;
  esac
}

apply_drift_gate_profile() {
  case "${DRIFT_GATE_PROFILE,,}" in
    ""|"none"|"off")
      ;;
    "bounded")
      if [[ "$REQUIRE_IDENTICAL_COMPARISONS_WAS_SET" == "0" ]]; then
        REQUIRE_IDENTICAL_COMPARISONS=0
      fi
      set_default_if_empty MAX_HIST_RAW_DELTA_PPM 5000
      set_default_if_empty MAX_HIST_UNIQUE_DELTA_PPM 5000
      set_default_if_empty MAX_RHIST_READS_DELTA_PPM 2000
      set_default_if_empty MAX_RHIST_BASES_DELTA_PPM 2000
      ;;
    "strict10k")
      if [[ "$REQUIRE_IDENTICAL_COMPARISONS_WAS_SET" == "0" ]]; then
        REQUIRE_IDENTICAL_COMPARISONS=0
      fi
      set_default_if_empty MAX_HIST_RAW_DELTA_PPM 700
      set_default_if_empty MAX_HIST_UNIQUE_DELTA_PPM 600
      set_default_if_empty MAX_RHIST_READS_DELTA_PPM 250
      set_default_if_empty MAX_RHIST_BASES_DELTA_PPM 250
      ;;
    *)
      printf 'Unknown DRIFT_GATE_PROFILE=%s; expected none, bounded, or strict10k.\n' \
        "$DRIFT_GATE_PROFILE" >&2
      exit 2
      ;;
  esac
}

apply_mode_profile
apply_drift_gate_profile

validate_optional_integer MAX_HIST_ABS_RAW_DELTA
validate_optional_integer MAX_HIST_ABS_UNIQUE_DELTA
validate_optional_integer MAX_RHIST_ABS_READS_DELTA
validate_optional_integer MAX_RHIST_ABS_BASES_DELTA
validate_optional_integer MAX_HIST_RAW_DELTA_PPM
validate_optional_integer MAX_HIST_UNIQUE_DELTA_PPM
validate_optional_integer MAX_RHIST_READS_DELTA_PPM
validate_optional_integer MAX_RHIST_BASES_DELTA_PPM
validate_optional_integer MIN_RUST_JAVA_SKETCH_CELL_PPM
validate_optional_integer MIN_RUST_JAVA_SKETCH_MEMORY_PPM
validate_optional_integer SKETCH_UNDERPROVISIONED_PPM
validate_optional_integer MAX_COUNTUP_SPILL_INITIAL_RUNS
validate_optional_integer MAX_COUNTUP_SPILL_MERGE_RUNS
validate_optional_integer MAX_COUNTUP_SPILL_FINAL_RUNS
validate_optional_integer MAX_COUNTUP_SPILL_BYTES_WRITTEN
validate_optional_integer MAX_COUNTUP_SPILL_PEAK_LIVE_BYTES
validate_optional_integer MAX_COUNTUP_SPILL_FINAL_LIVE_BYTES
if [[ "$SKIP_EXPECTED_FAILURE_JAVA" != "0" && "$SKIP_EXPECTED_FAILURE_JAVA" != "1" ]]; then
  printf 'SKIP_EXPECTED_FAILURE_JAVA must be 0 or 1, got: %s\n' \
    "$SKIP_EXPECTED_FAILURE_JAVA" >&2
  exit 2
fi

mode_args() {
  local mode="$1"
  local override_var="MODE_ARGS_${mode}"
  if [[ -n "${!override_var:-}" ]]; then
    printf '%s\n' "${!override_var}"
    return 0
  fi
  case "$mode" in
    default) printf '\n' ;;
    prefilter) printf 'prefilter=t\n' ;;
    keepall) printf 'keepall=t\n' ;;
    ecc_mark) printf 'ecc=t markuncorrectableerrors=t\n' ;;
    qtrim_right) printf 'qtrim=r trimq=10\n' ;;
    minlen) printf 'minlen=100\n' ;;
    k40) printf 'k=40\n' ;;
    k40_fixspikes) printf 'k=40 fixspikes=t\n' ;;
    passes2) printf 'passes=2\n' ;;
    passes2_ecc_mark) printf 'passes=2 ecc=t markuncorrectableerrors=t\n' ;;
    countup) printf 'countup=t\n' ;;
    countup_prefilter) printf 'countup=t prefilter=t\n' ;;
    countup_tossbadreads) printf 'countup=t tossbadreads=t\n' ;;
    countup_prefilter_tossbadreads) printf 'countup=t prefilter=t tossbadreads=t\n' ;;
    *)
      printf 'Unknown mode %s; set MODE_ARGS_%s to define it.\n' "$mode" "$mode" >&2
      return 2
      ;;
  esac
}

metric_field() {
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

comparison_field() {
  local file="$1"
  local key="$2"
  [[ -s "$file" ]] || return 0
  awk -v key="$key" '$1==key {print $2; exit}' "$file"
}

environment_field() {
  local file="$1"
  local key="$2"
  [[ -s "$file" ]] || return 0
  awk -F '\t' -v key="$key" '$1==key {print $2; exit}' "$file"
}

unique_field() {
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

histogram_diff_field() {
  local file="$1"
  local label="$2"
  local column="$3"
  [[ -s "$file" ]] || return 0
  awk -v label="$label" -v column="$column" '
    BEGIN { FS="\t" }
    NR==1 {
      for (i=1; i<=NF; i++) { idx[$i]=i }
      next
    }
    $1==label && (column in idx) { print $idx[column]; exit }
  ' "$file"
}

stage_field() {
  local file="$1"
  local tool="$2"
  local stage="$3"
  [[ -s "$file" ]] || return 0
  awk -v tool="$tool" -v stage="$stage" '
    BEGIN { FS="\t" }
    NR==1 {
      for (i=1; i<=NF; i++) { idx[$i]=i }
      next
    }
    $1==tool && $2==stage && ("seconds" in idx) { print $idx["seconds"]; exit }
  ' "$file"
}

numeric_delta() {
  local left="$1"
  local right="$2"
  if [[ "$left" =~ ^[0-9]+$ && "$right" =~ ^[0-9]+$ ]]; then
    printf '%s\n' "$((right - left))"
  fi
}

ratio_ppm() {
  local numerator="$1"
  local denominator="$2"
  if [[ "$numerator" =~ ^[0-9]+$ && "$denominator" =~ ^[0-9]+$ && "$denominator" -gt 0 ]]; then
    awk -v numerator="$numerator" -v denominator="$denominator" \
      'BEGIN { printf "%d\n", int((numerator * 1000000 / denominator) + 0.5) }'
  fi
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

min_present_integer() {
  local left="$1"
  local right="$2"
  if [[ "$left" =~ ^[0-9]+$ && "$right" =~ ^[0-9]+$ ]]; then
    if (( left < right )); then
      printf '%s\n' "$left"
    else
      printf '%s\n' "$right"
    fi
  elif [[ "$left" =~ ^[0-9]+$ ]]; then
    printf '%s\n' "$left"
  elif [[ "$right" =~ ^[0-9]+$ ]]; then
    printf '%s\n' "$right"
  fi
}

drift_exceeds_limit() {
  local value="$1"
  local limit="$2"
  [[ -n "$limit" ]] || return 1
  [[ "$value" =~ ^[0-9]+$ ]] || return 0
  (( value > limit ))
}

resource_exceeds_limit() {
  local value="$1"
  local limit="$2"
  [[ -n "$limit" ]] || return 1
  [[ "$value" =~ ^[0-9]+$ ]] || value=0
  (( value > limit ))
}

append_reason() {
  local current="$1"
  local reason="$2"
  if [[ -n "$current" ]]; then
    printf '%s,%s\n' "$current" "$reason"
  else
    printf '%s\n' "$reason"
  fi
}

mode_is_expected_failure() {
  local mode="$1"
  local expected
  for expected in $EXPECTED_FAILURE_MODES; do
    if [[ "$expected" == "$mode" ]]; then
      return 0
    fi
  done
  return 1
}

{
  printf 'timestamp\t%s\n' "$(date -Is)"
  printf 'r1\t%s\n' "$R1"
  printf 'r2\t%s\n' "${R2:-}"
  printf 'mode_profile\t%s\n' "$MODE_PROFILE"
  printf 'mode_cases\t%s\n' "$MODE_CASES"
  printf 'harness\t%s\n' "$HARNESS"
  printf 'reads\t%s\n' "${READS:-}"
  printf 'tablereads\t%s\n' "${TABLE_READS:-}"
  printf 'threads\t%s\n' "${THREADS:-}"
  printf 'zipthreads\t%s\n' "${ZIPTHREADS:-}"
  printf 'write_outputs\t%s\n' "${WRITE_OUTPUTS:-}"
  printf 'rust_mem_auto_from_java\t%s\n' "${RUST_MEM_AUTO_FROM_JAVA:-}"
  printf 'rust_mem_auto_max_bytes\t%s\n' "${RUST_MEM_AUTO_MAX_BYTES:-}"
  printf 'java_max_rss_kb\t%s\n' "${JAVA_MAX_RSS_KB:-${MAX_RSS_KB:-}}"
  printf 'rust_max_rss_kb\t%s\n' "${RUST_MAX_RSS_KB:-${MAX_RSS_KB:-}}"
  printf 'allow_mode_failures\t%s\n' "$ALLOW_MODE_FAILURES"
  printf 'drift_gate_profile\t%s\n' "$DRIFT_GATE_PROFILE"
  printf 'skip_expected_failure_java\t%s\n' "$SKIP_EXPECTED_FAILURE_JAVA"
  printf 'require_identical_comparisons\t%s\n' "$REQUIRE_IDENTICAL_COMPARISONS"
  printf 'expected_failure_modes\t%s\n' "$EXPECTED_FAILURE_MODES"
  printf 'max_hist_abs_raw_delta\t%s\n' "$MAX_HIST_ABS_RAW_DELTA"
  printf 'max_hist_abs_unique_delta\t%s\n' "$MAX_HIST_ABS_UNIQUE_DELTA"
  printf 'max_rhist_abs_reads_delta\t%s\n' "$MAX_RHIST_ABS_READS_DELTA"
  printf 'max_rhist_abs_bases_delta\t%s\n' "$MAX_RHIST_ABS_BASES_DELTA"
  printf 'max_hist_raw_delta_ppm\t%s\n' "$MAX_HIST_RAW_DELTA_PPM"
  printf 'max_hist_unique_delta_ppm\t%s\n' "$MAX_HIST_UNIQUE_DELTA_PPM"
  printf 'max_rhist_reads_delta_ppm\t%s\n' "$MAX_RHIST_READS_DELTA_PPM"
  printf 'max_rhist_bases_delta_ppm\t%s\n' "$MAX_RHIST_BASES_DELTA_PPM"
  printf 'min_rust_java_sketch_cell_ppm\t%s\n' "$MIN_RUST_JAVA_SKETCH_CELL_PPM"
  printf 'min_rust_java_sketch_memory_ppm\t%s\n' "$MIN_RUST_JAVA_SKETCH_MEMORY_PPM"
  printf 'sketch_underprovisioned_ppm\t%s\n' "$SKETCH_UNDERPROVISIONED_PPM"
  printf 'max_countup_spill_initial_runs\t%s\n' "$MAX_COUNTUP_SPILL_INITIAL_RUNS"
  printf 'max_countup_spill_merge_runs\t%s\n' "$MAX_COUNTUP_SPILL_MERGE_RUNS"
  printf 'max_countup_spill_final_runs\t%s\n' "$MAX_COUNTUP_SPILL_FINAL_RUNS"
  printf 'max_countup_spill_bytes_written\t%s\n' "$MAX_COUNTUP_SPILL_BYTES_WRITTEN"
  printf 'max_countup_spill_peak_live_bytes\t%s\n' "$MAX_COUNTUP_SPILL_PEAK_LIVE_BYTES"
  printf 'max_countup_spill_final_live_bytes\t%s\n' "$MAX_COUNTUP_SPILL_FINAL_LIVE_BYTES"
} > "$OUT/config.tsv"

summary="$OUT/summary.tsv"
printf 'mode\tstatus\texpected_failure\textra_args\tjava_status\trust_status\tjava_seconds\tjava_rss_kb\trust_seconds\trust_rss_kb\thist_cmp\trhist_cmp\thist_abs_raw_delta\thist_raw_delta_ppm\thist_abs_unique_delta\thist_unique_delta_ppm\trhist_abs_reads_delta\trhist_reads_delta_ppm\trhist_abs_bases_delta\trhist_bases_delta_ppm\tdrift_gate\tdrift_gate_reason\tdrift_classification\tjava_unique_kmers\trust_unique_kmers\tunique_delta\tjava_hist_unique_kmers\trust_hist_unique_kmers\thist_unique_delta\tjava_low_depth_kmers\trust_low_depth_kmers\tjava_high_depth_kmers\trust_high_depth_kmers\tjava_sketch_tables\trust_sketch_tables\tjava_sketch_total_cells\trust_sketch_total_cells\tjava_sketch_memory_bytes\trust_sketch_memory_bytes\tsketch_cell_ratio_ppm\tsketch_memory_ratio_ppm\trust_mem_for_java_sketch_bytes\trust_mem_for_java_sketch\trust_mem_auto_status\trust_mem_auto_java_sketch_bytes\trust_mem_auto_recommended_bytes\trust_mem_auto_recommended\tjava_table_creation_s\tjava_table_read_s\tjava_total_s\trust_input_counting_s\trust_input_exact_counting_s\trust_input_prefilter_counting_s\trust_input_main_counting_s\trust_input_hist_s\trust_input_rhist_s\trust_normalize_s\trust_summary_counts_s\trust_output_hist_s\trust_output_rhist_s\trust_countup_work_source_s\trust_countup_normalize_s\trust_countup_spill_initial_runs\trust_countup_spill_merge_runs\trust_countup_spill_final_runs\trust_countup_spill_bytes_written\trust_countup_spill_peak_live_bytes\trust_countup_spill_final_live_bytes\tartifact_dir\n' > "$summary"

failures=0
for mode in $MODE_CASES; do
  args="$(mode_args "$mode")"
  mode_out="$OUT/$mode"
  expected_failure=0
  if mode_is_expected_failure "$mode"; then
    expected_failure=1
  fi
  skip_java=0
  if [[ "$expected_failure" == "1" && "$SKIP_EXPECTED_FAILURE_JAVA" == "1" ]]; then
    skip_java=1
  fi
  printf 'Running Java/Rust mode %s with EXTRA_ARGS=%q...\n' "$mode" "$args"
  set +e
  harness_env=(
    "ALLOW_JAVA_FAILURE=$expected_failure"
    "SKIP_JAVA=$skip_java"
    "EXTRA_ARGS=$args"
    "MAX_COUNTUP_SPILL_INITIAL_RUNS=$MAX_COUNTUP_SPILL_INITIAL_RUNS"
    "MAX_COUNTUP_SPILL_MERGE_RUNS=$MAX_COUNTUP_SPILL_MERGE_RUNS"
    "MAX_COUNTUP_SPILL_FINAL_RUNS=$MAX_COUNTUP_SPILL_FINAL_RUNS"
    "MAX_COUNTUP_SPILL_BYTES_WRITTEN=$MAX_COUNTUP_SPILL_BYTES_WRITTEN"
    "MAX_COUNTUP_SPILL_PEAK_LIVE_BYTES=$MAX_COUNTUP_SPILL_PEAK_LIVE_BYTES"
    "MAX_COUNTUP_SPILL_FINAL_LIVE_BYTES=$MAX_COUNTUP_SPILL_FINAL_LIVE_BYTES"
  )
  if [[ -n "$R2" ]]; then
    env "${harness_env[@]}" "$HARNESS" "$R1" "$R2" "$mode_out"
  else
    env "${harness_env[@]}" "$HARNESS" "$R1" "" "$mode_out"
  fi
  status=$?
  set -e

  java_status="$(metric_field "$mode_out/results.tsv" java status)"
  rust_status="$(metric_field "$mode_out/results.tsv" rust status)"
  java_seconds="$(metric_field "$mode_out/results.tsv" java elapsed_seconds)"
  java_rss="$(metric_field "$mode_out/results.tsv" java max_rss_kb)"
  rust_seconds="$(metric_field "$mode_out/results.tsv" rust elapsed_seconds)"
  rust_rss="$(metric_field "$mode_out/results.tsv" rust max_rss_kb)"
  hist_cmp="$(comparison_field "$mode_out/comparison.tsv" hist_cmp)"
  rhist_cmp="$(comparison_field "$mode_out/comparison.tsv" rhist_cmp)"
  java_unique="$(unique_field "$mode_out/unique_kmers.tsv" java estimated_unique_kmers)"
  rust_unique="$(unique_field "$mode_out/unique_kmers.tsv" rust estimated_unique_kmers)"
  unique_delta="$(numeric_delta "${java_unique:-}" "${rust_unique:-}")"
  java_hist_unique="$(unique_field "$mode_out/unique_kmers.tsv" java hist_unique_kmers)"
  rust_hist_unique="$(unique_field "$mode_out/unique_kmers.tsv" rust hist_unique_kmers)"
  hist_unique_delta="$(numeric_delta "${java_hist_unique:-}" "${rust_hist_unique:-}")"
  java_low_unique="$(unique_field "$mode_out/unique_kmers.tsv" java low_depth_kmers)"
  rust_low_unique="$(unique_field "$mode_out/unique_kmers.tsv" rust low_depth_kmers)"
  java_high_unique="$(unique_field "$mode_out/unique_kmers.tsv" java high_depth_kmers)"
  rust_high_unique="$(unique_field "$mode_out/unique_kmers.tsv" rust high_depth_kmers)"
  java_sketch_tables="$(unique_field "$mode_out/unique_kmers.tsv" java sketch_tables)"
  rust_sketch_tables="$(unique_field "$mode_out/unique_kmers.tsv" rust sketch_tables)"
  java_sketch_total_cells="$(unique_field "$mode_out/unique_kmers.tsv" java sketch_total_cells)"
  rust_sketch_total_cells="$(unique_field "$mode_out/unique_kmers.tsv" rust sketch_total_cells)"
  java_sketch_memory_bytes="$(unique_field "$mode_out/unique_kmers.tsv" java sketch_memory_bytes)"
  rust_sketch_memory_bytes="$(unique_field "$mode_out/unique_kmers.tsv" rust sketch_memory_bytes)"
  sketch_cell_ratio_ppm="$(ratio_ppm "${rust_sketch_total_cells:-}" "${java_sketch_total_cells:-}")"
  sketch_memory_ratio_ppm="$(ratio_ppm "${rust_sketch_memory_bytes:-}" "${java_sketch_memory_bytes:-}")"
  sketch_min_ratio_ppm="$(min_present_integer "${sketch_cell_ratio_ppm:-}" "${sketch_memory_ratio_ppm:-}")"
  rust_mem_for_java_sketch_bytes="$(suggest_rust_mem_bytes_for_table "${java_sketch_memory_bytes:-}")"
  rust_mem_for_java_sketch="$(format_bytes_mib "${rust_mem_for_java_sketch_bytes:-}")"
  rust_mem_auto_status="$(environment_field "$mode_out/environment.tsv" rust_mem_auto_status)"
  rust_mem_auto_java_sketch_bytes="$(
    environment_field "$mode_out/environment.tsv" rust_mem_auto_java_sketch_bytes
  )"
  rust_mem_auto_recommended_bytes="$(
    environment_field "$mode_out/environment.tsv" rust_mem_auto_recommended_bytes
  )"
  rust_mem_auto_recommended="$(environment_field "$mode_out/environment.tsv" rust_mem_auto_recommended)"
  java_table_creation_s="$(stage_field "$mode_out/stage_timings.tsv" java table_creation)"
  java_table_read_s="$(stage_field "$mode_out/stage_timings.tsv" java table_read)"
  java_total_s="$(stage_field "$mode_out/stage_timings.tsv" java total)"
  rust_input_counting_s="$(stage_field "$mode_out/stage_timings.tsv" rust input_counting)"
  rust_input_exact_counting_s="$(stage_field "$mode_out/stage_timings.tsv" rust input_exact_counting)"
  rust_input_prefilter_counting_s="$(stage_field "$mode_out/stage_timings.tsv" rust input_prefilter_counting)"
  rust_input_main_counting_s="$(stage_field "$mode_out/stage_timings.tsv" rust input_main_counting)"
  rust_input_hist_s="$(stage_field "$mode_out/stage_timings.tsv" rust input_hist)"
  rust_input_rhist_s="$(stage_field "$mode_out/stage_timings.tsv" rust input_rhist)"
  rust_normalize_s="$(stage_field "$mode_out/stage_timings.tsv" rust normalize)"
  rust_summary_counts_s="$(stage_field "$mode_out/stage_timings.tsv" rust summary_counts)"
  rust_output_hist_s="$(stage_field "$mode_out/stage_timings.tsv" rust output_hist)"
  rust_output_rhist_s="$(stage_field "$mode_out/stage_timings.tsv" rust output_rhist)"
  rust_countup_work_source_s="$(stage_field "$mode_out/stage_timings.tsv" rust countup_work_source)"
  rust_countup_normalize_s="$(stage_field "$mode_out/stage_timings.tsv" rust countup_normalize)"
  rust_countup_spill_initial_runs="$(
    unique_field "$mode_out/unique_kmers.tsv" rust countup_spill_initial_runs
  )"
  rust_countup_spill_merge_runs="$(
    unique_field "$mode_out/unique_kmers.tsv" rust countup_spill_merge_runs
  )"
  rust_countup_spill_final_runs="$(
    unique_field "$mode_out/unique_kmers.tsv" rust countup_spill_final_runs
  )"
  rust_countup_spill_bytes_written="$(
    unique_field "$mode_out/unique_kmers.tsv" rust countup_spill_bytes_written
  )"
  rust_countup_spill_peak_live_bytes="$(
    unique_field "$mode_out/unique_kmers.tsv" rust countup_spill_peak_live_bytes
  )"
  rust_countup_spill_final_live_bytes="$(
    unique_field "$mode_out/unique_kmers.tsv" rust countup_spill_final_live_bytes
  )"
  hist_abs_raw_delta="$(histogram_diff_field "$mode_out/histogram_diffs.tsv" hist col2_abs_delta_sum)"
  hist_abs_unique_delta="$(histogram_diff_field "$mode_out/histogram_diffs.tsv" hist col3_abs_delta_sum)"
  rhist_abs_reads_delta="$(histogram_diff_field "$mode_out/histogram_diffs.tsv" rhist col2_abs_delta_sum)"
  rhist_abs_bases_delta="$(histogram_diff_field "$mode_out/histogram_diffs.tsv" rhist col3_abs_delta_sum)"
  hist_raw_delta_ppm="$(histogram_diff_field "$mode_out/histogram_diffs.tsv" hist col2_abs_delta_ppm)"
  hist_unique_delta_ppm="$(histogram_diff_field "$mode_out/histogram_diffs.tsv" hist col3_abs_delta_ppm)"
  rhist_reads_delta_ppm="$(histogram_diff_field "$mode_out/histogram_diffs.tsv" rhist col2_abs_delta_ppm)"
  rhist_bases_delta_ppm="$(histogram_diff_field "$mode_out/histogram_diffs.tsv" rhist col3_abs_delta_ppm)"
  drift_gate=ok
  drift_gate_reason=
  countup_spill_gate_failed=0
  if [[ "$expected_failure" == "1" && "${java_status:-}" != "0" ]]; then
    drift_gate=skipped_java_failed
  else
    if drift_exceeds_limit "${hist_abs_raw_delta:-}" "$MAX_HIST_ABS_RAW_DELTA"; then
      drift_gate=fail
      drift_gate_reason="$(append_reason "$drift_gate_reason" "hist_abs_raw_delta>${MAX_HIST_ABS_RAW_DELTA}")"
    fi
    if drift_exceeds_limit "${hist_abs_unique_delta:-}" "$MAX_HIST_ABS_UNIQUE_DELTA"; then
      drift_gate=fail
      drift_gate_reason="$(append_reason "$drift_gate_reason" "hist_abs_unique_delta>${MAX_HIST_ABS_UNIQUE_DELTA}")"
    fi
    if drift_exceeds_limit "${rhist_abs_reads_delta:-}" "$MAX_RHIST_ABS_READS_DELTA"; then
      drift_gate=fail
      drift_gate_reason="$(append_reason "$drift_gate_reason" "rhist_abs_reads_delta>${MAX_RHIST_ABS_READS_DELTA}")"
    fi
    if drift_exceeds_limit "${rhist_abs_bases_delta:-}" "$MAX_RHIST_ABS_BASES_DELTA"; then
      drift_gate=fail
      drift_gate_reason="$(append_reason "$drift_gate_reason" "rhist_abs_bases_delta>${MAX_RHIST_ABS_BASES_DELTA}")"
    fi
    if drift_exceeds_limit "${hist_raw_delta_ppm:-}" "$MAX_HIST_RAW_DELTA_PPM"; then
      drift_gate=fail
      drift_gate_reason="$(append_reason "$drift_gate_reason" "hist_raw_delta_ppm>${MAX_HIST_RAW_DELTA_PPM}")"
    fi
    if drift_exceeds_limit "${hist_unique_delta_ppm:-}" "$MAX_HIST_UNIQUE_DELTA_PPM"; then
      drift_gate=fail
      drift_gate_reason="$(append_reason "$drift_gate_reason" "hist_unique_delta_ppm>${MAX_HIST_UNIQUE_DELTA_PPM}")"
    fi
    if drift_exceeds_limit "${rhist_reads_delta_ppm:-}" "$MAX_RHIST_READS_DELTA_PPM"; then
      drift_gate=fail
      drift_gate_reason="$(append_reason "$drift_gate_reason" "rhist_reads_delta_ppm>${MAX_RHIST_READS_DELTA_PPM}")"
    fi
    if drift_exceeds_limit "${rhist_bases_delta_ppm:-}" "$MAX_RHIST_BASES_DELTA_PPM"; then
      drift_gate=fail
      drift_gate_reason="$(append_reason "$drift_gate_reason" "rhist_bases_delta_ppm>${MAX_RHIST_BASES_DELTA_PPM}")"
    fi
    if [[ -n "$MIN_RUST_JAVA_SKETCH_CELL_PPM" && "$sketch_cell_ratio_ppm" =~ ^[0-9]+$ ]] \
      && (( sketch_cell_ratio_ppm < MIN_RUST_JAVA_SKETCH_CELL_PPM )); then
      drift_gate=fail
      drift_gate_reason="$(append_reason "$drift_gate_reason" "sketch_cell_ratio_ppm<${MIN_RUST_JAVA_SKETCH_CELL_PPM}")"
    fi
    if [[ -n "$MIN_RUST_JAVA_SKETCH_MEMORY_PPM" && "$sketch_memory_ratio_ppm" =~ ^[0-9]+$ ]] \
      && (( sketch_memory_ratio_ppm < MIN_RUST_JAVA_SKETCH_MEMORY_PPM )); then
      drift_gate=fail
      drift_gate_reason="$(append_reason "$drift_gate_reason" "sketch_memory_ratio_ppm<${MIN_RUST_JAVA_SKETCH_MEMORY_PPM}")"
    fi
  fi
  if resource_exceeds_limit "${rust_countup_spill_initial_runs:-}" "$MAX_COUNTUP_SPILL_INITIAL_RUNS"; then
    drift_gate=fail
    countup_spill_gate_failed=1
    drift_gate_reason="$(append_reason "$drift_gate_reason" "countup_spill_initial_runs>${MAX_COUNTUP_SPILL_INITIAL_RUNS}")"
  fi
  if resource_exceeds_limit "${rust_countup_spill_merge_runs:-}" "$MAX_COUNTUP_SPILL_MERGE_RUNS"; then
    drift_gate=fail
    countup_spill_gate_failed=1
    drift_gate_reason="$(append_reason "$drift_gate_reason" "countup_spill_merge_runs>${MAX_COUNTUP_SPILL_MERGE_RUNS}")"
  fi
  if resource_exceeds_limit "${rust_countup_spill_final_runs:-}" "$MAX_COUNTUP_SPILL_FINAL_RUNS"; then
    drift_gate=fail
    countup_spill_gate_failed=1
    drift_gate_reason="$(append_reason "$drift_gate_reason" "countup_spill_final_runs>${MAX_COUNTUP_SPILL_FINAL_RUNS}")"
  fi
  if resource_exceeds_limit "${rust_countup_spill_bytes_written:-}" "$MAX_COUNTUP_SPILL_BYTES_WRITTEN"; then
    drift_gate=fail
    countup_spill_gate_failed=1
    drift_gate_reason="$(append_reason "$drift_gate_reason" "countup_spill_bytes_written>${MAX_COUNTUP_SPILL_BYTES_WRITTEN}")"
  fi
  if resource_exceeds_limit "${rust_countup_spill_peak_live_bytes:-}" "$MAX_COUNTUP_SPILL_PEAK_LIVE_BYTES"; then
    drift_gate=fail
    countup_spill_gate_failed=1
    drift_gate_reason="$(append_reason "$drift_gate_reason" "countup_spill_peak_live_bytes>${MAX_COUNTUP_SPILL_PEAK_LIVE_BYTES}")"
  fi
  if resource_exceeds_limit "${rust_countup_spill_final_live_bytes:-}" "$MAX_COUNTUP_SPILL_FINAL_LIVE_BYTES"; then
    drift_gate=fail
    countup_spill_gate_failed=1
    drift_gate_reason="$(append_reason "$drift_gate_reason" "countup_spill_final_live_bytes>${MAX_COUNTUP_SPILL_FINAL_LIVE_BYTES}")"
  fi
  drift_classification=ok
  if [[ "$countup_spill_gate_failed" == "1" ]]; then
    drift_classification=countup_spill_guard
  elif [[ "$drift_gate" == "skipped_java_failed" ]]; then
    drift_classification=skipped_java_failed
  elif [[ "${sketch_min_ratio_ppm:-}" =~ ^[0-9]+$ ]] \
    && (( sketch_min_ratio_ppm < SKETCH_UNDERPROVISIONED_PPM )); then
    if [[ "$drift_gate" == "fail" ]]; then
      drift_classification=underprovisioned_sketch_drift
    else
      drift_classification=underprovisioned_sketch_ok
    fi
  elif [[ "$drift_gate" == "fail" ]]; then
    drift_classification=quantified_drift
  fi

  row=(
    "$mode" "$status" "$expected_failure" "$args" "${java_status:-}" "${rust_status:-}"
    "${java_seconds:-}" "${java_rss:-}" "${rust_seconds:-}" "${rust_rss:-}"
    "${hist_cmp:-}" "${rhist_cmp:-}" "${hist_abs_raw_delta:-}" "${hist_raw_delta_ppm:-}"
    "${hist_abs_unique_delta:-}" "${hist_unique_delta_ppm:-}"
    "${rhist_abs_reads_delta:-}" "${rhist_reads_delta_ppm:-}"
    "${rhist_abs_bases_delta:-}" "${rhist_bases_delta_ppm:-}"
    "$drift_gate" "$drift_gate_reason" "$drift_classification"
    "${java_unique:-}" "${rust_unique:-}"
    "${unique_delta:-}" "${java_hist_unique:-}" "${rust_hist_unique:-}"
    "${hist_unique_delta:-}" "${java_low_unique:-}" "${rust_low_unique:-}"
    "${java_high_unique:-}" "${rust_high_unique:-}"
    "${java_sketch_tables:-}" "${rust_sketch_tables:-}"
    "${java_sketch_total_cells:-}" "${rust_sketch_total_cells:-}"
    "${java_sketch_memory_bytes:-}" "${rust_sketch_memory_bytes:-}"
    "${sketch_cell_ratio_ppm:-}" "${sketch_memory_ratio_ppm:-}"
    "${rust_mem_for_java_sketch_bytes:-}" "${rust_mem_for_java_sketch:-}"
    "${rust_mem_auto_status:-}" "${rust_mem_auto_java_sketch_bytes:-}"
    "${rust_mem_auto_recommended_bytes:-}" "${rust_mem_auto_recommended:-}"
    "${java_table_creation_s:-}" "${java_table_read_s:-}" "${java_total_s:-}"
    "${rust_input_counting_s:-}" "${rust_input_exact_counting_s:-}"
    "${rust_input_prefilter_counting_s:-}" "${rust_input_main_counting_s:-}"
    "${rust_input_hist_s:-}" "${rust_input_rhist_s:-}"
    "${rust_normalize_s:-}" "${rust_summary_counts_s:-}"
    "${rust_output_hist_s:-}" "${rust_output_rhist_s:-}"
    "${rust_countup_work_source_s:-}" "${rust_countup_normalize_s:-}"
    "${rust_countup_spill_initial_runs:-}" "${rust_countup_spill_merge_runs:-}"
    "${rust_countup_spill_final_runs:-}" "${rust_countup_spill_bytes_written:-}"
    "${rust_countup_spill_peak_live_bytes:-}" "${rust_countup_spill_final_live_bytes:-}"
    "$mode_out"
  )
  (
    IFS=$'\t'
    printf '%s\n' "${row[*]}"
  ) >> "$summary"

  if [[ "$expected_failure" == "1" && "${java_status:-}" != "0" ]]; then
    if [[ "${rust_status:-}" != "0" || "$drift_gate" == "fail" ]]; then
      failures=$((failures + 1))
    fi
  elif [[ "$status" -ne 0 ]]; then
    failures=$((failures + 1))
  elif [[ "$drift_gate" == "fail" ]]; then
    failures=$((failures + 1))
  elif [[ "$REQUIRE_IDENTICAL_COMPARISONS" != "0" ]]; then
    if [[ "${hist_cmp:-}" != "identical" || "${rhist_cmp:-}" != "identical" ]]; then
      failures=$((failures + 1))
    fi
  fi
done

cat "$summary"

if (( failures > 0 )) && [[ "$ALLOW_MODE_FAILURES" != "1" ]]; then
  printf '\n%s Java/Rust mode case(s) failed or produced non-identical hist/rhist comparisons. See %s\n' \
    "$failures" "$summary" >&2
  exit 1
fi

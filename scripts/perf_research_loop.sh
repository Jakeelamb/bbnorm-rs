#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT"
export PATH="$HOME/.local/bin:$PATH"

R1=${1:-$ROOT/tmp/human_benchmark_8threads/human_GRCh38_500k_R1.fq.gz}
R2=${2:-$ROOT/tmp/human_benchmark_8threads/human_GRCh38_500k_R2.fq.gz}
OUT=${3:-$ROOT/tmp/perf_research_$(date +%Y%m%d_%H%M%S)}
THREADS=${THREADS:-8}
ZIPTHREADS=${ZIPTHREADS:-1}
READS=${READS:-50000}
TABLE_READS=${TABLE_READS:-$READS}
JAVA_XMX=${JAVA_XMX:-4g}
MEM=${MEM:-768m}
TIMEOUT=${TIMEOUT:-10m}
WRITE_OUTPUTS=${WRITE_OUTPUTS:-0}
JAVA_MAX_RSS_KB=${JAVA_MAX_RSS_KB:-8000000}
RUST_MAX_RSS_KB=${RUST_MAX_RSS_KB:-6000000}
DRIFT_GATE_PROFILE=${DRIFT_GATE_PROFILE:-bounded}
RUST_MEM_AUTO_FROM_JAVA=${RUST_MEM_AUTO_FROM_JAVA:-1}
RUST_MEM_AUTO_MAX_BYTES=${RUST_MEM_AUTO_MAX_BYTES:-5000000000}
MODE_CASES=${MODE_CASES:-"default prefilter k40_fixspikes passes2"}
ALLOW_MODE_FAILURES=${ALLOW_MODE_FAILURES:-1}
REQUIRE_IDENTICAL_COMPARISONS=${REQUIRE_IDENTICAL_COMPARISONS:-0}
HARNESS=${HARNESS:-scripts/benchmark_java_rust_modes.sh}
PYTHON=${PYTHON:-python3}

mkdir -p "$OUT"

if [[ ! -r "$R1" || ! -r "$R2" ]]; then
  echo "Missing input dataset(s)." >&2
  exit 2
fi

MODE_CASES="$MODE_CASES" \
READS="$READS" \
TABLE_READS="$TABLE_READS" \
THREADS="$THREADS" \
ZIPTHREADS="$ZIPTHREADS" \
JAVA_XMX="$JAVA_XMX" \
MEM="$MEM" \
TIMEOUT="$TIMEOUT" \
WRITE_OUTPUTS="$WRITE_OUTPUTS" \
JAVA_MAX_RSS_KB="$JAVA_MAX_RSS_KB" \
RUST_MAX_RSS_KB="$RUST_MAX_RSS_KB" \
DRIFT_GATE_PROFILE="$DRIFT_GATE_PROFILE" \
RUST_MEM_AUTO_FROM_JAVA="$RUST_MEM_AUTO_FROM_JAVA" \
RUST_MEM_AUTO_MAX_BYTES="$RUST_MEM_AUTO_MAX_BYTES" \
ALLOW_MODE_FAILURES="$ALLOW_MODE_FAILURES" \
REQUIRE_IDENTICAL_COMPARISONS="$REQUIRE_IDENTICAL_COMPARISONS" \
"$HARNESS" "$R1" "$R2" "$OUT/modes"

"$PYTHON" scripts/summarize_perf_research.py "$OUT/modes/summary.tsv" "$OUT/research_summary.md"

cat > "$OUT/research_config.tsv" <<EOF
key	value
r1	$R1
r2	$R2
threads	$THREADS
zipthreads	$ZIPTHREADS
reads	$READS
tablereads	$TABLE_READS
java_xmx	$JAVA_XMX
mem	$MEM
timeout	$TIMEOUT
mode_cases	$MODE_CASES
drift_gate_profile	$DRIFT_GATE_PROFILE
rust_mem_auto_from_java	$RUST_MEM_AUTO_FROM_JAVA
rust_mem_auto_max_bytes	$RUST_MEM_AUTO_MAX_BYTES
require_identical_comparisons	$REQUIRE_IDENTICAL_COMPARISONS
allow_mode_failures	$ALLOW_MODE_FAILURES
EOF

echo "Performance research complete:"
echo "  summary: $OUT/research_summary.md"
echo "  mode matrix: $OUT/modes/summary.tsv"
echo "  config: $OUT/research_config.tsv"

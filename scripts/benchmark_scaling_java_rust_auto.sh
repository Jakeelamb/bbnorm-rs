#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT"
export PATH="$HOME/.local/bin:$PATH"

R1=${1:-$ROOT/tmp/human_benchmark_8threads/human_GRCh38_500k_R1.fq.gz}
R2=${2:-$ROOT/tmp/human_benchmark_8threads/human_GRCh38_500k_R2.fq.gz}
OUT=${3:-$ROOT/tmp/java_rust_auto_scaling_$(date +%Y%m%d_%H%M%S)}
THREADS=${THREADS:-8}
ZIPTHREADS=${ZIPTHREADS:-$THREADS}
JAVA_XMX=${JAVA_XMX:-4g}
MEM=${MEM:-768m}
TIMEOUT=${TIMEOUT:-20m}
WRITE_OUTPUTS=${WRITE_OUTPUTS:-0}
JAVA_MAX_RSS_KB=${JAVA_MAX_RSS_KB:-8000000}
RUST_MAX_RSS_KB=${RUST_MAX_RSS_KB:-2500000}
EXTRA_ARGS=${EXTRA_ARGS:-}
READ_POINTS=${READ_POINTS:-"1000 5000 10000 50000 100000 250000 500000"}
HARNESS=${HARNESS:-scripts/benchmark_java_rust_human.sh}
PYTHON=${PYTHON:-python3}

mkdir -p "$OUT/runs"

if [[ ! -r "$R1" ]]; then
  echo "Missing readable R1 input: $R1" >&2
  exit 2
fi
if [[ ! -r "$R2" ]]; then
  echo "Missing readable R2 input: $R2" >&2
  exit 2
fi
if [[ ! -x "$HARNESS" && ! -f "$HARNESS" ]]; then
  echo "Missing harness: $HARNESS" >&2
  exit 2
fi

printf 'reads\ttool\telapsed_seconds\tmax_rss_kb\tmax_rss_gib\tstatus\tspeed_vs_java\tmemory_vs_java\tartifact_dir\n' > "$OUT/scaling.tsv"

commit="$(git rev-parse HEAD 2>/dev/null || true)"
branch="$(git rev-parse --abbrev-ref HEAD 2>/dev/null || true)"
[[ -n "$commit" ]] || commit="unknown"
[[ -n "$branch" ]] || branch="unknown"
{
  printf 'key\tvalue\n'
  printf 'commit\t%s\n' "$commit"
  printf 'branch\t%s\n' "$branch"
  printf 'r1\t%s\n' "$R1"
  printf 'r2\t%s\n' "$R2"
  printf 'threads\t%s\n' "$THREADS"
  printf 'zipthreads\t%s\n' "$ZIPTHREADS"
  printf 'java_xmx\t%s\n' "$JAVA_XMX"
  printf 'mem\t%s\n' "$MEM"
  printf 'timeout\t%s\n' "$TIMEOUT"
  printf 'write_outputs\t%s\n' "$WRITE_OUTPUTS"
  printf 'java_max_rss_kb\t%s\n' "$JAVA_MAX_RSS_KB"
  printf 'rust_max_rss_kb\t%s\n' "$RUST_MAX_RSS_KB"
  printf 'extra_args\t%s\n' "$EXTRA_ARGS"
  printf 'read_points\t%s\n' "$READ_POINTS"
  printf 'harness\t%s\n' "$HARNESS"
} > "$OUT/environment.tsv"

for reads in $READ_POINTS; do
  run_out="$OUT/runs/reads_${reads}"
  echo "==> reads=$reads -> $run_out"
  READS="$reads" \
  TABLE_READS="$reads" \
  THREADS="$THREADS" \
  ZIPTHREADS="$ZIPTHREADS" \
  JAVA_XMX="$JAVA_XMX" \
  MEM="$MEM" \
  TIMEOUT="$TIMEOUT" \
  WRITE_OUTPUTS="$WRITE_OUTPUTS" \
  JAVA_MAX_RSS_KB="$JAVA_MAX_RSS_KB" \
  RUST_MAX_RSS_KB="$RUST_MAX_RSS_KB" \
  EXTRA_ARGS="$EXTRA_ARGS" \
  "$HARNESS" "$R1" "$R2" "$run_out"

  "$PYTHON" - "$reads" "$run_out/results.tsv" "$run_out" >> "$OUT/scaling.tsv" <<'PY'
import csv, sys
reads = sys.argv[1]
results_path = sys.argv[2]
artifact_dir = sys.argv[3]
with open(results_path, newline='') as fh:
    rows = list(csv.DictReader(fh, delimiter='\t'))
by_tool = {row['tool']: row for row in rows}
java_elapsed = float(by_tool['java']['elapsed_seconds'])
java_rss = float(by_tool['java']['max_rss_kb'])
for tool in ('java', 'rust'):
    row = by_tool[tool]
    elapsed = float(row['elapsed_seconds'])
    rss_kb = float(row['max_rss_kb'])
    speed_vs_java = elapsed / java_elapsed if java_elapsed else ''
    memory_vs_java = rss_kb / java_rss if java_rss else ''
    print('\t'.join([
        reads,
        tool,
        row['elapsed_seconds'],
        row['max_rss_kb'],
        f"{rss_kb / 1024 / 1024:.6f}",
        row['status'],
        f"{speed_vs_java:.6f}" if speed_vs_java != '' else '',
        f"{memory_vs_java:.6f}" if memory_vs_java != '' else '',
        artifact_dir,
    ]))
PY

done

"$PYTHON" scripts/plot_scaling_java_rust_auto_svg.py "$OUT/scaling.tsv" "$OUT"

echo
printf 'Scaling sweep complete. Raw data: %s\n' "$OUT/scaling.tsv"
printf 'Speed graph: %s\n' "$OUT/speed_vs_reads.svg"
printf 'Memory graph: %s\n' "$OUT/memory_vs_reads.svg"
printf 'Summary table: %s\n' "$OUT/summary.tsv"

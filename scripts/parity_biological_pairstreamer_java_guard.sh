#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT"

BIO_ROOT="${BIO_ROOT:-/home/jake/Projects/biological data}"
DATA1="${DATA1:-$BIO_ROOT/reads/short_reads/ecoli_mg1655_pe_drr023054/DRR023054_1.fastq.gz}"
DATA2="${DATA2:-$BIO_ROOT/reads/short_reads/ecoli_mg1655_pe_drr023054/DRR023054_2.fastq.gz}"
CP="vendor/BBTools-master/current"
OUT="${1:-tmp/biological_pairstreamer_java_guard}"
READS="${READS:-1000}"
TABLE_READS="${TABLE_READS:-$READS}"
PYTHON="${PYTHON:-python}"
JAVA_TIMEOUT_SECONDS="${JAVA_TIMEOUT_SECONDS:-30}"

rm -rf "$OUT"
mkdir -p "$OUT"
printf 'mode\tstatus\toutput\n' > "$OUT/summary.tsv"

if [[ ! -f "$DATA1" || ! -f "$DATA2" ]]; then
  printf 'Missing paired dataset files:\n  DATA1=%s\n  DATA2=%s\n' "$DATA1" "$DATA2" >&2
  exit 2
fi

COMMON=(
  "in=$DATA1"
  "in2=$DATA2"
  "reads=$READS"
  "tablereads=$TABLE_READS"
  "keepall=t"
  "lowbindepth=1"
  "highbindepth=1"
  "k=31"
  "minq=0"
  "minprob=0"
  "min=0"
  "minkmers=1"
  "target=999999999"
  "max=999999999"
  "threads=1"
  "overwrite=t"
  "bits=32"
)

measure() {
  local metrics_file="$1"
  shift
  "$PYTHON" - "$metrics_file" "$@" <<'PY'
import resource
import subprocess
import sys
import time

metrics_file = sys.argv[1]
cmd = sys.argv[2:]
start = time.perf_counter()
proc = subprocess.run(cmd)
elapsed = time.perf_counter() - start
usage = resource.getrusage(resource.RUSAGE_CHILDREN)
with open(metrics_file, "w", encoding="utf-8") as handle:
    handle.write(f"{elapsed:.6f}\t{usage.ru_maxrss}\t{proc.returncode}\n")
sys.exit(proc.returncode)
PY
}

run_case() {
  local mode_name="$1"
  shift
  local extra_args=("$@")
  local case_out="$OUT/$mode_name"

  mkdir -p "$case_out"
  printf '\nRunning biological Java guard for %s on %s and %s...\n' "$mode_name" "$DATA1" "$DATA2"

  set +e
  measure "$case_out/java.metrics.tsv" \
    timeout "$JAVA_TIMEOUT_SECONDS" \
    java -Xmx1g -cp "$CP" jgi.KmerNormalize \
    "${COMMON[@]}" \
    "${extra_args[@]}" \
    "out=$case_out/java.keep1.fq" "out2=$case_out/java.keep2.fq" \
    "outlow=$case_out/java.low1.fq" "outlow2=$case_out/java.low2.fq" \
    "outmid=$case_out/java.mid1.fq" "outmid2=$case_out/java.mid2.fq" \
    "outhigh=$case_out/java.high1.fq" "outhigh2=$case_out/java.high2.fq" \
    "hist=$case_out/java.hist.tsv" "rhist=$case_out/java.rhist.tsv" \
    >"$case_out/java.stdout.log" 2>"$case_out/java.stderr.log"
  java_status=$?
  set -e

  if [[ "$java_status" -eq 0 ]]; then
    echo "Expected vendored Java to fail for $mode_name on this dataset, but it succeeded." >&2
    exit 1
  fi

  grep -q 'IndexOutOfBoundsException' "$case_out/java.stderr.log"
  grep -q 'PairStreamer.nextList' "$case_out/java.stderr.log"

  cargo build --release --quiet
  measure "$case_out/rust.metrics.tsv" \
    target/release/bbnorm-rs \
    "${COMMON[@]}" \
    "${extra_args[@]}" \
    "out=$case_out/rust.keep1.fq" "out2=$case_out/rust.keep2.fq" \
    "outlow=$case_out/rust.low1.fq" "outlow2=$case_out/rust.low2.fq" \
    "outmid=$case_out/rust.mid1.fq" "outmid2=$case_out/rust.mid2.fq" \
    "outhigh=$case_out/rust.high1.fq" "outhigh2=$case_out/rust.high2.fq" \
    "hist=$case_out/rust.hist.tsv" "rhist=$case_out/rust.rhist.tsv" \
    >"$case_out/rust.stdout.log" 2>"$case_out/rust.stderr.log"

  test -s "$case_out/rust.hist.tsv"
  test -s "$case_out/rust.rhist.tsv"
  printf '%s\tguarded\t%s\n' "$mode_name" "$case_out" >> "$OUT/summary.tsv"
}

run_case default "passes=1"
run_case passes2 "passes=2"

printf '\nBiological PairStreamer Java guard summary:\n'
column -t -s $'\t' "$OUT/summary.tsv" || cat "$OUT/summary.tsv"
printf '\nBiological PairStreamer Java guard passed. Vendored Java fails or times out after the PairStreamer crash on this paired dataset while Rust completes. Outputs and logs: %s\n' "$OUT"

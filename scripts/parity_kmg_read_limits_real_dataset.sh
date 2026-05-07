#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT"

DATA1="vendor/BBTools-master/resources/sample1.fq.gz"
DATA2="vendor/BBTools-master/resources/sample2.fq.gz"
CP="vendor/BBTools-master/current"
OUT="${1:-tmp/real_kmg_read_limits_parity}"
mkdir -p "$OUT"
rm -f "$OUT"/*

COMMON=(
  "in=$DATA1"
  "in2=$DATA2"
  "passes=1"
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

cargo build --quiet

run_case() {
  local name="$1"
  local limit_arg="$2"

  printf 'Running Java BBNorm KMG read-limit case %s...\n' "$name"
  local java_start
  java_start=$(python - <<'PYTIME'
import time
print(time.perf_counter_ns())
PYTIME
)
  java -Xmx1g -cp "$CP" jgi.KmerNormalize \
    "${COMMON[@]}" \
    "$limit_arg" \
    "out=$OUT/java.$name.keep1.fq" \
    "out2=$OUT/java.$name.keep2.fq" \
    "hist=$OUT/java.$name.hist.tsv" \
    "rhist=$OUT/java.$name.rhist.tsv" \
    >"$OUT/java.$name.stdout.log" 2>"$OUT/java.$name.stderr.log"
  local java_end
  java_end=$(python - <<'PYTIME'
import time
print(time.perf_counter_ns())
PYTIME
)
  local java_ms=$(((java_end - java_start) / 1000000))

  printf 'Running Rust bbnorm-rs KMG read-limit case %s...\n' "$name"
  local rust_start
  rust_start=$(python - <<'PYTIME'
import time
print(time.perf_counter_ns())
PYTIME
)
  target/debug/bbnorm-rs \
    "${COMMON[@]}" \
    "$limit_arg" \
    "out=$OUT/rust.$name.keep1.fq" \
    "out2=$OUT/rust.$name.keep2.fq" \
    "hist=$OUT/rust.$name.hist.tsv" \
    "rhist=$OUT/rust.$name.rhist.tsv" \
    >"$OUT/rust.$name.stdout.log" 2>"$OUT/rust.$name.stderr.log"
  local rust_end
  rust_end=$(python - <<'PYTIME'
import time
print(time.perf_counter_ns())
PYTIME
)
  local rust_ms=$(((rust_end - rust_start) / 1000000))

  for suffix in keep1.fq keep2.fq hist.tsv rhist.tsv; do
    cmp "$OUT/java.$name.$suffix" "$OUT/rust.$name.$suffix"
  done

  printf 'Case %s passed. Java: %sms, Rust: %sms\n' "$name" "$java_ms" "$rust_ms"
}

run_case reads_0_01k "reads=0.01k"
run_case tablereads_0_01k "tablereads=0.01k"

printf 'KMG read-limit parity passed. Outputs and logs: %s\n' "$OUT"

#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT"

DATA1="vendor/BBTools-master/resources/sample1.fq.gz"
DATA2="vendor/BBTools-master/resources/sample2.fq.gz"
CP="vendor/BBTools-master/current"
OUT="${1:-tmp/real_threshold_clamps_parity}"
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
  shift
  local java_start java_end rust_start rust_end java_ms rust_ms

  printf 'Running Java BBNorm threshold-clamp case %s...\n' "$name"
  java_start=$(python - <<'PYTIME'
import time
print(time.perf_counter_ns())
PYTIME
)
  java -Xmx1g -cp "$CP" jgi.KmerNormalize \
    "${COMMON[@]}" "$@" \
    "out=$OUT/java.$name.keep1.fq" "out2=$OUT/java.$name.keep2.fq" \
    "hist=$OUT/java.$name.hist.tsv" "rhist=$OUT/java.$name.rhist.tsv" \
    "outlow=$OUT/java.$name.low1.fq" "outlow2=$OUT/java.$name.low2.fq" \
    "outmid=$OUT/java.$name.mid1.fq" "outmid2=$OUT/java.$name.mid2.fq" \
    "outhigh=$OUT/java.$name.high1.fq" "outhigh2=$OUT/java.$name.high2.fq" \
    >"$OUT/java.$name.stdout.log" 2>"$OUT/java.$name.stderr.log"
  java_end=$(python - <<'PYTIME'
import time
print(time.perf_counter_ns())
PYTIME
)
  java_ms=$(((java_end - java_start) / 1000000))

  printf 'Running Rust bbnorm-rs threshold-clamp case %s...\n' "$name"
  rust_start=$(python - <<'PYTIME'
import time
print(time.perf_counter_ns())
PYTIME
)
  target/debug/bbnorm-rs \
    "${COMMON[@]}" "$@" \
    "out=$OUT/rust.$name.keep1.fq" "out2=$OUT/rust.$name.keep2.fq" \
    "hist=$OUT/rust.$name.hist.tsv" "rhist=$OUT/rust.$name.rhist.tsv" \
    "outlow=$OUT/rust.$name.low1.fq" "outlow2=$OUT/rust.$name.low2.fq" \
    "outmid=$OUT/rust.$name.mid1.fq" "outmid2=$OUT/rust.$name.mid2.fq" \
    "outhigh=$OUT/rust.$name.high1.fq" "outhigh2=$OUT/rust.$name.high2.fq" \
    >"$OUT/rust.$name.stdout.log" 2>"$OUT/rust.$name.stderr.log"
  rust_end=$(python - <<'PYTIME'
import time
print(time.perf_counter_ns())
PYTIME
)
  rust_ms=$(((rust_end - rust_start) / 1000000))

  for suffix in keep1.fq keep2.fq hist.tsv rhist.tsv low1.fq low2.fq mid1.fq mid2.fq high1.fq high2.fq; do
    cmp "$OUT/java.$name.$suffix" "$OUT/rust.$name.$suffix"
  done

  printf 'Case %s passed. Java: %sms, Rust: %sms\n' "$name" "$java_ms" "$rust_ms"
}

run_case max_clamped_to_target "target=100" "max=50"
run_case zero_minkmers_clamped "minkmers=0"

printf 'Threshold clamp parity passed. Outputs and logs: %s\n' "$OUT"

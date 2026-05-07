#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT"

DATA1="vendor/BBTools-master/resources/sample1.fq.gz"
DATA2="vendor/BBTools-master/resources/sample2.fq.gz"
CP="vendor/BBTools-master/current"
OUT="${1:-tmp/real_stage_alias_parity}"
mkdir -p "$OUT"
rm -f "$OUT"/*

COMMON=(
  "in=$DATA1"
  "in2=$DATA2"
  "passes=1"
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

  printf 'Running Java BBNorm stage-alias case %s...\n' "$name"
  java_start=$(python - <<'PYTIME'
import time
print(time.perf_counter_ns())
PYTIME
)
  java -Xmx1g -cp "$CP" jgi.KmerNormalize \
    "${COMMON[@]}" "$@" \
    "out=$OUT/java.$name.keep1.fq" "out2=$OUT/java.$name.keep2.fq" \
    "outt=$OUT/java.$name.toss1.fq" "outt2=$OUT/java.$name.toss2.fq" \
    >"$OUT/java.$name.stdout.log" 2>"$OUT/java.$name.stderr.log"
  java_end=$(python - <<'PYTIME'
import time
print(time.perf_counter_ns())
PYTIME
)
  java_ms=$(((java_end - java_start) / 1000000))

  printf 'Running Rust bbnorm-rs stage-alias case %s...\n' "$name"
  rust_start=$(python - <<'PYTIME'
import time
print(time.perf_counter_ns())
PYTIME
)
  target/debug/bbnorm-rs \
    "${COMMON[@]}" "$@" \
    "out=$OUT/rust.$name.keep1.fq" "out2=$OUT/rust.$name.keep2.fq" \
    "outt=$OUT/rust.$name.toss1.fq" "outt2=$OUT/rust.$name.toss2.fq" \
    >"$OUT/rust.$name.stdout.log" 2>"$OUT/rust.$name.stderr.log"
  rust_end=$(python - <<'PYTIME'
import time
print(time.perf_counter_ns())
PYTIME
)
  rust_ms=$(((rust_end - rust_start) / 1000000))

  cmp "$OUT/java.$name.keep1.fq" "$OUT/rust.$name.keep1.fq"
  cmp "$OUT/java.$name.keep2.fq" "$OUT/rust.$name.keep2.fq"
  cmp "$OUT/java.$name.toss1.fq" "$OUT/rust.$name.toss1.fq"
  cmp "$OUT/java.$name.toss2.fq" "$OUT/rust.$name.toss2.fq"

  printf 'Stage-alias case %s parity passed. Java: %sms, Rust: %sms\n' "$name" "$java_ms" "$rust_ms"
}

run_case final_toss "tossbadreadsf=t"
run_case final_toss2_alias "tossbadreads2=t"
run_case final_ter2_alias "ter2=t"
run_case final_tbr2_alias "tbr2=t"
run_case first_toss_noop "tossbadreads1=t"

printf 'Stage-alias parity outputs and logs: %s\n' "$OUT"

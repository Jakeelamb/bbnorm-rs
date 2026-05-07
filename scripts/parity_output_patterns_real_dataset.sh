#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT"

DATA1="vendor/BBTools-master/resources/sample1.fq.gz"
DATA2="vendor/BBTools-master/resources/sample2.fq.gz"
CP="vendor/BBTools-master/current"
OUT="${1:-tmp/real_output_patterns_parity}"
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
  local java_args=("$@")
  local rust_args=()
  for arg in "${java_args[@]}"; do
    rust_args+=("${arg/java./rust.}")
  done

  printf 'Running Java BBNorm output-pattern case %s...\n' "$name"
  local java_start
  java_start=$(python - <<'PYTIME'
import time
print(time.perf_counter_ns())
PYTIME
)
  java -Xmx1g -cp "$CP" jgi.KmerNormalize \
    "${COMMON[@]}" \
    "${java_args[@]}" \
    >"$OUT/java.$name.stdout.log" 2>"$OUT/java.$name.stderr.log"
  local java_end
  java_end=$(python - <<'PYTIME'
import time
print(time.perf_counter_ns())
PYTIME
)
  local java_ms=$(((java_end - java_start) / 1000000))

  printf 'Running Rust bbnorm-rs output-pattern case %s...\n' "$name"
  local rust_start
  rust_start=$(python - <<'PYTIME'
import time
print(time.perf_counter_ns())
PYTIME
)
  target/debug/bbnorm-rs \
    "${COMMON[@]}" \
    "${rust_args[@]}" \
    >"$OUT/rust.$name.stdout.log" 2>"$OUT/rust.$name.stderr.log"
  local rust_end
  rust_end=$(python - <<'PYTIME'
import time
print(time.perf_counter_ns())
PYTIME
)
  local rust_ms=$(((rust_end - rust_start) / 1000000))

  printf 'Case %s ran. Java: %sms, Rust: %sms\n' "$name" "$java_ms" "$rust_ms"
}

run_case hash_split \
  "out=$OUT/java.hash.keep#.fq" \
  "outlow=$OUT/java.hash.low#.fq"
cmp "$OUT/java.hash.keep1.fq" "$OUT/rust.hash.keep1.fq"
cmp "$OUT/java.hash.keep2.fq" "$OUT/rust.hash.keep2.fq"
cmp "$OUT/java.hash.low1.fq" "$OUT/rust.hash.low1.fq"
cmp "$OUT/java.hash.low2.fq" "$OUT/rust.hash.low2.fq"

run_case implicit_interleaved \
  "out=$OUT/java.interleaved.keep.fq"
cmp "$OUT/java.interleaved.keep.fq" "$OUT/rust.interleaved.keep.fq"

printf 'Output-pattern parity passed. Outputs and logs: %s\n' "$OUT"

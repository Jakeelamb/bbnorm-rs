#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT"

DATA_HASH="vendor/BBTools-master/resources/sample#.fq.gz"
CP="vendor/BBTools-master/current"
OUT="${1:-tmp/real_hash_input_parity}"
mkdir -p "$OUT"
rm -f "$OUT"/*

COMMON=(
  "in=$DATA_HASH"
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

printf 'Running Java BBNorm paired # input expansion on real phiX fixture...\n'
JAVA_START=$(python - <<'PYTIME'
import time
print(time.perf_counter_ns())
PYTIME
)
java -Xmx1g -cp "$CP" jgi.KmerNormalize \
  "${COMMON[@]}" \
  "out=$OUT/java.keep1.fq" "out2=$OUT/java.keep2.fq" \
  "outt=$OUT/java.toss1.fq" "outt2=$OUT/java.toss2.fq" \
  >"$OUT/java.stdout.log" 2>"$OUT/java.stderr.log"
JAVA_END=$(python - <<'PYTIME'
import time
print(time.perf_counter_ns())
PYTIME
)
JAVA_MS=$(((JAVA_END - JAVA_START) / 1000000))

printf 'Running Rust bbnorm-rs paired # input expansion on same fixture...\n'
RUST_START=$(python - <<'PYTIME'
import time
print(time.perf_counter_ns())
PYTIME
)
target/debug/bbnorm-rs \
  "${COMMON[@]}" \
  "out=$OUT/rust.keep1.fq" "out2=$OUT/rust.keep2.fq" \
  "outt=$OUT/rust.toss1.fq" "outt2=$OUT/rust.toss2.fq" \
  >"$OUT/rust.stdout.log" 2>"$OUT/rust.stderr.log"
RUST_END=$(python - <<'PYTIME'
import time
print(time.perf_counter_ns())
PYTIME
)
RUST_MS=$(((RUST_END - RUST_START) / 1000000))

cmp "$OUT/java.keep1.fq" "$OUT/rust.keep1.fq"
cmp "$OUT/java.keep2.fq" "$OUT/rust.keep2.fq"
cmp "$OUT/java.toss1.fq" "$OUT/rust.toss1.fq"
cmp "$OUT/java.toss2.fq" "$OUT/rust.toss2.fq"

printf 'Paired # input expansion parity passed. Java: %sms, Rust: %sms\n' "$JAVA_MS" "$RUST_MS"
printf 'Outputs and logs: %s\n' "$OUT"

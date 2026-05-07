#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT"

DATA1="vendor/BBTools-master/resources/sample1.fq.gz"
DATA2="vendor/BBTools-master/resources/sample2.fq.gz"
CP="vendor/BBTools-master/current"
OUT="${1:-tmp/real_parity}"
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

printf 'Running Java BBNorm on real paired phiX fixture...\n'
JAVA_START=$(python - <<'PYTIME'
import time
print(time.perf_counter_ns())
PYTIME
)
java -Xmx1g -cp "$CP" jgi.KmerNormalize \
  "${COMMON[@]}" \
  "out=$OUT/java.keep1.fq" "out2=$OUT/java.keep2.fq" \
  "outlow=$OUT/java.low1.fq" "outlow2=$OUT/java.low2.fq" \
  "outmid=$OUT/java.mid1.fq" "outmid2=$OUT/java.mid2.fq" \
  "outhigh=$OUT/java.high1.fq" "outhigh2=$OUT/java.high2.fq" \
  "hist=$OUT/java.hist.tsv" "rhist=$OUT/java.rhist.tsv" \
  >"$OUT/java.stdout.log" 2>"$OUT/java.stderr.log"
JAVA_END=$(python - <<'PYTIME'
import time
print(time.perf_counter_ns())
PYTIME
)
JAVA_MS=$(((JAVA_END - JAVA_START) / 1000000))

printf 'Running Rust bbnorm-rs on same fixture...\n'
RUST_START=$(python - <<'PYTIME'
import time
print(time.perf_counter_ns())
PYTIME
)
target/debug/bbnorm-rs \
  "${COMMON[@]}" \
  "out=$OUT/rust.keep1.fq" "out2=$OUT/rust.keep2.fq" \
  "outlow=$OUT/rust.low1.fq" "outlow2=$OUT/rust.low2.fq" \
  "outmid=$OUT/rust.mid1.fq" "outmid2=$OUT/rust.mid2.fq" \
  "outhigh=$OUT/rust.high1.fq" "outhigh2=$OUT/rust.high2.fq" \
  "hist=$OUT/rust.hist.tsv" "rhist=$OUT/rust.rhist.tsv" \
  >"$OUT/rust.stdout.log" 2>"$OUT/rust.stderr.log"
RUST_END=$(python - <<'PYTIME'
import time
print(time.perf_counter_ns())
PYTIME
)
RUST_MS=$(((RUST_END - RUST_START) / 1000000))

cmp "$OUT/java.keep1.fq" "$OUT/rust.keep1.fq"
cmp "$OUT/java.keep2.fq" "$OUT/rust.keep2.fq"
cmp "$OUT/java.hist.tsv" "$OUT/rust.hist.tsv"
cmp "$OUT/java.rhist.tsv" "$OUT/rust.rhist.tsv"
cmp "$OUT/java.low1.fq" "$OUT/rust.low1.fq"
cmp "$OUT/java.low2.fq" "$OUT/rust.low2.fq"
cmp "$OUT/java.mid1.fq" "$OUT/rust.mid1.fq"
cmp "$OUT/java.mid2.fq" "$OUT/rust.mid2.fq"
cmp "$OUT/java.high1.fq" "$OUT/rust.high1.fq"
cmp "$OUT/java.high2.fq" "$OUT/rust.high2.fq"

printf 'Parity passed. Java: %sms, Rust: %sms\n' "$JAVA_MS" "$RUST_MS"
printf 'Outputs and logs: %s\n' "$OUT"

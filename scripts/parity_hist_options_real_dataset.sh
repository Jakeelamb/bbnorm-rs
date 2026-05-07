#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT"

DATA1="vendor/BBTools-master/resources/sample1.fq.gz"
DATA2="vendor/BBTools-master/resources/sample2.fq.gz"
CP="vendor/BBTools-master/current"
OUT="${1:-tmp/real_hist_options_parity}"
mkdir -p "$OUT"
rm -f "$OUT"/*

COMMON=(
  "in=$DATA1"
  "in2=$DATA2"
  "passes=1"
  "keepall=t"
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
  local extra=("$@")

  printf 'Running Java BBNorm histogram option case %s...\n' "$name"
  local java_start
  java_start=$(python - <<'PYTIME'
import time
print(time.perf_counter_ns())
PYTIME
)
  java -Xmx1g -cp "$CP" jgi.KmerNormalize \
    "${COMMON[@]}" \
    "${extra[@]}" \
    "out=$OUT/java.$name.keep1.fq" \
    "out2=$OUT/java.$name.keep2.fq" \
    "hist=$OUT/java.$name.hist.tsv" \
    >"$OUT/java.$name.stdout.log" 2>"$OUT/java.$name.stderr.log"
  local java_end
  java_end=$(python - <<'PYTIME'
import time
print(time.perf_counter_ns())
PYTIME
)
  local java_ms=$(((java_end - java_start) / 1000000))

  printf 'Running Rust bbnorm-rs histogram option case %s...\n' "$name"
  local rust_start
  rust_start=$(python - <<'PYTIME'
import time
print(time.perf_counter_ns())
PYTIME
)
  target/debug/bbnorm-rs \
    "${COMMON[@]}" \
    "${extra[@]}" \
    "out=$OUT/rust.$name.keep1.fq" \
    "out2=$OUT/rust.$name.keep2.fq" \
    "hist=$OUT/rust.$name.hist.tsv" \
    >"$OUT/rust.$name.stdout.log" 2>"$OUT/rust.$name.stderr.log"
  local rust_end
  rust_end=$(python - <<'PYTIME'
import time
print(time.perf_counter_ns())
PYTIME
)
  local rust_ms=$(((rust_end - rust_start) / 1000000))

  cmp "$OUT/java.$name.keep1.fq" "$OUT/rust.$name.keep1.fq"
  cmp "$OUT/java.$name.keep2.fq" "$OUT/rust.$name.keep2.fq"
  cmp "$OUT/java.$name.hist.tsv" "$OUT/rust.$name.hist.tsv"

  printf 'Case %s passed. Java: %sms, Rust: %sms\n' "$name" "$java_ms" "$rust_ms"
}

run_case histcol1 "histcol=1"
run_case histcol2 "histcol=2"
run_case zerobin "zerobin=t"
run_case print_zero_coverage "printzerocoverage=t" "histlen=5"

printf 'Histogram option parity passed. Outputs and logs: %s\n' "$OUT"

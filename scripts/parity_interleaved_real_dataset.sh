#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT"

DATA1="vendor/BBTools-master/resources/sample1.fq.gz"
DATA2="vendor/BBTools-master/resources/sample2.fq.gz"
CP="vendor/BBTools-master/current"
OUT="${1:-tmp/real_interleaved_parity}"
mkdir -p "$OUT"
rm -f "$OUT"/*

INTERLEAVED="$OUT/interleaved.fq"
python - "$DATA1" "$DATA2" "$INTERLEAVED" <<'PYBUILD'
import gzip
import sys

r1_path, r2_path, out_path = sys.argv[1:]
with gzip.open(r1_path, 'rt') as r1, gzip.open(r2_path, 'rt') as r2, open(out_path, 'w') as out:
    while True:
        rec1 = [r1.readline() for _ in range(4)]
        rec2 = [r2.readline() for _ in range(4)]
        if not rec1[0] and not rec2[0]:
            break
        if not all(rec1) or not all(rec2):
            raise SystemExit('paired phiX fixture is truncated or has uneven record counts')
        out.writelines(rec1)
        out.writelines(rec2)
PYBUILD

COMMON=(
  "in=$INTERLEAVED"
  "interleaved=t"
  "passes=1"
  "k=31"
  "minq=0"
  "minprob=0"
  "min=0"
  "minkmers=1"
  "target=1"
  "max=1"
  "threads=1"
  "overwrite=t"
  "bits=32"
)

cargo build --quiet

printf 'Running Java BBNorm normalization on interleaved real paired phiX fixture...\n'
JAVA_START=$(python - <<'PYTIME'
import time
print(time.perf_counter_ns())
PYTIME
)
java -Xmx1g -cp "$CP" jgi.KmerNormalize \
  "${COMMON[@]}" \
  "out=$OUT/java.keep.fq" \
  "outt=$OUT/java.toss.fq" \
  >"$OUT/java.stdout.log" 2>"$OUT/java.stderr.log"
JAVA_END=$(python - <<'PYTIME'
import time
print(time.perf_counter_ns())
PYTIME
)
JAVA_MS=$(((JAVA_END - JAVA_START) / 1000000))

printf 'Running Rust bbnorm-rs normalization on same interleaved fixture...\n'
RUST_START=$(python - <<'PYTIME'
import time
print(time.perf_counter_ns())
PYTIME
)
target/debug/bbnorm-rs \
  "${COMMON[@]}" \
  "out=$OUT/rust.keep.fq" \
  "outt=$OUT/rust.toss.fq" \
  >"$OUT/rust.stdout.log" 2>"$OUT/rust.stderr.log"
RUST_END=$(python - <<'PYTIME'
import time
print(time.perf_counter_ns())
PYTIME
)
RUST_MS=$(((RUST_END - RUST_START) / 1000000))

cmp "$OUT/java.keep.fq" "$OUT/rust.keep.fq"
cmp "$OUT/java.toss.fq" "$OUT/rust.toss.fq"

printf 'Interleaved normalization parity passed. Java: %sms, Rust: %sms\n' "$JAVA_MS" "$RUST_MS"
printf 'Outputs and logs: %s\n' "$OUT"

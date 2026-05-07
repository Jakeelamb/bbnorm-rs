#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT"

DATA1="vendor/BBTools-master/resources/sample1.fq.gz"
CP="vendor/BBTools-master/current"
OUT="${1:-tmp/real_qtrim_parity}"
mkdir -p "$OUT"
rm -f "$OUT"/*

FIXTURE="$OUT/real_qtrim_tail.fq"
python - "$DATA1" "$FIXTURE" <<'PYFIXTURE'
import gzip
import sys

source, dest = sys.argv[1:3]
low_tail = 12
chosen = None
with gzip.open(source, "rt") as reader:
    while True:
        record = [reader.readline(), reader.readline(), reader.readline(), reader.readline()]
        if not record[0]:
            break
        bases = record[1].strip()
        if "N" not in bases and "n" not in bases and len(bases) > low_tail + 31:
            chosen = record
            break

if chosen is None:
    raise SystemExit(f"no long N-free read found in {source}")

bases = chosen[1].strip()
qualities = "I" * (len(bases) - low_tail) + "!" * low_tail + "\n"
with open(dest, "w") as writer:
    for copy in range(4):
        writer.write(f"@real_qtrim_{copy}\n")
        writer.write(chosen[1])
        writer.write("+\n")
        writer.write(qualities)
PYFIXTURE

COMMON=(
  "in=$FIXTURE"
  "passes=1"
  "keepall=t"
  "qtrim=r"
  "trimq=10"
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

printf 'Running Java BBNorm qtrim real-derived keep-all case...\n'
JAVA_START=$(python - <<'PYTIME'
import time
print(time.perf_counter_ns())
PYTIME
)
java -Xmx1g -cp "$CP" jgi.KmerNormalize \
  "${COMMON[@]}" \
  "out=$OUT/java.keep.fq" \
  "hist=$OUT/java.hist.tsv" \
  "rhist=$OUT/java.rhist.tsv" \
  >"$OUT/java.stdout.log" 2>"$OUT/java.stderr.log"
JAVA_END=$(python - <<'PYTIME'
import time
print(time.perf_counter_ns())
PYTIME
)
JAVA_MS=$(((JAVA_END - JAVA_START) / 1000000))

printf 'Running Rust bbnorm-rs qtrim real-derived keep-all case...\n'
RUST_START=$(python - <<'PYTIME'
import time
print(time.perf_counter_ns())
PYTIME
)
target/debug/bbnorm-rs \
  "${COMMON[@]}" \
  "out=$OUT/rust.keep.fq" \
  "hist=$OUT/rust.hist.tsv" \
  "rhist=$OUT/rust.rhist.tsv" \
  >"$OUT/rust.stdout.log" 2>"$OUT/rust.stderr.log"
RUST_END=$(python - <<'PYTIME'
import time
print(time.perf_counter_ns())
PYTIME
)
RUST_MS=$(((RUST_END - RUST_START) / 1000000))

cmp "$OUT/java.keep.fq" "$OUT/rust.keep.fq"
cmp "$OUT/java.hist.tsv" "$OUT/rust.hist.tsv"
cmp "$OUT/java.rhist.tsv" "$OUT/rust.rhist.tsv"

printf 'qtrim parity passed. Java: %sms, Rust: %sms\n' "$JAVA_MS" "$RUST_MS"
printf 'Outputs and logs: %s\n' "$OUT"

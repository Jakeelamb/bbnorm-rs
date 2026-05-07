#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT"

CP="vendor/BBTools-master/current"
OUT="${1:-tmp/representative_long_kmer_fixspikes_parity}"
mkdir -p "$OUT"
rm -f "$OUT"/*

FIXTURE="$OUT/long_kmer_spike.fq"
python - "$FIXTURE" <<'PYFIXTURE'
import sys

dest = sys.argv[1]
center = "ACGTTGCAAGTCGATCGTAGCTAGGATCCGATGCTAGTCA"
assert len(center) == 40
records = [("target", "A" + center + "C")]
records.extend((f"center_dup_{idx}", center) for idx in range(3))
with open(dest, "w") as writer:
    for name, seq in records:
        writer.write(f"@{name}\n{seq}\n+\n{'I' * len(seq)}\n")
PYFIXTURE

COMMON=(
  "in=$FIXTURE"
  "passes=1"
  "keepall=t"
  "k=40"
  "fixspikes=t"
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

printf 'Running Java BBNorm long-kmer fixspikes representative spike case...\n'
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

printf 'Running Rust bbnorm-rs long-kmer fixspikes representative spike case...\n'
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

printf 'Long-kmer fixspikes parity passed. Java: %sms, Rust: %sms\n' "$JAVA_MS" "$RUST_MS"
printf 'Outputs and logs: %s\n' "$OUT"

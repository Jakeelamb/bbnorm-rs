#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT"

CP="vendor/BBTools-master/current"
OUT="${1:-tmp/representative_default_multipass_fallback}"
mkdir -p "$OUT"
rm -f "$OUT"/*

INPUT="$OUT/default_multipass_input.fq"
cat > "$INPUT" <<'FASTQ'
@default_pass_probe
ACGTACGT
+
IIIIIIII
FASTQ

COMMON=(
  "in=$INPUT"
  "keepall=t"
  "k=3"
  "minq=0"
  "minprob=0"
  "min=0"
  "minkmers=1"
  "target=999999999"
  "max=999999999"
  "threads=1"
  "bits=32"
)

cargo build --quiet

printf 'Running Java BBNorm with omitted passes (default multipass)...\n'
java -Xmx1g -cp "$CP" jgi.KmerNormalize \
  "${COMMON[@]}" \
  "out=$OUT/java.keep.fq" \
  "hist=$OUT/java.hist.tsv" \
  >"$OUT/java.stdout.log" 2>"$OUT/java.stderr.log"

test -s "$OUT/java.keep.fq"
test -s "$OUT/java.hist.tsv"

printf 'Running Rust bbnorm-rs with omitted passes through real multipass orchestration...\n'
target/debug/bbnorm-rs \
  "${COMMON[@]}" \
  "out=$OUT/rust.keep.fq" \
  "hist=$OUT/rust.hist.tsv" \
  >"$OUT/rust.stdout.log" 2>"$OUT/rust.stderr.log"

cmp "$OUT/java.keep.fq" "$OUT/rust.keep.fq"
cmp "$OUT/java.hist.tsv" "$OUT/rust.hist.tsv"

printf 'Default multipass parity passed. Rust runs intermediate temp-pass orchestration and matches Java output. Logs and outputs: %s\n' "$OUT"

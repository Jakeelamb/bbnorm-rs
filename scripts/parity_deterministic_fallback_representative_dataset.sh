#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT"

CP="vendor/BBTools-master/current"
OUT="${1:-tmp/representative_deterministic_mode}"
rm -rf "$OUT"
mkdir -p "$OUT"

INPUT="$OUT/deterministic_input.fq"
cat > "$INPUT" <<'FASTQ'
@deterministic_probe_1
ACGTACGT
+
IIIIIIII
@deterministic_probe_2
ACGTACGT
+
IIIIIIII
FASTQ

COMMON=(
  "in=$INPUT"
  "passes=1"
  "keepall=t"
  "k=3"
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

printf 'Running Java deterministic baseline...\n'
java -Xmx1g -cp "$CP" jgi.KmerNormalize \
  "${COMMON[@]}" \
  "out=$OUT/java.keep.fq" \
  "hist=$OUT/java.hist.tsv" \
  >"$OUT/java.stdout.log" 2>"$OUT/java.stderr.log"

printf 'Running Rust deterministic=f supported mode...\n'
target/debug/bbnorm-rs \
  "${COMMON[@]}" \
  "deterministic=f" \
  "out=$OUT/rust.keep.fq" \
  "hist=$OUT/rust.hist.tsv" \
  >"$OUT/rust.stdout.log" 2>"$OUT/rust.stderr.log"

if grep -q 'deterministic=f is not implemented yet' "$OUT/rust.stderr.log"; then
  echo "deterministic=f still emitted the old fallback note" >&2
  exit 1
fi
cmp "$OUT/java.keep.fq" "$OUT/rust.keep.fq"
cmp "$OUT/java.hist.tsv" "$OUT/rust.hist.tsv"

printf 'deterministic=f supported-mode smoke passed. The keepall fixture remains Java-identical while Rust keeps nondeterministic read selection enabled for non-keepall normalization. Logs and outputs: %s\n' "$OUT"

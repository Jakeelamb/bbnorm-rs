#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT"

CP="vendor/BBTools-master/current"
OUT="${1:-tmp/representative_append_parity}"
mkdir -p "$OUT"
rm -f "$OUT"/*

INPUT="$OUT/append_input.fq"
cat > "$INPUT" <<'FASTQ'
@append_probe
ACGTACGT
+
IIIIIIII
FASTQ

printf 'PREFASTQ\n' > "$OUT/java.keep.fq"
printf 'PREFASTQ\n' > "$OUT/rust.keep.fq"
printf 'PREHIST\n' > "$OUT/java.hist.tsv"
printf 'PREHIST\n' > "$OUT/rust.hist.tsv"

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
  "append=t"
  "bits=32"
)

cargo build --quiet

printf 'Running Java BBNorm append representative case...\n'
java -Xmx1g -cp "$CP" jgi.KmerNormalize \
  "${COMMON[@]}" \
  "out=$OUT/java.keep.fq" \
  "hist=$OUT/java.hist.tsv" \
  >"$OUT/java.stdout.log" 2>"$OUT/java.stderr.log"

printf 'Running Rust bbnorm-rs append representative case...\n'
target/debug/bbnorm-rs \
  "${COMMON[@]}" \
  "out=$OUT/rust.keep.fq" \
  "hist=$OUT/rust.hist.tsv" \
  >"$OUT/rust.stdout.log" 2>"$OUT/rust.stderr.log"

cmp "$OUT/java.keep.fq" "$OUT/rust.keep.fq"
cmp "$OUT/java.hist.tsv" "$OUT/rust.hist.tsv"
grep -q '^PREFASTQ$' "$OUT/rust.keep.fq"
if grep -q '^PREHIST$' "$OUT/rust.hist.tsv"; then
  echo 'Expected append=t to append read output but overwrite hist output like vendored BBNorm.' >&2
  exit 1
fi

printf 'Append representative parity passed. Logs and outputs: %s\n' "$OUT"

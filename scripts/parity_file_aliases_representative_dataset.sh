#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT"

OUT="${1:-tmp/representative_file_aliases}"
mkdir -p "$OUT"
rm -f "$OUT"/*
CP="vendor/BBTools-master/current"
INPUT1="$OUT/alias_r1.fq"
INPUT2="$OUT/alias_r2.fq"

cat > "$INPUT1" <<'FASTQ'
@r1/1
ACGTACGTACGTACGTACGTACGTACGTACGT
+
IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
@r2/1
TGCATGCATGCATGCATGCATGCATGCATGCA
+
IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
FASTQ
cat > "$INPUT2" <<'FASTQ'
@r1/2
ACGTACGTACGTACGTACGTACGTACGTACGT
+
IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
@r2/2
TGCATGCATGCATGCATGCATGCATGCATGCA
+
IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
FASTQ

COMMON=(
  "passes=1"
  "keepall=t"
  "k=15"
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

printf 'Running Java BBNorm baseline with canonical in/out names...\n'
java -Xmx1g -cp "$CP" jgi.KmerNormalize \
  "in=$INPUT1" "in2=$INPUT2" "${COMMON[@]}" \
  "out=$OUT/java.keep1.fq" "out2=$OUT/java.keep2.fq" \
  "hist=$OUT/java.hist.tsv" \
  >"$OUT/java.stdout.log" 2>"$OUT/java.stderr.log"

printf 'Confirming vendored KmerNormalize rejects shared input/output aliases...\n'
if java -Xmx1g -cp "$CP" jgi.KmerNormalize \
  "input=$INPUT1" "input2=$INPUT2" "${COMMON[@]}" \
  "output=$OUT/java.reject.keep1.fq" "output2=$OUT/java.reject.keep2.fq" \
  >"$OUT/java.reject.stdout.log" 2>"$OUT/java.reject.stderr.log"; then
  printf 'Expected vendored Java to reject input/output aliases, but it succeeded.\n' >&2
  exit 1
fi
grep -q 'Unknown parameter input=' "$OUT/java.reject.stderr.log"

printf 'Running Rust bbnorm-rs with shared input/output aliases...\n'
target/debug/bbnorm-rs \
  "input=$INPUT1" "input2=$INPUT2" "${COMMON[@]}" \
  "output=$OUT/rust.keep1.fq" "output2=$OUT/rust.keep2.fq" \
  "hist=$OUT/rust.hist.tsv" \
  >"$OUT/rust.stdout.log" 2>"$OUT/rust.stderr.log"

for suffix in keep1.fq keep2.fq hist.tsv; do
  cmp "$OUT/java.$suffix" "$OUT/rust.$suffix"
done

printf 'Shared file-alias fallback passed. Vendored KmerNormalize rejects input/output aliases; Rust accepts them and matches canonical Java output. Logs: %s\n' "$OUT"

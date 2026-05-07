#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT"

CP="vendor/BBTools-master/current"
OUT="${1:-tmp/representative_extra_parity}"
mkdir -p "$OUT"
rm -f "$OUT"/*

MAIN="$OUT/representative_extra_main.fq"
EXTRA1="$OUT/representative_extra_1.fq"
EXTRA2="$OUT/representative_extra_2.fq"
EXTRA_LITERAL="$OUT/representative,extra,literal.fq"
EXTRA_HASH="$OUT/representative_extra_hash#.fq"
EXTRA_HASH1="$OUT/representative_extra_hash1.fq"
EXTRA_HASH2="$OUT/representative_extra_hash2.fq"
cat > "$MAIN" <<'FASTQ'
@main
ACGTACGT
+
IIIIIIII
FASTQ
cat > "$EXTRA1" <<'FASTQ'
@extra1
ACGTACGT
+
IIIIIIII
@extra2
ACGTACGT
+
IIIIIIII
FASTQ
cat > "$EXTRA2" <<'FASTQ'
@extra3
TTTTTTTT
+
IIIIIIII
FASTQ
cat > "$EXTRA_LITERAL" <<'FASTQ'
@literal_extra1
ACGTACGT
+
IIIIIIII
@literal_extra2
ACGTACGT
+
IIIIIIII
FASTQ
cat > "$EXTRA_HASH1" <<'FASTQ'
@hash_extra1
ACGTACGT
+
IIIIIIII
FASTQ
cat > "$EXTRA_HASH2" <<'FASTQ'
@hash_extra2
ACGTACGT
+
IIIIIIII
FASTQ

COMMON=(
  "in=$MAIN"
  "extra=$EXTRA1,$EXTRA2"
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

printf 'Running Java BBNorm extra-input representative case...\n'
java -Xmx1g -cp "$CP" jgi.KmerNormalize \
  "${COMMON[@]}" \
  "out=$OUT/java.keep.fq" \
  "hist=$OUT/java.hist.tsv" \
  >"$OUT/java.stdout.log" 2>"$OUT/java.stderr.log"

printf 'Running Rust bbnorm-rs extra-input representative case...\n'
target/debug/bbnorm-rs \
  "${COMMON[@]}" \
  "out=$OUT/rust.keep.fq" \
  "hist=$OUT/rust.hist.tsv" \
  >"$OUT/rust.stdout.log" 2>"$OUT/rust.stderr.log"

cmp "$OUT/java.keep.fq" "$OUT/rust.keep.fq"
cmp "$OUT/java.hist.tsv" "$OUT/rust.hist.tsv"

LITERAL_COMMON=(
  "in=$MAIN"
  "extra=$EXTRA_LITERAL"
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

printf 'Running Java BBNorm literal-comma extra-input case...\n'
java -Xmx1g -cp "$CP" jgi.KmerNormalize \
  "${LITERAL_COMMON[@]}" \
  "out=$OUT/java.literal_comma.keep.fq" \
  "hist=$OUT/java.literal_comma.hist.tsv" \
  >"$OUT/java.literal_comma.stdout.log" 2>"$OUT/java.literal_comma.stderr.log"

printf 'Running Rust bbnorm-rs literal-comma extra-input case...\n'
target/debug/bbnorm-rs \
  "${LITERAL_COMMON[@]}" \
  "out=$OUT/rust.literal_comma.keep.fq" \
  "hist=$OUT/rust.literal_comma.hist.tsv" \
  >"$OUT/rust.literal_comma.stdout.log" 2>"$OUT/rust.literal_comma.stderr.log"

cmp "$OUT/java.literal_comma.keep.fq" "$OUT/rust.literal_comma.keep.fq"
cmp "$OUT/java.literal_comma.hist.tsv" "$OUT/rust.literal_comma.hist.tsv"

TABLEREADS_COMMON=(
  "in=$MAIN"
  "extra=$EXTRA1"
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
  "tablereads=1"
)

printf 'Running Java BBNorm extra-input tablereads-limit case...\n'
java -Xmx1g -cp "$CP" jgi.KmerNormalize \
  "${TABLEREADS_COMMON[@]}" \
  "out=$OUT/java.tablereads.keep.fq" \
  "hist=$OUT/java.tablereads.hist.tsv" \
  >"$OUT/java.tablereads.stdout.log" 2>"$OUT/java.tablereads.stderr.log"

printf 'Running Rust bbnorm-rs extra-input tablereads-limit case...\n'
target/debug/bbnorm-rs \
  "${TABLEREADS_COMMON[@]}" \
  "out=$OUT/rust.tablereads.keep.fq" \
  "hist=$OUT/rust.tablereads.hist.tsv" \
  >"$OUT/rust.tablereads.stdout.log" 2>"$OUT/rust.tablereads.stderr.log"

cmp "$OUT/java.tablereads.keep.fq" "$OUT/rust.tablereads.keep.fq"
cmp "$OUT/java.tablereads.hist.tsv" "$OUT/rust.tablereads.hist.tsv"

HASH_COMMON=(
  "in=$MAIN"
  "extra=$EXTRA_HASH"
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

printf 'Checking Java BBNorm rejection for missing literal hash-pattern extra...\n'
if java -Xmx1g -cp "$CP" jgi.KmerNormalize \
  "${HASH_COMMON[@]}" \
  "out=$OUT/java.hash_extra.keep.fq" \
  "hist=$OUT/java.hash_extra.hist.tsv" \
  >"$OUT/java.hash_extra.stdout.log" 2>"$OUT/java.hash_extra.stderr.log"; then
  printf 'Expected Java BBNorm to reject missing literal hash-pattern extra input.\n' >&2
  exit 1
fi

printf 'Checking Rust bbnorm-rs rejection for missing literal hash-pattern extra...\n'
if target/debug/bbnorm-rs \
  "${HASH_COMMON[@]}" \
  "out=$OUT/rust.hash_extra.keep.fq" \
  "hist=$OUT/rust.hash_extra.hist.tsv" \
  >"$OUT/rust.hash_extra.stdout.log" 2>"$OUT/rust.hash_extra.stderr.log"; then
  printf 'Expected Rust bbnorm-rs to reject missing literal hash-pattern extra input.\n' >&2
  exit 1
fi

printf 'Extra-input representative parity passed, including literal comma filenames, tablereads handling, and hash-pattern extra rejection. Logs and outputs: %s\n' "$OUT"

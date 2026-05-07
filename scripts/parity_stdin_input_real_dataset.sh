#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT"

DATA="vendor/BBTools-master/resources/sample1.fq.gz"
CP="vendor/BBTools-master/current"
OUT="${1:-tmp/real_stdin_input_parity}"
mkdir -p "$OUT"
rm -f "$OUT"/*

COMMON=(
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

printf 'Running Java BBNorm local baseline on bundled phiX read 1...\n'
java -Xmx1g -cp "$CP" jgi.KmerNormalize \
  "in=$DATA" "${COMMON[@]}" \
  "out=$OUT/java.keep.fq" \
  "outlow=$OUT/java.low.fq" \
  "outmid=$OUT/java.mid.fq" \
  "outhigh=$OUT/java.high.fq" \
  "hist=$OUT/java.hist.tsv" \
  "rhist=$OUT/java.rhist.tsv" \
  >"$OUT/java.stdout.log" 2>"$OUT/java.stderr.log"

printf 'Running Rust bbnorm-rs with plain FASTQ streamed through in=stdin...\n'
gzip -cd "$DATA" | target/debug/bbnorm-rs \
  "in=stdin" "${COMMON[@]}" \
  "out=$OUT/rust.stdin.keep.fq" \
  "outlow=$OUT/rust.stdin.low.fq" \
  "outmid=$OUT/rust.stdin.mid.fq" \
  "outhigh=$OUT/rust.stdin.high.fq" \
  "hist=$OUT/rust.stdin.hist.tsv" \
  "rhist=$OUT/rust.stdin.rhist.tsv" \
  >"$OUT/rust.stdin.stdout.log" 2>"$OUT/rust.stdin.stderr.log"

grep -q 'stdin was materialized from stdin into a temporary file' "$OUT/rust.stdin.stderr.log"

printf 'Running Rust bbnorm-rs with gzip FASTQ streamed through in=stdin.fq.gz...\n'
cat "$DATA" | target/debug/bbnorm-rs \
  "in=stdin.fq.gz" "${COMMON[@]}" \
  "out=$OUT/rust.stdin_gz.keep.fq" \
  "outlow=$OUT/rust.stdin_gz.low.fq" \
  "outmid=$OUT/rust.stdin_gz.mid.fq" \
  "outhigh=$OUT/rust.stdin_gz.high.fq" \
  "hist=$OUT/rust.stdin_gz.hist.tsv" \
  "rhist=$OUT/rust.stdin_gz.rhist.tsv" \
  >"$OUT/rust.stdin_gz.stdout.log" 2>"$OUT/rust.stdin_gz.stderr.log"

grep -q 'stdin.fq.gz was materialized from stdin into a temporary file' "$OUT/rust.stdin_gz.stderr.log"

for rust_prefix in rust.stdin rust.stdin_gz; do
  for suffix in keep.fq low.fq mid.fq high.fq hist.tsv rhist.tsv; do
    cmp "$OUT/java.$suffix" "$OUT/$rust_prefix.$suffix"
  done
done

printf 'Stdin input passed. Rust materializes plain and gzip stdin streams for safe rereads and matches the Java file-input baseline on bundled phiX. Logs: %s\n' "$OUT"

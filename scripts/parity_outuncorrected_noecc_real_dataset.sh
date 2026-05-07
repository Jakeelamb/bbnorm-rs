#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT"

DATA1="vendor/BBTools-master/resources/sample1.fq.gz"
DATA2="vendor/BBTools-master/resources/sample2.fq.gz"
CP="vendor/BBTools-master/current"
OUT="${1:-tmp/real_outuncorrected_noecc_parity}"
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

printf 'Running Java BBNorm outuncorrected without ecc on paired phiX...\n'
java -Xmx1g -cp "$CP" jgi.KmerNormalize \
  "${COMMON[@]}" \
  "out=$OUT/java.keep1.fq" "out2=$OUT/java.keep2.fq" \
  "outuncorrected=$OUT/java.unc1.fq" "outuncorrected2=$OUT/java.unc2.fq" \
  >"$OUT/java.stdout.log" 2>"$OUT/java.stderr.log"

printf 'Running Rust bbnorm-rs outuncorrected without ecc on paired phiX...\n'
target/debug/bbnorm-rs \
  "${COMMON[@]}" \
  "out=$OUT/rust.keep1.fq" "out2=$OUT/rust.keep2.fq" \
  "outuncorrected=$OUT/rust.unc1.fq" "outuncorrected2=$OUT/rust.unc2.fq" \
  >"$OUT/rust.stdout.log" 2>"$OUT/rust.stderr.log"

cmp "$OUT/java.keep1.fq" "$OUT/rust.keep1.fq"
cmp "$OUT/java.keep2.fq" "$OUT/rust.keep2.fq"
cmp "$OUT/java.unc1.fq" "$OUT/rust.unc1.fq"
cmp "$OUT/java.unc2.fq" "$OUT/rust.unc2.fq"
test ! -s "$OUT/rust.unc1.fq"
test ! -s "$OUT/rust.unc2.fq"

printf 'outuncorrected no-ecc parity passed. Logs and outputs: %s\n' "$OUT"

#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT"

DATA1="vendor/BBTools-master/resources/sample1.fq.gz"
DATA2="vendor/BBTools-master/resources/sample2.fq.gz"
CP="vendor/BBTools-master/current"
OUT="${1:-tmp/real_multipass_parity}"
mkdir -p "$OUT"
rm -f "$OUT"/*

COMMON=(
  "in=$DATA1"
  "in2=$DATA2"
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

for passes in 2 3 4; do
  printf 'Running Java BBNorm %s-pass baseline on paired phiX...\n' "$passes"
  java -Xmx1g -cp "$CP" jgi.KmerNormalize \
    "${COMMON[@]}" "passes=$passes" \
    "out=$OUT/java.p${passes}.keep1.fq" "out2=$OUT/java.p${passes}.keep2.fq" \
    "hist=$OUT/java.p${passes}.hist.tsv" \
    >"$OUT/java.p${passes}.stdout.log" 2>"$OUT/java.p${passes}.stderr.log"

  printf 'Running Rust bbnorm-rs %s-pass temp-file orchestration on paired phiX...\n' "$passes"
  target/debug/bbnorm-rs \
    "${COMMON[@]}" "passes=$passes" \
    "out=$OUT/rust.p${passes}.keep1.fq" "out2=$OUT/rust.p${passes}.keep2.fq" \
    "hist=$OUT/rust.p${passes}.hist.tsv" \
    >"$OUT/rust.p${passes}.stdout.log" 2>"$OUT/rust.p${passes}.stderr.log"

  cmp "$OUT/java.p${passes}.keep1.fq" "$OUT/rust.p${passes}.keep1.fq"
  cmp "$OUT/java.p${passes}.keep2.fq" "$OUT/rust.p${passes}.keep2.fq"
  cmp "$OUT/java.p${passes}.hist.tsv" "$OUT/rust.p${passes}.hist.tsv"
done

printf 'Paired phiX passes=2/3/4 parity passed. Rust runs intermediate temp-pass orchestration and matches Java output. Logs: %s\n' "$OUT"

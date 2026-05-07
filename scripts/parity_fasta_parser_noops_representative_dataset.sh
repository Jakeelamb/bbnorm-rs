#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT"

CP="vendor/BBTools-master/current"
OUT="${1:-tmp/representative_fasta_parser_noops_parity}"
mkdir -p "$OUT"
rm -f "$OUT"/*

FASTA="$OUT/representative.fa"
cat > "$FASTA" <<'FASTA'
>seq1
ACGTNN
>seq2
TTTTAC
FASTA

COMMON=(
  "in=$FASTA"
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

for case in \
  "fastaminread=1" \
  "fastaminlen=1" \
  "fastaminlength=1" \
  "forcesectionname=t" \
  "fastadump=f"; do
  stem="${case%%=*}"
  printf 'Running Java/Rust FASTA parser no-op case %s...\n' "$case"
  java -Xmx1g -cp "$CP" jgi.KmerNormalize \
    "${COMMON[@]}" "$case" \
    "out=$OUT/java.$stem.keep.fa" \
    >"$OUT/java.$stem.stdout.log" 2>"$OUT/java.$stem.stderr.log"

  target/debug/bbnorm-rs \
    "${COMMON[@]}" "$case" \
    "out=$OUT/rust.$stem.keep.fa" \
    >"$OUT/rust.$stem.stdout.log" 2>"$OUT/rust.$stem.stderr.log"

  cmp "$OUT/java.$stem.keep.fa" "$OUT/rust.$stem.keep.fa"
done

for option in "fastaminread=abc" "fastaminlen=abc" "fastaminlength=abc"; do
  stem="${option%%=*}"
  printf 'Checking malformed FASTA parser option %s is rejected by Java and Rust...\n' "$option"
  if java -Xmx1g -cp "$CP" jgi.KmerNormalize \
    "${COMMON[@]}" "$option" \
    "out=$OUT/java.$stem.malformed.keep.fa" \
    >"$OUT/java.$stem.malformed.stdout.log" 2>"$OUT/java.$stem.malformed.stderr.log"; then
    echo "Expected Java BBNorm to reject malformed $option" >&2
    exit 1
  fi

  if target/debug/bbnorm-rs \
    "${COMMON[@]}" "$option" \
    "out=$OUT/rust.$stem.malformed.keep.fa" \
    >"$OUT/rust.$stem.malformed.stdout.log" 2>"$OUT/rust.$stem.malformed.stderr.log"; then
    echo "Expected bbnorm-rs to reject malformed $option" >&2
    exit 1
  fi
  test ! -e "$OUT/rust.$stem.malformed.keep.fa"
done

printf 'FASTA parser no-op parity passed. Logs and outputs: %s\n' "$OUT"

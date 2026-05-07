#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT"

DATA1="vendor/BBTools-master/resources/sample1.fq.gz"
DATA2="vendor/BBTools-master/resources/sample2.fq.gz"
CP="vendor/BBTools-master/current"
OUT="${1:-tmp/real_buildstepsize_parity}"
mkdir -p "$OUT"
rm -f "$OUT"/*

COMMON=(
  "in=$DATA1"
  "in2=$DATA2"
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

for case in "buildstepsize=1" "stepsize=2"; do
  label="${case%%=*}"
  value="${case#*=}"
  stem="${label}${value}"

  printf 'Running Java BBNorm %s on paired phiX...\n' "$case"
  java -Xmx1g -cp "$CP" jgi.KmerNormalize \
    "${COMMON[@]}" "$case" \
    "out=$OUT/java.$stem.keep1.fq" "out2=$OUT/java.$stem.keep2.fq" \
    "outlow=$OUT/java.$stem.low1.fq" "outlow2=$OUT/java.$stem.low2.fq" \
    "outmid=$OUT/java.$stem.mid1.fq" "outmid2=$OUT/java.$stem.mid2.fq" \
    "outhigh=$OUT/java.$stem.high1.fq" "outhigh2=$OUT/java.$stem.high2.fq" \
    "hist=$OUT/java.$stem.hist.tsv" "rhist=$OUT/java.$stem.rhist.tsv" \
    >"$OUT/java.$stem.stdout.log" 2>"$OUT/java.$stem.stderr.log"

  printf 'Running Rust bbnorm-rs %s on paired phiX...\n' "$case"
  target/debug/bbnorm-rs \
    "${COMMON[@]}" "$case" \
    "out=$OUT/rust.$stem.keep1.fq" "out2=$OUT/rust.$stem.keep2.fq" \
    "outlow=$OUT/rust.$stem.low1.fq" "outlow2=$OUT/rust.$stem.low2.fq" \
    "outmid=$OUT/rust.$stem.mid1.fq" "outmid2=$OUT/rust.$stem.mid2.fq" \
    "outhigh=$OUT/rust.$stem.high1.fq" "outhigh2=$OUT/rust.$stem.high2.fq" \
    "hist=$OUT/rust.$stem.hist.tsv" "rhist=$OUT/rust.$stem.rhist.tsv" \
    >"$OUT/rust.$stem.stdout.log" 2>"$OUT/rust.$stem.stderr.log"

  for suffix in keep1.fq keep2.fq low1.fq low2.fq mid1.fq mid2.fq high1.fq high2.fq hist.tsv rhist.tsv; do
    cmp "$OUT/java.$stem.$suffix" "$OUT/rust.$stem.$suffix"
  done
done

for suffix in keep1.fq keep2.fq low1.fq low2.fq mid1.fq mid2.fq high1.fq high2.fq hist.tsv rhist.tsv; do
  cmp "$OUT/java.buildstepsize1.$suffix" "$OUT/java.stepsize2.$suffix"
done

printf 'buildstepsize/stepsize paired phiX parity passed. Logs and outputs: %s\n' "$OUT"

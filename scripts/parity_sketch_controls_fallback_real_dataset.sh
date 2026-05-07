#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT"

DATA1="vendor/BBTools-master/resources/sample1.fq.gz"
DATA2="vendor/BBTools-master/resources/sample2.fq.gz"
CP="vendor/BBTools-master/current"
OUT="${1:-tmp/real_sketch_controls_fallback}"
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

run_java() {
  local stem="$1"
  shift
  java -Xmx1g -cp "$CP" jgi.KmerNormalize \
    "${COMMON[@]}" "$@" \
    "out=$OUT/java.$stem.keep1.fq" "out2=$OUT/java.$stem.keep2.fq" \
    "outlow=$OUT/java.$stem.low1.fq" "outlow2=$OUT/java.$stem.low2.fq" \
    "outmid=$OUT/java.$stem.mid1.fq" "outmid2=$OUT/java.$stem.mid2.fq" \
    "outhigh=$OUT/java.$stem.high1.fq" "outhigh2=$OUT/java.$stem.high2.fq" \
    "hist=$OUT/java.$stem.hist.tsv" "rhist=$OUT/java.$stem.rhist.tsv" \
    >"$OUT/java.$stem.stdout.log" 2>"$OUT/java.$stem.stderr.log"
}

run_rust() {
  local stem="$1"
  shift
  target/debug/bbnorm-rs \
    "${COMMON[@]}" "$@" \
    "out=$OUT/rust.$stem.keep1.fq" "out2=$OUT/rust.$stem.keep2.fq" \
    "outlow=$OUT/rust.$stem.low1.fq" "outlow2=$OUT/rust.$stem.low2.fq" \
    "outmid=$OUT/rust.$stem.mid1.fq" "outmid2=$OUT/rust.$stem.mid2.fq" \
    "outhigh=$OUT/rust.$stem.high1.fq" "outhigh2=$OUT/rust.$stem.high2.fq" \
    "hist=$OUT/rust.$stem.hist.tsv" "rhist=$OUT/rust.$stem.rhist.tsv" \
    >"$OUT/rust.$stem.stdout.log" 2>"$OUT/rust.$stem.stderr.log"
}

cargo build --quiet

printf 'Running Java BBNorm baseline and default prefilter cases on paired phiX...\n'
run_java base
run_java prefilter prefilter=t

printf 'Running Rust bbnorm-rs default prefilter case on paired phiX...\n'
run_rust prefilter prefilter=t

for suffix in keep1.fq keep2.fq low1.fq low2.fq mid1.fq mid2.fq high1.fq high2.fq hist.tsv rhist.tsv; do
  cmp "$OUT/java.base.$suffix" "$OUT/java.prefilter.$suffix"
  cmp "$OUT/java.prefilter.$suffix" "$OUT/rust.prefilter.$suffix"
done

for case in "buildpasses=2"; do
  stem="${case%%=*}${case#*=}"
  printf 'Running Java BBNorm behavior-changing %s case on paired phiX...\n' "$case"
  run_java "$stem" "$case"
  if cmp -s "$OUT/java.base.hist.tsv" "$OUT/java.$stem.hist.tsv"; then
    echo "Expected Java histogram to differ for $case" >&2
    exit 1
  fi

  printf 'Running Rust bbnorm-rs behavior-changing %s case...\n' "$case"
  run_rust "$stem" "$case"
  grep -q 'trusted-kmer filtering' "$OUT/rust.$stem.stderr.log"
  if cmp -s "$OUT/rust.prefilter.hist.tsv" "$OUT/rust.$stem.hist.tsv"; then
    echo "Expected Rust build-pass histogram to differ from exact baseline for $case" >&2
    exit 1
  fi
  for suffix in keep1.fq keep2.fq low1.fq low2.fq mid1.fq mid2.fq high1.fq high2.fq rhist.tsv; do
    test -e "$OUT/rust.$stem.$suffix"
  done
done

for case in "prehashes=1" "prefilterhashes=1" "prefiltercells=1k" "precells=1k"; do
  stem="${case%%=*}${case#*=}"
  printf 'Running Java BBNorm accepted prefilter sketch %s case on paired phiX...\n' "$case"
  run_java "$stem" "$case"

  printf 'Running Rust bbnorm-rs prefilter sketch %s case...\n' "$case"
  run_rust "$stem" "$case"
  grep -q 'prefilter collision estimates' "$OUT/rust.$stem.stderr.log"
  if cmp -s "$OUT/rust.prefilter.hist.tsv" "$OUT/rust.$stem.hist.tsv"; then
    echo "Expected Rust prefilter histogram to differ from exact baseline for $case" >&2
    exit 1
  fi
  for suffix in keep1.fq keep2.fq low1.fq low2.fq mid1.fq mid2.fq high1.fq high2.fq rhist.tsv; do
    test -e "$OUT/rust.$stem.$suffix"
  done
done

for case in "cells=1k" "matrixbits=10"; do
  stem="${case%%=*}${case#*=}"
  printf 'Running Java BBNorm constrained sketch %s case on paired phiX...\n' "$case"
  run_java "$stem" "$case"
  if cmp -s "$OUT/java.base.hist.tsv" "$OUT/java.$stem.hist.tsv"; then
    echo "Expected Java histogram to differ for $case" >&2
    exit 1
  fi

  printf 'Running Rust bbnorm-rs constrained count-min %s case...\n' "$case"
  run_rust "$stem" "$case"
  grep -q 'fixed-memory count-min input sketch' "$OUT/rust.$stem.stderr.log"
  if cmp -s "$OUT/rust.prefilter.hist.tsv" "$OUT/rust.$stem.hist.tsv"; then
    echo "Expected Rust count-min histogram to differ from exact baseline for $case" >&2
    exit 1
  fi
  for suffix in keep1.fq keep2.fq low1.fq low2.fq mid1.fq mid2.fq high1.fq high2.fq rhist.tsv; do
    test -e "$OUT/rust.$stem.$suffix"
  done
done

printf 'Running Rust bbnorm-rs explicit sketchmemory budget case...\n'
run_rust sketchmemory sketchmemory=64 hashes=2 bits=8
grep -q 'count-min memory budget' "$OUT/rust.sketchmemory.stderr.log"
if cmp -s "$OUT/rust.prefilter.hist.tsv" "$OUT/rust.sketchmemory.hist.tsv"; then
  echo 'Expected Rust sketchmemory histogram to differ from exact baseline' >&2
  exit 1
fi
for suffix in keep1.fq keep2.fq low1.fq low2.fq mid1.fq mid2.fq high1.fq high2.fq rhist.tsv; do
  test -e "$OUT/rust.sketchmemory.$suffix"
done

for option in "bits=abc" "hashes=abc" "cells=abc" "matrixbits=abc" "sketchmemory=abc" "prefiltersize=abc" "prehashes=abc"; do
  stem="${option%%=*}"
  printf 'Checking malformed sketch/table option %s is rejected by Java and Rust...\n' "$option"
  if java -Xmx1g -cp "$CP" jgi.KmerNormalize \
    "${COMMON[@]}" "$option" \
    "out=null" "out2=null" \
    >"$OUT/java.reject.$stem.stdout.log" 2>"$OUT/java.reject.$stem.stderr.log"; then
    printf 'Expected vendored Java to reject malformed %s, but it succeeded.\n' "$option" >&2
    exit 1
  fi

  if target/debug/bbnorm-rs \
    "${COMMON[@]}" "$option" \
    "out=null" "out2=null" \
    >"$OUT/rust.reject.$stem.stdout.log" 2>"$OUT/rust.reject.$stem.stderr.log"; then
    printf 'Expected Rust to reject malformed %s, but it succeeded.\n' "$option" >&2
    exit 1
  fi
  grep -Eq 'expects|unsupported KMG suffix|too many suffix letters' "$OUT/rust.reject.$stem.stderr.log"
done

printf 'Sketch/table-sizing compatibility passed. Buildpasses, prefilter hash/cell, and cells/matrixbits requests exercise real Rust sketch/table behavior, including direct fixed-memory count-min input sketches; malformed prefilter fraction values are rejected. Logs and outputs: %s\n' "$OUT"

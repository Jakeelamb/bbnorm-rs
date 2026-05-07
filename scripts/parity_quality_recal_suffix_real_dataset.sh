#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT"

DATA1="vendor/BBTools-master/resources/sample1.fq.gz"
DATA2="vendor/BBTools-master/resources/sample2.fq.gz"
CP="vendor/BBTools-master/current"
OUT="${1:-tmp/real_quality_recal_suffix}"
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

RECAL_FLAGS=(
  "loadq102_p1=f"
  "loadq_p2=t"
  "observationcutoff_p1=1k"
  "recalpasses_p2=1"
  "recalqmax_p1=50"
  "recalqmin_p2=2"
  "recalwithposition_p1=t"
  "qmatrixmode_p2=max"
  "recaltile_p1=f"
)

cargo build --quiet

printf 'Running Java BBNorm with pass-suffixed quality recalibration controls on paired phiX...\n'
java -Xmx1g -cp "$CP" jgi.KmerNormalize \
  "${COMMON[@]}" "${RECAL_FLAGS[@]}" \
  "out=$OUT/java.keep1.fq" "out2=$OUT/java.keep2.fq" \
  "outlow=$OUT/java.low1.fq" "outlow2=$OUT/java.low2.fq" \
  "outmid=$OUT/java.mid1.fq" "outmid2=$OUT/java.mid2.fq" \
  "outhigh=$OUT/java.high1.fq" "outhigh2=$OUT/java.high2.fq" \
  "hist=$OUT/java.hist.tsv" "rhist=$OUT/java.rhist.tsv" \
  >"$OUT/java.stdout.log" 2>"$OUT/java.stderr.log"

printf 'Running Rust bbnorm-rs with pass-suffixed quality recalibration controls on paired phiX...\n'
target/debug/bbnorm-rs \
  "${COMMON[@]}" "${RECAL_FLAGS[@]}" \
  "out=$OUT/rust.keep1.fq" "out2=$OUT/rust.keep2.fq" \
  "outlow=$OUT/rust.low1.fq" "outlow2=$OUT/rust.low2.fq" \
  "outmid=$OUT/rust.mid1.fq" "outmid2=$OUT/rust.mid2.fq" \
  "outhigh=$OUT/rust.high1.fq" "outhigh2=$OUT/rust.high2.fq" \
  "hist=$OUT/rust.hist.tsv" "rhist=$OUT/rust.rhist.tsv" \
  >"$OUT/rust.stdout.log" 2>"$OUT/rust.stderr.log"

grep -q 'loadq102_p1=f is a BBTools quality-recalibration control' "$OUT/rust.stderr.log"
grep -q 'observationcutoff_p1=1k is a BBTools quality-recalibration control' "$OUT/rust.stderr.log"
grep -q 'qmatrixmode_p2=max is a BBTools quality-recalibration matrix mode' "$OUT/rust.stderr.log"

for suffix in keep1.fq keep2.fq low1.fq low2.fq mid1.fq mid2.fq high1.fq high2.fq hist.tsv rhist.tsv; do
  cmp "$OUT/java.$suffix" "$OUT/rust.$suffix"
done

printf 'Pass-suffixed quality recalibration compatibility passed. Java accepts _p1/_p2 controls; Rust accepts them as no-ops and matches Java FASTQ/hist output. Logs: %s\n' "$OUT"

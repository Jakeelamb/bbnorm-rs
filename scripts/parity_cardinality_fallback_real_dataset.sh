#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT"

DATA1="vendor/BBTools-master/resources/sample1.fq.gz"
DATA2="vendor/BBTools-master/resources/sample2.fq.gz"
CP="vendor/BBTools-master/current"
OUT="${1:-tmp/real_cardinality_fallback}"
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

CARDINALITY_FLAGS=(
  "cardinality=t"
  "loglog=31"
  "loglogin=t"
  "cardinalityout=t"
  "loglogout=t"
  "buckets=1k"
  "loglogcorrection=t"
  "loglogbits=16"
  "loglogk=31"
  "loglogklist=21,31"
  "loglogseed=42"
  "loglogminprob=0.5"
  "loglogtype=loglog2"
  "loglogmean=t"
  "loglogcounts=t"
)

cargo build --quiet

printf 'Running Java BBNorm local baseline on paired phiX...\n'
java -Xmx1g -cp "$CP" jgi.KmerNormalize \
  "${COMMON[@]}" \
  "out=$OUT/java.keep1.fq" "out2=$OUT/java.keep2.fq" \
  "outlow=$OUT/java.low1.fq" "outlow2=$OUT/java.low2.fq" \
  "outmid=$OUT/java.mid1.fq" "outmid2=$OUT/java.mid2.fq" \
  "outhigh=$OUT/java.high1.fq" "outhigh2=$OUT/java.high2.fq" \
  "hist=$OUT/java.hist.tsv" "rhist=$OUT/java.rhist.tsv" \
  >"$OUT/java.stdout.log" 2>"$OUT/java.stderr.log"

printf 'Confirming vendored KmerNormalize rejects cardinality/loglog controls...\n'
if java -Xmx1g -cp "$CP" jgi.KmerNormalize \
  "${COMMON[@]}" "${CARDINALITY_FLAGS[@]}" \
  "out=$OUT/java.reject.keep1.fq" "out2=$OUT/java.reject.keep2.fq" \
  >"$OUT/java.reject.stdout.log" 2>"$OUT/java.reject.stderr.log"; then
  printf 'Expected vendored Java to reject cardinality/loglog controls, but it succeeded.\n' >&2
  exit 1
fi
grep -q 'Unknown parameter cardinality=' "$OUT/java.reject.stderr.log"

printf 'Running Rust bbnorm-rs with cardinality/loglog controls on paired phiX...\n'
target/debug/bbnorm-rs \
  "${COMMON[@]}" "${CARDINALITY_FLAGS[@]}" \
  "out=$OUT/rust.keep1.fq" "out2=$OUT/rust.keep2.fq" \
  "outlow=$OUT/rust.low1.fq" "outlow2=$OUT/rust.low2.fq" \
  "outmid=$OUT/rust.mid1.fq" "outmid2=$OUT/rust.mid2.fq" \
  "outhigh=$OUT/rust.high1.fq" "outhigh2=$OUT/rust.high2.fq" \
  "hist=$OUT/rust.hist.tsv" "rhist=$OUT/rust.rhist.tsv" \
  >"$OUT/rust.stdout.log" 2>"$OUT/rust.stderr.log"

grep -q 'cardinality=t is a BBTools cardinality/loglog control' "$OUT/rust.stderr.log"
grep -q 'loglog=31 is a BBTools cardinality/loglog control' "$OUT/rust.stderr.log"
grep -q 'buckets=1k is a BBTools cardinality/loglog bucket control' "$OUT/rust.stderr.log"
grep -q 'loglogklist=21,31 is a BBTools cardinality/loglog k-list' "$OUT/rust.stderr.log"
grep -q 'Cardinality estimate: scope=input; k=21; buckets=1000;' "$OUT/rust.stderr.log"
grep -q 'Cardinality estimate: scope=output; k=21; buckets=1000;' "$OUT/rust.stderr.log"
grep -q 'Stage timing: name=input_cardinality;' "$OUT/rust.stderr.log"

for suffix in keep1.fq keep2.fq low1.fq low2.fq mid1.fq mid2.fq high1.fq high2.fq hist.tsv rhist.tsv; do
  cmp "$OUT/java.$suffix" "$OUT/rust.$suffix"
done

printf 'Cardinality/loglog bounded-estimate parity passed. Vendored KmerNormalize rejects these controls; Rust accepts them, emits input/output estimates, and matches local Java FASTQ/hist baseline. Logs: %s\n' "$OUT"

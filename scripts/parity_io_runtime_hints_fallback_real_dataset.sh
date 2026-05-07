#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT"

DATA1="vendor/BBTools-master/resources/sample1.fq.gz"
DATA2="vendor/BBTools-master/resources/sample2.fq.gz"
CP="vendor/BBTools-master/current"
OUT="${1:-tmp/real_io_runtime_hints_fallback}"
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

HINT_FLAGS=(
  "extin=.fq.gz"
  "extout=.fq"
  "workers=auto"
  "workerthreads=1"
  "wt=auto"
  "threadsin=1"
  "tin=auto"
  "threadsout=1"
  "tout=auto"
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

printf 'Confirming vendored KmerNormalize rejects shared I/O hint controls...\n'
if java -Xmx1g -cp "$CP" jgi.KmerNormalize \
  "${COMMON[@]}" "${HINT_FLAGS[@]}" \
  "out=$OUT/java.reject.keep1.fq" "out2=$OUT/java.reject.keep2.fq" \
  >"$OUT/java.reject.stdout.log" 2>"$OUT/java.reject.stderr.log"; then
  printf 'Expected vendored Java to reject shared I/O hint controls, but it succeeded.\n' >&2
  exit 1
fi
grep -Eq 'Unknown parameter (extin|workers|workerthreads|wt|threadsin|tin|threadsout|tout)=' "$OUT/java.reject.stderr.log"

printf 'Running Rust bbnorm-rs with shared I/O hint fallbacks on paired phiX...\n'
target/debug/bbnorm-rs \
  "${COMMON[@]}" "${HINT_FLAGS[@]}" \
  "out=$OUT/rust.keep1.fq" "out2=$OUT/rust.keep2.fq" \
  "outlow=$OUT/rust.low1.fq" "outlow2=$OUT/rust.low2.fq" \
  "outmid=$OUT/rust.mid1.fq" "outmid2=$OUT/rust.mid2.fq" \
  "outhigh=$OUT/rust.high1.fq" "outhigh2=$OUT/rust.high2.fq" \
  "hist=$OUT/rust.hist.tsv" "rhist=$OUT/rust.rhist.tsv" \
  >"$OUT/rust.stdout.log" 2>"$OUT/rust.stderr.log"

grep -q 'extin=.fq.gz is a BBTools file-extension hint' "$OUT/rust.stderr.log"
grep -q 'extout=.fq is a BBTools file-extension hint' "$OUT/rust.stderr.log"
grep -q 'workers=auto is a BBTools I/O worker control' "$OUT/rust.stderr.log"
grep -q 'threadsin=1 is a BBTools I/O worker control' "$OUT/rust.stderr.log"
grep -q 'threadsout=1 is a BBTools I/O worker control' "$OUT/rust.stderr.log"

for suffix in keep1.fq keep2.fq low1.fq low2.fq mid1.fq mid2.fq high1.fq high2.fq hist.tsv rhist.tsv; do
  cmp "$OUT/java.$suffix" "$OUT/rust.$suffix"
done

printf 'I/O hint fallback passed. Vendored KmerNormalize rejects these shared controls; Rust accepts them and matches local Java FASTQ/hist baseline. Logs: %s\n' "$OUT"

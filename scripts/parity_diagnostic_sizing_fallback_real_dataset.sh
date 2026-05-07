#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT"

DATA1="vendor/BBTools-master/resources/sample1.fq.gz"
DATA2="vendor/BBTools-master/resources/sample2.fq.gz"
CP="vendor/BBTools-master/current"
OUT="${1:-tmp/real_diagnostic_sizing_fallback}"
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

SIZING_FLAGS=(
  "testsize=t"
  "breaklen=0"
  "breaklength=-1"
  "recalibrate=f"
  "recalibratequality=f"
  "recal=f"
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

printf 'Confirming vendored KmerNormalize rejects shared diagnostic sizing controls...\n'
if java -Xmx1g -cp "$CP" jgi.KmerNormalize \
  "${COMMON[@]}" "${SIZING_FLAGS[@]}" \
  "out=$OUT/java.reject.keep1.fq" "out2=$OUT/java.reject.keep2.fq" \
  >"$OUT/java.reject.stdout.log" 2>"$OUT/java.reject.stderr.log"; then
  printf 'Expected vendored Java to reject diagnostic sizing controls, but it succeeded.\n' >&2
  exit 1
fi
grep -q 'Unknown parameter testsize=' "$OUT/java.reject.stderr.log"

printf 'Confirming Rust still rejects output-affecting break/recalibrate controls instead of silently ignoring them...\n'
if target/debug/bbnorm-rs "${COMMON[@]}" "breaklen=50" "out=$OUT/rust.break.reject.fq" \
  >"$OUT/rust.break.reject.stdout.log" 2>"$OUT/rust.break.reject.stderr.log"; then
  printf 'Expected Rust to reject breaklen, but it succeeded.\n' >&2
  exit 1
fi
grep -q 'breaklen=50 enables BBTools read breaking' "$OUT/rust.break.reject.stderr.log"
if target/debug/bbnorm-rs "${COMMON[@]}" "breaklen=abc" "out=$OUT/rust.break.malformed.reject.fq" \
  >"$OUT/rust.break.malformed.reject.stdout.log" 2>"$OUT/rust.break.malformed.reject.stderr.log"; then
  printf 'Expected Rust to reject malformed breaklen=abc, but it succeeded.\n' >&2
  exit 1
fi
grep -q 'breaklen expects' "$OUT/rust.break.malformed.reject.stderr.log"
if target/debug/bbnorm-rs "${COMMON[@]}" "recalibrate=t" "out=$OUT/rust.recal.reject.fq" \
  >"$OUT/rust.recal.reject.stdout.log" 2>"$OUT/rust.recal.reject.stderr.log"; then
  printf 'Expected Rust to reject recalibrate=t, but it succeeded.\n' >&2
  exit 1
fi
grep -q 'recalibrate=t enables BBTools quality recalibration' "$OUT/rust.recal.reject.stderr.log"
if target/debug/bbnorm-rs "${COMMON[@]}" "recalibrate=maybe" "out=$OUT/rust.recal.malformed.reject.fq" \
  >"$OUT/rust.recal.malformed.reject.stdout.log" 2>"$OUT/rust.recal.malformed.reject.stderr.log"; then
  printf 'Expected Rust to reject malformed recalibrate=maybe, but it succeeded.\n' >&2
  exit 1
fi
grep -q 'recalibrate expects a boolean value' "$OUT/rust.recal.malformed.reject.stderr.log"

printf 'Running Rust bbnorm-rs with diagnostic sizing fallback controls on paired phiX...\n'
target/debug/bbnorm-rs \
  "${COMMON[@]}" "${SIZING_FLAGS[@]}" \
  "out=$OUT/rust.keep1.fq" "out2=$OUT/rust.keep2.fq" \
  "outlow=$OUT/rust.low1.fq" "outlow2=$OUT/rust.low2.fq" \
  "outmid=$OUT/rust.mid1.fq" "outmid2=$OUT/rust.mid2.fq" \
  "outhigh=$OUT/rust.high1.fq" "outhigh2=$OUT/rust.high2.fq" \
  "hist=$OUT/rust.hist.tsv" "rhist=$OUT/rust.rhist.tsv" \
  >"$OUT/rust.stdout.log" 2>"$OUT/rust.stderr.log"

grep -q 'testsize=t is a BBTools diagnostic sizing control' "$OUT/rust.stderr.log"
grep -q 'breaklen=0 keeps BBTools read breaking disabled' "$OUT/rust.stderr.log"
grep -q 'breaklength=-1 keeps BBTools read breaking disabled' "$OUT/rust.stderr.log"
grep -q 'recalibrate=f keeps BBTools quality recalibration disabled' "$OUT/rust.stderr.log"
grep -q 'recalibratequality=f keeps BBTools quality recalibration disabled' "$OUT/rust.stderr.log"
grep -q 'recal=f keeps BBTools quality recalibration disabled' "$OUT/rust.stderr.log"

for suffix in keep1.fq keep2.fq low1.fq low2.fq mid1.fq mid2.fq high1.fq high2.fq hist.tsv rhist.tsv; do
  cmp "$OUT/java.$suffix" "$OUT/rust.$suffix"
done

printf 'Diagnostic sizing fallback passed. Vendored KmerNormalize rejects testsize; Rust accepts diagnostic sizing plus disabled break/recalibration controls as no-ops, rejects output-affecting controls, and matches local Java FASTQ/hist baseline. Logs: %s\n' "$OUT"

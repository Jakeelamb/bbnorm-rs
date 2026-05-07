#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT"

DATA1="vendor/BBTools-master/resources/sample1.fq.gz"
DATA2="vendor/BBTools-master/resources/sample2.fq.gz"
CP="vendor/BBTools-master/current"
OUT="${1:-tmp/real_tmpdir_controls_parity}"
mkdir -p "$OUT"
rm -rf "$OUT"/*
mkdir -p "$OUT/java_tmp" "$OUT/rust_tmp"

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

printf 'Running Java BBNorm with temporary-directory controls on paired phiX...\n'
java -Xmx1g -cp "$CP" jgi.KmerNormalize \
  "${COMMON[@]}" "tmpdir=$OUT/java_tmp" "usetmpdir=t" \
  "out=$OUT/java.keep1.fq" "out2=$OUT/java.keep2.fq" \
  "outlow=$OUT/java.low1.fq" "outlow2=$OUT/java.low2.fq" \
  "outmid=$OUT/java.mid1.fq" "outmid2=$OUT/java.mid2.fq" \
  "outhigh=$OUT/java.high1.fq" "outhigh2=$OUT/java.high2.fq" \
  "hist=$OUT/java.hist.tsv" "rhist=$OUT/java.rhist.tsv" \
  >"$OUT/java.stdout.log" 2>"$OUT/java.stderr.log"

if find "$OUT/java_tmp" -mindepth 1 -print -quit | grep -q .; then
  printf 'Vendored Java unexpectedly created single-pass temp files under %s.\n' "$OUT/java_tmp" >&2
  exit 1
fi

printf 'Running Rust bbnorm-rs with temporary-directory controls on paired phiX...\n'
target/debug/bbnorm-rs \
  "${COMMON[@]}" "tmpdir=$OUT/rust_tmp" "usetmpdir=t" "usetempdir=f" \
  "out=$OUT/rust.keep1.fq" "out2=$OUT/rust.keep2.fq" \
  "outlow=$OUT/rust.low1.fq" "outlow2=$OUT/rust.low2.fq" \
  "outmid=$OUT/rust.mid1.fq" "outmid2=$OUT/rust.mid2.fq" \
  "outhigh=$OUT/rust.high1.fq" "outhigh2=$OUT/rust.high2.fq" \
  "hist=$OUT/rust.hist.tsv" "rhist=$OUT/rust.rhist.tsv" \
  >"$OUT/rust.stdout.log" 2>"$OUT/rust.stderr.log"

grep -q 'tmpdir=.* is a BBTools temporary-directory control' "$OUT/rust.stderr.log"
grep -q 'usetmpdir=t is a BBTools temporary-directory control' "$OUT/rust.stderr.log"
grep -q 'usetempdir=f is a BBTools temporary-directory control' "$OUT/rust.stderr.log"

if find "$OUT/rust_tmp" -mindepth 1 -print -quit | grep -q .; then
  printf 'Rust unexpectedly created single-pass temp files under %s.\n' "$OUT/rust_tmp" >&2
  exit 1
fi

for suffix in keep1.fq keep2.fq low1.fq low2.fq mid1.fq mid2.fq high1.fq high2.fq hist.tsv rhist.tsv; do
  cmp "$OUT/java.$suffix" "$OUT/rust.$suffix"
done

printf 'Confirming Rust multipass uses the requested managed temp parent and cleans up after itself...\n'
rm -rf "$OUT/rust_multipass_tmp"
target/debug/bbnorm-rs \
  "in=$DATA1" \
  "in2=$DATA2" \
  "passes=2" \
  "keepall=t" \
  "lowbindepth=1" \
  "highbindepth=1" \
  "k=31" \
  "minq=0" \
  "minprob=0" \
  "min=0" \
  "minkmers=1" \
  "target=999999999" \
  "max=999999999" \
  "threads=1" \
  "overwrite=t" \
  "bits=32" \
  "tmpdir=$OUT/rust_multipass_tmp" \
  "usetmpdir=t" \
  "out=$OUT/rust.multipass.keep1.fq" \
  "out2=$OUT/rust.multipass.keep2.fq" \
  "hist=$OUT/rust.multipass.hist.tsv" \
  >"$OUT/rust.multipass.stdout.log" 2>"$OUT/rust.multipass.stderr.log"

grep -q 'tmpdir=.* is a BBTools temporary-directory control' "$OUT/rust.multipass.stderr.log"
[[ -d "$OUT/rust_multipass_tmp" ]]
if find "$OUT/rust_multipass_tmp" -mindepth 1 -print -quit | grep -q .; then
  printf 'Rust left managed multipass temp files under %s.\n' "$OUT/rust_multipass_tmp" >&2
  exit 1
fi

printf 'Temporary-directory control parity passed. Java and Rust both leave tmpdir empty in the covered single-pass paired phiX path, Rust uses and cleans the requested managed temp parent for multipass, and outputs match. Logs: %s\n' "$OUT"

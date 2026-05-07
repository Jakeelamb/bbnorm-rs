#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT"

DATA1="vendor/BBTools-master/resources/sample1.fq.gz"
DATA2="vendor/BBTools-master/resources/sample2.fq.gz"
CP="vendor/BBTools-master/current"
OUT="${1:-tmp/real_sam_runtime_noops}"
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

SAM_FLAGS=(
  "sam=1.4"
  "samtools=f"
  "nativebam=f"
  "attachedsamline=f"
  "streamerthreads=1"
  "fastqstreamerthreads=1"
  "fastqwriterthreads=1"
  "prefermd=f"
  "mdtag=f"
  "notags=f"
  "xstag=us"
  "readgroupid=rg1"
  "rgsm=sample1"
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

printf 'Confirming vendored KmerNormalize rejects shared SAM/readgroup controls...\n'
if java -Xmx1g -cp "$CP" jgi.KmerNormalize \
  "${COMMON[@]}" "${SAM_FLAGS[@]}" \
  "out=$OUT/java.reject.keep1.fq" "out2=$OUT/java.reject.keep2.fq" \
  >"$OUT/java.reject.stdout.log" 2>"$OUT/java.reject.stderr.log"; then
  printf 'Expected vendored Java to reject SAM/readgroup controls, but it succeeded.\n' >&2
  exit 1
fi
grep -q 'Unknown parameter sam=1.4' "$OUT/java.reject.stderr.log"

printf 'Running Rust bbnorm-rs with shared SAM/readgroup runtime controls on paired phiX...\n'
target/debug/bbnorm-rs \
  "${COMMON[@]}" "${SAM_FLAGS[@]}" \
  "out=$OUT/rust.keep1.fq" "out2=$OUT/rust.keep2.fq" \
  "outlow=$OUT/rust.low1.fq" "outlow2=$OUT/rust.low2.fq" \
  "outmid=$OUT/rust.mid1.fq" "outmid2=$OUT/rust.mid2.fq" \
  "outhigh=$OUT/rust.high1.fq" "outhigh2=$OUT/rust.high2.fq" \
  "hist=$OUT/rust.hist.tsv" "rhist=$OUT/rust.rhist.tsv" \
  >"$OUT/rust.stdout.log" 2>"$OUT/rust.stderr.log"

grep -q 'sam=1.4 is a BBTools SAM-version control' "$OUT/rust.stderr.log"
grep -q 'fastqwriterthreads=1 is a BBTools SAM/streamer threading control' "$OUT/rust.stderr.log"
grep -q 'readgroupid=rg1 is a BBTools read-group metadata control' "$OUT/rust.stderr.log"

for suffix in keep1.fq keep2.fq low1.fq low2.fq mid1.fq mid2.fq high1.fq high2.fq hist.tsv rhist.tsv; do
  cmp "$OUT/java.$suffix" "$OUT/rust.$suffix"
done

printf 'SAM/readgroup runtime no-op fallback passed. Vendored KmerNormalize rejects these controls; Rust accepts them and matches local Java FASTQ baseline. Logs: %s\n' "$OUT"

#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT"

DATA1="${DATA1:-vendor/BBTools-master/resources/sample1.fq.gz}"
DATA2="${DATA2:-vendor/BBTools-master/resources/sample2.fq.gz}"
CP="vendor/BBTools-master/current"
OUT="${1:-tmp/real_kmer_table_runtime_fallback}"
READS="${READS:-}"
TABLE_READS="${TABLE_READS:-$READS}"
THREADS="${THREADS:-1}"
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
  "threads=$THREADS"
  "overwrite=t"
  "bits=32"
)

if [[ -n "$READS" ]]; then
  COMMON+=("reads=$READS")
fi

if [[ -n "$TABLE_READS" ]]; then
  COMMON+=("tablereads=$TABLE_READS")
fi

TABLE_FLAGS=(
  "initialsize=1k"
  "ways=31"
  "buflen=64k"
  "bufflen=64k"
  "bufferlength=64k"
  "tabletype=2"
  "rcomp=t"
  "maskmiddle=f"
  "showstats=t"
  "stats=f"
  "showspeed=f"
  "ss=t"
  "verbose2=t"
  "prealloc=0.25"
  "preallocate=f"
  "minprobprefilter=f"
  "mpp=t"
  "minprobmain=t"
  "mpm=f"
  "prefilterpasses=auto"
  "prepasses=1"
  "onepass=t"
)

PREFILTER_MEMORY_FLAGS=(
  "filtermemory=1k"
  "prefiltermemory=1k"
  "filtermem=1k"
  "filtermemoryoverride=1k"
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

printf 'Confirming vendored KmerNormalize rejects shared kmer-table runtime controls...\n'
if java -Xmx1g -cp "$CP" jgi.KmerNormalize \
  "${COMMON[@]}" "${TABLE_FLAGS[@]}" "${PREFILTER_MEMORY_FLAGS[@]}" \
  "out=$OUT/java.reject.keep1.fq" "out2=$OUT/java.reject.keep2.fq" \
  >"$OUT/java.reject.stdout.log" 2>"$OUT/java.reject.stderr.log"; then
  printf 'Expected vendored Java to reject shared kmer-table runtime controls, but it succeeded.\n' >&2
  exit 1
fi
grep -Eq 'Unknown parameter (initialsize|ways|buflen|bufflen|bufferlength|tabletype|rcomp|maskmiddle|showstats|stats|showspeed|ss|verbose2|prealloc|preallocate|filtermemory|prefiltermemory|filtermem|filtermemoryoverride|minprobprefilter|mpp|minprobmain|mpm|prefilterpasses|prepasses|onepass)=' "$OUT/java.reject.stderr.log"

printf 'Running Rust bbnorm-rs with kmer-table runtime fallbacks on paired phiX...\n'
target/debug/bbnorm-rs \
  "${COMMON[@]}" "${TABLE_FLAGS[@]}" \
  "out=$OUT/rust.keep1.fq" "out2=$OUT/rust.keep2.fq" \
  "outlow=$OUT/rust.low1.fq" "outlow2=$OUT/rust.low2.fq" \
  "outmid=$OUT/rust.mid1.fq" "outmid2=$OUT/rust.mid2.fq" \
  "outhigh=$OUT/rust.high1.fq" "outhigh2=$OUT/rust.high2.fq" \
  "hist=$OUT/rust.hist.tsv" "rhist=$OUT/rust.rhist.tsv" \
  >"$OUT/rust.stdout.log" 2>"$OUT/rust.stderr.log"

grep -q 'initialsize=1k is a BBTools kmer-table runtime sizing control' "$OUT/rust.stderr.log"
grep -q 'buflen=64k is a BBTools kmer-table buffer-length control' "$OUT/rust.stderr.log"
grep -q 'tabletype=2 is a BBTools kmer-table implementation control' "$OUT/rust.stderr.log"
grep -q 'rcomp=t is a BBTools kmer-table matching control' "$OUT/rust.stderr.log"
grep -q 'showstats=t is a BBTools kmer-table reporting control' "$OUT/rust.stderr.log"
grep -q 'prealloc=0.25 is a BBTools kmer-table preallocation control' "$OUT/rust.stderr.log"
grep -q 'prefilterpasses=auto is a BBTools prefilter pass-count control' "$OUT/rust.stderr.log"
grep -q 'onepass=t is a BBTools kmer-table construction-mode control' "$OUT/rust.stderr.log"

for suffix in keep1.fq keep2.fq low1.fq low2.fq mid1.fq mid2.fq high1.fq high2.fq hist.tsv rhist.tsv; do
  cmp "$OUT/java.$suffix" "$OUT/rust.$suffix"
done

for option in "${PREFILTER_MEMORY_FLAGS[@]}"; do
  stem="${option//[^A-Za-z0-9]/_}"
  printf 'Running Rust bbnorm-rs real prefilter memory-sizing option %s on paired phiX...\n' "$option"
  target/debug/bbnorm-rs \
    "${COMMON[@]}" "$option" \
    "out=$OUT/rust.$stem.keep1.fq" "out2=$OUT/rust.$stem.keep2.fq" \
    "hist=$OUT/rust.$stem.hist.tsv" "rhist=$OUT/rust.$stem.rhist.tsv" \
    >"$OUT/rust.$stem.stdout.log" 2>"$OUT/rust.$stem.stderr.log"
  grep -q "${option} is a BBTools prefilter memory-sizing control" "$OUT/rust.$stem.stderr.log"
  if cmp -s "$OUT/java.hist.tsv" "$OUT/rust.$stem.hist.tsv"; then
    printf 'Expected prefilter memory-sizing option %s to exercise real Rust prefilter collision behavior.\n' "$option" >&2
    exit 1
  fi
done

for option in "initialsize=abc" "ways=abc" "buflen=abc" "tabletype=abc" "prealloc=0.abc" "prealloc=1.5" "filtermemory=abc" "prepasses=abc"; do
  stem="${option//[^A-Za-z0-9]/_}"
  printf 'Checking malformed kmer-table runtime option %s is rejected by Rust...\n' "$option"
  if target/debug/bbnorm-rs \
    "${COMMON[@]}" "$option" \
    "out=null" "out2=null" \
    >"$OUT/rust.reject.$stem.stdout.log" 2>"$OUT/rust.reject.$stem.stderr.log"; then
    printf 'Expected Rust to reject malformed %s, but it succeeded.\n' "$option" >&2
    exit 1
  fi
  grep -Eq 'expects|unsupported KMG suffix|too many suffix letters' "$OUT/rust.reject.$stem.stderr.log"
done

printf 'Kmer-table runtime fallback passed. Vendored KmerNormalize rejects these shared controls; Rust accepts pure runtime controls with Java-matching output, and prefilter memory controls exercise real Rust collision behavior. Logs: %s\n' "$OUT"

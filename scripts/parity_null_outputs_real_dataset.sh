#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
DATA1="$ROOT/vendor/BBTools-master/resources/sample1.fq.gz"
DATA2="$ROOT/vendor/BBTools-master/resources/sample2.fq.gz"
CP="$ROOT/vendor/BBTools-master/current"
OUT="${1:-$ROOT/tmp/real_null_outputs_parity}"
if [[ "$OUT" != /* ]]; then
  OUT="$ROOT/$OUT"
fi
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

cargo build --quiet --manifest-path "$ROOT/Cargo.toml"

pushd "$OUT" >/dev/null

printf 'Running Java BBNorm with null sequence outputs on paired phiX...\n'
java -Xmx1g -cp "$CP" jgi.KmerNormalize \
  "${COMMON[@]}" \
  "out=null" "out2=null" \
  "outlow=null" "outlow2=null" \
  "outmid=null" "outmid2=null" \
  "outhigh=null" "outhigh2=null" \
  "hist=$OUT/java.hist.tsv" "rhist=$OUT/java.rhist.tsv" \
  >"$OUT/java.stdout.log" 2>"$OUT/java.stderr.log"

if [[ -e null ]]; then
  printf 'Vendored Java unexpectedly created a literal null output in %s.\n' "$OUT" >&2
  exit 1
fi

printf 'Running Rust bbnorm-rs with null sequence outputs on paired phiX...\n'
"$ROOT/target/debug/bbnorm-rs" \
  "${COMMON[@]}" \
  "out=null" "out2=null" \
  "outlow=null" "outlow2=null" \
  "outmid=null" "outmid2=null" \
  "outhigh=null" "outhigh2=null" \
  "hist=$OUT/rust.hist.tsv" "rhist=$OUT/rust.rhist.tsv" \
  >"$OUT/rust.stdout.log" 2>"$OUT/rust.stderr.log"

if [[ -e null ]]; then
  printf 'Rust created a literal null output in %s instead of treating out=null as a sink.\n' "$OUT" >&2
  exit 1
fi

cmp "$OUT/java.hist.tsv" "$OUT/rust.hist.tsv"
cmp "$OUT/java.rhist.tsv" "$OUT/rust.rhist.tsv"

printf 'Running Java BBNorm with uppercase NULL sequence and hist outputs on paired phiX...\n'
java -Xmx1g -cp "$CP" jgi.KmerNormalize \
  "${COMMON[@]}" \
  "out=NULL" "out2=NULL" \
  "outlow=NULL" "outlow2=NULL" \
  "outmid=NULL" "outmid2=NULL" \
  "outhigh=NULL" "outhigh2=NULL" \
  "hist=NULL" "rhist=$OUT/java.upper.rhist.tsv" \
  >"$OUT/java.upper.stdout.log" 2>"$OUT/java.upper.stderr.log"

if [[ -e NULL ]]; then
  printf 'Vendored Java unexpectedly created a literal NULL output in %s.\n' "$OUT" >&2
  exit 1
fi

printf 'Running Rust bbnorm-rs with uppercase NULL sequence and hist outputs on paired phiX...\n'
"$ROOT/target/debug/bbnorm-rs" \
  "${COMMON[@]}" \
  "out=NULL" "out2=NULL" \
  "outlow=NULL" "outlow2=NULL" \
  "outmid=NULL" "outmid2=NULL" \
  "outhigh=NULL" "outhigh2=NULL" \
  "hist=NULL" "rhist=$OUT/rust.upper.rhist.tsv" \
  >"$OUT/rust.upper.stdout.log" 2>"$OUT/rust.upper.stderr.log"

if [[ -e NULL ]]; then
  printf 'Rust created a literal NULL output in %s instead of treating NULL as a case-insensitive sink.\n' "$OUT" >&2
  exit 1
fi

cmp "$OUT/java.upper.rhist.tsv" "$OUT/rust.upper.rhist.tsv"

popd >/dev/null

printf 'Null output parity passed. Java and Rust both suppress literal null/NULL sequence and hist outputs and retain matching paired phiX hist/rhist outputs. Logs: %s\n' "$OUT"

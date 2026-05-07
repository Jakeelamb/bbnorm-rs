#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT"

DATA1="vendor/BBTools-master/resources/sample1.fq.gz"
DATA2="vendor/BBTools-master/resources/sample2.fq.gz"
CP="vendor/BBTools-master/current"
OUT="${1:-tmp/real_runtime_noops_parity}"
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

run_pair_case() {
  local stem="$1"
  shift

  printf 'Running Java BBNorm runtime no-op case %s on paired phiX...\n' "$stem"
  java -Xmx1g -cp "$CP" jgi.KmerNormalize \
    "${COMMON[@]}" "$@" \
    "out=$OUT/java.$stem.keep1.fq" "out2=$OUT/java.$stem.keep2.fq" \
    "outlow=$OUT/java.$stem.low1.fq" "outlow2=$OUT/java.$stem.low2.fq" \
    "outmid=$OUT/java.$stem.mid1.fq" "outmid2=$OUT/java.$stem.mid2.fq" \
    "outhigh=$OUT/java.$stem.high1.fq" "outhigh2=$OUT/java.$stem.high2.fq" \
    "hist=$OUT/java.$stem.hist.tsv" "rhist=$OUT/java.$stem.rhist.tsv" \
    >"$OUT/java.$stem.stdout.log" 2>"$OUT/java.$stem.stderr.log"

  printf 'Running Rust bbnorm-rs runtime no-op case %s on paired phiX...\n' "$stem"
  target/debug/bbnorm-rs \
    "${COMMON[@]}" "$@" \
    "out=$OUT/rust.$stem.keep1.fq" "out2=$OUT/rust.$stem.keep2.fq" \
    "outlow=$OUT/rust.$stem.low1.fq" "outlow2=$OUT/rust.$stem.low2.fq" \
    "outmid=$OUT/rust.$stem.mid1.fq" "outmid2=$OUT/rust.$stem.mid2.fq" \
    "outhigh=$OUT/rust.$stem.high1.fq" "outhigh2=$OUT/rust.$stem.high2.fq" \
    "hist=$OUT/rust.$stem.hist.tsv" "rhist=$OUT/rust.$stem.rhist.tsv" \
    >"$OUT/rust.$stem.stdout.log" 2>"$OUT/rust.$stem.stderr.log"

  for suffix in keep1.fq keep2.fq low1.fq low2.fq mid1.fq mid2.fq high1.fq high2.fq hist.tsv rhist.tsv; do
    cmp "$OUT/java.$stem.$suffix" "$OUT/rust.$stem.$suffix"
  done
}

run_pair_case ordered_false ordered=f
run_pair_case verbose_true verbose=t
run_pair_case printcoverage_true printcoverage=t

FASTA="$OUT/representative.fa"
cat > "$FASTA" <<'FASTA'
>seq1
ACGTNN
>seq2
TTTTAC
FASTA

for case in "fastareadlen=4" "fastareadlength=4"; do
  stem="${case%%=*}"
  printf 'Running Java/Rust FASTA runtime no-op case %s...\n' "$case"
  java -Xmx1g -cp "$CP" jgi.KmerNormalize \
    "in=$FASTA" "passes=1" "keepall=t" "k=3" "minq=0" "minprob=0" \
    "min=0" "minkmers=1" "target=999999999" "max=999999999" \
    "threads=1" "overwrite=t" "bits=32" "$case" \
    "out=$OUT/java.$stem.keep.fa" \
    >"$OUT/java.$stem.stdout.log" 2>"$OUT/java.$stem.stderr.log"

  target/debug/bbnorm-rs \
    "in=$FASTA" "passes=1" "keepall=t" "k=3" "minq=0" "minprob=0" \
    "min=0" "minkmers=1" "target=999999999" "max=999999999" \
    "threads=1" "overwrite=t" "bits=32" "$case" \
    "out=$OUT/rust.$stem.keep.fa" \
    >"$OUT/rust.$stem.stdout.log" 2>"$OUT/rust.$stem.stderr.log"

  cmp "$OUT/java.$stem.keep.fa" "$OUT/rust.$stem.keep.fa"
done

printf 'Checking Rust auto=f exact-count fallback...\n'
run_pair_case auto_false auto=f
grep -q 'automatic count-table sizing control' "$OUT/rust.auto_false.stderr.log"

printf 'Runtime no-op parity passed. Logs and outputs: %s\n' "$OUT"

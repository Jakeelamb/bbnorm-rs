#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT"

OUT="${1:-tmp/component_smoke}"
THREADS="${RAYON_NUM_THREADS:-$(getconf _NPROCESSORS_ONLN 2>/dev/null || printf '1')}"
export RAYON_NUM_THREADS="$THREADS"

mkdir -p "$OUT"
rm -f "$OUT"/*

printf 'Building bbnorm-rs with RAYON_NUM_THREADS=%s...\n' "$RAYON_NUM_THREADS"
cargo build --quiet

printf 'Running fast Rust integration tests...\n'
cargo test --quiet --test basic

REP_FASTQ="$OUT/representative.fq"
cat > "$REP_FASTQ" <<'FASTQ'
@r1
ACGTACGTACGT
+
IIIIIIIIIIII
@r2
ACGTACGTACGT
+
IIIIIIIIIIII
@r3
TTTTACGTAAAA
+
IIIIIIIIIIII
@r4
NNNNACGTNNNN
+
############
FASTQ

printf 'Running generated single-end keepall smoke...\n'
target/debug/bbnorm-rs \
  "in=$REP_FASTQ" \
  "out=$OUT/generated.keep.fq" \
  "hist=$OUT/generated.hist.tsv" \
  "rhist=$OUT/generated.rhist.tsv" \
  "passes=1" \
  "keepall=t" \
  "k=4" \
  "minq=0" \
  "minprob=0" \
  "min=0" \
  "minkmers=1" \
  "target=999999999" \
  "max=999999999" \
  "threads=$THREADS" \
  "overwrite=t" \
  "bits=32" \
  >"$OUT/generated.stdout.log" 2>"$OUT/generated.stderr.log"

grep -q '^@r1$' "$OUT/generated.keep.fq"
grep -q '^#Depth[[:space:]]Raw_Count[[:space:]]Unique_Kmers$' "$OUT/generated.hist.tsv"
grep -q '^#Depth[[:space:]]Reads[[:space:]]Bases$' "$OUT/generated.rhist.tsv"

printf 'Running generated single-end keep/toss smoke...\n'
target/debug/bbnorm-rs \
  "in=$REP_FASTQ" \
  "out=$OUT/generated.norm.keep.fq" \
  "outt=$OUT/generated.norm.toss.fq" \
  "passes=1" \
  "k=4" \
  "minq=0" \
  "minprob=0" \
  "min=0" \
  "minkmers=1" \
  "target=1" \
  "max=1" \
  "threads=$THREADS" \
  "overwrite=t" \
  "bits=32" \
  >"$OUT/generated.norm.stdout.log" 2>"$OUT/generated.norm.stderr.log"

test -s "$OUT/generated.norm.keep.fq"
test -e "$OUT/generated.norm.toss.fq"

DATA1="vendor/BBTools-master/resources/sample1.fq.gz"
DATA2="vendor/BBTools-master/resources/sample2.fq.gz"
if [[ -f "$DATA1" && -f "$DATA2" ]]; then
  printf 'Running bundled paired phiX keep/bin/hist smoke...\n'
  target/debug/bbnorm-rs \
    "in=$DATA1" \
    "in2=$DATA2" \
    "out=$OUT/phix.keep1.fq" \
    "out2=$OUT/phix.keep2.fq" \
    "outlow=$OUT/phix.low1.fq" \
    "outlow2=$OUT/phix.low2.fq" \
    "outmid=$OUT/phix.mid1.fq" \
    "outmid2=$OUT/phix.mid2.fq" \
    "outhigh=$OUT/phix.high1.fq" \
    "outhigh2=$OUT/phix.high2.fq" \
    "hist=$OUT/phix.hist.tsv" \
    "rhist=$OUT/phix.rhist.tsv" \
    "passes=1" \
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
    "threads=$THREADS" \
    "overwrite=t" \
    "bits=32" \
    >"$OUT/phix.stdout.log" 2>"$OUT/phix.stderr.log"

  test -s "$OUT/phix.keep1.fq"
  test -s "$OUT/phix.keep2.fq"
  grep -q '^#Depth[[:space:]]Raw_Count[[:space:]]Unique_Kmers$' "$OUT/phix.hist.tsv"
  grep -q '^#Depth[[:space:]]Reads[[:space:]]Bases$' "$OUT/phix.rhist.tsv"
else
  printf 'Skipping bundled paired phiX smoke; fixture files are not present.\n'
fi

printf 'Component smoke passed. Outputs and logs: %s\n' "$OUT"

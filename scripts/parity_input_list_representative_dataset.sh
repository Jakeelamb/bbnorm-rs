#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT"

CP="vendor/BBTools-master/current"
OUT="${1:-tmp/representative_input_list_parity}"
mkdir -p "$OUT"
rm -f "$OUT"/*

INPUT1="$OUT/representative_list_a.fq"
INPUT2="$OUT/representative_list_b.fq"
INPUT3="$OUT/representative_list_c.fq"
R1A="$OUT/representative_pair_list_a_1.fq"
R2A="$OUT/representative_pair_list_a_2.fq"
R1B="$OUT/representative_pair_list_b_1.fq"
R2B="$OUT/representative_pair_list_b_2.fq"
R1C="$OUT/representative_pair_list_c_1.fq"
R2C="$OUT/representative_pair_list_c_2.fq"
cat > "$INPUT1" <<'FASTQ'
@a1
ACGTACGT
+
IIIIIIII
@a2
CCCCCCCC
+
IIIIIIII
FASTQ
cat > "$INPUT2" <<'FASTQ'
@b1
TTTTTTTT
+
IIIIIIII
@b2
GGGGGGGG
+
IIIIIIII
FASTQ
cat > "$INPUT3" <<'FASTQ'
@c1
ACGT
+
IIII
FASTQ
cat > "$R1A" <<'FASTQ'
@a1/1
AAA
+
III
@a2/1
CCC
+
III
FASTQ
cat > "$R2A" <<'FASTQ'
@a1/2
GGG
+
III
@a2/2
TTT
+
III
FASTQ
cat > "$R1B" <<'FASTQ'
@b1/1
ACG
+
III
@b2/1
CGT
+
III
FASTQ
cat > "$R2B" <<'FASTQ'
@b1/2
GTA
+
III
@b2/2
TAC
+
III
FASTQ
cat > "$R1C" <<'FASTQ'
@c1/1
AAA
+
III
FASTQ
cat > "$R2C" <<'FASTQ'
@c1/2
GGG
+
III
FASTQ

COMMON=(
  "in=$INPUT1,$INPUT2"
  "passes=1"
  "keepall=t"
  "k=3"
  "minq=0"
  "minprob=0"
  "min=0"
  "minkmers=1"
  "target=999999999"
  "max=999999999"
  "threads=1"
  "overwrite=t"
  "bits=32"
  "reads=1"
)

cargo build --quiet

printf 'Running Java BBNorm single-end input-list representative case...\n'
java -Xmx1g -cp "$CP" jgi.KmerNormalize \
  "${COMMON[@]}" \
  "out=$OUT/java.keep.fq" \
  "hist=$OUT/java.hist.tsv" \
  >"$OUT/java.stdout.log" 2>"$OUT/java.stderr.log"

printf 'Running Rust bbnorm-rs single-end input-list representative case...\n'
target/debug/bbnorm-rs \
  "${COMMON[@]}" \
  "out=$OUT/rust.keep.fq" \
  "hist=$OUT/rust.hist.tsv" \
  >"$OUT/rust.stdout.log" 2>"$OUT/rust.stderr.log"

cmp "$OUT/java.keep.fq" "$OUT/rust.keep.fq"
cmp "$OUT/java.hist.tsv" "$OUT/rust.hist.tsv"

printf 'Running Java BBNorm single-end input-list output-fanout representative case...\n'
java -Xmx1g -cp "$CP" jgi.KmerNormalize \
  "${COMMON[@]}" \
  "out=$OUT/java.list.keep.a.fq,$OUT/java.list.keep.b.fq" \
  "hist=$OUT/java.list.hist.tsv" \
  >"$OUT/java.list.stdout.log" 2>"$OUT/java.list.stderr.log"

printf 'Running Rust bbnorm-rs single-end input-list output-fanout representative case...\n'
target/debug/bbnorm-rs \
  "${COMMON[@]}" \
  "out=$OUT/rust.list.keep.a.fq,$OUT/rust.list.keep.b.fq" \
  "hist=$OUT/rust.list.hist.tsv" \
  >"$OUT/rust.list.stdout.log" 2>"$OUT/rust.list.stderr.log"

cmp "$OUT/java.list.keep.a.fq" "$OUT/rust.list.keep.a.fq"
cmp "$OUT/java.list.keep.b.fq" "$OUT/rust.list.keep.b.fq"
cmp "$OUT/java.list.hist.tsv" "$OUT/rust.list.hist.tsv"

printf 'Running Java BBNorm single-end input-list low-bin output-fanout representative case...\n'
java -Xmx1g -cp "$CP" jgi.KmerNormalize \
  "${COMMON[@]}" \
  "out=$OUT/java.low-list.keep.fq" \
  "outlow=$OUT/java.low-list.low.a.fq,$OUT/java.low-list.low.b.fq" \
  "hist=$OUT/java.low-list.hist.tsv" \
  >"$OUT/java.low-list.stdout.log" 2>"$OUT/java.low-list.stderr.log"

printf 'Running Rust bbnorm-rs single-end input-list low-bin output-fanout representative case...\n'
target/debug/bbnorm-rs \
  "${COMMON[@]}" \
  "out=$OUT/rust.low-list.keep.fq" \
  "outlow=$OUT/rust.low-list.low.a.fq,$OUT/rust.low-list.low.b.fq" \
  "hist=$OUT/rust.low-list.hist.tsv" \
  >"$OUT/rust.low-list.stdout.log" 2>"$OUT/rust.low-list.stderr.log"

cmp "$OUT/java.low-list.keep.fq" "$OUT/rust.low-list.keep.fq"
cmp "$OUT/java.low-list.low.a.fq" "$OUT/rust.low-list.low.a.fq"
cmp "$OUT/java.low-list.low.b.fq" "$OUT/rust.low-list.low.b.fq"
cmp "$OUT/java.low-list.hist.tsv" "$OUT/rust.low-list.hist.tsv"

for BIN_SPEC in "mid outmid 80" "high outhigh 0"; do
  read -r BIN_NAME BIN_OPT HIGH_BIN_DEPTH <<<"$BIN_SPEC"
  printf 'Running Java BBNorm single-end input-list %s-bin output-fanout representative case...\n' "$BIN_NAME"
  java -Xmx1g -cp "$CP" jgi.KmerNormalize \
    "${COMMON[@]}" \
    "lowbindepth=0" \
    "highbindepth=$HIGH_BIN_DEPTH" \
    "out=$OUT/java.$BIN_NAME-list.keep.fq" \
    "$BIN_OPT=$OUT/java.$BIN_NAME-list.$BIN_NAME.a.fq,$OUT/java.$BIN_NAME-list.$BIN_NAME.b.fq" \
    "hist=$OUT/java.$BIN_NAME-list.hist.tsv" \
    >"$OUT/java.$BIN_NAME-list.stdout.log" 2>"$OUT/java.$BIN_NAME-list.stderr.log"

  printf 'Running Rust bbnorm-rs single-end input-list %s-bin output-fanout representative case...\n' "$BIN_NAME"
  target/debug/bbnorm-rs \
    "${COMMON[@]}" \
    "lowbindepth=0" \
    "highbindepth=$HIGH_BIN_DEPTH" \
    "out=$OUT/rust.$BIN_NAME-list.keep.fq" \
    "$BIN_OPT=$OUT/rust.$BIN_NAME-list.$BIN_NAME.a.fq,$OUT/rust.$BIN_NAME-list.$BIN_NAME.b.fq" \
    "hist=$OUT/rust.$BIN_NAME-list.hist.tsv" \
    >"$OUT/rust.$BIN_NAME-list.stdout.log" 2>"$OUT/rust.$BIN_NAME-list.stderr.log"

  cmp "$OUT/java.$BIN_NAME-list.keep.fq" "$OUT/rust.$BIN_NAME-list.keep.fq"
  cmp "$OUT/java.$BIN_NAME-list.$BIN_NAME.a.fq" "$OUT/rust.$BIN_NAME-list.$BIN_NAME.a.fq"
  cmp "$OUT/java.$BIN_NAME-list.$BIN_NAME.b.fq" "$OUT/rust.$BIN_NAME-list.$BIN_NAME.b.fq"
  cmp "$OUT/java.$BIN_NAME-list.hist.tsv" "$OUT/rust.$BIN_NAME-list.hist.tsv"
done

TOSS_COMMON=(
  "in=$INPUT1,$INPUT2"
  "passes=1"
  "k=3"
  "minq=0"
  "minprob=0"
  "min=0"
  "minkmers=1"
  "target=999999999"
  "max=999999999"
  "threads=1"
  "overwrite=t"
  "bits=32"
  "reads=1"
  "minlen=9"
)

printf 'Running Java BBNorm single-end input-list toss output-fanout representative case...\n'
java -Xmx1g -cp "$CP" jgi.KmerNormalize \
  "${TOSS_COMMON[@]}" \
  "out=$OUT/java.toss-list.keep.fq" \
  "outt=$OUT/java.toss-list.toss.a.fq,$OUT/java.toss-list.toss.b.fq" \
  "hist=$OUT/java.toss-list.hist.tsv" \
  >"$OUT/java.toss-list.stdout.log" 2>"$OUT/java.toss-list.stderr.log"

printf 'Running Rust bbnorm-rs single-end input-list toss output-fanout representative case...\n'
target/debug/bbnorm-rs \
  "${TOSS_COMMON[@]}" \
  "out=$OUT/rust.toss-list.keep.fq" \
  "outt=$OUT/rust.toss-list.toss.a.fq,$OUT/rust.toss-list.toss.b.fq" \
  "hist=$OUT/rust.toss-list.hist.tsv" \
  >"$OUT/rust.toss-list.stdout.log" 2>"$OUT/rust.toss-list.stderr.log"

cmp "$OUT/java.toss-list.keep.fq" "$OUT/rust.toss-list.keep.fq"
cmp "$OUT/java.toss-list.toss.a.fq" "$OUT/rust.toss-list.toss.a.fq"
cmp "$OUT/java.toss-list.toss.b.fq" "$OUT/rust.toss-list.toss.b.fq"
cmp "$OUT/java.toss-list.hist.tsv" "$OUT/rust.toss-list.hist.tsv"

MISMATCH_COMMON=(
  "in=$INPUT1,$INPUT2,$INPUT3"
  "passes=1"
  "keepall=t"
  "k=3"
  "minq=0"
  "minprob=0"
  "min=0"
  "minkmers=1"
  "target=999999999"
  "max=999999999"
  "threads=1"
  "overwrite=t"
  "bits=32"
  "reads=1"
)

printf 'Running Java BBNorm short output-list negative representative case...\n'
set +e
java -Xmx1g -cp "$CP" jgi.KmerNormalize \
  "${MISMATCH_COMMON[@]}" \
  "out=$OUT/java.short-list.keep.a.fq,$OUT/java.short-list.keep.b.fq" \
  "hist=$OUT/java.short-list.hist.tsv" \
  >"$OUT/java.short-list.stdout.log" 2>"$OUT/java.short-list.stderr.log"
JAVA_SHORT_STATUS=$?
target/debug/bbnorm-rs \
  "${MISMATCH_COMMON[@]}" \
  "out=$OUT/rust.short-list.keep.a.fq,$OUT/rust.short-list.keep.b.fq" \
  "hist=$OUT/rust.short-list.hist.tsv" \
  >"$OUT/rust.short-list.stdout.log" 2>"$OUT/rust.short-list.stderr.log"
RUST_SHORT_STATUS=$?
set -e

test "$JAVA_SHORT_STATUS" -ne 0
test "$RUST_SHORT_STATUS" -ne 0
cmp "$OUT/java.short-list.keep.a.fq" "$OUT/rust.short-list.keep.a.fq"
cmp "$OUT/java.short-list.keep.b.fq" "$OUT/rust.short-list.keep.b.fq"

PAIRED_COMMON=(
  "in=$R1A,$R1B"
  "in2=$R2A,$R2B"
  "passes=1"
  "keepall=t"
  "k=3"
  "minq=0"
  "minprob=0"
  "min=0"
  "minkmers=1"
  "target=999999999"
  "max=999999999"
  "threads=1"
  "overwrite=t"
  "bits=32"
  "reads=1"
)

printf 'Running Java BBNorm paired input-list representative case...\n'
java -Xmx1g -cp "$CP" jgi.KmerNormalize \
  "${PAIRED_COMMON[@]}" \
  "out=$OUT/java.paired.keep1.fq" \
  "out2=$OUT/java.paired.keep2.fq" \
  "hist=$OUT/java.paired.hist.tsv" \
  >"$OUT/java.paired.stdout.log" 2>"$OUT/java.paired.stderr.log"

printf 'Running Rust bbnorm-rs paired input-list representative case...\n'
target/debug/bbnorm-rs \
  "${PAIRED_COMMON[@]}" \
  "out=$OUT/rust.paired.keep1.fq" \
  "out2=$OUT/rust.paired.keep2.fq" \
  "hist=$OUT/rust.paired.hist.tsv" \
  >"$OUT/rust.paired.stdout.log" 2>"$OUT/rust.paired.stderr.log"

cmp "$OUT/java.paired.keep1.fq" "$OUT/rust.paired.keep1.fq"
cmp "$OUT/java.paired.keep2.fq" "$OUT/rust.paired.keep2.fq"
cmp "$OUT/java.paired.hist.tsv" "$OUT/rust.paired.hist.tsv"

printf 'Running Java BBNorm paired input-list output-fanout representative case...\n'
java -Xmx1g -cp "$CP" jgi.KmerNormalize \
  "${PAIRED_COMMON[@]}" \
  "out=$OUT/java.paired.list.keep1.a.fq,$OUT/java.paired.list.keep1.b.fq" \
  "out2=$OUT/java.paired.list.keep2.a.fq,$OUT/java.paired.list.keep2.b.fq" \
  "hist=$OUT/java.paired.list.hist.tsv" \
  >"$OUT/java.paired.list.stdout.log" 2>"$OUT/java.paired.list.stderr.log"

printf 'Running Rust bbnorm-rs paired input-list output-fanout representative case...\n'
target/debug/bbnorm-rs \
  "${PAIRED_COMMON[@]}" \
  "out=$OUT/rust.paired.list.keep1.a.fq,$OUT/rust.paired.list.keep1.b.fq" \
  "out2=$OUT/rust.paired.list.keep2.a.fq,$OUT/rust.paired.list.keep2.b.fq" \
  "hist=$OUT/rust.paired.list.hist.tsv" \
  >"$OUT/rust.paired.list.stdout.log" 2>"$OUT/rust.paired.list.stderr.log"

cmp "$OUT/java.paired.list.keep1.a.fq" "$OUT/rust.paired.list.keep1.a.fq"
cmp "$OUT/java.paired.list.keep1.b.fq" "$OUT/rust.paired.list.keep1.b.fq"
cmp "$OUT/java.paired.list.keep2.a.fq" "$OUT/rust.paired.list.keep2.a.fq"
cmp "$OUT/java.paired.list.keep2.b.fq" "$OUT/rust.paired.list.keep2.b.fq"
cmp "$OUT/java.paired.list.hist.tsv" "$OUT/rust.paired.list.hist.tsv"

printf 'Running Java BBNorm paired input-list low-bin output-fanout representative case...\n'
java -Xmx1g -cp "$CP" jgi.KmerNormalize \
  "${PAIRED_COMMON[@]}" \
  "out=$OUT/java.paired.low-list.keep1.fq" \
  "out2=$OUT/java.paired.low-list.keep2.fq" \
  "outlow=$OUT/java.paired.low-list.low1.a.fq,$OUT/java.paired.low-list.low1.b.fq" \
  "outlow2=$OUT/java.paired.low-list.low2.a.fq,$OUT/java.paired.low-list.low2.b.fq" \
  "hist=$OUT/java.paired.low-list.hist.tsv" \
  >"$OUT/java.paired.low-list.stdout.log" 2>"$OUT/java.paired.low-list.stderr.log"

printf 'Running Rust bbnorm-rs paired input-list low-bin output-fanout representative case...\n'
target/debug/bbnorm-rs \
  "${PAIRED_COMMON[@]}" \
  "out=$OUT/rust.paired.low-list.keep1.fq" \
  "out2=$OUT/rust.paired.low-list.keep2.fq" \
  "outlow=$OUT/rust.paired.low-list.low1.a.fq,$OUT/rust.paired.low-list.low1.b.fq" \
  "outlow2=$OUT/rust.paired.low-list.low2.a.fq,$OUT/rust.paired.low-list.low2.b.fq" \
  "hist=$OUT/rust.paired.low-list.hist.tsv" \
  >"$OUT/rust.paired.low-list.stdout.log" 2>"$OUT/rust.paired.low-list.stderr.log"

cmp "$OUT/java.paired.low-list.keep1.fq" "$OUT/rust.paired.low-list.keep1.fq"
cmp "$OUT/java.paired.low-list.keep2.fq" "$OUT/rust.paired.low-list.keep2.fq"
cmp "$OUT/java.paired.low-list.low1.a.fq" "$OUT/rust.paired.low-list.low1.a.fq"
cmp "$OUT/java.paired.low-list.low1.b.fq" "$OUT/rust.paired.low-list.low1.b.fq"
cmp "$OUT/java.paired.low-list.low2.a.fq" "$OUT/rust.paired.low-list.low2.a.fq"
cmp "$OUT/java.paired.low-list.low2.b.fq" "$OUT/rust.paired.low-list.low2.b.fq"
cmp "$OUT/java.paired.low-list.hist.tsv" "$OUT/rust.paired.low-list.hist.tsv"

for BIN_SPEC in "mid outmid outmid2 80" "high outhigh outhigh2 0"; do
  read -r BIN_NAME BIN_OPT1 BIN_OPT2 HIGH_BIN_DEPTH <<<"$BIN_SPEC"
  printf 'Running Java BBNorm paired input-list %s-bin output-fanout representative case...\n' "$BIN_NAME"
  java -Xmx1g -cp "$CP" jgi.KmerNormalize \
    "${PAIRED_COMMON[@]}" \
    "lowbindepth=0" \
    "highbindepth=$HIGH_BIN_DEPTH" \
    "out=$OUT/java.paired.$BIN_NAME-list.keep1.fq" \
    "out2=$OUT/java.paired.$BIN_NAME-list.keep2.fq" \
    "$BIN_OPT1=$OUT/java.paired.$BIN_NAME-list.${BIN_NAME}1.a.fq,$OUT/java.paired.$BIN_NAME-list.${BIN_NAME}1.b.fq" \
    "$BIN_OPT2=$OUT/java.paired.$BIN_NAME-list.${BIN_NAME}2.a.fq,$OUT/java.paired.$BIN_NAME-list.${BIN_NAME}2.b.fq" \
    "hist=$OUT/java.paired.$BIN_NAME-list.hist.tsv" \
    >"$OUT/java.paired.$BIN_NAME-list.stdout.log" 2>"$OUT/java.paired.$BIN_NAME-list.stderr.log"

  printf 'Running Rust bbnorm-rs paired input-list %s-bin output-fanout representative case...\n' "$BIN_NAME"
  target/debug/bbnorm-rs \
    "${PAIRED_COMMON[@]}" \
    "lowbindepth=0" \
    "highbindepth=$HIGH_BIN_DEPTH" \
    "out=$OUT/rust.paired.$BIN_NAME-list.keep1.fq" \
    "out2=$OUT/rust.paired.$BIN_NAME-list.keep2.fq" \
    "$BIN_OPT1=$OUT/rust.paired.$BIN_NAME-list.${BIN_NAME}1.a.fq,$OUT/rust.paired.$BIN_NAME-list.${BIN_NAME}1.b.fq" \
    "$BIN_OPT2=$OUT/rust.paired.$BIN_NAME-list.${BIN_NAME}2.a.fq,$OUT/rust.paired.$BIN_NAME-list.${BIN_NAME}2.b.fq" \
    "hist=$OUT/rust.paired.$BIN_NAME-list.hist.tsv" \
    >"$OUT/rust.paired.$BIN_NAME-list.stdout.log" 2>"$OUT/rust.paired.$BIN_NAME-list.stderr.log"

  cmp "$OUT/java.paired.$BIN_NAME-list.keep1.fq" "$OUT/rust.paired.$BIN_NAME-list.keep1.fq"
  cmp "$OUT/java.paired.$BIN_NAME-list.keep2.fq" "$OUT/rust.paired.$BIN_NAME-list.keep2.fq"
  cmp "$OUT/java.paired.$BIN_NAME-list.${BIN_NAME}1.a.fq" "$OUT/rust.paired.$BIN_NAME-list.${BIN_NAME}1.a.fq"
  cmp "$OUT/java.paired.$BIN_NAME-list.${BIN_NAME}1.b.fq" "$OUT/rust.paired.$BIN_NAME-list.${BIN_NAME}1.b.fq"
  cmp "$OUT/java.paired.$BIN_NAME-list.${BIN_NAME}2.a.fq" "$OUT/rust.paired.$BIN_NAME-list.${BIN_NAME}2.a.fq"
  cmp "$OUT/java.paired.$BIN_NAME-list.${BIN_NAME}2.b.fq" "$OUT/rust.paired.$BIN_NAME-list.${BIN_NAME}2.b.fq"
  cmp "$OUT/java.paired.$BIN_NAME-list.hist.tsv" "$OUT/rust.paired.$BIN_NAME-list.hist.tsv"
done

PAIRED_TOSS_COMMON=(
  "in=$R1A,$R1B"
  "in2=$R2A,$R2B"
  "passes=1"
  "k=3"
  "minq=0"
  "minprob=0"
  "min=0"
  "minkmers=1"
  "target=999999999"
  "max=999999999"
  "threads=1"
  "overwrite=t"
  "bits=32"
  "reads=1"
  "minlen=4"
)

printf 'Running Java BBNorm paired input-list toss output-fanout representative case...\n'
java -Xmx1g -cp "$CP" jgi.KmerNormalize \
  "${PAIRED_TOSS_COMMON[@]}" \
  "out=$OUT/java.paired.toss-list.keep1.fq" \
  "out2=$OUT/java.paired.toss-list.keep2.fq" \
  "outt=$OUT/java.paired.toss-list.toss1.a.fq,$OUT/java.paired.toss-list.toss1.b.fq" \
  "outt2=$OUT/java.paired.toss-list.toss2.a.fq,$OUT/java.paired.toss-list.toss2.b.fq" \
  "hist=$OUT/java.paired.toss-list.hist.tsv" \
  >"$OUT/java.paired.toss-list.stdout.log" 2>"$OUT/java.paired.toss-list.stderr.log"

printf 'Running Rust bbnorm-rs paired input-list toss output-fanout representative case...\n'
target/debug/bbnorm-rs \
  "${PAIRED_TOSS_COMMON[@]}" \
  "out=$OUT/rust.paired.toss-list.keep1.fq" \
  "out2=$OUT/rust.paired.toss-list.keep2.fq" \
  "outt=$OUT/rust.paired.toss-list.toss1.a.fq,$OUT/rust.paired.toss-list.toss1.b.fq" \
  "outt2=$OUT/rust.paired.toss-list.toss2.a.fq,$OUT/rust.paired.toss-list.toss2.b.fq" \
  "hist=$OUT/rust.paired.toss-list.hist.tsv" \
  >"$OUT/rust.paired.toss-list.stdout.log" 2>"$OUT/rust.paired.toss-list.stderr.log"

cmp "$OUT/java.paired.toss-list.keep1.fq" "$OUT/rust.paired.toss-list.keep1.fq"
cmp "$OUT/java.paired.toss-list.keep2.fq" "$OUT/rust.paired.toss-list.keep2.fq"
cmp "$OUT/java.paired.toss-list.toss1.a.fq" "$OUT/rust.paired.toss-list.toss1.a.fq"
cmp "$OUT/java.paired.toss-list.toss1.b.fq" "$OUT/rust.paired.toss-list.toss1.b.fq"
cmp "$OUT/java.paired.toss-list.toss2.a.fq" "$OUT/rust.paired.toss-list.toss2.a.fq"
cmp "$OUT/java.paired.toss-list.toss2.b.fq" "$OUT/rust.paired.toss-list.toss2.b.fq"
cmp "$OUT/java.paired.toss-list.hist.tsv" "$OUT/rust.paired.toss-list.hist.tsv"

PAIRED_MISMATCH_COMMON=(
  "in=$R1A,$R1B,$R1C"
  "in2=$R2A,$R2B,$R2C"
  "passes=1"
  "keepall=t"
  "k=3"
  "minq=0"
  "minprob=0"
  "min=0"
  "minkmers=1"
  "target=999999999"
  "max=999999999"
  "threads=1"
  "overwrite=t"
  "bits=32"
  "reads=1"
)

printf 'Running Java BBNorm paired short out2-list negative representative case...\n'
set +e
java -Xmx1g -cp "$CP" jgi.KmerNormalize \
  "${PAIRED_MISMATCH_COMMON[@]}" \
  "out=$OUT/java.paired.short-second.keep1.a.fq,$OUT/java.paired.short-second.keep1.b.fq,$OUT/java.paired.short-second.keep1.c.fq" \
  "out2=$OUT/java.paired.short-second.keep2.a.fq,$OUT/java.paired.short-second.keep2.b.fq" \
  "hist=$OUT/java.paired.short-second.hist.tsv" \
  >"$OUT/java.paired.short-second.stdout.log" 2>"$OUT/java.paired.short-second.stderr.log"
JAVA_PAIRED_SHORT_STATUS=$?
target/debug/bbnorm-rs \
  "${PAIRED_MISMATCH_COMMON[@]}" \
  "out=$OUT/rust.paired.short-second.keep1.a.fq,$OUT/rust.paired.short-second.keep1.b.fq,$OUT/rust.paired.short-second.keep1.c.fq" \
  "out2=$OUT/rust.paired.short-second.keep2.a.fq,$OUT/rust.paired.short-second.keep2.b.fq" \
  "hist=$OUT/rust.paired.short-second.hist.tsv" \
  >"$OUT/rust.paired.short-second.stdout.log" 2>"$OUT/rust.paired.short-second.stderr.log"
RUST_PAIRED_SHORT_STATUS=$?
set -e

test "$JAVA_PAIRED_SHORT_STATUS" -ne 0
test "$RUST_PAIRED_SHORT_STATUS" -ne 0
cmp "$OUT/java.paired.short-second.keep1.a.fq" "$OUT/rust.paired.short-second.keep1.a.fq"
cmp "$OUT/java.paired.short-second.keep1.b.fq" "$OUT/rust.paired.short-second.keep1.b.fq"
test ! -e "$OUT/java.paired.short-second.keep1.c.fq"
test ! -e "$OUT/rust.paired.short-second.keep1.c.fq"
cmp "$OUT/java.paired.short-second.keep2.a.fq" "$OUT/rust.paired.short-second.keep2.a.fq"
cmp "$OUT/java.paired.short-second.keep2.b.fq" "$OUT/rust.paired.short-second.keep2.b.fq"

printf 'Running Java BBNorm paired input-list second-output-list-only representative case...\n'
java -Xmx1g -cp "$CP" jgi.KmerNormalize \
  "${PAIRED_COMMON[@]}" \
  "out=$OUT/java.paired.second-list.keep1.fq" \
  "out2=$OUT/java.paired.second-list.keep2.a.fq,$OUT/java.paired.second-list.keep2.b.fq" \
  "hist=$OUT/java.paired.second-list.hist.tsv" \
  >"$OUT/java.paired.second-list.stdout.log" 2>"$OUT/java.paired.second-list.stderr.log"

printf 'Running Rust bbnorm-rs paired input-list second-output-list-only representative case...\n'
target/debug/bbnorm-rs \
  "${PAIRED_COMMON[@]}" \
  "out=$OUT/rust.paired.second-list.keep1.fq" \
  "out2=$OUT/rust.paired.second-list.keep2.a.fq,$OUT/rust.paired.second-list.keep2.b.fq" \
  "hist=$OUT/rust.paired.second-list.hist.tsv" \
  >"$OUT/rust.paired.second-list.stdout.log" 2>"$OUT/rust.paired.second-list.stderr.log"

cmp "$OUT/java.paired.second-list.keep1.fq" "$OUT/rust.paired.second-list.keep1.fq"
cmp "$OUT/java.paired.second-list.keep2.a.fq" "$OUT/rust.paired.second-list.keep2.a.fq"
test ! -e "$OUT/java.paired.second-list.keep2.b.fq"
test ! -e "$OUT/rust.paired.second-list.keep2.b.fq"
cmp "$OUT/java.paired.second-list.hist.tsv" "$OUT/rust.paired.second-list.hist.tsv"

printf 'Input-list representative parity passed for single-end and paired lists, including output fanout edge cases. Logs and outputs: %s\n' "$OUT"

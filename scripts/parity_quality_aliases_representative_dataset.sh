#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT"

CP="vendor/BBTools-master/current"
OUT="${1:-tmp/representative_quality_aliases_parity}"
mkdir -p "$OUT"
rm -f "$OUT"/*

Q33="$OUT/representative_q33.fq"
Q64="$OUT/representative_q64.fq"
R1_Q33="$OUT/representative_qauto_r1.fq"
R2_Q64="$OUT/representative_qauto_r2.fq"
INTERLEAVED="$OUT/representative_qauto_interleaved.fq"
cat > "$Q33" <<'FASTQ'
@q33
ACGTNNACGT
+
!#I~IIIIII
FASTQ
cat > "$Q64" <<'FASTQ'
@q64
ACGTNNACGT
+
@Bh|hhhhhh
FASTQ
cat > "$R1_Q33" <<'FASTQ'
@pair0/1
ACGTNNACGT
+
!#I~IIIIII
@pair1/1
TGCATGCATG
+
IIIIIIIIII
FASTQ
cat > "$R2_Q64" <<'FASTQ'
@pair0/2
TGCANNACGT
+
@Bh|hhhhhh
@pair1/2
ACGTACGTAC
+
hhhhhhhhhh
FASTQ
cat > "$INTERLEAVED" <<'FASTQ'
@pair0/1
ACGTNNACGT
+
!#I~IIIIII
@pair0/2
TGCANNACGT
+
@Bh|hhhhhh
@pair1/1
TGCATGCATG
+
IIIIIIIIII
@pair1/2
ACGTACGTAC
+
hhhhhhhhhh
FASTQ

COMMON=(
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
)

cargo build --quiet

run_case() {
  local label="$1"
  local input="$2"
  shift 2
  local args=("$@")

  printf 'Running Java BBNorm quality alias case %s...\n' "$label"
  java -Xmx1g -cp "$CP" jgi.KmerNormalize \
    "in=$input" \
    "${COMMON[@]}" \
    "${args[@]}" \
    "out=$OUT/java.$label.keep.fq" \
    >"$OUT/java.$label.stdout.log" 2>"$OUT/java.$label.stderr.log"

  printf 'Running Rust bbnorm-rs quality alias case %s...\n' "$label"
  target/debug/bbnorm-rs \
    "in=$input" \
    "${COMMON[@]}" \
    "${args[@]}" \
    "out=$OUT/rust.$label.keep.fq" \
    >"$OUT/rust.$label.stdout.log" 2>"$OUT/rust.$label.stderr.log"

  cmp "$OUT/java.$label.keep.fq" "$OUT/rust.$label.keep.fq"
}

run_list_case() {
  local label="$1"
  local input_arg="$2"
  shift 2
  local args=("$@")

  printf 'Running Java BBNorm quality alias input-list case %s...\n' "$label"
  java -Xmx1g -cp "$CP" jgi.KmerNormalize \
    "in=$input_arg" \
    "${COMMON[@]}" \
    "${args[@]}" \
    "out=$OUT/java.$label.a.fq,$OUT/java.$label.b.fq" \
    "hist=$OUT/java.$label.hist.tsv" \
    >"$OUT/java.$label.stdout.log" 2>"$OUT/java.$label.stderr.log"

  printf 'Running Rust bbnorm-rs quality alias input-list case %s...\n' "$label"
  target/debug/bbnorm-rs \
    "in=$input_arg" \
    "${COMMON[@]}" \
    "${args[@]}" \
    "out=$OUT/rust.$label.a.fq,$OUT/rust.$label.b.fq" \
    "hist=$OUT/rust.$label.hist.tsv" \
    >"$OUT/rust.$label.stdout.log" 2>"$OUT/rust.$label.stderr.log"

  cmp "$OUT/java.$label.a.fq" "$OUT/rust.$label.a.fq"
  cmp "$OUT/java.$label.b.fq" "$OUT/rust.$label.b.fq"
  cmp "$OUT/java.$label.hist.tsv" "$OUT/rust.$label.hist.tsv"
}

run_paired_case() {
  local label="$1"
  shift
  local args=("$@")

  printf 'Running Java BBNorm quality alias paired case %s...\n' "$label"
  java -Xmx1g -cp "$CP" jgi.KmerNormalize \
    "in=$R1_Q33" \
    "in2=$R2_Q64" \
    "${COMMON[@]}" \
    "${args[@]}" \
    "out=$OUT/java.$label.1.fq" \
    "out2=$OUT/java.$label.2.fq" \
    "hist=$OUT/java.$label.hist.tsv" \
    >"$OUT/java.$label.stdout.log" 2>"$OUT/java.$label.stderr.log"

  printf 'Running Rust bbnorm-rs quality alias paired case %s...\n' "$label"
  target/debug/bbnorm-rs \
    "in=$R1_Q33" \
    "in2=$R2_Q64" \
    "${COMMON[@]}" \
    "${args[@]}" \
    "out=$OUT/rust.$label.1.fq" \
    "out2=$OUT/rust.$label.2.fq" \
    "hist=$OUT/rust.$label.hist.tsv" \
    >"$OUT/rust.$label.stdout.log" 2>"$OUT/rust.$label.stderr.log"

  cmp "$OUT/java.$label.1.fq" "$OUT/rust.$label.1.fq"
  cmp "$OUT/java.$label.2.fq" "$OUT/rust.$label.2.fq"
  cmp "$OUT/java.$label.hist.tsv" "$OUT/rust.$label.hist.tsv"
}

run_interleaved_case() {
  local label="$1"
  local interleaved_arg="$2"
  shift 2
  local args=("$@")

  printf 'Running Java BBNorm quality alias interleaved case %s...\n' "$label"
  java -Xmx1g -cp "$CP" jgi.KmerNormalize \
    "in=$INTERLEAVED" \
    "$interleaved_arg" \
    "${COMMON[@]}" \
    "${args[@]}" \
    "out=$OUT/java.$label.fq" \
    "hist=$OUT/java.$label.hist.tsv" \
    >"$OUT/java.$label.stdout.log" 2>"$OUT/java.$label.stderr.log"

  printf 'Running Rust bbnorm-rs quality alias interleaved case %s...\n' "$label"
  target/debug/bbnorm-rs \
    "in=$INTERLEAVED" \
    "$interleaved_arg" \
    "${COMMON[@]}" \
    "${args[@]}" \
    "out=$OUT/rust.$label.fq" \
    "hist=$OUT/rust.$label.hist.tsv" \
    >"$OUT/rust.$label.stdout.log" 2>"$OUT/rust.$label.stderr.log"

  cmp "$OUT/java.$label.fq" "$OUT/rust.$label.fq"
  cmp "$OUT/java.$label.hist.tsv" "$OUT/rust.$label.hist.tsv"
}

run_case ascii64 "$Q64" "ascii=64"
run_case asciiin64 "$Q64" "asciiin=64"
run_case asciiout64 "$Q33" "asciiout=64" "qin=33"
run_case qinauto "$Q64" "qin=auto"
run_case qauto "$Q64" "qauto=t"
run_list_case qinauto_list_q33_q64 "$Q33,$Q64" "qin=auto"
run_list_case qauto_list_q64_q33 "$Q64,$Q33" "qauto=t"
run_paired_case qinauto_paired_q33_q64 "qin=auto"
run_interleaved_case qauto_interleaved_q33_q64 "interleaved=t" "qauto=t"
run_interleaved_case qinauto_auto_interleaved_q33_q64 "interleaved=auto" "qin=auto"

printf 'Quality alias representative parity passed, including mixed q33/q64 input-list, paired, and interleaved auto-quality cases. Logs and outputs: %s\n' "$OUT"

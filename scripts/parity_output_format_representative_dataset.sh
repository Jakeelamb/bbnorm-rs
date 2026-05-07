#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT"

CP="vendor/BBTools-master/current"
OUT="${1:-tmp/representative_output_format_parity}"
mkdir -p "$OUT"
rm -f "$OUT"/*

FASTA="$OUT/representative_output.fa"
FASTQ="$OUT/representative_output.fq"
LONG_FASTQ="$OUT/representative_output_wrap.fq"
cat > "$FASTA" <<'FASTA'
>seq1
ACGTNN
>seq2
TTTTAC
FASTA
cat > "$FASTQ" <<'FASTQ'
@seq1
ACGTNN
+
IIIIII
FASTQ
LONG_BASES="$(printf 'ACGT%.0s' {1..25})"
LONG_QUALS="$(printf 'I%.0s' {1..100})"
cat > "$LONG_FASTQ" <<FASTQ
@long
$LONG_BASES
+
$LONG_QUALS
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
  local label="$1" input="$2" output_name="$3"
  shift 3
  local args=("$@")

  printf 'Running Java BBNorm output-format case %s...\n' "$label"
  java -Xmx1g -cp "$CP" jgi.KmerNormalize \
    "in=$input" \
    "${COMMON[@]}" \
    "${args[@]}" \
    "out=$OUT/java.$label.$output_name" \
    >"$OUT/java.$label.stdout.log" 2>"$OUT/java.$label.stderr.log"

  printf 'Running Rust bbnorm-rs output-format case %s...\n' "$label"
  target/debug/bbnorm-rs \
    "in=$input" \
    "${COMMON[@]}" \
    "${args[@]}" \
    "out=$OUT/rust.$label.$output_name" \
    >"$OUT/rust.$label.stdout.log" 2>"$OUT/rust.$label.stderr.log"

  cmp "$OUT/java.$label.$output_name" "$OUT/rust.$label.$output_name"
}

run_case fasta_to_fastq_default "$FASTA" keep.fq
run_case fasta_to_fastq_qout64 "$FASTA" keep.fq "qout=64"
run_case fasta_to_fastq_fakefasta20 "$FASTA" keep.fq "fakefastaquality=20"
run_case fasta_to_fastq_qfake15 "$FASTA" keep.fq "qfake=15"
run_case fasta_to_fasta "$FASTA" keep.fa
run_case fasta_to_unknown_txt "$FASTA" keep.txt
run_case fastq_to_fasta "$FASTQ" keep.fa
run_case fastq_to_fasta_default_wrap "$LONG_FASTQ" keep.fa
run_case fastq_to_fasta_fastawrap20 "$LONG_FASTQ" keep.fa "fastawrap=20"
run_case fastq_to_fasta_wrap_alias20 "$LONG_FASTQ" keep.fa "wrap=20"
run_case fastq_to_fasta_fastawrap0 "$LONG_FASTQ" keep.fa "fastawrap=0"

printf 'Output-format/wrap representative parity passed. Logs and outputs: %s\n' "$OUT"

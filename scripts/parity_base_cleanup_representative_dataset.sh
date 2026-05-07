#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT"

CP="vendor/BBTools-master/current"
OUT="${1:-tmp/representative_base_cleanup_parity}"
mkdir -p "$OUT"
rm -f "$OUT"/*

write_fastq() {
  local path="$1" id="$2" bases="$3"
  local qualities
  qualities=$(printf '%*s' "${#bases}" '' | tr ' ' 'I')
  printf '@%s\n%s\n+\n%s\n' "$id" "$bases" "$qualities" > "$path"
}

IUPAC="$OUT/representative_base_iupac.fq"
DDX="$OUT/representative_base_dotdashx.fq"
JUNK="$OUT/representative_base_junk.fq"
COMBO="$OUT/representative_base_combo.fq"
LOWER="$OUT/representative_base_lower.fq"
write_fastq "$IUPAC" iupac 'acgtuURYSWKMBDHVNn'
write_fastq "$DDX" ddx 'ACGT.-XxNn'
write_fastq "$JUNK" junk 'ACGT?ZN'
write_fastq "$COMBO" combo 'acgtuUnN.-XxRrYy'
write_fastq "$LOWER" lower 'acgtnu'

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
  "qin=33"
)

cargo build --quiet

run_case() {
  local label="$1" input="$2"
  shift 2
  local args=("$@")

  printf 'Running Java BBNorm base-cleanup case %s...\n' "$label"
  java -Xmx1g -cp "$CP" jgi.KmerNormalize \
    "in=$input" \
    "${COMMON[@]}" \
    "${args[@]}" \
    "out=$OUT/java.$label.keep.fq" \
    >"$OUT/java.$label.stdout.log" 2>"$OUT/java.$label.stderr.log"

  printf 'Running Rust bbnorm-rs base-cleanup case %s...\n' "$label"
  target/debug/bbnorm-rs \
    "in=$input" \
    "${COMMON[@]}" \
    "${args[@]}" \
    "out=$OUT/rust.$label.keep.fq" \
    >"$OUT/rust.$label.stdout.log" 2>"$OUT/rust.$label.stderr.log"

  cmp "$OUT/java.$label.keep.fq" "$OUT/rust.$label.keep.fq"
}

run_case utot "$IUPAC" "utot=t"
run_case touppercase "$IUPAC" "touppercase=t"
run_case tuc_alias "$IUPAC" "tuc=t"
run_case lowercaseton "$IUPAC" "lowercaseton=t"
run_case lctn_alias "$IUPAC" "lctn=t"
run_case iupacton "$IUPAC" "iupacton=t"
run_case undefinedton_alias "$IUPAC" "undefinedton=t"
run_case dotdashxton "$DDX" "dotdashxton=t"
run_case fixjunk "$JUNK" "fixjunk=t"
run_case ignorejunk "$JUNK" "ignorejunk=t"
run_case flagjunk "$JUNK" "flagjunk=t"
run_case tossjunk "$JUNK" "tossjunk=t"
run_case junk_flag "$JUNK" "junk=flag"
run_case junk_discard "$JUNK" "junk=discard"
run_case crashjunk_false "$JUNK" "crashjunk=f"
run_case failjunk_false "$JUNK" "failjunk=f"
run_case junk_crash_valid "$LOWER" "junk=crash"
run_case junk_fail_valid "$LOWER" "junk=fail"
run_case lowercaseton_changequality_false "$LOWER" "lowercaseton=t" "changequality=f"
run_case dotdashxton_iupacton "$COMBO" "dotdashxton=t" "iupacton=t"
run_case dotdashxton_iupacton_lctn "$COMBO" "dotdashxton=t" "iupacton=t" "lowercaseton=t"
run_case junk_iupacton "$COMBO" "junk=iupacton"

printf 'Base-cleanup representative parity passed. Logs and outputs: %s\n' "$OUT"

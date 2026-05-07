#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT"

CP="vendor/BBTools-master/current"
OUT="${1:-tmp/representative_quality_change_controls_parity}"
mkdir -p "$OUT"
rm -f "$OUT"/*

READS="$OUT/representative_quality_controls.fq"
cat > "$READS" <<'FASTQ'
@quality_controls
ACGTNNACGT
+
!#I~IIIIII
FASTQ

COMMON=(
  "in=$READS"
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
  local label="$1"
  shift
  local args=("$@")

  printf 'Running Java BBNorm quality-control case %s...\n' "$label"
  java -Xmx1g -cp "$CP" jgi.KmerNormalize \
    "${COMMON[@]}" \
    "${args[@]}" \
    "out=$OUT/java.$label.keep.fq" \
    >"$OUT/java.$label.stdout.log" 2>"$OUT/java.$label.stderr.log"

  printf 'Running Rust bbnorm-rs quality-control case %s...\n' "$label"
  target/debug/bbnorm-rs \
    "${COMMON[@]}" \
    "${args[@]}" \
    "out=$OUT/rust.$label.keep.fq" \
    >"$OUT/rust.$label.stdout.log" 2>"$OUT/rust.$label.stderr.log"

  cmp "$OUT/java.$label.keep.fq" "$OUT/rust.$label.keep.fq"
}

run_case changequality_false "changequality=f"
run_case cq_false_alias "cq=f"
run_case ignorebadquality_true "ignorebadquality=t"
run_case ibq_true "ibq=t"
run_case min5_max30 "mincalledquality=5" "maxcalledquality=30"

printf 'Quality change-control representative parity passed. Logs and outputs: %s\n' "$OUT"

#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT"

CP="vendor/BBTools-master/current"
OUT="${1:-tmp/representative_qin64_parity}"
mkdir -p "$OUT"
rm -f "$OUT"/*

READS="$OUT/representative_qin64.fq"
cat > "$READS" <<'FASTQ'
@qin64
ACGTNNACGT
+
@Bh|hhhhhh
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
  "qin=64"
)

cargo build --quiet

printf 'Running Java BBNorm qin=64 representative case...\n'
java -Xmx1g -cp "$CP" jgi.KmerNormalize \
  "${COMMON[@]}" \
  "out=$OUT/java.keep.fq" \
  >"$OUT/java.stdout.log" 2>"$OUT/java.stderr.log"

printf 'Running Rust bbnorm-rs qin=64 representative case...\n'
target/debug/bbnorm-rs \
  "${COMMON[@]}" \
  "out=$OUT/rust.keep.fq" \
  >"$OUT/rust.stdout.log" 2>"$OUT/rust.stderr.log"

cmp "$OUT/java.keep.fq" "$OUT/rust.keep.fq"
python - "$OUT/rust.keep.fq" <<'PY'
import pathlib
import sys
quality = pathlib.Path(sys.argv[1]).read_text().splitlines()[3]
expected = "##IS!!IIII"
if quality != expected:
    raise SystemExit(f"unexpected qin=64 normalized qualities: {quality!r} != {expected!r}")
PY

printf 'qin=64 representative parity passed. Logs and outputs: %s\n' "$OUT"

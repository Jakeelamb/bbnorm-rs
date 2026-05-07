#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT"

CP="vendor/BBTools-master/current"
OUT="${1:-tmp/representative_sampling_options_fallback}"
mkdir -p "$OUT"
rm -f "$OUT"/*

INPUT="$OUT/sampling_guard_input.fq"
cat > "$INPUT" <<'FASTQ'
@sampling_guard
ACGTACGT
+
IIIIIIII
FASTQ

COMMON=(
  "in=$INPUT"
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

java -Xmx1g -cp "$CP" jgi.KmerNormalize \
  "${COMMON[@]}" \
  "out=$OUT/java.base.keep.fq" \
  >"$OUT/java.base.stdout.log" 2>"$OUT/java.base.stderr.log"

for opt in sampleoutput=1 readsample=1 kmersample=1 samplerate=0.5 sample=0.5 sampleseed=1 seed=1; do
  safe_opt="${opt//=/_}"
  printf 'Confirming vendored Java BBNorm rejects %s...\n' "$opt"
  if java -Xmx1g -cp "$CP" jgi.KmerNormalize \
    "${COMMON[@]}" \
    "$opt" \
    "out=$OUT/java.$safe_opt.keep.fq" \
    >"$OUT/java.$safe_opt.stdout.log" 2>"$OUT/java.$safe_opt.stderr.log"; then
    echo "Expected vendored Java KmerNormalize to reject $opt, but it succeeded." >&2
    exit 1
  fi
  grep -q 'Unknown parameter' "$OUT/java.$safe_opt.stderr.log"
  test ! -e "$OUT/java.$safe_opt.keep.fq"

  printf 'Confirming Rust bbnorm-rs ignores %s and keeps working output...\n' "$opt"
  target/debug/bbnorm-rs \
    "${COMMON[@]}" \
    "$opt" \
    "out=$OUT/rust.$safe_opt.keep.fq" \
    >"$OUT/rust.$safe_opt.stdout.log" 2>"$OUT/rust.$safe_opt.stderr.log"
  grep -q 'Rust ignores it' "$OUT/rust.$safe_opt.stderr.log"
  cmp "$OUT/java.base.keep.fq" "$OUT/rust.$safe_opt.keep.fq"
done

printf 'Sampling option fallback parity passed. Logs: %s\n' "$OUT"

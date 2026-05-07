#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT"

CP="vendor/BBTools-master/current"
OUT="${1:-tmp/representative_peak_short_aliases}"
mkdir -p "$OUT"
rm -f "$OUT"/*

READS="$OUT/representative.fq"
python - <<'PY' > "$READS"
seq = "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT"
qual = "I" * len(seq)
for i in range(20):
    print(f"@rep{i}")
    print(seq)
    print("+")
    print(qual)
PY

BASE=(
  "in=$READS"
  "passes=1"
  "keepall=t"
  "k=5"
  "minq=0"
  "minprob=0"
  "min=0"
  "minkmers=1"
  "target=999999999"
  "max=999999999"
  "threads=1"
  "overwrite=t"
  "bits=32"
  "maxpeakcount=8"
)
LONG_PEAK=(
  "minheight=1"
  "minvolume=1"
  "minwidth=1"
  "minpeak=1"
  "maxpeak=100"
)
SHORT_PEAK=(
  "h=1"
  "v=1"
  "w=1"
  "minp=1"
  "maxp=100"
)

cargo build --quiet

printf 'Running Java BBNorm representative peak long-option baseline...\n'
java -Xmx1g -cp "$CP" jgi.KmerNormalize \
  "${BASE[@]}" \
  "${LONG_PEAK[@]}" \
  "out=$OUT/java.long.keep.fq" \
  "hist=$OUT/java.long.hist.tsv" \
  "peaks=$OUT/java.long.peaks.tsv" \
  "peaksout=$OUT/java.long.peaksout.tsv" \
  >"$OUT/java.long.stdout.log" 2>"$OUT/java.long.stderr.log"

printf 'Confirming vendored Java rejects CallPeaks short aliases in KmerNormalize...\n'
if java -Xmx1g -cp "$CP" jgi.KmerNormalize \
  "${BASE[@]}" \
  "${SHORT_PEAK[@]}" \
  "out=$OUT/java.short.keep.fq" \
  "hist=$OUT/java.short.hist.tsv" \
  "peaks=$OUT/java.short.peaks.tsv" \
  "peaksout=$OUT/java.short.peaksout.tsv" \
  >"$OUT/java.short.stdout.log" 2>"$OUT/java.short.stderr.log"; then
  printf 'Expected vendored Java KmerNormalize to reject short peak aliases, but it succeeded.\n' >&2
  exit 1
fi

grep -Eq 'Unknown parameter|Unknown argument|Could not find|Exception' "$OUT/java.short.stderr.log"

printf 'Running Rust bbnorm-rs representative peak long-option baseline...\n'
target/debug/bbnorm-rs \
  "${BASE[@]}" \
  "${LONG_PEAK[@]}" \
  "out=$OUT/rust.long.keep.fq" \
  "hist=$OUT/rust.long.hist.tsv" \
  "peaks=$OUT/rust.long.peaks.tsv" \
  "peaksout=$OUT/rust.long.peaksout.tsv" \
  >"$OUT/rust.long.stdout.log" 2>"$OUT/rust.long.stderr.log"

printf 'Running Rust bbnorm-rs representative peak short aliases...\n'
target/debug/bbnorm-rs \
  "${BASE[@]}" \
  "${SHORT_PEAK[@]}" \
  "out=$OUT/rust.short.keep.fq" \
  "hist=$OUT/rust.short.hist.tsv" \
  "peaks=$OUT/rust.short.peaks.tsv" \
  "peaksout=$OUT/rust.short.peaksout.tsv" \
  >"$OUT/rust.short.stdout.log" 2>"$OUT/rust.short.stderr.log"

cmp "$OUT/java.long.keep.fq" "$OUT/rust.long.keep.fq"
cmp "$OUT/java.long.hist.tsv" "$OUT/rust.long.hist.tsv"
cmp "$OUT/java.long.peaks.tsv" "$OUT/rust.long.peaks.tsv"
cmp "$OUT/java.long.peaksout.tsv" "$OUT/rust.long.peaksout.tsv"

cmp "$OUT/rust.long.keep.fq" "$OUT/rust.short.keep.fq"
cmp "$OUT/rust.long.hist.tsv" "$OUT/rust.short.hist.tsv"
cmp "$OUT/rust.long.peaks.tsv" "$OUT/rust.short.peaks.tsv"
cmp "$OUT/rust.long.peaksout.tsv" "$OUT/rust.short.peaksout.tsv"

grep -q $'17\t20\t23\t36\t36' "$OUT/rust.short.peaks.tsv"

printf 'Representative peak short-alias Rust support passed. Logs and outputs: %s\n' "$OUT"

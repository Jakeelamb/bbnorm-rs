#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT"

CP="vendor/BBTools-master/current"
OUT="${1:-tmp/representative_peaks_parity}"
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

COMMON=(
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
  "minheight=1"
  "minvolume=1"
  "minwidth=1"
  "minpeak=1"
  "maxpeak=100"
  "maxpeakcount=8"
)

cargo build --quiet

printf 'Running Java BBNorm representative peak case...\n'
java -Xmx1g -cp "$CP" jgi.KmerNormalize \
  "${COMMON[@]}" \
  "out=$OUT/java.keep.fq" \
  "hist=$OUT/java.hist.tsv" \
  "peaks=$OUT/java.peaks.tsv" \
  "peaksout=$OUT/java.peaksout.tsv" \
  >"$OUT/java.stdout.log" 2>"$OUT/java.stderr.log"

grep -q $'#start\tcenter\tstop\tmax\tvolume' "$OUT/java.peaks.tsv"
grep -q $'17\t20\t23\t36\t36' "$OUT/java.peaks.tsv"

printf 'Running Rust bbnorm-rs representative peak case...\n'
target/debug/bbnorm-rs \
  "${COMMON[@]}" \
  "out=$OUT/rust.keep.fq" \
  "hist=$OUT/rust.hist.tsv" \
  "peaks=$OUT/rust.peaks.tsv" \
  "peaksout=$OUT/rust.peaksout.tsv" \
  >"$OUT/rust.stdout.log" 2>"$OUT/rust.stderr.log"

cmp "$OUT/java.keep.fq" "$OUT/rust.keep.fq"
cmp "$OUT/java.hist.tsv" "$OUT/rust.hist.tsv"
cmp "$OUT/java.peaks.tsv" "$OUT/rust.peaks.tsv"
cmp "$OUT/java.peaksout.tsv" "$OUT/rust.peaksout.tsv"

printf 'Representative peak parity passed. Logs and outputs: %s\n' "$OUT"

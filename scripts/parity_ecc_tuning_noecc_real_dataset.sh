#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT"

DATA1="vendor/BBTools-master/resources/sample1.fq.gz"
DATA2="vendor/BBTools-master/resources/sample2.fq.gz"
CP="vendor/BBTools-master/current"
OUT="${1:-tmp/real_ecc_tuning_noecc_parity}"
mkdir -p "$OUT"
rm -f "$OUT"/*

COMMON=(
  "in=$DATA1"
  "in2=$DATA2"
  "passes=1"
  "keepall=t"
  "lowbindepth=1"
  "highbindepth=1"
  "k=31"
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

ECC_TUNING=(
  "ecclimit=3"
  "eccmaxqual=127"
  "errorcorrectratio=140"
  "echighthresh=22"
  "eclowthresh=2"
  "suflen=3"
  "prefixlen=3"
  "cfl=t"
  "cfr=f"
)

cargo build --quiet

printf 'Running Java BBNorm no-ECC tuning-control case on paired phiX...\n'
java -Xmx1g -cp "$CP" jgi.KmerNormalize \
  "${COMMON[@]}" "${ECC_TUNING[@]}" \
  "out=$OUT/java.keep1.fq" "out2=$OUT/java.keep2.fq" \
  "outlow=$OUT/java.low1.fq" "outlow2=$OUT/java.low2.fq" \
  "outmid=$OUT/java.mid1.fq" "outmid2=$OUT/java.mid2.fq" \
  "outhigh=$OUT/java.high1.fq" "outhigh2=$OUT/java.high2.fq" \
  "hist=$OUT/java.hist.tsv" "rhist=$OUT/java.rhist.tsv" \
  >"$OUT/java.stdout.log" 2>"$OUT/java.stderr.log"

printf 'Running Rust bbnorm-rs no-ECC tuning-control case on paired phiX...\n'
target/debug/bbnorm-rs \
  "${COMMON[@]}" "${ECC_TUNING[@]}" \
  "out=$OUT/rust.keep1.fq" "out2=$OUT/rust.keep2.fq" \
  "outlow=$OUT/rust.low1.fq" "outlow2=$OUT/rust.low2.fq" \
  "outmid=$OUT/rust.mid1.fq" "outmid2=$OUT/rust.mid2.fq" \
  "outhigh=$OUT/rust.high1.fq" "outhigh2=$OUT/rust.high2.fq" \
  "hist=$OUT/rust.hist.tsv" "rhist=$OUT/rust.rhist.tsv" \
  >"$OUT/rust.stdout.log" 2>"$OUT/rust.stderr.log"

for suffix in keep1.fq keep2.fq low1.fq low2.fq mid1.fq mid2.fq high1.fq high2.fq hist.tsv rhist.tsv; do
  cmp "$OUT/java.$suffix" "$OUT/rust.$suffix"
done

REP="$OUT/representative.fq"
cat > "$REP" <<'FASTQ'
@malformed
ACGTACGT
+
IIIIIIII
FASTQ

for option in "ecclimit=abc" "eccmaxqual=abc" "ecr=abc" "echthresh=abc" "eclt=abc" "suflen=abc" "prelen=abc"; do
  stem="${option%%=*}"
  printf 'Checking malformed ECC tuning option %s is rejected by Java and Rust...\n' "$option"
  if java -Xmx1g -cp "$CP" jgi.KmerNormalize \
    "in=$REP" "passes=1" "keepall=t" "k=3" "minq=0" "minprob=0" \
    "min=0" "minkmers=1" "target=999999999" "max=999999999" \
    "threads=1" "overwrite=t" "bits=32" "$option" \
    "out=$OUT/java.$stem.keep.fq" \
    >"$OUT/java.$stem.stdout.log" 2>"$OUT/java.$stem.stderr.log"; then
    echo "Expected Java BBNorm to reject malformed $option" >&2
    exit 1
  fi

  if target/debug/bbnorm-rs \
    "in=$REP" "passes=1" "keepall=t" "k=3" "minq=0" "minprob=0" \
    "min=0" "minkmers=1" "target=999999999" "max=999999999" \
    "threads=1" "overwrite=t" "bits=32" "$option" \
    "out=$OUT/rust.$stem.keep.fq" \
    >"$OUT/rust.$stem.stdout.log" 2>"$OUT/rust.$stem.stderr.log"; then
    echo "Expected bbnorm-rs to reject malformed $option" >&2
    exit 1
  fi
  test ! -e "$OUT/rust.$stem.keep.fq"
done

printf 'No-ECC tuning-control parity passed. Logs and outputs: %s\n' "$OUT"

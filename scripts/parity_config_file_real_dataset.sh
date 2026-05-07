#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT"

DATA1="vendor/BBTools-master/resources/sample1.fq.gz"
DATA2="vendor/BBTools-master/resources/sample2.fq.gz"
CP="vendor/BBTools-master/current"
OUT="${1:-tmp/real_config_file_parity}"
mkdir -p "$OUT"
rm -f "$OUT"/*

CONFIG1="$OUT/input.config"
CONFIG2="$OUT/normalization.config"
cat > "$CONFIG1" <<EOF_CFG
# BBTools-style one-argument-per-line config
in=$DATA1
in2=$DATA2
passes=1
keepall=t
lowbindepth=1
highbindepth=1
EOF_CFG
cat > "$CONFIG2" <<EOF_CFG
# A second file exercises comma-separated config lists.
k=31
minq=0
minprob=0
min=0
minkmers=1
target=999999999
max=999999999
threads=1
overwrite=t
bits=32
EOF_CFG

cargo build --quiet

printf 'Running Java BBNorm through comma-separated config files on paired phiX...\n'
java -Xmx1g -cp "$CP" jgi.KmerNormalize \
  "config=$CONFIG1,$CONFIG2" \
  "out=$OUT/java.keep1.fq" "out2=$OUT/java.keep2.fq" \
  "outlow=$OUT/java.low1.fq" "outlow2=$OUT/java.low2.fq" \
  "outmid=$OUT/java.mid1.fq" "outmid2=$OUT/java.mid2.fq" \
  "outhigh=$OUT/java.high1.fq" "outhigh2=$OUT/java.high2.fq" \
  "hist=$OUT/java.hist.tsv" "rhist=$OUT/java.rhist.tsv" \
  >"$OUT/java.stdout.log" 2>"$OUT/java.stderr.log"

printf 'Running Rust bbnorm-rs through comma-separated config files on paired phiX...\n'
target/debug/bbnorm-rs \
  "config=$CONFIG1,$CONFIG2" \
  "out=$OUT/rust.keep1.fq" "out2=$OUT/rust.keep2.fq" \
  "outlow=$OUT/rust.low1.fq" "outlow2=$OUT/rust.low2.fq" \
  "outmid=$OUT/rust.mid1.fq" "outmid2=$OUT/rust.mid2.fq" \
  "outhigh=$OUT/rust.high1.fq" "outhigh2=$OUT/rust.high2.fq" \
  "hist=$OUT/rust.hist.tsv" "rhist=$OUT/rust.rhist.tsv" \
  >"$OUT/rust.stdout.log" 2>"$OUT/rust.stderr.log"

grep -q 'config=.*expanded into 16 BBTools-style argument line' "$OUT/rust.stderr.log"

for suffix in keep1.fq keep2.fq low1.fq low2.fq mid1.fq mid2.fq high1.fq high2.fq hist.tsv rhist.tsv; do
  cmp "$OUT/java.$suffix" "$OUT/rust.$suffix"
done

printf 'Config-file parity passed. Rust expands BBTools config files and matches local Java baseline. Logs: %s\n' "$OUT"

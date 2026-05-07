#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT"

DATA1="vendor/BBTools-master/resources/sample1.fq.gz"
DATA2="vendor/BBTools-master/resources/sample2.fq.gz"
CP="vendor/BBTools-master/current"
OUT="${1:-tmp/real_rename_parity}"
mkdir -p "$OUT"
rm -f "$OUT"/*

COMMON=(
  "passes=1"
  "keepall=t"
  "rename=t"
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

cargo build --quiet

printf 'Running Java BBNorm renamed paired keep-all case...\n'
JAVA_START=$(python - <<'PYTIME'
import time
print(time.perf_counter_ns())
PYTIME
)
java -Xmx1g -cp "$CP" jgi.KmerNormalize \
  "in=$DATA1" "in2=$DATA2" "${COMMON[@]}" \
  "out=$OUT/java.paired.keep1.fq" "out2=$OUT/java.paired.keep2.fq" \
  >"$OUT/java.paired.stdout.log" 2>"$OUT/java.paired.stderr.log"
JAVA_END=$(python - <<'PYTIME'
import time
print(time.perf_counter_ns())
PYTIME
)
JAVA_PAIRED_MS=$(((JAVA_END - JAVA_START) / 1000000))

printf 'Running Rust bbnorm-rs renamed paired keep-all case...\n'
RUST_START=$(python - <<'PYTIME'
import time
print(time.perf_counter_ns())
PYTIME
)
target/debug/bbnorm-rs \
  "in=$DATA1" "in2=$DATA2" "${COMMON[@]}" \
  "out=$OUT/rust.paired.keep1.fq" "out2=$OUT/rust.paired.keep2.fq" \
  >"$OUT/rust.paired.stdout.log" 2>"$OUT/rust.paired.stderr.log"
RUST_END=$(python - <<'PYTIME'
import time
print(time.perf_counter_ns())
PYTIME
)
RUST_PAIRED_MS=$(((RUST_END - RUST_START) / 1000000))

cmp "$OUT/java.paired.keep1.fq" "$OUT/rust.paired.keep1.fq"
cmp "$OUT/java.paired.keep2.fq" "$OUT/rust.paired.keep2.fq"

INTERLEAVED="$OUT/interleaved.fq"
python - "$DATA1" "$DATA2" "$INTERLEAVED" <<'PYBUILD'
import gzip
import sys

r1_path, r2_path, out_path = sys.argv[1:]
with gzip.open(r1_path, 'rt') as r1, gzip.open(r2_path, 'rt') as r2, open(out_path, 'w') as out:
    while True:
        rec1 = [r1.readline() for _ in range(4)]
        rec2 = [r2.readline() for _ in range(4)]
        if not rec1[0] and not rec2[0]:
            break
        if not all(rec1) or not all(rec2):
            raise SystemExit('paired phiX fixture is truncated or has uneven record counts')
        out.writelines(rec1)
        out.writelines(rec2)
PYBUILD

printf 'Running Java BBNorm renamed interleaved keep-all case...\n'
JAVA_START=$(python - <<'PYTIME'
import time
print(time.perf_counter_ns())
PYTIME
)
java -Xmx1g -cp "$CP" jgi.KmerNormalize \
  "in=$INTERLEAVED" "interleaved=t" "${COMMON[@]}" \
  "out=$OUT/java.interleaved.keep.fq" \
  >"$OUT/java.interleaved.stdout.log" 2>"$OUT/java.interleaved.stderr.log"
JAVA_END=$(python - <<'PYTIME'
import time
print(time.perf_counter_ns())
PYTIME
)
JAVA_INTERLEAVED_MS=$(((JAVA_END - JAVA_START) / 1000000))

printf 'Running Rust bbnorm-rs renamed interleaved keep-all case...\n'
RUST_START=$(python - <<'PYTIME'
import time
print(time.perf_counter_ns())
PYTIME
)
target/debug/bbnorm-rs \
  "in=$INTERLEAVED" "interleaved=t" "${COMMON[@]}" \
  "out=$OUT/rust.interleaved.keep.fq" \
  >"$OUT/rust.interleaved.stdout.log" 2>"$OUT/rust.interleaved.stderr.log"
RUST_END=$(python - <<'PYTIME'
import time
print(time.perf_counter_ns())
PYTIME
)
RUST_INTERLEAVED_MS=$(((RUST_END - RUST_START) / 1000000))

cmp "$OUT/java.interleaved.keep.fq" "$OUT/rust.interleaved.keep.fq"

CLEAN="ACGTTGCATGTCAGTACCGTAACGTTGCA"
MUTANT="ACGTTGCATGTCAGAACCGTAACGTTGCA"
QUAL="IIIIIIIIIIIIIIIIIIIIIIIIIIIII"
ECC_RENAME="$OUT/ecc_rename.fq"
{
  for i in $(seq 1 30); do
    printf '@clean%s\n%s\n+\n%s\n' "$i" "$CLEAN" "$QUAL"
  done
  printf '@mutant\n%s\n+\n%s\n' "$MUTANT" "$QUAL"
} > "$ECC_RENAME"

printf 'Running Java BBNorm renamed ECC single-end case...\n'
java -Xmx1g -cp "$CP" jgi.KmerNormalize \
  "in=$ECC_RENAME" "passes=1" "keepall=t" "rename=t" "ecc=t" "ecco=f" \
  "k=7" "minq=0" "minprob=0" "min=0" "minkmers=1" \
  "target=999999999" "max=999999999" "threads=1" "overwrite=t" "bits=32" \
  "out=$OUT/java.ecc.rename.fq" \
  >"$OUT/java.ecc.rename.stdout.log" 2>"$OUT/java.ecc.rename.stderr.log"

printf 'Running Rust bbnorm-rs renamed ECC single-end case...\n'
target/debug/bbnorm-rs \
  "in=$ECC_RENAME" "passes=1" "keepall=t" "rename=t" "ecc=t" "ecco=f" \
  "k=7" "minq=0" "minprob=0" "min=0" "minkmers=1" \
  "target=999999999" "max=999999999" "threads=1" "overwrite=t" "bits=32" \
  "out=$OUT/rust.ecc.rename.fq" \
  >"$OUT/rust.ecc.rename.stdout.log" 2>"$OUT/rust.ecc.rename.stderr.log"

cmp "$OUT/java.ecc.rename.fq" "$OUT/rust.ecc.rename.fq"
grep -q 'e1=0' "$OUT/rust.ecc.rename.fq"

ECC_RENAME_R1="$OUT/ecc_rename_r1.fq"
ECC_RENAME_R2="$OUT/ecc_rename_r2.fq"
{
  for i in $(seq 1 30); do
    printf '@clean%s/1\n%s\n+\n%s\n' "$i" "$CLEAN" "$QUAL"
  done
  printf '@mutant/1\n%s\n+\n%s\n' "$MUTANT" "$QUAL"
} > "$ECC_RENAME_R1"
{
  for i in $(seq 1 30); do
    printf '@clean%s/2\n%s\n+\n%s\n' "$i" "$CLEAN" "$QUAL"
  done
  printf '@mutant/2\n%s\n+\n%s\n' "$CLEAN" "$QUAL"
} > "$ECC_RENAME_R2"

printf 'Running Java BBNorm renamed ECC paired case...\n'
java -Xmx1g -cp "$CP" jgi.KmerNormalize \
  "in=$ECC_RENAME_R1" "in2=$ECC_RENAME_R2" "passes=1" "keepall=t" "rename=t" "ecc=t" "ecco=f" \
  "k=7" "minq=0" "minprob=0" "min=0" "minkmers=1" \
  "target=999999999" "max=999999999" "threads=1" "overwrite=t" "bits=32" \
  "out=$OUT/java.ecc.paired.rename1.fq" "out2=$OUT/java.ecc.paired.rename2.fq" \
  >"$OUT/java.ecc.paired.rename.stdout.log" 2>"$OUT/java.ecc.paired.rename.stderr.log"

printf 'Running Rust bbnorm-rs renamed ECC paired case...\n'
target/debug/bbnorm-rs \
  "in=$ECC_RENAME_R1" "in2=$ECC_RENAME_R2" "passes=1" "keepall=t" "rename=t" "ecc=t" "ecco=f" \
  "k=7" "minq=0" "minprob=0" "min=0" "minkmers=1" \
  "target=999999999" "max=999999999" "threads=1" "overwrite=t" "bits=32" \
  "out=$OUT/rust.ecc.paired.rename1.fq" "out2=$OUT/rust.ecc.paired.rename2.fq" \
  >"$OUT/rust.ecc.paired.rename.stdout.log" 2>"$OUT/rust.ecc.paired.rename.stderr.log"

cmp "$OUT/java.ecc.paired.rename1.fq" "$OUT/rust.ecc.paired.rename1.fq"
cmp "$OUT/java.ecc.paired.rename2.fq" "$OUT/rust.ecc.paired.rename2.fq"
grep -q 'e1=0,e2=0' "$OUT/rust.ecc.paired.rename1.fq"
grep -q 'e1=0,e2=0' "$OUT/rust.ecc.paired.rename2.fq"

printf 'Rename paired parity passed. Java: %sms, Rust: %sms\n' "$JAVA_PAIRED_MS" "$RUST_PAIRED_MS"
printf 'Rename interleaved parity passed. Java: %sms, Rust: %sms\n' "$JAVA_INTERLEAVED_MS" "$RUST_INTERLEAVED_MS"
printf 'Rename ECC single-end and paired parity passed with Java-shaped e1/e2 fields.\n'
printf 'Outputs and logs: %s\n' "$OUT"

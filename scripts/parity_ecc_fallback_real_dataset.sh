#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT"

DATA1="vendor/BBTools-master/resources/sample1.fq.gz"
DATA2="vendor/BBTools-master/resources/sample2.fq.gz"
CP="vendor/BBTools-master/current"
OUT="${1:-tmp/real_ecc_behavior}"
mkdir -p "$OUT"
rm -f "$OUT"/*

COMMON=(
  "in=$DATA1"
  "in2=$DATA2"
  "passes=1"
  "keepall=t"
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

printf 'Running Java BBNorm no-ECC baseline on paired phiX...\n'
java -Xmx1g -cp "$CP" jgi.KmerNormalize \
  "${COMMON[@]}" \
  "out=$OUT/java.keep1.fq" "out2=$OUT/java.keep2.fq" \
  "outuncorrected=$OUT/java.unc1.fq" "outuncorrected2=$OUT/java.unc2.fq" \
  >"$OUT/java.stdout.log" 2>"$OUT/java.stderr.log"

printf 'Running Rust bbnorm-rs table-based ecc=t on paired phiX...\n'
target/debug/bbnorm-rs \
  "${COMMON[@]}" \
  "ecc=t" \
  "out=$OUT/rust.keep1.fq" "out2=$OUT/rust.keep2.fq" \
  "outuncorrected=$OUT/rust.unc1.fq" "outuncorrected2=$OUT/rust.unc2.fq" \
  >"$OUT/rust.stdout.log" 2>"$OUT/rust.stderr.log"

if grep -q 'error correction option ecc=t is not implemented yet' "$OUT/rust.stderr.log"; then
  echo "Rust still emitted the old ECC fallback note" >&2
  exit 1
fi
cmp "$OUT/java.keep1.fq" "$OUT/rust.keep1.fq"
cmp "$OUT/java.keep2.fq" "$OUT/rust.keep2.fq"
cmp "$OUT/java.unc1.fq" "$OUT/rust.unc1.fq"
cmp "$OUT/java.unc2.fq" "$OUT/rust.unc2.fq"
test ! -s "$OUT/rust.unc1.fq"
test ! -s "$OUT/rust.unc2.fq"

REP="$OUT/ecc_substitution.fq"
CLEAN="ACGTTGCATGTCAGTACCGTAACGTTGCA"
MUTANT="ACGTTGCATGTCAGAACCGTAACGTTGCA"
QUAL="IIIIIIIIIIIIIIIIIIIIIIIIIIIII"
{
  for i in $(seq 1 30); do
    printf '@clean%s\n%s\n+\n%s\n' "$i" "$CLEAN" "$QUAL"
  done
  printf '@mutant\n%s\n+\n%s\n' "$MUTANT" "$QUAL"
} > "$REP"

printf 'Running Rust bbnorm-rs ecc=t on a representative one-substitution fixture...\n'
target/debug/bbnorm-rs \
  "in=$REP" "passes=1" "keepall=t" "ecc=t" \
  "k=7" "minq=0" "minprob=0" "min=0" "minkmers=1" \
  "target=999999999" "max=999999999" "threads=1" "overwrite=t" "bits=32" \
  "out=$OUT/rust.corrected.fq" \
  >"$OUT/rust.corrected.stdout.log" 2>"$OUT/rust.corrected.stderr.log"

if grep -q "$MUTANT" "$OUT/rust.corrected.fq"; then
  echo "Rust ECC left the mutant read uncorrected" >&2
  exit 1
fi
grep -q "$CLEAN" "$OUT/rust.corrected.fq"
tail -4 "$OUT/rust.corrected.fq" | grep -q "$CLEAN"

printf 'Running Java and Rust bbnorm-rs on a representative uncorrectable ECC fixture...\n'
java -Xmx1g -cp "$CP" jgi.KmerNormalize \
  "in=$REP" "passes=1" "keepall=t" "ecc=t" "ecco=f" "eccmaxqual=0" \
  "k=7" "minq=0" "minprob=0" "min=0" "minkmers=1" \
  "target=999999999" "max=999999999" "threads=1" "overwrite=t" "bits=32" \
  "out=$OUT/java.uncorrectable.keep.fq" "outuncorrected=$OUT/java.uncorrectable.unc.fq" \
  >"$OUT/java.uncorrectable.stdout.log" 2>"$OUT/java.uncorrectable.stderr.log"

target/debug/bbnorm-rs \
  "in=$REP" "passes=1" "keepall=t" "ecc=t" "ecco=f" "eccmaxqual=0" \
  "k=7" "minq=0" "minprob=0" "min=0" "minkmers=1" \
  "target=999999999" "max=999999999" "threads=1" "overwrite=t" "bits=32" \
  "out=$OUT/rust.uncorrectable.keep.fq" "outuncorrected=$OUT/rust.uncorrectable.unc.fq" \
  >"$OUT/rust.uncorrectable.stdout.log" 2>"$OUT/rust.uncorrectable.stderr.log"

cmp "$OUT/java.uncorrectable.keep.fq" "$OUT/rust.uncorrectable.keep.fq"
cmp "$OUT/java.uncorrectable.unc.fq" "$OUT/rust.uncorrectable.unc.fq"
grep -q "$MUTANT" "$OUT/rust.uncorrectable.unc.fq"
tail -4 "$OUT/rust.uncorrectable.unc.fq" | grep -q "$MUTANT"

printf 'Running Java and Rust bbnorm-rs on markuncorrectableerrors=t for the uncorrectable fixture...\n'
java -Xmx1g -cp "$CP" jgi.KmerNormalize \
  "in=$REP" "passes=1" "keepall=t" "ecc=t" "ecco=f" "eccmaxqual=0" "markuncorrectableerrors=t" \
  "k=7" "minq=0" "minprob=0" "min=0" "minkmers=1" \
  "target=999999999" "max=999999999" "threads=1" "overwrite=t" "bits=32" \
  "out=$OUT/java.markuncorrectable.keep.fq" "outuncorrected=$OUT/java.markuncorrectable.unc.fq" \
  >"$OUT/java.markuncorrectable.stdout.log" 2>"$OUT/java.markuncorrectable.stderr.log"

target/debug/bbnorm-rs \
  "in=$REP" "passes=1" "keepall=t" "ecc=t" "ecco=f" "eccmaxqual=0" "markuncorrectableerrors=t" \
  "k=7" "minq=0" "minprob=0" "min=0" "minkmers=1" \
  "target=999999999" "max=999999999" "threads=1" "overwrite=t" "bits=32" \
  "out=$OUT/rust.markuncorrectable.keep.fq" "outuncorrected=$OUT/rust.markuncorrectable.unc.fq" \
  >"$OUT/rust.markuncorrectable.stdout.log" 2>"$OUT/rust.markuncorrectable.stderr.log"

cmp "$OUT/java.markuncorrectable.keep.fq" "$OUT/rust.markuncorrectable.keep.fq"
cmp "$OUT/java.markuncorrectable.unc.fq" "$OUT/rust.markuncorrectable.unc.fq"
tail -1 "$OUT/rust.markuncorrectable.keep.fq" | grep -q '2'
tail -1 "$OUT/rust.markuncorrectable.unc.fq" | grep -q '2'

MULTI_MARK_REP="$OUT/ecc_multi_mark.fq"
MULTI_MARK_MUTANT="ACGTTGCATGTCAAAACCGTAATGTTGCA"
{
  for i in $(seq 1 40); do
    printf '@clean%s\n%s\n+\n%s\n' "$i" "$CLEAN" "$QUAL"
  done
  printf '@multi_mark\n%s\n+\n%s\n' "$MULTI_MARK_MUTANT" "$QUAL"
} > "$MULTI_MARK_REP"

printf 'Running Java and Rust bbnorm-rs on multi-site markerrors ecclimit semantics...\n'
java -Xmx1g -cp "$CP" jgi.KmerNormalize \
  "in=$MULTI_MARK_REP" "passes=1" "keepall=t" "ecc=t" "ecco=f" \
  "markerrors=t" "ecclimit=1" "cfl=t" "cfr=f" \
  "k=7" "minq=0" "minprob=0" "min=0" "minkmers=1" \
  "target=999999999" "max=999999999" "threads=1" "overwrite=t" "bits=32" \
  "out=$OUT/java.multi_mark.fq" \
  >"$OUT/java.multi_mark.stdout.log" 2>"$OUT/java.multi_mark.stderr.log"

target/debug/bbnorm-rs \
  "in=$MULTI_MARK_REP" "passes=1" "keepall=t" "ecc=t" "ecco=f" \
  "markerrors=t" "ecclimit=1" "cfl=t" "cfr=f" \
  "k=7" "minq=0" "minprob=0" "min=0" "minkmers=1" \
  "target=999999999" "max=999999999" "threads=1" "overwrite=t" "bits=32" \
  "out=$OUT/rust.multi_mark.fq" \
  >"$OUT/rust.multi_mark.stdout.log" 2>"$OUT/rust.multi_mark.stderr.log"

cmp "$OUT/java.multi_mark.fq" "$OUT/rust.multi_mark.fq"
tail -4 "$OUT/rust.multi_mark.fq" | grep -q "$MULTI_MARK_MUTANT"
tail -1 "$OUT/rust.multi_mark.fq" | grep -q '2'

PAIR_R1="$OUT/ecc_uncorrectable_r1.fq"
PAIR_R2="$OUT/ecc_uncorrectable_r2.fq"
{
  for i in $(seq 1 30); do
    printf '@clean%s/1\n%s\n+\n%s\n' "$i" "$CLEAN" "$QUAL"
  done
  printf '@mutant/1\n%s\n+\n%s\n' "$MUTANT" "$QUAL"
} > "$PAIR_R1"
{
  for i in $(seq 1 30); do
    printf '@clean%s/2\n%s\n+\n%s\n' "$i" "$CLEAN" "$QUAL"
  done
  printf '@mutant/2\n%s\n+\n%s\n' "$CLEAN" "$QUAL"
} > "$PAIR_R2"

printf 'Running Java and Rust bbnorm-rs on a paired uncorrectable ECC fixture...\n'
java -Xmx1g -cp "$CP" jgi.KmerNormalize \
  "in=$PAIR_R1" "in2=$PAIR_R2" "passes=1" "keepall=t" "ecc=t" "ecco=f" "eccmaxqual=0" \
  "k=7" "minq=0" "minprob=0" "min=0" "minkmers=1" \
  "target=999999999" "max=999999999" "threads=1" "overwrite=t" "bits=32" \
  "out=$OUT/java.paired.uncorrectable.keep1.fq" "out2=$OUT/java.paired.uncorrectable.keep2.fq" \
  "outuncorrected=$OUT/java.paired.uncorrectable.unc1.fq" "outuncorrected2=$OUT/java.paired.uncorrectable.unc2.fq" \
  >"$OUT/java.paired.uncorrectable.stdout.log" 2>"$OUT/java.paired.uncorrectable.stderr.log"

target/debug/bbnorm-rs \
  "in=$PAIR_R1" "in2=$PAIR_R2" "passes=1" "keepall=t" "ecc=t" "ecco=f" "eccmaxqual=0" \
  "k=7" "minq=0" "minprob=0" "min=0" "minkmers=1" \
  "target=999999999" "max=999999999" "threads=1" "overwrite=t" "bits=32" \
  "out=$OUT/rust.paired.uncorrectable.keep1.fq" "out2=$OUT/rust.paired.uncorrectable.keep2.fq" \
  "outuncorrected=$OUT/rust.paired.uncorrectable.unc1.fq" "outuncorrected2=$OUT/rust.paired.uncorrectable.unc2.fq" \
  >"$OUT/rust.paired.uncorrectable.stdout.log" 2>"$OUT/rust.paired.uncorrectable.stderr.log"

cmp "$OUT/java.paired.uncorrectable.keep1.fq" "$OUT/rust.paired.uncorrectable.keep1.fq"
cmp "$OUT/java.paired.uncorrectable.keep2.fq" "$OUT/rust.paired.uncorrectable.keep2.fq"
cmp "$OUT/java.paired.uncorrectable.unc1.fq" "$OUT/rust.paired.uncorrectable.unc1.fq"
cmp "$OUT/java.paired.uncorrectable.unc2.fq" "$OUT/rust.paired.uncorrectable.unc2.fq"
grep -q "$MUTANT" "$OUT/rust.paired.uncorrectable.unc1.fq"
tail -4 "$OUT/rust.paired.uncorrectable.unc1.fq" | grep -q "$MUTANT"
tail -4 "$OUT/rust.paired.uncorrectable.unc2.fq" | grep -q "$CLEAN"

TAIL_MUTANT="ACGTTGCATGTCAGTACCGTAACGTTAC"
TAIL_UNTRIMMED="ACGTTGCATGTCAGTACCGTAACGTTACA"
TAIL_REP="$OUT/ecc_mark_trim.fq"
{
  for i in $(seq 1 30); do
    printf '@clean%s\n%s\n+\n%s\n' "$i" "$CLEAN" "$QUAL"
  done
  printf '@tail_mutant\n%s\n+\n%s\n' "$TAIL_UNTRIMMED" "$QUAL"
} > "$TAIL_REP"

printf 'Running Java and Rust bbnorm-rs on markerrors+trimaftermarking fixture...\n'
java -Xmx1g -cp "$CP" jgi.KmerNormalize \
  "in=$TAIL_REP" "passes=1" "keepall=t" "ecc=t" "ecco=f" \
  "markerrors=t" "trimaftermarking=t" "qtrim=r" "trimq=20" "optitrim=f" \
  "k=7" "minq=0" "minprob=0" "min=0" "minkmers=1" \
  "target=999999999" "max=999999999" "threads=1" "overwrite=t" "bits=32" \
  "out=$OUT/java.marktrim.fq" \
  >"$OUT/java.marktrim.stdout.log" 2>"$OUT/java.marktrim.stderr.log"

target/debug/bbnorm-rs \
  "in=$TAIL_REP" "passes=1" "keepall=t" "ecc=t" "ecco=f" \
  "markerrors=t" "trimaftermarking=t" "qtrim=r" "trimq=20" "optitrim=f" \
  "k=7" "minq=0" "minprob=0" "min=0" "minkmers=1" \
  "target=999999999" "max=999999999" "threads=1" "overwrite=t" "bits=32" \
  "out=$OUT/rust.marktrim.fq" \
  >"$OUT/rust.marktrim.stdout.log" 2>"$OUT/rust.marktrim.stderr.log"

cmp "$OUT/java.marktrim.fq" "$OUT/rust.marktrim.fq"
tail -4 "$OUT/rust.marktrim.fq" | grep -q "$TAIL_MUTANT"
if tail -4 "$OUT/rust.marktrim.fq" | grep -q "$TAIL_UNTRIMMED"; then
  echo "Rust markerrors+trimaftermarking left the tail mutant untrimmed" >&2
  exit 1
fi

printf 'Running Java and Rust bbnorm-rs on correction-staged multipass fixtures...\n'
for case_name in final_only first_only both_passes; do
  case "$case_name" in
    final_only) ECC_STAGE=("ecc1=f" "eccf=t") ;;
    first_only) ECC_STAGE=("ecc1=t" "eccf=f") ;;
    both_passes) ECC_STAGE=("ecc=t") ;;
  esac

  java -Xmx1g -cp "$CP" jgi.KmerNormalize \
    "in=$REP" "passes=2" "keepall=t" "ecco=f" "${ECC_STAGE[@]}" \
    "k=7" "minq=0" "minprob=0" "min=0" "minkmers=1" \
    "target=999999999" "max=999999999" "threads=1" "overwrite=t" "bits=32" \
    "out=$OUT/java.multipass.$case_name.fq" \
    >"$OUT/java.multipass.$case_name.stdout.log" 2>"$OUT/java.multipass.$case_name.stderr.log"

  target/debug/bbnorm-rs \
    "in=$REP" "passes=2" "keepall=t" "ecco=f" "${ECC_STAGE[@]}" \
    "k=7" "minq=0" "minprob=0" "min=0" "minkmers=1" \
    "target=999999999" "max=999999999" "threads=1" "overwrite=t" "bits=32" \
    "out=$OUT/rust.multipass.$case_name.fq" \
    >"$OUT/rust.multipass.$case_name.stdout.log" 2>"$OUT/rust.multipass.$case_name.stderr.log"

  cmp "$OUT/java.multipass.$case_name.fq" "$OUT/rust.multipass.$case_name.fq"
  tail -4 "$OUT/rust.multipass.$case_name.fq" | grep -q "$CLEAN"
done

NOISY_REP="$OUT/ecc_noisy_multipass.fq"
NOISY_MUTANT_A="ACGTTGCATGTCAGAACCGTAACGTTGCA"
NOISY_MUTANT_B="ACGTTGCATGTCAGTACCGTAATGTTGCA"
NOISY_DOUBLE_MUTANT="ACGTTGCATGTCAGAACCGTAATGTTGCA"
NOISY_QUAL="$(printf 'I%.0s' $(seq 1 ${#CLEAN}))"
{
  for i in $(seq 1 40); do
    printf '@clean%s\n%s\n+\n%s\n' "$i" "$CLEAN" "$NOISY_QUAL"
  done
  printf '@noisy_a\n%s\n+\n%s\n' "$NOISY_MUTANT_A" "$NOISY_QUAL"
  printf '@noisy_b\n%s\n+\n%s\n' "$NOISY_MUTANT_B" "$NOISY_QUAL"
  printf '@noisy_double\n%s\n+\n%s\n' "$NOISY_DOUBLE_MUTANT" "$NOISY_QUAL"
} > "$NOISY_REP"

printf 'Running Java and Rust bbnorm-rs on noisy correction-staged multipass fixture...\n'
for case_name in final_only first_only both_passes; do
  case "$case_name" in
    final_only) ECC_STAGE=("ecc1=f" "eccf=t") ;;
    first_only) ECC_STAGE=("ecc1=t" "eccf=f") ;;
    both_passes) ECC_STAGE=("ecc=t") ;;
  esac

  java -Xmx1g -cp "$CP" jgi.KmerNormalize \
    "in=$NOISY_REP" "passes=2" "keepall=t" "ecco=f" "${ECC_STAGE[@]}" \
    "k=7" "minq=0" "minprob=0" "min=0" "minkmers=1" \
    "target=999999999" "max=999999999" "threads=1" "overwrite=t" "bits=32" \
    "out=$OUT/java.noisy_multipass.$case_name.fq" \
    >"$OUT/java.noisy_multipass.$case_name.stdout.log" 2>"$OUT/java.noisy_multipass.$case_name.stderr.log"

  target/debug/bbnorm-rs \
    "in=$NOISY_REP" "passes=2" "keepall=t" "ecco=f" "${ECC_STAGE[@]}" \
    "k=7" "minq=0" "minprob=0" "min=0" "minkmers=1" \
    "target=999999999" "max=999999999" "threads=1" "overwrite=t" "bits=32" \
    "out=$OUT/rust.noisy_multipass.$case_name.fq" \
    >"$OUT/rust.noisy_multipass.$case_name.stdout.log" 2>"$OUT/rust.noisy_multipass.$case_name.stderr.log"

  cmp "$OUT/java.noisy_multipass.$case_name.fq" "$OUT/rust.noisy_multipass.$case_name.fq"
  if grep -q "$NOISY_MUTANT_A" "$OUT/rust.noisy_multipass.$case_name.fq" ||
    grep -q "$NOISY_MUTANT_B" "$OUT/rust.noisy_multipass.$case_name.fq" ||
    grep -q "$NOISY_DOUBLE_MUTANT" "$OUT/rust.noisy_multipass.$case_name.fq"; then
    echo "Rust noisy multipass ECC left a mutant read uncorrected in $case_name" >&2
    exit 1
  fi
done

NORM_REP="$OUT/ecc_normalization_toss.fq"
{
  for i in $(seq 1 30); do
    printf '@clean%s\n%s\n+\n%s\n' "$i" "$CLEAN" "$NOISY_QUAL"
  done
  printf '@noisy_a\n%s\n+\n%s\n' "$NOISY_MUTANT_A" "$NOISY_QUAL"
  printf '@noisy_b\n%s\n+\n%s\n' "$NOISY_MUTANT_B" "$NOISY_QUAL"
  printf '@noisy_double\n%s\n+\n%s\n' "$NOISY_DOUBLE_MUTANT" "$NOISY_QUAL"
} > "$NORM_REP"

printf 'Running Java and Rust bbnorm-rs on ECC-enabled normalization/toss fixture...\n'
for case_name in tossbadreads_false tossbadreads_true; do
  case "$case_name" in
    tossbadreads_false) TOSS_OPTION="tossbadreads=f" ;;
    tossbadreads_true) TOSS_OPTION="tossbadreads=t" ;;
  esac

  java -Xmx1g -cp "$CP" jgi.KmerNormalize \
    "in=$NORM_REP" "passes=1" "ecc=t" "ecco=f" "$TOSS_OPTION" \
    "target=2" "max=2" \
    "k=7" "minq=0" "minprob=0" "min=0" "minkmers=1" \
    "threads=1" "overwrite=t" "bits=32" \
    "out=$OUT/java.norm_toss.$case_name.keep.fq" \
    "outt=$OUT/java.norm_toss.$case_name.toss.fq" \
    >"$OUT/java.norm_toss.$case_name.stdout.log" 2>"$OUT/java.norm_toss.$case_name.stderr.log"

  target/debug/bbnorm-rs \
    "in=$NORM_REP" "passes=1" "ecc=t" "ecco=f" "$TOSS_OPTION" \
    "target=2" "max=2" \
    "k=7" "minq=0" "minprob=0" "min=0" "minkmers=1" \
    "threads=1" "overwrite=t" "bits=32" \
    "out=$OUT/rust.norm_toss.$case_name.keep.fq" \
    "outt=$OUT/rust.norm_toss.$case_name.toss.fq" \
    >"$OUT/rust.norm_toss.$case_name.stdout.log" 2>"$OUT/rust.norm_toss.$case_name.stderr.log"

  cmp "$OUT/java.norm_toss.$case_name.keep.fq" "$OUT/rust.norm_toss.$case_name.keep.fq"
  cmp "$OUT/java.norm_toss.$case_name.toss.fq" "$OUT/rust.norm_toss.$case_name.toss.fq"
  test -s "$OUT/rust.norm_toss.$case_name.keep.fq"
  test -s "$OUT/rust.norm_toss.$case_name.toss.fq"

  if [[ "$case_name" == "tossbadreads_false" ]]; then
    grep -q '^@noisy_double$' "$OUT/rust.norm_toss.$case_name.keep.fq"
    grep -q "$CLEAN" "$OUT/rust.norm_toss.$case_name.keep.fq"
    if grep -q "$NOISY_DOUBLE_MUTANT" "$OUT/rust.norm_toss.$case_name.keep.fq"; then
      echo "Rust ECC normalization kept the double mutant without correction" >&2
      exit 1
    fi
    grep -q '^@noisy_a$' "$OUT/rust.norm_toss.$case_name.toss.fq"
    grep -q '^@noisy_b$' "$OUT/rust.norm_toss.$case_name.toss.fq"
  else
    grep -q '^@noisy_a$' "$OUT/rust.norm_toss.$case_name.toss.fq"
    grep -q '^@noisy_b$' "$OUT/rust.norm_toss.$case_name.toss.fq"
    grep -q '^@noisy_double$' "$OUT/rust.norm_toss.$case_name.toss.fq"
    if grep -q '^@noisy_' "$OUT/rust.norm_toss.$case_name.keep.fq"; then
      echo "Rust ECC normalization kept noisy reads despite tossbadreads=t" >&2
      exit 1
    fi
  fi
done

OVERLAP_R1="$OUT/ecc_overlap_r1.fq"
OVERLAP_R2="$OUT/ecc_overlap_r2.fq"
OVERLAP_READ1="TTAGTTGTGCCGCAGCGAAGTAGTGCTTGAAATATGCGAC"
OVERLAP_READ2_CLEAN="GTCGCATATTTCAAGCACTACTTCGCTGCGGCACAACTAA"
OVERLAP_READ2_MUTANT="GTCGCATATTTCAAGCACTAATTCGCTGCGGCACAACTAA"
OVERLAP_QUAL1="IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII"
OVERLAP_QUAL2="IIIIIIIIIIIIIIIIIIII!IIIIIIIIIIIIIIIIIII"
{
  for i in $(seq 1 5); do
    printf '@overlap%s/1\n%s\n+\n%s\n' "$i" "$OVERLAP_READ1" "$OVERLAP_QUAL1"
  done
} > "$OVERLAP_R1"
{
  for i in $(seq 1 5); do
    printf '@overlap%s/2\n%s\n+\n%s\n' "$i" "$OVERLAP_READ2_MUTANT" "$OVERLAP_QUAL2"
  done
} > "$OVERLAP_R2"

printf 'Running Rust bbnorm-rs on overlap-only ECC fixture...\n'
for ecco in f t; do
  target/debug/bbnorm-rs \
    "in=$OVERLAP_R1" "in2=$OVERLAP_R2" "passes=1" "keepall=t" "ecc=t" "ecco=$ecco" \
    "k=62" "minq=0" "minprob=0" \
    "threads=1" "overwrite=t" "bits=32" \
    "out=$OUT/rust.overlap.ecco_$ecco.1.fq" \
    "out2=$OUT/rust.overlap.ecco_$ecco.2.fq" \
    >"$OUT/rust.overlap.ecco_$ecco.stdout.log" 2>"$OUT/rust.overlap.ecco_$ecco.stderr.log"
done

grep -q "$OVERLAP_READ2_MUTANT" "$OUT/rust.overlap.ecco_f.2.fq"
if grep -q "$OVERLAP_READ2_CLEAN" "$OUT/rust.overlap.ecco_f.2.fq"; then
  echo "Rust overlap ECC changed reads despite ecco=f" >&2
  exit 1
fi
grep -q "$OVERLAP_READ2_CLEAN" "$OUT/rust.overlap.ecco_t.2.fq"
if grep -q "$OVERLAP_READ2_MUTANT" "$OUT/rust.overlap.ecco_t.2.fq"; then
  echo "Rust overlap ECC left the lower-quality mate base uncorrected" >&2
  exit 1
fi
grep -q 'paired overlap repair before the table-based ECC path' "$OUT/rust.overlap.ecco_t.stderr.log"

COMPETING_R1="$OUT/ecc_overlap_competing_r1.fq"
COMPETING_R2="$OUT/ecc_overlap_competing_r2.fq"
COMPETING_READ1="CAGTAACCAATGCCTGTTGAGATGCCAGACGCGTAACCAAAA"
COMPETING_READ2_MUTANT="TTTTGCTAACGCGTCTGGCATCTCAACAGGCATTGGTTAC"
COMPETING_QUAL1="IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII"
COMPETING_QUAL2="IIIII!I'IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII"
{
  for i in $(seq 1 5); do
    printf '@competing%s/1\n%s\n+\n%s\n' "$i" "$COMPETING_READ1" "$COMPETING_QUAL1"
  done
} > "$COMPETING_R1"
{
  for i in $(seq 1 5); do
    printf '@competing%s/2\n%s\n+\n%s\n' "$i" "$COMPETING_READ2_MUTANT" "$COMPETING_QUAL2"
  done
} > "$COMPETING_R2"

printf 'Running Java and Rust on competing-overlap ambiguity fixture...\n'
java -Xmx1g -cp vendor/BBTools-master/current jgi.KmerNormalize \
  "in=$COMPETING_R1" "in2=$COMPETING_R2" "passes=1" "keepall=t" "ecc=t" "ecco=t" \
  "k=62" "minq=0" "minprob=0" \
  "threads=1" "overwrite=t" "bits=32" \
  "out=$OUT/java.competing_overlap.1.fq" \
  "out2=$OUT/java.competing_overlap.2.fq" \
  >"$OUT/java.competing_overlap.stdout.log" 2>"$OUT/java.competing_overlap.stderr.log"

target/debug/bbnorm-rs \
  "in=$COMPETING_R1" "in2=$COMPETING_R2" "passes=1" "keepall=t" "ecc=t" "ecco=t" \
  "k=62" "minq=0" "minprob=0" \
  "threads=1" "overwrite=t" "bits=32" \
  "out=$OUT/rust.competing_overlap.1.fq" \
  "out2=$OUT/rust.competing_overlap.2.fq" \
  >"$OUT/rust.competing_overlap.stdout.log" 2>"$OUT/rust.competing_overlap.stderr.log"

cmp "$OUT/java.competing_overlap.1.fq" "$OUT/rust.competing_overlap.1.fq"
cmp "$OUT/java.competing_overlap.2.fq" "$OUT/rust.competing_overlap.2.fq"

printf 'ECC behavior passed. Rust no longer falls back, preserves real phiX no-error output, corrects one representative substitution, matches Java single-end and paired uncorrectable-read routing, preserves Java markuncorrectableerrors quality marking, matches Java multi-site markerrors ecclimit semantics, matches markerrors+trimaftermarking output, honors Java-shaped ecc1/eccf multipass staging, matches Java on a noisy correction-staged multipass fixture, matches Java ECC-enabled normalization/toss routing, and now also matches compact Java overlap behavior for accepted low-confidence repairs, rejected high-confidence mismatches, and competing short-overlap ambiguity fixtures. Logs and outputs: %s\n' "$OUT"

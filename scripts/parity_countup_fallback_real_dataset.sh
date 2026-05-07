#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT"

DATA1="vendor/BBTools-master/resources/sample1.fq.gz"
DATA2="vendor/BBTools-master/resources/sample2.fq.gz"
CP="vendor/BBTools-master/current"
OUT="${1:-tmp/real_countup_behavior}"
mkdir -p "$OUT"
rm -f "$OUT"/*
MINLEN_FIXTURE="$OUT/countup_minlen_short.fq"
ECC_FIXTURE="$OUT/countup_ecc_substitution.fq"
MARK_TRIM_FIXTURE="$OUT/countup_ecc_mark_trim.fq"
OVERLAP_R1="$OUT/countup_overlap_r1.fq"
OVERLAP_R2="$OUT/countup_overlap_r2.fq"
REJECT_R1="$OUT/countup_overlap_reject_r1.fq"
REJECT_R2="$OUT/countup_overlap_reject_r2.fq"

COMMON=(
  "in=$DATA1"
  "in2=$DATA2"
  "passes=1"
  "countup=t"
  "k=31"
  "minq=0"
  "minprob=0"
  "min=1"
  "minkmers=1"
  "target=1"
  "max=1"
  "threads=1"
  "overwrite=t"
  "bits=32"
)

cargo build --quiet

python - "$DATA1" "$MINLEN_FIXTURE" <<'PYFIXTURE'
import gzip
import sys

source, dest = sys.argv[1:3]
with gzip.open(source, "rt") as reader:
    record = [reader.readline(), reader.readline(), reader.readline(), reader.readline()]
if not all(record):
    raise SystemExit(f"{source} did not contain a complete FASTQ record")
with open(dest, "w") as writer:
    writer.writelines(record)
PYFIXTURE

CLEAN="ACGTTGCATGTCAGTACCGTAACGTTGCA"
MUTANT="ACGTTGCATGTCAGAACCGTAACGTTGCA"
TAIL_MUTANT="ACGTTGCATGTCAGTACCGTAACGTTAC"
TAIL_UNTRIMMED="ACGTTGCATGTCAGTACCGTAACGTTACA"
QUAL="IIIIIIIIIIIIIIIIIIIIIIIIIIIII"
OVERLAP_READ1="TTAGTTGTGCCGCAGCGAAGTAGTGCTTGAAATATGCGAC"
OVERLAP_READ2_CLEAN="GTCGCATATTTCAAGCACTACTTCGCTGCGGCACAACTAA"
OVERLAP_READ2_MUTANT="GTCGCATATTTCAAGCACTAATTCGCTGCGGCACAACTAA"
OVERLAP_QUAL1="IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII"
OVERLAP_QUAL2="IIIIIIIIIIIIIIIIIIII!IIIIIIIIIIIIIIIIIII"
REJECT_READ1="CAGTAACCAATGCCTGTTGAGATGCCAGACGCGTAACCAAAA"
REJECT_READ2_MUTANT="TTTTGCTAACGCGTCTGGCATCTCAACAGGCATTGGTTAC"
REJECT_QUAL1="IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII"
REJECT_QUAL2="IIIII!I'IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII"
{
  for i in $(seq 1 30); do
    printf '@clean%s\n%s\n+\n%s\n' "$i" "$CLEAN" "$QUAL"
  done
  printf '@mutant\n%s\n+\n%s\n' "$MUTANT" "$QUAL"
} > "$ECC_FIXTURE"
{
  for i in $(seq 1 30); do
    printf '@clean%s\n%s\n+\n%s\n' "$i" "$CLEAN" "$QUAL"
  done
  printf '@tail_mutant\n%s\n+\n%s\n' "$TAIL_UNTRIMMED" "$QUAL"
} > "$MARK_TRIM_FIXTURE"
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
{
  for i in $(seq 1 5); do
    printf '@reject%s/1\n%s\n+\n%s\n' "$i" "$REJECT_READ1" "$REJECT_QUAL1"
  done
} > "$REJECT_R1"
{
  for i in $(seq 1 5); do
    printf '@reject%s/2\n%s\n+\n%s\n' "$i" "$REJECT_READ2_MUTANT" "$REJECT_QUAL2"
  done
} > "$REJECT_R2"

printf 'Confirming vendored Java BBNorm countup=t is not a normal comparable output path...\n'
if java -Xmx1g -cp "$CP" jgi.KmerNormalize \
  "${COMMON[@]}" \
  "out=$OUT/java.keep1.fq" "out2=$OUT/java.keep2.fq" \
  "outt=$OUT/java.toss1.fq" "outt2=$OUT/java.toss2.fq" \
  "hist=$OUT/java.hist.tsv" "rhist=$OUT/java.rhist.tsv" \
  >"$OUT/java.stdout.log" 2>"$OUT/java.stderr.log"; then
  echo 'Expected vendored Java countup=t to fail on this guard fixture, but it succeeded.' >&2
  exit 1
fi

grep -Eq 'NullPointerException|BBNorm terminated in an error state' "$OUT/java.stderr.log"

printf 'Confirming Rust countup=t runs exact count-up normalization and produces usable outputs...\n'
target/debug/bbnorm-rs \
  "${COMMON[@]}" \
  "out=$OUT/rust.keep1.fq" "out2=$OUT/rust.keep2.fq" \
  "outt=$OUT/rust.toss1.fq" "outt2=$OUT/rust.toss2.fq" \
  "hist=$OUT/rust.hist.tsv" "rhist=$OUT/rust.rhist.tsv" \
  >"$OUT/rust.stdout.log" 2>"$OUT/rust.stderr.log"

test -s "$OUT/rust.hist.tsv"
test -s "$OUT/rust.rhist.tsv"
test -e "$OUT/rust.keep1.fq"
test -e "$OUT/rust.keep2.fq"
test -e "$OUT/rust.toss1.fq"
test -e "$OUT/rust.toss2.fq"

printf 'Confirming Rust countup=t can use a bounded count-min sketch table...\n'
target/debug/bbnorm-rs \
  "${COMMON[@]}" \
  "hashes=2" \
  "bits=8" \
  "sketchmemory=1k" \
  "out=$OUT/rust.sketch.keep1.fq" "out2=$OUT/rust.sketch.keep2.fq" \
  "outt=$OUT/rust.sketch.toss1.fq" "outt2=$OUT/rust.sketch.toss2.fq" \
  "hist=$OUT/rust.sketch.hist.tsv" "rhist=$OUT/rust.sketch.rhist.tsv" \
  >"$OUT/rust.sketch.stdout.log" 2>"$OUT/rust.sketch.stderr.log"

grep -q 'count-min memory budget' "$OUT/rust.sketch.stderr.log"
test -s "$OUT/rust.sketch.hist.tsv"
test -s "$OUT/rust.sketch.rhist.tsv"
test -e "$OUT/rust.sketch.keep1.fq"
test -e "$OUT/rust.sketch.keep2.fq"
test -e "$OUT/rust.sketch.toss1.fq"
test -e "$OUT/rust.sketch.toss2.fq"

printf 'Confirming Rust countup=t addbadreadscountup=t is accepted and produces usable outputs...\n'
target/debug/bbnorm-rs \
  "${COMMON[@]}" \
  "addbadreadscountup=t" \
  "rename=t" \
  "out=$OUT/rust.abrc.keep1.fq" "out2=$OUT/rust.abrc.keep2.fq" \
  "outt=$OUT/rust.abrc.toss1.fq" "outt2=$OUT/rust.abrc.toss2.fq" \
  "outlow=$OUT/rust.abrc.low1.fq" "outlow2=$OUT/rust.abrc.low2.fq" \
  "outmid=$OUT/rust.abrc.mid1.fq" "outmid2=$OUT/rust.abrc.mid2.fq" \
  "outhigh=$OUT/rust.abrc.high1.fq" "outhigh2=$OUT/rust.abrc.high2.fq" \
  "outuncorrected=$OUT/rust.abrc.unc1.fq" "outuncorrected2=$OUT/rust.abrc.unc2.fq" \
  "hist=$OUT/rust.abrc.hist.tsv" "rhist=$OUT/rust.abrc.rhist.tsv" \
  >"$OUT/rust.abrc.stdout.log" 2>"$OUT/rust.abrc.stderr.log"

test -s "$OUT/rust.abrc.hist.tsv"
test -s "$OUT/rust.abrc.rhist.tsv"
test -e "$OUT/rust.abrc.keep1.fq"
test -e "$OUT/rust.abrc.keep2.fq"
test -e "$OUT/rust.abrc.toss1.fq"
test -e "$OUT/rust.abrc.toss2.fq"
test -e "$OUT/rust.abrc.low1.fq"
test -e "$OUT/rust.abrc.low2.fq"
test -e "$OUT/rust.abrc.mid1.fq"
test -e "$OUT/rust.abrc.mid2.fq"
test -e "$OUT/rust.abrc.high1.fq"
test -e "$OUT/rust.abrc.high2.fq"
test -e "$OUT/rust.abrc.unc1.fq"
test -e "$OUT/rust.abrc.unc2.fq"
if ! grep -hqm1 '^@id=.*d1=' "$OUT"/rust.abrc.keep*.fq "$OUT"/rust.abrc.toss*.fq "$OUT"/rust.abrc.low*.fq "$OUT"/rust.abrc.mid*.fq "$OUT"/rust.abrc.high*.fq; then
  echo 'Expected renamed count-up side outputs to contain BBTools-style id/d1 headers.' >&2
  exit 1
fi

printf 'Confirming Rust countup=t applies Java-shaped prepass minlen filtering...\n'
target/debug/bbnorm-rs \
  "in=$MINLEN_FIXTURE" \
  "passes=1" \
  "countup=t" \
  "k=31" \
  "minq=0" \
  "minprob=0" \
  "min=1" \
  "minkmers=1" \
  "target=1" \
  "max=1" \
  "minlen=101" \
  "threads=1" \
  "overwrite=t" \
  "bits=32" \
  "out=$OUT/rust.minlen.keep.fq" \
  "outt=$OUT/rust.minlen.toss.fq" \
  >"$OUT/rust.minlen.stdout.log" 2>"$OUT/rust.minlen.stderr.log"

test ! -s "$OUT/rust.minlen.keep.fq"
test ! -s "$OUT/rust.minlen.toss.fq"

printf 'Confirming Rust countup=t addbadreadscountup=t carries prepass-tossed reads forward...\n'
target/debug/bbnorm-rs \
  "in=$MINLEN_FIXTURE" \
  "passes=1" \
  "countup=t" \
  "addbadreadscountup=t" \
  "k=31" \
  "minq=0" \
  "minprob=0" \
  "min=1" \
  "minkmers=1" \
  "target=1" \
  "max=1" \
  "minlen=101" \
  "threads=1" \
  "overwrite=t" \
  "bits=32" \
  "out=$OUT/rust.minlen.abrc.keep.fq" \
  "outt=$OUT/rust.minlen.abrc.toss.fq" \
  >"$OUT/rust.minlen.abrc.stdout.log" 2>"$OUT/rust.minlen.abrc.stderr.log"

test ! -s "$OUT/rust.minlen.abrc.keep.fq"
test -s "$OUT/rust.minlen.abrc.toss.fq"

printf 'Confirming Rust countup=t ecc=t performs table-based correction...\n'
target/debug/bbnorm-rs \
  "in=$ECC_FIXTURE" \
  "passes=1" \
  "countup=t" \
  "keepall=t" \
  "ecc=t" \
  "k=7" \
  "minq=0" \
  "minprob=0" \
  "min=0" \
  "minkmers=1" \
  "target=999999999" \
  "max=999999999" \
  "threads=1" \
  "overwrite=t" \
  "bits=32" \
  "out=$OUT/rust.ecc.keep.fq" \
  "outuncorrected=$OUT/rust.ecc.unc.fq" \
  >"$OUT/rust.ecc.stdout.log" 2>"$OUT/rust.ecc.stderr.log"

if grep -q "$MUTANT" "$OUT/rust.ecc.keep.fq"; then
  echo "Rust countup ECC left the mutant read uncorrected" >&2
  exit 1
fi
grep -q "$CLEAN" "$OUT/rust.ecc.keep.fq"
tail -4 "$OUT/rust.ecc.keep.fq" | grep -q "$CLEAN"
test ! -s "$OUT/rust.ecc.unc.fq"

printf 'Confirming Rust countup=t ecc=t routes uncorrectable reads to outuncorrected...\n'
target/debug/bbnorm-rs \
  "in=$ECC_FIXTURE" \
  "passes=1" \
  "countup=t" \
  "keepall=t" \
  "ecc=t" \
  "ecco=f" \
  "eccmaxqual=0" \
  "k=7" \
  "minq=0" \
  "minprob=0" \
  "min=0" \
  "minkmers=1" \
  "target=999999999" \
  "max=999999999" \
  "threads=1" \
  "overwrite=t" \
  "bits=32" \
  "out=$OUT/rust.ecc_unc.keep.fq" \
  "outuncorrected=$OUT/rust.ecc_unc.unc.fq" \
  >"$OUT/rust.ecc_unc.stdout.log" 2>"$OUT/rust.ecc_unc.stderr.log"

grep -q "$MUTANT" "$OUT/rust.ecc_unc.keep.fq"
grep -q "$MUTANT" "$OUT/rust.ecc_unc.unc.fq"
tail -4 "$OUT/rust.ecc_unc.unc.fq" | grep -q "$MUTANT"

printf 'Confirming Rust countup=t markuncorrectableerrors=t preserves marked qualities on uncorrectable reads...\n'
target/debug/bbnorm-rs \
  "in=$ECC_FIXTURE" \
  "passes=1" \
  "countup=t" \
  "keepall=t" \
  "ecc=t" \
  "ecco=f" \
  "eccmaxqual=0" \
  "markuncorrectableerrors=t" \
  "k=7" \
  "minq=0" \
  "minprob=0" \
  "min=0" \
  "minkmers=1" \
  "target=999999999" \
  "max=999999999" \
  "threads=1" \
  "overwrite=t" \
  "bits=32" \
  "out=$OUT/rust.ecc_unc_mark.keep.fq" \
  "outuncorrected=$OUT/rust.ecc_unc_mark.unc.fq" \
  >"$OUT/rust.ecc_unc_mark.stdout.log" 2>"$OUT/rust.ecc_unc_mark.stderr.log"

tail -4 "$OUT/rust.ecc_unc_mark.keep.fq" > "$OUT/rust.ecc_unc_mark.keep.tail.fq"
cmp "$OUT/rust.ecc_unc_mark.keep.tail.fq" "$OUT/rust.ecc_unc_mark.unc.fq"
if tail -1 "$OUT/rust.ecc_unc_mark.keep.fq" | grep -q '^IIIIIIIIIIIIIIIIIIIIIIIIIIIII$'; then
  echo "Rust countup markuncorrectableerrors left the original qualities unchanged" >&2
  exit 1
fi

printf 'Confirming Rust countup=t markerrors+trimaftermarking works with ECC...\n'
target/debug/bbnorm-rs \
  "in=$MARK_TRIM_FIXTURE" \
  "passes=1" \
  "countup=t" \
  "keepall=t" \
  "ecc=t" \
  "ecco=f" \
  "markerrors=t" \
  "trimaftermarking=t" \
  "qtrim=r" \
  "trimq=20" \
  "optitrim=f" \
  "k=7" \
  "minq=0" \
  "minprob=0" \
  "min=0" \
  "minkmers=1" \
  "target=999999999" \
  "max=999999999" \
  "threads=1" \
  "overwrite=t" \
  "bits=32" \
  "out=$OUT/rust.marktrim.fq" \
  >"$OUT/rust.marktrim.stdout.log" 2>"$OUT/rust.marktrim.stderr.log"

tail -4 "$OUT/rust.marktrim.fq" | grep -q "$TAIL_MUTANT"
if tail -4 "$OUT/rust.marktrim.fq" | grep -q "$TAIL_UNTRIMMED"; then
  echo "Rust countup markerrors+trimaftermarking left the tail mutant untrimmed" >&2
  exit 1
fi

printf 'Confirming Rust countup=t ecco=t performs overlap-only mate repair...\n'
for ecco in f t; do
  target/debug/bbnorm-rs \
    "in=$OVERLAP_R1" \
    "in2=$OVERLAP_R2" \
    "passes=1" \
    "countup=t" \
    "keepall=t" \
    "ecc=t" \
    "ecco=$ecco" \
    "k=62" \
    "minq=0" \
    "minprob=0" \
    "min=0" \
    "minkmers=1" \
    "target=999999999" \
    "max=999999999" \
    "threads=1" \
    "overwrite=t" \
    "bits=32" \
    "out=$OUT/rust.overlap.ecco_$ecco.1.fq" \
    "out2=$OUT/rust.overlap.ecco_$ecco.2.fq" \
    >"$OUT/rust.overlap.ecco_$ecco.stdout.log" 2>"$OUT/rust.overlap.ecco_$ecco.stderr.log"
done

grep -q "$OVERLAP_READ2_MUTANT" "$OUT/rust.overlap.ecco_f.2.fq"
if grep -q "$OVERLAP_READ2_CLEAN" "$OUT/rust.overlap.ecco_f.2.fq"; then
  echo "Rust countup overlap ECC changed reads despite ecco=f" >&2
  exit 1
fi
grep -q "$OVERLAP_READ2_CLEAN" "$OUT/rust.overlap.ecco_t.2.fq"
if grep -q "$OVERLAP_READ2_MUTANT" "$OUT/rust.overlap.ecco_t.2.fq"; then
  echo "Rust countup overlap ECC left the lower-quality mate base uncorrected" >&2
  exit 1
fi
grep -q 'paired overlap repair before the table-based ECC path' "$OUT/rust.overlap.ecco_t.stderr.log"

printf 'Confirming Rust countup=t ecco=t also rejects competing short-overlap ambiguity...\n'
for ecco in f t; do
  target/debug/bbnorm-rs \
    "in=$REJECT_R1" \
    "in2=$REJECT_R2" \
    "passes=1" \
    "countup=t" \
    "keepall=t" \
    "ecc=t" \
    "ecco=$ecco" \
    "k=62" \
    "minq=0" \
    "minprob=0" \
    "min=0" \
    "minkmers=1" \
    "target=999999999" \
    "max=999999999" \
    "threads=1" \
    "overwrite=t" \
    "bits=32" \
    "out=$OUT/rust.reject.ecco_$ecco.1.fq" \
    "out2=$OUT/rust.reject.ecco_$ecco.2.fq" \
    >"$OUT/rust.reject.ecco_$ecco.stdout.log" 2>"$OUT/rust.reject.ecco_$ecco.stderr.log"
done

grep -q "$REJECT_READ2_MUTANT" "$OUT/rust.reject.ecco_f.2.fq"
grep -q "$REJECT_READ2_MUTANT" "$OUT/rust.reject.ecco_t.2.fq"
if grep -q 'TTTTGGTTACGCGTCTGGCATCTCAACAGGCATTGGTTAC' "$OUT/rust.reject.ecco_t.2.fq"; then
  echo "Rust countup overlap ECC over-corrected a competing short-overlap ambiguity fixture" >&2
  exit 1
fi

printf 'Countup Rust behavior passed. Vendored Java fails on this local guard fixture; Rust runs count-up keep/toss output, including bounded sketch kept-count tables via sketchmemory, addbadreadscountup=t table updates, renamed headers, side-output stream creation, minlen filtering, table-based ECC, count-up outuncorrected routing, markerrors+trimaftermarking, overlap-only mate repair for accepted ecco=t overlaps, and strict rejection of competing short-overlap ambiguity fixtures. Logs: %s\n' "$OUT"

#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT"

DATA1="vendor/BBTools-master/resources/sample1.fq.gz"
DATA2="vendor/BBTools-master/resources/sample2.fq.gz"
CP="vendor/BBTools-master/current"
OUT="${1:-tmp/real_side_output_stats_fallback}"
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

STATS_FLAGS=(
  "qhist=$OUT/aux.qhist.tsv"
  "bqhist=$OUT/aux.bqhist.tsv"
  "qchist=$OUT/aux.qchist.tsv"
  "aqhist=$OUT/aux.aqhist.tsv"
  "obqhist=$OUT/aux.obqhist.tsv"
  "mhist=$OUT/aux.mhist.tsv"
  "ihist=$OUT/aux.ihist.tsv"
  "qahist=$OUT/aux.qahist.tsv"
  "indelhist=$OUT/aux.indelhist.tsv"
  "ehist=$OUT/aux.ehist.tsv"
  "bhist=$OUT/aux.bhist.tsv"
  "lhist=$OUT/aux.lhist.tsv"
  "gchist=$OUT/aux.gchist.tsv"
  "enthist=$OUT/aux.enthist.tsv"
  "idhist=$OUT/aux.idhist.tsv"
  "gcbins=auto"
  "entropybins=100"
  "idbins=auto"
  "gcplot=f"
  "entropyns=t"
  "maxhistlen=1k"
)

cargo build --quiet

printf 'Running Java BBNorm local baseline on paired phiX...\n'
java -Xmx1g -cp "$CP" jgi.KmerNormalize \
  "${COMMON[@]}" \
  "out=$OUT/java.keep1.fq" "out2=$OUT/java.keep2.fq" \
  "outlow=$OUT/java.low1.fq" "outlow2=$OUT/java.low2.fq" \
  "outmid=$OUT/java.mid1.fq" "outmid2=$OUT/java.mid2.fq" \
  "outhigh=$OUT/java.high1.fq" "outhigh2=$OUT/java.high2.fq" \
  "hist=$OUT/java.hist.tsv" "rhist=$OUT/java.rhist.tsv" \
  >"$OUT/java.stdout.log" 2>"$OUT/java.stderr.log"

printf 'Confirming vendored KmerNormalize rejects shared side-output stats controls...\n'
if java -Xmx1g -cp "$CP" jgi.KmerNormalize \
  "${COMMON[@]}" "${STATS_FLAGS[@]}" \
  "out=$OUT/java.reject.keep1.fq" "out2=$OUT/java.reject.keep2.fq" \
  >"$OUT/java.reject.stdout.log" 2>"$OUT/java.reject.stderr.log"; then
  printf 'Expected vendored Java to reject side-output stats controls, but it succeeded.\n' >&2
  exit 1
fi
grep -q 'Unknown parameter qhist=' "$OUT/java.reject.stderr.log"

printf 'Running Rust bbnorm-rs with shared side-output stats controls on paired phiX...\n'
target/debug/bbnorm-rs \
  "${COMMON[@]}" "${STATS_FLAGS[@]}" \
  "out=$OUT/rust.keep1.fq" "out2=$OUT/rust.keep2.fq" \
  "outlow=$OUT/rust.low1.fq" "outlow2=$OUT/rust.low2.fq" \
  "outmid=$OUT/rust.mid1.fq" "outmid2=$OUT/rust.mid2.fq" \
  "outhigh=$OUT/rust.high1.fq" "outhigh2=$OUT/rust.high2.fq" \
  "hist=$OUT/rust.hist.tsv" "rhist=$OUT/rust.rhist.tsv" \
  >"$OUT/rust.stdout.log" 2>"$OUT/rust.stderr.log"

grep -q 'qhist=.*side-output quality histogram' "$OUT/rust.stderr.log"
grep -q 'bqhist=.*side-output base-quality histogram' "$OUT/rust.stderr.log"
grep -q 'qchist=.*side-output quality-count histogram' "$OUT/rust.stderr.log"
grep -q 'aqhist=.*side-output average-quality histogram' "$OUT/rust.stderr.log"
grep -q 'obqhist=.*side-output overall base-quality histogram' "$OUT/rust.stderr.log"
grep -q 'mhist=.*side-output match histogram' "$OUT/rust.stderr.log"
grep -q 'ihist=.*side-output insert histogram' "$OUT/rust.stderr.log"
grep -q 'qahist=.*side-output quality-accuracy histogram' "$OUT/rust.stderr.log"
grep -q 'indelhist=.*side-output indel histogram' "$OUT/rust.stderr.log"
grep -q 'ehist=.*side-output error histogram' "$OUT/rust.stderr.log"
grep -q 'lhist=.*side-output length histogram' "$OUT/rust.stderr.log"
grep -q 'gchist=.*side-output GC histogram' "$OUT/rust.stderr.log"
grep -q 'bhist=.*side-output base-content histogram' "$OUT/rust.stderr.log"
grep -q 'enthist=.*side-output entropy histogram' "$OUT/rust.stderr.log"
grep -q 'idhist=.*side-output identity histogram' "$OUT/rust.stderr.log"
grep -q 'gcbins=auto is a BBTools side-output GC histogram sizing control' "$OUT/rust.stderr.log"
grep -q 'entropybins=100 is a BBTools side-output entropy histogram sizing control' "$OUT/rust.stderr.log"
grep -q 'idbins=auto is a BBTools side-output identity histogram sizing control' "$OUT/rust.stderr.log"
grep -q 'entropyns=t is a BBTools side-output entropy control' "$OUT/rust.stderr.log"
grep -q 'maxhistlen=1k is a BBTools side-output histogram length control' "$OUT/rust.stderr.log"
test -s "$OUT/aux.qhist.tsv"
grep -q '^#Quality[[:space:]]Bases$' "$OUT/aux.qhist.tsv"
test -s "$OUT/aux.bqhist.tsv"
grep -q '^#BaseNum[[:space:]]count_1[[:space:]]min_1[[:space:]]max_1[[:space:]]mean_1[[:space:]]Q1_1[[:space:]]med_1[[:space:]]Q3_1[[:space:]]LW_1[[:space:]]RW_1[[:space:]]count_2[[:space:]]min_2[[:space:]]max_2[[:space:]]mean_2[[:space:]]Q1_2[[:space:]]med_2[[:space:]]Q3_2[[:space:]]LW_2[[:space:]]RW_2$' "$OUT/aux.bqhist.tsv"
test -s "$OUT/aux.qchist.tsv"
grep -q '^#Quality[[:space:]]count1[[:space:]]fraction1[[:space:]]count2[[:space:]]fraction2$' "$OUT/aux.qchist.tsv"
test -s "$OUT/aux.aqhist.tsv"
grep -q '^#Quality[[:space:]]count1[[:space:]]fraction1[[:space:]]count2[[:space:]]fraction2$' "$OUT/aux.aqhist.tsv"
test -s "$OUT/aux.obqhist.tsv"
grep -q '^#Quality[[:space:]]bases[[:space:]]fraction$' "$OUT/aux.obqhist.tsv"
test -s "$OUT/aux.lhist.tsv"
grep -q '^#Length[[:space:]]Reads[[:space:]]Bases$' "$OUT/aux.lhist.tsv"
test -s "$OUT/aux.gchist.tsv"
grep -q '^#GC_Bin[[:space:]]Reads[[:space:]]Bases$' "$OUT/aux.gchist.tsv"
test -s "$OUT/aux.bhist.tsv"
grep -q '^#Pos[[:space:]]A[[:space:]]C[[:space:]]G[[:space:]]T[[:space:]]N$' "$OUT/aux.bhist.tsv"
test -s "$OUT/aux.enthist.tsv"
grep -q '^#Value[[:space:]]Count$' "$OUT/aux.enthist.tsv"
test -s "$OUT/aux.idhist.tsv"
grep -q '^#Identity[[:space:]]Reads[[:space:]]Bases$' "$OUT/aux.idhist.tsv"
test -s "$OUT/aux.mhist.tsv"
grep -q '^#BaseNum[[:space:]]Match1[[:space:]]Sub1[[:space:]]Del1[[:space:]]Ins1[[:space:]]N1[[:space:]]Other1[[:space:]]Match2[[:space:]]Sub2[[:space:]]Del2[[:space:]]Ins2[[:space:]]N2[[:space:]]Other2$' "$OUT/aux.mhist.tsv"
test -s "$OUT/aux.ihist.tsv"
grep -q '^#InsertSize[[:space:]]Count$' "$OUT/aux.ihist.tsv"
test -s "$OUT/aux.qahist.tsv"
grep -q '^#Quality[[:space:]]Match[[:space:]]Sub[[:space:]]Ins[[:space:]]Del[[:space:]]TrueQuality[[:space:]]TrueQualitySub$' "$OUT/aux.qahist.tsv"
test -s "$OUT/aux.indelhist.tsv"
grep -q '^#Length[[:space:]]Deletions[[:space:]]Insertions$' "$OUT/aux.indelhist.tsv"
test -s "$OUT/aux.ehist.tsv"
grep -q '^#Errors[[:space:]]Count$' "$OUT/aux.ehist.tsv"

python3 - "$OUT/aux.qhist.tsv" "$OUT/aux.lhist.tsv" "$OUT/aux.gchist.tsv" "$OUT/aux.bhist.tsv" "$OUT/aux.qchist.tsv" "$OUT/aux.aqhist.tsv" "$OUT/aux.obqhist.tsv" "$OUT/aux.bqhist.tsv" "$OUT/aux.enthist.tsv" "$OUT/aux.idhist.tsv" "$OUT/aux.mhist.tsv" "$OUT/aux.qahist.tsv" "$OUT/aux.ehist.tsv" <<'PY'
import sys
from pathlib import Path

qhist = Path(sys.argv[1])
quality_bases = 0
for line in qhist.read_text().splitlines():
    if not line or line.startswith('#'):
        continue
    quality, base_count = map(int, line.split('\t'))
    if quality < 0:
        raise SystemExit(f'quality histogram row has a negative quality: {line}')
    quality_bases += base_count
if quality_bases != 20000:
    raise SystemExit(f'expected 20000 bases in qhist, saw {quality_bases}')

path = Path(sys.argv[2])
reads = 0
bases = 0
for line in path.read_text().splitlines():
    if not line or line.startswith('#'):
        continue
    length, count, base_count = map(int, line.split('\t'))
    if count:
        if base_count != length * count:
            raise SystemExit(f'length histogram row has inconsistent bases: {line}')
        reads += count
        bases += base_count
if reads != 200 or bases != 20000:
    raise SystemExit(f'expected 200 reads and 20000 bases in lhist, saw {reads} reads and {bases} bases')

path = Path(sys.argv[3])
reads = 0
bases = 0
for line in path.read_text().splitlines():
    if not line or line.startswith('#'):
        continue
    gc_bin, count, base_count = map(int, line.split('\t'))
    if gc_bin < 0:
        raise SystemExit(f'GC histogram row has a negative bin: {line}')
    reads += count
    bases += base_count
if reads != 200 or bases != 20000:
    raise SystemExit(f'expected 200 reads and 20000 bases in gchist, saw {reads} reads and {bases} bases')

path = Path(sys.argv[4])
rows = 0
last_pos = -1
for line in path.read_text().splitlines():
    if not line or line.startswith('#'):
        continue
    fields = line.split('\t')
    if len(fields) != 6:
        raise SystemExit(f'expected six base histogram fields, saw: {line}')
    pos = int(fields[0])
    if pos != last_pos + 1:
        raise SystemExit(f'base histogram positions are not contiguous at: {line}')
    values = [float(value) for value in fields[1:]]
    if not (0.9999 <= sum(values) <= 1.0001):
        raise SystemExit(f'base histogram row does not sum to 1.0: {line}')
    rows += 1
    last_pos = pos
if rows != 200:
    raise SystemExit(f'expected 200 base-content rows for paired 100bp phiX reads, saw {rows}')

path = Path(sys.argv[5])
read1_bases = 0
read2_bases = 0
for line in path.read_text().splitlines():
    if not line or line.startswith('#'):
        continue
    quality, count1, fraction1, count2, fraction2 = line.split('\t')
    read1_bases += int(count1)
    read2_bases += int(count2)
if read1_bases != 10000 or read2_bases != 10000:
    raise SystemExit(f'expected qchist to split 10000/10000 bases, saw {read1_bases}/{read2_bases}')

path = Path(sys.argv[6])
read1_reads = 0
read2_reads = 0
for line in path.read_text().splitlines():
    if not line or line.startswith('#'):
        continue
    quality, count1, fraction1, count2, fraction2 = line.split('\t')
    read1_reads += int(count1)
    read2_reads += int(count2)
if read1_reads != 100 or read2_reads != 100:
    raise SystemExit(f'expected aqhist to split 100/100 reads, saw {read1_reads}/{read2_reads}')

path = Path(sys.argv[7])
overall_bases = 0
for line in path.read_text().splitlines():
    if not line or line.startswith('#'):
        continue
    quality, bases, fraction = line.split('\t')
    overall_bases += int(bases)
if overall_bases != 20000:
    raise SystemExit(f'expected obqhist to count 20000 bases, saw {overall_bases}')

path = Path(sys.argv[8])
bq_rows = 0
for line in path.read_text().splitlines():
    if not line or line.startswith('#'):
        continue
    fields = line.split('\t')
    if len(fields) != 19:
        raise SystemExit(f'expected paired bqhist row to have 19 fields, saw: {line}')
    pos = int(fields[0])
    count1 = int(fields[1])
    count2 = int(fields[10])
    if count1 != 100 or count2 != 100:
        raise SystemExit(f'expected 100 read1/read2 bases per bqhist position, saw: {line}')
    bq_rows += 1
if bq_rows != 100:
    raise SystemExit(f'expected 100 bqhist rows for 100bp phiX reads, saw {bq_rows}')

path = Path(sys.argv[9])
entropy_reads = 0
for line in path.read_text().splitlines():
    if not line or line.startswith('#'):
        continue
    value, count = line.split('\t')
    value = float(value)
    if not (0.0 <= value <= 1.0):
        raise SystemExit(f'expected entropy value in [0, 1], saw: {line}')
    entropy_reads += int(count)
if entropy_reads != 200:
    raise SystemExit(f'expected enthist to count 200 reads, saw {entropy_reads}')

path = Path(sys.argv[10])
identity_reads = 0
identity_bases = 0
seen_full_identity = False
for line in path.read_text().splitlines():
    if not line or line.startswith('#'):
        continue
    identity, reads, bases = line.split('\t')
    identity = float(identity)
    reads = int(reads)
    bases = int(bases)
    if not (0.0 <= identity <= 100.0):
        raise SystemExit(f'expected identity value in [0, 100], saw: {line}')
    identity_reads += reads
    identity_bases += bases
    seen_full_identity |= identity == 100.0 and reads == 200 and bases == 20000
if identity_reads != 200 or identity_bases != 20000 or not seen_full_identity:
    raise SystemExit(f'expected idhist to put 200 reads/20000 bases at 100 identity, saw {identity_reads}/{identity_bases}')

path = Path(sys.argv[11])
mhist_rows = 0
for line in path.read_text().splitlines():
    if not line or line.startswith('#'):
        continue
    fields = line.split('\t')
    if len(fields) != 13:
        raise SystemExit(f'expected paired mhist row to have 13 fields, saw: {line}')
    first = [float(value) for value in fields[1:7]]
    second = [float(value) for value in fields[7:13]]
    if not (0.9999 <= sum(first) <= 1.0001 and 0.9999 <= sum(second) <= 1.0001):
        raise SystemExit(f'mhist fallback fractions should sum to 1.0 per mate: {line}')
    mhist_rows += 1
if mhist_rows != 100:
    raise SystemExit(f'expected 100 mhist rows for 100bp phiX reads, saw {mhist_rows}')

path = Path(sys.argv[12])
qahist_matches = 0
for line in path.read_text().splitlines():
    if not line or line.startswith('#'):
        continue
    fields = line.split('\t')
    if len(fields) < 5:
        raise SystemExit(f'expected qahist row to include quality and four counts, saw: {line}')
    qahist_matches += int(fields[1])
if qahist_matches != 20000:
    raise SystemExit(f'expected qahist fallback to mark 20000 matching quality observations, saw {qahist_matches}')

path = Path(sys.argv[13])
error_zero_reads = 0
for line in path.read_text().splitlines():
    if not line or line.startswith('#'):
        continue
    errors, count = map(int, line.split('\t'))
    if errors == 0:
        error_zero_reads += count
if error_zero_reads != 200:
    raise SystemExit(f'expected ehist fallback to put 200 reads at zero observed errors, saw {error_zero_reads}')
PY

for suffix in keep1.fq keep2.fq low1.fq low2.fq mid1.fq mid2.fq high1.fq high2.fq hist.tsv rhist.tsv; do
  cmp "$OUT/java.$suffix" "$OUT/rust.$suffix"
done

printf 'Side-output stats fallback passed. Vendored KmerNormalize rejects these controls; Rust accepts them, emits covered quality-family, read-length, GC, base-content, entropy, sequence-input identity, and no-alignment match/error histograms, and matches local Java FASTQ/hist baseline. Logs: %s\n' "$OUT"

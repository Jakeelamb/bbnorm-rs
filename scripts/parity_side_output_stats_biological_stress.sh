#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT"

BIO_ROOT="${BIO_ROOT:-/home/jake/Projects/biological data}"
DATA1="${DATA1:-$BIO_ROOT/reads/short_reads/scer_s288c_pe_srr23631023/SRR23631023_1.fastq.gz}"
DATA2="${DATA2:-$BIO_ROOT/reads/short_reads/scer_s288c_pe_srr23631023/SRR23631023_2.fastq.gz}"
OUT="${1:-tmp/side_output_stats_biological_stress}"
READS="${READS:-1000}"
TABLE_READS="${TABLE_READS:-$READS}"
THREAD_CASES="${THREAD_CASES:-${THREADS:-1 2 auto}}"
PYTHON="${PYTHON:-python}"
mkdir -p "$OUT"
rm -f "$OUT"/*

if [[ ! -f "$DATA1" || ! -f "$DATA2" ]]; then
  printf 'Missing paired dataset files:\n  DATA1=%s\n  DATA2=%s\n' "$DATA1" "$DATA2" >&2
  exit 2
fi
for value_name in READS TABLE_READS; do
  value="${!value_name}"
  if [[ ! "$value" =~ ^[0-9]+$ || "$value" -lt 1 ]]; then
    printf '%s must be a positive integer; got %s\n' "$value_name" "$value" >&2
    exit 2
  fi
done

label_for_threads() {
  local label="threads_$1"
  printf '%s\n' "${label//[^A-Za-z0-9_]/_}"
}

cargo build --quiet

COMMON=(
  "in=$DATA1"
  "in2=$DATA2"
  "passes=1"
  "reads=$READS"
  "tablereads=$TABLE_READS"
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
  "overwrite=t"
  "bits=32"
)

run_case() {
  local threads="$1"
  local label="$2"
  local prefix="$OUT/$label"

  printf 'Running paired biological side-output stats reads=%s threads=%s...\n' "$READS" "$threads"
  target/debug/bbnorm-rs \
    "${COMMON[@]}" \
    "threads=$threads" \
    "out=$prefix.keep1.fq" "out2=$prefix.keep2.fq" \
    "outlow=$prefix.low1.fq" "outlow2=$prefix.low2.fq" \
    "outmid=$prefix.mid1.fq" "outmid2=$prefix.mid2.fq" \
    "outhigh=$prefix.high1.fq" "outhigh2=$prefix.high2.fq" \
    "hist=$prefix.hist.tsv" "rhist=$prefix.rhist.tsv" \
    "qhist=$prefix.qhist.tsv" \
    "bqhist=$prefix.bqhist.tsv" \
    "qchist=$prefix.qchist.tsv" \
    "aqhist=$prefix.aqhist.tsv" \
    "obqhist=$prefix.obqhist.tsv" \
    "lhist=$prefix.lhist.tsv" \
    "gchist=$prefix.gchist.tsv" \
    "bhist=$prefix.bhist.tsv" \
    "enthist=$prefix.enthist.tsv" \
    "idhist=$prefix.idhist.tsv" \
    "mhist=$prefix.mhist.tsv" \
    "ihist=$prefix.ihist.tsv" \
    "qahist=$prefix.qahist.tsv" \
    "indelhist=$prefix.indelhist.tsv" \
    "ehist=$prefix.ehist.tsv" \
    "maxhistlen=1k" \
    "entropybins=100" \
    "idbins=auto" \
    >"$prefix.stdout.log" 2>"$prefix.stderr.log"

  for suffix in keep1.fq keep2.fq low1.fq low2.fq mid1.fq mid2.fq high1.fq high2.fq; do
    test -e "$prefix.$suffix"
  done
  for suffix in \
    hist.tsv rhist.tsv qhist.tsv bqhist.tsv qchist.tsv aqhist.tsv obqhist.tsv \
    lhist.tsv gchist.tsv bhist.tsv enthist.tsv idhist.tsv mhist.tsv ihist.tsv \
    qahist.tsv indelhist.tsv ehist.tsv; do
    test -s "$prefix.$suffix"
  done
  grep -q 'mhist=.*side-output match histogram' "$prefix.stderr.log"
  grep -q 'qahist=.*side-output quality-accuracy histogram' "$prefix.stderr.log"
  grep -q 'enthist=.*side-output entropy histogram' "$prefix.stderr.log"
}

labels=()
for threads in $THREAD_CASES; do
  label="$(label_for_threads "$threads")"
  labels+=("$label")
  run_case "$threads" "$label"
done

"$PYTHON" - "$OUT" "$READS" "${labels[@]}" <<'PY'
import sys
from pathlib import Path

out = Path(sys.argv[1])
expected_pairs = int(sys.argv[2])
labels = sys.argv[3:]
if not labels:
    raise SystemExit('No thread labels were supplied')

def path(label, suffix):
    return out / f'{label}.{suffix}'

def fastq_stats(path):
    text = Path(path).read_text()
    if not text:
        return 0, 0, 0
    lines = text.splitlines()
    if len(lines) % 4:
        raise SystemExit(f'{path} is not a complete FASTQ file')
    lengths = [len(lines[i + 1]) for i in range(0, len(lines), 4)]
    return len(lengths), sum(lengths), max(lengths, default=0)

def table_rows(path):
    for line in Path(path).read_text().splitlines():
        if line and not line.startswith('#'):
            yield line.split('\t')

def sum_two_column_hist(path):
    total = 0
    for fields in table_rows(path):
        total += int(fields[1])
    return total

def sum_three_column_hist(path):
    reads = 0
    bases = 0
    for fields in table_rows(path):
        reads += int(fields[1])
        bases += int(fields[2])
    return reads, bases

compare_suffixes = [
    'keep1.fq', 'keep2.fq',
    'low1.fq', 'low2.fq',
    'mid1.fq', 'mid2.fq',
    'high1.fq', 'high2.fq',
    'hist.tsv', 'rhist.tsv',
    'qhist.tsv', 'bqhist.tsv', 'qchist.tsv', 'aqhist.tsv', 'obqhist.tsv',
    'lhist.tsv', 'gchist.tsv', 'bhist.tsv', 'enthist.tsv', 'idhist.tsv',
    'mhist.tsv', 'ihist.tsv', 'qahist.tsv', 'indelhist.tsv', 'ehist.tsv',
]

summary = ['metric\tvalue']
for label in labels:
    r1_reads, r1_bases, r1_max = fastq_stats(path(label, 'keep1.fq'))
    r2_reads, r2_bases, r2_max = fastq_stats(path(label, 'keep2.fq'))
    if r1_reads != expected_pairs or r2_reads != expected_pairs:
        raise SystemExit(f'{label} kept {r1_reads}/{r2_reads} pairs, expected {expected_pairs}')
    total_reads = r1_reads + r2_reads
    total_bases = r1_bases + r2_bases

    if sum_two_column_hist(path(label, 'qhist.tsv')) != total_bases:
        raise SystemExit(f'{label} qhist does not sum to {total_bases} bases')
    lhist_reads, lhist_bases = sum_three_column_hist(path(label, 'lhist.tsv'))
    if (lhist_reads, lhist_bases) != (total_reads, total_bases):
        raise SystemExit(f'{label} lhist saw {lhist_reads}/{lhist_bases}, expected {total_reads}/{total_bases}')
    gchist_reads, gchist_bases = sum_three_column_hist(path(label, 'gchist.tsv'))
    if (gchist_reads, gchist_bases) != (total_reads, total_bases):
        raise SystemExit(f'{label} gchist saw {gchist_reads}/{gchist_bases}, expected {total_reads}/{total_bases}')

    qchist_r1 = qchist_r2 = 0
    for fields in table_rows(path(label, 'qchist.tsv')):
        if len(fields) != 5:
            raise SystemExit(f'{label} qchist should have paired rows: {fields}')
        qchist_r1 += int(fields[1])
        qchist_r2 += int(fields[3])
    if (qchist_r1, qchist_r2) != (r1_bases, r2_bases):
        raise SystemExit(f'{label} qchist saw {qchist_r1}/{qchist_r2}, expected {r1_bases}/{r2_bases}')

    aqhist_r1 = aqhist_r2 = 0
    for fields in table_rows(path(label, 'aqhist.tsv')):
        if len(fields) != 5:
            raise SystemExit(f'{label} aqhist should have paired rows: {fields}')
        aqhist_r1 += int(fields[1])
        aqhist_r2 += int(fields[3])
    if (aqhist_r1, aqhist_r2) != (r1_reads, r2_reads):
        raise SystemExit(f'{label} aqhist saw {aqhist_r1}/{aqhist_r2}, expected {r1_reads}/{r2_reads}')

    if sum(int(fields[1]) for fields in table_rows(path(label, 'obqhist.tsv'))) != total_bases:
        raise SystemExit(f'{label} obqhist does not sum to {total_bases} bases')

    bhist_rows = 0
    for fields in table_rows(path(label, 'bhist.tsv')):
        if len(fields) != 6:
            raise SystemExit(f'{label} bhist expected six fields, saw {fields}')
        fractions = [float(value) for value in fields[1:]]
        if not 0.9999 <= sum(fractions) <= 1.0001:
            raise SystemExit(f'{label} bhist fractions do not sum to 1.0: {fields}')
        bhist_rows += 1
    if bhist_rows == 0 or bhist_rows > r1_max + r2_max:
        raise SystemExit(f'{label} bhist row count {bhist_rows} is outside expected read-length bounds')

    bqhist_rows = 0
    for fields in table_rows(path(label, 'bqhist.tsv')):
        if len(fields) != 19:
            raise SystemExit(f'{label} bqhist expected 19 paired fields, saw {fields}')
        if int(fields[1]) < 1 or int(fields[10]) < 1:
            raise SystemExit(f'{label} bqhist paired row should have positive observations: {fields}')
        bqhist_rows += 1
    if bqhist_rows == 0 or bqhist_rows > max(r1_max, r2_max):
        raise SystemExit(f'{label} bqhist row count {bqhist_rows} is outside expected read-length bounds')

    entropy_reads = sum(int(fields[1]) for fields in table_rows(path(label, 'enthist.tsv')))
    if entropy_reads != total_reads:
        raise SystemExit(f'{label} enthist saw {entropy_reads}, expected {total_reads} reads')

    id_reads = id_bases = 0
    for fields in table_rows(path(label, 'idhist.tsv')):
        id_reads += int(fields[1])
        id_bases += int(fields[2])
    if (id_reads, id_bases) != (total_reads, total_bases):
        raise SystemExit(f'{label} idhist saw {id_reads}/{id_bases}, expected {total_reads}/{total_bases}')

    mhist_rows = 0
    for fields in table_rows(path(label, 'mhist.tsv')):
        if len(fields) != 13:
            raise SystemExit(f'{label} mhist expected 13 paired fields, saw {fields}')
        first = [float(value) for value in fields[1:7]]
        second = [float(value) for value in fields[7:13]]
        if not (0.9999 <= sum(first) <= 1.0001 and 0.9999 <= sum(second) <= 1.0001):
            raise SystemExit(f'{label} mhist fractions should sum to 1.0 per mate: {fields}')
        mhist_rows += 1
    if mhist_rows == 0 or mhist_rows > max(r1_max, r2_max):
        raise SystemExit(f'{label} mhist row count {mhist_rows} is outside expected read-length bounds')

    qahist_matches = sum(int(fields[1]) for fields in table_rows(path(label, 'qahist.tsv')))
    if qahist_matches != total_bases:
        raise SystemExit(f'{label} qahist saw {qahist_matches}, expected {total_bases} matching observations')

    ehist_zero = sum(int(fields[1]) for fields in table_rows(path(label, 'ehist.tsv')) if int(fields[0]) == 0)
    if ehist_zero != total_reads:
        raise SystemExit(f'{label} ehist zero-error count {ehist_zero}, expected {total_reads}')

    summary.extend([
        f'{label}_pairs\t{r1_reads}',
        f'{label}_bases\t{total_bases}',
        f'{label}_bqhist_rows\t{bqhist_rows}',
        f'{label}_mhist_rows\t{mhist_rows}',
    ])

baseline = labels[0]
for label in labels[1:]:
    for suffix in compare_suffixes:
        if path(baseline, suffix).read_bytes() != path(label, suffix).read_bytes():
            raise SystemExit(f'Side-output biological stress changed across threads for {suffix}: {baseline} vs {label}')

summary.append(f'thread_cases\t{" ".join(labels)}')
summary.append(f'thread_compare_baseline\t{baseline}')
(out / 'summary.tsv').write_text('\n'.join(summary) + '\n')
PY

printf 'Biological side-output stats stress passed. Summary: %s\n' "$OUT/summary.tsv"

#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT"

BIO_ROOT="${BIO_ROOT:-/home/jake/Projects/biological data}"
DATA1="${DATA1:-$BIO_ROOT/reads/short_reads/scer_s288c_pe_srr23631023/SRR23631023_1.fastq.gz}"
DATA2="${DATA2:-$BIO_ROOT/reads/short_reads/scer_s288c_pe_srr23631023/SRR23631023_2.fastq.gz}"
OUT="${1:-tmp/side_routing_biological_stress}"
READS="${READS:-1000}"
TABLE_READS="${TABLE_READS:-$READS}"
THREAD_CASES="${THREAD_CASES:-${THREADS:-1 2 auto}}"
UNCORRECTED_CASES="${UNCORRECTED_CASES:-1}"
MIN_READ_LEN="${MIN_READ_LEN:-40}"
PYTHON="${PYTHON:-python}"
mkdir -p "$OUT"
rm -f "$OUT"/*

if [[ ! -f "$DATA1" || ! -f "$DATA2" ]]; then
  printf 'Missing paired dataset files:\n  DATA1=%s\n  DATA2=%s\n' "$DATA1" "$DATA2" >&2
  exit 2
fi
for value_name in READS TABLE_READS UNCORRECTED_CASES; do
  value="${!value_name}"
  if [[ ! "$value" =~ ^[0-9]+$ || "$value" -lt 1 ]]; then
    printf '%s must be a positive integer; got %s\n' "$value_name" "$value" >&2
    exit 2
  fi
done
if [[ ! "$MIN_READ_LEN" =~ ^[0-9]+$ || "$MIN_READ_LEN" -lt 34 ]]; then
  printf 'MIN_READ_LEN must be an integer >=34 for k=31 ECC stress; got %s\n' "$MIN_READ_LEN" >&2
  exit 2
fi

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
  "k=31"
  "minq=0"
  "minprob=0"
  "min=1"
  "minkmers=1"
  "target=1"
  "max=1"
  "overwrite=t"
  "bits=32"
)

run_side_case() {
  local threads="$1"
  local label="$2"
  local prefix="$OUT/$label"

  printf 'Running paired side-routing biological stress on %s and %s with reads=%s threads=%s...\n' "$DATA1" "$DATA2" "$READS" "$threads"
  target/debug/bbnorm-rs \
    "${COMMON[@]}" \
    "threads=$threads" \
    "out=$prefix.keep1.fq" "out2=$prefix.keep2.fq" \
    "outt=$prefix.toss1.fq" "outt2=$prefix.toss2.fq" \
    "outlow=$prefix.low1.fq" "outlow2=$prefix.low2.fq" \
    "outmid=$prefix.mid1.fq" "outmid2=$prefix.mid2.fq" \
    "outhigh=$prefix.high1.fq" "outhigh2=$prefix.high2.fq" \
    "hist=$prefix.hist.tsv" "rhist=$prefix.rhist.tsv" \
    >"$prefix.stdout.log" 2>"$prefix.stderr.log"

  test -s "$prefix.hist.tsv"
  test -s "$prefix.rhist.tsv"
  for path in \
    "$prefix.keep1.fq" "$prefix.keep2.fq" \
    "$prefix.toss1.fq" "$prefix.toss2.fq" \
    "$prefix.low1.fq" "$prefix.low2.fq" \
    "$prefix.mid1.fq" "$prefix.mid2.fq" \
    "$prefix.high1.fq" "$prefix.high2.fq"; do
    test -e "$path"
  done
}

printf 'Building real-derived paired uncorrected-routing fixture from %s and %s...\n' "$DATA1" "$DATA2"
"$PYTHON" - "$DATA1" "$DATA2" "$OUT" "$UNCORRECTED_CASES" "$MIN_READ_LEN" <<'PY'
import gzip
import sys
from pathlib import Path

data1, data2, out, wanted, min_read_len = (
    sys.argv[1],
    sys.argv[2],
    Path(sys.argv[3]),
    int(sys.argv[4]),
    int(sys.argv[5]),
)

def open_text(path):
    if path.lower().endswith('.gz'):
        return gzip.open(path, 'rt')
    return open(path, 'rt', encoding='utf-8')

def records(path):
    with open_text(path) as handle:
        while True:
            header = handle.readline().rstrip('\n')
            if not header:
                return
            seq = handle.readline().rstrip('\n').upper()
            plus = handle.readline().rstrip('\n')
            qual = handle.readline().rstrip('\n')
            yield header, seq, plus, qual

def mutate(seq, pos):
    for base in 'ACGT':
        if base != seq[pos]:
            return seq[:pos] + base + seq[pos + 1:]
    raise AssertionError('unreachable')

selected = []
for idx, (r1, r2) in enumerate(zip(records(data1), records(data2)), start=1):
    _, seq1, _, _ = r1
    _, seq2, _, _ = r2
    if len(seq1) >= min_read_len and len(seq2) >= min_read_len and set(seq1) <= set('ACGT') and set(seq2) <= set('ACGT'):
        selected.append((idx, seq1, seq2))
        if len(selected) == wanted:
            break
if len(selected) < wanted:
    raise SystemExit(f'Could not find {wanted} no-N paired reads at least {min_read_len} bp long')

meta = ['case\tselected_record\tclean1\tclean2\tmutant1']
with (out / 'uncorrected.1.fq').open('w') as r1_out, (out / 'uncorrected.2.fq').open('w') as r2_out:
    for case, (idx, seq1, seq2) in enumerate(selected, start=1):
        pos = min(max(0, len(seq1) // 2), len(seq1) - 1)
        mutant1 = mutate(seq1, pos)
        qual1 = 'I' * len(seq1)
        qual2 = 'I' * len(seq2)
        for i in range(30):
            r1_out.write(f'@case{case}_clean{i + 1:03d}/1\n{seq1}\n+\n{qual1}\n')
            r2_out.write(f'@case{case}_clean{i + 1:03d}/2\n{seq2}\n+\n{qual2}\n')
        r1_out.write(f'@case{case}_uncorrectable/1\n{mutant1}\n+\n{qual1}\n')
        r2_out.write(f'@case{case}_uncorrectable/2\n{seq2}\n+\n{qual2}\n')
        meta.append(f'{case}\t{idx}\t{seq1}\t{seq2}\t{mutant1}')
(out / 'uncorrected_meta.tsv').write_text('\n'.join(meta) + '\n')
PY

run_uncorrected_case() {
  local threads="$1"
  local label="$2"
  local prefix="$OUT/$label.uncorrected"

  printf 'Running paired ECC outuncorrected biological-derived stress with threads=%s...\n' "$threads"
  target/debug/bbnorm-rs \
    "in=$OUT/uncorrected.1.fq" \
    "in2=$OUT/uncorrected.2.fq" \
    "passes=1" \
    "keepall=t" \
    "ecc=t" "ecco=f" "eccmaxqual=0" \
    "k=31" "minq=0" "minprob=0" "min=0" "minkmers=1" \
    "target=999999999" "max=999999999" \
    "threads=$threads" "overwrite=t" "bits=32" \
    "out=$prefix.keep1.fq" "out2=$prefix.keep2.fq" \
    "outuncorrected=$prefix.unc1.fq" "outuncorrected2=$prefix.unc2.fq" \
    >"$prefix.stdout.log" 2>"$prefix.stderr.log"

  if grep -q 'error correction option ecc=t is not implemented yet' "$prefix.stderr.log"; then
    echo 'Rust emitted the old ECC fallback note in outuncorrected biological stress' >&2
    exit 1
  fi
  test -s "$prefix.keep1.fq"
  test -s "$prefix.keep2.fq"
  test -s "$prefix.unc1.fq"
  test -s "$prefix.unc2.fq"
}

run_mark_uncorrected_case() {
  local threads="$1"
  local label="$2"
  local prefix="$OUT/$label.markuncorrectable"

  printf 'Running paired ECC markuncorrectable biological-derived stress with threads=%s...\n' "$threads"
  target/debug/bbnorm-rs \
    "in=$OUT/uncorrected.1.fq" \
    "in2=$OUT/uncorrected.2.fq" \
    "passes=1" \
    "keepall=t" \
    "ecc=t" "ecco=f" "eccmaxqual=0" "markuncorrectableerrors=t" \
    "k=31" "minq=0" "minprob=0" "min=0" "minkmers=1" \
    "target=999999999" "max=999999999" \
    "threads=$threads" "overwrite=t" "bits=32" \
    "out=$prefix.keep1.fq" "out2=$prefix.keep2.fq" \
    "outuncorrected=$prefix.unc1.fq" "outuncorrected2=$prefix.unc2.fq" \
    >"$prefix.stdout.log" 2>"$prefix.stderr.log"

  if grep -q 'error correction option ecc=t is not implemented yet' "$prefix.stderr.log"; then
    echo 'Rust emitted the old ECC fallback note in markuncorrectable biological stress' >&2
    exit 1
  fi
  test -s "$prefix.keep1.fq"
  test -s "$prefix.keep2.fq"
  test -s "$prefix.unc1.fq"
  test -s "$prefix.unc2.fq"
}

labels=()
for threads in $THREAD_CASES; do
  label="$(label_for_threads "$threads")"
  labels+=("$label")
  run_side_case "$threads" "$label"
  run_uncorrected_case "$threads" "$label"
  run_mark_uncorrected_case "$threads" "$label"
done

"$PYTHON" - "$OUT" "$READS" "$UNCORRECTED_CASES" "${labels[@]}" <<'PY'
import sys
from pathlib import Path

out = Path(sys.argv[1])
expected = int(sys.argv[2])
uncorrected_cases = int(sys.argv[3])
labels = sys.argv[4:]
if not labels:
    raise SystemExit('No thread labels were supplied')

def records(path):
    text = Path(path).read_text()
    if not text:
        return 0
    lines = text.splitlines()
    if len(lines) % 4:
        raise SystemExit(f'{path} is not a complete FASTQ file')
    return len(lines) // 4

def fastq_map(path):
    lines = Path(path).read_text().splitlines()
    if len(lines) % 4:
        raise SystemExit(f'{path} is not a complete FASTQ file')
    return {
        lines[i][1:]: (lines[i + 1], lines[i + 3])
        for i in range(0, len(lines), 4)
    }

def path(label, suffix):
    return out / f'{label}.{suffix}'

def uncorrected_path(label, suffix):
    return out / f'{label}.uncorrected.{suffix}'

def mark_path(label, suffix):
    return out / f'{label}.markuncorrectable.{suffix}'

compare_suffixes = [
    'keep1.fq', 'keep2.fq',
    'toss1.fq', 'toss2.fq',
    'low1.fq', 'low2.fq',
    'mid1.fq', 'mid2.fq',
    'high1.fq', 'high2.fq',
    'hist.tsv', 'rhist.tsv',
]
unc_compare_suffixes = ['keep1.fq', 'keep2.fq', 'unc1.fq', 'unc2.fq']
summary = ['metric\tvalue']
for label in labels:
    counts = {}
    for stem in ['keep', 'toss', 'low', 'mid', 'high']:
        left = records(path(label, f'{stem}1.fq'))
        right = records(path(label, f'{stem}2.fq'))
        if left != right:
            raise SystemExit(f'Paired count mismatch for {label}.{stem}: {left} vs {right}')
        counts[stem] = left
    if counts['keep'] + counts['toss'] != expected:
        raise SystemExit(f'{label} keep+toss={counts["keep"] + counts["toss"]}, expected {expected}')
    if counts['low'] + counts['mid'] + counts['high'] != expected:
        raise SystemExit(f'{label} low+mid+high={counts["low"] + counts["mid"] + counts["high"]}, expected {expected}')

    unc_keep1 = records(uncorrected_path(label, 'keep1.fq'))
    unc_keep2 = records(uncorrected_path(label, 'keep2.fq'))
    unc1 = records(uncorrected_path(label, 'unc1.fq'))
    unc2 = records(uncorrected_path(label, 'unc2.fq'))
    expected_unc_fixture = uncorrected_cases * 31
    if unc_keep1 != expected_unc_fixture or unc_keep2 != expected_unc_fixture:
        raise SystemExit(f'{label} ECC keepall retained {unc_keep1}/{unc_keep2}, expected {expected_unc_fixture}')
    if unc1 != uncorrected_cases or unc2 != uncorrected_cases:
        raise SystemExit(f'{label} outuncorrected wrote {unc1}/{unc2}, expected {uncorrected_cases}')

    mark_keep1 = records(mark_path(label, 'keep1.fq'))
    mark_keep2 = records(mark_path(label, 'keep2.fq'))
    mark_unc1 = records(mark_path(label, 'unc1.fq'))
    mark_unc2 = records(mark_path(label, 'unc2.fq'))
    if mark_keep1 != expected_unc_fixture or mark_keep2 != expected_unc_fixture:
        raise SystemExit(f'{label} markuncorrectable keepall retained {mark_keep1}/{mark_keep2}, expected {expected_unc_fixture}')
    if mark_unc1 != uncorrected_cases or mark_unc2 != uncorrected_cases:
        raise SystemExit(f'{label} markuncorrectable outuncorrected wrote {mark_unc1}/{mark_unc2}, expected {uncorrected_cases}')

    mark_keep1_map = fastq_map(mark_path(label, 'keep1.fq'))
    mark_keep2_map = fastq_map(mark_path(label, 'keep2.fq'))
    mark_unc1_map = fastq_map(mark_path(label, 'unc1.fq'))
    mark_unc2_map = fastq_map(mark_path(label, 'unc2.fq'))
    for case in range(1, uncorrected_cases + 1):
        r1_id = f'case{case}_uncorrectable/1'
        r2_id = f'case{case}_uncorrectable/2'
        if r1_id not in mark_keep1_map or r1_id not in mark_unc1_map:
            raise SystemExit(f'{label} missing marked uncorrectable read {r1_id}')
        if r2_id not in mark_keep2_map or r2_id not in mark_unc2_map:
            raise SystemExit(f'{label} missing mate read {r2_id}')
        _, keep1_qual = mark_keep1_map[r1_id]
        _, unc1_qual = mark_unc1_map[r1_id]
        _, keep2_qual = mark_keep2_map[r2_id]
        _, unc2_qual = mark_unc2_map[r2_id]
        if keep1_qual == 'I' * len(keep1_qual):
            raise SystemExit(f'{label} left marked uncorrectable qualities unchanged for {r1_id}')
        if unc1_qual != keep1_qual:
            raise SystemExit(f'{label} keep/uncorrected quality mismatch for {r1_id}')
        if keep2_qual != 'I' * len(keep2_qual) or unc2_qual != 'I' * len(unc2_qual):
            raise SystemExit(f'{label} unexpectedly changed clean mate qualities for {r2_id}')

    summary.extend([
        f'{label}_keep_pairs\t{counts["keep"]}',
        f'{label}_toss_pairs\t{counts["toss"]}',
        f'{label}_low_pairs\t{counts["low"]}',
        f'{label}_mid_pairs\t{counts["mid"]}',
        f'{label}_high_pairs\t{counts["high"]}',
        f'{label}_uncorrected_keep_pairs\t{unc_keep1}',
        f'{label}_uncorrected_pairs\t{unc1}',
        f'{label}_markuncorrectable_keep_pairs\t{mark_keep1}',
        f'{label}_markuncorrectable_pairs\t{mark_unc1}',
    ])

baseline = labels[0]
for label in labels[1:]:
    for suffix in compare_suffixes:
        if path(baseline, suffix).read_bytes() != path(label, suffix).read_bytes():
            raise SystemExit(f'Side-routing output changed across threads for {suffix}: {baseline} vs {label}')
    for suffix in unc_compare_suffixes:
        if uncorrected_path(baseline, suffix).read_bytes() != uncorrected_path(label, suffix).read_bytes():
            raise SystemExit(f'Uncorrected output changed across threads for {suffix}: {baseline} vs {label}')
    for suffix in unc_compare_suffixes:
        base = out / f'{baseline}.markuncorrectable.{suffix}'
        other = out / f'{label}.markuncorrectable.{suffix}'
        if base.read_bytes() != other.read_bytes():
            raise SystemExit(f'Markuncorrectable output changed across threads for {suffix}: {baseline} vs {label}')

summary.append(f'thread_cases\t{" ".join(labels)}')
summary.append(f'thread_compare_baseline\t{baseline}')
(out / 'summary.tsv').write_text('\n'.join(summary) + '\n')
PY

printf 'Biological side-routing stress passed. Summary: %s\n' "$OUT/summary.tsv"

#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT"

BIO_ROOT="${BIO_ROOT:-/home/jake/Projects/biological data}"
DATA1="${DATA1:-$BIO_ROOT/reads/short_reads/scer_s288c_pe_srr23631023/SRR23631023_1.fastq.gz}"
DATA2="${DATA2:-$BIO_ROOT/reads/short_reads/scer_s288c_pe_srr23631023/SRR23631023_2.fastq.gz}"
OUT="${1:-tmp/countup_biological_stress}"
READS="${READS:-1000}"
TABLE_READS="${TABLE_READS:-$READS}"
THREAD_CASES="${THREAD_CASES:-${THREADS:-auto}}"
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
  "countup=t"
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

printf 'Building real-derived paired countup markuncorrectable fixture from %s and %s...\n' "$DATA1" "$DATA2"
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
with (out / 'markuncorrectable.1.fq').open('w') as r1_out, (out / 'markuncorrectable.2.fq').open('w') as r2_out:
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
(out / 'markuncorrectable_meta.tsv').write_text('\n'.join(meta) + '\n')
PY

run_case() {
  local threads="$1"
  local label="$2"
  local prefix="$OUT/$label"

  printf 'Running Rust countup=t biological stress on %s and %s with reads=%s threads=%s...\n' "$DATA1" "$DATA2" "$READS" "$threads"
  target/debug/bbnorm-rs \
    "${COMMON[@]}" \
    "threads=$threads" \
    "addbadreadscountup=t" \
    "rename=t" \
    "out=$prefix.countup.keep1.fq" "out2=$prefix.countup.keep2.fq" \
    "outt=$prefix.countup.toss1.fq" "outt2=$prefix.countup.toss2.fq" \
    "outlow=$prefix.countup.low1.fq" "outlow2=$prefix.countup.low2.fq" \
    "outmid=$prefix.countup.mid1.fq" "outmid2=$prefix.countup.mid2.fq" \
    "outhigh=$prefix.countup.high1.fq" "outhigh2=$prefix.countup.high2.fq" \
    "outuncorrected=$prefix.countup.unc1.fq" "outuncorrected2=$prefix.countup.unc2.fq" \
    "hist=$prefix.countup.hist.tsv" "rhist=$prefix.countup.rhist.tsv" \
    >"$prefix.countup.stdout.log" 2>"$prefix.countup.stderr.log"

  if grep -qi 'countup.*not implemented\|not implemented.*countup' "$prefix.countup.stderr.log"; then
    echo 'Rust emitted an old count-up fallback note' >&2
    exit 1
  fi

  test -s "$prefix.countup.hist.tsv"
  test -s "$prefix.countup.rhist.tsv"
  for path in \
    "$prefix.countup.keep1.fq" "$prefix.countup.keep2.fq" \
    "$prefix.countup.toss1.fq" "$prefix.countup.toss2.fq" \
    "$prefix.countup.low1.fq" "$prefix.countup.low2.fq" \
    "$prefix.countup.mid1.fq" "$prefix.countup.mid2.fq" \
    "$prefix.countup.high1.fq" "$prefix.countup.high2.fq" \
    "$prefix.countup.unc1.fq" "$prefix.countup.unc2.fq"; do
    test -e "$path"
  done

  printf 'Running Rust countup=t ecc=t biological smoke with threads=%s...\n' "$threads"
  target/debug/bbnorm-rs \
    "${COMMON[@]}" \
    "threads=$threads" \
    "keepall=t" \
    "target=999999999" "max=999999999" \
    "ecc=t" "ecco=f" "eccmaxqual=99" \
    "out=$prefix.ecc.keep1.fq" "out2=$prefix.ecc.keep2.fq" \
    "outuncorrected=$prefix.ecc.unc1.fq" "outuncorrected2=$prefix.ecc.unc2.fq" \
    "hist=$prefix.ecc.hist.tsv" "rhist=$prefix.ecc.rhist.tsv" \
    >"$prefix.ecc.stdout.log" 2>"$prefix.ecc.stderr.log"

  if grep -q 'error correction option ecc=t is not implemented yet' "$prefix.ecc.stderr.log"; then
    echo 'Rust emitted the old ECC fallback note in count-up mode' >&2
    exit 1
  fi

  test -s "$prefix.ecc.keep1.fq"
  test -s "$prefix.ecc.keep2.fq"
  test -e "$prefix.ecc.unc1.fq"
  test -e "$prefix.ecc.unc2.fq"
  test -s "$prefix.ecc.hist.tsv"
  test -s "$prefix.ecc.rhist.tsv"

  printf 'Running Rust countup=t ecc=t markuncorrectable biological stress with threads=%s...\n' "$threads"
  target/debug/bbnorm-rs \
    "in=$OUT/markuncorrectable.1.fq" \
    "in2=$OUT/markuncorrectable.2.fq" \
    "passes=1" \
    "countup=t" \
    "keepall=t" \
    "addbadreadscountup=t" \
    "ecc=t" "ecco=f" "eccmaxqual=0" "markuncorrectableerrors=t" \
    "k=31" "minq=0" "minprob=0" "min=0" "minkmers=1" \
    "target=999999999" "max=999999999" \
    "threads=$threads" "overwrite=t" "bits=32" \
    "out=$prefix.mark.keep1.fq" "out2=$prefix.mark.keep2.fq" \
    "outuncorrected=$prefix.mark.unc1.fq" "outuncorrected2=$prefix.mark.unc2.fq" \
    "hist=$prefix.mark.hist.tsv" "rhist=$prefix.mark.rhist.tsv" \
    >"$prefix.mark.stdout.log" 2>"$prefix.mark.stderr.log"

  if grep -q 'error correction option ecc=t is not implemented yet' "$prefix.mark.stderr.log"; then
    echo 'Rust emitted the old ECC fallback note in count-up markuncorrectable biological stress' >&2
    exit 1
  fi

  test -s "$prefix.mark.keep1.fq"
  test -s "$prefix.mark.keep2.fq"
  test -s "$prefix.mark.unc1.fq"
  test -s "$prefix.mark.unc2.fq"
  test -s "$prefix.mark.hist.tsv"
  test -s "$prefix.mark.rhist.tsv"
}

labels=()
for threads in $THREAD_CASES; do
  label="$(label_for_threads "$threads")"
  labels+=("$label")
  run_case "$threads" "$label"
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

paired_outputs = [
    'countup.keep',
    'countup.toss',
    'countup.low',
    'countup.mid',
    'countup.high',
    'countup.unc',
    'ecc.keep',
    'ecc.unc',
]
compare_suffixes = [
    'countup.keep1.fq', 'countup.keep2.fq',
    'countup.toss1.fq', 'countup.toss2.fq',
    'countup.low1.fq', 'countup.low2.fq',
    'countup.mid1.fq', 'countup.mid2.fq',
    'countup.high1.fq', 'countup.high2.fq',
    'countup.unc1.fq', 'countup.unc2.fq',
    'countup.hist.tsv', 'countup.rhist.tsv',
    'ecc.keep1.fq', 'ecc.keep2.fq',
    'ecc.unc1.fq', 'ecc.unc2.fq',
    'ecc.hist.tsv', 'ecc.rhist.tsv',
]

def records(path):
    text = Path(path).read_text()
    if not text:
        return 0
    lines = text.splitlines()
    if len(lines) % 4:
        raise SystemExit(f'{path} is not a complete FASTQ file')
    return len(lines) // 4

def first_header(path):
    p = Path(path)
    if not p.exists() or not p.read_text():
        return None
    return p.read_text().splitlines()[0]

def output_path(label, suffix):
    return out / f'{label}.{suffix}'

def fastq_map(path):
    lines = Path(path).read_text().splitlines()
    if len(lines) % 4:
        raise SystemExit(f'{path} is not a complete FASTQ file')
    return {
        lines[i][1:]: (lines[i + 1], lines[i + 3])
        for i in range(0, len(lines), 4)
    }

summary = ['metric\tvalue']
for label in labels:
    counts = {}
    for stem in paired_outputs:
        left = records(output_path(label, f'{stem}1.fq'))
        right = records(output_path(label, f'{stem}2.fq'))
        if left != right:
            raise SystemExit(f'Paired count mismatch for {label}.{stem}: {left} vs {right}')
        counts[stem] = left

    if counts['countup.keep'] + counts['countup.toss'] != expected:
        raise SystemExit(
            f'Expected {label} keep+toss={expected}; got {counts["countup.keep"] + counts["countup.toss"]}'
        )
    if counts['countup.low'] + counts['countup.mid'] + counts['countup.high'] != expected:
        raise SystemExit(
            f'Expected {label} low+mid+high={expected}; '
            f'got {counts["countup.low"] + counts["countup.mid"] + counts["countup.high"]}'
        )
    if counts['ecc.keep'] != expected:
        raise SystemExit(f'Expected {label} countup ECC keepall to retain {expected} pairs; got {counts["ecc.keep"]}')

    mark_keep1 = records(output_path(label, 'mark.keep1.fq'))
    mark_keep2 = records(output_path(label, 'mark.keep2.fq'))
    mark_unc1 = records(output_path(label, 'mark.unc1.fq'))
    mark_unc2 = records(output_path(label, 'mark.unc2.fq'))
    expected_mark_fixture = uncorrected_cases * 31
    if mark_keep1 != expected_mark_fixture or mark_keep2 != expected_mark_fixture:
        raise SystemExit(
            f'Expected {label} countup markuncorrectable keepall to retain '
            f'{expected_mark_fixture} pairs; got {mark_keep1}/{mark_keep2}'
        )
    if mark_unc1 != uncorrected_cases or mark_unc2 != uncorrected_cases:
        raise SystemExit(
            f'Expected {label} countup markuncorrectable outuncorrected to contain '
            f'{uncorrected_cases} pairs; got {mark_unc1}/{mark_unc2}'
        )

    mark_keep1_map = fastq_map(output_path(label, 'mark.keep1.fq'))
    mark_keep2_map = fastq_map(output_path(label, 'mark.keep2.fq'))
    mark_unc1_map = fastq_map(output_path(label, 'mark.unc1.fq'))
    mark_unc2_map = fastq_map(output_path(label, 'mark.unc2.fq'))
    for case in range(1, uncorrected_cases + 1):
        r1_id = f'case{case}_uncorrectable/1'
        r2_id = f'case{case}_uncorrectable/2'
        if r1_id not in mark_keep1_map or r1_id not in mark_unc1_map:
            raise SystemExit(f'{label} missing marked countup uncorrectable read {r1_id}')
        if r2_id not in mark_keep2_map or r2_id not in mark_unc2_map:
            raise SystemExit(f'{label} missing marked countup mate read {r2_id}')
        _, keep1_qual = mark_keep1_map[r1_id]
        _, unc1_qual = mark_unc1_map[r1_id]
        _, keep2_qual = mark_keep2_map[r2_id]
        _, unc2_qual = mark_unc2_map[r2_id]
        if keep1_qual == 'I' * len(keep1_qual):
            raise SystemExit(f'{label} left countup markuncorrectable qualities unchanged for {r1_id}')
        if unc1_qual != keep1_qual:
            raise SystemExit(f'{label} countup keep/uncorrected quality mismatch for {r1_id}')
        if keep2_qual != 'I' * len(keep2_qual) or unc2_qual != 'I' * len(unc2_qual):
            raise SystemExit(f'{label} unexpectedly changed clean mate qualities for {r2_id}')

    renamed_header = None
    for suffix in [
        'countup.keep1.fq',
        'countup.toss1.fq',
        'countup.low1.fq',
        'countup.mid1.fq',
        'countup.high1.fq',
    ]:
        renamed_header = first_header(output_path(label, suffix))
        if renamed_header:
            break
    if not renamed_header or not renamed_header.startswith('@id=') or 'd1=' not in renamed_header:
        raise SystemExit(f'Expected a renamed count-up output header with id/d1 fields for {label}')

    for key in ['countup.keep', 'countup.toss', 'countup.low', 'countup.mid', 'countup.high', 'countup.unc', 'ecc.keep', 'ecc.unc']:
        summary.append(f'{label}_{key.replace(".", "_")}_pairs\t{counts[key]}')
    summary.append(f'{label}_mark_keep_pairs\t{mark_keep1}')
    summary.append(f'{label}_mark_unc_pairs\t{mark_unc1}')
    summary.append(f'{label}_renamed_header\t{renamed_header}')

baseline = labels[0]
for label in labels[1:]:
    for suffix in compare_suffixes:
        left = output_path(baseline, suffix).read_bytes()
        right = output_path(label, suffix).read_bytes()
        if left != right:
            raise SystemExit(f'Output changed across thread cases for {suffix}: {baseline} vs {label}')
    for suffix in ['mark.keep1.fq', 'mark.keep2.fq', 'mark.unc1.fq', 'mark.unc2.fq', 'mark.hist.tsv', 'mark.rhist.tsv']:
        left = output_path(baseline, suffix).read_bytes()
        right = output_path(label, suffix).read_bytes()
        if left != right:
            raise SystemExit(f'Countup markuncorrectable output changed across thread cases for {suffix}: {baseline} vs {label}')
summary.append(f'thread_cases\t{" ".join(labels)}')
summary.append(f'thread_compare_baseline\t{baseline}')
(out / 'summary.tsv').write_text('\n'.join(summary) + '\n')
PY

printf 'Count-up biological stress passed. Summary: %s\n' "$OUT/summary.tsv"

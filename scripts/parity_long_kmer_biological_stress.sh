#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT"

BIO_ROOT="${BIO_ROOT:-/home/jake/Projects/biological data}"
DATA1="${DATA1:-$BIO_ROOT/reads/short_reads/scer_s288c_pe_srr23631023/SRR23631023_1.fastq.gz}"
DATA2="${DATA2:-$BIO_ROOT/reads/short_reads/scer_s288c_pe_srr23631023/SRR23631023_2.fastq.gz}"
OUT="${1:-tmp/long_kmer_biological_stress}"
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
  "k=40"
  "minq=0"
  "minprob=0"
  "min=1"
  "minkmers=1"
  "target=1"
  "max=1"
  "overwrite=t"
  "bits=32"
)

run_mode() {
  local mode="$1"
  local threads="$2"
  local label="$3"
  local prefix="$OUT/$label.$mode"
  local extra=()
  if [[ "$mode" == "fixspikes" ]]; then
    extra=("fixspikes=t")
  fi

  printf 'Running paired biological long-kmer mode=%s reads=%s threads=%s...\n' "$mode" "$READS" "$threads"
  target/debug/bbnorm-rs \
    "${COMMON[@]}" \
    "${extra[@]}" \
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

labels=()
for threads in $THREAD_CASES; do
  label="$(label_for_threads "$threads")"
  labels+=("$label")
  run_mode "plain" "$threads" "$label"
  run_mode "fixspikes" "$threads" "$label"
done

"$PYTHON" - "$OUT" "$READS" "${labels[@]}" <<'PY'
import sys
from pathlib import Path

out = Path(sys.argv[1])
expected = int(sys.argv[2])
labels = sys.argv[3:]
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

def path(label, mode, suffix):
    return out / f'{label}.{mode}.{suffix}'

modes = ['plain', 'fixspikes']
compare_suffixes = [
    'keep1.fq', 'keep2.fq',
    'toss1.fq', 'toss2.fq',
    'low1.fq', 'low2.fq',
    'mid1.fq', 'mid2.fq',
    'high1.fq', 'high2.fq',
    'hist.tsv', 'rhist.tsv',
]
summary = ['metric\tvalue']
for label in labels:
    for mode in modes:
        counts = {}
        for stem in ['keep', 'toss', 'low', 'mid', 'high']:
            left = records(path(label, mode, f'{stem}1.fq'))
            right = records(path(label, mode, f'{stem}2.fq'))
            if left != right:
                raise SystemExit(f'Paired count mismatch for {label}.{mode}.{stem}: {left} vs {right}')
            counts[stem] = left
        if counts['keep'] + counts['toss'] != expected:
            raise SystemExit(f'{label}.{mode} keep+toss={counts["keep"] + counts["toss"]}, expected {expected}')
        if counts['low'] + counts['mid'] + counts['high'] != expected:
            raise SystemExit(f'{label}.{mode} low+mid+high={counts["low"] + counts["mid"] + counts["high"]}, expected {expected}')
        hist = path(label, mode, 'hist.tsv').read_text()
        rhist = path(label, mode, 'rhist.tsv').read_text()
        if '#Depth\tRaw_Count\tUnique_Kmers' not in hist:
            raise SystemExit(f'{label}.{mode} missing k-mer histogram header')
        if '#Depth\tReads\tBases' not in rhist:
            raise SystemExit(f'{label}.{mode} missing read-depth histogram header')
        summary.extend([
            f'{label}_{mode}_keep_pairs\t{counts["keep"]}',
            f'{label}_{mode}_toss_pairs\t{counts["toss"]}',
            f'{label}_{mode}_low_pairs\t{counts["low"]}',
            f'{label}_{mode}_mid_pairs\t{counts["mid"]}',
            f'{label}_{mode}_high_pairs\t{counts["high"]}',
        ])

baseline = labels[0]
for label in labels[1:]:
    for mode in modes:
        for suffix in compare_suffixes:
            if path(baseline, mode, suffix).read_bytes() != path(label, mode, suffix).read_bytes():
                raise SystemExit(f'Long-kmer output changed across threads for {mode}.{suffix}: {baseline} vs {label}')

summary.append(f'thread_cases\t{" ".join(labels)}')
summary.append(f'thread_compare_baseline\t{baseline}')
(out / 'summary.tsv').write_text('\n'.join(summary) + '\n')
PY

printf 'Biological long-kmer stress passed. Summary: %s\n' "$OUT/summary.tsv"

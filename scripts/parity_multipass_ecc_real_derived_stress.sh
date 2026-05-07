#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT"

DATA1="${DATA1:-vendor/BBTools-master/resources/sample1.fq.gz}"
DATA2="${DATA2:-vendor/BBTools-master/resources/sample2.fq.gz}"
OUT="${1:-tmp/multipass_ecc_real_derived_stress}"
STRESS_RECORDS="${STRESS_RECORDS:-1}"
THREAD_CASES="${THREAD_CASES:-${THREADS:-1}}"
MIN_READ_LEN="${MIN_READ_LEN:-40}"
mkdir -p "$OUT"
rm -f "$OUT"/*

if [[ ! -f "$DATA1" || ! -f "$DATA2" ]]; then
  printf 'Missing paired dataset files:\n  DATA1=%s\n  DATA2=%s\n' "$DATA1" "$DATA2" >&2
  exit 2
fi
if [[ ! "$STRESS_RECORDS" =~ ^[0-9]+$ || "$STRESS_RECORDS" -lt 1 ]]; then
  printf 'STRESS_RECORDS must be a positive integer; got %s\n' "$STRESS_RECORDS" >&2
  exit 2
fi
if [[ ! "$MIN_READ_LEN" =~ ^[0-9]+$ || "$MIN_READ_LEN" -lt 34 ]]; then
  printf 'MIN_READ_LEN must be an integer >=34 for k=31 ECC stress; got %s\n' "$MIN_READ_LEN" >&2
  exit 2
fi

label_for_threads() {
  printf '%s' "$1" | tr -c 'A-Za-z0-9_-' '_'
}

cargo build --quiet

printf 'Building a no-N paired real-derived multipass ECC stress fixture from %s and %s...\n' "$DATA1" "$DATA2"
python - "$DATA1" "$DATA2" "$OUT" "$STRESS_RECORDS" "$MIN_READ_LEN" <<'PY'
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
            seq = handle.readline().rstrip('\n')
            plus = handle.readline().rstrip('\n')
            qual = handle.readline().rstrip('\n')
            yield header, seq.upper(), plus, qual

def mutate(seq, pos):
    for base in 'ACGT':
        if base != seq[pos]:
            return seq[:pos] + base + seq[pos + 1:]
    raise AssertionError('unreachable')

selected = []
for idx, (r1, r2) in enumerate(zip(records(data1), records(data2)), start=1):
    _, seq1, _, qual1 = r1
    _, seq2, _, qual2 = r2
    if len(seq1) >= min_read_len and len(seq2) >= min_read_len and set(seq1) <= set('ACGT') and set(seq2) <= set('ACGT'):
        selected.append((idx, seq1, seq2, qual1, qual2))
        if len(selected) == wanted:
            break
if len(selected) < wanted:
    raise SystemExit(
        f'Could not find {wanted} no-N paired reads at least {min_read_len} bp long for the multipass ECC fixture'
    )

meta_lines = ['case\tselected_record\tpos1\tpos2\tclean1\tclean2\tmutants1\tmutants2\tuncorrectable1']
with (
    (out / 'stress.1.fq').open('w') as r1_out,
    (out / 'stress.2.fq').open('w') as r2_out,
    (out / 'uncorrectable.1.fq').open('w') as unc1_out,
    (out / 'uncorrectable.2.fq').open('w') as unc2_out,
):
    for case, (idx, seq1, seq2, qual1, qual2) in enumerate(selected, start=1):
        pos1 = min(50, max(0, len(seq1) // 2))
        pos2 = min(50, max(0, len(seq2) // 2))
        pos1 = min(pos1, len(seq1) - 3)
        pos2 = min(pos2, len(seq2) - 3)
        mut1 = [mutate(seq1, pos1 + offset) for offset in range(3)]
        mut2 = [mutate(seq2, pos2 + offset) for offset in range(3)]
        uncorrectable1 = mutate(seq1, pos1)
        for i in range(40):
            r1_out.write(f'@case{case}_clean{i + 1:03d}/1\n{seq1}\n+\n{qual1}\n')
            r2_out.write(f'@case{case}_clean{i + 1:03d}/2\n{seq2}\n+\n{qual2}\n')
            unc1_out.write(f'@case{case}_clean{i + 1:03d}/1\n{seq1}\n+\n{"I" * len(seq1)}\n')
            unc2_out.write(f'@case{case}_clean{i + 1:03d}/2\n{seq2}\n+\n{"I" * len(seq2)}\n')
        for i, seq in enumerate(mut1, start=1):
            r1_out.write(f'@case{case}_mut_r1_{i}/1\n{seq}\n+\n{qual1}\n')
            r2_out.write(f'@case{case}_mut_r1_{i}/2\n{seq2}\n+\n{qual2}\n')
        for i, seq in enumerate(mut2, start=1):
            r1_out.write(f'@case{case}_mut_r2_{i}/1\n{seq1}\n+\n{qual1}\n')
            r2_out.write(f'@case{case}_mut_r2_{i}/2\n{seq}\n+\n{qual2}\n')
        unc1_out.write(f'@case{case}_uncorrectable/1\n{uncorrectable1}\n+\n{"I" * len(seq1)}\n')
        unc2_out.write(f'@case{case}_uncorrectable/2\n{seq2}\n+\n{"I" * len(seq2)}\n')
        meta_lines.append(
            f'{case}\t{idx}\t{pos1}\t{pos2}\t{seq1}\t{seq2}\t{",".join(mut1)}\t{",".join(mut2)}\t{uncorrectable1}'
        )

(out / 'fixture_meta.tsv').write_text('\n'.join(meta_lines) + '\n')
PY

COMMON=(
  "in=$OUT/stress.1.fq"
  "in2=$OUT/stress.2.fq"
  "passes=2"
  "ecco=f"
  "eccmaxqual=99"
  "k=31"
  "minq=0"
  "minprob=0"
  "min=0"
  "minkmers=1"
  "target=999999999"
  "max=999999999"
  "keepall=t"
  "overwrite=t"
  "bits=32"
)

printf 'case\tthread_label\tthreads\toutput1\toutput2\n' > "$OUT/run_matrix.tsv"
base_label=""

for thread_case in $THREAD_CASES; do
  thread_label="threads_$(label_for_threads "$thread_case")"
  [[ -n "$base_label" ]] || base_label="$thread_label"

  for case_name in final_only first_only both_passes; do
    case "$case_name" in
      final_only) ECC_STAGE=("ecc1=f" "eccf=t") ;;
      first_only) ECC_STAGE=("ecc1=t" "eccf=f") ;;
      both_passes) ECC_STAGE=("ecc=t") ;;
      *) echo "unknown multipass ECC case: $case_name" >&2; exit 2 ;;
    esac

    printf 'Running Rust multipass ECC case=%s threads=%s...\n' "$case_name" "$thread_case"
    target/debug/bbnorm-rs \
      "${COMMON[@]}" \
      "${ECC_STAGE[@]}" \
      "threads=$thread_case" \
      "out=$OUT/rust.$case_name.$thread_label.1.fq" \
      "out2=$OUT/rust.$case_name.$thread_label.2.fq" \
      >"$OUT/rust.$case_name.$thread_label.stdout.log" \
      2>"$OUT/rust.$case_name.$thread_label.stderr.log"

    if grep -q 'error correction option ecc=t is not implemented yet' "$OUT/rust.$case_name.$thread_label.stderr.log"; then
      echo "Rust emitted the old ECC fallback note for $case_name/$thread_case" >&2
      exit 1
    fi
    printf '%s\t%s\t%s\t%s\t%s\n' \
      "$case_name" "$thread_label" "$thread_case" \
      "$OUT/rust.$case_name.$thread_label.1.fq" \
      "$OUT/rust.$case_name.$thread_label.2.fq" >> "$OUT/run_matrix.tsv"
  done
done

python - "$OUT" "$base_label" <<'PY'
import sys
from pathlib import Path

out = Path(sys.argv[1])
base_label = sys.argv[2]
meta_rows = []
for line in (out / 'fixture_meta.tsv').read_text().splitlines()[1:]:
    case, selected_record, pos1, pos2, clean1, clean2, mutants1, mutants2, uncorrectable1 = line.split('\t')
    meta_rows.append({
        'case': case,
        'selected_record': selected_record,
        'clean1': clean1,
        'clean2': clean2,
        'mutants1': mutants1.split(','),
        'mutants2': mutants2.split(','),
        'uncorrectable1': uncorrectable1,
    })

runs = []
for line in (out / 'run_matrix.tsv').read_text().splitlines()[1:]:
    case_name, thread_label, threads, output1, output2 = line.split('\t')
    runs.append((case_name, thread_label, threads, Path(output1), Path(output2)))

if not runs:
    raise SystemExit('No multipass ECC stress runs were recorded')

def count_records(text):
    lines = text.splitlines()
    if len(lines) % 4 != 0:
        raise SystemExit('FASTQ output line count is not divisible by 4')
    return len(lines) // 4

def count_any(text, seqs):
    return sum(text.count(seq) for seq in seqs)

summary = ['metric\tvalue']
expected_records = len(meta_rows) * 46
base_outputs = {}

for case_name, thread_label, threads, output1, output2 in runs:
    text1 = output1.read_text()
    text2 = output2.read_text()
    records1 = count_records(text1)
    records2 = count_records(text2)
    if records1 != expected_records or records2 != expected_records:
        raise SystemExit(
            f'{case_name}/{threads} wrote {records1}/{records2} mate records, expected {expected_records}'
        )

    total_mutants = 0
    for row in meta_rows:
        clean1 = row['clean1']
        clean2 = row['clean2']
        mutants1 = row['mutants1']
        mutants2 = row['mutants2']
        mutant_count = count_any(text1, mutants1) + count_any(text2, mutants2)
        total_mutants += mutant_count
        if mutant_count != 0:
            raise SystemExit(f'{case_name}/{threads} left mutant sequences in case {row["case"]}')
        clean_count1 = text1.count(clean1)
        clean_count2 = text2.count(clean2)
        if clean_count1 != 46 or clean_count2 != 46:
            raise SystemExit(
                f'{case_name}/{threads} did not restore all reads to clean sequences in case {row["case"]}: '
                f'clean1={clean_count1} clean2={clean_count2}'
            )

    summary.extend([
        f'{case_name}_{thread_label}_records1\t{records1}',
        f'{case_name}_{thread_label}_records2\t{records2}',
        f'{case_name}_{thread_label}_mutants\t{total_mutants}',
    ])

    key = case_name
    if thread_label == base_label:
        base_outputs[key] = (text1, text2)
    elif key in base_outputs and base_outputs[key] != (text1, text2):
        raise SystemExit(f'{case_name}/{threads} output differs from base thread case {base_label}')

summary.extend([
    f'total_cases\t{len(meta_rows)}',
    f'expected_records_per_mate\t{expected_records}',
    f'base_thread_label\t{base_label}',
])
(out / 'summary.tsv').write_text('\n'.join(summary) + '\n')
PY

printf 'Real-derived multipass ECC stress passed. Summary: %s\n' "$OUT/summary.tsv"

TOSS_COMMON=(
  "in=$OUT/stress.1.fq"
  "in2=$OUT/stress.2.fq"
  "passes=2"
  "ecco=f"
  "eccmaxqual=99"
  "k=31"
  "minq=0"
  "minprob=0"
  "min=0"
  "minkmers=1"
  "target=2"
  "max=2"
  "tossbadreads=t"
  "overwrite=t"
  "bits=32"
)

printf 'case\tthread_label\tthreads\tkeep1\tkeep2\ttoss1\ttoss2\n' > "$OUT/toss_run_matrix.tsv"

for thread_case in $THREAD_CASES; do
  thread_label="threads_$(label_for_threads "$thread_case")"

  for case_name in final_only first_only both_passes; do
    case "$case_name" in
      final_only) ECC_STAGE=("ecc1=f" "eccf=t") ;;
      first_only) ECC_STAGE=("ecc1=t" "eccf=f") ;;
      both_passes) ECC_STAGE=("ecc=t") ;;
      *) echo "unknown multipass ECC toss case: $case_name" >&2; exit 2 ;;
    esac

    printf 'Running Rust multipass ECC toss-routing case=%s threads=%s...\n' "$case_name" "$thread_case"
    target/debug/bbnorm-rs \
      "${TOSS_COMMON[@]}" \
      "${ECC_STAGE[@]}" \
      "threads=$thread_case" \
      "out=$OUT/rust.$case_name.$thread_label.tossmode.keep1.fq" \
      "out2=$OUT/rust.$case_name.$thread_label.tossmode.keep2.fq" \
      "outt=$OUT/rust.$case_name.$thread_label.tossmode.toss1.fq" \
      "outt2=$OUT/rust.$case_name.$thread_label.tossmode.toss2.fq" \
      >"$OUT/rust.$case_name.$thread_label.tossmode.stdout.log" \
      2>"$OUT/rust.$case_name.$thread_label.tossmode.stderr.log"

    printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
      "$case_name" "$thread_label" "$thread_case" \
      "$OUT/rust.$case_name.$thread_label.tossmode.keep1.fq" \
      "$OUT/rust.$case_name.$thread_label.tossmode.keep2.fq" \
      "$OUT/rust.$case_name.$thread_label.tossmode.toss1.fq" \
      "$OUT/rust.$case_name.$thread_label.tossmode.toss2.fq" >> "$OUT/toss_run_matrix.tsv"
  done
done

python - "$OUT" "$base_label" <<'PY'
import sys
from pathlib import Path

out = Path(sys.argv[1])
base_label = sys.argv[2]
meta_rows = []
for line in (out / 'fixture_meta.tsv').read_text().splitlines()[1:]:
    case, selected_record, pos1, pos2, clean1, clean2, mutants1, mutants2, uncorrectable1 = line.split('\t')
    meta_rows.append({
        'case': case,
        'clean1': clean1,
        'clean2': clean2,
        'mutants1': mutants1.split(','),
        'mutants2': mutants2.split(','),
    })

runs = []
for line in (out / 'toss_run_matrix.tsv').read_text().splitlines()[1:]:
    case_name, thread_label, threads, keep1, keep2, toss1, toss2 = line.split('\t')
    runs.append((case_name, thread_label, threads, Path(keep1), Path(keep2), Path(toss1), Path(toss2)))

if not runs:
    raise SystemExit('No multipass ECC toss-routing stress runs were recorded')

def count_records(text):
    lines = text.splitlines()
    if len(lines) % 4 != 0:
        raise SystemExit('FASTQ output line count is not divisible by 4')
    return len(lines) // 4

def count_any(text, seqs):
    return sum(text.count(seq) for seq in seqs)

summary = (out / 'summary.tsv').read_text().splitlines()
expected_records = len(meta_rows) * 46
base_outputs = {}

for case_name, thread_label, threads, keep1, keep2, toss1, toss2 in runs:
    keep_text1 = keep1.read_text()
    keep_text2 = keep2.read_text()
    toss_text1 = toss1.read_text()
    toss_text2 = toss2.read_text()
    keep_records1 = count_records(keep_text1)
    keep_records2 = count_records(keep_text2)
    toss_records1 = count_records(toss_text1)
    toss_records2 = count_records(toss_text2)

    if keep_records1 + toss_records1 != expected_records or keep_records2 + toss_records2 != expected_records:
        raise SystemExit(
            f'{case_name}/{threads} keep+toss wrote {keep_records1 + toss_records1}/'
            f'{keep_records2 + toss_records2} mate records, expected {expected_records}'
        )
    if toss_records1 == 0 or toss_records2 == 0:
        raise SystemExit(f'{case_name}/{threads} did not exercise paired multipass toss routing')

    kept_mutants = 0
    routed_mutants = 0
    for row in meta_rows:
        keep_mutant_count = count_any(keep_text1, row['mutants1']) + count_any(keep_text2, row['mutants2'])
        routed_mutant_count = (
            keep_mutant_count
            + count_any(toss_text1, row['mutants1'])
            + count_any(toss_text2, row['mutants2'])
        )
        kept_mutants += keep_mutant_count
        routed_mutants += routed_mutant_count
        if keep_mutant_count != 0:
            raise SystemExit(f'{case_name}/{threads} left mutant sequences in kept outputs for case {row["case"]}')

    summary.extend([
        f'{case_name}_{thread_label}_tossmode_keep1\t{keep_records1}',
        f'{case_name}_{thread_label}_tossmode_keep2\t{keep_records2}',
        f'{case_name}_{thread_label}_tossmode_toss1\t{toss_records1}',
        f'{case_name}_{thread_label}_tossmode_toss2\t{toss_records2}',
        f'{case_name}_{thread_label}_tossmode_kept_mutants\t{kept_mutants}',
        f'{case_name}_{thread_label}_tossmode_routed_mutants\t{routed_mutants}',
    ])

    key = case_name
    output_tuple = (keep_text1, keep_text2, toss_text1, toss_text2)
    if thread_label == base_label:
        base_outputs[key] = output_tuple
    elif key in base_outputs and base_outputs[key] != output_tuple:
        raise SystemExit(f'{case_name}/{threads} toss-routing output differs from base thread case {base_label}')

(out / 'summary.tsv').write_text('\n'.join(summary) + '\n')
PY

printf 'Real-derived multipass ECC toss-routing stress passed. Summary: %s\n' "$OUT/summary.tsv"

MARK_COMMON=(
  "in=$OUT/uncorrectable.1.fq"
  "in2=$OUT/uncorrectable.2.fq"
  "passes=2"
  "ecco=f"
  "eccmaxqual=0"
  "markuncorrectableerrors=t"
  "k=31"
  "minq=0"
  "minprob=0"
  "min=0"
  "minkmers=1"
  "target=999999999"
  "max=999999999"
  "keepall=t"
  "overwrite=t"
  "bits=32"
)

printf 'case\tthread_label\tthreads\tkeep1\tkeep2\tunc1\tunc2\n' > "$OUT/mark_run_matrix.tsv"

for thread_case in $THREAD_CASES; do
  thread_label="threads_$(label_for_threads "$thread_case")"

  for case_name in final_only first_only both_passes; do
    case "$case_name" in
      final_only) ECC_STAGE=("ecc1=f" "eccf=t") ;;
      first_only) ECC_STAGE=("ecc1=t" "eccf=f") ;;
      both_passes) ECC_STAGE=("ecc=t") ;;
      *) echo "unknown multipass ECC mark case: $case_name" >&2; exit 2 ;;
    esac

    printf 'Running Rust multipass ECC marked-uncorrectable case=%s threads=%s...\n' "$case_name" "$thread_case"
    target/debug/bbnorm-rs \
      "${MARK_COMMON[@]}" \
      "${ECC_STAGE[@]}" \
      "threads=$thread_case" \
      "out=$OUT/rust.$case_name.$thread_label.mark.keep1.fq" \
      "out2=$OUT/rust.$case_name.$thread_label.mark.keep2.fq" \
      "outuncorrected=$OUT/rust.$case_name.$thread_label.mark.unc1.fq" \
      "outuncorrected2=$OUT/rust.$case_name.$thread_label.mark.unc2.fq" \
      >"$OUT/rust.$case_name.$thread_label.mark.stdout.log" \
      2>"$OUT/rust.$case_name.$thread_label.mark.stderr.log"

    if grep -q 'error correction option ecc=t is not implemented yet' "$OUT/rust.$case_name.$thread_label.mark.stderr.log"; then
      echo "Rust emitted the old ECC fallback note for marked multipass $case_name/$thread_case" >&2
      exit 1
    fi

    printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
      "$case_name" "$thread_label" "$thread_case" \
      "$OUT/rust.$case_name.$thread_label.mark.keep1.fq" \
      "$OUT/rust.$case_name.$thread_label.mark.keep2.fq" \
      "$OUT/rust.$case_name.$thread_label.mark.unc1.fq" \
      "$OUT/rust.$case_name.$thread_label.mark.unc2.fq" >> "$OUT/mark_run_matrix.tsv"
  done
done

python - "$OUT" "$base_label" <<'PY'
import sys
from pathlib import Path

out = Path(sys.argv[1])
base_label = sys.argv[2]
meta_rows = []
for line in (out / 'fixture_meta.tsv').read_text().splitlines()[1:]:
    case, selected_record, pos1, pos2, clean1, clean2, mutants1, mutants2, uncorrectable1 = line.split('\t')
    meta_rows.append({
        'case': case,
        'clean1': clean1,
        'clean2': clean2,
        'uncorrectable1': uncorrectable1,
    })

runs = []
for line in (out / 'mark_run_matrix.tsv').read_text().splitlines()[1:]:
    case_name, thread_label, threads, keep1, keep2, unc1, unc2 = line.split('\t')
    runs.append((case_name, thread_label, threads, Path(keep1), Path(keep2), Path(unc1), Path(unc2)))

if not runs:
    raise SystemExit('No multipass ECC marked-uncorrectable stress runs were recorded')

def count_records(path):
    lines = Path(path).read_text().splitlines()
    if len(lines) % 4 != 0:
        raise SystemExit(f'{path} is not a complete FASTQ file')
    return len(lines) // 4

def fastq_map(path):
    lines = Path(path).read_text().splitlines()
    if len(lines) % 4 != 0:
        raise SystemExit(f'{path} is not a complete FASTQ file')
    return {
        lines[i][1:]: (lines[i + 1], lines[i + 3])
        for i in range(0, len(lines), 4)
    }

summary = (out / 'summary.tsv').read_text().splitlines()
expected_records = len(meta_rows) * 41
base_outputs = {}

for case_name, thread_label, threads, keep1, keep2, unc1, unc2 in runs:
    keep_records1 = count_records(keep1)
    keep_records2 = count_records(keep2)
    unc_records1 = count_records(unc1)
    unc_records2 = count_records(unc2)
    if keep_records1 != expected_records or keep_records2 != expected_records:
        raise SystemExit(
            f'{case_name}/{threads} marked keep outputs wrote {keep_records1}/{keep_records2}, expected {expected_records}'
        )
    expected_uncorrected = len(meta_rows) * (0 if case_name == 'first_only' else 1)
    if unc_records1 != expected_uncorrected or unc_records2 != expected_uncorrected:
        raise SystemExit(
            f'{case_name}/{threads} marked outuncorrected outputs wrote {unc_records1}/{unc_records2}, expected {expected_uncorrected}'
        )

    keep1_map = fastq_map(keep1)
    keep2_map = fastq_map(keep2)
    unc1_map = fastq_map(unc1)
    unc2_map = fastq_map(unc2)
    for row in meta_rows:
        r1_id = f'case{row["case"]}_uncorrectable/1'
        r2_id = f'case{row["case"]}_uncorrectable/2'
        if r1_id not in keep1_map:
            raise SystemExit(f'{case_name}/{threads} missing kept marked uncorrectable read {r1_id}')
        if r2_id not in keep2_map:
            raise SystemExit(f'{case_name}/{threads} missing kept clean mate read {r2_id}')
        keep1_seq, keep1_qual = keep1_map[r1_id]
        keep2_seq, keep2_qual = keep2_map[r2_id]
        if keep1_seq != row['uncorrectable1']:
            raise SystemExit(f'{case_name}/{threads} changed the forced uncorrectable keep sequence for {r1_id}')
        if keep2_seq != row['clean2']:
            raise SystemExit(f'{case_name}/{threads} changed the clean mate keep sequence for {r2_id}')
        if keep1_qual == 'I' * len(keep1_qual):
            raise SystemExit(f'{case_name}/{threads} left marked qualities unchanged for {r1_id}')
        if keep2_qual != 'I' * len(keep2_qual):
            raise SystemExit(f'{case_name}/{threads} unexpectedly changed clean mate qualities for {r2_id}')
        if expected_uncorrected:
            if r1_id not in unc1_map:
                raise SystemExit(f'{case_name}/{threads} missing uncorrected marked read {r1_id}')
            if r2_id not in unc2_map:
                raise SystemExit(f'{case_name}/{threads} missing uncorrected clean mate {r2_id}')
            unc1_seq, unc1_qual = unc1_map[r1_id]
            unc2_seq, unc2_qual = unc2_map[r2_id]
            if unc1_seq != row['uncorrectable1']:
                raise SystemExit(f'{case_name}/{threads} changed the forced uncorrectable side-output sequence for {r1_id}')
            if unc2_seq != row['clean2']:
                raise SystemExit(f'{case_name}/{threads} changed the clean mate side-output sequence for {r2_id}')
            if unc1_qual != keep1_qual:
                raise SystemExit(f'{case_name}/{threads} keep/uncorrected quality mismatch for {r1_id}')
            if unc2_qual != 'I' * len(unc2_qual):
                raise SystemExit(f'{case_name}/{threads} unexpectedly changed clean mate side-output qualities for {r2_id}')
        else:
            if r1_id in unc1_map or r2_id in unc2_map:
                raise SystemExit(
                    f'{case_name}/{threads} unexpectedly emitted staged outuncorrected reads before the final ECC pass'
                )

    summary.extend([
        f'{case_name}_{thread_label}_mark_keep1\t{keep_records1}',
        f'{case_name}_{thread_label}_mark_keep2\t{keep_records2}',
        f'{case_name}_{thread_label}_mark_unc1\t{unc_records1}',
        f'{case_name}_{thread_label}_mark_unc2\t{unc_records2}',
    ])

    key = case_name
    output_tuple = (keep1.read_bytes(), keep2.read_bytes(), unc1.read_bytes(), unc2.read_bytes())
    if thread_label == base_label:
        base_outputs[key] = output_tuple
    elif key in base_outputs and base_outputs[key] != output_tuple:
        raise SystemExit(f'{case_name}/{threads} marked output differs from base thread case {base_label}')

(out / 'summary.tsv').write_text('\n'.join(summary) + '\n')
PY

printf 'Real-derived multipass ECC marked-uncorrectable stress passed. Summary: %s\n' "$OUT/summary.tsv"

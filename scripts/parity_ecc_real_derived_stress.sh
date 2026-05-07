#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT"

DATA1="${DATA1:-vendor/BBTools-master/resources/sample1.fq.gz}"
DATA2="${DATA2:-vendor/BBTools-master/resources/sample2.fq.gz}"
OUT="${1:-tmp/ecc_real_derived_stress}"
STRESS_RECORDS="${STRESS_RECORDS:-1}"
THREADS="${THREADS:-1}"
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

cargo build --quiet

printf 'Building a no-N paired real-derived ECC stress fixture from %s and %s...\n' "$DATA1" "$DATA2"
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
        f'Could not find {wanted} no-N paired reads at least {min_read_len} bp long for the ECC fixture'
    )

meta_lines = ['case\tselected_record\tpos1\tpos2\tclean1\tclean2\tmutants1\tmutants2\tnoisy1\tnoisy2']
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
        noisy1 = mutate(mut1[0], pos1 + 1)
        noisy2 = mutate(mut2[0], pos2 + 1)
        for i in range(40):
            r1_out.write(f'@case{case}_clean{i + 1:03d}/1\n{seq1}\n+\n{qual1}\n')
            r2_out.write(f'@case{case}_clean{i + 1:03d}/2\n{seq2}\n+\n{qual2}\n')
        for i in range(30):
            unc1_out.write(f'@case{case}_clean{i + 1:03d}/1\n{seq1}\n+\n{qual1}\n')
            unc2_out.write(f'@case{case}_clean{i + 1:03d}/2\n{seq2}\n+\n{qual2}\n')
        for i, seq in enumerate(mut1, start=1):
            r1_out.write(f'@case{case}_mut_r1_{i}/1\n{seq}\n+\n{qual1}\n')
            r2_out.write(f'@case{case}_mut_r1_{i}/2\n{seq2}\n+\n{qual2}\n')
        for i, seq in enumerate(mut2, start=1):
            r1_out.write(f'@case{case}_mut_r2_{i}/1\n{seq1}\n+\n{qual1}\n')
            r2_out.write(f'@case{case}_mut_r2_{i}/2\n{seq}\n+\n{qual2}\n')
        unc1_out.write(f'@case{case}_noisy/1\n{noisy1}\n+\n{qual1}\n')
        unc2_out.write(f'@case{case}_noisy/2\n{noisy2}\n+\n{qual2}\n')
        meta_lines.append(
            f'{case}\t{idx}\t{pos1}\t{pos2}\t{seq1}\t{seq2}\t{",".join(mut1)}\t{",".join(mut2)}\t{noisy1}\t{noisy2}'
        )

(out / 'fixture_meta.tsv').write_text('\n'.join(meta_lines) + '\n')
PY

COMMON=(
  "in=$OUT/stress.1.fq"
  "in2=$OUT/stress.2.fq"
  "passes=1"
  "ecc=t"
  "ecco=f"
  "eccmaxqual=99"
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

printf 'Running Rust keep-all ECC on the real-derived stress fixture with threads=%s...\n' "$THREADS"
target/debug/bbnorm-rs \
  "${COMMON[@]}" \
  "keepall=t" \
  "threads=$THREADS" \
  "out=$OUT/rust.keepall.1.fq" "out2=$OUT/rust.keepall.2.fq" \
  >"$OUT/rust.keepall.stdout.log" 2>"$OUT/rust.keepall.stderr.log"

if grep -q 'error correction option ecc=t is not implemented yet' "$OUT/rust.keepall.stderr.log"; then
  echo "Rust emitted the old ECC fallback note" >&2
  exit 1
fi

printf 'Running Rust normalization/toss ECC on the real-derived stress fixture with threads=%s...\n' "$THREADS"
target/debug/bbnorm-rs \
  "${COMMON[@]}" \
  "target=2" "max=2" "tossbadreads=t" \
  "threads=$THREADS" \
  "out=$OUT/rust.norm.1.fq" "out2=$OUT/rust.norm.2.fq" \
  "outt=$OUT/rust.norm.toss1.fq" "outt2=$OUT/rust.norm.toss2.fq" \
  >"$OUT/rust.norm.stdout.log" 2>"$OUT/rust.norm.stderr.log"

printf 'Running Rust keep-all ECC outuncorrected stress on the noisier real-derived fixture with threads=%s...\n' "$THREADS"
target/debug/bbnorm-rs \
  "in=$OUT/uncorrectable.1.fq" \
  "in2=$OUT/uncorrectable.2.fq" \
  "passes=1" \
  "keepall=t" \
  "ecc=t" \
  "ecco=f" \
  "eccmaxqual=0" \
  "k=31" \
  "minq=0" \
  "minprob=0" \
  "min=0" \
  "minkmers=1" \
  "target=999999999" \
  "max=999999999" \
  "threads=$THREADS" \
  "overwrite=t" \
  "bits=32" \
  "out=$OUT/rust.uncorrectable.keep1.fq" \
  "out2=$OUT/rust.uncorrectable.keep2.fq" \
  "outuncorrected=$OUT/rust.uncorrectable.unc1.fq" \
  "outuncorrected2=$OUT/rust.uncorrectable.unc2.fq" \
  >"$OUT/rust.uncorrectable.stdout.log" 2>"$OUT/rust.uncorrectable.stderr.log"

python - "$OUT" <<'PY'
import sys
from pathlib import Path

out = Path(sys.argv[1])
meta_rows = []
for line in (out / 'fixture_meta.tsv').read_text().splitlines()[1:]:
    case, selected_record, pos1, pos2, clean1, clean2, mutants1, mutants2, noisy1, noisy2 = line.split('\t')
    meta_rows.append({
        'case': case,
        'selected_record': selected_record,
        'clean1': clean1,
        'clean2': clean2,
        'mutants1': mutants1.split(','),
        'mutants2': mutants2.split(','),
        'noisy1': noisy1,
        'noisy2': noisy2,
    })

input_r1 = (out / 'stress.1.fq').read_text()
input_r2 = (out / 'stress.2.fq').read_text()
keep1 = (out / 'rust.keepall.1.fq').read_text()
keep2 = (out / 'rust.keepall.2.fq').read_text()
norm1 = (out / 'rust.norm.1.fq').read_text()
norm2 = (out / 'rust.norm.2.fq').read_text()
toss1 = (out / 'rust.norm.toss1.fq').read_text()
toss2 = (out / 'rust.norm.toss2.fq').read_text()
unc_in1 = (out / 'uncorrectable.1.fq').read_text()
unc_in2 = (out / 'uncorrectable.2.fq').read_text()
unc_keep1 = (out / 'rust.uncorrectable.keep1.fq').read_text()
unc_keep2 = (out / 'rust.uncorrectable.keep2.fq').read_text()
unc_out1 = (out / 'rust.uncorrectable.unc1.fq').read_text()
unc_out2 = (out / 'rust.uncorrectable.unc2.fq').read_text()

def count_any(text, seqs):
    return sum(text.count(seq) for seq in seqs)

summary = ['metric\tvalue']
total_input_mutants = 0
total_keepall_mutants = 0
total_norm_keep_mutants = 0
for row in meta_rows:
    case = row['case']
    clean1 = row['clean1']
    clean2 = row['clean2']
    mutants1 = row['mutants1']
    mutants2 = row['mutants2']
    noisy1 = row['noisy1']
    noisy2 = row['noisy2']
    input_m1 = count_any(input_r1, mutants1)
    input_m2 = count_any(input_r2, mutants2)
    keep_m1 = count_any(keep1, mutants1)
    keep_m2 = count_any(keep2, mutants2)
    norm_m1 = count_any(norm1, mutants1)
    norm_m2 = count_any(norm2, mutants2)
    unc_input_noisy1 = unc_in1.count(noisy1)
    unc_input_noisy2 = unc_in2.count(noisy2)
    unc_keep_noisy1 = unc_keep1.count(noisy1)
    unc_keep_noisy2 = unc_keep2.count(noisy2)
    unc_out_noisy1 = unc_out1.count(noisy1)
    unc_out_noisy2 = unc_out2.count(noisy2)
    clean_count1 = keep1.count(clean1)
    clean_count2 = keep2.count(clean2)
    total_input_mutants += input_m1 + input_m2
    total_keepall_mutants += keep_m1 + keep_m2
    total_norm_keep_mutants += norm_m1 + norm_m2
    if input_m1 != 3 or input_m2 != 3:
        raise SystemExit(f'Expected exactly three mutant records per mate in input fixture case {case}')
    if keep_m1 != 0 or keep_m2 != 0:
        raise SystemExit(f'Rust keep-all ECC left real-derived mutant sequences in case {case}')
    if clean_count1 != 46 or clean_count2 != 46:
        raise SystemExit(f'Rust keep-all ECC did not restore all reads to clean sequences in case {case}')
    if norm_m1 != 0 or norm_m2 != 0:
        raise SystemExit(f'Rust normalization keep output retained real-derived mutant sequences in case {case}')
    if unc_input_noisy1 != 1 or unc_input_noisy2 != 1:
        raise SystemExit(f'Expected exactly one noisier multi-mutation input pair in case {case}')
    if unc_keep_noisy1 != 1 or unc_keep_noisy2 != 1:
        raise SystemExit(f'Rust keep-all uncorrectable stress did not retain the noisy pair in keep output for case {case}')
    if unc_out_noisy1 != 1 or unc_out_noisy2 != 1:
        raise SystemExit(f'Rust keep-all uncorrectable stress did not route the noisy pair to outuncorrected for case {case}')
    summary.extend([
        f'case{case}_input_mutants\t{input_m1 + input_m2}',
        f'case{case}_keepall_mutants\t{keep_m1 + keep_m2}',
        f'case{case}_keepall_clean1\t{clean_count1}',
        f'case{case}_keepall_clean2\t{clean_count2}',
        f'case{case}_norm_keep_mutants\t{norm_m1 + norm_m2}',
        f'case{case}_uncorrectable_keep_noisy\t{unc_keep_noisy1 + unc_keep_noisy2}',
        f'case{case}_uncorrectable_out_noisy\t{unc_out_noisy1 + unc_out_noisy2}',
    ])

if not toss1 or not toss2:
    raise SystemExit('Rust normalization stress run did not exercise paired toss routing')

if not unc_out1 or not unc_out2:
    raise SystemExit('Rust noisier uncorrectable stress did not exercise paired outuncorrected routing')

summary.extend([
    f'total_cases\t{len(meta_rows)}',
    f'total_input_mutants\t{total_input_mutants}',
    f'total_keepall_mutants\t{total_keepall_mutants}',
    f'total_norm_keep_mutants\t{total_norm_keep_mutants}',
    f'norm_toss1_records\t{len(toss1.splitlines()) // 4}',
    f'norm_toss2_records\t{len(toss2.splitlines()) // 4}',
    f'uncorrectable_out1_records\t{len(unc_out1.splitlines()) // 4}',
    f'uncorrectable_out2_records\t{len(unc_out2.splitlines()) // 4}',
])
(out / 'summary.tsv').write_text('\n'.join(summary) + '\n')
PY

printf 'Real-derived ECC stress passed. Summary: %s\n' "$OUT/summary.tsv"

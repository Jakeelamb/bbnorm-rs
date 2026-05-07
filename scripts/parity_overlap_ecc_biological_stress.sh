#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT"

BIO_ROOT="${BIO_ROOT:-/home/jake/Projects/biological data}"
DATA1="${DATA1:-$BIO_ROOT/reads/short_reads/scer_s288c_pe_srr23631023/SRR23631023_1.fastq.gz}"
DATA2="${DATA2:-$BIO_ROOT/reads/short_reads/scer_s288c_pe_srr23631023/SRR23631023_2.fastq.gz}"
OUT="${1:-tmp/overlap_ecc_biological_stress}"
THREAD_CASES="${THREAD_CASES:-${THREADS:-1 2 auto}}"
OVERLAP_PAIRS="${OVERLAP_PAIRS:-5}"
OVERLAP_LEN="${OVERLAP_LEN:-40}"
PYTHON="${PYTHON:-python}"
mkdir -p "$OUT"
rm -f "$OUT"/*
ACCEPT_PREFIX="accept"
REJECT_PREFIX="reject"

if [[ ! -f "$DATA1" || ! -f "$DATA2" ]]; then
  printf 'Missing paired dataset files:\n  DATA1=%s\n  DATA2=%s\n' "$DATA1" "$DATA2" >&2
  exit 2
fi
for value_name in OVERLAP_PAIRS OVERLAP_LEN; do
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

printf 'Building real-derived overlap-only ECC fixture from %s and %s...\n' "$DATA1" "$DATA2"
"$PYTHON" - "$DATA1" "$DATA2" "$OUT" "$OVERLAP_PAIRS" "$OVERLAP_LEN" <<'PY'
import gzip
import sys
from pathlib import Path

data1, data2, out, pairs, overlap_len = (
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

def revcomp(seq):
    table = str.maketrans('ACGT', 'TGCA')
    return seq.translate(table)[::-1]

def mutate(seq, pos):
    for base in 'ACGT':
        if base != seq[pos]:
            return seq[:pos] + base + seq[pos + 1:]
    raise AssertionError('unreachable')

selected = None
for idx, (r1, r2) in enumerate(zip(records(data1), records(data2)), start=1):
    _, seq1, _, _ = r1
    _, seq2, _, _ = r2
    if len(seq1) >= overlap_len and len(seq2) >= overlap_len and set(seq1) <= set('ACGT') and set(seq2) <= set('ACGT'):
        selected = (idx, seq1[:overlap_len])
        break
if selected is None:
    raise SystemExit(f'Could not find a no-N paired read long enough for overlap length {overlap_len}')

idx, read1 = selected
read2_clean = revcomp(read1)
mut_pos = min(max(1, overlap_len // 2), overlap_len - 2)
read2_mutant = mutate(read2_clean, mut_pos)
qual1 = 'I' * overlap_len
qual2 = ''.join('!' if i == mut_pos else 'I' for i in range(overlap_len))

with (out / 'accept.1.fq').open('w') as r1_out, (out / 'accept.2.fq').open('w') as r2_out:
    for i in range(1, pairs + 1):
        r1_out.write(f'@accept{i}/1\n{read1}\n+\n{qual1}\n')
        r2_out.write(f'@accept{i}/2\n{read2_mutant}\n+\n{qual2}\n')

(out / 'accept_fixture_meta.tsv').write_text(
    '\n'.join([
        'selected_record\tread1\tread2_clean\tread2_mutant\tqual1\tqual2',
        f'{idx}\t{read1}\t{read2_clean}\t{read2_mutant}\t{qual1}\t{qual2}',
    ]) + '\n'
)
PY

cat > "$OUT/$REJECT_PREFIX.1.fq" <<'EOF'
@reject1/1
CAGTAACCAATGCCTGTTGAGATGCCAGACGCGTAACCAAAA
+
IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
@reject2/1
CAGTAACCAATGCCTGTTGAGATGCCAGACGCGTAACCAAAA
+
IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
@reject3/1
CAGTAACCAATGCCTGTTGAGATGCCAGACGCGTAACCAAAA
+
IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
@reject4/1
CAGTAACCAATGCCTGTTGAGATGCCAGACGCGTAACCAAAA
+
IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
@reject5/1
CAGTAACCAATGCCTGTTGAGATGCCAGACGCGTAACCAAAA
+
IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
EOF
cat > "$OUT/$REJECT_PREFIX.2.fq" <<'EOF'
@reject1/2
TTTTGCTAACGCGTCTGGCATCTCAACAGGCATTGGTTAC
+
IIIII!I'IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
@reject2/2
TTTTGCTAACGCGTCTGGCATCTCAACAGGCATTGGTTAC
+
IIIII!I'IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
@reject3/2
TTTTGCTAACGCGTCTGGCATCTCAACAGGCATTGGTTAC
+
IIIII!I'IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
@reject4/2
TTTTGCTAACGCGTCTGGCATCTCAACAGGCATTGGTTAC
+
IIIII!I'IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
@reject5/2
TTTTGCTAACGCGTCTGGCATCTCAACAGGCATTGGTTAC
+
IIIII!I'IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
EOF
cat > "$OUT/reject_fixture_meta.tsv" <<'EOF'
read1	read2_mutant	qual1	qual2
CAGTAACCAATGCCTGTTGAGATGCCAGACGCGTAACCAAAA	TTTTGCTAACGCGTCTGGCATCTCAACAGGCATTGGTTAC	IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII	IIIII!I'IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
EOF

run_single_case() {
  local threads="$1"
  local label="$2"

  for fixture in "$ACCEPT_PREFIX" "$REJECT_PREFIX"; do
    for ecco in f t; do
      local prefix="$OUT/$label.$fixture.single.ecco_$ecco"
      printf 'Running overlap-only ECC biological stress fixture=%s threads=%s ecco=%s...\n' "$fixture" "$threads" "$ecco"
      target/debug/bbnorm-rs \
      "in=$OUT/$fixture.1.fq" \
      "in2=$OUT/$fixture.2.fq" \
      "passes=1" \
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
      "threads=$threads" \
      "overwrite=t" \
      "bits=32" \
        "out=$prefix.1.fq" \
        "out2=$prefix.2.fq" \
        >"$prefix.stdout.log" 2>"$prefix.stderr.log"
    done
  done
}

run_countup_case() {
  local threads="$1"
  local label="$2"

  for fixture in "$ACCEPT_PREFIX" "$REJECT_PREFIX"; do
    for ecco in f t; do
      local prefix="$OUT/$label.$fixture.countup.ecco_$ecco"
      printf 'Running overlap-only countup ECC biological stress fixture=%s threads=%s ecco=%s...\n' "$fixture" "$threads" "$ecco"
      target/debug/bbnorm-rs \
      "in=$OUT/$fixture.1.fq" \
      "in2=$OUT/$fixture.2.fq" \
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
      "threads=$threads" \
      "overwrite=t" \
      "bits=32" \
        "out=$prefix.1.fq" \
        "out2=$prefix.2.fq" \
        >"$prefix.stdout.log" 2>"$prefix.stderr.log"
    done
  done
}

labels=()
for threads in $THREAD_CASES; do
  label="$(label_for_threads "$threads")"
  labels+=("$label")
  run_single_case "$threads" "$label"
  run_countup_case "$threads" "$label"
done

"$PYTHON" - "$OUT" "$OVERLAP_PAIRS" "${labels[@]}" <<'PY'
import sys
from pathlib import Path

out = Path(sys.argv[1])
expected_pairs = int(sys.argv[2])
labels = sys.argv[3:]
if not labels:
    raise SystemExit('No thread labels were supplied')

accept_meta = (out / 'accept_fixture_meta.tsv').read_text().splitlines()[1].split('\t')
_, accept_read1, accept_read2_clean, accept_read2_mutant, accept_qual1, _ = accept_meta
reject_meta = (out / 'reject_fixture_meta.tsv').read_text().splitlines()[1].split('\t')
reject_read1, reject_read2_mutant, reject_qual1, _ = reject_meta

def records(path):
    lines = Path(path).read_text().splitlines()
    if len(lines) % 4:
        raise SystemExit(f'{path} is not a complete FASTQ file')
    return len(lines) // 4

def read_fastq(path):
    lines = Path(path).read_text().splitlines()
    if len(lines) % 4:
        raise SystemExit(f'{path} is not a complete FASTQ file')
    return [tuple(lines[i:i + 4]) for i in range(0, len(lines), 4)]

def output_path(label, fixture, mode, ecco, mate):
    return out / f'{label}.{fixture}.{mode}.ecco_{ecco}.{mate}.fq'

summary = ['metric\tvalue']
compare_suffixes = [
    'accept.single.ecco_f.1.fq', 'accept.single.ecco_f.2.fq',
    'accept.single.ecco_t.1.fq', 'accept.single.ecco_t.2.fq',
    'accept.countup.ecco_f.1.fq', 'accept.countup.ecco_f.2.fq',
    'accept.countup.ecco_t.1.fq', 'accept.countup.ecco_t.2.fq',
    'reject.single.ecco_f.1.fq', 'reject.single.ecco_f.2.fq',
    'reject.single.ecco_t.1.fq', 'reject.single.ecco_t.2.fq',
    'reject.countup.ecco_f.1.fq', 'reject.countup.ecco_f.2.fq',
    'reject.countup.ecco_t.1.fq', 'reject.countup.ecco_t.2.fq',
]

for label in labels:
    for mode in ['single', 'countup']:
        accept_unchanged = read_fastq(output_path(label, 'accept', mode, 'f', 1))
        accept_mutant = read_fastq(output_path(label, 'accept', mode, 'f', 2))
        accept_repaired = read_fastq(output_path(label, 'accept', mode, 't', 2))
        reject_unchanged = read_fastq(output_path(label, 'reject', mode, 'f', 1))
        reject_mutant = read_fastq(output_path(label, 'reject', mode, 'f', 2))
        reject_rejected = read_fastq(output_path(label, 'reject', mode, 't', 2))
        if records(output_path(label, 'accept', mode, 'f', 1)) != expected_pairs:
            raise SystemExit(f'{label} accept {mode} ecco=f mate1 wrote the wrong record count')
        if records(output_path(label, 'accept', mode, 'f', 2)) != expected_pairs:
            raise SystemExit(f'{label} accept {mode} ecco=f mate2 wrote the wrong record count')
        if records(output_path(label, 'accept', mode, 't', 2)) != expected_pairs:
            raise SystemExit(f'{label} accept {mode} ecco=t mate2 wrote the wrong record count')
        if records(output_path(label, 'reject', mode, 'f', 1)) != expected_pairs:
            raise SystemExit(f'{label} reject {mode} ecco=f mate1 wrote the wrong record count')
        if records(output_path(label, 'reject', mode, 'f', 2)) != expected_pairs:
            raise SystemExit(f'{label} reject {mode} ecco=f mate2 wrote the wrong record count')
        if records(output_path(label, 'reject', mode, 't', 2)) != expected_pairs:
            raise SystemExit(f'{label} reject {mode} ecco=t mate2 wrote the wrong record count')
        if any(seq != accept_read1 or q != accept_qual1 for _, seq, _, q in accept_unchanged):
            raise SystemExit(f'{label} accept {mode} changed mate1 unexpectedly')
        if any(seq != accept_read2_mutant for _, seq, _, _ in accept_mutant):
            raise SystemExit(f'{label} accept {mode} ecco=f changed the lower-quality mate unexpectedly')
        if any(seq != accept_read2_clean for _, seq, _, _ in accept_repaired):
            raise SystemExit(f'{label} accept {mode} ecco=t failed to repair the lower-quality mate')
        if any(seq != reject_read1 or q != reject_qual1 for _, seq, _, q in reject_unchanged):
            raise SystemExit(f'{label} reject {mode} changed mate1 unexpectedly')
        if any(seq != reject_read2_mutant for _, seq, _, _ in reject_mutant):
            raise SystemExit(f'{label} reject {mode} ecco=f changed the competing-overlap mate unexpectedly')
        if any(seq != reject_read2_mutant for _, seq, _, _ in reject_rejected):
            raise SystemExit(f'{label} reject {mode} ecco=t failed to reject the competing-overlap ambiguity')
        if Path(out / f'{label}.accept.{mode}.ecco_t.stderr.log').read_text().find('paired overlap repair before the table-based ECC path') == -1:
            raise SystemExit(f'{label} accept {mode} missing overlap-repair stderr note for ecco=t')
        if Path(out / f'{label}.reject.{mode}.ecco_t.stderr.log').read_text().find('paired overlap repair before the table-based ECC path') == -1:
            raise SystemExit(f'{label} reject {mode} missing overlap-repair stderr note for ecco=t')
        summary.extend([
            f'{label}_{mode}_pairs\t{expected_pairs}',
            f'{label}_{mode}_accept_mutant_seq\t{accept_read2_mutant}',
            f'{label}_{mode}_accept_clean_seq\t{accept_read2_clean}',
            f'{label}_{mode}_reject_seq\t{reject_read2_mutant}',
        ])

baseline = labels[0]
for label in labels[1:]:
    for suffix in compare_suffixes:
        base = out / f'{baseline}.{suffix}'
        other = out / f'{label}.{suffix}'
        if base.read_bytes() != other.read_bytes():
            raise SystemExit(f'Overlap ECC output changed across threads for {suffix}: {baseline} vs {label}')

summary.append(f'thread_cases\t{" ".join(labels)}')
summary.append(f'thread_compare_baseline\t{baseline}')
(out / 'summary.tsv').write_text('\n'.join(summary) + '\n')
PY

printf 'Overlap ECC biological stress passed. Summary: %s\n' "$OUT/summary.tsv"

#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT"

OUT="${1:-tmp/biological_working_modes_truncated_parity}"
BIO_ROOT="${BIO_ROOT:-/home/jake/Projects/biological data}"
DATA1="${DATA1:-$BIO_ROOT/reads/short_reads/scer_s288c_pe_srr23631023/SRR23631023_1.fastq.gz}"
DATA2="${DATA2:-$BIO_ROOT/reads/short_reads/scer_s288c_pe_srr23631023/SRR23631023_2.fastq.gz}"
READS="${READS:-1000}"
TABLE_READS="${TABLE_READS:-$READS}"

rm -rf "$OUT"
mkdir -p "$OUT"
printf 'mode\tstatus\toutput\n' > "$OUT/summary.tsv"

run_mode() {
  local mode_name="$1"
  local extra_args="$2"
  local mode_out="$OUT/$mode_name"

  printf '\nRunning truncated biological Java parity for %s...\n' "$mode_name"
  BIO_ROOT="$BIO_ROOT" \
    DATA1="$DATA1" \
    DATA2="$DATA2" \
    READS="$READS" \
    TABLE_READS="$TABLE_READS" \
    EXTRA_ARGS="$extra_args" \
    scripts/parity_biological_dataset_truncated.sh "$mode_out"
  printf '%s\tpassed\t%s\n' "$mode_name" "$mode_out" >> "$OUT/summary.tsv"
}

probe_mode_divergence() {
  local mode_name="$1"
  local extra_args="$2"
  local mode_out="$OUT/$mode_name"

  printf '\nRunning truncated biological Java parity probe for %s...\n' "$mode_name"
  set +e
  BIO_ROOT="$BIO_ROOT" \
    DATA1="$DATA1" \
    DATA2="$DATA2" \
    READS="$READS" \
    TABLE_READS="$TABLE_READS" \
    EXTRA_ARGS="$extra_args" \
    scripts/parity_biological_dataset_truncated.sh "$mode_out" \
    >"$mode_out.probe.stdout.log" 2>"$mode_out.probe.stderr.log"
  status=$?
  set -e

  if [[ "$status" -eq 0 ]]; then
    printf '%s\tmatched\t%s\n' "$mode_name" "$mode_out" >> "$OUT/summary.tsv"
    return 0
  fi

  python - "$mode_out" <<'PY'
from pathlib import Path
import sys

base = Path(sys.argv[1])

def record_count(path: Path) -> int:
    if not path.exists():
        return 0
    return len(path.read_text().splitlines()) // 4

summary = [
    ("java_mid_pairs", record_count(base / "java.mid1.fq")),
    ("rust_mid_pairs", record_count(base / "rust.mid1.fq")),
    ("java_high_pairs", record_count(base / "java.high1.fq")),
    ("rust_high_pairs", record_count(base / "rust.high1.fq")),
]

def record_ids(path: Path):
    if not path.exists():
        return []
    rows = []
    with path.open(encoding="utf-8") as handle:
        while True:
            header = handle.readline().rstrip("\n")
            if not header:
                break
            handle.readline()
            handle.readline()
            handle.readline()
            rows.append(header[1:].split()[0])
    return rows

java_high = set(record_ids(base / "java.high1.fq"))
rust_high = set(record_ids(base / "rust.high1.fq"))
shifted_to_mid = sorted(java_high - rust_high)

(base / "divergence.tsv").write_text(
    "metric\tvalue\n"
    + "\n".join(f"{key}\t{value}" for key, value in summary)
    + ("\nshifted_high_to_mid_ids\t" + ",".join(shifted_to_mid) if shifted_to_mid else "\nshifted_high_to_mid_ids\tnone")
    + "\n",
    encoding="utf-8",
)
PY

  printf '%s\tdivergence_confirmed\t%s\n' "$mode_name" "$mode_out" >> "$OUT/summary.tsv"
}

run_mode default ""
probe_mode_divergence long_kmer "k=40"
probe_mode_divergence long_kmer_fixspikes "k=40 fixspikes=t"

printf '\nTruncated biological working-mode parity summary:\n'
column -t -s $'\t' "$OUT/summary.tsv" || cat "$OUT/summary.tsv"
printf '\nTruncated biological working-mode parity probe completed. Outputs and logs: %s\n' "$OUT"

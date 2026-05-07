#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT"

BIO_ROOT="${BIO_ROOT:-/home/jake/Projects/biological data}"
DATA1="${DATA1:-$BIO_ROOT/reads/short_reads/scer_s288c_pe_srr23631023/SRR23631023_1.fastq.gz}"
DATA2="${DATA2:-$BIO_ROOT/reads/short_reads/scer_s288c_pe_srr23631023/SRR23631023_2.fastq.gz}"
CP="vendor/BBTools-master/current"
OUT="${1:-tmp/biological_dataset_truncated_parity}"
READS="${READS:-1000}"
TABLE_READS="${TABLE_READS:-$READS}"
EXTRA_ARGS="${EXTRA_ARGS:-}"
PYTHON="${PYTHON:-python}"

rm -rf "$OUT"
mkdir -p "$OUT"

if [[ ! -f "$DATA1" || ! -f "$DATA2" ]]; then
  printf 'Missing paired dataset files:\n  DATA1=%s\n  DATA2=%s\n' "$DATA1" "$DATA2" >&2
  exit 2
fi

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
  "threads=1"
  "overwrite=t"
  "bits=32"
)

EXTRA_ARGS_ARRAY=()
if [[ -n "$EXTRA_ARGS" ]]; then
  read -r -a EXTRA_ARGS_ARRAY <<< "$EXTRA_ARGS"
fi

measure() {
  local metrics_file="$1"
  shift
  "$PYTHON" - "$metrics_file" "$@" <<'PY'
import resource
import subprocess
import sys
import time

metrics_file = sys.argv[1]
cmd = sys.argv[2:]
start = time.perf_counter()
proc = subprocess.run(cmd)
elapsed = time.perf_counter() - start
usage = resource.getrusage(resource.RUSAGE_CHILDREN)
with open(metrics_file, "w", encoding="utf-8") as handle:
    handle.write(f"{elapsed:.6f}\t{usage.ru_maxrss}\t{proc.returncode}\n")
sys.exit(proc.returncode)
PY
}

write_result() {
  local tool="$1"
  local metrics_file="$2"
  local elapsed rss status
  read -r elapsed rss status < "$metrics_file"
  printf '%s\t%s\t%s\t%s\n' "$tool" "$elapsed" "$rss" "$status" >> "$OUT/results.tsv"
}

cargo build --release --quiet

printf 'dataset1\t%s\n' "$DATA1" > "$OUT/dataset.tsv"
printf 'dataset2\t%s\n' "$DATA2" >> "$OUT/dataset.tsv"
printf 'reads\t%s\n' "$READS" >> "$OUT/dataset.tsv"
printf 'tablereads\t%s\n' "$TABLE_READS" >> "$OUT/dataset.tsv"
printf 'extra_args\t%s\n' "${EXTRA_ARGS:-none}" >> "$OUT/dataset.tsv"
printf 'tool\telapsed_seconds\tmax_rss_kb\tstatus\n' > "$OUT/results.tsv"

printf 'Running Java BBNorm on truncated biological dataset...\n'
measure "$OUT/java.metrics.tsv" \
  java -Xmx1g -cp "$CP" jgi.KmerNormalize \
  "${COMMON[@]}" \
  "${EXTRA_ARGS_ARRAY[@]}" \
  "out=$OUT/java.keep1.fq" "out2=$OUT/java.keep2.fq" \
  "outlow=$OUT/java.low1.fq" "outlow2=$OUT/java.low2.fq" \
  "outmid=$OUT/java.mid1.fq" "outmid2=$OUT/java.mid2.fq" \
  "outhigh=$OUT/java.high1.fq" "outhigh2=$OUT/java.high2.fq" \
  "hist=$OUT/java.hist.tsv" "rhist=$OUT/java.rhist.tsv" \
  >"$OUT/java.stdout.log" 2>"$OUT/java.stderr.log"
write_result java "$OUT/java.metrics.tsv"

printf 'Running Rust bbnorm-rs on same truncated biological dataset...\n'
measure "$OUT/rust.metrics.tsv" \
  target/release/bbnorm-rs \
  "${COMMON[@]}" \
  "${EXTRA_ARGS_ARRAY[@]}" \
  "out=$OUT/rust.keep1.fq" "out2=$OUT/rust.keep2.fq" \
  "outlow=$OUT/rust.low1.fq" "outlow2=$OUT/rust.low2.fq" \
  "outmid=$OUT/rust.mid1.fq" "outmid2=$OUT/rust.mid2.fq" \
  "outhigh=$OUT/rust.high1.fq" "outhigh2=$OUT/rust.high2.fq" \
  "hist=$OUT/rust.hist.tsv" "rhist=$OUT/rust.rhist.tsv" \
  >"$OUT/rust.stdout.log" 2>"$OUT/rust.stderr.log"
write_result rust "$OUT/rust.metrics.tsv"

for suffix in keep1.fq keep2.fq low1.fq low2.fq mid1.fq mid2.fq high1.fq high2.fq hist.tsv rhist.tsv; do
  cmp "$OUT/java.$suffix" "$OUT/rust.$suffix"
done

"$PYTHON" - "$OUT/results.tsv" <<'PY'
import csv
import sys

rows = list(csv.DictReader(open(sys.argv[1], encoding="utf-8"), delimiter="\t"))
print("\nTruncated biological parity summary:")
print("tool\telapsed_seconds\tmax_rss_kb")
for row in rows:
    print(f"{row['tool']}\t{float(row['elapsed_seconds']):.6f}\t{row['max_rss_kb']}")
PY

printf '\nTruncated biological parity passed. Results: %s\n' "$OUT/results.tsv"
printf 'Dataset metadata: %s\n' "$OUT/dataset.tsv"
printf 'Outputs and logs: %s\n' "$OUT"

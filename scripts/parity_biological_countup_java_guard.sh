#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT"

BIO_ROOT="${BIO_ROOT:-/home/jake/Projects/biological data}"
DATA1="${DATA1:-$BIO_ROOT/reads/short_reads/scer_s288c_pe_srr23631023/SRR23631023_1.fastq.gz}"
DATA2="${DATA2:-$BIO_ROOT/reads/short_reads/scer_s288c_pe_srr23631023/SRR23631023_2.fastq.gz}"
CP="vendor/BBTools-master/current"
OUT="${1:-tmp/biological_countup_java_guard}"
READS="${READS:-1000}"
TABLE_READS="${TABLE_READS:-$READS}"
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
  "countup=t"
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

cargo build --release --quiet

printf 'dataset1\t%s\n' "$DATA1" > "$OUT/dataset.tsv"
printf 'dataset2\t%s\n' "$DATA2" >> "$OUT/dataset.tsv"
printf 'reads\t%s\n' "$READS" >> "$OUT/dataset.tsv"
printf 'tablereads\t%s\n' "$TABLE_READS" >> "$OUT/dataset.tsv"

printf 'Running Java BBNorm countup=t on truncated biological dataset; expecting the vendored crash...\n'
set +e
measure "$OUT/java.metrics.tsv" \
  java -Xmx1g -cp "$CP" jgi.KmerNormalize \
  "${COMMON[@]}" \
  "out=$OUT/java.keep1.fq" "out2=$OUT/java.keep2.fq" \
  "outlow=$OUT/java.low1.fq" "outlow2=$OUT/java.low2.fq" \
  "outmid=$OUT/java.mid1.fq" "outmid2=$OUT/java.mid2.fq" \
  "outhigh=$OUT/java.high1.fq" "outhigh2=$OUT/java.high2.fq" \
  "hist=$OUT/java.hist.tsv" "rhist=$OUT/java.rhist.tsv" \
  >"$OUT/java.stdout.log" 2>"$OUT/java.stderr.log"
java_status=$?
set -e

if [[ "$java_status" -eq 0 ]]; then
  echo 'Expected vendored Java countup=t biological run to fail, but it succeeded.' >&2
  exit 1
fi

grep -q 'NullPointerException' "$OUT/java.stderr.log"
grep -q 'normalizeInThread' "$OUT/java.stderr.log"
grep -q 'BBNorm terminated in an error state' "$OUT/java.stderr.log"

printf 'Running Rust bbnorm-rs countup=t on the same truncated biological dataset...\n'
measure "$OUT/rust.metrics.tsv" \
  target/release/bbnorm-rs \
  "${COMMON[@]}" \
  "out=$OUT/rust.keep1.fq" "out2=$OUT/rust.keep2.fq" \
  "outlow=$OUT/rust.low1.fq" "outlow2=$OUT/rust.low2.fq" \
  "outmid=$OUT/rust.mid1.fq" "outmid2=$OUT/rust.mid2.fq" \
  "outhigh=$OUT/rust.high1.fq" "outhigh2=$OUT/rust.high2.fq" \
  "hist=$OUT/rust.hist.tsv" "rhist=$OUT/rust.rhist.tsv" \
  >"$OUT/rust.stdout.log" 2>"$OUT/rust.stderr.log"

for suffix in keep1.fq keep2.fq low1.fq low2.fq mid1.fq mid2.fq high1.fq high2.fq hist.tsv rhist.tsv; do
  test -e "$OUT/rust.$suffix"
done
test -s "$OUT/rust.hist.tsv"
test -s "$OUT/rust.rhist.tsv"

"$PYTHON" - "$OUT/java.metrics.tsv" "$OUT/rust.metrics.tsv" <<'PY'
import sys

def read_metrics(path):
    elapsed, rss, status = open(path, encoding="utf-8").read().strip().split("\t")
    return float(elapsed), int(rss), int(status)

java_elapsed, java_rss, java_status = read_metrics(sys.argv[1])
rust_elapsed, rust_rss, rust_status = read_metrics(sys.argv[2])

print("\nTruncated biological countup Java guard summary:")
print("tool\telapsed_seconds\tmax_rss_kb\tstatus")
print(f"java\t{java_elapsed:.6f}\t{java_rss}\t{java_status}")
print(f"rust\t{rust_elapsed:.6f}\t{rust_rss}\t{rust_status}")
PY

printf '\nTruncated biological countup Java guard passed. Vendored Java crashed in normalizeInThread while Rust produced complete countup outputs.\n'
printf 'Dataset metadata: %s\n' "$OUT/dataset.tsv"
printf 'Outputs and logs: %s\n' "$OUT"

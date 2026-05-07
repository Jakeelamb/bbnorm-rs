#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT"

DATA1="vendor/BBTools-master/resources/sample1.fq.gz"
DATA2="vendor/BBTools-master/resources/sample2.fq.gz"
CP="vendor/BBTools-master/current"
OUT="${1:-tmp/real_use_lower_depth_parity}"
mkdir -p "$OUT"
rm -f "$OUT"/*

R1_FIXTURE="$OUT/mixed_depth.r1.fq"
R2_FIXTURE="$OUT/mixed_depth.r2.fq"
python - "$DATA1" "$DATA2" "$R1_FIXTURE" "$R2_FIXTURE" <<'PYFIXTURE'
import gzip
import sys

source1, source2, dest1, dest2 = sys.argv[1:5]

def n_free_records(path, count):
    records = []
    with gzip.open(path, "rt") as reader:
        while len(records) < count:
            record = [reader.readline(), reader.readline(), reader.readline(), reader.readline()]
            if not record[0]:
                break
            if "N" not in record[1] and "n" not in record[1] and len(record[1].strip()) >= 80:
                records.append(record)
    if len(records) != count:
        raise SystemExit(f"not enough N-free reads in {path}")
    return records

def mutate_middle(seq):
    bases = list(seq.strip())
    middle = len(bases) // 2
    bases[middle] = {"A": "C", "C": "G", "G": "T", "T": "A"}.get(bases[middle].upper(), "A")
    return "".join(bases) + "\n"

r1 = n_free_records(source1, 1)[0]
r2_records = n_free_records(source2, 8)
mutated_r1 = mutate_middle(r1[1])

with open(dest1, "w") as writer1, open(dest2, "w") as writer2:
    for pair in range(8):
        writer1.write(f"@mixed_r1_{pair}\n")
        writer1.write(mutated_r1 if pair == 0 else r1[1])
        writer1.write("+\n")
        writer1.write(r1[3])

        r2 = r2_records[pair]
        writer2.write(f"@mixed_r2_{pair}\n")
        writer2.write(r2[1])
        writer2.write("+\n")
        writer2.write(r2[3])
PYFIXTURE

COMMON=(
  "in=$R1_FIXTURE"
  "in2=$R2_FIXTURE"
  "passes=1"
  "k=31"
  "minq=0"
  "minprob=0"
  "min=0"
  "minkmers=1"
  "target=1"
  "max=999999999"
  "tossbadreads=t"
  "saverarereads=t"
  "lowthresh=1"
  "highthresh=1"
  "errordetectratio=2"
  "threads=1"
  "overwrite=t"
  "bits=32"
)

cargo build --quiet

run_case() {
  local name="$1"
  local use_lower_depth="$2"

  printf 'Running Java BBNorm paired use-lower-depth case %s...\n' "$name"
  local java_start
  java_start=$(python - <<'PYTIME'
import time
print(time.perf_counter_ns())
PYTIME
)
  java -Xmx1g -cp "$CP" jgi.KmerNormalize \
    "${COMMON[@]}" \
    "uselowerdepth=$use_lower_depth" \
    "out=$OUT/java.$name.keep1.fq" \
    "out2=$OUT/java.$name.keep2.fq" \
    "outt=$OUT/java.$name.toss1.fq" \
    "outt2=$OUT/java.$name.toss2.fq" \
    >"$OUT/java.$name.stdout.log" 2>"$OUT/java.$name.stderr.log"
  local java_end
  java_end=$(python - <<'PYTIME'
import time
print(time.perf_counter_ns())
PYTIME
)
  local java_ms=$(((java_end - java_start) / 1000000))

  printf 'Running Rust bbnorm-rs paired use-lower-depth case %s...\n' "$name"
  local rust_start
  rust_start=$(python - <<'PYTIME'
import time
print(time.perf_counter_ns())
PYTIME
)
  target/debug/bbnorm-rs \
    "${COMMON[@]}" \
    "uselowerdepth=$use_lower_depth" \
    "out=$OUT/rust.$name.keep1.fq" \
    "out2=$OUT/rust.$name.keep2.fq" \
    "outt=$OUT/rust.$name.toss1.fq" \
    "outt2=$OUT/rust.$name.toss2.fq" \
    >"$OUT/rust.$name.stdout.log" 2>"$OUT/rust.$name.stderr.log"
  local rust_end
  rust_end=$(python - <<'PYTIME'
import time
print(time.perf_counter_ns())
PYTIME
)
  local rust_ms=$(((rust_end - rust_start) / 1000000))

  cmp "$OUT/java.$name.keep1.fq" "$OUT/rust.$name.keep1.fq"
  cmp "$OUT/java.$name.keep2.fq" "$OUT/rust.$name.keep2.fq"
  cmp "$OUT/java.$name.toss1.fq" "$OUT/rust.$name.toss1.fq"
  cmp "$OUT/java.$name.toss2.fq" "$OUT/rust.$name.toss2.fq"

  printf 'Case %s passed. Java: %sms, Rust: %sms\n' "$name" "$java_ms" "$rust_ms"
}

run_case lower_depth t
run_case higher_depth f

printf 'Paired use-lower-depth parity passed. Outputs and logs: %s\n' "$OUT"

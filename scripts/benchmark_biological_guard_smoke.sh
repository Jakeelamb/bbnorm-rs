#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT"

OUT="${1:-tmp/biological_guard_smoke}"
PASS_OUT="$OUT/pass"
FAIL_OUT="$OUT/fail"
STDOUT_LOG="$OUT/fail.stdout.log"
STDERR_LOG="$OUT/fail.stderr.log"
PASS_THREAD_CASES="${PASS_THREAD_CASES:-${THREAD_CASES:-1 2 auto}}"
FAIL_THREAD_CASES="${FAIL_THREAD_CASES:-1}"

rm -rf "$OUT"
mkdir -p "$OUT"

READS=1000 THREAD_CASES="$PASS_THREAD_CASES" KEEP_OUTPUTS=0 MAX_RSS_KB=1000000 \
  scripts/benchmark_biological_dataset.sh "$PASS_OUT"

set +e
READS=1000 THREAD_CASES="$FAIL_THREAD_CASES" KEEP_OUTPUTS=0 MAX_RSS_KB=1 \
  scripts/benchmark_biological_dataset.sh "$FAIL_OUT" \
  >"$STDOUT_LOG" 2>"$STDERR_LOG"
status=$?
set -e

if [[ "$status" -ne 3 ]]; then
  printf 'Expected MAX_RSS_KB guard failure status 3, got %s\n' "$status" >&2
  exit 1
fi

grep -q 'Memory guard failed' "$STDERR_LOG"

if find "$PASS_OUT" "$FAIL_OUT" -maxdepth 1 -type f \( -name '*.fq' -o -name '*.fastq' \) | grep -q .; then
  printf 'Expected KEEP_OUTPUTS=0 to remove FASTQ outputs in pass and fail smoke outputs.\n' >&2
  exit 1
fi

printf 'Biological benchmark RSS guard smoke passed. Outputs and logs: %s\n' "$OUT"

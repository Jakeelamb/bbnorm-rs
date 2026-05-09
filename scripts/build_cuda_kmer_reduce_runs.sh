#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT"

NVCC="${NVCC:-nvcc}"
OUT="${1:-tmp/cuda_kmer_reduce_runs}"
PERSISTENT_OUT="${2:-tmp/cuda_kmer_reduce_runs_persistent}"
mkdir -p "$(dirname "$OUT")" "$(dirname "$PERSISTENT_OUT")"
"$NVCC" -O3 -std=c++17 scripts/cuda_kmer_reduce_runs.cu -o "$OUT"
"$NVCC" -O3 -std=c++17 scripts/cuda_kmer_reduce_runs_persistent.cu -o "$PERSISTENT_OUT"
printf '%s\n%s\n' "$OUT" "$PERSISTENT_OUT"

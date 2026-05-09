#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT"

NVCC="${NVCC:-nvcc}"
OUT="${1:-tmp/cuda_kmer_reduce_runs}"
mkdir -p "$(dirname "$OUT")"
"$NVCC" -O3 -std=c++17 scripts/cuda_kmer_reduce_runs.cu -o "$OUT"
printf '%s\n' "$OUT"

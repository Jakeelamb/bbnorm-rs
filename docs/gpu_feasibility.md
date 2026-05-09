# GPU Feasibility Notes

This repo does not have a production GPU backend. The first useful test is a
standalone feasibility probe for the part of bounded counting most likely to
benefit from the GPU: sort/reduce of packed canonical short k-mers.

## Probe

Use:

```bash
scripts/cuda_kmer_sort_reduce_probe.py \
  --reads 500000 \
  --outdir tmp/cuda_probe_500k_$(date +%Y%m%d) \
  --force-rebuild
```

The probe:

- Extracts canonical `k=31` short k-mers from paired FASTQ/FASTQ.GZ into a
  binary `u64` stream. The default `--extractor rust` path uses the crate's
  real parser and k-mer generation code; `--extractor python` is retained as a
  simple control.
- Compiles a tiny CUDA/Thrust program with `nvcc`.
- Compares CPU `std::sort` plus adjacent reduce against GPU
  host-to-device copy, `thrust::sort`, `thrust::reduce_by_key`, and
  device-to-host copy.
- Checks that CPU and GPU reduced `(key, count)` streams match by unique count
  and checksum.

Extraction time is reported separately from GPU timing. The Rust extractor is
still a helper-process feasibility path, but it is close enough to the real
production shape to make the result useful.

## Local Results

Machine:

- GPU: NVIDIA GeForce RTX 5070 Laptop GPU
- VRAM: 8151 MiB
- Driver: 595.71.05
- CUDA: 13.2

Input:

- `tmp/human_benchmark_8threads/human_GRCh38_500k_R1.fq.gz`
- `tmp/human_benchmark_8threads/human_GRCh38_500k_R2.fq.gz`
- `k=31`

Rust extractor:

| Read pairs | K-mers | Extract | CPU sort/reduce | GPU total | GPU sort | Match |
| ---: | ---: | ---: | ---: | ---: | ---: | --- |
| 1,000 | 240,000 | 0.008 s | 12.677 ms | 273.973 ms | 273.394 ms | true |
| 50,000 | 12,000,000 | 0.250 s | 682.257 ms | 27.748 ms | 8.846 ms | true |
| 500,000 | 120,000,000 | 2.367 s | 7660.800 ms | 235.400 ms | 66.522 ms | true |

Earlier Python extractor control:

| Read pairs | K-mers | CPU sort/reduce | GPU total | GPU sort | Match |
| ---: | ---: | ---: | ---: | ---: | --- |
| 1,000 | 240,000 | 16.803 ms | 337.756 ms | 337.156 ms | true |
| 50,000 | 12,000,000 | 733.955 ms | 321.658 ms | 302.764 ms | true |
| 100,000 | 24,000,000 | 1451.000 ms | 321.488 ms | 284.518 ms | true |
| 500,000 | 120,000,000 | 7917.670 ms | 517.186 ms | 338.936 ms | true |

The 500k Rust-extractor slice shows that GPU sort/reduce is feasible at the
same scale as the current publish benchmark: 120M canonical k-mers fit in VRAM
and reduce correctly. The small 1k run loses badly, so a production GPU path
should only activate above a large batch threshold.

## Interpretation

The promising path is not direct GPU count-min atomics. It is:

1. CPU/Rust extracts packed canonical k-mers in large batches.
2. GPU sorts and reduces them into `(key, count)` runs.
3. Rust replays reduced counts into the existing deterministic conservative
   sketch path, preserving the parity-safe update semantics.

This preserves the semantics that made deterministic sorted replay publishable,
while moving the expensive key ordering/reduction stage to the GPU.

The next production-grade experiment should stream Rust-side k-mer batches to
either:

- a CUDA helper process that consumes binary `u64` batches, or
- a feature-gated CUDA FFI path.

The helper-process route is the lower-risk first integration because it keeps
CUDA build requirements out of the default crates.io package.

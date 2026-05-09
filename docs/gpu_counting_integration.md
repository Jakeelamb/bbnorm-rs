# GPU Counting Integration Notes

`gpucounting=t gpuhelper=<path>` is an experimental full-tool integration path.
It is parity-safe on the tested bounded no-prefilter lanes, but it is not the
default path because CPU counting is still faster.

The current implementation:

1. Writes primary short k-mers with the same Rust parser and per-pair duplicate
   removal rules as CPU counting.
2. Preserves the CPU deterministic chunk boundaries.
3. Runs either a one-shot external CUDA/CUB helper built from
   `scripts/cuda_kmer_reduce_runs.cu` for each chunk, or a persistent helper
   built from `scripts/cuda_kmer_reduce_runs_persistent.cu` when
   `gpupersistent=t` is set.
4. Replays each chunk's reduced `(u64 kmer, u32 count)` runs into the existing
   Rust count-min sketch before advancing to the next chunk.

Build the helper with:

```bash
scripts/build_cuda_kmer_reduce_runs.sh
```

Then run an explicit experiment:

```bash
target/release/bbnorm-rs \
  in=tmp/human_benchmark_8threads/human_GRCh38_500k_R1.fq.gz \
  in2=tmp/human_benchmark_8threads/human_GRCh38_500k_R2.fq.gz \
  reads=50000 tablereads=50000 passes=1 threads=8 zipthreads=1 \
  autocountmin=t autocountminreads=1 gpucounting=t \
  gpuhelper=tmp/cuda_kmer_reduce_runs \
  out=null out2=null hist=tmp/gpu.hist.tsv rhist=tmp/gpu.rhist.tsv
```

To use the persistent helper protocol:

```bash
target/release/bbnorm-rs \
  in=tmp/human_benchmark_8threads/human_GRCh38_500k_R1.fq.gz \
  in2=tmp/human_benchmark_8threads/human_GRCh38_500k_R2.fq.gz \
  reads=50000 tablereads=50000 passes=1 threads=8 zipthreads=1 \
  autocountmin=t autocountminreads=1 gpucounting=t gpupersistent=t \
  gpuhelper=tmp/cuda_kmer_reduce_runs_persistent \
  out=null out2=null hist=tmp/gpu.hist.tsv rhist=tmp/gpu.rhist.tsv
```

## Current Finding

The first naive global GPU reduce path was semantically close but not parity-safe for
deterministic conservative count-min sketches. Conservative sketch updates are
collision-order-sensitive. CPU deterministic counting replays reduced keys in
bounded chunks; the naive GPU helper globally sorts and reduces the whole input
before replay. That changes collision order.

Observed on the 50k paired-human lane:

| Metric | CPU deterministic | Global GPU reduce |
| --- | ---: | ---: |
| Wall | 2.741 s | 3.439 s |
| Input counting stage | 2.107 s | 2.908 s |
| Kept reads | 3,510 | 3,512 |
| Hist absolute raw delta | 314 |  |
| Hist absolute unique delta | 117 |  |
| Rhist absolute read delta | 12 |  |
| Rhist absolute base delta | 1,800 |  |

The 1k smoke lane was byte-identical, but the 50k lane exposes real collision
order drift.

The chunk-preserving one-shot path fixes that drift on the same 50k lane:

| Metric | CPU deterministic | Chunked GPU reduce |
| --- | ---: | ---: |
| Wall | 2.371 s | 5.564 s |
| Input counting stage | 1.708 s | 4.994 s |
| Kept reads | 3,510 | 3,510 |
| Hist comparison | identical | identical |
| Rhist comparison | identical | identical |

This is now a parity-safe integration probe for the tested bounded no-prefilter
lane, but it is not fast yet. The slowdown comes from launching the external
CUDA helper once per deterministic chunk.

The persistent helper removes that repeated process/CUDA context startup cost
while preserving the same chunk replay contract:

| Dataset | Metric | CPU deterministic | One-shot GPU reduce | Persistent GPU reduce |
| --- | --- | ---: | ---: | ---: |
| 50k read pairs | Wall | 2.142 s | 5.544 s | 2.569 s |
| 50k read pairs | Input counting stage | 1.589 s | 4.991 s | 2.027 s |
| 50k read pairs | Kept reads | 3,510 | 3,510 | 3,510 |
| 50k read pairs | Input unique kmers | 11,140,933 | 11,140,933 | 11,140,933 |
| 50k read pairs | Hist/rhist comparison | baseline | identical | identical |
| 500k read pairs | Wall | 21.499 s | 41.137 s | 25.616 s |
| 500k read pairs | Input counting stage | 16.029 s | 35.642 s | 19.705 s |
| 500k read pairs | Kept reads | 54,798 | 54,798 | 54,798 |
| 500k read pairs | Input unique kmers | 108,712,290 | 108,712,290 | 108,712,290 |
| 500k read pairs | Hist/rhist comparison | baseline | identical | identical |

The persistent helper is a real structural improvement over launching a helper
per chunk, but it still loses to the CPU implementation on this machine. The
remaining overhead is likely dominated by host/device copies, binary pipe I/O,
and per-chunk CUDA allocations inside the helper.

## Next Correct Target

The next performance target is keeping chunk-preserving semantics while removing
the remaining per-chunk overheads:

1. Keep the same deterministic chunk boundaries as CPU counting.
2. Reuse device buffers inside the persistent helper or move to a feature-gated
   FFI backend.
3. Replay each chunk's sorted reduced runs in chunk order.

That keeps the existing deterministic replay semantics while allowing CUB to
replace the expensive per-chunk sort/reduce work without changing conservative
count-min collision order.

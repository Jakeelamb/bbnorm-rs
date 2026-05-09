# Performance Frontier

This note records the current CPU hot-path evidence for the packed 16-bit
bounded human benchmark lane. It is intended to prevent future optimization
passes from repeating changes that looked plausible but did not improve the
measured result.

## Current Profile

Profile artifact:

```text
tmp/profile_500k_bits16_cpu_hotpath_20260508/
```

Command shape:

```text
target/release/bbnorm-rs
  in=tmp/human_benchmark_8threads/human_GRCh38_500k_R1.fq.gz
  in2=tmp/human_benchmark_8threads/human_GRCh38_500k_R2.fq.gz
  reads=500000 tablereads=500000 passes=1
  target=40 max=80 min=5 k=31
  threads=8 zipthreads=1 bits=16 mem=4279m
  autocountmin=t autocountminreads=1
  out=null outt=null out2=null outt2=null
```

The profile captured about 29k samples. The largest self-time symbols were:

| Symbol | Self cycles |
| --- | ---: |
| `PackedCountMinSketch::depth` via `CountLookup` | 35.42% |
| `PackedCountMinSketch::increment_16bit_3hash_conservative_and_return_unincremented` | 21.69% |
| `increment_pair_counts_with_prefilter` | 7.38% |
| `rayon::slice::sort::recurse` | 5.27% |
| `unfiltered_kmer_windows` | 4.09% |

That means the Java gap is mostly in random-access packed sketch lookup and
conservative update cost, not gzip, process overhead, or obvious allocation
churn.

## Negative Results

These were tested after the profile and rejected:

| Attempt | 500k bits16 median wall | Result |
| --- | ---: | --- |
| Baseline final refresh | 19.786s | Current README baseline |
| Direct large-sketch 16-bit cell update | 20.476s | Rejected; no input-counting win |
| Raw short-key `depth` trait fast path | 20.503s | Rejected; made lookup path slower |
| Packed 16-bit word-scan sparse histogram | 20.336s | Rejected; histogram change was noise/regression |
| `COUNT_PARALLEL_CHUNK_SIZE=16384` | 20.848s | Rejected; worse than 8192 |
| Raw short-kmer vector sort/reduce | 20.640s | Rejected; lower RSS but slower in the publish harness |
| Sharded raw short-kmer exact counter | 25.040s | Rejected; byte-identical output, lower RSS, but much slower input counting |

All rejected attempts preserved output shape in focused checks, but none beat
the existing benchmark median. They should not be reintroduced without a new
profile showing a different bottleneck.

## Next Serious Moves

Small generic micro-specialization is not enough, and a plain raw-vector
comparison sort is not enough either. The next plausible work needs a stronger
algorithmic change:

1. Replace comparison sorting with radix/counting-style short-kmer run
   formation, or use a sharded exact counter that preserves deterministic
   replay order without one global comparison sort per chunk.
2. Keep deterministic chunk replay order exactly the same.
3. Feed reduced `(raw_short_kmer, count)` runs directly into a raw packed
   16-bit/3-hash replay function.
4. Compare against Java and the current Rust baseline only with
   `scripts/benchmark_trustworthy_baseline.py`; ad hoc timing loops were too
   optimistic for this change.

The useful lesson from the raw short-kmer and sharded exact-counter
experiments is that removing `FxHashMap<KmerKey, u64>` can reduce memory, but
extra run formation, per-shard maps, and merge work are expensive enough to
lose badly in the publish lane.

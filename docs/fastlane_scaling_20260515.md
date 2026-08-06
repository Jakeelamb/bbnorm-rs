# BBNorm Java vs bbnorm-rs Fast-Lane Scaling

Benchmark lane: paired human GRCh38 slices, `bits=16`, `threads=8`, `zipthreads=1`, `passes=1`, `target=40`, `max=80`, `min=5`, null read outputs, `hist` and `rhist` enabled. Rust uses `deterministic=f` bounded approximate atomic packed sketch updates.

Each point is the median of 3 repeats from `scripts/benchmark_trustworthy_baseline.py --variants java,rust_nondeterministic`.

| Read pairs | Java s | Rust s | Speedup | Wall reduction | Java RSS GiB | Rust RSS GiB | RSS reduction |
| ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| 1,000 | 0.612 | 0.152 | 4.02x | 75.1% | 3.20 | 2.78 | 12.9% |
| 5,000 | 0.769 | 0.204 | 3.77x | 73.5% | 3.22 | 2.78 | 13.6% |
| 10,000 | 0.870 | 0.305 | 2.85x | 64.9% | 3.27 | 2.79 | 14.6% |
| 50,000 | 1.535 | 0.813 | 1.89x | 47.1% | 3.39 | 2.79 | 17.6% |
| 100,000 | 2.412 | 1.472 | 1.64x | 39.0% | 3.39 | 2.79 | 17.8% |
| 250,000 | 4.774 | 3.737 | 1.28x | 21.7% | 3.41 | 2.79 | 18.2% |
| 500,000 | 8.468 | 6.613 | 1.28x | 21.9% | 3.41 | 2.79 | 18.2% |

At 500k read pairs, Rust is 1.28x faster (21.9% lower wall time) and uses 18.2% less peak RSS than Java on this lane.

Artifacts:
- `fastlane_scaling_20260515.tsv`: aggregate table
- `assets/fastlane_wall_time.svg`: wall-clock scaling
- `assets/fastlane_speedup.svg`: Java/Rust speedup by read count
- `assets/fastlane_rss.svg`: peak RSS scaling

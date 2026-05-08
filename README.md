# bbnorm-rs

`bbnorm-rs` is a Rust implementation of the practical BBNorm read-depth
normalization workflow from BBTools. It focuses on local FASTA/FASTQ
normalization, histogram output, paired/interleaved routing, bounded memory
counting, and reproducible Java-parity behavior for covered modes.

This is an early working release, not a complete BBTools replacement. The Git
repository includes a vendored BBTools snapshot for parity tests; crates.io
packages intentionally exclude that snapshot to keep the package small.

## Install

From crates.io, once published:

```bash
cargo install bbnorm-rs
```

From a checkout:

```bash
cargo install --path .
```

## Basic Usage

```bash
bbnorm-rs in=reads_R1.fq.gz in2=reads_R2.fq.gz \
  out=normalized_R1.fq.gz out2=normalized_R2.fq.gz \
  target=40 max=80 min=5 k=31 threads=8
```

Common outputs:

```bash
bbnorm-rs in=reads.fq.gz out=keep.fq.gz outt=toss.fq.gz \
  hist=depth.tsv rhist=read_depth.tsv peaksout=peaks.tsv
```

Supported inputs and outputs include plain and gzip FASTA/FASTQ, single-end,
paired two-file, explicit interleaved, auto-interleaved, BBTools-style `#`
paired filename expansion, `null` output sinks, and common BBNorm-style
`key=value` aliases.

## Current Status

Verified locally:

- `cargo fmt --all -- --check`
- `cargo clippy --all-targets --all-features -- -D warnings`
- `cargo test --all`

Current tests cover 242 library tests, 8 basic integration tests, and 108
Java-parity tests against the vendored BBTools snapshot.

Implemented working areas include:

- BBNorm-style CLI parsing for common normalization options and aliases.
- Plain and gzip FASTA/FASTQ I/O.
- Single-end, paired, and interleaved read routing.
- Short canonical k-mers and long BBTools-shaped hashed k-mers.
- Exact and bounded count-min counting paths.
- Depth histograms, read-depth histograms, and peak output.
- Deterministic normalization decisions for covered modes.
- Multipass, count-up, and ECC behavior for tested subsets.
- Guarded benchmark and parity harnesses in the source repository.

Known gaps remain:

- Full BBTools sketch, prefilter, and cardinality/loglog collision parity is
  not complete.
- ECC and overlap behavior is covered by compact and biological stress tests,
  but not every BBMerge/BBNorm edge case.
- Large human-read benchmarks show excellent memory usage, but input counting
  remains the main speed bottleneck in some comparable modes.

See [docs/parity.md](docs/parity.md) and
[docs/component_buildout.md](docs/component_buildout.md) for the detailed
ledger. The acceptance matrix in
[docs/parity_matrix.md](docs/parity_matrix.md) is the current front-door
workflow for deciding whether a mode is exact parity, bounded approximate
parity, accepted Rust-over-Java divergence, or still a gap.

## Benchmark Snapshot

The acceptance matrix is the publishable benchmark source of truth. The latest
local matrix run at
`tmp/parity_acceptance_publish_ready_20260508/acceptance_summary.tsv` verified
9 bundled phiX exact-output modes and 6 local human bounded-sketch modes.

Exact bundled rows:

- `default`, `k=40`, `k=40 fixspikes=t`
- `passes=2`
- `keepall=t`
- `ecc=t markuncorrectableerrors=t`
- `qtrim=r trimq=10`
- `minlen=100`
- `passes=2 ecc=t markuncorrectableerrors=t`

Local human bounded-sketch rows:

| Row | Mode | Verdict | Java time | Rust time | Java RSS | Rust RSS | Hist drift ppm | Rhist drift ppm |
| --- | --- | --- | ---: | ---: | ---: | ---: | ---: | ---: |
| 50k | default | bounded_drift | 1.54s | 2.57s | 3.35 GiB | 3.30 GiB | 4 | 840 |
| 50k | prefilter | bounded_drift | 2.05s | 2.68s | 3.35 GiB | 3.10 GiB | 4 | 840 |
| 50k | k40_fixspikes | bounded_drift | 1.64s | 2.47s | 3.34 GiB | 3.12 GiB | 2 | 140 |
| 500k | default | bounded_drift | 9.03s | 30.12s | 3.38 GiB | 3.38 GiB | 1227 | 1492 |
| 500k | prefilter | bounded_drift | 10.53s | 31.81s | 3.26 GiB | 3.14 GiB | 49 | 998 |
| 500k | k40_fixspikes | bounded_drift | 10.82s | 28.60s | 3.39 GiB | 3.25 GiB | 245 | 982 |

The conservative publishable claim is: exact covered fixture modes match the
vendored Java oracle byte-for-byte, bounded human rows stay within the matrix
drift gate, and large-slice Rust speed still needs work in the input-counting
hot path. `countup=t` is tracked separately as an accepted Rust-over-Java
divergence guard rather than normal Java parity.

For high-throughput bounded approximate runs where byte-stable collision order
is less important than speed, `deterministic=f` enables direct parallel sketch
updates. On the local 50k human-pair no-output benchmark, default deterministic
mode measured 3.18s / 0.87 GiB RSS after chunk-size tuning, while
`deterministic=f` measured 1.42s / 0.40 GiB RSS with the same read limit and
memory settings.

For repeatable current baselines, use
[`scripts/benchmark_trustworthy_baseline.py`](scripts/benchmark_trustworthy_baseline.py).
It records git/tool/input metadata, command lines, raw run data, stage timings,
Java/Rust histogram drift, and aggregate p10/median/p90 summaries. See
[`docs/trustworthy_benchmarking.md`](docs/trustworthy_benchmarking.md) for the
Java-default and packed 16-bit benchmark lanes.

## Repository Layout

- `src/`: Rust library and CLI implementation.
- `tests/basic.rs`: package-friendly integration tests.
- `tests/java_parity.rs`: repository-only Java parity tests requiring
  `vendor/BBTools-master`.
- `docs/`: implementation status and parity notes.
- `scripts/`: local parity, benchmark, and stress harnesses.
- `vendor/`: BBTools reference snapshot for repository testing only.

## License

`bbnorm-rs` is licensed under the BSD 3-Clause License. The vendored BBTools
reference snapshot in the source repository is distributed under its own license
at `vendor/BBTools-master/license.txt` and is not included in crates.io
packages.

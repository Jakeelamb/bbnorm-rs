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

Current tests cover 242 library tests, 8 basic integration tests, and 106
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
ledger.

## Benchmark Snapshot

Latest local human-slice scaling results show the tradeoff clearly:

| Read pairs | Java time | Rust time | Java RSS | Rust RSS |
| ---: | ---: | ---: | ---: | ---: |
| 1k | 0.71s | 0.10s | 3.28 GiB | 0.42 GiB |
| 10k | 1.02s | 0.66s | 3.30 GiB | 0.67 GiB |
| 50k | 1.53s | 3.04s | 3.37 GiB | 0.81 GiB |
| 500k | 8.12s | 30.97s | 3.37 GiB | 0.81 GiB |

Other bounded/no-output and mode-specific runs are faster, but the conservative
publishable claim is: memory use is already much lower than vendored Java in
these local runs, while large-slice speed still needs work in the input-counting
hot path.

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

# Trustworthy Benchmarking

Use `scripts/benchmark_trustworthy_baseline.py` when benchmark numbers need to
be published, compared across commits, or used to choose the next optimization.
It runs repeated Java/Rust benchmarks, records metadata, and writes aggregate
p10/median/p90 summaries.

The harness intentionally records:

- Git commit, branch, status, host, tool versions, command lines, and input
  sizes/checksums in `metadata.json` and `commands.tsv`.
- One run directory per repeat and variant under `runs/`.
- Raw wall-clock/RSS and selected stage timings in `raw_runs.tsv`.
- Complete long-form stage timings in `stage_timings.tsv`.
- Java-vs-Rust histogram drift in `comparisons.tsv`.
- Aggregate min/p10/median/p90/max/mean values in `summary.tsv`.

Default variants are:

- `java`
- `rust_deterministic`
- `rust_nondeterministic`

The Rust commands include `autocountmin=t autocountminreads=1` by default so
bounded-counting probes actually exercise the bounded sketch path at small and
medium read limits.

## Java-default bounded lane

This compares Rust's default bounded-sketch shape against Java's default
bounded-sketch shape.

```bash
scripts/benchmark_trustworthy_baseline.py \
  --repeats 3 \
  --reads 100000 \
  --table-reads 100000 \
  --outdir tmp/trustworthy_baseline_100k_$(date +%Y%m%d)
```

## Packed 16-bit bounded lane

This is the lane that exercises the packed 16-bit Rust sketch hot path.

```bash
scripts/benchmark_trustworthy_baseline.py \
  --repeats 3 \
  --reads 100000 \
  --table-reads 100000 \
  --bits 16 \
  --outdir tmp/trustworthy_baseline_100k_bits16_$(date +%Y%m%d)
```

## Faster smoke check

Use a one-repeat 1k run to validate the harness after edits. This is not a
publishable benchmark.

```bash
scripts/benchmark_trustworthy_baseline.py \
  --repeats 1 \
  --reads 1000 \
  --table-reads 1000 \
  --skip-input-sha256 \
  --outdir tmp/trustworthy_baseline_smoke
```

## Interpreting Results

Use `summary.tsv` for headline numbers. Prefer median stage timings over single
runs. For bounded approximate modes, use `comparisons.tsv` to report drift
instead of expecting byte-identical histograms.

Key columns:

- `elapsed_seconds`: measured wall time for the full command.
- `max_rss_kb`: peak sampled resident memory for the process tree.
- `stage_input_counting`: Rust input table/sketch build time.
- `stage_normalize`: Rust normalization time.
- `stage_table_creation` and `stage_table_read`: Java table construction and
  read/normalization stages.

Do not mix lanes when making performance claims. Java-default `bits=32` atomic
Rust runs and packed `bits=16` Rust runs exercise different data structures.

# Parity Acceptance Matrix

This repo treats BBNorm parity as an acceptance matrix rather than a vague
"complete replacement" claim. Each row names a dataset, a set of modes, and the
comparison policy required for those modes.

## Verdict Tiers

- `exact_match`: Java and Rust both completed, `hist` and `rhist` match
  byte-for-byte, and sequence outputs match when the row enables
  `require_sequence_outputs`.
- `bounded_drift`: Java and Rust both completed, and histogram/read-depth drift
  stays under the configured gate. This is the expected class for bounded
  approximate count-min/prefilter rows.
- `rust_better_java_crashes`: Rust completed but the vendored Java oracle failed
  or was intentionally skipped for a known-broken mode.
- `fail_*`: the row violated its stated policy.
- `skipped_missing_dataset`: the row references a local dataset that is not
  present on this machine.

The default matrix is defined in
[`docs/parity_acceptance_matrix.tsv`](parity_acceptance_matrix.tsv).

## Running

Fast bundled fixture:

```bash
ROW_CASES=phix_exact_core scripts/parity_acceptance_matrix.py
```

Core local matrix:

```bash
scripts/parity_acceptance_matrix.py
```

Disabled guard rows are excluded from the default pass. Run them explicitly
when checking accepted Rust-over-Java divergences:

```bash
ROW_CASES=human_countup_guard_10k MATRIX_INCLUDE_DISABLED=1 scripts/parity_acceptance_matrix.py
```

The driver writes `acceptance_summary.tsv` plus the underlying Java/Rust
artifact directories under `tmp/parity_acceptance_matrix_<timestamp>/`.

## Current Scope

The enabled matrix covers bundled exact-output probes for default, `k=40`,
`fixspikes`, `passes=2`, `keepall`, `ecc=t markuncorrectableerrors=t`, right
quality trimming, `minlen`, and multipass ECC. It also covers local human
bounded-sketch rows at 50k and 500k read pairs when the local dataset is
present.

`countup=t` is intentionally not treated as normal parity in this matrix. The
guard row classifies it as `rust_better_java_crashes`: Rust must complete, but
the vendored Java oracle is skipped or expected to fail for that probe.

The 500k row is deliberately included as a publish-readiness check rather than
a speed victory lap. Current local evidence shows bounded drift within gate, but
Rust remains slower than Java there because input counting dominates wall time.

## Development Rule

New performance claims should be made only against rows that are already
`exact_match` or accepted `bounded_drift`. New feature work should add or update
rows in the TSV before being treated as covered parity.

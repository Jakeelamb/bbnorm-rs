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

The driver writes `acceptance_summary.tsv` plus the underlying Java/Rust
artifact directories under `tmp/parity_acceptance_matrix_<timestamp>/`.

## Development Rule

New performance claims should be made only against rows that are already
`exact_match` or accepted `bounded_drift`. New feature work should add or update
rows in the TSV before being treated as covered parity.

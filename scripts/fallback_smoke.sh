#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT"

OUT="${1:-tmp/fallback_smoke}"
rm -rf "$OUT"
mkdir -p "$OUT"
printf 'stage\tstatus\toutput\n' > "$OUT/summary.tsv"

printf 'Running default multipass fallback smoke...\n'
scripts/parity_default_multipass_fallback_representative_dataset.sh "$OUT/default_multipass"
printf 'default_multipass\tpassed\t%s\n' "$OUT/default_multipass" >> "$OUT/summary.tsv"

printf '\nRunning paired multipass smoke...\n'
scripts/parity_multipass_real_dataset.sh "$OUT/multipass"
printf 'paired_multipass\tpassed\t%s\n' "$OUT/multipass" >> "$OUT/summary.tsv"

printf '\nRunning countup behavior smoke...\n'
scripts/parity_countup_fallback_real_dataset.sh "$OUT/countup"
printf 'countup_behavior\tpassed\t%s\n' "$OUT/countup" >> "$OUT/summary.tsv"

printf '\nRunning ECC behavior smoke...\n'
scripts/parity_ecc_fallback_real_dataset.sh "$OUT/ecc"
printf 'ecc_behavior\tpassed\t%s\n' "$OUT/ecc" >> "$OUT/summary.tsv"

printf '\nRunning real-derived ECC stress smoke...\n'
scripts/parity_ecc_real_derived_stress.sh "$OUT/ecc_real_derived_stress"
printf 'ecc_real_derived_stress\tpassed\t%s\n' "$OUT/ecc_real_derived_stress" >> "$OUT/summary.tsv"

printf '\nRunning real-derived multipass ECC stress smoke...\n'
scripts/parity_multipass_ecc_real_derived_stress.sh "$OUT/multipass_ecc_real_derived_stress"
printf 'multipass_ecc_real_derived_stress\tpassed\t%s\n' "$OUT/multipass_ecc_real_derived_stress" >> "$OUT/summary.tsv"

printf '\nRunning staged multipass marked-uncorrectable ECC smoke...\n'
scripts/parity_multipass_ecc_real_derived_stress.sh "$OUT/multipass_ecc_markuncorrectable"
printf 'multipass_ecc_markuncorrectable\tpassed\t%s\n' "$OUT/multipass_ecc_markuncorrectable" >> "$OUT/summary.tsv"

printf '\nRunning sketch-control compatibility smoke...\n'
scripts/parity_sketch_controls_fallback_real_dataset.sh "$OUT/sketch_controls"
printf 'sketch_controls_compatibility\tpassed\t%s\n' "$OUT/sketch_controls" >> "$OUT/summary.tsv"

printf '\nRunning kmer-table runtime fallback smoke...\n'
scripts/parity_kmer_table_runtime_fallback_real_dataset.sh "$OUT/kmer_table_runtime"
printf 'kmer_table_runtime_fallback\tpassed\t%s\n' "$OUT/kmer_table_runtime" >> "$OUT/summary.tsv"

printf '\nRunning wrapper-sampling fallback smoke...\n'
scripts/parity_sampling_options_fallback_representative_dataset.sh "$OUT/sampling_options"
printf 'sampling_options_fallback\tpassed\t%s\n' "$OUT/sampling_options" >> "$OUT/summary.tsv"

printf '\nRunning deterministic=f supported-mode smoke...\n'
scripts/parity_deterministic_fallback_representative_dataset.sh "$OUT/deterministic"
printf 'deterministic_mode\tpassed\t%s\n' "$OUT/deterministic" >> "$OUT/summary.tsv"

printf '\nRunning peak short-alias compatibility smoke...\n'
scripts/parity_peak_short_aliases_representative_dataset.sh "$OUT/peak_short_aliases"
printf 'peak_short_aliases\tpassed\t%s\n' "$OUT/peak_short_aliases" >> "$OUT/summary.tsv"

printf '\nRunning trimq comma fallback smoke...\n'
scripts/parity_trimq_comma_fallback_real_dataset.sh "$OUT/trimq_comma"
printf 'trimq_comma_fallback\tpassed\t%s\n' "$OUT/trimq_comma" >> "$OUT/summary.tsv"

printf '\nRunning MPI local fallback smoke...\n'
scripts/parity_mpi_fallback_real_dataset.sh "$OUT/mpi"
printf 'mpi_fallback\tpassed\t%s\n' "$OUT/mpi" >> "$OUT/summary.tsv"

printf '\nRunning pairing runtime fallback smoke...\n'
scripts/parity_pairing_runtime_fallback_real_dataset.sh "$OUT/pairing_runtime"
printf 'pairing_runtime_fallback\tpassed\t%s\n' "$OUT/pairing_runtime" >> "$OUT/summary.tsv"

printf '\nRunning preparser runtime no-op smoke...\n'
scripts/parity_preparser_runtime_noops_real_dataset.sh "$OUT/preparser_runtime"
printf 'preparser_runtime_noops\tpassed\t%s\n' "$OUT/preparser_runtime" >> "$OUT/summary.tsv"

printf '\nRunning pass-suffixed quality recalibration compatibility smoke...\n'
scripts/parity_quality_recal_suffix_real_dataset.sh "$OUT/quality_recal_suffix"
printf 'quality_recal_suffix\tpassed\t%s\n' "$OUT/quality_recal_suffix" >> "$OUT/summary.tsv"

printf '\nRunning config-file expansion smoke...\n'
scripts/parity_config_file_real_dataset.sh "$OUT/config_file"
printf 'config_file_expansion\tpassed\t%s\n' "$OUT/config_file" >> "$OUT/summary.tsv"

printf '\nRunning SAM/readgroup runtime no-op smoke...\n'
scripts/parity_sam_runtime_noops_real_dataset.sh "$OUT/sam_runtime"
printf 'sam_runtime_noops\tpassed\t%s\n' "$OUT/sam_runtime" >> "$OUT/summary.tsv"

printf '\nRunning shared file-alias smoke...\n'
scripts/parity_file_aliases_representative_dataset.sh "$OUT/file_aliases"
printf 'file_aliases\tpassed\t%s\n' "$OUT/file_aliases" >> "$OUT/summary.tsv"

printf '\nRunning stdin input smoke...\n'
scripts/parity_stdin_input_real_dataset.sh "$OUT/stdin_input"
printf 'stdin_input\tpassed\t%s\n' "$OUT/stdin_input" >> "$OUT/summary.tsv"

printf '\nRunning null-output smoke...\n'
scripts/parity_null_outputs_real_dataset.sh "$OUT/null_outputs"
printf 'null_outputs\tpassed\t%s\n' "$OUT/null_outputs" >> "$OUT/summary.tsv"

printf '\nRunning temporary-directory control smoke...\n'
scripts/parity_tmpdir_controls_real_dataset.sh "$OUT/tmpdir_controls"
printf 'tmpdir_controls\tpassed\t%s\n' "$OUT/tmpdir_controls" >> "$OUT/summary.tsv"

printf '\nRunning side-output stats fallback smoke...\n'
scripts/parity_side_output_stats_fallback_real_dataset.sh "$OUT/side_output_stats"
printf 'side_output_stats_fallback\tpassed\t%s\n' "$OUT/side_output_stats" >> "$OUT/summary.tsv"

printf '\nRunning cardinality/loglog bounded-estimate smoke...\n'
scripts/parity_cardinality_fallback_real_dataset.sh "$OUT/cardinality"
printf 'cardinality_estimate\tpassed\t%s\n' "$OUT/cardinality" >> "$OUT/summary.tsv"

printf '\nRunning I/O hint fallback smoke...\n'
scripts/parity_io_runtime_hints_fallback_real_dataset.sh "$OUT/io_runtime_hints"
printf 'io_runtime_hints_fallback\tpassed\t%s\n' "$OUT/io_runtime_hints" >> "$OUT/summary.tsv"

printf '\nRunning genome-build context fallback smoke...\n'
scripts/parity_genome_context_fallback_real_dataset.sh "$OUT/genome_context"
printf 'genome_context_fallback\tpassed\t%s\n' "$OUT/genome_context" >> "$OUT/summary.tsv"

printf '\nRunning diagnostic sizing fallback smoke...\n'
scripts/parity_diagnostic_sizing_fallback_real_dataset.sh "$OUT/diagnostic_sizing"
printf 'diagnostic_sizing_fallback\tpassed\t%s\n' "$OUT/diagnostic_sizing" >> "$OUT/summary.tsv"

printf '\nFallback smoke summary:\n'
column -t -s $'\t' "$OUT/summary.tsv" || cat "$OUT/summary.tsv"
printf '\nFallback smoke passed. Outputs and logs: %s\n' "$OUT"

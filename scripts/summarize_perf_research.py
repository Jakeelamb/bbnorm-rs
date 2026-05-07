#!/usr/bin/env python3
import csv
import sys
from collections import defaultdict
from pathlib import Path

STAGE_COLUMNS = [
    ("rust_input_counting_s", "input_counting"),
    ("rust_input_exact_counting_s", "input_exact_counting"),
    ("rust_input_prefilter_counting_s", "input_prefilter_counting"),
    ("rust_input_main_counting_s", "input_main_counting"),
    ("rust_input_hist_s", "input_hist"),
    ("rust_input_rhist_s", "input_rhist"),
    ("rust_normalize_s", "normalize"),
    ("rust_summary_counts_s", "summary_counts"),
    ("rust_output_hist_s", "output_hist"),
    ("rust_output_rhist_s", "output_rhist"),
    ("rust_countup_work_source_s", "countup_work_source"),
    ("rust_countup_normalize_s", "countup_normalize"),
]

HYPOTHESES = {
    "input_counting": "Primary hotspot is counting overall, but use the largest counting substage below it to choose the actual target. Usually this means main-counting or prefilter-counting hot loops, not orchestration.",
    "input_main_counting": "Primary hotspot is main counting. Attack hot-loop update path, hash/key materialization, duplicate suppression, and atomic/packed table update shape before touching higher-level logic.",
    "input_prefilter_counting": "Primary hotspot is prefilter counting. Attack prefilter update path and gate placement; avoid paying for keys that the prefilter will reject.",
    "summary_counts": "Primary hotspot is summary scanning. Attack full-table scans with occupied-cell tracking, sparse summaries, or fused summary generation.",
    "input_hist": "Primary hotspot is input histogram generation. Fuse passes and keep histogram state sparse/compact where exact dense bins are unnecessary.",
    "input_rhist": "Primary hotspot is read-depth histogram generation. Prefer sparse reducers and combined collectors instead of a second full read pass.",
    "normalize": "Primary hotspot is normalization. Attack repeated depth lookups, pair analysis, output routing decisions, and avoid post-decision re-analysis.",
    "output_hist": "Primary hotspot is output histogram generation. Favor sparse kept-count reporting, combined hist/rhist collection, and bounded output side sketches.",
    "output_rhist": "Primary hotspot is output read-depth histogram generation. Same story: combine passes and avoid dense output-side buffers.",
    "countup_work_source": "Primary hotspot is count-up work-source construction. Attack presort analysis, chunking, temp-run payload size, and duplicated metadata in sort keys.",
    "countup_normalize": "Primary hotspot is count-up normalization. Fuse decision planning with kept-table updates and skip optional post-analysis when outputs do not need it.",
}


def f(x):
    try:
        return float(x)
    except Exception:
        return 0.0


def main():
    if len(sys.argv) != 3:
        raise SystemExit("Usage: summarize_perf_research.py <summary.tsv> <out.md>")
    summary_path = Path(sys.argv[1])
    out_path = Path(sys.argv[2])
    with summary_path.open(newline="") as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))

    by_mode = []
    stage_totals = defaultdict(float)
    mode_lines = []
    for row in rows:
        rust_seconds = f(row.get("rust_seconds"))
        java_seconds = f(row.get("java_seconds"))
        stage_pairs = []
        for col, label in STAGE_COLUMNS:
            val = f(row.get(col))
            if val > 0:
                stage_pairs.append((label, val))
                stage_totals[label] += val
        stage_pairs.sort(key=lambda kv: kv[1], reverse=True)
        top = ", ".join(f"{k}={v:.3f}s" for k, v in stage_pairs[:3]) if stage_pairs else "no stage timings"
        ratio = (rust_seconds / java_seconds) if java_seconds else 0.0
        by_mode.append((row.get("mode", "?"), ratio, top, row))
        mode_lines.append(f"- `{row.get('mode','?')}`: rust={rust_seconds:.3f}s, java={java_seconds:.3f}s, ratio={ratio:.2f}x, top stages: {top}")

    ranked = sorted(stage_totals.items(), key=lambda kv: kv[1], reverse=True)
    bottleneck_lines = [f"- `{name}` total={total:.3f}s" for name, total in ranked[:6]]
    lead = ranked[0][0] if ranked else None
    actionable_lead = lead
    if lead == "input_counting":
        for candidate in ("input_main_counting", "input_prefilter_counting", "input_exact_counting"):
            if stage_totals.get(candidate, 0.0) > 0:
                actionable_lead = candidate
                break
    hypothesis = HYPOTHESES.get(actionable_lead or lead, "No strong bottleneck inferred.")

    text = []
    text.append("# Performance Research Summary\n")
    text.append("This artifact is the repo's local 'deep research' equivalent for performance work: capture comparable measurements, preserve raw evidence, and force the next optimization target to be chosen from measured stage timings instead of vibes.\n")
    text.append("## Mode results\n")
    text.extend(mode_lines)
    text.append("\n## Aggregated Rust bottlenecks\n")
    text.extend(bottleneck_lines if bottleneck_lines else ["- none"])
    text.append("\n## Lead hypothesis\n")
    text.append(f"- Aggregate lead stage: `{lead}`" if lead else "- Aggregate lead stage: none")
    text.append(f"- Actionable lead stage: `{actionable_lead}`" if actionable_lead else "- Actionable lead stage: none")
    text.append(f"- Recommendation: {hypothesis}")
    text.append("\n## Next experiment rule\n")
    text.append("- Change exactly one hotspot mechanism at a time.")
    text.append("- Re-run the same benchmark profile.")
    text.append("- Keep the change only if elapsed time or peak RSS improves without breaking the intended drift/parity gate.")
    out_path.write_text("\n".join(text) + "\n")


if __name__ == "__main__":
    main()

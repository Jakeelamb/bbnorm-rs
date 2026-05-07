use anyhow::Result;
use bbnorm_rs::{parse_args, run};
use std::fs;
use std::io::{self, Write};
use std::path::Path;
use tempfile::{Builder, NamedTempFile};

fn main() {
    if let Err(err) = try_main() {
        eprintln!("error: {err:#}");
        std::process::exit(1);
    }
}

fn try_main() -> Result<()> {
    let mut config = parse_args(std::env::args_os().skip(1))?;
    let _stdin_spool = materialize_stdin_input(&mut config)?;
    if let Some(threads) = config.threads {
        rayon::ThreadPoolBuilder::new()
            .num_threads(threads)
            .build_global()?;
    }
    for note in &config.notes {
        eprintln!("note: {note}");
    }

    let summary = run(&config)?;
    let input_unique_split = summary
        .unique_kmers_in_split
        .map(|split| {
            format!(
                "; input unique kmers depth 1-{}: {}; input unique kmers depth {}+: {}",
                split.low_depth_max,
                split.low_depth_kmers,
                split.high_depth_min,
                split.high_depth_kmers
            )
        })
        .unwrap_or_default();
    let output_unique = summary
        .unique_kmers_out
        .map(|n| format!("; output unique kmers: {n}"))
        .unwrap_or_default();
    for layout in &summary.sketch_layouts {
        let prefilter_limit = layout
            .prefilter_limit
            .map(|limit| limit.to_string())
            .unwrap_or_default();
        eprintln!(
            "Sketch layout: table={}; kind={}; cells={}; hashes={}; bits={}; arrays={}; cells_per_array={}; mask_seed={}; update={}; max_count={}; memory_bytes={}; prefilter_limit={}",
            layout.table,
            layout.kind,
            layout.cells,
            layout.hashes,
            layout.bits,
            layout.arrays,
            layout.cells_per_array,
            layout.mask_seed,
            layout.update_mode,
            layout.max_count,
            layout.memory_bytes,
            prefilter_limit,
        );
    }
    for timing in &summary.stage_timings {
        eprintln!(
            "Stage timing: name={}; seconds={:.6}",
            timing.name,
            timing.elapsed_micros as f64 / 1_000_000.0
        );
    }
    if let Some(cardinality) = &summary.cardinality_in {
        eprintln!(
            "Cardinality estimate: scope=input; k={}; buckets={}; unique_kmers={}",
            cardinality.k, cardinality.buckets, cardinality.estimated_unique_kmers
        );
    }
    if let Some(cardinality) = &summary.cardinality_out {
        eprintln!(
            "Cardinality estimate: scope=output; k={}; buckets={}; unique_kmers={}",
            cardinality.k, cardinality.buckets, cardinality.estimated_unique_kmers
        );
    }
    if summary.countup_spill.has_spills() {
        eprintln!(
            "Count-up spill: initial_runs={}; merge_runs={}; final_runs={}; bytes_written={}; peak_live_bytes={}; final_live_bytes={}",
            summary.countup_spill.initial_runs,
            summary.countup_spill.merge_runs,
            summary.countup_spill.final_runs,
            summary.countup_spill.bytes_written,
            summary.countup_spill.peak_live_bytes,
            summary.countup_spill.final_live_bytes,
        );
    }
    eprintln!(
        "Processed {} reads ({} bases); kept {} reads, tossed {} reads; input unique kmers: {}{}{}",
        summary.reads_in,
        summary.bases_in,
        summary.reads_kept,
        summary.reads_tossed,
        summary.unique_kmers_in,
        input_unique_split,
        output_unique
    );
    Ok(())
}

fn materialize_stdin_input(config: &mut bbnorm_rs::Config) -> Result<Option<NamedTempFile>> {
    let in1_stdin = config.in1.as_deref().is_some_and(is_stdin_input);
    let in2_stdin = config.in2.as_deref().is_some_and(is_stdin_input);
    if !in1_stdin && !in2_stdin {
        return Ok(None);
    }
    if in1_stdin && in2_stdin {
        anyhow::bail!(
            "in=stdin and in2=stdin cannot both read the same stream; use interleaved=t for paired reads in one stdin stream"
        );
    }

    let original = if in1_stdin {
        config.in1.as_ref().expect("checked in1").clone()
    } else {
        config.in2.as_ref().expect("checked in2").clone()
    };
    let mut builder = Builder::new();
    builder
        .prefix("bbnorm-rs-stdin-")
        .suffix(stdin_temp_suffix(&original));
    let mut spool = if config.use_temp_dir {
        if let Some(dir) = config.temp_dir.as_deref() {
            fs::create_dir_all(dir)?;
            builder.tempfile_in(dir)?
        } else {
            builder.tempfile()?
        }
    } else {
        builder.tempfile()?
    };
    let bytes = io::copy(&mut io::stdin().lock(), spool.as_file_mut())?;
    spool.as_file_mut().flush()?;
    let path = spool.path().to_path_buf();
    if in1_stdin {
        config.in1 = Some(path);
    } else {
        config.in2 = Some(path);
    }
    config.notes.push(format!(
        "{} was materialized from stdin into a temporary file ({bytes} bytes) so the Rust count/normalize/histogram passes can reread it safely",
        original.display()
    ));
    Ok(Some(spool))
}

fn is_stdin_input(path: &Path) -> bool {
    path == Path::new("-")
        || path.to_str().is_some_and(|text| {
            text.eq_ignore_ascii_case("stdin") || text.to_ascii_lowercase().starts_with("stdin.")
        })
}

fn stdin_temp_suffix(path: &Path) -> &'static str {
    let lower = path.to_string_lossy().to_ascii_lowercase();
    if lower.ends_with(".gz") {
        if lower.contains(".fa") || lower.contains(".fasta") {
            ".fa.gz"
        } else {
            ".fq.gz"
        }
    } else if lower.contains(".fa") || lower.contains(".fasta") {
        ".fa"
    } else {
        ".fq"
    }
}

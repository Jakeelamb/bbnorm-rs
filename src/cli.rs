use crate::seqio::JunkMode;
use anyhow::{Context, Result, bail};
use std::collections::VecDeque;
use std::ffi::OsString;
use std::fs;
use std::path::PathBuf;

pub const USAGE: &str = "bbnorm-rs: Rust BBNorm compatibility port\n\nUsage:\n  bbnorm-rs in=<reads.fq> out=<kept.fq> outt=<tossed.fq> hist=<hist.tsv> [passes=1]\n\nThis working Rust slice supports exact k-mer counting for small inputs, automatic bounded count-min input sketches for large inputs, explicit bounded sketches via cells/matrixbits/sketchmemory, conservative atomic bits=32 sketch insertion with packed small-bit fallbacks, constrained and memory-sized prefilter sketch collision behavior, deterministic normalization, managed multipass temp-file orchestration, count-up mode with bounded kept-count sketches when requested, table-based ECC for covered paths, hist/rhist/peaks output, low/mid/high depth bins, zlib-rs gzip, BBTools-style pigz/unpigz hooks when available, bounded cardinality/loglog estimates when requested, and Rayon worker controls including threads=auto/max/all. Wrapper-sampling requests fall back to the supported engine with notes.";

pub const CARDINALITY_DEFAULT_BUCKETS: usize = 2048;
pub const CARDINALITY_MAX_BUCKETS: usize = 1 << 26;

#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub struct CountMinSettings {
    pub cells: Option<usize>,
    pub hashes: Option<usize>,
    pub bits: Option<u8>,
    pub memory_bytes: Option<usize>,
}

#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub struct PrefilterSettings {
    pub enabled: bool,
    pub force_disabled: bool,
    pub cells: Option<usize>,
    pub hashes: Option<usize>,
    pub bits: Option<u8>,
    pub memory_bytes: Option<usize>,
    pub memory_fraction_micros: Option<u32>,
}

#[derive(Debug, Clone, PartialEq)]
pub struct CardinalitySettings {
    pub input: bool,
    pub output: bool,
    pub buckets: usize,
    pub k: Option<usize>,
    pub seed: u64,
    pub min_probability: f64,
}

impl Default for CardinalitySettings {
    fn default() -> Self {
        Self {
            input: false,
            output: false,
            buckets: CARDINALITY_DEFAULT_BUCKETS,
            k: None,
            seed: 0,
            min_probability: 0.0,
        }
    }
}

#[derive(Debug, Clone)]
pub struct Config {
    pub in1: Option<PathBuf>,
    pub in2: Option<PathBuf>,
    pub extra: Vec<PathBuf>,
    pub out1: Option<PathBuf>,
    pub out2: Option<PathBuf>,
    pub out_toss1: Option<PathBuf>,
    pub out_toss2: Option<PathBuf>,
    pub out_low1: Option<PathBuf>,
    pub out_low2: Option<PathBuf>,
    pub out_mid1: Option<PathBuf>,
    pub out_mid2: Option<PathBuf>,
    pub out_high1: Option<PathBuf>,
    pub out_high2: Option<PathBuf>,
    pub out_uncorrected1: Option<PathBuf>,
    pub out_uncorrected2: Option<PathBuf>,
    pub hist_in: Option<PathBuf>,
    pub hist_out: Option<PathBuf>,
    pub rhist_in: Option<PathBuf>,
    pub rhist_out: Option<PathBuf>,
    pub peaks_in: Option<PathBuf>,
    pub peaks_out: Option<PathBuf>,
    pub match_hist_out: Option<PathBuf>,
    pub insert_hist_out: Option<PathBuf>,
    pub quality_accuracy_hist_out: Option<PathBuf>,
    pub indel_hist_out: Option<PathBuf>,
    pub error_hist_out: Option<PathBuf>,
    pub quality_hist_out: Option<PathBuf>,
    pub base_quality_hist_out: Option<PathBuf>,
    pub quality_count_hist_out: Option<PathBuf>,
    pub average_quality_hist_out: Option<PathBuf>,
    pub overall_base_quality_hist_out: Option<PathBuf>,
    pub length_hist_out: Option<PathBuf>,
    pub gc_hist_out: Option<PathBuf>,
    pub base_hist_out: Option<PathBuf>,
    pub entropy_hist_out: Option<PathBuf>,
    pub identity_hist_out: Option<PathBuf>,
    pub k: usize,
    pub min_quality: u8,
    pub quality_in_offset: u8,
    pub quality_out_offset: u8,
    pub change_quality: bool,
    pub min_called_quality: u8,
    pub max_called_quality: u8,
    pub fake_quality: u8,
    pub fasta_wrap: usize,
    pub u_to_t: bool,
    pub to_upper_case: bool,
    pub lower_case_to_n: bool,
    pub dot_dash_x_to_n: bool,
    pub iupac_to_n: bool,
    pub fix_junk_and_iupac: bool,
    pub junk_mode: JunkMode,
    pub min_prob: f64,
    pub max_reads: Option<u64>,
    pub table_reads: Option<u64>,
    pub min_length: usize,
    pub trim_left: bool,
    pub trim_right: bool,
    pub trim_quality: f64,
    pub trim_optimal: bool,
    pub trim_optimal_bias: Option<f64>,
    pub trim_window: bool,
    pub trim_window_length: usize,
    pub trim_min_good_interval: usize,
    pub interleaved: bool,
    pub test_interleaved: bool,
    pub keep_all: bool,
    pub zero_bin: bool,
    pub deterministic: bool,
    pub rename_reads: bool,
    pub canonical: bool,
    pub remove_duplicate_kmers: bool,
    pub fix_spikes: bool,
    pub target_depth: u64,
    pub target_depth_first: Option<u64>,
    pub target_bad_percent_low: f64,
    pub target_bad_percent_high: f64,
    pub max_depth: Option<u64>,
    pub min_depth: u64,
    pub min_kmers_over_min_depth: usize,
    pub depth_percentile: f64,
    pub high_percentile: f64,
    pub low_percentile: f64,
    pub error_detect_ratio: u64,
    pub high_thresh: u64,
    pub low_thresh: u64,
    pub toss_error_reads: bool,
    pub toss_error_reads_first: bool,
    pub require_both_bad: bool,
    pub save_rare_reads: bool,
    pub discard_bad_only: bool,
    pub discard_bad_only_first: bool,
    pub error_correct: bool,
    pub error_correct_first: bool,
    pub error_correct_final: bool,
    pub overlap_error_correct: bool,
    pub overlap_error_correct_auto: bool,
    pub mark_errors_only: bool,
    pub mark_uncorrectable_errors: bool,
    pub trim_after_marking: bool,
    pub mark_with_one: bool,
    pub error_correct_ratio: u64,
    pub error_correct_high_thresh: u64,
    pub error_correct_low_thresh: u64,
    pub max_errors_to_correct: usize,
    pub max_quality_to_correct: u8,
    pub correct_from_left: bool,
    pub correct_from_right: bool,
    pub suffix_len: usize,
    pub prefix_len: usize,
    pub count_up: bool,
    pub add_bad_reads_countup: bool,
    pub use_lower_depth: bool,
    pub toss_by_low_true_depth: bool,
    pub low_bin_depth: i64,
    pub high_bin_depth: i64,
    pub hist_len: usize,
    pub side_hist_len: Option<usize>,
    pub gc_bins: Option<usize>,
    pub entropy_bins: usize,
    pub entropy_k: usize,
    pub entropy_window: usize,
    pub allow_entropy_ns: bool,
    pub identity_bins: usize,
    pub cardinality: CardinalitySettings,
    pub hist_columns: u8,
    pub print_zero_coverage: bool,
    pub peak_min_height: u64,
    pub peak_min_volume: u64,
    pub peak_min_width: usize,
    pub peak_min_peak: usize,
    pub peak_max_peak: usize,
    pub peak_max_count: usize,
    pub peak_ploidy: i32,
    pub overwrite: bool,
    pub append: bool,
    pub passes: usize,
    pub threads: Option<usize>,
    pub gzip_threads: Option<usize>,
    pub temp_dir: Option<PathBuf>,
    pub use_temp_dir: bool,
    pub max_countup_spill_initial_runs: Option<usize>,
    pub max_countup_spill_merge_runs: Option<usize>,
    pub max_countup_spill_final_runs: Option<usize>,
    pub max_countup_spill_live_bytes: Option<u64>,
    pub max_countup_spill_final_live_bytes: Option<u64>,
    pub max_countup_spill_write_bytes: Option<u64>,
    pub table_initial_size: Option<usize>,
    pub table_prealloc_fraction: Option<f64>,
    pub build_passes: usize,
    pub auto_count_min: bool,
    pub force_exact_counts: bool,
    pub auto_count_min_input_bytes: usize,
    pub auto_count_min_read_threshold: u64,
    pub auto_count_min_memory_bytes: Option<usize>,
    pub count_min: CountMinSettings,
    pub count_min_bits_first: Option<u8>,
    pub prefilter: PrefilterSettings,
    pub locked_increment: Option<bool>,
    pub notes: Vec<String>,
}

impl Default for Config {
    fn default() -> Self {
        Self {
            in1: None,
            in2: None,
            extra: Vec::new(),
            out1: None,
            out2: None,
            out_toss1: None,
            out_toss2: None,
            out_low1: None,
            out_low2: None,
            out_mid1: None,
            out_mid2: None,
            out_high1: None,
            out_high2: None,
            out_uncorrected1: None,
            out_uncorrected2: None,
            hist_in: None,
            hist_out: None,
            rhist_in: None,
            rhist_out: None,
            peaks_in: None,
            peaks_out: None,
            match_hist_out: None,
            insert_hist_out: None,
            quality_accuracy_hist_out: None,
            indel_hist_out: None,
            error_hist_out: None,
            quality_hist_out: None,
            base_quality_hist_out: None,
            quality_count_hist_out: None,
            average_quality_hist_out: None,
            overall_base_quality_hist_out: None,
            length_hist_out: None,
            gc_hist_out: None,
            base_hist_out: None,
            entropy_hist_out: None,
            identity_hist_out: None,
            k: 31,
            min_quality: 5,
            quality_in_offset: 33,
            quality_out_offset: 33,
            change_quality: true,
            min_called_quality: 2,
            max_called_quality: 50,
            fake_quality: 30,
            fasta_wrap: 70,
            u_to_t: false,
            to_upper_case: false,
            lower_case_to_n: false,
            dot_dash_x_to_n: false,
            iupac_to_n: false,
            fix_junk_and_iupac: false,
            junk_mode: JunkMode::Crash,
            min_prob: 0.5,
            max_reads: None,
            table_reads: None,
            min_length: 1,
            trim_left: false,
            trim_right: false,
            trim_quality: 5.0,
            trim_optimal: true,
            trim_optimal_bias: None,
            trim_window: false,
            trim_window_length: 4,
            trim_min_good_interval: 2,
            interleaved: false,
            test_interleaved: true,
            keep_all: false,
            zero_bin: false,
            deterministic: true,
            rename_reads: false,
            canonical: true,
            remove_duplicate_kmers: true,
            fix_spikes: false,
            target_depth: 100,
            target_depth_first: None,
            target_bad_percent_low: 0.85,
            target_bad_percent_high: 1.5,
            max_depth: None,
            min_depth: 5,
            min_kmers_over_min_depth: 15,
            depth_percentile: 0.54,
            high_percentile: 0.90,
            low_percentile: 0.25,
            error_detect_ratio: 125,
            high_thresh: 12,
            low_thresh: 3,
            toss_error_reads: false,
            toss_error_reads_first: false,
            require_both_bad: false,
            save_rare_reads: false,
            discard_bad_only: false,
            discard_bad_only_first: false,
            error_correct: false,
            error_correct_first: false,
            error_correct_final: false,
            overlap_error_correct: false,
            overlap_error_correct_auto: false,
            mark_errors_only: false,
            mark_uncorrectable_errors: false,
            trim_after_marking: false,
            mark_with_one: false,
            error_correct_ratio: 140,
            error_correct_high_thresh: 22,
            error_correct_low_thresh: 2,
            max_errors_to_correct: 3,
            max_quality_to_correct: 127,
            correct_from_left: true,
            correct_from_right: true,
            suffix_len: 3,
            prefix_len: 3,
            count_up: false,
            add_bad_reads_countup: false,
            use_lower_depth: true,
            toss_by_low_true_depth: true,
            low_bin_depth: 10,
            high_bin_depth: 80,
            hist_len: (1 << 20) + 1,
            side_hist_len: None,
            gc_bins: None,
            entropy_bins: 1000,
            entropy_k: 5,
            entropy_window: 50,
            allow_entropy_ns: true,
            identity_bins: 750,
            cardinality: CardinalitySettings::default(),
            hist_columns: 3,
            print_zero_coverage: false,
            peak_min_height: 2,
            peak_min_volume: 5,
            peak_min_width: 3,
            peak_min_peak: 2,
            peak_max_peak: i32::MAX as usize,
            peak_max_count: 10,
            peak_ploidy: -1,
            overwrite: false,
            append: false,
            passes: 2,
            threads: None,
            gzip_threads: None,
            temp_dir: None,
            use_temp_dir: false,
            max_countup_spill_initial_runs: None,
            max_countup_spill_merge_runs: None,
            max_countup_spill_final_runs: None,
            max_countup_spill_live_bytes: None,
            max_countup_spill_final_live_bytes: None,
            max_countup_spill_write_bytes: None,
            table_initial_size: None,
            table_prealloc_fraction: None,
            build_passes: 1,
            auto_count_min: true,
            force_exact_counts: false,
            auto_count_min_input_bytes: 32 * 1024 * 1024,
            auto_count_min_read_threshold: 250_000,
            auto_count_min_memory_bytes: None,
            count_min: CountMinSettings::default(),
            count_min_bits_first: None,
            prefilter: PrefilterSettings::default(),
            locked_increment: None,
            notes: Vec::new(),
        }
    }
}

pub fn parse_args<I>(args: I) -> Result<Config>
where
    I: IntoIterator<Item = OsString>,
{
    let mut config = Config::default();
    let mut positional = Vec::new();
    let mut saw_arg = false;
    let mut pending: VecDeque<OsString> = args.into_iter().collect();

    while let Some(raw) = pending.pop_front() {
        saw_arg = true;
        let arg = raw.to_string_lossy().into_owned();
        if arg == "-h" || arg == "--help" || arg.eq_ignore_ascii_case("help") {
            bail!(USAGE);
        }
        if let Some((key, value)) = arg.split_once('=') {
            let key = key.to_ascii_lowercase();
            if key == "config" {
                let expanded = read_config_args(value)?;
                config.notes.push(format!(
                    "config={value} expanded into {} BBTools-style argument line(s)",
                    expanded.len()
                ));
                for item in expanded.into_iter().rev() {
                    pending.push_front(OsString::from(item));
                }
            } else {
                handle_key_value(&mut config, &key, value)?;
            }
        } else if arg.eq_ignore_ascii_case("null") {
            // BBTools treats a literal "null" argument as an inert placeholder.
        } else if arg.eq_ignore_ascii_case("1pass") || arg.eq_ignore_ascii_case("1p") {
            config.passes = 1;
            config.notes.push("single-pass mode selected".to_string());
        } else if arg.eq_ignore_ascii_case("2pass") || arg.eq_ignore_ascii_case("2p") {
            config.passes = 2;
        } else {
            let key = arg.to_ascii_lowercase();
            if is_bare_boolean_key(&key) {
                handle_key_value(&mut config, &key, "t")?;
            } else {
                positional.push(PathBuf::from(arg));
            }
        }
    }

    if !saw_arg {
        bail!(USAGE);
    }

    if positional.len() > 2 {
        bail!(
            "expected at most two positional inputs; use in=<file> and in2=<file> for paired input"
        );
    }
    if config.in1.is_none() {
        config.in1 = positional.first().cloned();
    }
    if config.in2.is_none() {
        config.in2 = positional.get(1).cloned();
    }

    fill_default_gzip_threads(&mut config);
    validate(&mut config)?;
    Ok(config)
}

fn is_bare_boolean_key(key: &str) -> bool {
    matches!(
        key,
        "keepall"
            | "zerobin"
            | "deterministic"
            | "dr"
            | "det"
            | "rn"
            | "rename"
            | "renamereads"
            | "canonical"
            | "removeduplicatekmers"
            | "rdk"
            | "fixspikes"
            | "fs"
            | "tossbadreads"
            | "tosserrorreads"
            | "tbr"
            | "ter"
            | "requirebothbad"
            | "rbb"
            | "removeifeitherbad"
            | "rieb"
            | "saverarereads"
            | "srr"
            | "discardbadonly"
            | "dbo"
            | "uselowerdepth"
            | "uld"
            | "printzerocoverage"
            | "pzc"
            | "overwrite"
            | "ow"
            | "ignorebadquality"
            | "ibq"
            | "changequality"
            | "cq"
            | "utot"
            | "tuc"
            | "touppercase"
            | "lctn"
            | "lowercaseton"
            | "dotdashxton"
            | "undefinedton"
            | "iupacton"
            | "itn"
            | "fixjunk"
            | "ignorejunk"
            | "usebgzip"
            | "bgzip"
            | "usepigz"
            | "pigz"
            | "usegunzip"
            | "gunzip"
            | "ungzip"
            | "useunpigz"
            | "unpigz"
            | "useunbgzip"
            | "unbgzip"
            | "usegzip"
            | "gzip"
            | "usebgzf"
            | "bgzf"
            | "ordered"
            | "ord"
            | "verbose"
            | "printcoverage"
            | "append"
            | "app"
            | "interleaved"
            | "int"
            | "testinterleaved"
            | "forceinterleaved"
            | "prefilter"
            | "autocountmin"
            | "autosketch"
            | "autosketchtable"
            | "autosketchtables"
            | "exact"
            | "exactcount"
            | "exactcounts"
            | "useexact"
            | "sketchexact"
            | "auto"
            | "automatic"
            | "countup"
            | "abrc"
            | "addbadreadscountup"
            | "markerrors"
            | "markonly"
            | "meo"
            | "markuncorrectableerrors"
            | "markuncorrectable"
            | "mue"
            | "tam"
            | "trimaftermarking"
            | "markwith1"
            | "markwithone"
            | "mw1"
            | "aec"
            | "aecc"
            | "aggressiveerrorcorrection"
            | "cec"
            | "cecc"
            | "conservativeerrorcorrection"
            | "ecc"
            | "ecc1"
            | "ecc2"
            | "eccf"
            | "eccbyoverlap"
            | "ecco"
            | "overlap"
            | "cfl"
            | "cfr"
            | "cardinality"
            | "loglog"
            | "loglogin"
            | "cardinalityout"
            | "loglogout"
    )
}

fn fill_default_gzip_threads(config: &mut Config) {
    if config.gzip_threads.is_some() {
        return;
    }
    let Some(threads) = config.threads.filter(|threads| *threads > 1) else {
        return;
    };
    config.gzip_threads = Some(threads);
    config.notes.push(format!(
        "threads={threads} also enables gzip input/output workers up to {threads}; use zipthreads=1 to force single-thread gzip I/O"
    ));
}

fn read_config_args(value: &str) -> Result<Vec<String>> {
    let mut args = Vec::new();
    for file in value.split(',').filter(|part| !part.trim().is_empty()) {
        let path = PathBuf::from(file.trim());
        let text = fs::read_to_string(&path)
            .with_context(|| format!("could not process config file {}", path.display()))?;
        for line in text.lines() {
            let trimmed = line.trim();
            if !trimmed.is_empty() && !trimmed.starts_with('#') {
                args.push(trimmed.to_string());
            }
        }
    }
    Ok(args)
}

fn handle_key_value(config: &mut Config, key: &str, value: &str) -> Result<()> {
    match key {
        "in" | "input" | "in1" | "input1" => config.in1 = Some(path(value)),
        "in2" | "input2" => config.in2 = Some(path(value)),
        "extra" => config.extra.extend(extra_paths(value)),
        "out" | "output" | "out1" | "output1" | "outk" | "outkeep" | "outgood" => {
            config.out1 = Some(path(value))
        }
        "out2" | "output2" | "outk2" | "outkeep2" | "outgood2" => config.out2 = Some(path(value)),
        "outt" | "outt1" | "outtoss" | "outoss" | "outbad" => {
            config.out_toss1 = Some(path(value));
        }
        "outt2" | "outtoss2" | "outoss2" | "outbad2" => config.out_toss2 = Some(path(value)),
        "outl" | "outl1" | "outlow" | "outlow1" => config.out_low1 = Some(path(value)),
        "outl2" | "outlow2" => config.out_low2 = Some(path(value)),
        "outm" | "outm1" | "outmid" | "outmid1" | "outmiddle" => {
            config.out_mid1 = Some(path(value));
        }
        "outm2" | "outmid2" | "outmiddle2" => config.out_mid2 = Some(path(value)),
        "outh" | "outh1" | "outhigh" | "outhigh1" => config.out_high1 = Some(path(value)),
        "outh2" | "outhigh2" => config.out_high2 = Some(path(value)),
        "outu" | "outu1" | "outuncorrected" => config.out_uncorrected1 = Some(path(value)),
        "outu2" | "outuncorrected2" => config.out_uncorrected2 = Some(path(value)),
        "hist" | "histin" | "inhist" | "khist" => config.hist_in = Some(path(value)),
        "histout" | "outhist" | "hist2" | "khistout" => config.hist_out = Some(path(value)),
        "rhist" => config.rhist_in = Some(path(value)),
        "rhistout" => config.rhist_out = Some(path(value)),
        "peaks" => config.peaks_in = Some(path(value)),
        "peaksout" => config.peaks_out = Some(path(value)),
        "extin" | "extout" => {
            config.notes.push(format!(
                "{key}={value} is a BBTools file-extension hint; covered Rust paths infer FASTA/FASTQ format from explicit filenames"
            ));
        }
        "k" | "kmer" => config.k = parse_usize(value, key)?,
        "minq" | "minqual" => config.min_quality = parse_u8(value, key)?,
        "minprob" => config.min_prob = parse_f64(value, key)?,
        "reads" | "maxreads" => config.max_reads = parse_limit(value, key)?,
        "tablereads" | "buildreads" => config.table_reads = parse_limit(value, key)?,
        "ml" | "minlen" | "minlength" => config.min_length = parse_kmg_usize(value, key)?,
        "maxlength" | "maxreadlength" | "maxreadlen" | "maxlen" => {
            let _ = parse_kmg_usize(value, key)?;
            config.notes.push(format!(
                "{key}={value} is parsed by BBNorm but not used by KmerNormalize"
            ));
        }
        "mingc" | "maxgc" | "mlf" | "minlenfrac" | "minlenfraction" | "minlengthfraction" => {
            let _ = parse_f64(value, key)?;
            config.notes.push(format!(
                "{key}={value} is parsed by BBNorm but not used by KmerNormalize"
            ));
        }
        "usepairgc"
        | "pairgc"
        | "trimbadsequence"
        | "chastityfilter"
        | "cf"
        | "failnobarcode"
        | "averagequalitybyprobability"
        | "aqbp"
        | "untrim" => {
            let _ = parse_bool(value, key)?;
            config.notes.push(format!(
                "{key}={value} is parsed by BBNorm but not used by KmerNormalize"
            ));
        }
        "badbarcodes" | "barcodefilter" => {
            if !value.eq_ignore_ascii_case("crash") && !value.eq_ignore_ascii_case("fail") {
                let _ = parse_bool(value, key)?;
            }
            config.notes.push(format!(
                "{key}={value} is parsed by BBNorm but not used by KmerNormalize"
            ));
        }
        "barcodes" | "barcode" => {
            config.notes.push(format!(
                "{key}={value} is parsed by BBNorm but not used by KmerNormalize"
            ));
        }
        "maxns" => {
            let _ = parse_i64(value, key)?;
            config.notes.push(format!(
                "{key}={value} is parsed by BBNorm but not used by KmerNormalize"
            ));
        }
        "minconsecutivebases"
        | "mcb"
        | "minavgqualitybases"
        | "maqb"
        | "mintl"
        | "mintrimlen"
        | "mintrimlength" => {
            let _ = parse_usize(value, key)?;
            config.notes.push(format!(
                "{key}={value} is parsed by BBNorm but not used by KmerNormalize"
            ));
        }
        "minavgquality" | "minaveragequality" | "maq" => {
            parse_min_average_quality(value, key)?;
            config.notes.push(format!(
                "{key}={value} is parsed by BBNorm but not used by KmerNormalize"
            ));
        }
        "minbasequality" | "mbq" => {
            let _ = parse_i8(value, key)?;
            config.notes.push(format!(
                "{key}={value} is parsed by BBNorm but not used by KmerNormalize"
            ));
        }
        "build" | "genome" => {
            let _ = parse_i32(value, key)?;
            config.notes.push(format!(
                "{key}={value} is a BBTools genome-build context control; covered Rust FASTA/FASTQ normalization does not use reference build metadata"
            ));
        }
        "qtrim" | "qtrim1" | "qtrim2" => parse_qtrim(config, value, key)?,
        "trimq" | "trimquality" | "trimq2" => {
            config.trim_quality = parse_trim_quality(config, value, key)?
        }
        "trimleft" | "qtrimleft" => config.trim_left = parse_bool(value, key)?,
        "trimright" | "qtrimright" => config.trim_right = parse_bool(value, key)?,
        "optitrim" | "otf" | "otm" => parse_optitrim(config, value, key)?,
        "trimgoodinterval" => config.trim_min_good_interval = parse_usize(value, key)?,
        "trimclip" => {
            let _ = parse_bool(value, key)?;
            config.notes.push(format!(
                "{key}={value} is parsed by BBNorm but not used by its trimFast call"
            ));
        }
        "trimpolya" | "trimpolyg" | "trimpolygleft" | "trimpolygright" | "filterpolyg"
        | "trimpolyc" | "trimpolycleft" | "trimpolycright" | "filterpolyc" | "maxnonpoly" => {
            let _ = parse_poly(value, key)?;
            config.notes.push(format!(
                "{key}={value} is parsed by BBNorm but not used by KmerNormalize"
            ));
        }
        "forcetrimmod" | "forcemrimmodulo" | "ftm" | "ftl" | "forcetrimleft" | "ftr"
        | "forcetrimright" | "ftr2" | "forcetrimright2" => {
            let _ = parse_i64(value, key)?;
            config.notes.push(format!(
                "{key}={value} is parsed by BBNorm but not used by KmerNormalize"
            ));
        }
        "keepall" => config.keep_all = parse_bool(value, key)?,
        "zerobin" => config.zero_bin = parse_bool(value, key)?,
        "deterministic" | "dr" | "det" => {
            config.deterministic = parse_bool(value, key)?;
            if !config.deterministic {
                config.notes.push(format!(
                    "{key}={value} enables nondeterministic read selection and faster parallel replay for bounded approximate sketches"
                ));
            }
        }
        "rn" | "rename" | "renamereads" => config.rename_reads = parse_bool(value, key)?,
        "canonical" => config.canonical = parse_bool(value, key)?,
        "removeduplicatekmers" | "rdk" => config.remove_duplicate_kmers = parse_bool(value, key)?,
        "fixspikes" | "fs" => config.fix_spikes = parse_bool(value, key)?,
        "target" | "targetdepth" | "tgt" => config.target_depth = parse_u64(value, key)?,
        "max" | "maxdepth" => config.max_depth = Some(parse_u64(value, key)?),
        "min" | "mindepth" => config.min_depth = parse_u64(value, key)?,
        "minkmers" | "minkmersovermindepth" | "mingoodkmersperread" | "mgkpr" => {
            config.min_kmers_over_min_depth = parse_usize(value, key)?.max(1);
        }
        "percentile" | "depthpercentile" | "dp" => {
            config.depth_percentile = parse_percent(value, key)?;
        }
        "highdepthpercentile" | "highpercentile" | "hdp" => {
            config.high_percentile = parse_percent(value, key)?;
        }
        "lowdepthpercentile" | "lowpercentile" | "ldp" => {
            config.low_percentile = parse_percent(value, key)?;
        }
        "errordetectratio" | "edr" => config.error_detect_ratio = parse_u64(value, key)?,
        "highthresh" | "hthresh" | "ht" => config.high_thresh = parse_u64(value, key)?,
        "lowthresh" | "lthresh" | "lt" => config.low_thresh = parse_u64(value, key)?,
        "tossbadreads" | "tosserrorreads" | "tbr" | "ter" => {
            let enabled = parse_bool(value, key)?;
            config.toss_error_reads = enabled;
            config.toss_error_reads_first = enabled;
        }
        "tossbadreads2" | "tosserrorreads2" | "tbr2" | "ter2" | "tossbadreadsf"
        | "tosserrorreadsf" | "tbrf" | "terf" => {
            config.toss_error_reads = parse_bool(value, key)?;
        }
        "tossbadreads1" | "tosserrorreads1" | "tbr1" | "ter1" => {
            config.toss_error_reads_first = parse_bool(value, key)?;
        }
        "requirebothbad" | "rbb" => config.require_both_bad = parse_bool(value, key)?,
        "removeifeitherbad" | "rieb" => config.require_both_bad = !parse_bool(value, key)?,
        "saverarereads" | "srr" => config.save_rare_reads = parse_bool(value, key)?,
        "discardbadonly" | "dbo" | "discardbadonlyf" | "dbof" | "discardbadonly2" | "dbo2" => {
            let enabled = parse_bool(value, key)?;
            config.discard_bad_only = enabled;
            config.discard_bad_only_first = enabled;
        }
        "discardbadonly1" | "dbo1" => {
            config.discard_bad_only_first = parse_bool(value, key)?;
        }
        "uselowerdepth" | "uld" => config.use_lower_depth = parse_bool(value, key)?,
        "lbd" | "lowbindepth" | "lowerlimit" => config.low_bin_depth = parse_i64(value, key)?,
        "hbd" | "highbindepth" | "upperlimit" => config.high_bin_depth = parse_i64(value, key)?,
        "histlen" | "histogramlen" => config.hist_len = parse_usize(value, key)?.saturating_add(1),
        "histcol" | "histcolumns" | "histogramcolumns" => {
            config.hist_columns = parse_u8(value, key)?
        }
        "printzerocoverage" | "pzc" => config.print_zero_coverage = parse_bool(value, key)?,
        "minheight" | "h" => config.peak_min_height = parse_u64(value, key)?,
        "minvolume" | "v" => config.peak_min_volume = parse_u64(value, key)?,
        "minwidth" | "w" => config.peak_min_width = parse_usize(value, key)?,
        "minpeak" | "minp" => config.peak_min_peak = parse_usize(value, key)?,
        "maxpeak" | "maxp" => config.peak_max_peak = parse_usize(value, key)?,
        "ploidy" => config.peak_ploidy = parse_i32(value, key)?,
        "maxpeakcount" | "maxpc" | "maxpeaks" => {
            config.peak_max_count = parse_usize(value, key)?.max(1)
        }
        "overwrite" | "ow" => config.overwrite = parse_bool(value, key)?,
        "passes" | "p" => {
            config.passes = parse_usize(value, key)?;
        }
        "1pass" | "1p" => {
            config.passes = 1;
            config.notes.push("single-pass mode selected".to_string());
        }
        "2pass" | "2p" => {
            config.passes = 2;
        }
        "ascii" | "asciioffset" | "quality" | "qual" => {
            let offset = parse_quality_offset(value, key)?.unwrap_or(33);
            config.quality_in_offset = offset;
            config.quality_out_offset = offset;
        }
        "qin" | "asciiin" | "qualityin" | "qualin" => {
            config.quality_in_offset = parse_quality_offset(value, key)?.unwrap_or(33);
        }
        "qout" | "asciiout" | "qualityout" | "qualout" => {
            config.quality_out_offset = parse_quality_offset(value, key)?.unwrap_or(33);
        }
        "qauto" => config
            .notes
            .push("qauto accepted for BBTools-compatible quality alias handling".into()),
        key if matches!(
            quality_recal_base_key(key),
            "recalibrate" | "recalibratequality" | "recal"
        ) =>
        {
            if parse_bool(value, key)? {
                bail!(
                    "{key}={value} enables BBTools quality recalibration; Rust does not implement output-affecting recalibration yet"
                );
            }
            config.notes.push(format!(
                "{key}={value} keeps BBTools quality recalibration disabled in the supported Rust path"
            ));
        }
        key if is_quality_recal_bool_key(key) => {
            let _ = parse_java_bool(value);
            config.notes.push(format!(
                "{key}={value} is a BBTools quality-recalibration control; covered Rust output is unchanged"
            ));
        }
        key if quality_recal_base_key(key) == "observationcutoff" => {
            let _ = parse_kmg_i64(value, key)?;
            config.notes.push(format!(
                "{key}={value} is a BBTools quality-recalibration control; covered Rust output is unchanged"
            ));
        }
        key if matches!(
            quality_recal_base_key(key),
            "recalpasses" | "recalqmax" | "recalqmin"
        ) =>
        {
            let _ = parse_i32(value, key)?;
            config.notes.push(format!(
                "{key}={value} is a BBTools quality-recalibration control; covered Rust output is unchanged"
            ));
        }
        key if quality_recal_base_key(key) == "qmatrixmode" => {
            config.notes.push(format!(
                "{key}={value} is a BBTools quality-recalibration matrix mode; covered Rust output is unchanged"
            ));
        }
        "ignorebadquality" | "ibq" => {
            if parse_bool(value, key)? {
                config.change_quality = false;
            }
        }
        "changequality" | "cq" => config.change_quality = parse_bool(value, key)?,
        "mincalledquality" => {
            config.min_called_quality = parse_i32_clamped(value, key, 0, 93)? as u8
        }
        "maxcalledquality" => {
            config.max_called_quality = parse_i32_clamped(value, key, 1, 93)? as u8
        }
        "fakequality" | "qfake" => {
            config.fake_quality = parse_i32_clamped(value, key, 0, 93)? as u8
        }
        "fakefastaqual" | "fakefastaquality" | "ffq" => parse_fake_fasta_quality(config, value)?,
        "fastawrap" | "wrap" => config.fasta_wrap = parse_fasta_wrap(value, key)?,
        "trd" | "trc" | "trimreaddescription" | "trimreaddescriptions" => {
            let _ = parse_java_bool(value);
            config.notes.push(format!(
                "{key}={value} is accepted for KmerNormalize compatibility; covered FASTA/FASTQ read output keeps full headers like Java"
            ));
        }
        "trimrefdescription" | "trimrefdescriptions" | "trimrname" => {
            let _ = parse_java_bool(value);
            config.notes.push(format!(
                "{key}={value} is a BBTools reference-name trimming control; covered FASTA/FASTQ read output is unchanged"
            ));
        }
        "utot" => config.u_to_t = parse_bool(value, key)?,
        "tuc" | "touppercase" => config.to_upper_case = parse_bool(value, key)?,
        "lctn" | "lowercaseton" => config.lower_case_to_n = parse_bool(value, key)?,
        "dotdashxton" => config.dot_dash_x_to_n = parse_bool(value, key)?,
        "undefinedton" | "iupacton" | "itn" => config.iupac_to_n = parse_bool(value, key)?,
        "fixjunk" => {
            if parse_bool(value, key)? {
                config.junk_mode = JunkMode::Fix;
            } else if config.junk_mode == JunkMode::Fix {
                config.junk_mode = JunkMode::Crash;
            }
        }
        "ignorejunk" => {
            if parse_bool(value, key)? {
                config.junk_mode = JunkMode::Ignore;
            } else if config.junk_mode == JunkMode::Ignore {
                config.junk_mode = JunkMode::Crash;
            }
        }
        "flagjunk" => {
            if parse_bool(value, key)? {
                config.junk_mode = JunkMode::Flag;
            } else if config.junk_mode == JunkMode::Flag {
                config.junk_mode = JunkMode::Crash;
            }
        }
        "tossjunk" => {
            if parse_bool(value, key)? {
                config.junk_mode = JunkMode::Flag;
            }
        }
        "crashjunk" | "failjunk" => {
            if parse_bool(value, key)? {
                config.junk_mode = JunkMode::Crash;
            } else if config.junk_mode == JunkMode::Crash {
                config.junk_mode = JunkMode::Ignore;
            }
        }
        "junk" => parse_junk_mode(config, value)?,
        "threads" | "t" => {
            let threads = value.to_ascii_lowercase();
            if threads == "auto" {
                config
                    .notes
                    .push("threads=auto accepted; Rayon will use its default worker count".into());
            } else if matches!(threads.as_str(), "max" | "all") {
                let workers = std::thread::available_parallelism()
                    .map(|threads| threads.get())
                    .unwrap_or(1);
                config.threads = Some(workers);
                config.notes.push(format!(
                    "threads={threads} accepted; Rayon worker count will use all {workers} available workers"
                ));
            } else {
                let threads = parse_i64(value, key)?;
                if threads > 1 {
                    config.threads = Some(threads as usize);
                    config.notes.push(format!(
                        "threads={threads} accepted; Rayon worker count will be capped to {threads}"
                    ));
                } else if threads == 1 {
                    config.threads = Some(1);
                }
            }
        }
        "null" => {}
        "monitor" | "killswitch" => {
            parse_monitor(value, key)?;
            config.notes.push(format!(
                "{key}={value} is a BBTools watchdog runtime control; the Rust CLI accepts it as a no-op"
            ));
        }
        "outstream" | "proxyhost" | "proxyport" | "metadatafile" => {
            config.notes.push(format!(
                "{key}={value} is a BBTools preparser runtime control; covered Rust output records are unchanged"
            ));
        }
        "json" | "silent" | "printexecuting" | "bufferbf" | "bufferbf1" => {
            let _ = parse_java_bool(value);
            config.notes.push(format!(
                "{key}={value} is a BBTools preparser runtime control; covered Rust output records are unchanged"
            ));
        }
        "testsize" => {
            let _ = parse_java_bool(value);
            config.notes.push(format!(
                "{key}={value} is a BBTools diagnostic sizing control; covered Rust output records are unchanged"
            ));
        }
        "breaklen" | "breaklength" => {
            let break_len = parse_i32(value, key)?;
            if break_len > 0 {
                bail!(
                    "{key}={value} enables BBTools read breaking; Rust does not implement output-affecting read splitting yet"
                );
            }
            config.notes.push(format!(
                "{key}={value} keeps BBTools read breaking disabled in the supported Rust path"
            ));
        }
        "usejni" | "jni" | "skipvalidation" | "validate" | "validateinconstructor" | "vic" => {
            let _ = parse_java_bool(value);
            config.notes.push(format!(
                "{key}={value} is a BBTools shared runtime/validation control; covered Rust output is unchanged"
            ));
        }
        "usempi" | "mpi" => {
            let enabled = parse_mpi_enabled(value, key)?;
            config.notes.push(format!(
                "{key}={value} is a BBTools MPI execution control; Rust runs locally and ignores MPI mode{}",
                if enabled { " for ASAP output" } else { "" }
            ));
        }
        "crismpi" | "mpikeepall" => {
            let enabled = parse_java_bool(value);
            config.notes.push(format!(
                "{key}={value} is a BBTools MPI stream control; Rust runs locally and ignores MPI stream mode{}",
                if enabled { " for ASAP output" } else { "" }
            ));
        }
        "bf1" | "bytefile1" | "bf2" | "bytefile2" | "bf3" | "bytefile3" | "bf4" | "bytefile4" => {
            let _ = parse_java_bool(value);
            config.notes.push(format!(
                "{key}={value} is a BBTools byte-file runtime control; covered Rust output is unchanged"
            ));
        }
        "bf1bufferlen" | "readbufferlength" | "readbufferlen" | "readbufferdata" => {
            let _ = parse_kmg_i64(value, key)?;
            config.notes.push(format!(
                "{key}={value} is a BBTools buffer-sizing control; covered Rust output is unchanged"
            ));
        }
        "bf4threads" | "bfthreads" | "readbuffers" => {
            let _ = parse_i32(value, key)?;
            config.notes.push(format!(
                "{key}={value} is a BBTools I/O threading control; current Rust engine manages I/O internally"
            ));
        }
        "workers" | "workerthreads" | "wt" | "threadsin" | "tin" | "threadsout" | "tout" => {
            parse_auto_or_i32(value, key)?;
            config.notes.push(format!(
                "{key}={value} is a BBTools I/O worker control; current Rust engine manages I/O internally"
            ));
        }
        "zipthreads" | "bgzfthreadsin" | "bgzftin" | "bgzfreadthreads" | "bgzfthreadsout"
        | "bgzftout" | "bgzfwritethreads" => {
            let threads = parse_i32(value, key)?;
            if threads > 0 {
                config.gzip_threads = Some(threads as usize);
            }
            config.notes.push(format!(
                "{key}={value} is a BBTools compression/threading control; Rust uses gzip input/output worker settings for .gz files when threads > 1"
            ));
        }
        "ziplevel" | "zl" | "bziplevel" | "bzl" | "blocksize" | "pigziterations" | "pigziters" => {
            let _ = parse_i32(value, key)?;
            config.notes.push(format!(
                "{key}={value} is a BBTools compression/threading control; covered Rust output records are unchanged"
            ));
        }
        "zipthreaddivisor" | "ztd" => {
            let _ = parse_f64(value, key)?;
            config.notes.push(format!(
                "{key}={value} is a BBTools compression/threading control; covered Rust output records are unchanged"
            ));
        }
        "usebgzip" | "bgzip" | "usepigz" | "pigz" => {
            if value
                .as_bytes()
                .first()
                .is_some_and(|byte| byte.is_ascii_digit())
            {
                let threads = parse_i32(value, key)?;
                if threads > 0 {
                    config.gzip_threads = Some(threads as usize);
                }
            } else if parse_java_bool(value) {
                let workers = config
                    .threads
                    .unwrap_or_else(|| std::thread::available_parallelism().map_or(1, |n| n.get()));
                if workers > 1 {
                    config.gzip_threads = Some(workers);
                }
            } else {
                config.gzip_threads = Some(1);
            }
            config.notes.push(format!(
                "{key}={value} is a BBTools compression control; Rust uses zlib-rs gzip plus pigz/unpigz hooks for .gz input/output when enabled and available"
            ));
        }
        "usegunzip" | "gunzip" | "ungzip" | "useunpigz" | "unpigz" | "useunbgzip" | "unbgzip" => {
            if value
                .as_bytes()
                .first()
                .is_some_and(|byte| byte.is_ascii_digit())
            {
                let threads = parse_i32(value, key)?;
                if threads > 0 {
                    config.gzip_threads = Some(threads as usize);
                }
            } else if parse_java_bool(value) {
                let workers = config
                    .threads
                    .unwrap_or_else(|| std::thread::available_parallelism().map_or(1, |n| n.get()));
                if workers > 1 {
                    config.gzip_threads = Some(workers);
                }
            } else {
                config.gzip_threads = Some(1);
            }
            config.notes.push(format!(
                "{key}={value} is a BBTools gzip-input control; Rust uses zlib-rs and tries pigz/unpigz for .gz input when worker count is >1"
            ));
        }
        "allowziplevelchange"
        | "usegzip"
        | "gzip"
        | "usebgzf"
        | "bgzf"
        | "forcepigz"
        | "forcebgzip"
        | "preferbgzip"
        | "nativebgzip"
        | "nativebgzf"
        | "usenativebgzip"
        | "usenativebgzf"
        | "allownativebgzip"
        | "allownativebgzf"
        | "nativebgzipin"
        | "nativebgzfin"
        | "nativebgzipout"
        | "nativebgzfout"
        | "prefernativebgzip"
        | "prefernativebgzf"
        | "nativebgzipmt"
        | "nativebgzfmt"
        | "multithreadedbgzf"
        | "bgzfosmt2"
        | "filteredbgzf"
        | "preferunbgzip"
        | "usebzip2"
        | "bzip2"
        | "usepbzip2"
        | "pbzip2"
        | "uselbzip2"
        | "lbzip2" => {
            let _ = parse_java_bool(value);
            config.notes.push(format!(
                "{key}={value} is a BBTools compression/runtime control; covered Rust output records are unchanged"
            ));
        }
        "samversion" | "samv" | "sam" => {
            let _ = parse_f64(value, key)?;
            config.notes.push(format!(
                "{key}={value} is a BBTools SAM-version control; covered Rust FASTA/FASTQ output is unchanged"
            ));
        }
        "streamerthreads"
        | "ssthreads"
        | "bsthreads"
        | "fastqstreamerthreads"
        | "fqsthreads"
        | "fastastreamerthreads"
        | "fasthreads"
        | "samwriterthreads"
        | "swthreads"
        | "bamwriterthreads"
        | "bwthreads"
        | "fastqwriterthreads"
        | "fqwthreads"
        | "intronlen"
        | "intronlength" => {
            let _ = parse_i32(value, key)?;
            config.notes.push(format!(
                "{key}={value} is a BBTools SAM/streamer threading control; current Rust engine manages FASTA/FASTQ I/O internally"
            ));
        }
        "sambamba"
        | "samtools"
        | "printheaderwait"
        | "nativebam"
        | "usenativebam"
        | "allownativebam"
        | "nativebamout"
        | "usenativebamout"
        | "nativebamin"
        | "usenativebamin"
        | "prefernativebamout"
        | "prefernativebamin"
        | "prefernativebam"
        | "userssw"
        | "attachedsamline"
        | "useattachedsamline"
        | "fastastreamer2"
        | "prefermd"
        | "prefermdtag"
        | "notags"
        | "mdtag"
        | "md"
        | "idtag"
        | "mateqtag"
        | "xmtag"
        | "xm"
        | "smtag"
        | "amtag"
        | "nmtag"
        | "xttag"
        | "stoptag"
        | "lengthtag"
        | "boundstag"
        | "scoretag"
        | "sortscaffolds"
        | "customtag"
        | "nhtag"
        | "keepnames"
        | "saa"
        | "secondaryalignmentasterisks"
        | "inserttag"
        | "correctnesstag"
        | "suppressheader"
        | "noheader"
        | "noheadersequences"
        | "nhs"
        | "suppressheadersequences"
        | "tophat"
        | "flipsam" => {
            let _ = parse_java_bool(value);
            config.notes.push(format!(
                "{key}={value} is a BBTools SAM/BAM runtime control; covered Rust FASTA/FASTQ output is unchanged"
            ));
        }
        "xstag" | "xs" => {
            let lower = value.to_ascii_lowercase();
            if !matches!(
                lower.strip_prefix("fr-").unwrap_or(&lower),
                "ss" | "secondstrand" | "fs" | "firststrand" | "us" | "unstranded"
            ) {
                let _ = parse_java_bool(value);
            }
            config.notes.push(format!(
                "{key}={value} is a BBTools SAM XS-tag control; covered Rust FASTA/FASTQ output is unchanged"
            ));
        }
        "readgroup" | "readgroupid" | "rgid" | "readgroupcn" | "rgcn" | "readgroupds" | "rgds"
        | "readgroupdt" | "rgdt" | "readgroupfo" | "rgfo" | "readgroupks" | "rgks"
        | "readgrouplb" | "rglb" | "readgrouppg" | "rgpg" | "readgrouppi" | "rgpi"
        | "readgrouppl" | "rgpl" | "readgrouppu" | "rgpu" | "readgroupsm" | "rgsm" => {
            config.notes.push(format!(
                "{key}={value} is a BBTools read-group metadata control; covered Rust FASTA/FASTQ output is unchanged"
            ));
        }
        "tossbrokenreads"
        | "nullifybrokenquality"
        | "nbq"
        | "rbm"
        | "renamebymapping"
        | "don"
        | "deleteoldname"
        | "assertcigar"
        | "verbosesamline"
        | "parsecustom"
        | "fastqparsecustom"
        | "shrinkheaders"
        | "fixheader"
        | "fixheaders"
        | "allownullheader"
        | "allownullheaders"
        | "recalpairnum"
        | "recalibratepairnum" => {
            let _ = parse_java_bool(value);
            config.notes.push(format!(
                "{key}={value} is a BBTools shared read/header runtime control; covered Rust FASTA/FASTQ output is unchanged"
            ));
        }
        "pairreads" | "flipr2" => {
            let enabled = parse_java_bool(value);
            config.notes.push(format!(
                "{key}={value} is a BBTools global pairing behavior control; Rust pairing uses explicit in2=, interleaved=, and # routing{}",
                if enabled { " for ASAP output" } else { "" }
            ));
        }
        "aminoin" | "amino" | "amino8" => {
            if parse_java_bool(value) {
                bail!(
                    "{key}={value} enables BBTools amino-acid kmer mode; the Rust engine currently supports nucleotide BBNorm only"
                );
            }
            config.notes.push(format!(
                "{key}={value} keeps BBTools amino-acid kmer mode disabled in the supported Rust path"
            ));
        }
        "validatebranchless"
        | "fairqueues"
        | "fixextensions"
        | "fixextension"
        | "tryallextensions"
        | "2passresize"
        | "twopassresize"
        | "parallelsort"
        | "paralellsort"
        | "gcbeforemem"
        | "warnifnosequence"
        | "warnfirsttimeonly"
        | "kmg"
        | "outputkmg"
        | "forcejavaparsedouble"
        | "simdsparse"
        | "simdmultsparse"
        | "simdfmasparse"
        | "simdcopy"
        | "awsservers"
        | "aws"
        | "nerscservers"
        | "nersc"
        | "lowmem"
        | "lowram"
        | "lowmemory"
        | "buffer"
        | "buffered"
        | "sidechannelstats"
        | "comment"
        | "taxpath"
        | "silva"
        | "unite"
        | "imghq"
        | "callins"
        | "callinss"
        | "calldel"
        | "calldels"
        | "callsub"
        | "callsubs"
        | "callsnp"
        | "callsnps"
        | "callindel"
        | "callindels"
        | "calljunct"
        | "calljunction"
        | "calljunctions"
        | "callnocall"
        | "callnocalls"
        | "protfull" => {
            let _ = parse_java_bool(value);
            config.notes.push(format!(
                "{key}={value} is a BBTools shared environment/performance control; covered Rust output is unchanged"
            ));
        }
        "lockedincrement" | "symmetricwrite" | "symmetric" | "sw" => {
            if value.eq_ignore_ascii_case("auto") {
                config.locked_increment = None;
            } else {
                config.locked_increment = Some(parse_java_bool(value));
            }
            config.notes.push(format!(
                "{key}={value} is a BBTools KCountArray write-symmetry control; bounded Rust sketches use the matching locked/conservative update mode when applicable"
            ));
        }
        "simd" => {
            if !value.eq_ignore_ascii_case("auto") {
                let _ = parse_java_bool(value);
            }
            config.notes.push(format!(
                "{key}={value} is a BBTools SIMD runtime control; covered Rust output is unchanged"
            ));
        }
        "entropyk" | "ek" | "entropywindow" | "ew" => {
            let parsed = parse_i32(value, key)?;
            if parsed <= 0 {
                bail!("{key} expects a positive integer, got {value}");
            }
            if matches!(key, "entropyk" | "ek") {
                config.entropy_k = parsed as usize;
            } else {
                config.entropy_window = parsed as usize;
            }
            config.notes.push(format!(
                "{key}={value} is a BBTools entropy-stat runtime control; Rust applies it to emitted entropy histograms"
            ));
        }
        "barcodestats" | "barcodecounts" | "timehistogram" | "thist" => {
            config.notes.push(format!(
                "{key}={value} is a BBTools side-output stats histogram; Rust does not emit this auxiliary file yet and keeps the supported normalization path"
            ));
        }
        "matchhistogram" | "matchhist" | "mhist" => {
            config.match_hist_out = Some(path(value));
            config.notes.push(format!(
                "{key}={value} is a BBTools side-output match histogram; Rust emits a covered no-alignment sequence-match fallback histogram"
            ));
        }
        "inserthistogram" | "inserthist" | "ihist" => {
            config.insert_hist_out = Some(path(value));
            config.notes.push(format!(
                "{key}={value} is a BBTools side-output insert histogram; Rust emits a covered no-alignment insert-size fallback histogram"
            ));
        }
        "qualityaccuracyhistogram" | "qahist" => {
            config.quality_accuracy_hist_out = Some(path(value));
            config.notes.push(format!(
                "{key}={value} is a BBTools side-output quality-accuracy histogram; Rust emits a covered no-alignment quality-accuracy fallback histogram"
            ));
        }
        "indelhistogram" | "indelhist" => {
            config.indel_hist_out = Some(path(value));
            config.notes.push(format!(
                "{key}={value} is a BBTools side-output indel histogram; Rust emits a covered no-alignment indel fallback histogram"
            ));
        }
        "errorhistogram" | "ehist" => {
            config.error_hist_out = Some(path(value));
            config.notes.push(format!(
                "{key}={value} is a BBTools side-output error histogram; Rust emits a covered no-alignment error-count fallback histogram"
            ));
        }
        "gchistogram" | "gchist" => {
            config.gc_hist_out = Some(path(value));
            config.notes.push(format!(
                "{key}={value} is a BBTools side-output GC histogram; Rust emits a covered primary input GC-bin histogram"
            ));
        }
        "qualityhistogram" | "qualityhist" | "qhist" => {
            config.quality_hist_out = Some(path(value));
            config.notes.push(format!(
                "{key}={value} is a BBTools side-output quality histogram; Rust emits a covered primary input quality histogram"
            ));
        }
        "basequalityhistogram" | "basequalityhist" | "bqhist" => {
            config.base_quality_hist_out = Some(path(value));
            config.notes.push(format!(
                "{key}={value} is a BBTools side-output base-quality histogram; Rust emits a covered primary input base-quality histogram"
            ));
        }
        "qualitycounthistogram" | "qualitycounthist" | "qchist" | "qdhist" | "qfhist" => {
            config.quality_count_hist_out = Some(path(value));
            config.notes.push(format!(
                "{key}={value} is a BBTools side-output quality-count histogram; Rust emits a covered primary input quality-count histogram"
            ));
        }
        "averagequalityhistogram" | "aqhist" => {
            config.average_quality_hist_out = Some(path(value));
            config.notes.push(format!(
                "{key}={value} is a BBTools side-output average-quality histogram; Rust emits a covered primary input average-quality histogram"
            ));
        }
        "overallbasequalityhistogram" | "overallbasequalityhist" | "obqhist" => {
            config.overall_base_quality_hist_out = Some(path(value));
            config.notes.push(format!(
                "{key}={value} is a BBTools side-output overall base-quality histogram; Rust emits a covered primary input overall base-quality histogram"
            ));
        }
        "lengthhistogram" | "lhist" => {
            config.length_hist_out = Some(path(value));
            config.notes.push(format!(
                "{key}={value} is a BBTools side-output length histogram; Rust emits a covered read-length histogram for the primary input"
            ));
        }
        "basehistogram" | "basehist" | "bhist" => {
            config.base_hist_out = Some(path(value));
            config.notes.push(format!(
                "{key}={value} is a BBTools side-output base-content histogram; Rust emits a covered primary input base-content histogram"
            ));
        }
        "entropyhistogram" | "entropyhist" | "enhist" | "enthist" => {
            config.entropy_hist_out = Some(path(value));
            config.notes.push(format!(
                "{key}={value} is a BBTools side-output entropy histogram; Rust emits a covered primary input entropy histogram"
            ));
        }
        "identityhistogram" | "idhist" => {
            config.identity_hist_out = Some(path(value));
            config.notes.push(format!(
                "{key}={value} is a BBTools side-output identity histogram; Rust emits a covered sequence-input identity fallback histogram because this BBNorm path has no aligner"
            ));
        }
        "gcbins" | "gchistbins" => {
            if !value.eq_ignore_ascii_case("auto") {
                let bins = parse_i32(value, key)?;
                if bins <= 0 {
                    bail!("{key} expects a positive integer or auto, got {value}");
                }
                config.gc_bins = Some(bins as usize);
            }
            config.notes.push(format!(
                "{key}={value} is a BBTools side-output GC histogram sizing control; Rust applies it to emitted GC histograms"
            ));
        }
        "entropybins" | "entropyhistbins" | "entbins" | "enthistbins" => {
            if !value.eq_ignore_ascii_case("auto") {
                let bins = parse_i32(value, key)?;
                if bins <= 0 {
                    bail!("{key} expects a positive integer or auto, got {value}");
                }
                config.entropy_bins = bins as usize;
            }
            config.notes.push(format!(
                "{key}={value} is a BBTools side-output entropy histogram sizing control; Rust applies it to emitted entropy histograms"
            ));
        }
        "idhistlen" | "idhistlength" | "idhistbins" | "idbins" => {
            if !value.eq_ignore_ascii_case("auto") {
                let bins = parse_i32(value, key)?;
                if bins <= 0 {
                    bail!("{key} expects a positive integer or auto, got {value}");
                }
                config.identity_bins = bins as usize;
            }
            config.notes.push(format!(
                "{key}={value} is a BBTools side-output identity histogram sizing control; Rust applies it to emitted identity histograms"
            ));
        }
        "entropyns" | "entropyhistns" => {
            config.allow_entropy_ns = parse_java_bool(value);
            config.notes.push(format!(
                "{key}={value} is a BBTools side-output entropy control; Rust applies it to emitted entropy histograms"
            ));
        }
        "gcchart" | "gcplot" | "fixindels" | "ignorevcfindels" => {
            let _ = parse_java_bool(value);
            config.notes.push(format!(
                "{key}={value} is a BBTools side-output stats control; covered Rust FASTA/FASTQ output is unchanged"
            ));
        }
        "maxhistlen" => {
            let len = parse_kmg_i64(value, key)?;
            if len <= 0 {
                bail!("{key} expects a positive KMG value, got {value}");
            }
            config.side_hist_len = Some(
                usize::try_from(len)
                    .map_err(|_| anyhow::anyhow!("{key} value is out of range: {value}"))?,
            );
            config.notes.push(format!(
                "{key}={value} is a BBTools side-output histogram length control; Rust applies it to emitted side histograms"
            ));
        }
        "cardinality" | "loglog" => {
            match parse_cardinality_bool_or_int(value, key)? {
                CardinalityToggle::Bool(enabled) => config.cardinality.input = enabled,
                CardinalityToggle::Int(k) => {
                    config.cardinality.input = true;
                    config.cardinality.k = Some(k);
                }
            }
            config.notes.push(format!(
                "{key}={value} is a BBTools cardinality/loglog control; Rust emits a bounded input estimate when enabled"
            ));
        }
        "loglogin" => {
            match parse_cardinality_bool_or_int(value, key)? {
                CardinalityToggle::Bool(enabled) => config.cardinality.input = enabled,
                CardinalityToggle::Int(k) => {
                    config.cardinality.input = true;
                    config.cardinality.k = Some(k);
                }
            }
            config.notes.push(format!(
                "{key}={value} is a BBTools cardinality/loglog input control; Rust emits a bounded input estimate when enabled"
            ));
        }
        "cardinalityout" | "loglogout" => {
            match parse_cardinality_bool_or_int(value, key)? {
                CardinalityToggle::Bool(enabled) => config.cardinality.output = enabled,
                CardinalityToggle::Int(k) => {
                    config.cardinality.output = true;
                    config.cardinality.k = Some(k);
                }
            }
            config.notes.push(format!(
                "{key}={value} is a BBTools cardinality/loglog control; Rust emits a bounded output estimate when enabled"
            ));
        }
        "buckets" | "loglogbuckets" => {
            let buckets = parse_cardinality_buckets(value, key)?;
            config.cardinality.buckets = buckets;
            config.notes.push(format!(
                "{key}={value} is a BBTools cardinality/loglog bucket control; Rust applies it to bounded cardinality estimates"
            ));
        }
        "loglogk" | "cardinalityk" | "kcardinality" => {
            config.cardinality.k = Some(parse_cardinality_k(value, key)?);
            config.notes.push(format!(
                "{key}={value} is a BBTools cardinality/loglog numeric control; Rust applies it to bounded cardinality estimates"
            ));
        }
        "loglogbits" | "loglogmantissa" => {
            let _ = parse_i32(value, key)?;
            config.notes.push(format!(
                "{key}={value} is a BBTools cardinality/loglog numeric control; Rust accepts it while using compact byte registers"
            ));
        }
        "loglogklist" => {
            let mut first_k = None;
            for part in value.split(',') {
                let trimmed = part.trim();
                if trimmed.is_empty() {
                    bail!("{key} expects a comma-separated integer list, got {value}");
                }
                let k = parse_cardinality_k(trimmed, key)?;
                first_k.get_or_insert(k);
            }
            config.cardinality.k = first_k;
            config.notes.push(format!(
                "{key}={value} is a BBTools cardinality/loglog k-list; Rust uses the first k for bounded cardinality estimates"
            ));
        }
        "loglogseed" => {
            config.cardinality.seed = parse_cardinality_seed(value, key)?;
            config.notes.push(format!(
                "{key}={value} is a BBTools cardinality/loglog seed; Rust applies it to bounded cardinality estimates"
            ));
        }
        "loglogminprob" => {
            let min_probability = parse_f64(value, key)?;
            if !(0.0..=1.0).contains(&min_probability) {
                bail!("{key} expects a probability between 0 and 1, got {value}");
            }
            config.cardinality.min_probability = min_probability;
            config.notes.push(format!(
                "{key}={value} is a BBTools cardinality/loglog probability threshold; Rust records it for bounded cardinality estimates"
            ));
        }
        "loglogtype" => {
            config.notes.push(format!(
                "{key}={value} is a BBTools cardinality/loglog estimator type; Rust uses its compact bounded estimator"
            ));
        }
        "loglogcorrection" | "loglogcf" | "loglogmean" | "loglogmedian" | "loglogmwa"
        | "logloghmean" | "logloggmean" | "loglogcounts" | "loglogcount" => {
            let _ = parse_java_bool(value);
            config.notes.push(format!(
                "{key}={value} is a BBTools cardinality/loglog output-control toggle; Rust emits compact summary estimates"
            ));
        }
        "countup" => {
            config.count_up = parse_bool(value, key)?;
            if !config.count_up {
                config.notes.push(
                    "countup=f selected; standard single-pass normalization remains active"
                        .to_string(),
                );
            }
        }
        "bits" | "cbits" | "cellbits" => {
            let bits = parse_kcount_cell_bits(value, key)?;
            config.count_min.bits = Some(bits);
            config.notes.push(format!(
                "{key}={bits} is a BBTools count-min cell-width control; constrained Rust count-min tables use it for saturation"
            ));
        }
        "bits1" | "cbits1" | "cellbits1" => {
            let bits = parse_kcount_cell_bits(value, key)?;
            config.count_min_bits_first = Some(bits);
            config.notes.push(format!(
                "{key}={bits} is a BBTools first/intermediate-pass sketch width control; Rust uses it for multipass bounded sketches"
            ));
        }
        "hashes" => {
            let hashes = parse_kcount_hashes(value, key)?;
            config.count_min.hashes = Some(hashes);
            config.notes.push(format!(
                "hashes={hashes} is a BBTools count-min hashing control; constrained Rust count-min tables use it for collision estimates"
            ));
        }
        "cells" | "matrixbits" => {
            let cells = if key == "matrixbits" {
                parse_matrixbits_cells(value, key)?
            } else {
                parse_positive_kmg_usize(value, key)?
            };
            config.count_min.cells = Some(cells.max(1));
            config.notes.push(format!(
                "{key}={value} is a BBTools count-min table-sizing control; Rust treats it as a total-cell budget and builds a fixed-memory count-min input sketch"
            ));
        }
        "sketchmemory" | "sketchmem" | "countminmemory" | "countminmem" | "cmem" => {
            let bytes = parse_positive_kmg_usize(value, key)?;
            config.count_min.memory_bytes = Some(bytes);
            config.notes.push(format!(
                "{key}={value} is a Rust count-min memory budget; Rust sizes the fixed-memory input sketch from this budget when cells/matrixbits are not set"
            ));
        }
        "maxcountupspillbytes"
        | "maxcountupspilllivebytes"
        | "countupspillbytes"
        | "countupspilllimit" => {
            let bytes = parse_kmg_usize(value, key)?;
            config.max_countup_spill_live_bytes = Some(bytes as u64);
            config.notes.push(format!(
                "{key}={value} is a Rust count-up temp-spill safety cap; Rust aborts count-up if peak live spill bytes exceed {bytes}"
            ));
        }
        "maxcountupspillfinallivebytes"
        | "maxcountupspillfinalbytes"
        | "countupspillfinallivebytes" => {
            let bytes = parse_kmg_usize(value, key)?;
            config.max_countup_spill_final_live_bytes = Some(bytes as u64);
            config.notes.push(format!(
                "{key}={value} is a Rust count-up temp-spill safety cap; Rust aborts count-up if current/final live spill bytes exceed {bytes}"
            ));
        }
        "maxcountupspillinitialruns" | "countupspillinitialruns" => {
            let runs = parse_kmg_usize(value, key)?;
            config.max_countup_spill_initial_runs = Some(runs);
            config.notes.push(format!(
                "{key}={value} is a Rust count-up temp-spill safety cap; Rust aborts count-up if initial spill run count exceeds {runs}"
            ));
        }
        "maxcountupspillmergeruns" | "countupspillmergeruns" => {
            let runs = parse_kmg_usize(value, key)?;
            config.max_countup_spill_merge_runs = Some(runs);
            config.notes.push(format!(
                "{key}={value} is a Rust count-up temp-spill safety cap; Rust aborts count-up if merge spill run count exceeds {runs}"
            ));
        }
        "maxcountupspillfinalruns" | "maxcountupspillruns" | "countupspillfinalruns" => {
            let runs = parse_kmg_usize(value, key)?;
            config.max_countup_spill_final_runs = Some(runs);
            config.notes.push(format!(
                "{key}={value} is a Rust count-up temp-spill safety cap; Rust aborts count-up if live/final spill run count exceeds {runs}"
            ));
        }
        "maxcountupspillwritebytes" | "maxcountupspillwrittenbytes" | "countupspillwritebytes" => {
            let bytes = parse_kmg_usize(value, key)?;
            config.max_countup_spill_write_bytes = Some(bytes as u64);
            config.notes.push(format!(
                "{key}={value} is a Rust count-up temp-spill I/O safety cap; Rust aborts count-up if cumulative spill bytes written exceed {bytes}"
            ));
        }
        "memory" | "mem" | "ram" | "maxmemory" | "maxmem" | "xmx" => {
            let bytes = parse_positive_kmg_usize(value, key)?;
            config.auto_count_min_memory_bytes = Some(bytes);
            config.auto_count_min = true;
            config.notes.push(format!(
                "{key}={value} is a BBTools-style memory budget; automatic Rust count-min sizing uses it for large inputs"
            ));
        }
        "autocountmin" | "autosketch" | "autosketchtable" | "autosketchtables" => {
            config.auto_count_min = parse_bool(value, key)?;
            config.notes.push(format!(
                "{key}={} controls Rust's large-input automatic bounded count-min table selection",
                config.auto_count_min
            ));
        }
        "exact" | "exactcount" | "exactcounts" | "useexact" | "sketchexact" => {
            config.force_exact_counts = parse_bool(value, key)?;
            config.notes.push(format!(
                "{key}={} forces Rust exact-count maps and disables automatic/explicit count-min sketches",
                config.force_exact_counts
            ));
        }
        "autosketchbytes" | "autosketchminbytes" | "autocountminbytes" | "autocountminminbytes" => {
            config.auto_count_min_input_bytes = parse_positive_kmg_usize(value, key)?;
            config.notes.push(format!(
                "{key}={value} sets the compressed/uncompressed input-size trigger for automatic Rust count-min tables"
            ));
        }
        "autosketchtablereads" | "autocountminreads" | "autosketchtablereadthreshold" => {
            config.auto_count_min_read_threshold = parse_u64(value, key)?.max(1);
            config.notes.push(format!(
                "{key}={value} sets the read-limit trigger for automatic Rust count-min tables"
            ));
        }
        "precells" | "prefiltercells" => {
            let cells = parse_kmg_usize(value, key)?;
            config.prefilter.cells = (cells > 0).then_some(cells);
            if cells == 0 {
                config.notes.push(format!(
                    "{key}=0 is a BBTools prefilter sketch control; Rust leaves prefilter cells unset unless prefiltering is otherwise requested"
                ));
            } else {
                config.prefilter.enabled = true;
                config.prefilter.force_disabled = false;
                config.notes.push(format!(
                    "{key}={value} is a BBTools prefilter sketch control; Rust applies deterministic prefilter collision estimates when prefilter cells are constrained"
                ));
            }
        }
        "prefiltersize" | "prefilterfraction" => {
            let fraction = parse_fraction_micros(value, key)?;
            config.prefilter.memory_fraction_micros = (fraction > 0).then_some(fraction);
            config.prefilter.enabled = fraction > 0;
            config.prefilter.force_disabled = fraction == 0;
            if fraction == 0 {
                config.notes.push(format!(
                    "{key}=0 is a BBTools prefilter sketch control; Rust disables fraction-derived prefilter sizing"
                ));
            } else {
                config.notes.push(format!(
                    "{key}={value} is a BBTools prefilter sketch control; Rust derives deterministic prefilter collision memory from the configured table memory budget"
                ));
            }
        }
        "prefilterbits" | "prebits" | "pbits" => {
            let bits = parse_kcount_cell_bits(value, key)?;
            config.prefilter.bits = Some(bits);
            config.notes.push(format!(
                "{key}={value} is a BBTools prefilter sketch control; Rust uses it with constrained prefilter cells"
            ));
        }
        "prehashes" | "prefilterhashes" => {
            let hashes = parse_prefilter_hashes(value, key)?;
            config.prefilter.hashes = (hashes > 0).then_some(hashes);
            if hashes == 0 {
                config.notes.push(format!(
                    "{key}=0 is a BBTools prefilter sketch control; Rust leaves prefilter hashes unset unless prefiltering is otherwise requested"
                ));
            } else {
                config.prefilter.enabled = true;
                config.prefilter.force_disabled = false;
                config.notes.push(format!(
                    "{key}={value} is a BBTools prefilter sketch control; Rust applies deterministic prefilter collision estimates with explicit or implicit prefilter cells"
                ));
            }
        }
        "buildpasses" => {
            let build_passes = parse_i64(value, key)?;
            if build_passes <= 0 {
                bail!("{key} expects a positive integer, got {value}");
            }
            config.build_passes = usize::try_from(build_passes)
                .map_err(|_| anyhow::anyhow!("{key} value is out of range: {value}"))?;
            config.notes.push(format!(
                "{key}={build_passes} is a BBTools table-construction pass control; Rust applies deterministic trusted-kmer filtering when buildpasses is greater than 1"
            ));
        }
        "initialsize" => {
            let initial_size = parse_positive_kmg_usize(value, key)?;
            config.table_initial_size = Some(initial_size);
            config.notes.push(format!(
                "{key}={value} is a BBTools kmer-table runtime sizing control; Rust pre-reserves exact-count table capacity when practical"
            ));
        }
        "ways" => {
            let _ = parse_kmg_i64(value, key)?;
            config.notes.push(format!(
                "{key}={value} is a BBTools kmer-table runtime sizing control; exact Rust counting keeps native map sharding"
            ));
        }
        "buflen" | "bufflen" | "bufferlength" => {
            let _ = parse_kmg_i64(value, key)?;
            config.notes.push(format!(
                "{key}={value} is a BBTools kmer-table buffer-length control; covered Rust output records are unchanged"
            ));
        }
        "tabletype" => {
            let _ = parse_i32(value, key)?;
            config.notes.push(format!(
                "{key}={value} is a BBTools kmer-table implementation control; exact Rust counting uses its native map"
            ));
        }
        "rcomp" | "maskmiddle" => {
            let _ = parse_java_bool(value);
            config.notes.push(format!(
                "{key}={value} is a BBTools kmer-table matching control; covered Rust BBNorm canonical/exact-count behavior is unchanged"
            ));
        }
        "showstats" | "stats" | "showspeed" | "ss" | "verbose2" => {
            let _ = parse_java_bool(value);
            config.notes.push(format!(
                "{key}={value} is a BBTools kmer-table reporting control; covered Rust output records are unchanged"
            ));
        }
        "prealloc" | "preallocate" => {
            config.table_prealloc_fraction = parse_preallocation_fraction(value, key)?;
            config.notes.push(format!(
                "{key}={value} is a BBTools kmer-table preallocation control; Rust pre-reserves exact-count table capacity when practical"
            ));
        }
        "filtermemory" | "prefiltermemory" | "filtermem" | "filtermemoryoverride" => {
            let bytes = parse_positive_kmg_usize(value, key)?;
            config.prefilter.memory_bytes = Some(bytes);
            config.prefilter.enabled = true;
            config.prefilter.force_disabled = false;
            config.notes.push(format!(
                "{key}={value} is a BBTools prefilter memory-sizing control; Rust sizes deterministic prefilter collision estimates from this budget when prefilter cells are not set"
            ));
        }
        "minprobprefilter" | "mpp" | "minprobmain" | "mpm" => {
            let _ = parse_java_bool(value);
            config.notes.push(format!(
                "{key}={value} is a BBTools kmer-table minprob routing control; covered Rust minprob behavior is unchanged"
            ));
        }
        "prefilterpasses" | "prepasses" => {
            parse_auto_or_kmg_i64(value, key)?;
            config.notes.push(format!(
                "{key}={value} is a BBTools prefilter pass-count control; exact Rust counting uses one deterministic table build"
            ));
        }
        "onepass" => {
            let _ = parse_java_bool(value);
            config.notes.push(format!(
                "{key}={value} is a BBTools kmer-table construction-mode control; covered Rust output remains single-pass"
            ));
        }
        "stepsize" | "buildstepsize" => {
            let _ = parse_i32(value, key)?;
            config.notes.push(format!(
                "{key}={value} is a BBTools trusted-kmer sampling control; the covered no-ECC single-pass path ignores it"
            ));
        }
        "prefilter" => {
            config.prefilter.enabled = parse_bool(value, key)?;
            if config.prefilter.enabled {
                config.prefilter.force_disabled = false;
                config.notes.push(
                    "prefilter=t requested; Rust applies BBTools-style default prefilter partitioning when bounded count-min counting is selected"
                        .to_string(),
                );
            } else {
                config.prefilter.force_disabled = true;
                config.notes.push(
                    "prefilter=f requested; Rust disables prefilter sketch construction unless a later prefilter control re-enables it"
                        .to_string(),
                );
            }
        }
        "auto" | "automatic" => {
            let enabled = parse_bool(value, key)?;
            config.auto_count_min = enabled;
            config.notes.push(format!(
                "{key}={enabled} is a BBTools automatic count-table sizing control; Rust uses it to select bounded count-min tables for large inputs"
            ));
        }
        "tmpdir" => {
            config.temp_dir = Some(PathBuf::from(value));
            config.use_temp_dir = true;
            config.notes.push(format!(
                "{key}={value} is a BBTools temporary-directory control; covered Rust multipass and stdin paths use managed temp files there when enabled"
            ));
        }
        "usetmpdir" | "usetempdir" => {
            config.use_temp_dir = parse_java_bool(value);
            config.notes.push(format!(
                "{key}={value} is a BBTools temporary-directory control; covered Rust multipass and stdin paths use managed temp files there when enabled"
            ));
        }
        "ordered" | "ord" | "verbose" | "printcoverage" => {
            config.notes.push(format!(
                "{key}={value} is accepted as a no-op in this Rust parity slice"
            ));
        }
        "append" | "app" => {
            config.append = parse_bool(value, key)?;
        }
        "interleaved" | "int" => {
            let lower = value.to_ascii_lowercase();
            if lower == "auto" {
                config.interleaved = false;
                config.test_interleaved = true;
            } else {
                config.interleaved = parse_bool(value, key)?;
                config.test_interleaved = false;
            }
        }
        "testinterleaved" => {
            config.test_interleaved = parse_bool(value, key)?;
        }
        "forceinterleaved" => {
            config.interleaved = parse_bool(value, key)?;
            config.test_interleaved = false;
        }
        "overrideinterleaved" => {
            let _ = parse_bool(value, key)?;
            config.notes.push(format!(
                "{key}={value} is a BBTools paired-output assertion override; covered Rust paired output is unchanged"
            ));
        }
        "fastareadlen" | "fastareadlength" => {
            if parse_u64(value, key)? != u64::MAX && value != "2147483647" {
                config.notes.push(
                    "fastareadlen is accepted for KmerNormalize parity; covered FASTA records are processed as-is".to_string(),
                );
            }
        }
        "fastaminread" | "fastaminlen" | "fastaminlength" => {
            let _ = parse_i32(value, key)?;
            config.notes.push(format!(
                "{key}={value} is a BBTools FASTA parser control; covered KmerNormalize FASTA records are processed as-is"
            ));
        }
        "forcesectionname" | "fastadump" => {
            let _ = parse_java_bool(value);
            config.notes.push(format!(
                "{key}={value} is a BBTools FASTA parser control; covered Rust output is unchanged"
            ));
        }
        "sampleoutput" | "readsample" | "kmersample" => {
            config.notes.push(format!(
                "{key}={value} is advertised in bbnorm.sh but rejected by vendored KmerNormalize; Rust ignores it and keeps the supported normalization path"
            ));
        }
        "samplerate" | "sample" | "sampleseed" | "seed" => {
            config.notes.push(format!(
                "{key}={value} is a BBTools stream-wrapper sampling option; Rust ignores it and keeps the supported normalization path"
            ));
        }
        "markerrors" | "markonly" | "meo" => {
            config.mark_errors_only = parse_bool(value, key)?;
            if config.mark_errors_only {
                enable_error_correction_if_unset(config);
            }
        }
        "markuncorrectableerrors" | "markuncorrectable" | "mue" => {
            config.mark_uncorrectable_errors = parse_bool(value, key)?;
        }
        "tam" | "trimaftermarking" => {
            config.trim_after_marking = parse_bool(value, key)?;
        }
        "markwith1" | "markwithone" | "mw1" => {
            config.mark_with_one = parse_bool(value, key)?;
        }
        "aec" | "aecc" | "aggressiveerrorcorrection" => {
            let enabled = parse_bool(value, key)?;
            if enabled {
                config.error_correct = true;
                config.error_correct_first = true;
                config.error_correct_final = true;
                config.error_correct_high_thresh = config.error_correct_high_thresh.min(16);
                config.error_correct_low_thresh = config.error_correct_low_thresh.max(3);
                config.error_correct_ratio = config.error_correct_ratio.min(100);
                config.max_errors_to_correct = config.max_errors_to_correct.max(7);
                config.suffix_len = config.suffix_len.min(3);
                config.prefix_len = config.prefix_len.min(2);
            }
        }
        "cec" | "cecc" | "conservativeerrorcorrection" => {
            let enabled = parse_bool(value, key)?;
            if enabled {
                config.error_correct = true;
                config.error_correct_first = true;
                config.error_correct_final = true;
                config.error_correct_high_thresh = config.error_correct_high_thresh.max(30);
                config.error_correct_low_thresh = config.error_correct_low_thresh.min(1);
                config.error_correct_ratio = config.error_correct_ratio.max(170);
                config.max_errors_to_correct = config.max_errors_to_correct.min(2);
                config.max_quality_to_correct = config.max_quality_to_correct.min(25);
                config.suffix_len = config.suffix_len.max(4);
                config.prefix_len = config.prefix_len.max(4);
            }
        }
        "ecc" => {
            let enabled = parse_bool(value, key)?;
            config.error_correct = enabled;
            config.error_correct_first = enabled;
            config.error_correct_final = enabled;
            config.overlap_error_correct &= enabled;
            config.overlap_error_correct_auto &= enabled;
        }
        "ecc1" => {
            config.error_correct_first = parse_bool(value, key)?;
            config.error_correct = config.error_correct_first || config.error_correct_final;
        }
        "ecc2" | "eccf" => {
            config.error_correct_final = parse_bool(value, key)?;
            config.error_correct = config.error_correct_first || config.error_correct_final;
        }
        "eccbyoverlap" | "ecco" | "overlap" => {
            if value.eq_ignore_ascii_case("auto") {
                config.notes.push(format!(
                    "{key}=auto requests automatic overlap-based error correction; Rust samples paired reads and enables paired overlap repair when the overlap fraction is high"
                ));
                config.error_correct = true;
                config.error_correct_first = true;
                config.error_correct_final = true;
                config.overlap_error_correct = false;
                config.overlap_error_correct_auto = true;
            } else if parse_bool(value, key)? {
                config.notes.push(format!(
                    "{key}={value} requests overlap-based error correction; Rust uses paired overlap repair before the table-based ECC path"
                ));
                config.error_correct = true;
                config.error_correct_first = true;
                config.error_correct_final = true;
                config.overlap_error_correct = true;
                config.overlap_error_correct_auto = false;
            } else {
                config.overlap_error_correct = false;
                config.overlap_error_correct_auto = false;
            }
        }
        "ecclimit" => config.max_errors_to_correct = parse_usize(value, key)?,
        "eccmaxqual" => config.max_quality_to_correct = parse_u8(value, key)?,
        "errorcorrectratio" | "ecr" => config.error_correct_ratio = parse_u64(value, key)?,
        "echighthresh" | "echthresh" | "echt" => {
            config.error_correct_high_thresh = parse_u64(value, key)?
        }
        "eclowthresh" | "eclthresh" | "eclt" => {
            config.error_correct_low_thresh = parse_u64(value, key)?
        }
        "sl" | "suflen" | "suffixlen" => config.suffix_len = parse_usize(value, key)?,
        "pl" | "prelen" | "prefixlen" => config.prefix_len = parse_usize(value, key)?,
        "cfl" => config.correct_from_left = parse_bool(value, key)?,
        "cfr" => config.correct_from_right = parse_bool(value, key)?,
        "target1" | "targetdepth1" | "tgt1" => {
            config.target_depth_first = Some(parse_u64(value, key)?);
        }
        "targetbadpercentilelow" | "tbpl" => {
            let value = parse_percent(value, key)?;
            config.target_bad_percent_low = value;
            config.target_bad_percent_high = config.target_bad_percent_high.max(value);
        }
        "targetbadpercentilehigh" | "tbph" => {
            let value = parse_percent(value, key)?;
            config.target_bad_percent_high = value;
            config.target_bad_percent_low = config.target_bad_percent_low.min(value);
        }
        "abrc" | "addbadreadscountup" => {
            config.add_bad_reads_countup = parse_bool(value, key)?;
        }
        _ => bail!("unknown or unsupported BBNorm option: {key}={value}"),
    }
    Ok(())
}

fn validate(config: &mut Config) -> Result<()> {
    if config.in1.is_none() {
        bail!("missing input: provide in=<reads.fq>");
    }
    if !(1..=4).contains(&config.passes) {
        bail!("passes should be in range 1 through 4");
    }
    expand_hash_paired_input(config);
    validate_extra_inputs(config)?;
    if config.k == 0 {
        bail!("k must be greater than zero");
    }
    if !(0.0..1.0).contains(&config.min_prob) && (config.min_prob - 1.0).abs() > f64::EPSILON {
        bail!("minprob must be between 0 and 1");
    }
    if config.target_depth == 0 {
        bail!("target depth must be greater than zero");
    }
    if config.passes == 1 {
        config.target_bad_percent_low = 1.0;
        config.target_bad_percent_high = 1.0;
    }
    config.max_depth = Some(
        config
            .max_depth
            .unwrap_or(config.target_depth)
            .max(config.target_depth),
    );
    if config.error_detect_ratio == 0 {
        bail!("errordetectratio must be greater than zero");
    }
    if config.hist_columns == 0 || config.hist_columns > 3 {
        bail!("histcol must be 1, 2, or 3");
    }
    if config.hist_len < 2 {
        bail!("histlen must be at least 1");
    }
    if config.in2.is_some() {
        if config.out2.is_some() && config.out1.is_none() {
            bail!("out2 requires out=<file> for paired input");
        }
        if config.out_toss2.is_some() && config.out_toss1.is_none() {
            bail!("outt2 requires outt=<file> for paired input");
        }
        if config.out_low2.is_some() && config.out_low1.is_none() {
            bail!("outlow2 requires outlow=<file> for paired input");
        }
        if config.out_mid2.is_some() && config.out_mid1.is_none() {
            bail!("outmid2 requires outmid=<file> for paired input");
        }
        if config.out_high2.is_some() && config.out_high1.is_none() {
            bail!("outhigh2 requires outhigh=<file> for paired input");
        }
        if config.out_uncorrected2.is_some() && config.out_uncorrected1.is_none() {
            bail!("outuncorrected2 requires outuncorrected=<file> for paired input");
        }
    } else if config.interleaved {
        if config.out2.is_some() && config.out1.is_none() {
            bail!("out2 requires out=<file> for interleaved input");
        }
        if config.out_toss2.is_some() && config.out_toss1.is_none() {
            bail!("outt2 requires outt=<file> for interleaved input");
        }
        if config.out_low2.is_some() && config.out_low1.is_none() {
            bail!("outlow2 requires outlow=<file> for interleaved input");
        }
        if config.out_mid2.is_some() && config.out_mid1.is_none() {
            bail!("outmid2 requires outmid=<file> for interleaved input");
        }
        if config.out_high2.is_some() && config.out_high1.is_none() {
            bail!("outhigh2 requires outhigh=<file> for interleaved input");
        }
        if config.out_uncorrected2.is_some() && config.out_uncorrected1.is_none() {
            bail!("outuncorrected2 requires outuncorrected=<file> for interleaved input");
        }
    } else if !config.test_interleaved && (config.out2.is_some() || config.out_toss2.is_some()) {
        bail!("out2/outt2 require paired input with in2=<file> or interleaved=t");
    } else if !config.test_interleaved
        && (config.out_low2.is_some()
            || config.out_mid2.is_some()
            || config.out_high2.is_some()
            || config.out_uncorrected2.is_some())
    {
        bail!(
            "outlow2/outmid2/outhigh2/outuncorrected2 require paired input with in2=<file> or interleaved=t"
        );
    }
    Ok(())
}

fn validate_extra_inputs(config: &Config) -> Result<()> {
    for extra in &config.extra {
        if !extra.exists() || !extra.is_file() {
            bail!("extra input {} does not exist", extra.display());
        }
    }
    Ok(())
}

fn expand_hash_paired_input(config: &mut Config) {
    let Some(input) = config.in1.as_ref() else {
        return;
    };
    if input.exists() {
        return;
    }
    let text = input.to_string_lossy().into_owned();
    if !text.contains('#') {
        return;
    }

    config.in1 = Some(PathBuf::from(text.replacen('#', "1", 1)));
    config.in2 = Some(PathBuf::from(text.replacen('#', "2", 1)));
}

fn path(value: &str) -> PathBuf {
    PathBuf::from(value)
}

fn split_paths(value: &str) -> Vec<PathBuf> {
    value
        .split(',')
        .filter(|part| !part.trim().is_empty())
        .map(|part| PathBuf::from(part.trim()))
        .collect()
}

fn extra_paths(value: &str) -> Vec<PathBuf> {
    let trimmed = value.trim();
    if trimmed.is_empty() || trimmed.eq_ignore_ascii_case("null") {
        return Vec::new();
    }
    let literal = PathBuf::from(trimmed);
    if literal.exists() {
        vec![literal]
    } else {
        split_paths(trimmed)
    }
}

fn parse_bool(value: &str, key: &str) -> Result<bool> {
    match value.to_ascii_lowercase().as_str() {
        "t" | "true" | "1" | "y" | "yes" => Ok(true),
        "f" | "false" | "0" | "n" | "no" => Ok(false),
        _ => bail!("{key} expects a boolean value, got {value}"),
    }
}

fn quality_recal_base_key(key: &str) -> &str {
    key.strip_suffix("_p1")
        .or_else(|| key.strip_suffix("_p2"))
        .unwrap_or(key)
}

fn is_quality_recal_bool_key(key: &str) -> bool {
    matches!(
        quality_recal_base_key(key),
        "trackall"
            | "clearmatrices"
            | "loadq102"
            | "loadqap"
            | "loadqbp"
            | "loadqpt"
            | "loadqbt"
            | "loadq10"
            | "loadq12"
            | "loadqb12"
            | "loadqb012"
            | "loadqb123"
            | "loadqb234"
            | "loadq12b12"
            | "loadqp"
            | "loadq"
            | "recalwithposition"
            | "recalwithpos"
            | "recalusepos"
            | "recaltile"
            | "recaltiles"
            | "usetiles"
    )
}

fn parse_java_bool(value: &str) -> bool {
    if value.is_empty() {
        return true;
    }
    if value.len() == 1 {
        let byte = value.as_bytes()[0].to_ascii_lowercase();
        return byte == b't' || byte == b'1';
    }
    if value.eq_ignore_ascii_case("null") || value.eq_ignore_ascii_case("none") {
        return false;
    }
    value.eq_ignore_ascii_case("true")
}

fn parse_mpi_enabled(value: &str, key: &str) -> Result<bool> {
    if value
        .as_bytes()
        .first()
        .is_some_and(|byte| byte.is_ascii_digit())
    {
        Ok(parse_i32(value, key)? > 0)
    } else {
        Ok(parse_java_bool(value))
    }
}

enum CardinalityToggle {
    Bool(bool),
    Int(usize),
}

fn parse_cardinality_bool_or_int(value: &str, key: &str) -> Result<CardinalityToggle> {
    if value
        .as_bytes()
        .first()
        .is_some_and(|byte| byte.is_ascii_digit())
    {
        Ok(CardinalityToggle::Int(parse_cardinality_k(value, key)?))
    } else {
        Ok(CardinalityToggle::Bool(parse_bool(value, key)?))
    }
}

fn parse_cardinality_k(value: &str, key: &str) -> Result<usize> {
    let parsed = parse_i32(value, key)?;
    if parsed <= 0 {
        bail!("{key} expects a positive integer, got {value}");
    }
    usize::try_from(parsed).map_err(|_| anyhow::anyhow!("{key} value is out of range: {value}"))
}

fn parse_cardinality_buckets(value: &str, key: &str) -> Result<usize> {
    let buckets = parse_kmg_i64(value, key)?;
    if buckets <= 0 {
        bail!("{key} expects a positive KMG value, got {value}");
    }
    let buckets = usize::try_from(buckets)
        .map_err(|_| anyhow::anyhow!("{key} value is out of range: {value}"))?;
    if buckets > CARDINALITY_MAX_BUCKETS {
        bail!(
            "{key} requests {buckets} cardinality buckets, above the Rust safety cap of {CARDINALITY_MAX_BUCKETS}"
        );
    }
    Ok(buckets)
}

fn parse_cardinality_seed(value: &str, key: &str) -> Result<u64> {
    let parsed = parse_i64(value, key)?;
    if parsed < 0 {
        Ok(parsed as u64)
    } else {
        Ok(u64::try_from(parsed)
            .map_err(|_| anyhow::anyhow!("{key} value is out of range: {value}"))?)
    }
}

fn parse_kcount_cell_bits(value: &str, key: &str) -> Result<u8> {
    let bits = parse_i64(value, key)?;
    if bits <= 0 || bits > 32 || !(bits as u64).is_power_of_two() {
        bail!("{key} expects a power-of-two integer from 1 to 32, got {value}");
    }
    Ok(bits as u8)
}

fn parse_kcount_hashes(value: &str, key: &str) -> Result<usize> {
    let hashes = parse_i64(value, key)?;
    if !(1..=8).contains(&hashes) {
        bail!("{key} expects an integer from 1 to 8, got {value}");
    }
    Ok(hashes as usize)
}

fn parse_prefilter_hashes(value: &str, key: &str) -> Result<usize> {
    let hashes = parse_i64(value, key)?;
    if !(0..=8).contains(&hashes) {
        bail!("{key} expects an integer from 0 to 8, got {value}");
    }
    Ok(hashes as usize)
}

fn parse_matrixbits_cells(value: &str, key: &str) -> Result<usize> {
    let bits = parse_i64(value, key)?;
    if !(1..63).contains(&bits) {
        bail!("{key} expects an integer exponent from 1 to 62, got {value}");
    }
    1usize
        .checked_shl(bits as u32)
        .with_context(|| format!("{key} exponent is too large for this platform: {value}"))
}

fn parse_auto_or_i32(value: &str, key: &str) -> Result<()> {
    if !value.eq_ignore_ascii_case("auto") {
        let _ = parse_i32(value, key)?;
    }
    Ok(())
}

fn parse_auto_or_kmg_i64(value: &str, key: &str) -> Result<()> {
    if !value.eq_ignore_ascii_case("auto") {
        let _ = parse_kmg_i64(value, key)?;
    }
    Ok(())
}

fn parse_preallocation_fraction(value: &str, key: &str) -> Result<Option<f64>> {
    if value
        .as_bytes()
        .first()
        .is_some_and(|byte| byte.is_ascii_digit() || *byte == b'.')
    {
        let fraction = parse_f64(value, key)?;
        if !(0.0..=1.0).contains(&fraction) {
            bail!("{key} expects a fraction between 0 and 1 or a boolean value, got {value}");
        }
        Ok((fraction > 0.0).then_some(fraction))
    } else if parse_java_bool(value) {
        Ok(Some(1.0))
    } else {
        Ok(None)
    }
}

fn parse_fraction_micros(value: &str, key: &str) -> Result<u32> {
    let fraction = parse_f64(value, key)?;
    if !(0.0..=1.0).contains(&fraction) {
        bail!("{key} expects a fraction between 0 and 1, got {value}");
    }
    Ok((fraction * 1_000_000.0).round() as u32)
}

fn parse_monitor(value: &str, key: &str) -> Result<()> {
    if value
        .as_bytes()
        .first()
        .is_some_and(|byte| byte.is_ascii_digit() || *byte == b'.')
    {
        let mut parts = value.split(',');
        let first = parts.next().unwrap_or_default();
        parse_f64(first, key)?;
        if let Some(second) = parts.next() {
            parse_f64(second, key)?;
        }
        if parts.next().is_some() {
            bail!("{key} expects one or two numeric watchdog values, got {value}");
        }
    } else {
        let _ = parse_java_bool(value);
    }
    Ok(())
}

fn parse_qtrim(config: &mut Config, value: &str, key: &str) -> Result<()> {
    let lower = value.to_ascii_lowercase();
    match lower.as_str() {
        "" => {
            config.trim_left = true;
            config.trim_right = true;
        }
        "left" | "l" => {
            config.trim_left = true;
            config.trim_right = false;
        }
        "right" | "r" => {
            config.trim_left = false;
            config.trim_right = true;
        }
        "both" | "rl" | "lr" => {
            config.trim_left = true;
            config.trim_right = true;
        }
        "window" | "w" => {
            config.trim_left = false;
            config.trim_right = true;
            config.trim_window = true;
            config.trim_optimal = false;
            config.trim_optimal_bias = None;
        }
        _ if lower.starts_with("window,") || lower.starts_with("w,") => {
            let Some((_, length)) = value.split_once(',') else {
                unreachable!("guard requires a comma");
            };
            config.trim_window_length = parse_usize(length, key)?;
            config.trim_left = false;
            config.trim_right = true;
            config.trim_window = true;
            config.trim_optimal = false;
            config.trim_optimal_bias = None;
        }
        _ if value
            .as_bytes()
            .first()
            .is_some_and(|byte| byte.is_ascii_digit()) =>
        {
            config.trim_quality = parse_trim_quality(config, value, key)?;
            config.trim_right = true;
        }
        _ => {
            let enabled = parse_bool(value, key)?;
            config.trim_left = enabled;
            config.trim_right = enabled;
        }
    }
    Ok(())
}

fn parse_trim_quality(config: &mut Config, value: &str, key: &str) -> Result<f64> {
    if value.contains(',') {
        let mut parts = value.split(',');
        let first = parts.next().unwrap_or_default();
        let trim_quality = parse_f64(first, key)?;
        for part in parts {
            parse_f64(part, key)?;
        }
        config.notes.push(format!(
            "{key}={value} requests position-specific trim qualities; Rust uses the first threshold {trim_quality} for the supported trimming path"
        ));
        return Ok(trim_quality);
    }
    parse_f64(value, key)
}

fn parse_poly(value: &str, key: &str) -> Result<usize> {
    if value.is_empty() {
        bail!("{key} expects a polymer threshold or boolean value");
    }
    if value
        .as_bytes()
        .first()
        .is_some_and(|byte| byte.is_ascii_digit())
    {
        parse_usize(value, key)
    } else {
        Ok(if parse_bool(value, key)? { 2 } else { 0 })
    }
}

fn parse_optitrim(config: &mut Config, value: &str, key: &str) -> Result<()> {
    if value
        .as_bytes()
        .first()
        .is_some_and(|byte| *byte == b'.' || byte.is_ascii_digit())
    {
        let bias = parse_f64(value, key)?;
        if !(0.0..1.0).contains(&bias) {
            bail!("{key} bias must be greater than or equal to 0 and less than 1");
        }
        config.trim_optimal = true;
        config.trim_optimal_bias = Some(bias);
    } else {
        config.trim_optimal = parse_bool(value, key)?;
        config.trim_optimal_bias = None;
    }
    Ok(())
}

fn enable_error_correction_if_unset(config: &mut Config) {
    if !config.error_correct_first && !config.error_correct_final {
        config.error_correct_first = true;
    }
    config.error_correct = config.error_correct_first || config.error_correct_final;
}

fn parse_u8(value: &str, key: &str) -> Result<u8> {
    value
        .parse::<u8>()
        .map_err(|_| anyhow::anyhow!("{key} expects an integer, got {value}"))
}

fn parse_i8(value: &str, key: &str) -> Result<i8> {
    value
        .parse::<i8>()
        .map_err(|_| anyhow::anyhow!("{key} expects a byte integer, got {value}"))
}

fn parse_usize(value: &str, key: &str) -> Result<usize> {
    value
        .parse::<usize>()
        .map_err(|_| anyhow::anyhow!("{key} expects a non-negative integer, got {value}"))
}

fn parse_u64(value: &str, key: &str) -> Result<u64> {
    value
        .parse::<u64>()
        .map_err(|_| anyhow::anyhow!("{key} expects a non-negative integer, got {value}"))
}

fn parse_i64(value: &str, key: &str) -> Result<i64> {
    value
        .parse::<i64>()
        .map_err(|_| anyhow::anyhow!("{key} expects an integer, got {value}"))
}

fn parse_i32(value: &str, key: &str) -> Result<i32> {
    value
        .parse::<i32>()
        .map_err(|_| anyhow::anyhow!("{key} expects an integer, got {value}"))
}

fn parse_i32_clamped(value: &str, key: &str, min: i32, max: i32) -> Result<i32> {
    Ok(parse_i32(value, key)?.clamp(min, max))
}

fn parse_f64(value: &str, key: &str) -> Result<f64> {
    value
        .parse::<f64>()
        .map_err(|_| anyhow::anyhow!("{key} expects a number, got {value}"))
}

fn parse_min_average_quality(value: &str, key: &str) -> Result<()> {
    let mut parts = value.split(',');
    let quality = parts.next().unwrap_or_default();
    parse_f64(quality, key)?;
    if let Some(bases) = parts.next() {
        parse_usize(bases, key)?;
    }
    if parts.next().is_some() {
        bail!("{key} expects quality or quality,bases, got {value}");
    }
    Ok(())
}

fn parse_quality_offset(value: &str, key: &str) -> Result<Option<u8>> {
    match value.to_ascii_lowercase().as_str() {
        "auto" => Ok(None),
        "sanger" => Ok(Some(33)),
        "illumina" => Ok(Some(64)),
        "33" => Ok(Some(33)),
        "64" => Ok(Some(64)),
        _ => bail!("{key} expects auto, sanger, illumina, 33, or 64, got {value}"),
    }
}

fn parse_fake_fasta_quality(config: &mut Config, value: &str) -> Result<()> {
    if value.is_empty() {
        return Ok(());
    }
    if value.as_bytes()[0].is_ascii_alphabetic() {
        let _ = parse_bool(value, "fakefastaquality")?;
        return Ok(());
    }

    let parsed = parse_i32(value, "fakefastaquality")?;
    if parsed > 0 {
        config.fake_quality = parsed.min(50) as u8;
    }
    Ok(())
}

fn parse_fasta_wrap(value: &str, key: &str) -> Result<usize> {
    let parsed = parse_kmg_i64(value, key)?;
    if parsed < 1 {
        Ok(0)
    } else {
        usize::try_from(parsed).map_err(|_| anyhow::anyhow!("{key} value is out of range: {value}"))
    }
}

fn parse_junk_mode(config: &mut Config, value: &str) -> Result<()> {
    match value.to_ascii_lowercase().as_str() {
        "ignore" => {
            config.fix_junk_and_iupac = false;
            config.junk_mode = JunkMode::Ignore;
        }
        "crash" | "fail" => {
            config.fix_junk_and_iupac = false;
            config.junk_mode = JunkMode::Crash;
        }
        "fix" => {
            config.fix_junk_and_iupac = false;
            config.junk_mode = JunkMode::Fix;
        }
        "flag" | "discard" => {
            config.fix_junk_and_iupac = false;
            config.junk_mode = JunkMode::Flag;
        }
        "iupacton" => {
            config.fix_junk_and_iupac = true;
            config.junk_mode = JunkMode::Fix;
        }
        _ => {
            bail!("junk expects ignore, crash, fail, fix, flag, discard, or iupacton, got {value}")
        }
    }
    Ok(())
}

fn parse_percent(value: &str, key: &str) -> Result<f64> {
    let mut parsed = parse_f64(value, key)?;
    if parsed > 1.0 && parsed <= 100.0 {
        parsed /= 100.0;
    }
    if !(0.0..=1.0).contains(&parsed) {
        bail!("{key} must be between 0 and 100");
    }
    Ok(parsed)
}

fn parse_limit(value: &str, key: &str) -> Result<Option<u64>> {
    let parsed = parse_kmg_i64(value, key)?;
    if parsed < 0 {
        Ok(None)
    } else {
        Ok(Some(parsed as u64))
    }
}

fn parse_kmg_i64(value: &str, key: &str) -> Result<i64> {
    let lower = value.to_ascii_lowercase();
    if matches!(lower.as_str(), "big" | "inf" | "infinity" | "max" | "huge") {
        return Ok(i64::MAX);
    }

    let Some(last) = lower.chars().last() else {
        bail!("{key} expects an integer or KMG value, got {value}");
    };
    let (number, multiplier) = match last {
        'k' => (&value[..value.len() - 1], 1_000_f64),
        'm' => (&value[..value.len() - 1], 1_000_000_f64),
        'g' | 'b' => (&value[..value.len() - 1], 1_000_000_000_f64),
        't' => (&value[..value.len() - 1], 1_000_000_000_000_f64),
        'p' | 'q' => (&value[..value.len() - 1], 1_000_000_000_000_000_f64),
        'e' => (&value[..value.len() - 1], 1_000_000_000_000_000_000_f64),
        'c' | 'h' => (&value[..value.len() - 1], 100_f64),
        'd' => (&value[..value.len() - 1], 10_f64),
        _ if last.is_ascii_alphabetic() => {
            bail!("{key} has an unsupported KMG suffix in {value}");
        }
        _ => (value, 1_f64),
    };

    if number
        .chars()
        .last()
        .is_some_and(|char| char.is_ascii_alphabetic())
    {
        bail!("{key} has too many suffix letters in {value}");
    }

    let parsed = if number.contains('.') || multiplier != 1.0 {
        let scaled = number
            .parse::<f64>()
            .map_err(|_| anyhow::anyhow!("{key} expects an integer or KMG value, got {value}"))?
            * multiplier;
        if scaled > i64::MAX as f64 || scaled < i64::MIN as f64 {
            bail!("{key} value is out of range: {value}");
        }
        scaled as i64
    } else {
        number
            .parse::<i64>()
            .map_err(|_| anyhow::anyhow!("{key} expects an integer or KMG value, got {value}"))?
    };
    Ok(parsed)
}

fn parse_kmg_usize(value: &str, key: &str) -> Result<usize> {
    let parsed = parse_kmg_i64(value, key)?;
    if parsed < 0 {
        bail!("{key} expects a non-negative KMG value, got {value}");
    }
    usize::try_from(parsed).map_err(|_| anyhow::anyhow!("{key} value is out of range: {value}"))
}

fn parse_positive_kmg_usize(value: &str, key: &str) -> Result<usize> {
    let parsed = parse_kmg_usize(value, key)?;
    if parsed == 0 {
        bail!("{key} expects a positive KMG value, got {value}");
    }
    Ok(parsed)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn parse(values: &[&str]) -> Config {
        let mut args: Vec<OsString> = values.iter().map(OsString::from).collect();
        if !values.iter().any(|value| is_pass_selector(value)) {
            args.push(OsString::from("passes=1"));
        }
        parse_args(args).unwrap()
    }

    fn is_pass_selector(value: &str) -> bool {
        let lower = value.to_ascii_lowercase();
        matches!(lower.as_str(), "1pass" | "1p" | "2pass" | "2p")
            || lower.split_once('=').is_some_and(|(key, _)| {
                matches!(key, "passes" | "p" | "1pass" | "1p" | "2pass" | "2p")
            })
    }

    #[test]
    fn implicit_bbnorm_default_keeps_two_pass_mode() {
        let cfg = parse_args(["in=reads.fq"].into_iter().map(OsString::from)).unwrap();
        assert_eq!(cfg.passes, 2);
    }

    #[test]
    fn one_pass_aliases_select_supported_single_pass_like_bbnorm() {
        let cfg = parse_args(["in=reads.fq", "1pass"].into_iter().map(OsString::from)).unwrap();
        assert_eq!(cfg.passes, 1);

        let cfg = parse_args(["in=reads.fq", "1pass=f"].into_iter().map(OsString::from)).unwrap();
        assert_eq!(cfg.passes, 1);
    }

    #[test]
    fn two_pass_aliases_select_multipass_like_bbnorm() {
        let cfg = parse_args(["in=reads.fq", "2pass=f"].into_iter().map(OsString::from)).unwrap();
        assert_eq!(cfg.passes, 2);
    }

    #[test]
    fn parses_core_aliases() {
        let cfg = parse(&[
            "reads.fq",
            "out=keep.fq",
            "outt=toss.fq",
            "hist=hist.tsv",
            "k=21",
            "min=3",
            "max=9",
            "minkmers=2",
            "ml=42",
            "dp=60",
            "tbr=t",
            "rbb=t",
            "srr=t",
            "overwrite=t",
            "append=t",
        ]);
        assert_eq!(cfg.in1.unwrap(), PathBuf::from("reads.fq"));
        assert_eq!(cfg.out1.unwrap(), PathBuf::from("keep.fq"));
        assert_eq!(cfg.out_toss1.unwrap(), PathBuf::from("toss.fq"));
        assert_eq!(cfg.hist_in.unwrap(), PathBuf::from("hist.tsv"));
        assert_eq!(cfg.k, 21);
        assert_eq!(cfg.min_depth, 3);
        assert_eq!(cfg.max_depth, Some(100));
        assert_eq!(cfg.min_kmers_over_min_depth, 2);
        assert_eq!(cfg.min_length, 42);
        assert!((cfg.depth_percentile - 0.60).abs() < f64::EPSILON);
        assert!(cfg.toss_error_reads);
        assert!(cfg.require_both_bad);
        assert!(cfg.save_rare_reads);
        assert!(cfg.overwrite);
        assert!(cfg.append);
    }

    #[test]
    fn accepts_shared_input_output_file_aliases() {
        let cfg = parse(&[
            "input=reads1.fq",
            "input2=reads2.fq",
            "output=keep1.fq",
            "output2=keep2.fq",
        ]);

        assert_eq!(cfg.in1.unwrap(), PathBuf::from("reads1.fq"));
        assert_eq!(cfg.in2.unwrap(), PathBuf::from("reads2.fq"));
        assert_eq!(cfg.out1.unwrap(), PathBuf::from("keep1.fq"));
        assert_eq!(cfg.out2.unwrap(), PathBuf::from("keep2.fq"));
    }

    #[test]
    fn parses_bare_boolean_flags_like_bbnorm() {
        let cfg = parse_args(
            [
                "reads.fq",
                "prefilter",
                "countup",
                "keepall",
                "ecc",
                "ecco",
                "ow",
            ]
            .into_iter()
            .map(OsString::from),
        )
        .unwrap();
        assert_eq!(cfg.in1.unwrap(), PathBuf::from("reads.fq"));
        assert!(cfg.in2.is_none());
        assert!(cfg.prefilter.enabled);
        assert!(!cfg.prefilter.force_disabled);
        assert!(cfg.count_up);
        assert!(cfg.keep_all);
        assert!(cfg.error_correct);
        assert!(cfg.overlap_error_correct);
        assert!(cfg.overwrite);

        let cfg = parse_args(
            ["in=x.fq", "prefilter", "prefilter=f"]
                .into_iter()
                .map(OsString::from),
        )
        .unwrap();
        assert!(!cfg.prefilter.enabled);
        assert!(cfg.prefilter.force_disabled);

        let cfg = parse_args(
            ["in=x.fq", "prefilter=f", "prefilter"]
                .into_iter()
                .map(OsString::from),
        )
        .unwrap();
        assert!(cfg.prefilter.enabled);
        assert!(!cfg.prefilter.force_disabled);
    }

    #[test]
    fn clamps_max_depth_and_minkmers_like_bbnorm() {
        let cfg = parse(&["in=reads.fq", "target=100", "max=50", "minkmers=0"]);
        assert_eq!(cfg.target_depth, 100);
        assert_eq!(cfg.max_depth, Some(100));
        assert_eq!(cfg.min_kmers_over_min_depth, 1);

        let cfg = parse(&["in=reads.fq", "max=150", "target=100"]);
        assert_eq!(cfg.max_depth, Some(150));
    }

    #[test]
    fn parses_fixspikes_aliases() {
        let cfg = parse(&["in=reads.fq", "fixspikes=t"]);
        assert!(cfg.fix_spikes);

        let cfg = parse(&["in=reads.fq", "fs=f"]);
        assert!(!cfg.fix_spikes);
    }

    #[test]
    fn parses_kmg_read_limits_like_bbnorm() {
        let cfg = parse(&["in=reads.fq", "reads=0.01k", "tablereads=1d"]);
        assert_eq!(cfg.max_reads, Some(10));
        assert_eq!(cfg.table_reads, Some(10));

        let cfg = parse(&["in=reads.fq", "reads=-1", "tablereads=max"]);
        assert_eq!(cfg.max_reads, None);
        assert_eq!(cfg.table_reads, Some(i64::MAX as u64));
    }

    #[test]
    fn parses_kmg_min_length_like_bbnorm() {
        let cfg = parse(&["in=reads.fq", "minlen=0.101k"]);
        assert_eq!(cfg.min_length, 101);
    }

    #[test]
    fn parses_quality_trimming_like_bbnorm() {
        let cfg = parse(&["in=reads.fq", "qtrim=r", "trimq=10"]);
        assert!(!cfg.trim_left);
        assert!(cfg.trim_right);
        assert!((cfg.trim_quality - 10.0).abs() < f64::EPSILON);

        let cfg = parse(&["in=reads.fq", "qtrim=12"]);
        assert!(!cfg.trim_left);
        assert!(cfg.trim_right);
        assert!((cfg.trim_quality - 12.0).abs() < f64::EPSILON);

        let cfg = parse(&["in=reads.fq", "qtrim=r", "trimq=10,20"]);
        assert!(!cfg.trim_left);
        assert!(cfg.trim_right);
        assert!((cfg.trim_quality - 10.0).abs() < f64::EPSILON);
        assert!(cfg.notes.iter().any(|note| note.contains("trimq=10,20")));

        let cfg = parse(&["in=reads.fq", "qtrim=12,20"]);
        assert!(!cfg.trim_left);
        assert!(cfg.trim_right);
        assert!((cfg.trim_quality - 12.0).abs() < f64::EPSILON);
        assert!(cfg.notes.iter().any(|note| note.contains("qtrim=12,20")));

        let cfg = parse(&["in=reads.fq", "qtrim=t", "optitrim=f", "trimgoodinterval=3"]);
        assert!(cfg.trim_left);
        assert!(cfg.trim_right);
        assert!(!cfg.trim_optimal);
        assert_eq!(cfg.trim_min_good_interval, 3);

        let cfg = parse(&["in=reads.fq", "qtrim=w,5"]);
        assert!(!cfg.trim_left);
        assert!(cfg.trim_right);
        assert!(cfg.trim_window);
        assert!(!cfg.trim_optimal);
        assert_eq!(cfg.trim_window_length, 5);
    }

    #[test]
    fn parses_quality_output_offset_like_bbnorm() {
        let cfg = parse(&["in=reads.fq", "qin=64", "qout=64"]);
        assert_eq!(cfg.quality_in_offset, 64);
        assert_eq!(cfg.quality_out_offset, 64);

        let cfg = parse(&["in=reads.fq", "qout=auto", "qin=sanger"]);
        assert_eq!(cfg.quality_in_offset, 33);
        assert_eq!(cfg.quality_out_offset, 33);

        let cfg = parse(&["in=reads.fq", "qual=illumina"]);
        assert_eq!(cfg.quality_in_offset, 64);
        assert_eq!(cfg.quality_out_offset, 64);

        let cfg = parse(&["in=reads.fq", "asciiin=64", "qualityout=64"]);
        assert_eq!(cfg.quality_in_offset, 64);
        assert_eq!(cfg.quality_out_offset, 64);

        let cfg = parse(&["in=reads.fq", "qauto=t"]);
        assert_eq!(cfg.quality_in_offset, 33);
        assert_eq!(cfg.quality_out_offset, 33);
        assert!(cfg.notes.iter().any(|note| note.contains("qauto")));

        let cfg = parse(&["in=reads.fq", "qin=64", "qauto=f", "qout=64"]);
        assert_eq!(cfg.quality_in_offset, 64);
        assert_eq!(cfg.quality_out_offset, 64);
    }

    #[test]
    fn parses_quality_change_controls_like_bbnorm() {
        let cfg = parse(&[
            "in=reads.fq",
            "changequality=f",
            "mincalledquality=5",
            "maxcalledquality=30",
        ]);
        assert!(!cfg.change_quality);
        assert_eq!(cfg.min_called_quality, 5);
        assert_eq!(cfg.max_called_quality, 30);

        let cfg = parse(&[
            "in=reads.fq",
            "cq=t",
            "mincalledquality=-5",
            "maxcalledquality=200",
        ]);
        assert!(cfg.change_quality);
        assert_eq!(cfg.min_called_quality, 0);
        assert_eq!(cfg.max_called_quality, 93);

        let cfg = parse(&["in=reads.fq", "ignorebadquality=t"]);
        assert!(!cfg.change_quality);

        let cfg = parse(&["in=reads.fq", "ibq=t"]);
        assert!(!cfg.change_quality);

        let cfg = parse(&["in=reads.fq", "changequality=f", "ignorebadquality=f"]);
        assert!(!cfg.change_quality);

        let cfg = parse(&["in=reads.fq", "ignorebadquality=t", "changequality=t"]);
        assert!(cfg.change_quality);
    }

    #[test]
    fn parses_fake_quality_controls_like_bbnorm() {
        let cfg = parse(&["in=reads.fa", "fakequality=20"]);
        assert_eq!(cfg.fake_quality, 20);

        let cfg = parse(&["in=reads.fa", "qfake=15"]);
        assert_eq!(cfg.fake_quality, 15);

        let cfg = parse(&["in=reads.fa", "fakefastaquality=80"]);
        assert_eq!(cfg.fake_quality, 50);

        let cfg = parse(&["in=reads.fa", "fakefastaquality=0"]);
        assert_eq!(cfg.fake_quality, 30);

        let cfg = parse(&["in=reads.fa", "ffq=t"]);
        assert_eq!(cfg.fake_quality, 30);
    }

    #[test]
    fn parses_fasta_wrap_like_bbnorm() {
        let cfg = parse(&["in=reads.fq"]);
        assert_eq!(cfg.fasta_wrap, 70);

        let cfg = parse(&["in=reads.fq", "fastawrap=20"]);
        assert_eq!(cfg.fasta_wrap, 20);

        let cfg = parse(&["in=reads.fq", "wrap=0"]);
        assert_eq!(cfg.fasta_wrap, 0);

        let cfg = parse(&["in=reads.fq", "wrap=-1"]);
        assert_eq!(cfg.fasta_wrap, 0);
    }

    #[test]
    fn accepts_thread_counts_like_bbnorm_as_rayon_controls() {
        let cfg = parse(&["in=reads.fq", "threads=2"]);
        assert_eq!(cfg.threads, Some(2));
        assert_eq!(cfg.gzip_threads, Some(2));
        assert!(
            cfg.notes
                .iter()
                .any(|note| note.contains("threads=2 accepted"))
        );
        assert!(
            cfg.notes
                .iter()
                .any(|note| note.contains("also enables gzip input/output workers"))
        );

        let cfg = parse(&["in=reads.fq", "threads=2", "zipthreads=1"]);
        assert_eq!(cfg.threads, Some(2));
        assert_eq!(cfg.gzip_threads, Some(1));

        let cfg = parse(&["in=reads.fq", "threads=2", "useunpigz=t"]);
        assert_eq!(cfg.gzip_threads, Some(2));

        let cfg = parse(&["in=reads.fq", "t=-1"]);
        assert_eq!(cfg.threads, None);
        assert!(cfg.notes.is_empty());

        let cfg = parse(&["in=reads.fq", "threads=auto"]);
        assert_eq!(cfg.threads, None);
        assert!(
            cfg.notes
                .iter()
                .any(|note| note.contains("threads=auto accepted"))
        );

        let cfg = parse(&["in=reads.fq", "threads=max"]);
        assert_eq!(
            cfg.threads,
            Some(
                std::thread::available_parallelism()
                    .map(|threads| threads.get())
                    .unwrap_or(1)
            )
        );
        assert!(
            cfg.notes
                .iter()
                .any(|note| note.contains("threads=max accepted"))
        );
    }

    #[test]
    fn accepts_build_step_size_controls_as_covered_noops() {
        for case in ["stepsize=2", "buildstepsize=4"] {
            let cfg = parse(&["in=reads.fq", case]);
            assert!(
                cfg.notes
                    .iter()
                    .any(|note| note.contains("trusted-kmer sampling control")),
                "missing trusted-kmer note for {case}: {:?}",
                cfg.notes
            );
        }

        for case in ["stepsize=abc", "buildstepsize=abc"] {
            let err = parse_args(
                ["in=reads.fq", "passes=1", case]
                    .into_iter()
                    .map(OsString::from),
            )
            .unwrap_err()
            .to_string();
            assert!(
                err.contains("expects"),
                "unexpected error for malformed {case}: {err}"
            );
        }
    }

    #[test]
    fn accepts_default_equivalent_sketch_controls_as_noops() {
        for case in [
            "bits=32",
            "bits1=16",
            "cbits1=16",
            "cellbits1=16",
            "hashes=3",
            "buildpasses=1",
            "prefilter=t",
        ] {
            let cfg = parse(&["in=reads.fq", case]);
            assert!(
                !cfg.notes.is_empty(),
                "expected an explanatory no-op note for {case}"
            );
            if case.contains("bits1") || case.contains("cbits1") || case.contains("cellbits1") {
                assert_eq!(
                    cfg.count_min_bits_first,
                    Some(16),
                    "expected first-pass bit width for {case}"
                );
            }
        }

        for case in ["bits1=abc", "cbits1=abc", "cellbits1=abc"] {
            let err = parse_args(
                ["in=reads.fq", "passes=1", case]
                    .into_iter()
                    .map(OsString::from),
            )
            .unwrap_err()
            .to_string();
            assert!(
                err.contains("expects"),
                "unexpected error for malformed {case}: {err}"
            );
        }
    }

    #[test]
    fn accepts_prefilter_controls_with_constrained_sketch_settings() {
        let cfg = parse_args(
            [
                "in=x.fq",
                "passes=1",
                "prefiltercells=1k",
                "prehashes=2",
                "pbits=8",
            ]
            .into_iter()
            .map(OsString::from),
        )
        .unwrap();
        assert_eq!(cfg.prefilter.cells, Some(1000));
        assert_eq!(cfg.prefilter.hashes, Some(2));
        assert_eq!(cfg.prefilter.bits, Some(8));
        assert_eq!(cfg.prefilter.memory_bytes, None);
        assert_eq!(cfg.prefilter.memory_fraction_micros, None);
        assert!(cfg.prefilter.enabled);
        assert!(!cfg.prefilter.force_disabled);
        assert!(
            cfg.notes
                .iter()
                .any(|note| note.contains("prefilter collision estimates")),
            "expected constrained prefilter note: {:?}",
            cfg.notes
        );

        let cfg = parse_args(
            ["in=x.fq", "passes=1", "precells=1k"]
                .into_iter()
                .map(OsString::from),
        )
        .unwrap();
        assert_eq!(cfg.prefilter.cells, Some(1000));
        assert!(cfg.prefilter.enabled);
        assert!(!cfg.prefilter.force_disabled);

        let cfg = parse_args(
            ["in=x.fq", "passes=1", "precells=0"]
                .into_iter()
                .map(OsString::from),
        )
        .unwrap();
        assert_eq!(cfg.prefilter.cells, None);
        assert!(!cfg.prefilter.enabled);
        assert!(!cfg.prefilter.force_disabled);

        let cfg = parse_args(
            ["in=x.fq", "passes=1", "prefilter=t", "prefiltercells=0"]
                .into_iter()
                .map(OsString::from),
        )
        .unwrap();
        assert_eq!(cfg.prefilter.cells, None);
        assert!(cfg.prefilter.enabled);
        assert!(!cfg.prefilter.force_disabled);

        let cfg = parse_args(
            ["in=x.fq", "passes=1", "prefilterhashes=1"]
                .into_iter()
                .map(OsString::from),
        )
        .unwrap();
        assert_eq!(cfg.prefilter.hashes, Some(1));
        assert!(cfg.prefilter.enabled);
        assert!(!cfg.prefilter.force_disabled);
        assert!(
            cfg.notes
                .iter()
                .any(|note| note.contains("prefilter collision estimates")),
            "expected implicit prefilter note: {:?}",
            cfg.notes
        );

        let cfg = parse_args(
            ["in=x.fq", "passes=1", "prehashes=0"]
                .into_iter()
                .map(OsString::from),
        )
        .unwrap();
        assert_eq!(cfg.prefilter.hashes, None);
        assert!(!cfg.prefilter.enabled);
        assert!(!cfg.prefilter.force_disabled);

        let cfg = parse_args(
            ["in=x.fq", "passes=1", "prefilter=t", "prehashes=0"]
                .into_iter()
                .map(OsString::from),
        )
        .unwrap();
        assert_eq!(cfg.prefilter.hashes, None);
        assert!(cfg.prefilter.enabled);
        assert!(!cfg.prefilter.force_disabled);

        let cfg = parse_args(
            ["in=x.fq", "passes=1", "prefiltermemory=1k"]
                .into_iter()
                .map(OsString::from),
        )
        .unwrap();
        assert_eq!(cfg.prefilter.memory_bytes, Some(1000));
        assert!(cfg.prefilter.enabled);
        assert!(!cfg.prefilter.force_disabled);
        assert!(
            cfg.notes
                .iter()
                .any(|note| note.contains("prefilter memory-sizing")),
            "expected memory-backed prefilter note: {:?}",
            cfg.notes
        );

        for case in ["prefiltersize=0.1", "prefilterfraction=0.1"] {
            let cfg = parse_args(
                ["in=x.fq", "passes=1", case]
                    .into_iter()
                    .map(OsString::from),
            )
            .unwrap();
            assert!(cfg.prefilter.enabled);
            assert!(!cfg.prefilter.force_disabled);
            assert_eq!(cfg.prefilter.memory_fraction_micros, Some(100_000));
            assert!(
                cfg.notes
                    .iter()
                    .any(|note| note.contains("prefilter collision memory")),
                "expected prefilter fraction note for {case}"
            );
        }

        for case in ["prefiltersize=0", "prefilterfraction=0"] {
            let cfg = parse_args(
                ["in=x.fq", "passes=1", case]
                    .into_iter()
                    .map(OsString::from),
            )
            .unwrap();
            assert!(!cfg.prefilter.enabled);
            assert!(cfg.prefilter.force_disabled);
            assert_eq!(cfg.prefilter.memory_fraction_micros, None);
            assert!(
                cfg.notes
                    .iter()
                    .any(|note| note.contains("disables fraction-derived")),
                "expected zero-fraction note for {case}"
            );
        }

        let cfg = parse_args(
            ["in=x.fq", "passes=1", "prefilter=t"]
                .into_iter()
                .map(OsString::from),
        )
        .unwrap();
        assert!(cfg.prefilter.enabled);
        assert!(!cfg.prefilter.force_disabled);
        assert!(
            cfg.notes
                .iter()
                .any(|note| note.contains("default prefilter partitioning")),
            "expected enabled prefilter note: {:?}",
            cfg.notes
        );

        let cfg = parse_args(
            ["in=x.fq", "passes=1", "prehashes=1", "prefilter=f"]
                .into_iter()
                .map(OsString::from),
        )
        .unwrap();
        assert_eq!(cfg.prefilter.hashes, Some(1));
        assert!(!cfg.prefilter.enabled);
        assert!(cfg.prefilter.force_disabled);

        let cfg = parse_args(
            [
                "in=x.fq",
                "passes=1",
                "prehashes=1",
                "prefilter=f",
                "prefilter=t",
            ]
            .into_iter()
            .map(OsString::from),
        )
        .unwrap();
        assert_eq!(cfg.prefilter.hashes, Some(1));
        assert!(cfg.prefilter.enabled);
        assert!(!cfg.prefilter.force_disabled);

        let cfg = parse_args(
            ["in=x.fq", "passes=1", "prefilter=f", "prehashes=1"]
                .into_iter()
                .map(OsString::from),
        )
        .unwrap();
        assert_eq!(cfg.prefilter.hashes, Some(1));
        assert!(cfg.prefilter.enabled);
        assert!(!cfg.prefilter.force_disabled);

        let cfg = parse_args(
            ["in=x.fq", "passes=1", "prefiltercells=1k", "prefilter=f"]
                .into_iter()
                .map(OsString::from),
        )
        .unwrap();
        assert_eq!(cfg.prefilter.cells, Some(1000));
        assert!(!cfg.prefilter.enabled);
        assert!(cfg.prefilter.force_disabled);

        let cfg = parse_args(
            [
                "in=x.fq",
                "passes=1",
                "prefilterfraction=0.1",
                "prefilter=f",
            ]
            .into_iter()
            .map(OsString::from),
        )
        .unwrap();
        assert_eq!(cfg.prefilter.memory_fraction_micros, Some(100_000));
        assert!(!cfg.prefilter.enabled);
        assert!(cfg.prefilter.force_disabled);

        let cfg = parse_args(
            [
                "in=x.fq",
                "passes=1",
                "prefilter=f",
                "prefilterfraction=0.1",
            ]
            .into_iter()
            .map(OsString::from),
        )
        .unwrap();
        assert_eq!(cfg.prefilter.memory_fraction_micros, Some(100_000));
        assert!(cfg.prefilter.enabled);
        assert!(!cfg.prefilter.force_disabled);

        let cfg = parse_args(
            ["in=x.fq", "passes=1", "buildpasses=2"]
                .into_iter()
                .map(OsString::from),
        )
        .unwrap();
        assert_eq!(cfg.build_passes, 2);
        assert!(
            cfg.notes
                .iter()
                .any(|note| note.contains("trusted-kmer filtering")),
            "expected build-pass trusted-filter note: {:?}",
            cfg.notes
        );
    }

    #[test]
    fn accepts_constrained_count_min_controls_as_real_sketch_settings() {
        let cfg = parse_args(
            ["in=x.fq", "passes=1", "bits=16", "hashes=2", "cells=1k"]
                .into_iter()
                .map(OsString::from),
        )
        .unwrap();
        assert_eq!(cfg.count_min.bits, Some(16));
        assert_eq!(cfg.count_min.hashes, Some(2));
        assert_eq!(cfg.count_min.cells, Some(1000));
        assert!(
            cfg.notes
                .iter()
                .any(|note| note.contains("fixed-memory count-min input sketch")),
            "expected fixed-memory count-min sketch note: {:?}",
            cfg.notes
        );

        let cfg = parse_args(
            ["in=x.fq", "passes=1", "matrixbits=10"]
                .into_iter()
                .map(OsString::from),
        )
        .unwrap();
        assert_eq!(cfg.count_min.cells, Some(1024));

        let cfg = parse_args(
            [
                "in=x.fq",
                "passes=1",
                "bits=8",
                "hashes=2",
                "sketchmemory=1k",
            ]
            .into_iter()
            .map(OsString::from),
        )
        .unwrap();
        assert_eq!(cfg.count_min.memory_bytes, Some(1000));
        assert!(
            cfg.notes
                .iter()
                .any(|note| note.contains("count-min memory budget")),
            "expected count-min memory-budget note: {:?}",
            cfg.notes
        );

        let cfg = parse_args(
            [
                "in=x.fq",
                "passes=1",
                "maxcountupspillbytes=64m",
                "maxcountupspillfinallivebytes=96m",
                "maxcountupspillwritebytes=128m",
                "maxcountupspillinitialruns=10",
                "maxcountupspillmergeruns=2",
                "maxcountupspillfinalruns=4",
            ]
            .into_iter()
            .map(OsString::from),
        )
        .unwrap();
        assert_eq!(cfg.max_countup_spill_live_bytes, Some(64_000_000));
        assert_eq!(cfg.max_countup_spill_final_live_bytes, Some(96_000_000));
        assert_eq!(cfg.max_countup_spill_write_bytes, Some(128_000_000));
        assert_eq!(cfg.max_countup_spill_initial_runs, Some(10));
        assert_eq!(cfg.max_countup_spill_merge_runs, Some(2));
        assert_eq!(cfg.max_countup_spill_final_runs, Some(4));
        assert!(
            cfg.notes
                .iter()
                .any(|note| note.contains("count-up temp-spill safety cap")),
            "expected count-up spill live cap note: {:?}",
            cfg.notes
        );
        assert!(
            cfg.notes
                .iter()
                .any(|note| note.contains("count-up temp-spill I/O safety cap")),
            "expected count-up spill write cap note: {:?}",
            cfg.notes
        );

        let cfg = parse_args(
            [
                "in=x.fq",
                "passes=1",
                "mem=2g",
                "autocountmin=f",
                "exact=t",
                "autosketchbytes=4m",
                "autocountminreads=500",
            ]
            .into_iter()
            .map(OsString::from),
        )
        .unwrap();
        assert_eq!(cfg.auto_count_min_memory_bytes, Some(2_000_000_000));
        assert!(!cfg.auto_count_min);
        assert!(cfg.force_exact_counts);
        assert_eq!(cfg.auto_count_min_input_bytes, 4_000_000);
        assert_eq!(cfg.auto_count_min_read_threshold, 500);

        for case in [
            "bits=abc",
            "bits=0",
            "bits=3",
            "bits=64",
            "hashes=abc",
            "hashes=0",
            "hashes=9",
            "cells=abc",
            "cells=0",
            "matrixbits=abc",
            "matrixbits=0",
            "matrixbits=64",
            "sketchmemory=abc",
            "sketchmemory=0",
            "maxcountupspillbytes=abc",
            "maxcountupspillfinallivebytes=abc",
            "maxcountupspillwritebytes=-1",
            "maxcountupspillinitialruns=abc",
            "maxcountupspillmergeruns=-1",
            "maxcountupspillfinalruns=abc",
            "mem=abc",
            "autosketchbytes=0",
            "autocountminreads=0x",
            "buildpasses=abc",
            "prehashes=abc",
            "prefilterhashes=abc",
            "prefiltercells=abc",
            "precells=abc",
            "prefiltersize=abc",
            "prefilterfraction=abc",
            "prefilterbits=abc",
            "prefilterbits=64",
            "prebits=abc",
            "prebits=3",
            "pbits=3",
            "prehashes=9",
        ] {
            let err = parse_args(
                ["in=x.fq", "passes=1", case]
                    .into_iter()
                    .map(OsString::from),
            )
            .unwrap_err()
            .to_string();
            assert!(
                err.contains("expects")
                    || err.contains("unsupported KMG suffix")
                    || err.contains("too many suffix letters"),
                "unexpected error for malformed {case}: {err}"
            );
        }
    }

    #[test]
    fn accepts_kmer_table_runtime_controls_as_working_fallbacks() {
        for case in [
            "initialsize=1k",
            "ways=31",
            "buflen=64k",
            "bufflen=64k",
            "bufferlength=64k",
            "tabletype=2",
            "rcomp=t",
            "maskmiddle=f",
            "showstats=t",
            "stats=f",
            "showspeed=f",
            "ss=t",
            "verbose2=t",
            "prealloc=0.25",
            "preallocate=f",
            "filtermemory=1k",
            "prefiltermemory=1k",
            "filtermem=1k",
            "filtermemoryoverride=1k",
            "minprobprefilter=f",
            "mpp=t",
            "minprobmain=t",
            "mpm=f",
            "prefilterpasses=auto",
            "prepasses=1",
            "onepass=t",
        ] {
            let cfg = parse_args(
                ["in=x.fq", "passes=1", case]
                    .into_iter()
                    .map(OsString::from),
            )
            .unwrap();
            assert!(
                cfg.notes.iter().any(|note| {
                    note.contains("kmer-table")
                        || note.contains("prefilter memory-sizing")
                        || note.contains("prefilter pass-count")
                }),
                "expected kmer-table fallback note for {case}: {:?}",
                cfg.notes
            );
        }

        let cfg = parse(&["in=x.fq", "passes=1", "initialsize=1k", "prealloc=0.25"]);
        assert_eq!(cfg.table_initial_size, Some(1000));
        assert_eq!(cfg.table_prealloc_fraction, Some(0.25));

        let cfg = parse(&["in=x.fq", "passes=1", "preallocate=t"]);
        assert_eq!(cfg.table_prealloc_fraction, Some(1.0));

        let cfg = parse(&["in=x.fq", "passes=1", "preallocate=f"]);
        assert_eq!(cfg.table_prealloc_fraction, None);

        for case in [
            "initialsize=abc",
            "ways=abc",
            "buflen=abc",
            "tabletype=abc",
            "prealloc=0.abc",
            "prealloc=1.5",
            "filtermemory=abc",
            "prepasses=abc",
        ] {
            let err = parse_args(
                ["in=x.fq", "passes=1", case]
                    .into_iter()
                    .map(OsString::from),
            )
            .unwrap_err()
            .to_string();
            assert!(
                err.contains("expects")
                    || err.contains("unsupported KMG suffix")
                    || err.contains("too many suffix letters"),
                "unexpected error for malformed {case}: {err}"
            );
        }
    }

    #[test]
    fn accepts_covered_runtime_noops_and_manual_auto_sizing_fallback() {
        for case in [
            "auto=t",
            "auto=f",
            "ordered=f",
            "verbose=t",
            "printcoverage=t",
            "tmpdir=/tmp",
            "usetmpdir=t",
            "usetmpdir=f",
            "usetempdir=f",
            "fastareadlen=4",
            "fastareadlength=4",
            "fastaminread=1",
            "fastaminlen=1",
            "fastaminlength=1",
            "forcesectionname=t",
            "fastadump=f",
        ] {
            let cfg = parse(&["in=reads.fq", case]);
            assert!(
                !cfg.notes.is_empty(),
                "expected an explanatory no-op note for {case}"
            );
        }

        for case in ["fastaminread=abc", "fastaminlen=abc", "fastaminlength=abc"] {
            let err = parse_args(
                ["in=reads.fq", "passes=1", case]
                    .into_iter()
                    .map(OsString::from),
            )
            .unwrap_err()
            .to_string();
            assert!(
                err.contains("expects"),
                "unexpected error for malformed {case}: {err}"
            );
        }
    }

    #[test]
    fn accepts_temporary_directory_controls_for_managed_temp_paths() {
        for case in [
            "tmpdir=/tmp/bbnorm",
            "usetmpdir=t",
            "usetmpdir=f",
            "usetempdir=t",
        ] {
            let cfg = parse(&["in=reads.fq", "passes=1", case]);
            assert!(
                cfg.notes
                    .iter()
                    .any(|note| note.contains("temporary-directory control")),
                "expected temporary-directory note for {case}: {:?}",
                cfg.notes
            );
        }
        let enabled = parse(&["in=reads.fq", "tmpdir=/tmp/bbnorm"]);
        assert_eq!(enabled.temp_dir, Some(PathBuf::from("/tmp/bbnorm")));
        assert!(enabled.use_temp_dir);

        let disabled = parse(&["in=reads.fq", "tmpdir=/tmp/bbnorm", "usetmpdir=f"]);
        assert_eq!(disabled.temp_dir, Some(PathBuf::from("/tmp/bbnorm")));
        assert!(!disabled.use_temp_dir);
    }

    #[test]
    fn parses_header_trimming_controls_like_bbnorm() {
        for case in ["trd=t", "trc=t", "trimreaddescriptions=f", "trimrname=t"] {
            let cfg = parse(&["in=reads.fq", case]);
            assert!(
                !cfg.notes.is_empty(),
                "expected an explanatory no-op note for {case}"
            );
        }
    }

    #[test]
    fn accepts_shared_io_runtime_controls_as_noops_and_validates_values() {
        let cfg = parse(&[
            "in=reads.fq",
            "null",
            "monitor=f",
            "killswitch=600,0.002",
            "json=t",
            "silent=t",
            "printexecuting=f",
            "proxyhost=localhost",
            "proxyport=8080",
            "metadatafile=metadata.json",
            "testsize=t",
            "extin=.fq.gz",
            "extout=.fq",
            "bufferbf=f",
            "bufferbf1=f",
            "usejni=f",
            "bytefile1=t",
            "bytefile2=maybe",
            "bf1bufferlen=64k",
            "bfthreads=1",
            "readbufferlength=64k",
            "readbufferdata=1m",
            "readbuffers=1",
            "workers=auto",
            "workerthreads=1",
            "wt=auto",
            "threadsin=1",
            "tin=auto",
            "threadsout=1",
            "tout=auto",
            "ziplevel=2",
            "pigz=2",
            "bgzip=f",
            "zipthreads=1",
            "ztd=2.0",
            "blocksize=128",
            "nativebgzip=f",
            "usebzip2=f",
            "skipvalidation=t",
            "validate=maybe",
            "vic=f",
            "usempi=f",
            "mpi=0",
            "crismpi=f",
            "mpikeepall=f",
            "tossbrokenreads=f",
            "nullifybrokenquality=f",
            "deleteoldname=f",
            "renamebymapping=f",
            "assertcigar=f",
            "parsecustom=f",
            "shrinkheaders=f",
            "fixheader=f",
            "allownullheader=f",
            "recalpairnum=f",
            "pairreads=f",
            "flipr2=f",
            "int=f",
            "testinterleaved=f",
            "forceinterleaved=f",
            "overrideinterleaved=t",
        ]);
        assert_eq!(cfg.notes.len(), 56);
        assert_eq!(cfg.gzip_threads, Some(1));

        for case in [
            "monitor=1,2,3",
            "bf1bufferlen=abc",
            "bfthreads=abc",
            "readbufferlength=abc",
            "readbuffers=abc",
            "workers=abc",
            "threadsin=abc",
            "threadsout=abc",
            "mpi=2k",
            "ziplevel=abc",
            "pigz=2k",
            "zipthreads=abc",
            "ztd=abc",
            "blocksize=abc",
        ] {
            let err = parse_args(
                ["in=reads.fq", "passes=1", case]
                    .into_iter()
                    .map(OsString::from),
            )
            .unwrap_err()
            .to_string();
            assert!(
                err.contains("expects") || err.contains("suffix"),
                "unexpected error for malformed {case}: {err}"
            );
        }

        for case in ["usempi=t", "mpi=2", "crismpi=t", "mpikeepall=t"] {
            let cfg = parse(&["in=reads.fq", "passes=1", case]);
            assert!(
                cfg.notes.iter().any(|note| note.contains("MPI")),
                "missing MPI fallback note for {case}: {:?}",
                cfg.notes
            );
        }

        for case in ["pairreads=t", "flipr2=t"] {
            let cfg = parse(&["in=reads.fq", "passes=1", case]);
            assert!(
                cfg.notes.iter().any(|note| note.contains("pairing")),
                "missing pairing fallback note for {case}: {:?}",
                cfg.notes
            );
        }
    }

    #[test]
    fn accepts_shared_sam_runtime_controls_as_fastq_noops_and_validates_values() {
        for case in [
            "sam=1.4",
            "samv=1.6",
            "samtools=f",
            "sambamba=f",
            "printHeaderWait=f",
            "nativebam=f",
            "prefernativebam=f",
            "userssw=f",
            "attachedsamline=f",
            "streamerthreads=1",
            "fastqstreamerthreads=1",
            "fastastreamerthreads=1",
            "samwriterthreads=1",
            "bamwriterthreads=1",
            "fastqwriterthreads=1",
            "fastastreamer2=f",
            "prefermd=f",
            "notags=f",
            "mdtag=f",
            "idtag=f",
            "mateqtag=f",
            "xmtag=f",
            "smtag=f",
            "amtag=f",
            "nmtag=f",
            "xttag=f",
            "stoptag=f",
            "lengthtag=f",
            "boundstag=f",
            "scoretag=f",
            "sortscaffolds=f",
            "customtag=f",
            "nhtag=f",
            "keepnames=f",
            "saa=f",
            "inserttag=f",
            "correctnesstag=f",
            "intronlen=10",
            "suppressheader=f",
            "noheadersequences=f",
            "tophat=f",
            "xs=us",
            "xstag=fr-ss",
            "flipsam=f",
            "readgroupid=rg1",
            "rgsm=sample",
        ] {
            let cfg = parse_args(
                ["in=reads.fq", "passes=1", case]
                    .into_iter()
                    .map(OsString::from),
            )
            .unwrap();
            assert!(
                cfg.notes
                    .iter()
                    .any(|note| note.contains("SAM") || note.contains("read-group")),
                "expected SAM/read-group no-op note for {case}: {:?}",
                cfg.notes
            );
        }

        for case in [
            "sam=abc",
            "streamerthreads=abc",
            "fastqwriterthreads=abc",
            "intronlen=abc",
        ] {
            let err = parse_args(
                ["in=reads.fq", "passes=1", case]
                    .into_iter()
                    .map(OsString::from),
            )
            .unwrap_err()
            .to_string();
            assert!(
                err.contains("expects") || err.contains("invalid float"),
                "unexpected error for malformed {case}: {err}"
            );
        }
    }

    #[test]
    fn accepts_side_output_stats_histograms_and_emits_quality_length_gc_and_base_histograms() {
        for case in [
            "qhist=qual.tsv",
            "bqhist=basequal.tsv",
            "qchist=qcount.tsv",
            "aqhist=avg.tsv",
            "obqhist=overall.tsv",
            "mhist=match.tsv",
            "ihist=insert.tsv",
            "bhist=base.tsv",
            "qahist=qacc.tsv",
            "indelhist=indel.tsv",
            "ehist=error.tsv",
            "lhist=length.tsv",
            "gchist=gc.tsv",
            "enthist=entropy.tsv",
            "barcodestats=barcode.tsv",
            "thist=time.tsv",
            "idhist=id.tsv",
            "gcbins=auto",
            "gchistbins=100",
            "entropybins=auto",
            "enthistbins=100",
            "idhistbins=auto",
            "idbins=100",
            "gcplot=f",
            "entropyns=t",
            "maxhistlen=1k",
            "fixindels=f",
        ] {
            let cfg = parse_args(
                ["in=reads.fq", "passes=1", case]
                    .into_iter()
                    .map(OsString::from),
            )
            .unwrap();
            assert!(
                cfg.notes.iter().any(|note| note.contains("side-output")),
                "expected side-output fallback note for {case}: {:?}",
                cfg.notes
            );
        }

        let cfg = parse_args(
            [
                "in=reads.fq",
                "passes=1",
                "qhist=quality.tsv",
                "bqhist=basequal.tsv",
                "qchist=qcount.tsv",
                "aqhist=avg.tsv",
                "obqhist=overall.tsv",
                "mhist=match.tsv",
                "ihist=insert.tsv",
                "qahist=qacc.tsv",
                "indelhist=indel.tsv",
                "ehist=error.tsv",
                "lhist=length.tsv",
                "gchist=gc.tsv",
                "bhist=base.tsv",
                "enthist=entropy.tsv",
                "idhist=id.tsv",
                "gcbins=100",
                "entropybins=100",
                "idbins=100",
                "maxhistlen=1k",
            ]
            .into_iter()
            .map(OsString::from),
        )
        .unwrap();
        assert_eq!(cfg.quality_hist_out, Some(PathBuf::from("quality.tsv")));
        assert_eq!(cfg.match_hist_out, Some(PathBuf::from("match.tsv")));
        assert_eq!(cfg.insert_hist_out, Some(PathBuf::from("insert.tsv")));
        assert_eq!(
            cfg.quality_accuracy_hist_out,
            Some(PathBuf::from("qacc.tsv"))
        );
        assert_eq!(cfg.indel_hist_out, Some(PathBuf::from("indel.tsv")));
        assert_eq!(cfg.error_hist_out, Some(PathBuf::from("error.tsv")));
        assert_eq!(
            cfg.base_quality_hist_out,
            Some(PathBuf::from("basequal.tsv"))
        );
        assert_eq!(
            cfg.quality_count_hist_out,
            Some(PathBuf::from("qcount.tsv"))
        );
        assert_eq!(cfg.average_quality_hist_out, Some(PathBuf::from("avg.tsv")));
        assert_eq!(
            cfg.overall_base_quality_hist_out,
            Some(PathBuf::from("overall.tsv"))
        );
        assert_eq!(cfg.length_hist_out, Some(PathBuf::from("length.tsv")));
        assert_eq!(cfg.gc_hist_out, Some(PathBuf::from("gc.tsv")));
        assert_eq!(cfg.base_hist_out, Some(PathBuf::from("base.tsv")));
        assert_eq!(cfg.entropy_hist_out, Some(PathBuf::from("entropy.tsv")));
        assert_eq!(cfg.identity_hist_out, Some(PathBuf::from("id.tsv")));
        assert_eq!(cfg.gc_bins, Some(100));
        assert_eq!(cfg.entropy_bins, 100);
        assert_eq!(cfg.identity_bins, 100);
        assert_eq!(cfg.side_hist_len, Some(1000));

        for case in [
            "gcbins=abc",
            "entropybins=abc",
            "idhistbins=abc",
            "maxhistlen=abc",
            "maxhistlen=0",
        ] {
            let err = parse_args(
                ["in=reads.fq", "passes=1", case]
                    .into_iter()
                    .map(OsString::from),
            )
            .unwrap_err()
            .to_string();
            assert!(
                err.contains("expects") || err.contains("suffix"),
                "unexpected error for malformed {case}: {err}"
            );
        }
    }

    #[test]
    fn accepts_cardinality_loglog_controls_as_bounded_estimates_and_validates_values() {
        for case in [
            "cardinality=t",
            "cardinality=31",
            "loglog=f",
            "loglogin=t",
            "cardinalityout=t",
            "loglogout=f",
            "buckets=1k",
            "loglogbuckets=100",
            "loglogcorrection=t",
            "loglogcf=f",
            "loglogbits=16",
            "loglogk=31",
            "cardinalityk=31",
            "kcardinality=31",
            "loglogklist=21,31",
            "loglogseed=42",
            "loglogminprob=0.5",
            "loglogtype=loglog2",
            "loglogmean=t",
            "loglogmedian=t",
            "loglogmwa=t",
            "logloghmean=t",
            "logloggmean=t",
            "loglogmantissa=8",
            "loglogcounts=t",
            "loglogcount=f",
        ] {
            let cfg = parse_args(
                ["in=reads.fq", "passes=1", case]
                    .into_iter()
                    .map(OsString::from),
            )
            .unwrap();
            assert!(
                cfg.notes
                    .iter()
                    .any(|note| note.contains("cardinality/loglog")),
                "expected cardinality/loglog fallback note for {case}: {:?}",
                cfg.notes
            );
        }

        let cfg = parse(&[
            "in=reads.fq",
            "passes=1",
            "cardinality=t",
            "cardinalityout=t",
            "buckets=1k",
            "loglogseed=42",
            "loglogk=25",
            "loglogminprob=0.25",
        ]);
        assert!(cfg.cardinality.input);
        assert!(cfg.cardinality.output);
        assert_eq!(cfg.cardinality.buckets, 1000);
        assert_eq!(cfg.cardinality.seed, 42);
        assert_eq!(cfg.cardinality.k, Some(25));
        assert_eq!(cfg.cardinality.min_probability, 0.25);

        let cfg = parse(&[
            "in=reads.fq",
            "passes=1",
            "cardinality=t",
            "cardinality=f",
            "cardinalityout=t",
            "loglogout=f",
        ]);
        assert!(!cfg.cardinality.input);
        assert!(!cfg.cardinality.output);

        for case in [
            "cardinality=maybe",
            "buckets=0",
            "buckets=100g",
            "loglogbits=abc",
            "loglogklist=21,abc",
            "loglogseed=abc",
            "loglogminprob=abc",
            "loglogminprob=2",
        ] {
            let err = parse_args(
                ["in=reads.fq", "passes=1", case]
                    .into_iter()
                    .map(OsString::from),
            )
            .unwrap_err()
            .to_string();
            assert!(
                err.contains("expects") || err.contains("above the Rust safety cap"),
                "unexpected error for malformed {case}: {err}"
            );
        }
    }

    #[test]
    fn accepts_quality_recalibration_controls_as_noops_and_validates_values() {
        let cfg = parse(&[
            "in=reads.fq",
            "trackall=f",
            "clearmatrices=f",
            "loadq=f",
            "loadq102=f",
            "loadqap=f",
            "loadqbp=f",
            "loadqpt=f",
            "loadqbt=f",
            "loadq10=f",
            "loadq12=f",
            "loadqb12=f",
            "loadqb012=f",
            "loadqb123=f",
            "loadqb234=f",
            "loadq12b12=f",
            "loadqp=f",
            "observationcutoff=1k",
            "recalpasses=1",
            "recalqmax=50",
            "recalqmin=2",
            "recalwithposition=t",
            "qmatrixmode=max",
            "recaltile=f",
        ]);
        assert_eq!(cfg.notes.len(), 23);

        let cfg = parse(&[
            "in=reads.fq",
            "loadq102_p1=f",
            "loadq_p2=t",
            "observationcutoff_p1=1k",
            "recalpasses_p2=1",
            "recalqmax_p1=50",
            "recalqmin_p2=2",
            "recalwithposition_p1=t",
            "qmatrixmode_p2=max",
            "recaltile_p1=f",
        ]);
        assert_eq!(cfg.notes.len(), 9);

        for case in [
            "observationcutoff=abc",
            "recalpasses=abc",
            "recalqmax=abc",
            "observationcutoff_p1=abc",
            "recalpasses_p2=abc",
            "recalqmax_p1=abc",
        ] {
            let err = parse_args(
                ["in=reads.fq", "passes=1", case]
                    .into_iter()
                    .map(OsString::from),
            )
            .unwrap_err()
            .to_string();
            assert!(
                err.contains("expects") || err.contains("suffix"),
                "unexpected error for malformed {case}: {err}"
            );
        }
    }

    #[test]
    fn accepts_disabled_recalibrate_controls_and_rejects_enabled_recalibration() {
        let cfg = parse(&[
            "in=reads.fq",
            "recalibrate=f",
            "recalibratequality=f",
            "recal=f",
            "recalibrate_p1=f",
        ]);
        assert_eq!(cfg.notes.len(), 4);
        assert!(
            cfg.notes
                .iter()
                .all(|note| note.contains("keeps BBTools quality recalibration disabled"))
        );

        for case in ["recalibrate=t", "recalibratequality=t", "recal=t"] {
            let err = parse_args(
                ["in=reads.fq", "passes=1", case]
                    .into_iter()
                    .map(OsString::from),
            )
            .unwrap_err()
            .to_string();
            assert!(
                err.contains("enables BBTools quality recalibration"),
                "unexpected error for enabled {case}: {err}"
            );
        }

        let err = parse_args(
            ["in=reads.fq", "passes=1", "recalibrate=maybe"]
                .into_iter()
                .map(OsString::from),
        )
        .unwrap_err()
        .to_string();
        assert!(
            err.contains("recalibrate expects a boolean value"),
            "unexpected error for malformed recalibrate: {err}"
        );
    }

    #[test]
    fn accepts_disabled_break_length_controls_and_rejects_read_splitting() {
        let cfg = parse(&["in=reads.fq", "breaklen=0", "breaklength=-1"]);
        assert_eq!(cfg.notes.len(), 2);
        assert!(
            cfg.notes
                .iter()
                .all(|note| note.contains("keeps BBTools read breaking disabled"))
        );

        for case in ["breaklen=50", "breaklength=1"] {
            let err = parse_args(
                ["in=reads.fq", "passes=1", case]
                    .into_iter()
                    .map(OsString::from),
            )
            .unwrap_err()
            .to_string();
            assert!(
                err.contains("enables BBTools read breaking"),
                "unexpected error for enabled {case}: {err}"
            );
        }

        let err = parse_args(
            ["in=reads.fq", "passes=1", "breaklen=abc"]
                .into_iter()
                .map(OsString::from),
        )
        .unwrap_err()
        .to_string();
        assert!(
            err.contains("breaklen expects"),
            "unexpected error for malformed breaklen: {err}"
        );
    }

    #[test]
    fn accepts_shared_environment_runtime_controls_as_noops_and_validates_values() {
        let cfg = parse(&[
            "in=reads.fq",
            "amino=f",
            "amino8=f",
            "validatebranchless=maybe",
            "fairqueues=t",
            "fixextensions=f",
            "2passresize=f",
            "parallelsort=f",
            "gcbeforemem=t",
            "warnifnosequence=f",
            "warnfirsttimeonly=f",
            "kmg=t",
            "forceJavaParseDouble=f",
            "simd=auto",
            "simdsparse=f",
            "simdmultsparse=f",
            "simdfmasparse=f",
            "simdcopy=f",
            "aws=f",
            "nersc=t",
            "lowmem=f",
            "lockedincrement=auto",
            "symmetricwrite=f",
            "buffer=10",
            "buffered=f",
            "sidechannelstats=f",
            "silva=f",
            "unite=f",
            "imghq=f",
            "callins=f",
            "calldel=f",
            "callsub=f",
            "callindel=f",
            "calljunct=f",
            "callnocall=f",
            "protFull=t",
            "entropyk=3",
            "entropywindow=50",
        ]);
        assert_eq!(cfg.notes.len(), 37);
        assert_eq!(cfg.locked_increment, Some(false));

        for case in ["entropyk=abc", "entropywindow=abc"] {
            let err = parse_args(
                ["in=reads.fq", "passes=1", case]
                    .into_iter()
                    .map(OsString::from),
            )
            .unwrap_err()
            .to_string();
            assert!(
                err.contains("expects"),
                "unexpected error for malformed {case}: {err}"
            );
        }

        for case in ["amino=t", "amino8=t"] {
            let err = parse_args(
                ["in=reads.fq", "passes=1", case]
                    .into_iter()
                    .map(OsString::from),
            )
            .unwrap_err()
            .to_string();
            assert!(
                err.contains("amino-acid kmer mode"),
                "unexpected error for enabled {case}: {err}"
            );
        }
    }

    #[test]
    fn parses_base_cleanup_controls_like_bbnorm() {
        let cfg = parse(&[
            "in=reads.fq",
            "utot=t",
            "tuc=t",
            "lctn=t",
            "dotdashxton=t",
            "itn=t",
            "fixjunk=t",
        ]);
        assert!(cfg.u_to_t);
        assert!(cfg.to_upper_case);
        assert!(cfg.lower_case_to_n);
        assert!(cfg.dot_dash_x_to_n);
        assert!(cfg.iupac_to_n);
        assert_eq!(cfg.junk_mode, JunkMode::Fix);

        let cfg = parse(&["in=reads.fq", "ignorejunk=t"]);
        assert_eq!(cfg.junk_mode, JunkMode::Ignore);

        let cfg = parse(&["in=reads.fq", "flagjunk=t"]);
        assert_eq!(cfg.junk_mode, JunkMode::Flag);

        let cfg = parse(&["in=reads.fq", "tossjunk=t"]);
        assert_eq!(cfg.junk_mode, JunkMode::Flag);

        let cfg = parse(&["in=reads.fq", "junk=discard"]);
        assert_eq!(cfg.junk_mode, JunkMode::Flag);

        let cfg = parse(&["in=reads.fq", "crashjunk=f"]);
        assert_eq!(cfg.junk_mode, JunkMode::Ignore);

        let cfg = parse(&["in=reads.fq", "failjunk=f"]);
        assert_eq!(cfg.junk_mode, JunkMode::Ignore);

        let cfg = parse(&["in=reads.fq", "ignorejunk=t", "crashjunk=t"]);
        assert_eq!(cfg.junk_mode, JunkMode::Crash);

        let cfg = parse(&["in=reads.fq", "junk=fail"]);
        assert_eq!(cfg.junk_mode, JunkMode::Crash);

        let cfg = parse(&["in=reads.fq", "junk=iupacton"]);
        assert!(cfg.fix_junk_and_iupac);
        assert_eq!(cfg.junk_mode, JunkMode::Fix);
    }

    #[test]
    fn false_flagjunk_alias_resets_to_crash_like_bbnorm() {
        let cfg = parse(&["in=reads.fq", "flagjunk=t", "flagjunk=f"]);
        assert_eq!(cfg.junk_mode, JunkMode::Crash);

        let cfg = parse(&["in=reads.fq", "tossjunk=t", "tossjunk=f"]);
        assert_eq!(cfg.junk_mode, JunkMode::Flag);
    }

    #[test]
    fn accepts_bbnorm_inactive_trim_parser_options_as_noops() {
        let cfg = parse(&[
            "in=reads.fq",
            "trimclip=t",
            "trimpolya=t",
            "trimpolyg=10",
            "trimpolygleft=f",
            "trimpolycright=2",
            "maxnonpoly=3",
            "ftr=10",
            "ftl=2",
            "ftm=4",
            "ftr2=7",
        ]);
        assert_eq!(cfg.notes.len(), 10);
    }

    #[test]
    fn accepts_bbnorm_inactive_read_filter_parser_options_as_noops() {
        let cfg = parse(&[
            "in=reads.fq",
            "maxlen=50",
            "minlenfraction=0.8",
            "maxns=0",
            "mingc=0.9",
            "maxgc=0.1",
            "usepairgc=t",
            "minconsecutivebases=200",
            "maq=40,20",
            "maqb=20",
            "mbq=30",
            "chastityfilter=t",
            "trimbadsequence=t",
            "failnobarcode=f",
            "badbarcodes=fail",
            "barcodefilter=f",
            "barcodes=ACGT,TGCA",
            "aqbp=t",
            "mintrimlen=10",
            "untrim=f",
        ]);
        assert_eq!(cfg.notes.len(), 19);

        for case in ["mintrimlen=abc", "badbarcodes=maybe"] {
            let err = parse_args(
                ["in=reads.fq", "passes=1", case]
                    .into_iter()
                    .map(OsString::from),
            )
            .unwrap_err()
            .to_string();
            assert!(
                err.contains("expects"),
                "unexpected error for malformed {case}: {err}"
            );
        }
    }

    #[test]
    fn accepts_genome_build_context_controls_as_normalization_noops() {
        for case in ["build=1", "genome=1"] {
            let cfg = parse_args(
                ["in=reads.fq", "passes=1", case]
                    .into_iter()
                    .map(OsString::from),
            )
            .unwrap();
            assert!(
                cfg.notes
                    .iter()
                    .any(|note| note.contains("genome-build context")),
                "expected genome-build context no-op note for {case}: {:?}",
                cfg.notes
            );
        }

        for case in ["genome=abc", "idfilter=0.9", "subfilter=1"] {
            let err = parse_args(
                ["in=reads.fq", "passes=1", case]
                    .into_iter()
                    .map(OsString::from),
            )
            .unwrap_err()
            .to_string();
            assert!(
                err.contains("expects") || err.contains("unknown or unsupported"),
                "unexpected error for malformed {case}: {err}"
            );
        }
    }

    #[test]
    fn parses_explicit_interleaved_single_stream_outputs() {
        let cfg = parse(&[
            "in=reads.fq",
            "interleaved=t",
            "out=keep.fq",
            "outt=toss.fq",
        ]);
        assert!(cfg.interleaved);
        assert_eq!(cfg.in1.unwrap(), PathBuf::from("reads.fq"));
        assert_eq!(cfg.out1.unwrap(), PathBuf::from("keep.fq"));
        assert_eq!(cfg.out_toss1.unwrap(), PathBuf::from("toss.fq"));
        assert!(cfg.out2.is_none());
        assert!(cfg.out_toss2.is_none());

        let cfg = parse(&["in=reads.fq", "int=t"]);
        assert!(cfg.interleaved);
        assert!(!cfg.test_interleaved);

        let cfg = parse(&["in=reads.fq", "forceinterleaved=t"]);
        assert!(cfg.interleaved);
        assert!(!cfg.test_interleaved);

        let cfg = parse(&["in=reads.fq", "testinterleaved=f"]);
        assert!(!cfg.interleaved);
        assert!(!cfg.test_interleaved);

        let cfg = parse(&["in=reads.fq", "overrideinterleaved=t"]);
        assert!(!cfg.notes.is_empty());
    }

    #[test]
    fn defaults_to_auto_interleaved_detection() {
        let cfg = parse(&[
            "in=reads.fq",
            "out=keep1.fq",
            "out2=keep2.fq",
            "outt=toss1.fq",
            "outt2=toss2.fq",
        ]);
        assert!(!cfg.interleaved);
        assert!(cfg.test_interleaved);
    }

    #[test]
    fn paired_input_allows_bbnorm_single_stream_or_hash_pattern_outputs() {
        let cfg = parse(&[
            "in=reads1.fq",
            "in2=reads2.fq",
            "out=keep#.fq",
            "outt=toss.fq",
        ]);
        assert_eq!(cfg.out1.unwrap(), PathBuf::from("keep#.fq"));
        assert!(cfg.out2.is_none());
        assert_eq!(cfg.out_toss1.unwrap(), PathBuf::from("toss.fq"));
        assert!(cfg.out_toss2.is_none());
    }

    #[test]
    fn interleaved_true_with_in2_remains_two_file_paired_like_bbnorm() {
        let cfg = parse(&["in=reads1.fq", "in2=reads2.fq", "interleaved=t"]);
        assert!(cfg.interleaved);
        assert_eq!(cfg.in1.unwrap(), PathBuf::from("reads1.fq"));
        assert_eq!(cfg.in2.unwrap(), PathBuf::from("reads2.fq"));
    }

    #[test]
    fn expands_missing_hash_input_pattern_like_bbnorm() {
        let cfg = parse(&["in=reads#.fq"]);
        assert_eq!(cfg.in1.unwrap(), PathBuf::from("reads1.fq"));
        assert_eq!(cfg.in2.unwrap(), PathBuf::from("reads2.fq"));
    }

    #[test]
    fn keeps_literal_hash_input_when_file_exists_like_bbnorm() {
        let dir = tempfile::tempdir().unwrap();
        let literal = dir.path().join("reads#.fq");
        std::fs::write(&literal, b"@r1\nACGT\n+\nIIII\n").unwrap();

        let cfg = parse_args(
            [format!("in={}", literal.display()), "passes=1".to_string()]
                .into_iter()
                .map(OsString::from),
        )
        .unwrap();
        assert_eq!(cfg.in1.unwrap(), literal);
        assert!(cfg.in2.is_none());
    }

    #[test]
    fn keeps_literal_comma_extra_when_file_exists_like_bbnorm() {
        let dir = tempfile::tempdir().unwrap();
        let input = dir.path().join("main.fq");
        let literal = dir.path().join("extra,with,commas.fq");
        std::fs::write(&input, b"@r1\nACGT\n+\nIIII\n").unwrap();
        std::fs::write(&literal, b"@r2\nACGT\n+\nIIII\n").unwrap();

        let cfg = parse_args(
            [
                format!("in={}", input.display()),
                format!("extra={}", literal.display()),
                "extra=null".to_string(),
                "passes=1".to_string(),
            ]
            .into_iter()
            .map(OsString::from),
        )
        .unwrap();
        assert_eq!(cfg.extra, vec![literal]);
    }

    #[test]
    fn expands_config_files_like_bbnorm() {
        let dir = tempfile::tempdir().unwrap();
        let cfg1 = dir.path().join("a.config");
        let cfg2 = dir.path().join("b.config");
        std::fs::write(
            &cfg1,
            "\n# comment\nin=reads.fq\npasses=1\nkeepall=t\nk=21\n",
        )
        .unwrap();
        std::fs::write(&cfg2, "target=7\nout=keep.fq\n").unwrap();

        let cfg = parse_args(
            [
                format!("config={},{}", cfg1.display(), cfg2.display()),
                "target=9".to_string(),
            ]
            .into_iter()
            .map(OsString::from),
        )
        .unwrap();

        assert_eq!(cfg.in1.unwrap(), PathBuf::from("reads.fq"));
        assert_eq!(cfg.k, 21);
        assert_eq!(cfg.target_depth, 9);
        assert_eq!(cfg.out1.unwrap(), PathBuf::from("keep.fq"));
        assert!(cfg.keep_all);
        assert!(
            cfg.notes
                .iter()
                .any(|note| note.contains("expanded into 6 BBTools-style argument line"))
        );
    }

    #[test]
    fn reports_missing_config_files_like_bbnorm() {
        let dir = tempfile::tempdir().unwrap();
        let missing = dir.path().join("missing.config");
        let err = parse_args(
            [format!("config={}", missing.display())]
                .into_iter()
                .map(OsString::from),
        )
        .unwrap_err()
        .to_string();

        assert!(err.contains("could not process config file"));
    }

    #[test]
    fn rejects_missing_extra_inputs_like_bbnorm() {
        let err = parse_args(
            ["in=reads.fq", "extra=missing#.fq", "passes=1"]
                .into_iter()
                .map(OsString::from),
        )
        .unwrap_err();
        assert!(err.to_string().contains("extra input missing#.fq"));
    }

    #[test]
    fn parses_single_pass_final_stage_aliases() {
        let cfg = parse(&[
            "in=reads.fq",
            "tbrf=t",
            "dbo2=t",
            "tossbadreads1=t",
            "dbo1=t",
        ]);
        assert!(cfg.toss_error_reads);
        assert!(cfg.toss_error_reads_first);
        assert!(cfg.discard_bad_only);
        assert!(cfg.discard_bad_only_first);
    }

    #[test]
    fn parses_multipass_and_countup_controls() {
        let cfg = parse(&[
            "in=reads.fq",
            "passes=2",
            "target1=7",
            "targetbadpercentilelow=20",
            "tbph=0.8",
            "abrc=t",
        ]);
        assert_eq!(cfg.target_depth_first, Some(7));
        assert_eq!(cfg.target_bad_percent_low, 0.2);
        assert_eq!(cfg.target_bad_percent_high, 0.8);
        assert!(cfg.add_bad_reads_countup);

        for case in ["target1=abc", "targetbadpercentilelow=abc", "tbph=abc"] {
            let err = parse_args(
                ["in=reads.fq", "passes=1", case]
                    .into_iter()
                    .map(OsString::from),
            )
            .unwrap_err()
            .to_string();
            assert!(
                err.contains("expects"),
                "unexpected error for malformed {case}: {err}"
            );
        }
    }

    #[test]
    fn allows_outuncorrected_in_multipass_runs() {
        let cfg = parse(&[
            "in=reads_1.fq",
            "in2=reads_2.fq",
            "passes=2",
            "out=keep_1.fq",
            "out2=keep_2.fq",
            "outuncorrected=unc_1.fq",
            "outuncorrected2=unc_2.fq",
        ]);
        assert_eq!(cfg.passes, 2);
        assert_eq!(
            cfg.out_uncorrected1.as_deref(),
            Some(std::path::Path::new("unc_1.fq"))
        );
        assert_eq!(
            cfg.out_uncorrected2.as_deref(),
            Some(std::path::Path::new("unc_2.fq"))
        );
    }

    #[test]
    fn final_stage_alias_can_override_conflated_alias() {
        let cfg = parse(&["in=reads.fq", "tossbadreads=t", "tossbadreadsf=f"]);
        assert!(!cfg.toss_error_reads);
    }

    #[test]
    fn remove_if_either_bad_alias_inverts_require_both_bad() {
        let cfg = parse(&["in=reads.fq", "requirebothbad=t", "removeifeitherbad=t"]);
        assert!(!cfg.require_both_bad);

        let cfg = parse(&["in=reads.fq", "rieb=f"]);
        assert!(cfg.require_both_bad);
    }

    #[test]
    fn explicit_interleaved_false_rejects_second_outputs_without_in2() {
        let err = parse_args(
            [
                "in=reads.fq",
                "interleaved=f",
                "out=keep1.fq",
                "out2=keep2.fq",
                "passes=1",
            ]
            .into_iter()
            .map(OsString::from),
        )
        .unwrap_err();
        assert!(err.to_string().contains("out2"));
    }

    #[test]
    fn enabled_ecc_sets_real_correction_fields() {
        let cfg = parse_args(["in=x.fq", "ecc=t"].into_iter().map(OsString::from)).unwrap();
        assert_eq!(cfg.passes, 2);
        assert!(cfg.error_correct);
        assert!(cfg.error_correct_first);
        assert!(cfg.error_correct_final);
        assert!(!cfg.overlap_error_correct);
        assert!(!cfg.mark_errors_only);
        assert!(cfg.notes.is_empty());

        let cfg = parse(&["in=x.fq", "ecc=f"]);
        assert!(!cfg.error_correct);
        assert!(!cfg.error_correct_first);
        assert!(!cfg.error_correct_final);
        assert!(!cfg.overlap_error_correct);

        let cfg = parse(&["in=x.fq", "ecc1=t", "ecc2=f"]);
        assert!(cfg.error_correct);
        assert!(cfg.error_correct_first);
        assert!(!cfg.error_correct_final);

        let cfg = parse(&["in=x.fq", "ecc1=f", "eccf=t"]);
        assert!(cfg.error_correct);
        assert!(!cfg.error_correct_first);
        assert!(cfg.error_correct_final);

        let cfg = parse(&["in=x.fq", "markerrors=t"]);
        assert!(cfg.error_correct);
        assert!(cfg.error_correct_first);
        assert!(!cfg.error_correct_final);

        let cfg = parse(&["in=x.fq", "ecco=t"]);
        assert!(cfg.error_correct);
        assert!(cfg.error_correct_first);
        assert!(cfg.error_correct_final);
        assert!(cfg.overlap_error_correct);
        assert!(!cfg.overlap_error_correct_auto);
        assert!(cfg.notes[0].contains("paired overlap repair"));

        let cfg = parse(&["in=x.fq", "ecco=auto"]);
        assert!(cfg.error_correct);
        assert!(cfg.error_correct_first);
        assert!(cfg.error_correct_final);
        assert!(!cfg.overlap_error_correct);
        assert!(cfg.overlap_error_correct_auto);
        assert!(cfg.notes[0].contains("automatic overlap"));

        let cfg = parse(&["in=x.fq", "ecco=t", "ecco=f"]);
        assert!(cfg.error_correct);
        assert!(!cfg.overlap_error_correct);
        assert!(!cfg.overlap_error_correct_auto);
    }

    #[test]
    fn accepts_ecc_tuning_controls_and_validates_integers() {
        let cfg = parse(&[
            "in=reads.fq",
            "ecclimit=3",
            "eccmaxqual=127",
            "errorcorrectratio=140",
            "echighthresh=22",
            "eclowthresh=2",
            "suflen=3",
            "prefixlen=3",
            "cfl=t",
            "cfr=f",
        ]);
        assert_eq!(cfg.max_errors_to_correct, 3);
        assert_eq!(cfg.max_quality_to_correct, 127);
        assert_eq!(cfg.error_correct_ratio, 140);
        assert_eq!(cfg.error_correct_high_thresh, 22);
        assert_eq!(cfg.error_correct_low_thresh, 2);
        assert_eq!(cfg.suffix_len, 3);
        assert_eq!(cfg.prefix_len, 3);
        assert!(cfg.correct_from_left);
        assert!(!cfg.correct_from_right);
        assert!(cfg.notes.is_empty());

        for case in [
            "ecclimit=abc",
            "eccmaxqual=abc",
            "ecr=abc",
            "echthresh=abc",
            "eclt=abc",
            "suflen=abc",
            "prelen=abc",
        ] {
            let err = parse_args(
                ["in=reads.fq", "passes=1", case]
                    .into_iter()
                    .map(OsString::from),
            )
            .unwrap_err()
            .to_string();
            assert!(
                err.contains("expects"),
                "unexpected error for malformed {case}: {err}"
            );
        }
    }

    #[test]
    fn parses_countup_mode() {
        let cfg = parse_args(["in=x.fq", "countup=t"].into_iter().map(OsString::from)).unwrap();
        assert!(cfg.count_up);

        let cfg = parse(&["in=x.fq", "countup=f"]);
        assert!(!cfg.count_up);
        assert!(cfg.notes.iter().any(|note| note.contains("countup=f")));
    }

    #[test]
    fn wrapper_sampling_options_fall_back_to_supported_normalization() {
        for case in [
            "sampleoutput=1",
            "readsample=1",
            "kmersample=1",
            "samplerate=0.5",
            "sample=0.5",
            "sampleseed=1",
            "seed=1",
        ] {
            let cfg = parse_args(
                ["in=x.fq", "passes=1", case]
                    .into_iter()
                    .map(OsString::from),
            )
            .unwrap();
            assert!(
                cfg.notes
                    .iter()
                    .any(|note| note.contains("Rust ignores it")),
                "expected sampling fallback note for {case}"
            );
        }
    }

    #[test]
    fn nondeterministic_mode_stays_enabled_for_random_selection() {
        for case in ["deterministic=t", "dr=t", "det=t"] {
            let cfg = parse(&["in=reads.fq", case]);
            assert!(cfg.deterministic, "expected deterministic mode for {case}");
        }

        let cfg = parse_args(
            ["in=reads.fq", "passes=1", "deterministic=f"]
                .into_iter()
                .map(OsString::from),
        )
        .unwrap();
        assert!(!cfg.deterministic);
        assert!(
            cfg.notes
                .iter()
                .all(|note| !note.contains("deterministic=f is not implemented yet"))
        );
        assert!(
            cfg.notes
                .iter()
                .any(|note| note.contains("faster parallel replay"))
        );
    }
}

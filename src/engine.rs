use crate::cli::{CARDINALITY_MAX_BUCKETS, Config};
use crate::kmer::{
    KmerKey, canonical_short_code, for_each_kmer_for_record, unfiltered_kmer_windows_for_record,
};
use crate::peaks::write_peaks;
use crate::seqio::{
    BaseSettings, QualitySettings, SeqFormat, SequenceReader, SequenceRecord, SequenceSettings,
    SequenceWriter, create_output_with_append, detect_interleaved_input_with_gzip_threads,
};
use anyhow::{Context, Result, bail, ensure};
use rayon::prelude::*;
use rustc_hash::FxHashMap;
use std::alloc::{Layout, alloc_zeroed};
use std::cmp::Ordering as CmpOrdering;
use std::collections::{BTreeMap, BinaryHeap};
use std::fs;
use std::io::{BufReader, BufWriter, ErrorKind, Read, Write};
use std::path::{Path, PathBuf};
use std::process::{Child, ChildStdin, ChildStdout, Command, Stdio};
use std::sync::{
    Mutex, OnceLock,
    atomic::{AtomicU32, AtomicU64, AtomicUsize, Ordering},
};
use std::time::{Instant, SystemTime, UNIX_EPOCH};

pub type CountMap = FxHashMap<KmerKey, u64>;

#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub struct CardinalityEstimate {
    pub k: usize,
    pub buckets: usize,
    pub estimated_unique_kmers: u64,
}

struct KmerCardinalityEstimator {
    k: usize,
    buckets: usize,
    seed: u64,
    registers: Vec<u8>,
}

impl KmerCardinalityEstimator {
    fn from_config(config: &Config) -> Self {
        let buckets = config.cardinality.buckets.clamp(1, CARDINALITY_MAX_BUCKETS);
        Self {
            k: config.cardinality.k.unwrap_or(config.k),
            buckets,
            seed: config.cardinality.seed,
            registers: vec![0; buckets],
        }
    }

    fn observe_pair(&mut self, config: &Config, r1: &SequenceRecord, r2: Option<&SequenceRecord>) {
        self.observe_record(config, r1);
        if let Some(mate) = r2 {
            self.observe_record(config, mate);
        }
    }

    fn observe_record(&mut self, config: &Config, record: &SequenceRecord) {
        for_each_kmer_for_record(record, config, |kmer| self.observe_key(&kmer));
    }

    fn observe_key(&mut self, key: &KmerKey) {
        let raw = raw_kmer_key(key);
        let kind_salt = match key {
            KmerKey::Short(_) => 0x9E37_79B9_7F4A_7C15,
            KmerKey::LongHash(_) => 0xD1B5_4A32_D192_ED03,
        };
        let hash = mix_seed(raw ^ self.seed ^ kind_salt);
        let bucket = (((hash as u128) * (self.buckets as u128)) >> 64) as usize;
        let rank_hash = mix_seed(hash ^ 0x94D0_49BB_1331_11EB);
        let rank = rank_hash.leading_zeros().saturating_add(1).min(64) as u8;
        if let Some(slot) = self.registers.get_mut(bucket) {
            *slot = (*slot).max(rank);
        }
    }

    fn estimate(&self) -> CardinalityEstimate {
        let m = self.buckets as f64;
        let zero_count = self
            .registers
            .iter()
            .filter(|&&register| register == 0)
            .count();
        let inverse_sum: f64 = self
            .registers
            .iter()
            .map(|&register| 2f64.powi(-(i32::from(register))))
            .sum();
        let raw_estimate = hll_alpha(self.buckets) * m * m / inverse_sum.max(f64::MIN_POSITIVE);
        let corrected = if raw_estimate <= 2.5 * m && zero_count > 0 {
            m * (m / zero_count as f64).ln()
        } else {
            raw_estimate
        };
        CardinalityEstimate {
            k: self.k,
            buckets: self.buckets,
            estimated_unique_kmers: corrected.round().max(0.0) as u64,
        }
    }
}

fn hll_alpha(buckets: usize) -> f64 {
    match buckets {
        16 => 0.673,
        32 => 0.697,
        64 => 0.709,
        _ => 0.7213 / (1.0 + 1.079 / buckets as f64),
    }
}

trait CountLookup: Sync {
    fn depth(&self, key: &KmerKey) -> u64;
    fn unique_kmers(&self) -> usize;
    fn unique_kmers_at_least(&self, min_depth: u64) -> usize;
}

impl CountLookup for CountMap {
    fn depth(&self, key: &KmerKey) -> u64 {
        self.get(key).copied().unwrap_or(0)
    }

    fn unique_kmers(&self) -> usize {
        self.len()
    }

    fn unique_kmers_at_least(&self, min_depth: u64) -> usize {
        if min_depth <= 1 {
            return self.len();
        }
        self.values().filter(|&&depth| depth >= min_depth).count()
    }
}

enum InputCounts {
    Exact(CountMap),
    Sketch(PackedCountMinSketch),
    AtomicSketch(AtomicCountMinSketch),
    PrefilteredSketch {
        prefilter: PrefilterCountMinSketch,
        limit: u64,
        main: Box<InputCounts>,
    },
}

#[derive(Clone, Copy)]
struct PrefilterGate<'a> {
    sketch: &'a PrefilterCountMinSketch,
    limit: u64,
}

impl<'a> PrefilterGate<'a> {
    fn new(sketch: &'a PrefilterCountMinSketch, limit: u64) -> Self {
        Self {
            sketch,
            limit: limit.min(sketch.max_count()),
        }
    }

    fn should_count_in_main(&self, key: &KmerKey) -> bool {
        self.sketch.depth(key) >= self.limit
    }
}

impl CountLookup for InputCounts {
    fn depth(&self, key: &KmerKey) -> u64 {
        match self {
            Self::Exact(counts) => counts.depth(key),
            Self::Sketch(sketch) => sketch.depth(key),
            Self::AtomicSketch(sketch) => sketch.depth(key),
            Self::PrefilteredSketch {
                prefilter,
                limit,
                main,
            } => {
                let prefilter_depth = prefilter.depth(key);
                if prefilter_depth < *limit {
                    prefilter_depth
                } else {
                    main.depth(key)
                }
            }
        }
    }

    fn unique_kmers(&self) -> usize {
        match self {
            Self::Exact(counts) => counts.unique_kmers(),
            Self::Sketch(sketch) => sketch.unique_kmers(),
            Self::AtomicSketch(sketch) => sketch.unique_kmers(),
            Self::PrefilteredSketch { prefilter, .. } => prefilter.unique_kmers(),
        }
    }

    fn unique_kmers_at_least(&self, min_depth: u64) -> usize {
        match self {
            Self::Exact(counts) => counts.unique_kmers_at_least(min_depth),
            Self::Sketch(sketch) => sketch.unique_kmers_at_least(min_depth),
            Self::AtomicSketch(sketch) => sketch.unique_kmers_at_least(min_depth),
            Self::PrefilteredSketch {
                prefilter,
                limit,
                main,
            } => {
                if min_depth < *limit {
                    prefilter.unique_kmers_at_least(min_depth)
                } else {
                    main.unique_kmers_at_least(min_depth)
                }
            }
        }
    }
}

impl InputCounts {
    #[cfg(test)]
    fn unique_kmer_estimate_split(&self) -> Option<UniqueKmerEstimateSplit> {
        self.unique_kmer_estimate().1
    }

    fn unique_kmer_estimate(&self) -> (usize, Option<UniqueKmerEstimateSplit>) {
        match self {
            Self::PrefilteredSketch {
                prefilter, main, ..
            } => {
                let low_depth_max = prefilter.max_count();
                let high_depth_min = low_depth_max.saturating_add(1);
                let total = prefilter.unique_kmers();
                let high_depth_kmers = main.unique_kmers_at_least(high_depth_min);
                (
                    total,
                    Some(UniqueKmerEstimateSplit {
                        low_depth_max,
                        low_depth_kmers: total.saturating_sub(high_depth_kmers),
                        high_depth_min,
                        high_depth_kmers,
                    }),
                )
            }
            _ => (self.unique_kmers(), None),
        }
    }

    fn sketch_layouts(&self) -> Vec<SketchLayoutSummary> {
        let mut layouts = Vec::new();
        self.append_sketch_layouts(&mut layouts, "input_main");
        layouts
    }

    fn append_sketch_layouts(&self, layouts: &mut Vec<SketchLayoutSummary>, table: &'static str) {
        match self {
            Self::Exact(_) => {}
            Self::Sketch(sketch) => layouts.push(sketch.layout_summary(table, None)),
            Self::AtomicSketch(sketch) => layouts.push(sketch.layout_summary(table, None)),
            Self::PrefilteredSketch {
                prefilter,
                limit,
                main,
            } => {
                layouts.push(prefilter.layout_summary("input_prefilter", Some(*limit)));
                main.append_sketch_layouts(layouts, "input_main");
            }
        }
    }
}

enum OutputCounts {
    Exact(CountMap),
    Sketch(PackedCountMinSketch),
    AtomicSketch(AtomicCountMinSketch),
}

impl CountLookup for OutputCounts {
    fn depth(&self, key: &KmerKey) -> u64 {
        match self {
            Self::Exact(counts) => counts.depth(key),
            Self::Sketch(sketch) => sketch.depth(key),
            Self::AtomicSketch(sketch) => sketch.depth(key),
        }
    }

    fn unique_kmers(&self) -> usize {
        match self {
            Self::Exact(counts) => counts.unique_kmers(),
            Self::Sketch(sketch) => sketch.unique_kmers(),
            Self::AtomicSketch(sketch) => sketch.unique_kmers(),
        }
    }

    fn unique_kmers_at_least(&self, min_depth: u64) -> usize {
        match self {
            Self::Exact(counts) => counts.unique_kmers_at_least(min_depth),
            Self::Sketch(sketch) => sketch.unique_kmers_at_least(min_depth),
            Self::AtomicSketch(sketch) => sketch.unique_kmers_at_least(min_depth),
        }
    }
}

impl OutputCounts {
    #[cfg(test)]
    #[cfg(test)]
    #[cfg(test)]
    fn depth_hist(&self, hist_len: usize) -> Vec<u64> {
        match self {
            Self::Exact(counts) => count_map_depth_hist(counts, hist_len),
            Self::Sketch(sketch) => sketch.depth_hist(hist_len),
            Self::AtomicSketch(sketch) => sketch.depth_hist(hist_len),
        }
    }

    fn sparse_depth_hist(&self, hist_len: usize) -> SparseHist {
        match self {
            Self::Exact(counts) => count_map_sparse_depth_hist(counts, hist_len),
            Self::Sketch(sketch) => sketch.sparse_depth_hist(hist_len),
            Self::AtomicSketch(sketch) => sketch.sparse_depth_hist(hist_len),
        }
    }

    fn append_sketch_layouts(&self, layouts: &mut Vec<SketchLayoutSummary>, table: &'static str) {
        match self {
            Self::Exact(_) => {}
            Self::Sketch(sketch) => layouts.push(sketch.layout_summary(table, None)),
            Self::AtomicSketch(sketch) => layouts.push(sketch.layout_summary(table, None)),
        }
    }
}

#[derive(Debug, Clone)]
struct PackedCountMinSketch {
    cells: usize,
    hashes: usize,
    bits: u8,
    max_count: u64,
    layout: KCountArrayLayout,
    update_mode: CountMinUpdateMode,
    words: Vec<u64>,
    increments: u64,
    occupied_slots: usize,
    tracked_slots: Option<Vec<usize>>,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum CountMinUpdateMode {
    Conservative,
    Independent,
}

impl CountMinUpdateMode {
    fn as_str(self) -> &'static str {
        match self {
            Self::Conservative => "conservative",
            Self::Independent => "independent",
        }
    }
}

struct AtomicCountMinSketch {
    cells: usize,
    hashes: usize,
    max_count: u32,
    layout: KCountArrayLayout,
    update_mode: CountMinUpdateMode,
    parallel_replay: bool,
    cells_by_hash: Vec<AtomicU32>,
    locks: Vec<Mutex<()>>,
    increments: AtomicU64,
    occupied_slots: AtomicUsize,
}

enum PrefilterCountMinSketch {
    Packed(PackedCountMinSketch),
    AtomicPacked(AtomicPackedCountMinSketch),
}

struct AtomicPackedCountMinSketch {
    cells: usize,
    hashes: usize,
    bits: u8,
    max_count: u64,
    layout: KCountArrayLayout,
    update_mode: CountMinUpdateMode,
    words: Vec<AtomicU64>,
    locks: Vec<Mutex<()>>,
    increments: AtomicU64,
    occupied_slots: AtomicUsize,
}

const BBTOOLS_HASH_BITS: u32 = 6;
const BBTOOLS_HASH_ARRAY_LENGTH: usize = 1 << BBTOOLS_HASH_BITS;
const BBTOOLS_HASH_CELL_MASK: u64 = (BBTOOLS_HASH_ARRAY_LENGTH as u64) - 1;
const BBTOOLS_LONG_MAX_VALUE: u64 = i64::MAX as u64;
type BbtoolsHashMaskTable = [[u64; BBTOOLS_HASH_ARRAY_LENGTH]; 8];
type BbtoolsHashMaskRef = &'static BbtoolsHashMaskTable;
type BbtoolsHashMaskCache = FxHashMap<u64, BbtoolsHashMaskRef>;

#[derive(Debug, Clone, Copy)]
struct KCountArrayLayout {
    array_mask: u64,
    array_bits: u32,
    cells_per_array: usize,
    mask_seed: u64,
    masks: BbtoolsHashMaskRef,
}

const COUNT_PARALLEL_CHUNK_SIZE: usize = 8192;
const COUNT_CHUNK_LOCAL_MAP_MAX_CAPACITY: usize = 131_072;
const COUNTUP_SORT_RUN_PAIR_LIMIT: usize = 65_536;
const COUNTUP_SORT_RUN_BYTE_LIMIT: usize = 64 * 1024 * 1024;
const COUNTUP_SORT_MERGE_FANIN: usize = 128;
const COUNTUP_RUN_IO_BUFFER_CAPACITY: usize = 1024 * 1024;
const COUNTUP_PREPASS_CHUNK_PAIR_LIMIT: usize = 1024;
const COUNTUP_PREPASS_CHUNK_BYTE_LIMIT: usize = 16 * 1024 * 1024;
const HIST_PARALLEL_CHUNK_SIZE: usize = 1024;
const NORMALIZE_PARALLEL_CHUNK_SIZE: usize = 1024;
const PAIRED_ANALYSIS_JOIN_MIN_BASES: usize = 1024;
const COVERAGE_PAR_SORT_MIN_WINDOWS: usize = 4096;
const OVERLAP_AUTO_SAMPLE_PAIRS: u64 = 1_000_000;
const ATOMIC_SKETCH_PAR_REPLAY_MIN_KEYS: usize = 16_384;
const PACKED_SKETCH_TRACKED_SLOT_LIMIT: usize = 8_000_000;
const OVERLAP_AUTO_SAMPLE_INTERVAL: u64 = 100;
const OVERLAP_AUTO_ENABLE_FRACTION: f64 = 0.25;
const DEFAULT_PREFILTER_CELLS: usize = 1 << 20;
const DEFAULT_PREFILTER_BITS: u8 = 2;
const DEFAULT_PREFILTER_FRACTION_MICROS: u32 = 350_000;
const OUTPUT_COUNT_MIN_AUTO_FRACTION_MICROS: u32 = 250_000;
const OUTPUT_COUNT_MIN_AUTO_MIN_MEMORY_BYTES: usize = 64 * 1024 * 1024;
const AUTO_COUNT_MIN_FALLBACK_MEMORY_BYTES: usize = 2 * 1024 * 1024 * 1024;
const AUTO_COUNT_MIN_MAX_MEMORY_BYTES: usize = 2 * 1024 * 1024 * 1024;
const AUTO_COUNT_MIN_MIN_MEMORY_BYTES: usize = 256 * 1024 * 1024;
const BBTOOLS_MEMORY_HEADROOM_BYTES: usize = 96_000_000;
const EXPLICIT_COUNT_MIN_SAFE_MEMORY_PERCENT: usize = 85;
const BBTOOLS_KCOUNT_ARRAY_MIN_ARRAYS: usize = 2;
const BBTOOLS_KCOUNT_ARRAY_SHARD_MIN_CELLS: usize = 64;
const BBTOOLS_KCOUNT_ARRAY_MAX_HASHES: usize = 8;
const BBTOOLS_KCOUNT_ARRAY_LOCKS: usize = 1999;
const BBTOOLS_KCOUNT_ARRAY_FIRST_MASK_SEED: u64 = 0;
const BBTOOLS_KCOUNT_ARRAY_SECOND_MASK_SEED: u64 = 7;
const BBTOOLS_KCOUNT_ARRAY_MASK_SEED_STEP: u64 = 7;
const BBTOOLS_KCOUNT_ARRAY_THIRD_MASK_SEED: u64 =
    BBTOOLS_KCOUNT_ARRAY_SECOND_MASK_SEED + BBTOOLS_KCOUNT_ARRAY_MASK_SEED_STEP;
const PEAK_COMPACT_ZERO_TAIL: usize = 32;
static NONDETERMINISTIC_SEED_COUNTER: AtomicU64 = AtomicU64::new(0);
type AnalysisPair = (SequenceRecord, Option<SequenceRecord>, Option<f64>);
type NormalizationInput = (usize, SequenceRecord, Option<SequenceRecord>, f64);
type SparseHist = FxHashMap<usize, u64>;
type SparseReadDepthHist = FxHashMap<usize, (u64, u64)>;

#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub struct UniqueKmerEstimateSplit {
    pub low_depth_max: u64,
    pub low_depth_kmers: usize,
    pub high_depth_min: u64,
    pub high_depth_kmers: usize,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct SketchLayoutSummary {
    pub table: &'static str,
    pub kind: &'static str,
    pub cells: usize,
    pub hashes: usize,
    pub bits: u8,
    pub arrays: usize,
    pub cells_per_array: usize,
    pub mask_seed: u64,
    pub update_mode: &'static str,
    pub max_count: u64,
    pub memory_bytes: usize,
    pub prefilter_limit: Option<u64>,
}

#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub struct StageTiming {
    pub name: &'static str,
    pub elapsed_micros: u128,
}

#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub struct CountupSpillSummary {
    pub initial_runs: usize,
    pub merge_runs: usize,
    pub final_runs: usize,
    pub bytes_written: u64,
    pub peak_live_bytes: u64,
    pub final_live_bytes: u64,
}

impl CountupSpillSummary {
    pub fn has_spills(&self) -> bool {
        self.initial_runs > 0 || self.merge_runs > 0 || self.bytes_written > 0
    }

    fn note_initial_run(&mut self, bytes: u64) {
        self.initial_runs = self.initial_runs.saturating_add(1);
        self.note_written(bytes);
    }

    fn note_merge_run(&mut self, bytes: u64) {
        self.merge_runs = self.merge_runs.saturating_add(1);
        self.note_written(bytes);
    }

    fn note_written(&mut self, bytes: u64) {
        self.bytes_written = self.bytes_written.saturating_add(bytes);
        self.final_live_bytes = self.final_live_bytes.saturating_add(bytes);
        self.peak_live_bytes = self.peak_live_bytes.max(self.final_live_bytes);
    }

    fn note_removed(&mut self, bytes: u64) {
        self.final_live_bytes = self.final_live_bytes.saturating_sub(bytes);
    }
}

#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub struct RunSummary {
    pub reads_in: u64,
    pub bases_in: u64,
    pub reads_kept: u64,
    pub reads_tossed: u64,
    pub bases_kept: u64,
    pub bases_tossed: u64,
    pub unique_kmers_in: usize,
    pub unique_kmers_in_split: Option<UniqueKmerEstimateSplit>,
    pub unique_kmers_out: Option<usize>,
    pub cardinality_in: Option<CardinalityEstimate>,
    pub cardinality_out: Option<CardinalityEstimate>,
    pub sketch_layouts: Vec<SketchLayoutSummary>,
    pub stage_timings: Vec<StageTiming>,
    pub countup_spill: CountupSpillSummary,
}

#[derive(Debug, Clone, Default)]
struct ReadAnalysis {
    depth_al: Option<u64>,
    true_depth: Option<u64>,
    min_true_depth: Option<u64>,
    low_kmer_count: usize,
    total_kmer_count: usize,
    error: bool,
    had_kmer_windows: bool,
    coverage_desc: Vec<i64>,
}

#[derive(Debug, Clone, Default)]
struct PairAnalysis {
    read1: ReadAnalysis,
    read2: Option<ReadAnalysis>,
    depth_proxy_al: Option<u64>,
    max_true_depth: Option<u64>,
    low_kmer_count: usize,
    total_kmer_count: usize,
    error1: bool,
    error2: bool,
}

#[derive(Debug, Clone, Default)]
struct PairDecision {
    toss: bool,
    analysis: PairAnalysis,
}

#[derive(Debug, Clone, Default, PartialEq, Eq)]
struct CountupDecisionPlan {
    toss: bool,
    eligible_key_indices: Vec<usize>,
}

#[derive(Debug, Clone)]
struct NormalizedPair {
    input_list_index: usize,
    r1: SequenceRecord,
    r2: Option<SequenceRecord>,
    out_r1: SequenceRecord,
    out_r2: Option<SequenceRecord>,
    decision: PairDecision,
    uncorrectable: bool,
    read_count: u64,
    base_count: u64,
}

#[derive(Debug, Clone)]
struct CountupWorkPair {
    input_list_index: usize,
    sort_key: CountupSortKey,
    r1: SequenceRecord,
    r2: Option<SequenceRecord>,
}

#[derive(Debug, Clone)]
struct CountupWorkCandidate {
    input_list_index: usize,
    original_index: usize,
    rand: f64,
    r1: SequenceRecord,
    r2: Option<SequenceRecord>,
}

struct CountupWorkBuild {
    source: CountupWorkSource,
    input_hist: Option<SparseHist>,
    input_read_hist: Option<SparseReadDepthHist>,
    input_hist_elapsed_micros: u128,
    format1: SeqFormat,
    format2: Option<SeqFormat>,
    spill_summary: CountupSpillSummary,
}

struct CountupChunkBuild {
    work_pairs: Vec<CountupWorkPair>,
    depth_hist: SparseHist,
    read_hist: SparseReadDepthHist,
}

struct CountupInputHistAccumulator<'a> {
    wants_depth_hist: bool,
    wants_read_hist: bool,
    depth_hist: &'a mut SparseHist,
    read_hist: &'a mut SparseReadDepthHist,
}

#[derive(Debug, Clone)]
struct CountupSortKey {
    errors: usize,
    total_len: usize,
    expected_errors: f64,
    numeric_id: u64,
    original_index: usize,
}

struct CountupPrepassResult {
    include: bool,
    sort_analysis: Option<PairAnalysis>,
}

struct CountupWorkSource {
    temp_dir: Option<tempfile::TempDir>,
    inner: CountupWorkSourceInner,
}

enum CountupWorkSourceInner {
    Memory(Vec<CountupWorkPair>),
    Spilled(Vec<PathBuf>),
}

struct CountupWorkIter {
    _temp_dir: Option<tempfile::TempDir>,
    inner: CountupWorkIterInner,
}

enum CountupWorkIterInner {
    Memory(std::vec::IntoIter<CountupWorkPair>),
    Spilled(CountupRunMerger),
}

struct CountupRunMerger {
    readers: Vec<CountupRunReader>,
    heap: BinaryHeap<CountupRunHead>,
}

struct CountupRunReader {
    reader: BufReader<fs::File>,
}

struct CountupRunHead {
    pair: CountupWorkPair,
    run_index: usize,
}

#[derive(Debug, Clone, Copy, Default)]
struct CorrectionResult {
    corrected: usize,
    marked: usize,
    uncorrectable: bool,
}

#[derive(Debug, Clone, Copy)]
struct CorrectionTarget {
    low: i64,
    lower_bound: i64,
    upper_bound: i64,
    mult: i64,
}

#[derive(Debug, Clone)]
struct InputLists {
    first: Vec<PathBuf>,
    second: Option<Vec<PathBuf>>,
}

#[derive(Debug, Clone, Default)]
struct ReadDepthHistogram {
    reads: Vec<u64>,
    bases: Vec<u64>,
}

impl ReadDepthHistogram {
    fn new(len: usize) -> Self {
        Self {
            reads: vec![0; len],
            bases: vec![0; len],
        }
    }
}

#[derive(Debug, Clone, Copy, Default)]
struct BaseCounts {
    a: u64,
    c: u64,
    g: u64,
    t: u64,
    n: u64,
}

impl BaseCounts {
    fn total(self) -> u64 {
        self.a + self.c + self.g + self.t + self.n
    }
}

#[derive(Debug, Clone, Default)]
struct BaseContentHistogram {
    first: Vec<BaseCounts>,
    second: Vec<BaseCounts>,
}

#[derive(Debug, Clone, Copy, Default)]
struct MatchCounts {
    matches: u64,
    n: u64,
}

#[derive(Debug, Clone, Default)]
struct AlignmentFallbackHistograms {
    first_match: Vec<MatchCounts>,
    second_match: Vec<MatchCounts>,
    quality_match: Vec<u64>,
    read_count: u64,
    base_count: u64,
    pair_count: u64,
    paired: bool,
}

#[derive(Debug, Clone, Default)]
struct QualitySideHistograms {
    overall: Vec<u64>,
    first_counts: Vec<u64>,
    second_counts: Vec<u64>,
    first_avg: Vec<u64>,
    second_avg: Vec<u64>,
    first_by_pos: Vec<Vec<u64>>,
    second_by_pos: Vec<Vec<u64>>,
    paired: bool,
}

#[derive(Debug, Clone, Default)]
struct ReadLocalSideHistograms {
    quality: Option<QualitySideHistograms>,
    length: Option<ReadDepthHistogram>,
    gc: Option<ReadDepthHistogram>,
    base: Option<BaseContentHistogram>,
    entropy: Option<Vec<u64>>,
    identity: Option<ReadDepthHistogram>,
    alignment: Option<AlignmentFallbackHistograms>,
    barcodes: Option<BTreeMap<String, u64>>,
}

#[derive(Debug, Clone)]
struct JavaXoshiro {
    s0: u64,
    s1: u64,
    s2: u64,
    s3: u64,
}

impl JavaXoshiro {
    fn new(seed: u64) -> Self {
        let mut rng = Self {
            s0: seed,
            s1: mix_seed(seed),
            s2: 0,
            s3: 0,
        };
        rng.s2 = mix_seed(rng.s1);
        rng.s3 = mix_seed(rng.s2);
        if rng.s0 == 0 && rng.s1 == 0 && rng.s2 == 0 && rng.s3 == 0 {
            rng.s0 = 0x5DEECE66D;
            rng.s1 = 0xB;
            rng.s2 = 0xCCA;
            rng.s3 = 0xF00;
        }
        for _ in 0..4 {
            rng.next_long();
        }
        rng
    }

    fn next_long(&mut self) -> u64 {
        let result = self.s0.wrapping_add(self.s3);
        let t = self.s1 << 17;

        self.s2 ^= self.s0;
        self.s3 ^= self.s1;
        self.s1 ^= self.s2;
        self.s0 ^= self.s3;

        self.s2 ^= t;
        self.s3 = self.s3.rotate_left(45);

        result
    }

    fn next_double(&mut self) -> f64 {
        ((self.next_long() >> 11) as f64) * (1.0 / ((1u64 << 53) as f64))
    }
}

fn run_random_seed(config: &Config) -> u64 {
    if config.deterministic {
        0
    } else {
        nondeterministic_seed()
    }
}

fn nondeterministic_seed() -> u64 {
    let nanos = SystemTime::now()
        .duration_since(UNIX_EPOCH)
        .map(|duration| duration.as_nanos() as u64)
        .unwrap_or(0);
    let counter = NONDETERMINISTIC_SEED_COUNTER.fetch_add(1, Ordering::Relaxed);
    nanos ^ ((std::process::id() as u64) << 32) ^ mix_seed(counter.wrapping_add(0x9E37_79B9))
}

fn mix_seed(mut x: u64) -> u64 {
    x = x.wrapping_add(0x9E37_79B9_7F4A_7C15);
    x = (x ^ (x >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
    x = (x ^ (x >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
    x ^ (x >> 31)
}

struct PrimaryReaders {
    r1: SequenceReader,
    r2: Option<SequenceReader>,
    interleaved: bool,
    input_list1: Vec<PathBuf>,
    input_list2: Option<Vec<PathBuf>>,
    input_list_index: usize,
    settings: SequenceSettings,
    limit_per_file: Option<u64>,
    pairs_seen_in_file: u64,
    format1: SeqFormat,
    format2: Option<SeqFormat>,
    next_pair_numeric_id: u64,
    gzip_threads: Option<usize>,
}

impl PrimaryReaders {
    fn open(config: &Config, limit_per_file: Option<u64>) -> Result<Self> {
        let in1 = config.in1.as_ref().context("missing in1")?;
        let sequence_settings = sequence_settings(config);
        let input_list = primary_input_lists(config);
        let first_path = input_list
            .as_ref()
            .and_then(|paths| paths.first.first())
            .unwrap_or(in1);
        let r2_path = input_list
            .as_ref()
            .and_then(|paths| paths.second.as_ref())
            .and_then(|paths| paths.first())
            .or(config.in2.as_ref());
        let gzip_threads = gzip_threads_for_paths(
            config.gzip_threads,
            [Some(first_path.as_path()), r2_path.map(PathBuf::as_path)],
        );
        let r1 =
            open_sequence_reader_with_gzip_threads(first_path, sequence_settings, gzip_threads)?;
        let interleaved = input_list.is_none()
            && config.in2.is_none()
            && (config.interleaved
                || (config.test_interleaved
                    && detect_interleaved_input_with_gzip_threads(
                        first_path,
                        sequence_settings,
                        config.gzip_threads,
                    )?));
        let r2 = r2_path
            .map(|path| {
                open_sequence_reader_with_gzip_threads(path, sequence_settings, gzip_threads)
            })
            .transpose()?;
        if let Some(r2_ref) = &r2
            && r1.format() != r2_ref.format()
        {
            bail!("paired inputs must use the same FASTA/FASTQ format");
        }
        let format1 = r1.format();
        let format2 = if interleaved {
            Some(format1)
        } else {
            r2.as_ref().map(SequenceReader::format)
        };
        Ok(Self {
            r1,
            r2,
            interleaved,
            input_list1: input_list
                .as_ref()
                .map(|paths| paths.first.clone())
                .unwrap_or_default(),
            input_list2: input_list.and_then(|paths| paths.second),
            input_list_index: 0,
            settings: sequence_settings,
            limit_per_file,
            pairs_seen_in_file: 0,
            format1,
            format2,
            next_pair_numeric_id: 0,
            gzip_threads: config.gzip_threads,
        })
    }

    fn format1(&self) -> SeqFormat {
        self.format1
    }

    fn format2(&self) -> Option<SeqFormat> {
        self.format2
    }

    fn input_list_index(&self) -> usize {
        self.input_list_index
    }

    fn next_pair(&mut self) -> Result<Option<(SequenceRecord, Option<SequenceRecord>)>> {
        if !self.input_list1.is_empty() {
            return self.next_list_record();
        }
        if limit_reached(self.limit_per_file, self.pairs_seen_in_file) {
            return Ok(None);
        }
        let r1 = self.r1.next_record()?;
        if self.interleaved {
            return match r1 {
                Some(mut record) => {
                    let mut mate = self
                        .r1
                        .next_record()?
                        .context("interleaved input ended after an unmatched first mate record")?;
                    record.numeric_id = self.next_pair_numeric_id;
                    mate.numeric_id = self.next_pair_numeric_id;
                    self.next_pair_numeric_id += 1;
                    self.pairs_seen_in_file += 1;
                    Ok(Some((record, Some(mate))))
                }
                None => Ok(None),
            };
        }

        let r2 = match &mut self.r2 {
            Some(reader) => reader.next_record()?,
            None => None,
        };

        match (r1, r2) {
            (None, None) => Ok(None),
            (Some(record), mate) => {
                self.pairs_seen_in_file += 1;
                Ok(Some((record, mate)))
            }
            (None, Some(_)) => bail!("in2 has more records than in1"),
        }
    }

    fn next_list_record(&mut self) -> Result<Option<(SequenceRecord, Option<SequenceRecord>)>> {
        loop {
            if limit_reached(self.limit_per_file, self.pairs_seen_in_file) {
                if !self.advance_list_reader()? {
                    return Ok(None);
                }
                continue;
            }
            let had_r2 = self.r2.is_some();
            let r1 = self.r1.next_record()?;
            let r2 = match &mut self.r2 {
                Some(reader) => reader.next_record()?,
                None => None,
            };
            match (r1, r2) {
                (Some(record), Some(mate)) => {
                    self.pairs_seen_in_file += 1;
                    return Ok(Some((record, Some(mate))));
                }
                (Some(record), None) if !had_r2 => {
                    self.pairs_seen_in_file += 1;
                    return Ok(Some((record, None)));
                }
                (Some(_), None) => bail!("in2 has fewer records than in1"),
                (None, Some(_)) => bail!("in2 has more records than in1"),
                (None, None) => {
                    if !self.advance_list_reader()? {
                        return Ok(None);
                    }
                }
            }
        }
    }

    fn advance_list_reader(&mut self) -> Result<bool> {
        if self.input_list_index + 1 >= self.input_list1.len() {
            return Ok(false);
        }
        self.input_list_index += 1;
        let path = &self.input_list1[self.input_list_index];
        let second_path = self
            .input_list2
            .as_ref()
            .and_then(|paths| paths.get(self.input_list_index));
        let gzip_threads = gzip_threads_for_paths(
            self.gzip_threads,
            [Some(path.as_path()), second_path.map(PathBuf::as_path)],
        );
        let reader =
            SequenceReader::from_path_with_gzip_threads(path, self.settings, gzip_threads)?;
        if reader.format() != self.format1 {
            bail!("comma-separated input list entries must use the same FASTA/FASTQ format");
        }
        self.r2 = self
            .input_list2
            .as_ref()
            .and_then(|paths| paths.get(self.input_list_index))
            .map(|path| {
                SequenceReader::from_path_with_gzip_threads(path, self.settings, gzip_threads)
            })
            .transpose()?;
        if let Some(r2_ref) = &self.r2
            && Some(r2_ref.format()) != self.format2
        {
            bail!("comma-separated paired input list entries must use the same FASTA/FASTQ format");
        }
        self.r1 = reader;
        self.pairs_seen_in_file = 0;
        Ok(true)
    }
}

fn open_sequence_reader(
    config: &Config,
    path: &Path,
    settings: SequenceSettings,
) -> Result<SequenceReader> {
    SequenceReader::from_path_with_gzip_threads(path, settings, config.gzip_threads)
}

fn open_sequence_reader_with_gzip_threads(
    path: &Path,
    settings: SequenceSettings,
    gzip_threads: Option<usize>,
) -> Result<SequenceReader> {
    SequenceReader::from_path_with_gzip_threads(path, settings, gzip_threads)
}

fn open_paired_sequence_readers(
    config: &Config,
    path1: &Path,
    path2: &Path,
    settings: SequenceSettings,
) -> Result<(SequenceReader, SequenceReader)> {
    let gzip_threads = gzip_threads_for_paths(config.gzip_threads, [Some(path1), Some(path2)]);
    let reader1 = open_sequence_reader_with_gzip_threads(path1, settings, gzip_threads)?;
    let reader2 = open_sequence_reader_with_gzip_threads(path2, settings, gzip_threads)?;
    Ok((reader1, reader2))
}

fn gzip_threads_for_paths<'a>(
    gzip_threads: Option<usize>,
    paths: impl IntoIterator<Item = Option<&'a Path>>,
) -> Option<usize> {
    let gzip_streams = paths
        .into_iter()
        .flatten()
        .filter(|path| path_uses_gzip(path))
        .count();
    gzip_threads_for_streams(gzip_threads, gzip_streams)
}

fn gzip_threads_for_streams(gzip_threads: Option<usize>, gzip_streams: usize) -> Option<usize> {
    gzip_threads.map(|threads| {
        if threads <= 1 || gzip_streams <= 1 {
            threads
        } else {
            (threads / gzip_streams).max(1)
        }
    })
}

fn path_uses_gzip(path: &Path) -> bool {
    path.extension()
        .and_then(|ext| ext.to_str())
        .is_some_and(|ext| ext.eq_ignore_ascii_case("gz"))
}

struct OptionalWriters {
    interleaved_output: bool,
    current_output_list_index: usize,
    keep_plan: OutputPathPlan,
    toss_plan: OutputPathPlan,
    low_plan: OutputPathPlan,
    mid_plan: OutputPathPlan,
    high_plan: OutputPathPlan,
    uncorrected_plan: OutputPathPlan,
    keep1: Option<SequenceWriter>,
    keep2: Option<SequenceWriter>,
    toss1: Option<SequenceWriter>,
    toss2: Option<SequenceWriter>,
    low1: Option<SequenceWriter>,
    low2: Option<SequenceWriter>,
    mid1: Option<SequenceWriter>,
    mid2: Option<SequenceWriter>,
    high1: Option<SequenceWriter>,
    high2: Option<SequenceWriter>,
    uncorrected1: Option<SequenceWriter>,
    uncorrected2: Option<SequenceWriter>,
}

impl OptionalWriters {
    fn open(config: &Config, _format1: SeqFormat, format2: Option<SeqFormat>) -> Result<Self> {
        if format2.is_none() && has_second_output(config) {
            bail!(
                "second-output paths require paired input; interleaved auto-detection did not detect paired records"
            );
        }
        let paired = format2.is_some();
        let input_list_len = primary_input_lists(config)
            .map(|paths| paths.first.len())
            .unwrap_or(1);
        let keep_plan = prepare_output_path_plan(
            config.out1.as_deref(),
            config.out2.as_deref(),
            paired,
            input_list_len,
        )?;
        let toss_plan = prepare_output_path_plan(
            config.out_toss1.as_deref(),
            config.out_toss2.as_deref(),
            paired,
            input_list_len,
        )?;
        let low_plan = prepare_output_path_plan(
            config.out_low1.as_deref(),
            config.out_low2.as_deref(),
            paired,
            input_list_len,
        )?;
        let mid_plan = prepare_output_path_plan(
            config.out_mid1.as_deref(),
            config.out_mid2.as_deref(),
            paired,
            input_list_len,
        )?;
        let high_plan = prepare_output_path_plan(
            config.out_high1.as_deref(),
            config.out_high2.as_deref(),
            paired,
            input_list_len,
        )?;
        let uncorrected_plan = prepare_output_path_plan(
            config.out_uncorrected1.as_deref(),
            config.out_uncorrected2.as_deref(),
            paired,
            input_list_len,
        )?;
        let output_gzip_threads = output_gzip_threads_for_plans(
            config.gzip_threads,
            [
                &keep_plan,
                &toss_plan,
                &low_plan,
                &mid_plan,
                &high_plan,
                &uncorrected_plan,
            ],
            0,
        )?;
        let (keep1, keep2) = open_output_pair(
            keep_plan.pair_for_index(0)?,
            config.overwrite,
            config.append,
            config.quality_out_offset,
            config.fake_quality,
            config.fasta_wrap,
            output_gzip_threads,
        )?;
        let (toss1, toss2) = open_output_pair(
            toss_plan.pair_for_index(0)?,
            config.overwrite,
            config.append,
            config.quality_out_offset,
            config.fake_quality,
            config.fasta_wrap,
            output_gzip_threads,
        )?;
        let (low1, low2) = open_output_pair(
            low_plan.pair_for_index(0)?,
            config.overwrite,
            config.append,
            config.quality_out_offset,
            config.fake_quality,
            config.fasta_wrap,
            output_gzip_threads,
        )?;
        let (mid1, mid2) = open_output_pair(
            mid_plan.pair_for_index(0)?,
            config.overwrite,
            config.append,
            config.quality_out_offset,
            config.fake_quality,
            config.fasta_wrap,
            output_gzip_threads,
        )?;
        let (high1, high2) = open_output_pair(
            high_plan.pair_for_index(0)?,
            config.overwrite,
            config.append,
            config.quality_out_offset,
            config.fake_quality,
            config.fasta_wrap,
            output_gzip_threads,
        )?;
        let (uncorrected1, uncorrected2) = open_output_pair(
            uncorrected_plan.pair_for_index(0)?,
            config.overwrite,
            config.append,
            config.quality_out_offset,
            config.fake_quality,
            config.fasta_wrap,
            output_gzip_threads,
        )?;
        Ok(Self {
            interleaved_output: paired,
            current_output_list_index: 0,
            keep_plan,
            toss_plan,
            low_plan,
            mid_plan,
            high_plan,
            uncorrected_plan,
            keep1,
            keep2,
            toss1,
            toss2,
            low1,
            low2,
            mid1,
            mid2,
            high1,
            high2,
            uncorrected1,
            uncorrected2,
        })
    }

    fn sync_to_input_list_index(&mut self, config: &Config, index: usize) -> Result<()> {
        if self.current_output_list_index == index {
            return Ok(());
        }
        self.flush()?;
        let output_gzip_threads = output_gzip_threads_for_plans(
            config.gzip_threads,
            [
                &self.keep_plan,
                &self.toss_plan,
                &self.low_plan,
                &self.mid_plan,
                &self.high_plan,
                &self.uncorrected_plan,
            ],
            index,
        )?;
        reopen_output_pair_if_fanout(
            &self.keep_plan,
            index,
            &mut self.keep1,
            &mut self.keep2,
            config,
            output_gzip_threads,
        )?;
        reopen_output_pair_if_fanout(
            &self.toss_plan,
            index,
            &mut self.toss1,
            &mut self.toss2,
            config,
            output_gzip_threads,
        )?;
        reopen_output_pair_if_fanout(
            &self.low_plan,
            index,
            &mut self.low1,
            &mut self.low2,
            config,
            output_gzip_threads,
        )?;
        reopen_output_pair_if_fanout(
            &self.mid_plan,
            index,
            &mut self.mid1,
            &mut self.mid2,
            config,
            output_gzip_threads,
        )?;
        reopen_output_pair_if_fanout(
            &self.high_plan,
            index,
            &mut self.high1,
            &mut self.high2,
            config,
            output_gzip_threads,
        )?;
        reopen_output_pair_if_fanout(
            &self.uncorrected_plan,
            index,
            &mut self.uncorrected1,
            &mut self.uncorrected2,
            config,
            output_gzip_threads,
        )?;
        self.current_output_list_index = index;
        Ok(())
    }

    fn write_pair(
        &mut self,
        toss: bool,
        r1: &SequenceRecord,
        r2: Option<&SequenceRecord>,
    ) -> Result<()> {
        if toss {
            write_to_optional_pair(
                &mut self.toss1,
                &mut self.toss2,
                self.interleaved_output,
                r1,
                r2,
            )?;
        } else {
            write_to_optional_pair(
                &mut self.keep1,
                &mut self.keep2,
                self.interleaved_output,
                r1,
                r2,
            )?;
        }
        Ok(())
    }

    fn write_depth_bin(
        &mut self,
        config: &Config,
        analysis: &PairAnalysis,
        r1: &SequenceRecord,
        r2: Option<&SequenceRecord>,
    ) -> Result<()> {
        let d1 = bin_depth(analysis.read1.depth_al);
        let d2 = analysis
            .read2
            .as_ref()
            .map(|read| bin_depth(read.depth_al))
            .unwrap_or(-1);
        let target = if d1 < config.low_bin_depth && d2 < config.low_bin_depth {
            DepthBin::Low
        } else if (d1 < config.low_bin_depth || d1 > config.high_bin_depth)
            && (d2 < config.low_bin_depth || d2 >= config.high_bin_depth)
        {
            DepthBin::High
        } else {
            DepthBin::Mid
        };

        match target {
            DepthBin::Low => write_to_optional_pair(
                &mut self.low1,
                &mut self.low2,
                self.interleaved_output,
                r1,
                r2,
            )?,
            DepthBin::Mid => write_to_optional_pair(
                &mut self.mid1,
                &mut self.mid2,
                self.interleaved_output,
                r1,
                r2,
            )?,
            DepthBin::High => write_to_optional_pair(
                &mut self.high1,
                &mut self.high2,
                self.interleaved_output,
                r1,
                r2,
            )?,
        }
        Ok(())
    }

    fn write_uncorrected(
        &mut self,
        r1: &SequenceRecord,
        r2: Option<&SequenceRecord>,
    ) -> Result<()> {
        write_to_optional_pair(
            &mut self.uncorrected1,
            &mut self.uncorrected2,
            self.interleaved_output,
            r1,
            r2,
        )
    }

    fn flush(&mut self) -> Result<()> {
        for writer in [
            self.keep1.as_mut(),
            self.keep2.as_mut(),
            self.toss1.as_mut(),
            self.toss2.as_mut(),
            self.low1.as_mut(),
            self.low2.as_mut(),
            self.mid1.as_mut(),
            self.mid2.as_mut(),
            self.high1.as_mut(),
            self.high2.as_mut(),
            self.uncorrected1.as_mut(),
            self.uncorrected2.as_mut(),
        ]
        .into_iter()
        .flatten()
        {
            writer.flush()?;
        }
        Ok(())
    }
}

#[derive(Debug, Clone, Copy)]
enum DepthBin {
    Low,
    Mid,
    High,
}

fn write_to_optional_pair(
    writer1: &mut Option<SequenceWriter>,
    writer2: &mut Option<SequenceWriter>,
    interleaved_output: bool,
    r1: &SequenceRecord,
    r2: Option<&SequenceRecord>,
) -> Result<()> {
    if let Some(writer) = writer1.as_mut() {
        writer.write_record(r1)?;
        if interleaved_output && writer2.is_none() {
            if let Some(mate) = r2 {
                writer.write_record(mate)?;
            }
            return Ok(());
        }
    }
    if let (Some(writer), Some(mate)) = (writer2.as_mut(), r2) {
        writer.write_record(mate)?;
    }
    Ok(())
}

fn has_second_output(config: &Config) -> bool {
    config.out2.is_some()
        || config.out_toss2.is_some()
        || config.out_low2.is_some()
        || config.out_mid2.is_some()
        || config.out_high2.is_some()
        || config.out_uncorrected2.is_some()
}

fn depth_bin_outputs_enabled(config: &Config) -> bool {
    config.out_low1.is_some()
        || config.out_low2.is_some()
        || config.out_mid1.is_some()
        || config.out_mid2.is_some()
        || config.out_high1.is_some()
        || config.out_high2.is_some()
}

fn needs_output_pair_analysis(config: &Config) -> bool {
    config.rename_reads || depth_bin_outputs_enabled(config)
}

#[derive(Debug, Clone)]
struct OutputPathPair {
    first: Option<PathBuf>,
    second: Option<PathBuf>,
}

#[derive(Debug, Clone)]
struct OutputPathPlan {
    pairs: Vec<OutputPathPair>,
    fanout: bool,
}

impl OutputPathPlan {
    fn pair_for_index(&self, index: usize) -> Result<&OutputPathPair> {
        if self.fanout {
            self.pairs
                .get(index)
                .with_context(|| format!("missing output path list entry for input {}", index + 1))
        } else {
            self.pairs.first().context("missing output path plan entry")
        }
    }
}

fn prepare_output_paths(
    first: Option<&Path>,
    second: Option<&Path>,
    paired: bool,
) -> OutputPathPair {
    let second = match second {
        Some(path) => Some(path.to_path_buf()),
        None if paired => first.and_then(|path| replace_hash_in_path(path, "2")),
        None => None,
    };
    let first =
        first.map(|path| replace_hash_in_path(path, "1").unwrap_or_else(|| path.to_path_buf()));
    OutputPathPair { first, second }
}

fn prepare_output_path_plan(
    first: Option<&Path>,
    second: Option<&Path>,
    paired: bool,
    input_list_len: usize,
) -> Result<OutputPathPlan> {
    if input_list_len > 1
        && let Some(first_values) = output_path_values(first)
        && first_values.len() > 1
    {
        let second_values = output_path_values(second);
        let fanout_len = second_values
            .as_ref()
            .map(|values| first_values.len().min(values.len()))
            .unwrap_or(first_values.len());
        let mut pairs = Vec::with_capacity(fanout_len);
        for index in 0..fanout_len {
            let mut first_path = first_values[index].clone();
            let second_path = if let Some(values) = &second_values {
                Some(values[index].clone())
            } else if paired {
                if let Some(second_path) = replace_hash_in_path(&first_path, "2") {
                    first_path = replace_hash_in_path(&first_path, "1").unwrap_or(first_path);
                    Some(second_path)
                } else {
                    None
                }
            } else {
                None
            };
            pairs.push(OutputPathPair {
                first: Some(first_path),
                second: second_path,
            });
        }
        return Ok(OutputPathPlan {
            pairs,
            fanout: true,
        });
    }

    if input_list_len > 1
        && let Some(second_values) = output_path_values(second)
        && second_values.len() > 1
    {
        let first_path =
            first.map(|path| replace_hash_in_path(path, "1").unwrap_or_else(|| path.to_path_buf()));
        return Ok(OutputPathPlan {
            pairs: vec![OutputPathPair {
                first: first_path,
                second: Some(second_values[0].clone()),
            }],
            fanout: false,
        });
    }

    Ok(OutputPathPlan {
        pairs: vec![prepare_output_paths(first, second, paired)],
        fanout: false,
    })
}

fn output_path_values(path: Option<&Path>) -> Option<Vec<PathBuf>> {
    let path = path?;
    if path.exists() {
        return Some(vec![path.to_path_buf()]);
    }
    let text = path.to_string_lossy();
    if text.contains(',') {
        let paths = split_path_list(&text);
        if paths.len() > 1 {
            return Some(paths);
        }
    }
    Some(vec![path.to_path_buf()])
}

fn reopen_output_pair_if_fanout(
    plan: &OutputPathPlan,
    index: usize,
    first: &mut Option<SequenceWriter>,
    second: &mut Option<SequenceWriter>,
    config: &Config,
    gzip_threads: Option<usize>,
) -> Result<()> {
    if !plan.fanout {
        return Ok(());
    }
    *first = None;
    *second = None;
    let (new_first, new_second) = open_output_pair(
        plan.pair_for_index(index)?,
        config.overwrite,
        config.append,
        config.quality_out_offset,
        config.fake_quality,
        config.fasta_wrap,
        gzip_threads,
    )?;
    *first = new_first;
    *second = new_second;
    Ok(())
}

fn output_gzip_threads_for_plans<'a>(
    gzip_threads: Option<usize>,
    plans: impl IntoIterator<Item = &'a OutputPathPlan>,
    index: usize,
) -> Result<Option<usize>> {
    let mut gzip_streams = 0usize;
    for plan in plans {
        gzip_streams =
            gzip_streams.saturating_add(output_pair_gzip_streams(plan.pair_for_index(index)?));
    }
    Ok(gzip_threads_for_streams(gzip_threads, gzip_streams))
}

fn output_pair_gzip_streams(pair: &OutputPathPair) -> usize {
    [pair.first.as_deref(), pair.second.as_deref()]
        .into_iter()
        .flatten()
        .filter(|path| path_uses_gzip(path))
        .count()
}

fn open_output_pair(
    pair: &OutputPathPair,
    overwrite: bool,
    append: bool,
    quality_out_offset: u8,
    fake_quality: u8,
    fasta_wrap: usize,
    gzip_threads: Option<usize>,
) -> Result<(Option<SequenceWriter>, Option<SequenceWriter>)> {
    let first = open_sequence_writer(
        pair.first.as_deref(),
        overwrite,
        append,
        quality_out_offset,
        fake_quality,
        fasta_wrap,
        gzip_threads,
    )?;
    let second = open_sequence_writer(
        pair.second.as_deref(),
        overwrite,
        append,
        quality_out_offset,
        fake_quality,
        fasta_wrap,
        gzip_threads,
    )?;
    Ok((first, second))
}

fn replace_hash_in_path(path: &Path, replacement: &str) -> Option<PathBuf> {
    let text = path.to_string_lossy();
    if text.contains('#') {
        Some(PathBuf::from(text.replacen('#', replacement, 1)))
    } else {
        None
    }
}

fn bin_depth(depth: Option<u64>) -> i64 {
    depth
        .and_then(|value| i64::try_from(value).ok())
        .unwrap_or(-1)
}

pub fn run(config: &Config) -> Result<RunSummary> {
    let resolved_config;
    let config = if config.overlap_error_correct_auto {
        resolved_config = resolve_overlap_error_correct_auto(config)?;
        &resolved_config
    } else {
        config
    };
    if config.passes > 1 {
        return run_multipass(config);
    }
    run_single_pass(config)
}

fn resolve_overlap_error_correct_auto(config: &Config) -> Result<Config> {
    let mut resolved = config.clone();
    resolved.overlap_error_correct_auto = false;
    resolved.overlap_error_correct = sampled_overlap_fraction(config)?
        .is_some_and(|fraction| fraction > OVERLAP_AUTO_ENABLE_FRACTION);
    Ok(resolved)
}

fn sampled_overlap_fraction(config: &Config) -> Result<Option<f64>> {
    let mut readers = PrimaryReaders::open(config, Some(OVERLAP_AUTO_SAMPLE_PAIRS))?;
    let mut sampled = 0u64;
    let mut seen = 0u64;
    let mut mergeable = 0u64;
    while let Some((r1, r2)) = readers.next_pair()? {
        let Some(r2) = r2 else {
            return Ok(None);
        };
        seen += 1;
        if !seen.is_multiple_of(OVERLAP_AUTO_SAMPLE_INTERVAL) {
            continue;
        }
        sampled += 1;
        if best_pair_overlap(&r1, &r2).is_some() {
            mergeable += 1;
        }
    }
    if sampled == 0 {
        Ok(None)
    } else {
        Ok(Some(mergeable as f64 / sampled as f64))
    }
}

fn run_multipass(config: &Config) -> Result<RunSummary> {
    let mut multipass_config = config.clone();
    apply_bbtools_multipass_cell_bits_cap(&mut multipass_config);
    let config = &multipass_config;
    let temp_dir = managed_temp_dir(config, "bbnorm-rs-multipass-")?;
    let paired = config.in2.is_some() || config.interleaved;
    let separate_pair_outputs = paired && config.out2.is_some();
    let temp_ext = temp_sequence_extension(config);
    let mut last_in1 = config.in1.clone().context("missing in1")?;
    let mut last_in2 = config.in2.clone();
    let mut last_interleaved = config.interleaved;

    for pass in 1..config.passes {
        let temp1 = temp_dir.path().join(format!("pass{pass}.r1.{temp_ext}"));
        let temp2 = separate_pair_outputs
            .then(|| temp_dir.path().join(format!("pass{pass}.r2.{temp_ext}")));
        let mut pass_config = pass_config_for_intermediate(
            config,
            pass,
            &last_in1,
            last_in2.as_deref(),
            last_interleaved,
            temp1.clone(),
            temp2.clone(),
            None,
            None,
        );
        run_single_pass(&pass_config)
            .with_context(|| format!("running Rust multipass intermediate pass {pass}"))?;

        last_in1 = temp1;
        last_in2 = temp2;
        last_interleaved = paired && last_in2.is_none();
        pass_config.notes.clear();
    }

    let mut final_config = config.clone();
    final_config.in1 = Some(last_in1);
    final_config.in2 = last_in2;
    final_config.interleaved = last_interleaved;
    final_config.test_interleaved = !last_interleaved && final_config.in2.is_none();
    final_config.extra.clear();
    final_config.hist_in = None;
    final_config.rhist_in = None;
    final_config.peaks_in = None;
    final_config.match_hist_out = None;
    final_config.insert_hist_out = None;
    final_config.quality_accuracy_hist_out = None;
    final_config.indel_hist_out = None;
    final_config.error_hist_out = None;
    final_config.quality_hist_out = None;
    final_config.base_quality_hist_out = None;
    final_config.quality_count_hist_out = None;
    final_config.average_quality_hist_out = None;
    final_config.overall_base_quality_hist_out = None;
    final_config.length_hist_out = None;
    final_config.gc_hist_out = None;
    final_config.base_hist_out = None;
    final_config.entropy_hist_out = None;
    final_config.identity_hist_out = None;
    final_config.target_bad_percent_low = 1.0;
    final_config.target_bad_percent_high = 1.0;
    final_config.error_correct = config.error_correct_final;
    final_config.overlap_error_correct = config.overlap_error_correct && config.error_correct_final;
    final_config.passes = 1;
    let final_toss1 = config.out_toss1.as_ref().map(|_| {
        temp_dir
            .path()
            .join(format!("pass{}.final.toss1.{temp_ext}", config.passes))
    });
    let final_toss2 = config.out_toss2.as_ref().map(|_| {
        temp_dir
            .path()
            .join(format!("pass{}.final.toss2.{temp_ext}", config.passes))
    });
    final_config.out_toss1 = final_toss1.clone();
    final_config.out_toss2 = final_toss2.clone();

    let summary = run_single_pass(&final_config).context("running Rust multipass final pass")?;

    if let Some(path) = final_toss1
        && let Some(output) = config.out_toss1.as_deref()
    {
        write_multipass_fragments(
            &[path],
            output,
            config.overwrite,
            config.append,
            "multipass toss output",
        )?;
    }
    if let Some(path) = final_toss2
        && let Some(output) = config.out_toss2.as_deref()
    {
        write_multipass_fragments(
            &[path],
            output,
            config.overwrite,
            config.append,
            "multipass paired toss output",
        )?;
    }
    Ok(summary)
}

fn apply_bbtools_multipass_cell_bits_cap(config: &mut Config) {
    if config.passes > 1 && config.count_min.bits.unwrap_or(32) > 16 {
        config.count_min.bits = Some(16);
    }
}

fn managed_temp_dir(config: &Config, prefix: &str) -> Result<tempfile::TempDir> {
    let mut builder = tempfile::Builder::new();
    builder.prefix(prefix);
    if config.use_temp_dir
        && let Some(dir) = config.temp_dir.as_deref()
    {
        fs::create_dir_all(dir)
            .with_context(|| format!("creating temporary directory parent {}", dir.display()))?;
        return builder
            .tempdir_in(dir)
            .with_context(|| format!("creating managed temporary directory in {}", dir.display()));
    }
    builder
        .tempdir()
        .context("creating managed temporary directory")
}

fn write_multipass_fragments(
    fragments: &[PathBuf],
    output: &Path,
    overwrite: bool,
    append: bool,
    label: &str,
) -> Result<()> {
    let mut writer = create_output_with_append(output, overwrite, append)
        .with_context(|| format!("opening {label} {}", output.display()))?;
    for fragment in fragments {
        if fragment.exists() {
            let mut input = std::fs::File::open(fragment)
                .with_context(|| format!("opening multipass fragment {}", fragment.display()))?;
            std::io::copy(&mut input, &mut writer)
                .with_context(|| format!("copying multipass fragment {}", fragment.display()))?;
        }
    }
    writer
        .flush()
        .with_context(|| format!("flushing {label} {}", output.display()))?;
    Ok(())
}

#[allow(clippy::too_many_arguments)]
fn pass_config_for_intermediate(
    config: &Config,
    pass: usize,
    in1: &Path,
    in2: Option<&Path>,
    interleaved: bool,
    out1: PathBuf,
    out2: Option<PathBuf>,
    out_toss1: Option<PathBuf>,
    out_toss2: Option<PathBuf>,
) -> Config {
    let mut pass_config = config.clone();
    let target = intermediate_target_depth(config, pass);
    let (target_bad_low, target_bad_high) = intermediate_bad_depth_targets(config, pass, target);
    pass_config.in1 = Some(in1.to_path_buf());
    pass_config.in2 = in2.map(Path::to_path_buf);
    pass_config.interleaved = interleaved;
    pass_config.test_interleaved = !interleaved && pass_config.in2.is_none();
    pass_config.extra = if pass == 1 {
        config.extra.clone()
    } else {
        Vec::new()
    };
    pass_config.out1 = Some(out1);
    pass_config.out2 = out2;
    pass_config.out_toss1 = out_toss1;
    pass_config.out_toss2 = out_toss2;
    pass_config.out_low1 = None;
    pass_config.out_low2 = None;
    pass_config.out_mid1 = None;
    pass_config.out_mid2 = None;
    pass_config.out_high1 = None;
    pass_config.out_high2 = None;
    pass_config.out_uncorrected1 = None;
    pass_config.out_uncorrected2 = None;
    pass_config.hist_in = (pass == 1).then(|| config.hist_in.clone()).flatten();
    pass_config.rhist_in = (pass == 1).then(|| config.rhist_in.clone()).flatten();
    pass_config.peaks_in = (pass == 1).then(|| config.peaks_in.clone()).flatten();
    pass_config.match_hist_out = (pass == 1).then(|| config.match_hist_out.clone()).flatten();
    pass_config.insert_hist_out = (pass == 1)
        .then(|| config.insert_hist_out.clone())
        .flatten();
    pass_config.quality_accuracy_hist_out = (pass == 1)
        .then(|| config.quality_accuracy_hist_out.clone())
        .flatten();
    pass_config.indel_hist_out = (pass == 1).then(|| config.indel_hist_out.clone()).flatten();
    pass_config.error_hist_out = (pass == 1).then(|| config.error_hist_out.clone()).flatten();
    pass_config.quality_hist_out = (pass == 1)
        .then(|| config.quality_hist_out.clone())
        .flatten();
    pass_config.base_quality_hist_out = (pass == 1)
        .then(|| config.base_quality_hist_out.clone())
        .flatten();
    pass_config.quality_count_hist_out = (pass == 1)
        .then(|| config.quality_count_hist_out.clone())
        .flatten();
    pass_config.average_quality_hist_out = (pass == 1)
        .then(|| config.average_quality_hist_out.clone())
        .flatten();
    pass_config.overall_base_quality_hist_out = (pass == 1)
        .then(|| config.overall_base_quality_hist_out.clone())
        .flatten();
    pass_config.length_hist_out = (pass == 1)
        .then(|| config.length_hist_out.clone())
        .flatten();
    pass_config.gc_hist_out = (pass == 1).then(|| config.gc_hist_out.clone()).flatten();
    pass_config.base_hist_out = (pass == 1).then(|| config.base_hist_out.clone()).flatten();
    pass_config.entropy_hist_out = (pass == 1)
        .then(|| config.entropy_hist_out.clone())
        .flatten();
    pass_config.identity_hist_out = (pass == 1)
        .then(|| config.identity_hist_out.clone())
        .flatten();
    pass_config.hist_out = None;
    pass_config.rhist_out = None;
    pass_config.peaks_out = None;
    if let Some(bits) = config.count_min_bits_first {
        pass_config.count_min.bits = Some(bits);
    }
    pass_config.target_depth = target;
    pass_config.target_bad_percent_low = target_bad_low as f64 / target as f64;
    pass_config.target_bad_percent_high = target_bad_high as f64 / target as f64;
    pass_config.max_depth = Some(target + target / 4);
    pass_config.min_depth =
        config
            .min_depth
            .min(if config.passes > 2 && pass < config.passes - 1 {
                2
            } else {
                3
            });
    pass_config.min_kmers_over_min_depth = if config.passes > 2 && pass < config.passes - 1 {
        config.min_kmers_over_min_depth.min(5)
    } else {
        config.min_kmers_over_min_depth
    };
    pass_config.depth_percentile = (config.depth_percentile.max(0.4) * 1.2).min(0.8);
    pass_config.toss_error_reads = if config.passes > 2 && pass < config.passes - 1 {
        false
    } else {
        config.toss_error_reads_first
    };
    pass_config.discard_bad_only = if config.passes > 2 && pass < config.passes - 1 {
        true
    } else {
        config.discard_bad_only_first
    };
    pass_config.low_percentile = if config.passes > 2 && pass < config.passes - 1 {
        0.0
    } else {
        config.low_percentile
    };
    pass_config.error_detect_ratio = if config.passes > 2 && pass < config.passes - 1 {
        if config.error_detect_ratio > 100 {
            100 + (config.error_detect_ratio - 100) / 2
        } else {
            config.error_detect_ratio
        }
    } else {
        config.error_detect_ratio
    };
    pass_config.fix_spikes = false;
    pass_config.count_up = false;
    pass_config.error_correct = config.error_correct_first;
    pass_config.overlap_error_correct = config.overlap_error_correct && config.error_correct_first;
    pass_config.rename_reads = false;
    pass_config.overwrite = true;
    pass_config.append = false;
    pass_config.passes = 1;
    pass_config.notes.clear();
    pass_config
}

fn intermediate_target_depth(config: &Config, pass: usize) -> u64 {
    if config.passes > 2 && pass == config.passes - 1 {
        config
            .target_depth_first
            .unwrap_or_else(|| config.target_depth.saturating_mul(2))
    } else if config.passes > 2 {
        config
            .target_depth_first
            .map(|target| target.saturating_mul(2))
            .unwrap_or_else(|| config.target_depth.saturating_mul(4))
    } else {
        config
            .target_depth_first
            .unwrap_or_else(|| config.target_depth.saturating_mul(4))
    }
}

fn intermediate_bad_depth_targets(config: &Config, pass: usize, target: u64) -> (u64, u64) {
    let early_multiplier = if config.passes > 2 && pass < config.passes - 1 {
        1.5
    } else {
        1.0
    };
    let target_f = config.target_depth as f64;
    let low = (target_f * config.target_bad_percent_low * early_multiplier)
        .ceil()
        .max(1.0) as u64;
    let high = (target_f * config.target_bad_percent_high * early_multiplier)
        .ceil()
        .max(1.0) as u64;
    let low = low.min(target);
    let high = high.min(target).max(low);
    (low, high)
}

fn temp_sequence_extension(config: &Config) -> &'static str {
    for path in [
        config.out1.as_ref(),
        config.in1.as_ref(),
        config.out2.as_ref(),
        config.in2.as_ref(),
    ]
    .into_iter()
    .flatten()
    {
        let text = path.to_string_lossy().to_ascii_lowercase();
        if text.ends_with(".fa")
            || text.ends_with(".fasta")
            || text.ends_with(".fna")
            || text.ends_with(".fa.gz")
            || text.ends_with(".fasta.gz")
            || text.ends_with(".fna.gz")
        {
            return "fa";
        }
    }
    "fq"
}

fn cardinality_kmer_config(config: &Config) -> Config {
    let mut cardinality_config = config.clone();
    if let Some(k) = config.cardinality.k {
        cardinality_config.k = k;
    }
    if config.cardinality.min_probability > 0.0 {
        cardinality_config.min_prob = config.cardinality.min_probability;
    }
    cardinality_config
}

fn estimate_primary_cardinality(
    config: &Config,
    cardinality_config: &Config,
) -> Result<CardinalityEstimate> {
    let mut estimator = KmerCardinalityEstimator::from_config(config);
    let mut readers = PrimaryReaders::open(config, config.table_reads)?;
    let mut chunk = Vec::with_capacity(HIST_PARALLEL_CHUNK_SIZE);

    while let Some((r1, r2)) = readers.next_pair()? {
        chunk.push((r1, r2));
        if chunk.len() >= HIST_PARALLEL_CHUNK_SIZE {
            observe_cardinality_chunk(&mut estimator, cardinality_config, &chunk);
            chunk.clear();
        }
    }
    if !chunk.is_empty() {
        observe_cardinality_chunk(&mut estimator, cardinality_config, &chunk);
    }
    Ok(estimator.estimate())
}

fn observe_cardinality_chunk(
    estimator: &mut KmerCardinalityEstimator,
    config: &Config,
    pairs: &[(SequenceRecord, Option<SequenceRecord>)],
) {
    for (r1, r2) in pairs {
        estimator.observe_pair(config, r1, r2.as_ref());
    }
}

fn run_single_pass(config: &Config) -> Result<RunSummary> {
    if config.count_up {
        return run_countup(config);
    }
    let mut stage_timings = Vec::new();
    let cardinality_config = cardinality_kmer_config(config);
    let random_seed = run_random_seed(config);
    let input_counts = build_input_counts_with_stage_timings(config, &mut stage_timings)?;
    let input_cardinality = if config.cardinality.input {
        let started = Instant::now();
        let estimate = estimate_primary_cardinality(config, &cardinality_config)?;
        record_stage_timing(&mut stage_timings, "input_cardinality", started);
        Some(estimate)
    } else {
        None
    };

    let wants_input_hist = config.hist_in.is_some() || config.peaks_in.is_some();
    let wants_input_rhist = config.rhist_in.is_some();
    let mut input_rhist_written_with_hist = false;
    let started = Instant::now();
    if wants_input_hist && wants_input_rhist {
        let (hist, read_hist) =
            collect_primary_sparse_hist_and_read_hist(config, &input_counts, None, random_seed)?;
        if let Some(path) = &config.hist_in {
            write_sparse_depth_hist(path, &hist, config.hist_len, config)?;
        }
        if let Some(path) = &config.peaks_in {
            let dense_hist = sparse_hist_to_peak_dense(&hist, config.hist_len);
            write_peaks(path, &dense_hist, config)?;
        }
        if let Some(path) = &config.rhist_in {
            write_sparse_read_depth_hist(path, &read_hist, config.hist_len, config)?;
            input_rhist_written_with_hist = true;
        }
    } else if wants_input_hist {
        let hist = collect_primary_sparse_hist(config, &input_counts, None, random_seed)?;
        if let Some(path) = &config.hist_in {
            write_sparse_depth_hist(path, &hist, config.hist_len, config)?;
        }
        if let Some(path) = &config.peaks_in {
            let dense_hist = sparse_hist_to_peak_dense(&hist, config.hist_len);
            write_peaks(path, &dense_hist, config)?;
        }
    }
    record_stage_timing(&mut stage_timings, "input_hist", started);

    if input_rhist_written_with_hist {
        record_stage_timing(&mut stage_timings, "input_rhist", Instant::now());
    } else if let Some(path) = &config.rhist_in {
        let started = Instant::now();
        let hist = collect_primary_sparse_read_hist(config, &input_counts, None, random_seed)?;
        write_sparse_read_depth_hist(path, &hist, config.hist_len, config)?;
        record_stage_timing(&mut stage_timings, "input_rhist", started);
    }

    let started = Instant::now();
    emit_read_local_side_outputs(config)?;
    record_stage_timing(&mut stage_timings, "side_outputs", started);

    let started = Instant::now();
    let mut output_counts =
        if config.hist_out.is_some() || config.rhist_out.is_some() || config.peaks_out.is_some() {
            Some(new_output_counts(config)?)
        } else {
            None
        };
    let mut output_cardinality = config
        .cardinality
        .output
        .then(|| KmerCardinalityEstimator::from_config(config));
    record_stage_timing(&mut stage_timings, "output_count_init", started);

    let started = Instant::now();
    let mut summary = normalize_primary(
        config,
        &input_counts,
        output_counts.as_mut(),
        output_cardinality.as_mut(),
        &cardinality_config,
        random_seed,
    )?;
    record_stage_timing(&mut stage_timings, "normalize", started);

    let started = Instant::now();
    (summary.unique_kmers_in, summary.unique_kmers_in_split) = input_counts.unique_kmer_estimate();
    summary.cardinality_in = input_cardinality;
    summary.cardinality_out = output_cardinality
        .as_ref()
        .map(KmerCardinalityEstimator::estimate);
    summary.sketch_layouts = input_counts.sketch_layouts();
    if let Some(counts) = output_counts.as_mut() {
        apply_output_count_adjustments(config, counts);
    }
    summary.unique_kmers_out = output_counts.as_ref().map(CountLookup::unique_kmers);
    if let Some(counts) = output_counts.as_ref() {
        counts.append_sketch_layouts(&mut summary.sketch_layouts, "output_kept");
    }
    record_stage_timing(&mut stage_timings, "summary_counts", started);

    let wants_output_hist = config.hist_out.is_some() || config.peaks_out.is_some();
    let wants_output_rhist = config.rhist_out.is_some();
    let mut output_rhist_written_with_hist = false;
    let started = Instant::now();
    if let Some(counts) = &output_counts {
        if wants_output_hist && wants_output_rhist {
            let (hist, read_hist) = collect_primary_sparse_hist_and_read_hist(
                config,
                counts,
                Some(&input_counts),
                random_seed,
            )?;
            if let Some(path) = &config.hist_out {
                write_sparse_depth_hist(path, &hist, config.hist_len, config)?;
            }
            if let Some(path) = &config.peaks_out {
                let dense_hist = sparse_hist_to_peak_dense(&hist, config.hist_len);
                write_peaks(path, &dense_hist, config)?;
            }
            if let Some(path) = &config.rhist_out {
                write_sparse_read_depth_hist(path, &read_hist, config.hist_len, config)?;
                output_rhist_written_with_hist = true;
            }
        } else if wants_output_hist {
            let hist =
                collect_primary_sparse_hist(config, counts, Some(&input_counts), random_seed)?;
            if let Some(path) = &config.hist_out {
                write_sparse_depth_hist(path, &hist, config.hist_len, config)?;
            }
            if let Some(path) = &config.peaks_out {
                let dense_hist = sparse_hist_to_peak_dense(&hist, config.hist_len);
                write_peaks(path, &dense_hist, config)?;
            }
        }
    }
    record_stage_timing(&mut stage_timings, "output_hist", started);

    if output_rhist_written_with_hist {
        record_stage_timing(&mut stage_timings, "output_rhist", Instant::now());
    } else if let (Some(path), Some(counts)) = (&config.rhist_out, &output_counts) {
        let started = Instant::now();
        let hist =
            collect_primary_sparse_read_hist(config, counts, Some(&input_counts), random_seed)?;
        write_sparse_read_depth_hist(path, &hist, config.hist_len, config)?;
        record_stage_timing(&mut stage_timings, "output_rhist", started);
    }

    summary.stage_timings = stage_timings;
    Ok(summary)
}

fn run_countup(config: &Config) -> Result<RunSummary> {
    let mut stage_timings = Vec::new();
    let cardinality_config = cardinality_kmer_config(config);
    let input_counts = build_input_counts_with_stage_timings(config, &mut stage_timings)?;
    let input_cardinality = if config.cardinality.input {
        let started = Instant::now();
        let estimate = estimate_primary_cardinality(config, &cardinality_config)?;
        record_stage_timing(&mut stage_timings, "input_cardinality", started);
        Some(estimate)
    } else {
        None
    };

    let wants_input_hist = config.hist_in.is_some() || config.peaks_in.is_some();
    let wants_input_rhist = config.rhist_in.is_some();

    let started = Instant::now();
    emit_read_local_side_outputs(config)?;
    record_stage_timing(&mut stage_timings, "side_outputs", started);

    let random_seed = run_random_seed(config);
    let started = Instant::now();
    let work_build = collect_countup_work_source(
        config,
        &input_counts,
        random_seed,
        wants_input_hist,
        wants_input_rhist,
    )?;
    let countup_work_elapsed = started
        .elapsed()
        .as_micros()
        .saturating_sub(work_build.input_hist_elapsed_micros);
    record_stage_timing_micros(
        &mut stage_timings,
        "input_hist",
        work_build.input_hist_elapsed_micros,
    );
    if let (Some(path), Some(hist)) = (&config.hist_in, &work_build.input_hist) {
        write_sparse_depth_hist(path, hist, config.hist_len, config)?;
    }
    if let (Some(path), Some(hist)) = (&config.peaks_in, &work_build.input_hist) {
        let dense_hist = sparse_hist_to_peak_dense(hist, config.hist_len);
        write_peaks(path, &dense_hist, config)?;
    }
    if let (Some(path), Some(hist)) = (&config.rhist_in, &work_build.input_read_hist) {
        write_sparse_read_depth_hist(path, hist, config.hist_len, config)?;
    }
    record_stage_timing_micros(&mut stage_timings, "input_rhist", 0);
    record_stage_timing_micros(
        &mut stage_timings,
        "countup_work_source",
        countup_work_elapsed,
    );
    let format1 = work_build.format1;
    let format2 = work_build.format2;
    let countup_spill = work_build.spill_summary;
    let mut work_pairs = work_build.source.into_iter()?;
    let mut writers = OptionalWriters::open(config, format1, format2)?;
    let mut summary = RunSummary {
        cardinality_in: input_cardinality,
        countup_spill,
        ..RunSummary::default()
    };
    let mut kept_counts = new_output_counts(config)?;
    let mut output_cardinality = config
        .cardinality
        .output
        .then(|| KmerCardinalityEstimator::from_config(config));
    let adjusted_target = ((config.target_depth as f64) * 0.95).round().max(1.0) as u64;
    let started = Instant::now();

    while let Some(CountupWorkPair {
        input_list_index,
        mut r1,
        mut r2,
        ..
    }) = work_pairs.next_pair()?
    {
        writers.sync_to_input_list_index(config, input_list_index)?;
        let keys = unique_pair_kmers(config, &r1, r2.as_ref());
        let mut decision_plan =
            countup_decision_plan(config, &input_counts, &kept_counts, &keys, adjusted_target);
        if countup_length_toss(config, &r1, r2.as_ref()) {
            decision_plan.toss = true;
        }
        update_countup_kept_counts_for_plan(config, &mut kept_counts, &keys, &decision_plan);

        let output_analysis = needs_output_pair_analysis(config)
            .then(|| analyze_pair(config, &input_counts, &r1, r2.as_ref()));
        let mut correction = CorrectionResult::default();
        if config.error_correct && !decision_plan.toss {
            correction =
                correct_pair_errors_with_rollback(config, &input_counts, &mut r1, r2.as_mut());
        }
        if config.trim_after_marking && config.error_correct {
            trim_pair(config, &mut r1, r2.as_mut());
        }
        let (out_r1, out_r2) = match output_analysis.as_ref() {
            Some(analysis) => maybe_rename_pair(config, &r1, r2.as_ref(), analysis),
            None => (r1.clone(), r2.clone()),
        };
        let read_count = 1 + u64::from(r2.is_some());
        let base_count = r1.len() as u64 + r2.as_ref().map(|r| r.len() as u64).unwrap_or(0);

        summary.reads_in += read_count;
        summary.bases_in += base_count;
        if decision_plan.toss {
            summary.reads_tossed += read_count;
            summary.bases_tossed += base_count;
        } else {
            summary.reads_kept += read_count;
            summary.bases_kept += base_count;
            if let Some(estimator) = output_cardinality.as_mut() {
                estimator.observe_pair(&cardinality_config, &r1, r2.as_ref());
            }
        }
        writers.write_pair(decision_plan.toss, &out_r1, out_r2.as_ref())?;
        if correction.uncorrectable {
            writers.write_uncorrected(&r1, r2.as_ref())?;
        }
        if let Some(analysis) = output_analysis.as_ref()
            && depth_bin_outputs_enabled(config)
        {
            writers.write_depth_bin(config, analysis, &out_r1, out_r2.as_ref())?;
        }
    }
    writers.flush()?;
    record_stage_timing(&mut stage_timings, "countup_normalize", started);

    let started = Instant::now();
    if config.hist_out.is_some() || config.peaks_out.is_some() || config.rhist_out.is_some() {
        apply_output_count_adjustments(config, &mut kept_counts);
    }
    record_stage_timing(&mut stage_timings, "output_count_adjust", started);

    let started = Instant::now();
    let output_hist = if config.hist_out.is_some() || config.peaks_out.is_some() {
        Some(kept_counts.sparse_depth_hist(config.hist_len))
    } else {
        None
    };
    if let (Some(path), Some(hist)) = (&config.hist_out, &output_hist) {
        write_sparse_depth_hist(path, hist, config.hist_len, config)?;
    }
    if let (Some(path), Some(hist)) = (&config.peaks_out, &output_hist) {
        let dense_hist = sparse_hist_to_peak_dense(hist, config.hist_len);
        write_peaks(path, &dense_hist, config)?;
    }
    record_stage_timing(&mut stage_timings, "output_hist", started);

    if let Some(path) = &config.rhist_out {
        let started = Instant::now();
        let hist = collect_primary_sparse_read_hist(config, &kept_counts, Some(&input_counts), 0)?;
        write_sparse_read_depth_hist(path, &hist, config.hist_len, config)?;
        record_stage_timing(&mut stage_timings, "output_rhist", started);
    }

    let started = Instant::now();
    (summary.unique_kmers_in, summary.unique_kmers_in_split) = input_counts.unique_kmer_estimate();
    summary.unique_kmers_out = Some(kept_counts.unique_kmers());
    summary.cardinality_out = output_cardinality
        .as_ref()
        .map(KmerCardinalityEstimator::estimate);
    summary.sketch_layouts = input_counts.sketch_layouts();
    kept_counts.append_sketch_layouts(&mut summary.sketch_layouts, "countup_kept");
    record_stage_timing(&mut stage_timings, "summary_counts", started);
    summary.stage_timings = stage_timings;
    Ok(summary)
}

fn record_stage_timing(timings: &mut Vec<StageTiming>, name: &'static str, started: Instant) {
    timings.push(StageTiming {
        name,
        elapsed_micros: started.elapsed().as_micros(),
    });
}

fn record_stage_timing_micros(
    timings: &mut Vec<StageTiming>,
    name: &'static str,
    elapsed_micros: u128,
) {
    timings.push(StageTiming {
        name,
        elapsed_micros,
    });
}

fn collect_countup_work_source(
    config: &Config,
    input_counts: &dyn CountLookup,
    random_seed: u64,
    wants_input_hist: bool,
    wants_input_rhist: bool,
) -> Result<CountupWorkBuild> {
    let mut readers = PrimaryReaders::open(config, config.max_reads)?;
    let format1 = readers.format1();
    let format2 = readers.format2();
    let presort_config = countup_prepass_config(config);
    let mut rng = JavaXoshiro::new(random_seed);
    let mut work_pairs = Vec::new();
    let mut work_pair_bytes = 0usize;
    let mut run_paths = Vec::new();
    let mut temp_dir = None;
    let mut spill_summary = CountupSpillSummary::default();
    let mut input_hist = wants_input_hist.then(SparseHist::default);
    let mut input_read_hist = wants_input_rhist.then(SparseReadDepthHist::default);
    let mut input_hist_elapsed_micros = 0u128;
    let mut candidates = Vec::with_capacity(COUNTUP_PREPASS_CHUNK_PAIR_LIMIT);
    let mut candidate_bytes = 0usize;
    let mut original_index = 0usize;
    while let Some((r1, r2)) = readers.next_pair()? {
        let candidate = CountupWorkCandidate {
            input_list_index: readers.input_list_index(),
            original_index,
            rand: rng.next_double(),
            r1,
            r2,
        };
        candidate_bytes =
            candidate_bytes.saturating_add(countup_work_candidate_memory_hint(&candidate));
        candidates.push(candidate);
        if countup_prepass_chunk_ready(candidates.len(), candidate_bytes) {
            let chunk = std::mem::take(&mut candidates);
            let chunk_build = process_countup_work_candidate_chunk(
                config,
                &presort_config,
                input_counts,
                wants_input_hist,
                wants_input_rhist,
                chunk,
            );
            let hist_started = Instant::now();
            if let Some(input_hist) = input_hist.as_mut() {
                merge_sparse_hist(input_hist, chunk_build.depth_hist);
            }
            if let Some(input_read_hist) = input_read_hist.as_mut() {
                merge_sparse_read_depth_hist(input_read_hist, chunk_build.read_hist);
            }
            input_hist_elapsed_micros =
                input_hist_elapsed_micros.saturating_add(hist_started.elapsed().as_micros());
            append_countup_work_pairs(
                config,
                &mut temp_dir,
                &mut run_paths,
                &mut spill_summary,
                &mut work_pairs,
                &mut work_pair_bytes,
                chunk_build.work_pairs,
            )?;
            candidates = Vec::with_capacity(COUNTUP_PREPASS_CHUNK_PAIR_LIMIT);
            candidate_bytes = 0;
        }
        original_index += 1;
    }
    if !candidates.is_empty() {
        let chunk_build = process_countup_work_candidate_chunk(
            config,
            &presort_config,
            input_counts,
            wants_input_hist,
            wants_input_rhist,
            candidates,
        );
        let hist_started = Instant::now();
        if let Some(input_hist) = input_hist.as_mut() {
            merge_sparse_hist(input_hist, chunk_build.depth_hist);
        }
        if let Some(input_read_hist) = input_read_hist.as_mut() {
            merge_sparse_read_depth_hist(input_read_hist, chunk_build.read_hist);
        }
        input_hist_elapsed_micros =
            input_hist_elapsed_micros.saturating_add(hist_started.elapsed().as_micros());
        append_countup_work_pairs(
            config,
            &mut temp_dir,
            &mut run_paths,
            &mut spill_summary,
            &mut work_pairs,
            &mut work_pair_bytes,
            chunk_build.work_pairs,
        )?;
    }
    let source = if run_paths.is_empty() {
        work_pairs.sort_by(compare_countup_work_pairs);
        CountupWorkSource {
            temp_dir: None,
            inner: CountupWorkSourceInner::Memory(work_pairs),
        }
    } else {
        if !work_pairs.is_empty() {
            spill_countup_run(
                config,
                &mut temp_dir,
                &mut run_paths,
                &mut spill_summary,
                &mut work_pairs,
            )?;
        }
        compact_countup_runs(config, &mut run_paths, &mut spill_summary)?;
        spill_summary.final_runs = run_paths.len();
        enforce_countup_spill_limits(config, &spill_summary, run_paths.len())?;
        CountupWorkSource {
            temp_dir,
            inner: CountupWorkSourceInner::Spilled(run_paths),
        }
    };
    Ok(CountupWorkBuild {
        source,
        input_hist,
        input_read_hist,
        input_hist_elapsed_micros,
        format1,
        format2,
        spill_summary,
    })
}

fn countup_prepass_chunk_ready(candidate_count: usize, candidate_bytes: usize) -> bool {
    candidate_count >= COUNTUP_PREPASS_CHUNK_PAIR_LIMIT
        || candidate_bytes >= COUNTUP_PREPASS_CHUNK_BYTE_LIMIT
}

fn process_countup_work_candidates(
    config: &Config,
    presort_config: &Config,
    input_counts: &dyn CountLookup,
    candidates: Vec<CountupWorkCandidate>,
) -> Vec<CountupWorkPair> {
    candidates
        .into_par_iter()
        .filter_map(|candidate| {
            countup_work_pair_from_candidate(config, presort_config, input_counts, candidate)
        })
        .collect()
}

fn process_countup_work_candidate_chunk(
    config: &Config,
    presort_config: &Config,
    input_counts: &dyn CountLookup,
    wants_depth_hist: bool,
    wants_read_hist: bool,
    candidates: Vec<CountupWorkCandidate>,
) -> CountupChunkBuild {
    if !wants_depth_hist && !wants_read_hist {
        return CountupChunkBuild {
            work_pairs: process_countup_work_candidates(
                config,
                presort_config,
                input_counts,
                candidates,
            ),
            depth_hist: SparseHist::default(),
            read_hist: SparseReadDepthHist::default(),
        };
    }

    candidates
        .into_par_iter()
        .fold(
            || CountupChunkBuild {
                work_pairs: Vec::new(),
                depth_hist: SparseHist::default(),
                read_hist: SparseReadDepthHist::default(),
            },
            |mut local, candidate| {
                let mut hist = CountupInputHistAccumulator {
                    wants_depth_hist,
                    wants_read_hist,
                    depth_hist: &mut local.depth_hist,
                    read_hist: &mut local.read_hist,
                };
                if let Some(work_pair) = countup_work_pair_from_candidate_with_input_hists(
                    config,
                    presort_config,
                    input_counts,
                    candidate,
                    &mut hist,
                ) {
                    local.work_pairs.push(work_pair);
                }
                local
            },
        )
        .reduce(
            || CountupChunkBuild {
                work_pairs: Vec::new(),
                depth_hist: SparseHist::default(),
                read_hist: SparseReadDepthHist::default(),
            },
            |mut left, mut right| {
                left.work_pairs.append(&mut right.work_pairs);
                merge_sparse_hist(&mut left.depth_hist, right.depth_hist);
                merge_sparse_read_depth_hist(&mut left.read_hist, right.read_hist);
                left
            },
        )
}

fn countup_work_pair_from_candidate(
    config: &Config,
    presort_config: &Config,
    input_counts: &dyn CountLookup,
    mut candidate: CountupWorkCandidate,
) -> Option<CountupWorkPair> {
    if !config.trim_after_marking {
        trim_pair(config, &mut candidate.r1, candidate.r2.as_mut());
    }
    let prepass_result = countup_prepass_pair(
        presort_config,
        config.add_bad_reads_countup,
        input_counts,
        &mut candidate.r1,
        candidate.r2.as_mut(),
        candidate.rand,
    );
    countup_work_pair_from_prepass_result(presort_config, input_counts, candidate, prepass_result)
}

fn countup_work_pair_from_candidate_with_input_hists(
    config: &Config,
    presort_config: &Config,
    input_counts: &dyn CountLookup,
    mut candidate: CountupWorkCandidate,
    hist: &mut CountupInputHistAccumulator<'_>,
) -> Option<CountupWorkPair> {
    if config.trim_after_marking {
        let mut hist_r1 = candidate.r1.clone();
        let mut hist_r2 = candidate.r2.clone();
        trim_pair(config, &mut hist_r1, hist_r2.as_mut());
        let hist_analysis = analyze_pair(config, input_counts, &hist_r1, hist_r2.as_ref());
        increment_countup_input_hists_from_analysis(
            config,
            hist,
            &hist_r1,
            hist_r2.as_ref(),
            &hist_analysis,
        );
    } else {
        trim_pair(config, &mut candidate.r1, candidate.r2.as_mut());
        let (hist_analysis, prepass_analysis) = analyze_pair_for_two_configs(
            config,
            presort_config,
            input_counts,
            &candidate.r1,
            candidate.r2.as_ref(),
        );
        increment_countup_input_hists_from_analysis(
            config,
            hist,
            &candidate.r1,
            candidate.r2.as_ref(),
            &hist_analysis,
        );
        let prepass_result = countup_prepass_pair_from_analysis(
            presort_config,
            config.add_bad_reads_countup,
            input_counts,
            &mut candidate.r1,
            candidate.r2.as_mut(),
            candidate.rand,
            prepass_analysis,
        );
        return countup_work_pair_from_prepass_result(
            presort_config,
            input_counts,
            candidate,
            prepass_result,
        );
    }

    if !config.trim_after_marking {
        trim_pair(config, &mut candidate.r1, candidate.r2.as_mut());
    }
    let prepass_result = countup_prepass_pair(
        presort_config,
        config.add_bad_reads_countup,
        input_counts,
        &mut candidate.r1,
        candidate.r2.as_mut(),
        candidate.rand,
    );
    countup_work_pair_from_prepass_result(presort_config, input_counts, candidate, prepass_result)
}

fn increment_countup_input_hists_from_analysis(
    config: &Config,
    hist: &mut CountupInputHistAccumulator<'_>,
    r1: &SequenceRecord,
    r2: Option<&SequenceRecord>,
    analysis: &PairAnalysis,
) {
    if hist.wants_depth_hist {
        increment_sparse_hist_from_analysis(hist.depth_hist, &analysis.read1, config.hist_len);
        if let Some(read2_analysis) = &analysis.read2 {
            increment_sparse_hist_from_analysis(hist.depth_hist, read2_analysis, config.hist_len);
        }
    }
    if hist.wants_read_hist {
        increment_sparse_read_hist(hist.read_hist, &analysis.read1, r1.len(), config.hist_len);
        if let (Some(read2_analysis), Some(read2)) = (&analysis.read2, r2) {
            increment_sparse_read_hist(
                hist.read_hist,
                read2_analysis,
                read2.len(),
                config.hist_len,
            );
        }
    }
}

fn countup_work_pair_from_prepass_result(
    presort_config: &Config,
    input_counts: &dyn CountLookup,
    candidate: CountupWorkCandidate,
    prepass_result: CountupPrepassResult,
) -> Option<CountupWorkPair> {
    if !prepass_result.include {
        return None;
    }
    let sort_key = prepass_result.sort_analysis.as_ref().map_or_else(
        || {
            countup_sort_key(
                presort_config,
                input_counts,
                &candidate.r1,
                candidate.r2.as_ref(),
                candidate.original_index,
            )
        },
        |analysis| {
            countup_sort_key_from_analysis(
                &candidate.r1,
                candidate.r2.as_ref(),
                candidate.original_index,
                analysis,
            )
        },
    );
    Some(CountupWorkPair {
        input_list_index: candidate.input_list_index,
        sort_key,
        r1: candidate.r1,
        r2: candidate.r2,
    })
}

fn append_countup_work_pairs(
    config: &Config,
    temp_dir: &mut Option<tempfile::TempDir>,
    run_paths: &mut Vec<PathBuf>,
    spill_summary: &mut CountupSpillSummary,
    work_pairs: &mut Vec<CountupWorkPair>,
    work_pair_bytes: &mut usize,
    new_pairs: Vec<CountupWorkPair>,
) -> Result<()> {
    for work_pair in new_pairs {
        *work_pair_bytes =
            (*work_pair_bytes).saturating_add(countup_work_pair_memory_hint(&work_pair));
        work_pairs.push(work_pair);
        if work_pairs.len() >= COUNTUP_SORT_RUN_PAIR_LIMIT
            || *work_pair_bytes >= COUNTUP_SORT_RUN_BYTE_LIMIT
        {
            spill_countup_run(config, temp_dir, run_paths, spill_summary, work_pairs)?;
            *work_pair_bytes = 0;
        }
    }
    Ok(())
}

fn countup_work_pair_memory_hint(pair: &CountupWorkPair) -> usize {
    std::mem::size_of::<CountupWorkPair>()
        .saturating_add(countup_sort_key_memory_hint(&pair.sort_key))
        .saturating_add(sequence_record_memory_hint(&pair.r1))
        .saturating_add(pair.r2.as_ref().map_or(0, sequence_record_memory_hint))
}

fn countup_work_candidate_memory_hint(candidate: &CountupWorkCandidate) -> usize {
    std::mem::size_of::<CountupWorkCandidate>()
        .saturating_add(sequence_record_memory_hint(&candidate.r1))
        .saturating_add(candidate.r2.as_ref().map_or(0, sequence_record_memory_hint))
}

fn countup_sort_key_memory_hint(key: &CountupSortKey) -> usize {
    let _ = key;
    std::mem::size_of::<CountupSortKey>()
}

fn sequence_record_memory_hint(record: &SequenceRecord) -> usize {
    std::mem::size_of::<SequenceRecord>()
        .saturating_add(record.id.capacity())
        .saturating_add(record.bases.capacity())
        .saturating_add(record.qualities.as_ref().map_or(0, Vec::capacity))
}

fn spill_countup_run(
    config: &Config,
    temp_dir: &mut Option<tempfile::TempDir>,
    run_paths: &mut Vec<PathBuf>,
    spill_summary: &mut CountupSpillSummary,
    work_pairs: &mut Vec<CountupWorkPair>,
) -> Result<()> {
    if work_pairs.is_empty() {
        return Ok(());
    }
    let dir = match temp_dir {
        Some(dir) => dir,
        None => temp_dir.insert(managed_temp_dir(config, "bbnorm-rs-countup-")?),
    };
    work_pairs.sort_by(compare_countup_work_pairs);
    let path = dir
        .path()
        .join(format!("countup-run-{:06}.bin", run_paths.len()));
    let bytes = write_countup_run(&path, work_pairs)?;
    spill_summary.note_initial_run(bytes);
    run_paths.push(path);
    enforce_countup_spill_limits(config, spill_summary, run_paths.len())?;
    work_pairs.clear();
    Ok(())
}

fn compact_countup_runs(
    config: &Config,
    run_paths: &mut Vec<PathBuf>,
    spill_summary: &mut CountupSpillSummary,
) -> Result<()> {
    if run_paths.len() <= COUNTUP_SORT_MERGE_FANIN {
        return Ok(());
    }
    let run_dir = run_paths
        .first()
        .and_then(|path| path.parent())
        .context("count-up spill runs had no parent directory")?
        .to_path_buf();
    let mut round = 0usize;
    while run_paths.len() > COUNTUP_SORT_MERGE_FANIN {
        let old_paths = std::mem::take(run_paths);
        for (group_index, group) in old_paths.chunks(COUNTUP_SORT_MERGE_FANIN).enumerate() {
            let merged_path =
                run_dir.join(format!("countup-merge-{round:03}-{group_index:06}.bin"));
            let merged_bytes = merge_countup_run_group(group, &merged_path)?;
            spill_summary.note_merge_run(merged_bytes);
            run_paths.push(merged_path);
            enforce_countup_spill_limits(config, spill_summary, run_paths.len())?;
        }
        for path in old_paths {
            let removed_bytes = path.metadata().map(|metadata| metadata.len()).unwrap_or(0);
            match fs::remove_file(&path) {
                Ok(()) => spill_summary.note_removed(removed_bytes),
                Err(err) if err.kind() == ErrorKind::NotFound => {}
                Err(err) => {
                    return Err(err).with_context(|| {
                        format!("removing compacted count-up run {}", path.display())
                    });
                }
            }
        }
        round += 1;
    }
    Ok(())
}

fn enforce_countup_spill_limits(
    config: &Config,
    spill_summary: &CountupSpillSummary,
    live_run_count: usize,
) -> Result<()> {
    if let Some(limit) = config.max_countup_spill_initial_runs
        && spill_summary.initial_runs > limit
    {
        bail!(
            "count-up spill exceeded maxcountupspillinitialruns: initial spill runs {} > limit {}",
            spill_summary.initial_runs,
            limit
        );
    }
    if let Some(limit) = config.max_countup_spill_merge_runs
        && spill_summary.merge_runs > limit
    {
        bail!(
            "count-up spill exceeded maxcountupspillmergeruns: merge spill runs {} > limit {}",
            spill_summary.merge_runs,
            limit
        );
    }
    if let Some(limit) = config.max_countup_spill_final_runs
        && live_run_count > limit
    {
        bail!(
            "count-up spill exceeded maxcountupspillfinalruns: live spill runs {} > limit {}",
            live_run_count,
            limit
        );
    }
    if let Some(limit) = config.max_countup_spill_live_bytes
        && spill_summary.peak_live_bytes > limit
    {
        bail!(
            "count-up spill exceeded maxcountupspillbytes: peak live spill bytes {} > limit {}",
            spill_summary.peak_live_bytes,
            limit
        );
    }
    if let Some(limit) = config.max_countup_spill_final_live_bytes
        && spill_summary.final_live_bytes > limit
    {
        bail!(
            "count-up spill exceeded maxcountupspillfinallivebytes: current/final live spill bytes {} > limit {}",
            spill_summary.final_live_bytes,
            limit
        );
    }
    if let Some(limit) = config.max_countup_spill_write_bytes
        && spill_summary.bytes_written > limit
    {
        bail!(
            "count-up spill exceeded maxcountupspillwritebytes: cumulative spill bytes written {} > limit {}",
            spill_summary.bytes_written,
            limit
        );
    }
    Ok(())
}

fn merge_countup_run_group(paths: &[PathBuf], output_path: &Path) -> Result<u64> {
    let mut merger = CountupRunMerger::new(paths)?;
    let file = fs::File::create(output_path)
        .with_context(|| format!("creating compacted count-up run {}", output_path.display()))?;
    let mut writer = BufWriter::with_capacity(COUNTUP_RUN_IO_BUFFER_CAPACITY, file);
    while let Some(pair) = merger.next_pair()? {
        write_countup_work_pair(&mut writer, &pair)?;
    }
    writer
        .flush()
        .with_context(|| format!("flushing compacted count-up run {}", output_path.display()))?;
    output_path
        .metadata()
        .map(|metadata| metadata.len())
        .with_context(|| format!("checking compacted count-up run {}", output_path.display()))
}

impl CountupWorkSource {
    fn into_iter(self) -> Result<CountupWorkIter> {
        let CountupWorkSource { temp_dir, inner } = self;
        let inner = match inner {
            CountupWorkSourceInner::Memory(work_pairs) => {
                CountupWorkIterInner::Memory(work_pairs.into_iter())
            }
            CountupWorkSourceInner::Spilled(paths) => {
                CountupWorkIterInner::Spilled(CountupRunMerger::new(&paths)?)
            }
        };
        Ok(CountupWorkIter {
            _temp_dir: temp_dir,
            inner,
        })
    }
}

impl CountupWorkIter {
    fn next_pair(&mut self) -> Result<Option<CountupWorkPair>> {
        match &mut self.inner {
            CountupWorkIterInner::Memory(iter) => Ok(iter.next()),
            CountupWorkIterInner::Spilled(merger) => merger.next_pair(),
        }
    }
}

impl CountupRunMerger {
    fn new(paths: &[PathBuf]) -> Result<Self> {
        let mut readers = Vec::with_capacity(paths.len());
        let mut heap = BinaryHeap::new();
        for path in paths {
            let mut reader = CountupRunReader::open(path)?;
            if let Some(pair) = reader.next_pair()? {
                heap.push(CountupRunHead {
                    pair,
                    run_index: readers.len(),
                });
            }
            readers.push(reader);
        }
        Ok(Self { readers, heap })
    }

    fn next_pair(&mut self) -> Result<Option<CountupWorkPair>> {
        let Some(head) = self.heap.pop() else {
            return Ok(None);
        };
        let pair = head.pair;
        if let Some(next) = self.readers[head.run_index].next_pair()? {
            self.heap.push(CountupRunHead {
                pair: next,
                run_index: head.run_index,
            });
        }
        Ok(Some(pair))
    }
}

impl CountupRunReader {
    fn open(path: &Path) -> Result<Self> {
        let file = fs::File::open(path)
            .with_context(|| format!("opening count-up run {}", path.display()))?;
        Ok(Self {
            reader: BufReader::with_capacity(COUNTUP_RUN_IO_BUFFER_CAPACITY, file),
        })
    }

    fn next_pair(&mut self) -> Result<Option<CountupWorkPair>> {
        read_countup_work_pair(&mut self.reader)
    }
}

impl PartialEq for CountupRunHead {
    fn eq(&self, other: &Self) -> bool {
        compare_countup_work_pairs(&self.pair, &other.pair) == CmpOrdering::Equal
            && self.run_index == other.run_index
    }
}

impl Eq for CountupRunHead {}

impl PartialOrd for CountupRunHead {
    fn partial_cmp(&self, other: &Self) -> Option<CmpOrdering> {
        Some(self.cmp(other))
    }
}

impl Ord for CountupRunHead {
    fn cmp(&self, other: &Self) -> CmpOrdering {
        compare_countup_work_pairs(&other.pair, &self.pair)
            .then_with(|| other.run_index.cmp(&self.run_index))
    }
}

fn write_countup_run(path: &Path, work_pairs: &[CountupWorkPair]) -> Result<u64> {
    let file = fs::File::create(path)
        .with_context(|| format!("creating count-up run {}", path.display()))?;
    let mut writer = BufWriter::with_capacity(COUNTUP_RUN_IO_BUFFER_CAPACITY, file);
    for pair in work_pairs {
        write_countup_work_pair(&mut writer, pair)?;
    }
    writer
        .flush()
        .with_context(|| format!("flushing count-up run {}", path.display()))?;
    path.metadata()
        .map(|metadata| metadata.len())
        .with_context(|| format!("checking count-up run {}", path.display()))
}

fn write_countup_work_pair(writer: &mut impl Write, pair: &CountupWorkPair) -> Result<()> {
    write_usize(writer, pair.input_list_index)?;
    write_usize(writer, pair.sort_key.errors)?;
    write_usize(writer, pair.sort_key.total_len)?;
    writer.write_all(&pair.sort_key.expected_errors.to_le_bytes())?;
    writer.write_all(&pair.sort_key.numeric_id.to_le_bytes())?;
    write_usize(writer, pair.sort_key.original_index)?;
    write_sequence_record(writer, &pair.r1)?;
    write_bool(writer, pair.r2.is_some())?;
    if let Some(r2) = &pair.r2 {
        write_sequence_record(writer, r2)?;
    }
    Ok(())
}

fn read_countup_work_pair(reader: &mut impl Read) -> Result<Option<CountupWorkPair>> {
    let Some(input_list_index) = read_usize_opt(reader)? else {
        return Ok(None);
    };
    let errors = read_usize(reader)?;
    let total_len = read_usize(reader)?;
    let expected_errors = read_f64(reader)?;
    let numeric_id = read_u64(reader)?;
    let original_index = read_usize(reader)?;
    let r1 = read_sequence_record(reader)?;
    let has_r2 = read_bool(reader)?;
    let r2 = has_r2.then(|| read_sequence_record(reader)).transpose()?;
    Ok(Some(CountupWorkPair {
        input_list_index,
        sort_key: CountupSortKey {
            errors,
            total_len,
            expected_errors,
            numeric_id,
            original_index,
        },
        r1,
        r2,
    }))
}

fn write_sequence_record(writer: &mut impl Write, record: &SequenceRecord) -> Result<()> {
    write_string(writer, &record.id)?;
    writer.write_all(&record.numeric_id.to_le_bytes())?;
    write_bytes(writer, &record.bases)?;
    write_bool(writer, record.qualities.is_some())?;
    if let Some(qualities) = &record.qualities {
        write_bytes(writer, qualities)?;
    }
    Ok(())
}

fn read_sequence_record(reader: &mut impl Read) -> Result<SequenceRecord> {
    let id = read_string(reader)?;
    let numeric_id = read_u64(reader)?;
    let bases = read_bytes(reader)?;
    let has_qualities = read_bool(reader)?;
    let qualities = has_qualities.then(|| read_bytes(reader)).transpose()?;
    Ok(SequenceRecord {
        id,
        numeric_id,
        bases,
        qualities,
    })
}

fn write_string(writer: &mut impl Write, value: &str) -> Result<()> {
    write_bytes(writer, value.as_bytes())
}

fn read_string(reader: &mut impl Read) -> Result<String> {
    let bytes = read_bytes(reader)?;
    String::from_utf8(bytes).context("count-up run contained invalid UTF-8 id")
}

fn write_bytes(writer: &mut impl Write, bytes: &[u8]) -> Result<()> {
    write_usize(writer, bytes.len())?;
    writer.write_all(bytes)?;
    Ok(())
}

fn read_bytes(reader: &mut impl Read) -> Result<Vec<u8>> {
    let len = read_usize(reader)?;
    let mut bytes = vec![0; len];
    reader.read_exact(&mut bytes)?;
    Ok(bytes)
}

fn write_bool(writer: &mut impl Write, value: bool) -> Result<()> {
    writer.write_all(&[u8::from(value)])?;
    Ok(())
}

fn read_bool(reader: &mut impl Read) -> Result<bool> {
    let mut buf = [0; 1];
    reader.read_exact(&mut buf)?;
    Ok(buf[0] != 0)
}

fn write_usize(writer: &mut impl Write, value: usize) -> Result<()> {
    writer.write_all(&(value as u64).to_le_bytes())?;
    Ok(())
}

fn read_usize(reader: &mut impl Read) -> Result<usize> {
    let value = read_u64(reader)?;
    usize::try_from(value).context("count-up run usize field exceeded this platform")
}

fn read_usize_opt(reader: &mut impl Read) -> Result<Option<usize>> {
    let Some(value) = read_u64_opt(reader)? else {
        return Ok(None);
    };
    Ok(Some(usize::try_from(value).context(
        "count-up run usize field exceeded this platform",
    )?))
}

fn read_u64(reader: &mut impl Read) -> Result<u64> {
    let mut buf = [0; 8];
    reader.read_exact(&mut buf)?;
    Ok(u64::from_le_bytes(buf))
}

fn read_u64_opt(reader: &mut impl Read) -> Result<Option<u64>> {
    let mut buf = [0; 8];
    match reader.read_exact(&mut buf) {
        Ok(()) => Ok(Some(u64::from_le_bytes(buf))),
        Err(err) if err.kind() == ErrorKind::UnexpectedEof => Ok(None),
        Err(err) => Err(err.into()),
    }
}

fn read_f64(reader: &mut impl Read) -> Result<f64> {
    let mut buf = [0; 8];
    reader.read_exact(&mut buf)?;
    Ok(f64::from_le_bytes(buf))
}

fn countup_prepass_config(config: &Config) -> Config {
    let mut prepass = config.clone();
    prepass.count_up = false;
    // Java sets REQUIRE_BOTH_BAD=(rbb || COUNTUP) before temporarily
    // disabling COUNTUP for the relaxed count-up prepass.
    prepass.require_both_bad = true;
    prepass.target_depth = config.target_depth.saturating_mul(4).max(1);
    prepass.target_bad_percent_low = config.target_bad_percent_low / 4.0;
    prepass.target_bad_percent_high = config.target_bad_percent_high / 4.0;
    prepass.max_depth = config.max_depth.map(|depth| depth.saturating_mul(4).max(1));
    prepass.min_depth = config.min_depth / 2;
    prepass.min_kmers_over_min_depth = config.min_kmers_over_min_depth / 2;
    prepass.low_percentile = 0.20;
    prepass
}

fn countup_prepass_pair(
    prepass_config: &Config,
    add_bad_reads_countup: bool,
    input_counts: &dyn CountLookup,
    r1: &mut SequenceRecord,
    r2: Option<&mut SequenceRecord>,
    rand: f64,
) -> CountupPrepassResult {
    let analysis = analyze_pair(prepass_config, input_counts, r1, r2.as_deref());
    countup_prepass_pair_from_analysis(
        prepass_config,
        add_bad_reads_countup,
        input_counts,
        r1,
        r2,
        rand,
        analysis,
    )
}

fn countup_prepass_pair_from_analysis(
    prepass_config: &Config,
    add_bad_reads_countup: bool,
    input_counts: &dyn CountLookup,
    r1: &mut SequenceRecord,
    mut r2: Option<&mut SequenceRecord>,
    rand: f64,
    analysis: PairAnalysis,
) -> CountupPrepassResult {
    let decision =
        decide_pair_from_analysis(prepass_config, r1, r2.as_deref(), analysis, Some(rand));
    let include = !decision.toss || add_bad_reads_countup;
    if prepass_config.error_correct && !decision.toss {
        let correction =
            correct_pair_errors_with_rollback(prepass_config, input_counts, r1, r2.as_deref_mut());
        if (!correction.uncorrectable || prepass_config.mark_uncorrectable_errors)
            && prepass_config.trim_after_marking
        {
            trim_pair(prepass_config, r1, r2);
        }
        return CountupPrepassResult {
            include,
            sort_analysis: None,
        };
    }
    CountupPrepassResult {
        include,
        sort_analysis: include.then_some(decision.analysis),
    }
}

fn compare_countup_work_pairs(left: &CountupWorkPair, right: &CountupWorkPair) -> CmpOrdering {
    left.input_list_index
        .cmp(&right.input_list_index)
        .then_with(|| compare_countup_sort_key(&left.sort_key, &right.sort_key))
        .then_with(|| left.r1.id.cmp(&right.r1.id))
        .then_with(|| {
            left.sort_key
                .original_index
                .cmp(&right.sort_key.original_index)
        })
}

fn compare_countup_sort_key(left: &CountupSortKey, right: &CountupSortKey) -> CmpOrdering {
    left.errors
        .cmp(&right.errors)
        .then_with(|| right.total_len.cmp(&left.total_len))
        .then_with(|| {
            left.expected_errors
                .partial_cmp(&right.expected_errors)
                .unwrap_or(CmpOrdering::Equal)
        })
        .then_with(|| left.numeric_id.cmp(&right.numeric_id))
}

fn countup_sort_key(
    config: &Config,
    input_counts: &dyn CountLookup,
    r1: &SequenceRecord,
    r2: Option<&SequenceRecord>,
    original_index: usize,
) -> CountupSortKey {
    let analysis = analyze_pair(config, input_counts, r1, r2);
    countup_sort_key_from_analysis(r1, r2, original_index, &analysis)
}

fn countup_sort_key_from_analysis(
    r1: &SequenceRecord,
    r2: Option<&SequenceRecord>,
    original_index: usize,
    analysis: &PairAnalysis,
) -> CountupSortKey {
    CountupSortKey {
        errors: analysis.low_kmer_count,
        total_len: r1.len() + r2.map(SequenceRecord::len).unwrap_or(0),
        expected_errors: expected_errors(r1) + r2.map(expected_errors).unwrap_or(0.0),
        numeric_id: r1.numeric_id,
        original_index,
    }
}

fn expected_errors(record: &SequenceRecord) -> f64 {
    let Some(qualities) = &record.qualities else {
        return 0.0;
    };
    record
        .bases
        .iter()
        .zip(qualities)
        .map(|(&base, &quality)| {
            let q = if is_defined_base(base) {
                quality.saturating_sub(33)
            } else {
                0
            };
            phred_error_probability(q)
        })
        .sum()
}

fn phred_error_probability(q: u8) -> f64 {
    match q {
        0 => 0.75,
        1 => 0.70,
        _ => 10f64.powf(-0.1 * f64::from(q)),
    }
}

fn unique_pair_kmers(
    config: &Config,
    r1: &SequenceRecord,
    r2: Option<&SequenceRecord>,
) -> Vec<KmerKey> {
    let mut keys = Vec::with_capacity(pair_kmer_window_capacity(config, r1, r2));
    fill_unique_pair_kmers(config, r1, r2, &mut keys);
    keys
}

fn fill_unique_pair_kmers(
    config: &Config,
    r1: &SequenceRecord,
    r2: Option<&SequenceRecord>,
    keys: &mut Vec<KmerKey>,
) {
    keys.clear();
    let required = pair_kmer_window_capacity(config, r1, r2);
    if keys.capacity() < required {
        keys.reserve(required - keys.capacity());
    }
    for_each_kmer_for_record(r1, config, |key| keys.push(key));
    if let Some(mate) = r2 {
        for_each_kmer_for_record(mate, config, |key| keys.push(key));
    }
    keys.sort_unstable();
    keys.dedup();
}

fn pair_kmer_window_capacity(
    config: &Config,
    r1: &SequenceRecord,
    r2: Option<&SequenceRecord>,
) -> usize {
    record_kmer_window_capacity(config.k, r1)
        .saturating_add(r2.map_or(0, |mate| record_kmer_window_capacity(config.k, mate)))
}

fn record_kmer_window_capacity(k: usize, record: &SequenceRecord) -> usize {
    if k == 0 {
        0
    } else {
        record.bases.len().saturating_sub(k).saturating_add(1)
    }
}

#[cfg(test)]
fn decide_countup_pair(
    config: &Config,
    input_counts: &dyn CountLookup,
    kept_counts: &dyn CountLookup,
    keys: &[KmerKey],
    target_depth: u64,
) -> bool {
    countup_decision_plan(config, input_counts, kept_counts, keys, target_depth).toss
}

fn countup_decision_plan(
    config: &Config,
    input_counts: &dyn CountLookup,
    kept_counts: &dyn CountLookup,
    keys: &[KmerKey],
    target_depth: u64,
) -> CountupDecisionPlan {
    let unique = keys.len();
    if unique == 0 {
        return CountupDecisionPlan {
            toss: !config.keep_all,
            eligible_key_indices: Vec::new(),
        };
    }

    let mut desired = 0usize;
    let mut needed = 0usize;
    let mut badly_needed = 0usize;
    let mut input_depths = config.toss_error_reads.then(Vec::new);
    let mut eligible_key_indices = Vec::with_capacity(keys.len());
    for (index, key) in keys.iter().enumerate() {
        let input_depth = input_counts.depth(key);
        if let Some(depths) = &mut input_depths {
            depths.push(input_depth);
        }
        if input_depth >= config.min_depth {
            desired += 1;
            eligible_key_indices.push(index);
            let kept_depth = kept_counts.depth(key);
            if kept_depth < target_depth {
                needed += 1;
                if kept_depth < target_depth.min(input_depth).saturating_mul(3) / 4 {
                    badly_needed += 1;
                }
            }
        }
    }

    let threshold_needed = 8usize.max(unique.div_ceil(6));
    let threshold_badly_needed = 2usize.max(unique.div_ceil(24));
    let keep = (needed >= threshold_needed || badly_needed >= threshold_badly_needed)
        && (desired >= config.min_kmers_over_min_depth || unique < config.min_kmers_over_min_depth);
    let mut toss = !keep;
    if config.toss_error_reads
        && let Some(mut depths) = input_depths
    {
        let errors = countup_error_count(&mut depths, config);
        if errors > 8 && needed < 2 * threshold_needed && badly_needed < 2 * threshold_badly_needed
        {
            toss = true;
        }
        if errors > unique / 2
            && needed < 3 * threshold_needed
            && badly_needed < 4 * threshold_badly_needed
        {
            toss = true;
        }
    }
    CountupDecisionPlan {
        toss: if config.keep_all { false } else { toss },
        eligible_key_indices,
    }
}

fn countup_error_count(depths: &mut [u64], config: &Config) -> usize {
    depths.sort_unstable();
    let mut previous: Option<u64> = None;
    for (index, &depth) in depths.iter().enumerate() {
        if let Some(prev) = previous
            && ((depth >= config.high_thresh && prev <= config.low_thresh)
                || depth >= prev.saturating_mul(config.error_detect_ratio))
        {
            return depths.len() - index;
        }
        previous = Some(depth);
    }
    0
}

#[cfg(test)]
fn increment_countup_kept_counts(
    config: &Config,
    kept_counts: &mut OutputCounts,
    input_counts: &dyn CountLookup,
    keys: &[KmerKey],
) {
    let mut atomic_increments = 0u64;
    for key in keys {
        if input_counts.depth(key) >= config.min_depth {
            match kept_counts {
                OutputCounts::Exact(counts) => {
                    *counts.entry(key.clone()).or_insert(0) += 1;
                }
                OutputCounts::Sketch(sketch) => sketch.increment(key),
                OutputCounts::AtomicSketch(sketch) => {
                    sketch.increment_key(key);
                    atomic_increments = atomic_increments.saturating_add(1);
                }
            }
        }
    }
    if let OutputCounts::AtomicSketch(sketch) = kept_counts {
        sketch.add_key_increments(atomic_increments);
    }
}

#[cfg(test)]
fn update_countup_kept_counts_for_decision(
    config: &Config,
    kept_counts: &mut OutputCounts,
    input_counts: &dyn CountLookup,
    keys: &[KmerKey],
    toss: bool,
) {
    if !toss || config.add_bad_reads_countup {
        increment_countup_kept_counts(config, kept_counts, input_counts, keys);
    }
}

fn update_countup_kept_counts_for_plan(
    config: &Config,
    kept_counts: &mut OutputCounts,
    keys: &[KmerKey],
    plan: &CountupDecisionPlan,
) {
    if plan.toss && !config.add_bad_reads_countup {
        return;
    }
    let mut atomic_increments = 0u64;
    for &index in &plan.eligible_key_indices {
        let Some(key) = keys.get(index) else {
            continue;
        };
        match kept_counts {
            OutputCounts::Exact(counts) => {
                *counts.entry(key.clone()).or_insert(0) += 1;
            }
            OutputCounts::Sketch(sketch) => sketch.increment(key),
            OutputCounts::AtomicSketch(sketch) => {
                sketch.increment_key(key);
                atomic_increments = atomic_increments.saturating_add(1);
            }
        }
    }
    if let OutputCounts::AtomicSketch(sketch) = kept_counts {
        sketch.add_key_increments(atomic_increments);
    }
}

fn countup_length_toss(config: &Config, r1: &SequenceRecord, r2: Option<&SequenceRecord>) -> bool {
    !config.keep_all
        && (r1.len() < config.min_length || r2.is_some_and(|mate| mate.len() < config.min_length))
}

#[cfg(test)]
fn count_map_depth_hist(counts: &CountMap, hist_len: usize) -> Vec<u64> {
    let mut hist = vec![0; hist_len];
    for &depth in counts.values() {
        let idx = (depth as usize).min(hist_len.saturating_sub(1));
        hist[idx] += depth;
    }
    hist
}

fn count_map_sparse_depth_hist(counts: &CountMap, hist_len: usize) -> SparseHist {
    let Some(last_index) = hist_len.checked_sub(1) else {
        return SparseHist::default();
    };
    let mut hist = SparseHist::default();
    for &depth in counts.values() {
        add_depth_to_sparse_hist(&mut hist, depth, last_index);
    }
    hist
}

#[cfg(test)]
fn add_depth_to_dynamic_hist(local: &mut Vec<u64>, depth: u64, last_index: usize) {
    if depth == 0 {
        return;
    }
    let idx = usize_from_u64_saturating(depth).min(last_index);
    if idx >= local.len() {
        local.resize(idx + 1, 0);
    }
    local[idx] = local[idx].saturating_add(depth);
}

#[cfg(test)]
fn merge_dynamic_depth_hist(mut left: Vec<u64>, right: Vec<u64>) -> Vec<u64> {
    if right.len() > left.len() {
        left.resize(right.len(), 0);
    }
    for (index, value) in right.into_iter().enumerate() {
        left[index] = left[index].saturating_add(value);
    }
    left
}

fn add_depth_to_sparse_hist(local: &mut SparseHist, depth: u64, last_index: usize) {
    if depth == 0 {
        return;
    }
    let idx = usize_from_u64_saturating(depth).min(last_index);
    let entry = local.entry(idx).or_insert(0);
    *entry = entry.saturating_add(depth);
}

fn merge_sparse_depth_hist(mut left: SparseHist, right: SparseHist) -> SparseHist {
    merge_sparse_hist(&mut left, right);
    left
}

#[cfg(test)]
fn build_input_counts(config: &Config) -> Result<InputCounts> {
    let mut stage_timings = Vec::new();
    build_input_counts_with_stage_timings(config, &mut stage_timings)
}

fn build_input_counts_with_stage_timings(
    config: &Config,
    stage_timings: &mut Vec<StageTiming>,
) -> Result<InputCounts> {
    let started = Instant::now();
    let counts = build_input_counts_inner(config, stage_timings)?;
    record_stage_timing(stage_timings, "input_counting", started);
    Ok(counts)
}

fn build_input_counts_inner(
    config: &Config,
    stage_timings: &mut Vec<StageTiming>,
) -> Result<InputCounts> {
    if use_bounded_input_sketch(config) {
        return build_sketch_input_counts(config, stage_timings);
    }
    let started = Instant::now();
    let mut counts = new_count_map(config);
    count_primary(config, &mut counts)?;
    for extra in &config.extra {
        count_single_file(config, extra, &mut counts, None)?;
    }
    apply_trusted_build_pass_filter(config, &mut counts);
    apply_prefilter_collision_estimates(config, &mut counts);
    apply_count_min_collision_estimates(config, &mut counts);
    record_stage_timing(stage_timings, "input_exact_counting", started);
    Ok(InputCounts::Exact(counts))
}

fn use_bounded_input_sketch(config: &Config) -> bool {
    if config.force_exact_counts {
        return false;
    }
    config.count_min.cells.is_some()
        || config.count_min.memory_bytes.is_some()
        || automatic_count_min_should_use(config)
}

fn gpu_counting_supported(config: &Config) -> bool {
    config.gpu_counting
        && config.gpu_helper.is_some()
        && config.k <= 31
        && !use_prefilter_collision_estimates(config)
}

fn build_sketch_input_counts(
    config: &Config,
    stage_timings: &mut Vec<StageTiming>,
) -> Result<InputCounts> {
    validate_gpu_counting_request(config)?;
    if use_prefilter_collision_estimates(config) {
        let started = Instant::now();
        let mut prefilter = new_input_prefilter_count_min_sketch(config)?;
        count_primary_prefilter_sketch(config, &mut prefilter)?;
        for extra in &config.extra {
            count_single_file_prefilter_sketch(config, extra, &mut prefilter, None)?;
        }
        let prefilter_limit = prefilter.max_count();
        record_stage_timing(stage_timings, "input_prefilter_counting", started);

        if use_atomic_count_min_sketch(config) {
            let started = Instant::now();
            let sketch = new_atomic_count_min_sketch_with_mask_seed(
                config,
                BBTOOLS_KCOUNT_ARRAY_SECOND_MASK_SEED,
            )?;
            count_primary_atomic_sketch(
                config,
                &sketch,
                Some(PrefilterGate::new(&prefilter, prefilter_limit)),
            )?;
            for extra in &config.extra {
                count_single_file_atomic_sketch(
                    config,
                    extra,
                    &sketch,
                    None,
                    Some(PrefilterGate::new(&prefilter, prefilter_limit)),
                )?;
            }
            record_stage_timing(stage_timings, "input_main_counting", started);
            return Ok(InputCounts::PrefilteredSketch {
                prefilter,
                limit: prefilter_limit,
                main: Box::new(InputCounts::AtomicSketch(sketch)),
            });
        }

        let started = Instant::now();
        let mut sketch = new_bounded_count_min_sketch_with_mask_seed(
            config,
            BBTOOLS_KCOUNT_ARRAY_SECOND_MASK_SEED,
        )?;
        count_primary_sketch(
            config,
            &mut sketch,
            Some(PrefilterGate::new(&prefilter, prefilter_limit)),
        )?;
        for extra in &config.extra {
            count_single_file_sketch(
                config,
                extra,
                &mut sketch,
                None,
                Some(PrefilterGate::new(&prefilter, prefilter_limit)),
            )?;
        }
        record_stage_timing(stage_timings, "input_main_counting", started);
        return Ok(InputCounts::PrefilteredSketch {
            prefilter,
            limit: prefilter_limit,
            main: Box::new(InputCounts::Sketch(sketch)),
        });
    }

    if use_atomic_count_min_sketch(config) {
        let started = Instant::now();
        let sketch = new_atomic_count_min_sketch(config)?;
        if gpu_counting_supported(config) {
            count_primary_gpu_reduced_runs_atomic_sketch(config, &sketch)?;
        } else {
            count_primary_atomic_sketch(config, &sketch, None)?;
        }
        for extra in &config.extra {
            count_single_file_atomic_sketch(config, extra, &sketch, None, None)?;
        }
        record_stage_timing(stage_timings, "input_main_counting", started);
        return Ok(InputCounts::AtomicSketch(sketch));
    }
    let started = Instant::now();
    let mut sketch = new_bounded_count_min_sketch(config)?;
    if gpu_counting_supported(config) {
        count_primary_gpu_reduced_runs_sketch(config, &mut sketch)?;
    } else {
        count_primary_sketch(config, &mut sketch, None)?;
    }
    for extra in &config.extra {
        count_single_file_sketch(config, extra, &mut sketch, None, None)?;
    }
    record_stage_timing(stage_timings, "input_main_counting", started);
    Ok(InputCounts::Sketch(sketch))
}

fn validate_gpu_counting_request(config: &Config) -> Result<()> {
    if !config.gpu_counting {
        return Ok(());
    }
    ensure!(
        config.gpu_helper.is_some(),
        "gpucounting=t requires gpuhelper=<cuda_kmer_reduce_runs binary>"
    );
    ensure!(
        config.k <= 31,
        "gpucounting=t currently supports short k-mers only (k<=31)"
    );
    ensure!(
        !use_prefilter_collision_estimates(config),
        "gpucounting=t currently supports the main bounded sketch without prefilter=t"
    );
    Ok(())
}

fn new_output_counts(config: &Config) -> Result<OutputCounts> {
    if use_bounded_input_sketch(config) {
        if config.count_up {
            return new_countup_output_counts(config);
        }
        if use_atomic_count_min_sketch(config) {
            new_atomic_output_count_min_sketch(config).map(OutputCounts::AtomicSketch)
        } else {
            new_bounded_output_count_min_sketch(config).map(OutputCounts::Sketch)
        }
    } else {
        Ok(OutputCounts::Exact(new_count_map(config)))
    }
}

fn new_atomic_output_count_min_sketch(config: &Config) -> Result<AtomicCountMinSketch> {
    let hashes = config
        .count_min
        .hashes
        .unwrap_or(3)
        .clamp(1, BBTOOLS_KCOUNT_ARRAY_MAX_HASHES);
    let total_cells = output_count_min_total_cells(config, 32);
    ensure_count_min_budget_fits_memory(
        "output_kept",
        total_cells,
        32,
        output_count_min_memory_bytes(config, 32),
    )?;
    let min_arrays = kcount_array_min_arrays(config);
    let cells = count_min_table_cells_from_total_bits_with_min_arrays(total_cells, 32, min_arrays);
    let update_mode = count_min_update_mode(config, 32, hashes);
    AtomicCountMinSketch::new_with_min_arrays_and_update_mode(
        cells,
        hashes,
        min_arrays,
        update_mode,
        kept_output_mask_seed(config),
    )
    .map(|sketch| sketch.with_parallel_replay(!config.deterministic))
}

fn new_bounded_output_count_min_sketch(config: &Config) -> Result<PackedCountMinSketch> {
    let hashes = config
        .count_min
        .hashes
        .unwrap_or(3)
        .clamp(1, BBTOOLS_KCOUNT_ARRAY_MAX_HASHES);
    let bits = config.count_min.bits.unwrap_or(32);
    let total_cells = output_count_min_total_cells(config, bits);
    ensure_count_min_budget_fits_memory(
        "output_kept",
        total_cells,
        bits,
        output_count_min_memory_bytes(config, bits),
    )?;
    let min_arrays = kcount_array_min_arrays(config);
    let cells =
        count_min_table_cells_from_total_bits_with_min_arrays(total_cells, bits, min_arrays);
    PackedCountMinSketch::new_with_min_arrays_and_mask_seed(
        cells,
        hashes,
        bits,
        min_arrays,
        kept_output_mask_seed(config),
    )
    .map(|sketch| sketch.with_update_mode(count_min_update_mode(config, bits, hashes)))
}

fn new_countup_output_counts(config: &Config) -> Result<OutputCounts> {
    let bits = countup_output_count_bits(config);
    let hashes = 3;
    let total_cells = countup_output_total_cells(config, bits);
    ensure_count_min_budget_fits_memory(
        "count-up output",
        total_cells,
        bits,
        config
            .count_min
            .memory_bytes
            .or_else(|| automatic_count_min_memory_bytes(config)),
    )?;
    let min_arrays = kcount_array_min_arrays(config);
    let cells =
        count_min_table_cells_from_total_bits_with_min_arrays(total_cells, bits, min_arrays);
    PackedCountMinSketch::new_with_min_arrays_and_mask_seed(
        cells,
        hashes,
        bits,
        min_arrays,
        countup_output_mask_seed(config),
    )
    .map(|sketch| sketch.with_update_mode(count_min_update_mode(config, bits, hashes)))
    .map(OutputCounts::Sketch)
}

fn countup_output_count_bits(config: &Config) -> u8 {
    let target = countup_adjusted_target_depth(config);
    if target <= 15 {
        4
    } else if target <= 255 {
        8
    } else {
        16
    }
}

fn countup_adjusted_target_depth(config: &Config) -> u64 {
    ((config.target_depth as f64) * 0.95).round().max(1.0) as u64
}

fn countup_output_total_cells(config: &Config, bits: u8) -> usize {
    config
        .count_min
        .cells
        .unwrap_or_else(|| count_min_cells_from_memory(count_min_memory_bytes(config), bits))
        .max(1)
}

fn countup_output_mask_seed(config: &Config) -> u64 {
    kept_output_mask_seed(config)
}

fn kept_output_mask_seed(config: &Config) -> u64 {
    let preceding_tables = if use_prefilter_collision_estimates(config) {
        2
    } else {
        1
    };
    BBTOOLS_KCOUNT_ARRAY_FIRST_MASK_SEED
        .saturating_add(BBTOOLS_KCOUNT_ARRAY_MASK_SEED_STEP.saturating_mul(preceding_tables))
}

fn output_count_min_total_cells(config: &Config, bits: u8) -> usize {
    let base = config
        .count_min
        .cells
        .unwrap_or_else(|| {
            count_min_cells_from_memory(output_count_min_memory_bytes(config, bits), bits)
        })
        .max(1);
    let Some(fraction_micros) = prefilter_memory_fraction_micros(config) else {
        return cap_main_cells_to_short_kmer_space(config, base);
    };
    let main_fraction = 1_000_000usize.saturating_sub(fraction_micros as usize);
    base.saturating_mul(main_fraction)
        .checked_div(1_000_000)
        .unwrap_or(0)
        .max(1)
}

fn output_count_min_memory_bytes(config: &Config, _bits: u8) -> Option<usize> {
    if config.count_min.cells.is_some() {
        return config
            .count_min
            .memory_bytes
            .or(config.auto_count_min_memory_bytes)
            .or_else(|| automatic_count_min_memory_bytes(config));
    }
    if config.count_min.memory_bytes.is_some() {
        return count_min_memory_bytes(config);
    }
    automatic_count_min_memory_bytes(config).map(output_count_min_auto_memory_bytes)
}

fn output_count_min_auto_memory_bytes(memory_bytes: usize) -> usize {
    let min_memory = OUTPUT_COUNT_MIN_AUTO_MIN_MEMORY_BYTES.min(memory_bytes);
    scale_by_micros(memory_bytes, OUTPUT_COUNT_MIN_AUTO_FRACTION_MICROS)
        .max(min_memory)
        .min(memory_bytes)
        .max(1)
}

fn use_atomic_count_min_sketch(config: &Config) -> bool {
    config.count_min.bits.unwrap_or(32) == 32
}

fn new_atomic_count_min_sketch(config: &Config) -> Result<AtomicCountMinSketch> {
    new_atomic_count_min_sketch_with_mask_seed(config, BBTOOLS_KCOUNT_ARRAY_FIRST_MASK_SEED)
}

fn new_atomic_count_min_sketch_with_mask_seed(
    config: &Config,
    mask_seed: u64,
) -> Result<AtomicCountMinSketch> {
    let hashes = config
        .count_min
        .hashes
        .unwrap_or(3)
        .clamp(1, BBTOOLS_KCOUNT_ARRAY_MAX_HASHES);
    let total_cells = main_count_min_total_cells(config, 32);
    ensure_count_min_budget_fits_memory(
        "main",
        total_cells,
        32,
        config
            .count_min
            .memory_bytes
            .or(config.auto_count_min_memory_bytes),
    )?;
    let min_arrays = kcount_array_min_arrays(config);
    let cells = count_min_table_cells_from_total_bits_with_min_arrays(total_cells, 32, min_arrays);
    let update_mode = count_min_update_mode(config, 32, hashes);
    AtomicCountMinSketch::new_with_min_arrays_and_update_mode(
        cells,
        hashes,
        min_arrays,
        update_mode,
        mask_seed,
    )
    .map(|sketch| sketch.with_parallel_replay(!config.deterministic))
}

fn new_bounded_count_min_sketch(config: &Config) -> Result<PackedCountMinSketch> {
    new_bounded_count_min_sketch_with_mask_seed(config, BBTOOLS_KCOUNT_ARRAY_FIRST_MASK_SEED)
}

fn new_bounded_count_min_sketch_with_mask_seed(
    config: &Config,
    mask_seed: u64,
) -> Result<PackedCountMinSketch> {
    let hashes = config
        .count_min
        .hashes
        .unwrap_or(3)
        .clamp(1, BBTOOLS_KCOUNT_ARRAY_MAX_HASHES);
    let bits = config.count_min.bits.unwrap_or(32);
    let total_cells = main_count_min_total_cells(config, bits);
    ensure_count_min_budget_fits_memory(
        "main",
        total_cells,
        bits,
        config
            .count_min
            .memory_bytes
            .or(config.auto_count_min_memory_bytes),
    )?;
    let min_arrays = kcount_array_min_arrays(config);
    let cells =
        count_min_table_cells_from_total_bits_with_min_arrays(total_cells, bits, min_arrays);
    PackedCountMinSketch::new_with_min_arrays_and_mask_seed(
        cells, hashes, bits, min_arrays, mask_seed,
    )
    .map(|sketch| sketch.with_update_mode(count_min_update_mode(config, bits, hashes)))
}

fn new_prefilter_count_min_sketch(config: &Config) -> Result<PackedCountMinSketch> {
    let hashes = config
        .prefilter
        .hashes
        .unwrap_or_else(|| default_prefilter_hashes(config))
        .clamp(1, BBTOOLS_KCOUNT_ARRAY_MAX_HASHES);
    let bits = config.prefilter.bits.unwrap_or(DEFAULT_PREFILTER_BITS);
    let total_cells = prefilter_total_cells(config, bits).max(1);
    ensure_count_min_budget_fits_memory(
        "prefilter",
        total_cells,
        bits,
        config
            .prefilter
            .memory_bytes
            .or(config.count_min.memory_bytes)
            .or(config.auto_count_min_memory_bytes),
    )?;
    let min_arrays = kcount_array_min_arrays(config);
    let cells =
        count_min_table_cells_from_total_bits_with_min_arrays(total_cells, bits, min_arrays);
    PackedCountMinSketch::new_with_min_arrays(cells, hashes, bits, min_arrays)
        .map(|sketch| sketch.with_update_mode(count_min_update_mode(config, bits, hashes)))
}

fn new_input_prefilter_count_min_sketch(config: &Config) -> Result<PrefilterCountMinSketch> {
    if config.deterministic {
        return new_prefilter_count_min_sketch(config).map(PrefilterCountMinSketch::Packed);
    }
    new_atomic_packed_prefilter_count_min_sketch(config).map(PrefilterCountMinSketch::AtomicPacked)
}

fn new_atomic_packed_prefilter_count_min_sketch(
    config: &Config,
) -> Result<AtomicPackedCountMinSketch> {
    let hashes = config
        .prefilter
        .hashes
        .unwrap_or_else(|| default_prefilter_hashes(config))
        .clamp(1, BBTOOLS_KCOUNT_ARRAY_MAX_HASHES);
    let bits = config.prefilter.bits.unwrap_or(DEFAULT_PREFILTER_BITS);
    let total_cells = prefilter_total_cells(config, bits).max(1);
    ensure_count_min_budget_fits_memory(
        "prefilter",
        total_cells,
        bits,
        config
            .prefilter
            .memory_bytes
            .or(config.count_min.memory_bytes)
            .or(config.auto_count_min_memory_bytes),
    )?;
    let min_arrays = kcount_array_min_arrays(config);
    let cells =
        count_min_table_cells_from_total_bits_with_min_arrays(total_cells, bits, min_arrays);
    AtomicPackedCountMinSketch::new_with_min_arrays_and_update_mode(
        cells,
        hashes,
        bits,
        min_arrays,
        count_min_update_mode(config, bits, hashes),
        BBTOOLS_KCOUNT_ARRAY_FIRST_MASK_SEED,
    )
}

fn default_prefilter_hashes(config: &Config) -> usize {
    let main_hashes = config
        .count_min
        .hashes
        .unwrap_or(3)
        .clamp(1, BBTOOLS_KCOUNT_ARRAY_MAX_HASHES);
    main_hashes.div_ceil(2)
}

fn count_min_update_mode(config: &Config, bits: u8, hashes: usize) -> CountMinUpdateMode {
    // KCountArray7MTA defaults to locked/symmetric writes for multi-hash,
    // multi-bit cells; lockedincrement=f / symmetricwrite=f opts into the
    // faster independent-row update path.
    if bits > 1 && hashes > 1 && config.locked_increment.unwrap_or(true) {
        CountMinUpdateMode::Conservative
    } else {
        CountMinUpdateMode::Independent
    }
}

fn count_min_memory_bytes(config: &Config) -> Option<usize> {
    config
        .count_min
        .memory_bytes
        .or_else(|| automatic_count_min_memory_bytes(config))
}

fn main_count_min_total_cells(config: &Config, bits: u8) -> usize {
    let base = config
        .count_min
        .cells
        .unwrap_or_else(|| count_min_cells_from_memory(count_min_memory_bytes(config), bits))
        .max(1);
    let Some(fraction_micros) = prefilter_memory_fraction_micros(config) else {
        return cap_main_cells_to_short_kmer_space(config, base);
    };
    let main_fraction = 1_000_000usize.saturating_sub(fraction_micros as usize);
    base.saturating_mul(main_fraction)
        .checked_div(1_000_000)
        .unwrap_or(0)
        .max(1)
}

fn cap_main_cells_to_short_kmer_space(config: &Config, cells: usize) -> usize {
    if use_prefilter_collision_estimates(config) {
        return cells;
    }
    short_kmer_space_cells(config.k)
        .map(|cap| cells.min(cap))
        .unwrap_or(cells)
        .max(1)
}

fn short_kmer_space_cells(k: usize) -> Option<usize> {
    if k >= 32 {
        return None;
    }
    1usize.checked_shl((2 * k) as u32)
}

fn prefilter_memory_fraction_micros(config: &Config) -> Option<u32> {
    if config.prefilter.force_disabled {
        return None;
    }
    if config.prefilter.cells.is_some() || config.prefilter.memory_bytes.is_some() {
        return None;
    }
    if let Some(fraction) = config
        .prefilter
        .memory_fraction_micros
        .filter(|fraction| *fraction > 0)
    {
        return Some(fraction);
    }
    if config.prefilter.enabled && use_bounded_input_sketch(config) {
        return Some(DEFAULT_PREFILTER_FRACTION_MICROS);
    }
    None
}

fn scale_by_micros(value: usize, micros: u32) -> usize {
    value
        .saturating_mul(micros as usize)
        .checked_div(1_000_000)
        .unwrap_or(0)
}

fn zeroed_u64_vec(len: usize) -> Result<Vec<u64>> {
    // SAFETY: an all-zero bit pattern is a valid initialized `u64`.
    unsafe { zeroed_vec_with_layout::<u64>(len, "u64") }
}

fn zeroed_atomic_u32_vec(len: usize) -> Result<Vec<AtomicU32>> {
    // SAFETY: an all-zero bit pattern is a valid initialized `AtomicU32`.
    unsafe { zeroed_vec_with_layout::<AtomicU32>(len, "AtomicU32") }
}

fn zeroed_atomic_u64_vec(len: usize) -> Result<Vec<AtomicU64>> {
    // SAFETY: an all-zero bit pattern is a valid initialized `AtomicU64`.
    unsafe { zeroed_vec_with_layout::<AtomicU64>(len, "AtomicU64") }
}

unsafe fn zeroed_vec_with_layout<T>(len: usize, type_name: &str) -> Result<Vec<T>> {
    if len == 0 {
        return Ok(Vec::new());
    }
    let layout = Layout::array::<T>(len)
        .with_context(|| format!("allocating zeroed {type_name} vector layout"))?;
    // SAFETY: The caller guarantees an all-zero bit pattern is valid for `T`.
    // The requested layout is for exactly `len` values of `T`, and the returned
    // pointer is converted back into a Vec with the same length and capacity.
    let ptr = unsafe { alloc_zeroed(layout) };
    if ptr.is_null() {
        bail!("allocating zeroed {type_name} vector failed for {len} elements");
    }
    // SAFETY: `ptr` came from the global allocator with `Layout::array::<T>(len)`,
    // is non-null, and points to zero-initialized memory valid for `T` by the
    // caller contract. Vec will deallocate it with the matching layout.
    Ok(unsafe { Vec::from_raw_parts(ptr.cast::<T>(), len, len) })
}

fn count_min_cells_from_memory(memory_bytes: Option<usize>, bits: u8) -> usize {
    let Some(memory_bytes) = memory_bytes else {
        return DEFAULT_PREFILTER_CELLS;
    };
    let bits_total = memory_bytes.saturating_mul(8);
    let bits_per_cell = bits.max(1) as usize;
    (bits_total / bits_per_cell).max(1)
}

fn count_min_total_bytes(total_cells: usize, bits: u8) -> Result<usize> {
    let total_cells = total_cells.max(1);
    let bits = bits.max(1) as usize;
    let total_bits = total_cells
        .checked_mul(bits)
        .context("bounded count-min sketch size overflowed")?;
    Ok(total_bits.div_ceil(8).max(1))
}

fn packed_sketch_should_track_slots(cells: usize) -> bool {
    cells <= PACKED_SKETCH_TRACKED_SLOT_LIMIT
}

fn safe_explicit_count_min_bytes(available: usize) -> usize {
    available
        .saturating_mul(EXPLICIT_COUNT_MIN_SAFE_MEMORY_PERCENT)
        .checked_div(100)
        .unwrap_or(0)
        .max(1)
}

fn count_min_safe_budget_bytes(
    configured_memory_bytes: Option<usize>,
    available_memory_bytes: Option<usize>,
) -> Option<usize> {
    let safe_available = available_memory_bytes.map(safe_explicit_count_min_bytes);
    match (configured_memory_bytes, safe_available) {
        (Some(configured), Some(available)) => Some(configured.min(available)),
        (Some(configured), None) => Some(configured),
        (None, Some(available)) => Some(available),
        (None, None) => None,
    }
}

fn ensure_count_min_budget_fits_ceiling(
    label: &str,
    total_cells: usize,
    bits: u8,
    safe_budget: usize,
) -> Result<()> {
    let requested = count_min_total_bytes(total_cells, bits)?;
    if requested > safe_budget {
        bail!(
            "{label} count-min table requests {requested} bytes ({total_cells} cells x {} bits), above safe memory budget {safe_budget} bytes; reduce cells/matrixbits/sketchmemory/mem",
            bits.max(1)
        );
    }
    Ok(())
}

fn ensure_count_min_budget_fits_memory(
    label: &str,
    total_cells: usize,
    bits: u8,
    configured_memory_bytes: Option<usize>,
) -> Result<()> {
    if let Some(safe_budget) =
        count_min_safe_budget_bytes(configured_memory_bytes, system_available_memory_bytes())
    {
        ensure_count_min_budget_fits_ceiling(label, total_cells, bits, safe_budget)
    } else {
        count_min_total_bytes(total_cells, bits).map(|_| ())
    }
}

#[cfg(test)]
fn count_min_table_cells_from_total(total_cells: usize, hashes: usize) -> usize {
    let _ = hashes;
    count_min_table_cells_from_total_bits(total_cells, 32)
}

#[cfg(test)]
fn count_min_table_cells_from_total_bits(total_cells: usize, bits: u8) -> usize {
    count_min_table_cells_from_total_bits_with_min_arrays(
        total_cells,
        bits,
        BBTOOLS_KCOUNT_ARRAY_MIN_ARRAYS,
    )
}

fn count_min_table_cells_from_total_bits_with_min_arrays(
    total_cells: usize,
    bits: u8,
    min_arrays: usize,
) -> usize {
    let total_cells = total_cells.max(1);
    let arrays = kcount_array_count(total_cells, bits, min_arrays);
    if arrays <= 1 {
        return prime_at_most(total_cells);
    }
    prime_at_most(total_cells.div_ceil(arrays)).saturating_mul(arrays)
}

fn kcount_array_min_arrays(config: &Config) -> usize {
    kcount_array_min_arrays_for_threads(config.threads.unwrap_or_else(rayon::current_num_threads))
}

fn kcount_array_min_arrays_for_threads(threads: usize) -> usize {
    let target = threads.max(BBTOOLS_KCOUNT_ARRAY_MIN_ARRAYS);
    let mut arrays = BBTOOLS_KCOUNT_ARRAY_MIN_ARRAYS;
    while arrays < target {
        let next = arrays.saturating_mul(2);
        if next == arrays {
            break;
        }
        arrays = next;
    }
    arrays
}

fn kcount_array_lock_index(key: &KmerKey) -> usize {
    let raw = match key {
        KmerKey::Short(raw) | KmerKey::LongHash(raw) => *raw,
    };
    ((raw & (i64::MAX as u64)) % BBTOOLS_KCOUNT_ARRAY_LOCKS as u64) as usize
}

fn kcount_array_count(desired_cells: usize, bits: u8, min_arrays: usize) -> usize {
    if desired_cells < BBTOOLS_KCOUNT_ARRAY_SHARD_MIN_CELLS {
        return 1;
    }
    let bits = bits.clamp(1, 64) as usize;
    let min_arrays = kcount_array_min_arrays_for_threads(min_arrays);
    let words = desired_cells
        .saturating_mul(bits)
        .saturating_add(31)
        .checked_div(32)
        .unwrap_or(usize::MAX)
        .max(min_arrays);
    let mut arrays = min_arrays;
    while words / arrays >= i32::MAX as usize {
        arrays = arrays.saturating_mul(2);
    }
    while arrays > desired_cells {
        arrays /= 2;
    }
    arrays.max(1)
}

fn prime_at_most(value: usize) -> usize {
    if value <= 2 {
        return value.max(1);
    }

    let mut candidate = if value.is_multiple_of(2) {
        value - 1
    } else {
        value
    };
    while candidate > 2 {
        if is_prime(candidate) {
            return candidate;
        }
        candidate -= 2;
    }
    2
}

fn is_prime(value: usize) -> bool {
    if value <= 3 {
        return value > 1;
    }
    if value.is_multiple_of(2) || value.is_multiple_of(3) {
        return false;
    }

    let mut divisor = 5usize;
    while divisor <= value / divisor {
        if value.is_multiple_of(divisor) || value.is_multiple_of(divisor + 2) {
            return false;
        }
        divisor += 6;
    }
    true
}

fn automatic_count_min_should_use(config: &Config) -> bool {
    if !config.auto_count_min || config.force_exact_counts {
        return false;
    }
    if config
        .table_reads
        .or(config.max_reads)
        .is_some_and(|reads| reads >= config.auto_count_min_read_threshold)
    {
        return true;
    }
    input_metadata_bytes(config)
        .is_some_and(|bytes| bytes >= config.auto_count_min_input_bytes as u64)
}

fn automatic_count_min_memory_bytes(config: &Config) -> Option<usize> {
    if !automatic_count_min_should_use(config) {
        return None;
    }
    let raw_memory = config
        .auto_count_min_memory_bytes
        .unwrap_or_else(default_auto_count_min_memory_bytes);
    Some(automatic_count_min_filter_memory_bytes(config, raw_memory))
}

fn automatic_count_min_filter_memory_bytes(config: &Config, raw_memory: usize) -> usize {
    let usable = bbtools_usable_table_memory_bytes(config, raw_memory).max(1);
    if config.count_up {
        (usable / 2).max(1)
    } else {
        usable
    }
}

fn default_auto_count_min_memory_bytes() -> usize {
    system_available_memory_bytes()
        .map(|bytes| {
            (bytes / 4).clamp(
                AUTO_COUNT_MIN_MIN_MEMORY_BYTES,
                AUTO_COUNT_MIN_MAX_MEMORY_BYTES,
            )
        })
        .unwrap_or(AUTO_COUNT_MIN_FALLBACK_MEMORY_BYTES)
}

fn bbtools_usable_table_memory_bytes(config: &Config, memory_bytes: usize) -> usize {
    let after_headroom = memory_bytes.saturating_sub(BBTOOLS_MEMORY_HEADROOM_BYTES) as f64 * 0.73;
    let fraction = memory_bytes as f64 * 0.45;
    let mut usable = after_headroom.max(fraction).max(1.0) as usize;
    if histogram_memory_is_reserved(config) {
        let threads = config
            .threads
            .unwrap_or_else(rayon::current_num_threads)
            .max(1);
        let hist_bytes = config
            .hist_len
            .saturating_mul(8)
            .saturating_mul(threads.saturating_add(1));
        usable = usable.saturating_sub(hist_bytes);
    }
    if config.build_passes > 1 {
        usable /= 2;
    }
    usable.max(1)
}

fn histogram_memory_is_reserved(config: &Config) -> bool {
    config.hist_in.is_some()
        || config.hist_out.is_some()
        || config.peaks_in.is_some()
        || config.peaks_out.is_some()
}

fn system_available_memory_bytes() -> Option<usize> {
    let text = fs::read_to_string("/proc/meminfo").ok()?;
    for line in text.lines() {
        if let Some(rest) = line.strip_prefix("MemAvailable:") {
            let kb = rest.split_whitespace().next()?.parse::<usize>().ok()?;
            return kb.checked_mul(1024);
        }
    }
    None
}

fn input_metadata_bytes(config: &Config) -> Option<u64> {
    let mut total = 0u64;
    let mut found = false;
    for path in input_metadata_paths(config) {
        let Ok(metadata) = fs::metadata(path) else {
            continue;
        };
        if metadata.is_file() {
            total = total.saturating_add(metadata.len());
            found = true;
        }
    }
    found.then_some(total)
}

fn input_metadata_paths(config: &Config) -> Vec<PathBuf> {
    let mut paths = Vec::new();
    if let Some(path) = &config.in1 {
        paths.extend(metadata_path_expansion(path));
    }
    if let Some(path) = &config.in2 {
        paths.extend(metadata_path_expansion(path));
    }
    for path in &config.extra {
        paths.extend(metadata_path_expansion(path));
    }
    paths
}

fn metadata_path_expansion(path: &Path) -> Vec<PathBuf> {
    if path.exists() {
        return vec![path.to_path_buf()];
    }
    let text = path.to_string_lossy();
    if text.contains(',') {
        split_path_list(&text)
    } else {
        vec![path.to_path_buf()]
    }
}

fn apply_output_count_adjustments(config: &Config, counts: &mut OutputCounts) {
    let OutputCounts::Exact(counts) = counts else {
        return;
    };
    apply_trusted_build_pass_filter(config, counts);
    apply_prefilter_collision_estimates(config, counts);
    apply_count_min_collision_estimates(config, counts);
}

fn apply_trusted_build_pass_filter(config: &Config, counts: &mut CountMap) {
    if config.build_passes <= 1 || counts.len() < 2 {
        return;
    }
    let decrement = (config.build_passes as u64).saturating_sub(1);
    for count in counts.values_mut() {
        if *count > 1 {
            *count = count.saturating_sub(decrement).max(1);
        }
    }
}

fn apply_prefilter_collision_estimates(config: &Config, counts: &mut CountMap) {
    if config.force_exact_counts {
        return;
    }
    if !use_prefilter_collision_estimates(config) {
        return;
    };
    if counts.len() < 2 {
        return;
    }
    let entries = sorted_count_entries(counts);
    let Ok(mut sketch) = new_prefilter_count_min_sketch(config) else {
        return;
    };
    sketch.add_key_counts(counts);

    for (key, exact) in entries {
        let estimate = sketch.depth(&key);
        if estimate < sketch.max_count {
            counts.insert(key, estimate);
        } else {
            counts.insert(key, exact);
        }
    }
}

fn use_prefilter_collision_estimates(config: &Config) -> bool {
    if config.prefilter.force_disabled {
        return false;
    }
    config.prefilter.cells.is_some()
        || config.prefilter.hashes.is_some()
        || config.prefilter.memory_bytes.is_some()
        || config
            .prefilter
            .memory_fraction_micros
            .is_some_and(|fraction| fraction > 0)
        || (config.prefilter.enabled && use_bounded_input_sketch(config))
}

fn prefilter_total_cells(config: &Config, bits: u8) -> usize {
    if let Some(cells) = config.prefilter.cells {
        return cells.max(1);
    }
    if let Some(memory_bytes) = config.prefilter.memory_bytes {
        return count_min_cells_from_memory(Some(memory_bytes), bits);
    }
    if let Some(fraction_micros) = prefilter_memory_fraction_micros(config) {
        if let Some(total_cells) = config.count_min.cells {
            let main_bits = config.count_min.bits.unwrap_or(32).max(1) as usize;
            let prefilter_bits = scale_by_micros(
                total_cells.max(1).saturating_mul(main_bits),
                fraction_micros,
            )
            .max(bits.max(1) as usize);
            return (prefilter_bits / bits.max(1) as usize).max(1);
        }
        if let Some(memory_bytes) =
            count_min_memory_bytes(config).or(config.auto_count_min_memory_bytes)
        {
            let prefilter_memory = scale_by_micros(memory_bytes, fraction_micros).max(1);
            return count_min_cells_from_memory(Some(prefilter_memory), bits);
        }
    }
    DEFAULT_PREFILTER_CELLS
}

fn apply_count_min_collision_estimates(config: &Config, counts: &mut CountMap) {
    if config.force_exact_counts {
        return;
    }
    let Some(cells) = config.count_min.cells else {
        return;
    };
    if cells == 0 || counts.len() < 2 {
        return;
    }
    let entries = sorted_count_entries(counts);
    let Ok(mut sketch) = new_bounded_count_min_sketch(config) else {
        return;
    };
    sketch.add_key_counts(counts);

    for (key, exact) in entries {
        let exact = exact.min(sketch.max_count);
        let estimate = sketch.depth(&key).max(exact).min(sketch.max_count);
        counts.insert(key, estimate);
    }
}

fn sorted_count_entries(counts: &CountMap) -> Vec<(KmerKey, u64)> {
    let mut entries: Vec<_> = counts
        .iter()
        .map(|(key, &count)| (key.clone(), count))
        .collect();
    entries.sort_unstable_by(|(left, _), (right, _)| left.cmp(right));
    entries
}

impl PackedCountMinSketch {
    #[cfg(test)]
    fn new(cells: usize, hashes: usize, bits: u8) -> Result<Self> {
        Self::new_with_min_arrays(cells, hashes, bits, BBTOOLS_KCOUNT_ARRAY_MIN_ARRAYS)
    }

    fn new_with_min_arrays(
        cells: usize,
        hashes: usize,
        bits: u8,
        min_arrays: usize,
    ) -> Result<Self> {
        Self::new_with_min_arrays_and_mask_seed(
            cells,
            hashes,
            bits,
            min_arrays,
            BBTOOLS_KCOUNT_ARRAY_FIRST_MASK_SEED,
        )
    }

    fn new_with_min_arrays_and_mask_seed(
        cells: usize,
        hashes: usize,
        bits: u8,
        min_arrays: usize,
        mask_seed: u64,
    ) -> Result<Self> {
        let cells = cells.max(1);
        let hashes = hashes.clamp(1, BBTOOLS_KCOUNT_ARRAY_MAX_HASHES);
        let bits = bits.clamp(1, 64);
        let layout = KCountArrayLayout::new_with_min_arrays_and_mask_seed(
            cells, bits, min_arrays, mask_seed,
        );
        let word_count = if bits == 64 {
            cells
        } else {
            let total_bits = cells
                .checked_mul(bits as usize)
                .context("bounded sketch bit count overflowed")?;
            total_bits.div_ceil(64)
        };
        let words = zeroed_u64_vec(word_count).context("allocating bounded count-min sketch")?;
        Ok(Self {
            cells,
            hashes,
            bits,
            max_count: count_min_max_count(bits),
            layout,
            update_mode: CountMinUpdateMode::Conservative,
            words,
            increments: 0,
            occupied_slots: 0,
            tracked_slots: packed_sketch_should_track_slots(cells).then(Vec::new),
        })
    }

    fn with_update_mode(mut self, update_mode: CountMinUpdateMode) -> Self {
        self.update_mode = update_mode;
        self
    }

    fn layout_summary(
        &self,
        table: &'static str,
        prefilter_limit: Option<u64>,
    ) -> SketchLayoutSummary {
        SketchLayoutSummary {
            table,
            kind: "packed",
            cells: self.cells,
            hashes: self.hashes,
            bits: self.bits,
            arrays: self.layout.array_count(),
            cells_per_array: self.layout.cells_per_array,
            mask_seed: self.layout.mask_seed,
            update_mode: self.update_mode.as_str(),
            max_count: self.max_count,
            memory_bytes: self.estimated_memory_bytes(),
            prefilter_limit,
        }
    }

    fn estimated_memory_bytes(&self) -> usize {
        self.words
            .len()
            .saturating_mul(std::mem::size_of::<u64>())
            .saturating_add(self.tracked_slot_memory_bytes())
    }

    fn tracked_slot_memory_bytes(&self) -> usize {
        self.tracked_slots.as_ref().map_or(0, |slots| {
            slots
                .capacity()
                .saturating_mul(std::mem::size_of::<usize>())
        })
    }

    fn increment(&mut self, key: &KmerKey) {
        self.add_key_count(key, 1);
        self.increments = self.increments.saturating_add(1);
    }

    fn add_key_count(&mut self, key: &KmerKey, count: u64) {
        let _ = self.increment_and_return_unincremented(key, count);
    }

    fn increment_and_return_unincremented(&mut self, key: &KmerKey, count: u64) -> u64 {
        if count == 0 {
            return self.depth(key);
        }
        if self.update_mode == CountMinUpdateMode::Independent {
            return self.increment_independent_and_return_unincremented(key, count);
        }
        if self.bits == 2 && self.hashes == 2 {
            return self.increment_2bit_2hash_conservative_and_return_unincremented(key, count);
        }
        if self.bits == 16 && self.hashes == 3 {
            return self.increment_16bit_3hash_conservative_and_return_unincremented(key, count);
        }
        let target_increment = count.min(self.max_count);
        let mut slots = [0usize; 16];
        let mut min_depth = self.max_count;
        fill_count_min_buckets(key, self.hashes, self.layout, &mut slots);
        for slot in slots.iter().take(self.hashes) {
            min_depth = min_depth.min(self.cell(*slot));
        }
        if min_depth >= self.max_count {
            return min_depth;
        }
        let target = min_depth
            .saturating_add(target_increment)
            .min(self.max_count);
        let mut previous_min = self.max_count;
        for slot in slots.iter().take(self.hashes) {
            let previous = self.cell(*slot);
            previous_min = previous_min.min(previous);
            if previous < target {
                self.set_cell_with_previous(*slot, previous, target);
            }
        }
        previous_min
    }

    fn increment_16bit_3hash_conservative_and_return_unincremented(
        &mut self,
        key: &KmerKey,
        count: u64,
    ) -> u64 {
        let [first, second, third] = count_min_three_buckets_raw(raw_kmer_key(key), self.layout);
        let first_depth = self.cell_16bit(first);
        let second_depth = self.cell_16bit(second);
        let third_depth = self.cell_16bit(third);
        let min_depth = first_depth.min(second_depth).min(third_depth);
        if min_depth >= self.max_count {
            return min_depth;
        }
        let target = min_depth
            .saturating_add(count.min(self.max_count))
            .min(self.max_count);
        if first_depth < target {
            self.set_cell_16bit_with_previous(first, first_depth, target);
        }
        if second_depth < target {
            self.set_cell_16bit_with_previous(second, second_depth, target);
        }
        if third_depth < target {
            self.set_cell_16bit_with_previous(third, third_depth, target);
        }
        min_depth
    }

    fn increment_2bit_2hash_conservative_and_return_unincremented(
        &mut self,
        key: &KmerKey,
        count: u64,
    ) -> u64 {
        let [first, second] = count_min_two_buckets(key, self.layout);
        let first_depth = self.cell_2bit(first);
        let second_depth = self.cell_2bit(second);
        let min_depth = first_depth.min(second_depth);
        if min_depth >= self.max_count {
            return min_depth;
        }
        let target = min_depth
            .saturating_add(count.min(self.max_count))
            .min(self.max_count);
        if first_depth < target {
            self.set_cell_2bit_with_previous(first, first_depth, target);
        }
        if second_depth < target {
            self.set_cell_2bit_with_previous(second, second_depth, target);
        }
        min_depth
    }

    fn increment_independent_and_return_unincremented(&mut self, key: &KmerKey, count: u64) -> u64 {
        if count == 0 {
            return self.depth(key);
        }
        let increment = count.min(self.max_count);
        let mut previous_min = self.max_count;
        let mut slots = [0usize; 16];
        fill_count_min_buckets(key, self.hashes, self.layout, &mut slots);
        for slot in slots.iter().take(self.hashes) {
            let previous = self.cell(*slot);
            previous_min = previous_min.min(previous);
            let next = previous.saturating_add(increment).min(self.max_count);
            self.set_cell_with_previous(*slot, previous, next);
        }
        previous_min
    }

    fn add_key_counts(&mut self, counts: &CountMap) {
        if self.update_mode == CountMinUpdateMode::Conservative
            && self.bits == 16
            && self.hashes == 3
        {
            for (key, count) in counts {
                let _ =
                    self.increment_16bit_3hash_conservative_and_return_unincremented(key, *count);
            }
            return;
        }
        for (key, count) in counts {
            self.add_key_count(key, *count);
        }
    }

    fn add_key_increments(&mut self, key_increments: u64) {
        self.increments = self.increments.saturating_add(key_increments);
    }

    fn depth_16bit_3hash(&self, key: &KmerKey) -> u64 {
        let [first, second, third] = count_min_three_buckets_raw(raw_kmer_key(key), self.layout);
        self.cell_16bit(first)
            .min(self.cell_16bit(second))
            .min(self.cell_16bit(third))
    }

    fn occupied_slots_at_least(&self, min_depth: u64) -> usize {
        if min_depth > self.max_count {
            return 0;
        }
        if min_depth <= 1 {
            return self.occupied_slots;
        }
        let min_depth = min_depth.max(1);
        if let Some(slots) = &self.tracked_slots {
            return slots
                .par_iter()
                .filter(|&&slot| self.cell(slot) >= min_depth)
                .count();
        }
        (0..self.cells)
            .into_par_iter()
            .filter(|&slot| self.cell(slot) >= min_depth)
            .count()
    }

    fn cell(&self, slot: usize) -> u64 {
        if self.bits == 64 {
            return self.words[slot];
        }
        if self.bits == 16 {
            return self.cell_16bit(slot);
        }
        if self.bits == 2 {
            return self.cell_2bit(slot);
        }
        let bit = slot * self.bits as usize;
        let word = bit / 64;
        let offset = bit % 64;
        let mask = (1u64 << self.bits) - 1;
        if offset + self.bits as usize <= 64 {
            (self.words[word] >> offset) & mask
        } else {
            let low_bits = 64 - offset;
            let high_bits = self.bits as usize - low_bits;
            let low = self.words[word] >> offset;
            let high = self.words[word + 1] & ((1u64 << high_bits) - 1);
            ((high << low_bits) | low) & mask
        }
    }

    fn cell_16bit(&self, slot: usize) -> u64 {
        let word = slot >> 2;
        let offset = (slot & 3) << 4;
        (self.words[word] >> offset) & 0xffff
    }

    fn cell_2bit(&self, slot: usize) -> u64 {
        let word = slot >> 5;
        let offset = (slot & 31) << 1;
        (self.words[word] >> offset) & 3
    }

    #[cfg(test)]
    fn set_cell(&mut self, slot: usize, value: u64) {
        let previous = self.cell(slot);
        self.set_cell_with_previous(slot, previous, value);
    }

    fn set_cell_with_previous(&mut self, slot: usize, previous: u64, value: u64) {
        let value = value.min(self.max_count);
        self.set_cell_raw(slot, value);
        self.note_cell_transition(previous, value, slot);
    }

    fn set_cell_raw(&mut self, slot: usize, value: u64) {
        if self.bits == 64 {
            self.words[slot] = value;
            return;
        }
        if self.bits == 16 {
            self.set_cell_16bit_raw(slot, value);
            return;
        }
        if self.bits == 2 {
            self.set_cell_2bit_raw(slot, value);
            return;
        }
        let bit = slot * self.bits as usize;
        let word = bit / 64;
        let offset = bit % 64;
        let mask = (1u64 << self.bits) - 1;
        if offset + self.bits as usize <= 64 {
            let shifted_mask = mask << offset;
            self.words[word] = (self.words[word] & !shifted_mask) | ((value & mask) << offset);
        } else {
            let low_bits = 64 - offset;
            let high_bits = self.bits as usize - low_bits;
            let low_mask = ((1u64 << low_bits) - 1) << offset;
            self.words[word] =
                (self.words[word] & !low_mask) | ((value & ((1u64 << low_bits) - 1)) << offset);
            let high_mask = (1u64 << high_bits) - 1;
            self.words[word + 1] =
                (self.words[word + 1] & !high_mask) | ((value >> low_bits) & high_mask);
        }
    }

    fn set_cell_16bit_raw(&mut self, slot: usize, value: u64) {
        let word = slot >> 2;
        let offset = (slot & 3) << 4;
        let shifted_mask = 0xffffu64 << offset;
        self.words[word] = (self.words[word] & !shifted_mask) | ((value & 0xffff) << offset);
    }

    fn set_cell_16bit_with_previous(&mut self, slot: usize, previous: u64, value: u64) {
        let value = value.min(self.max_count);
        self.set_cell_16bit_raw(slot, value);
        self.note_cell_transition(previous, value, slot);
    }

    fn set_cell_2bit_with_previous(&mut self, slot: usize, previous: u64, value: u64) {
        let value = value.min(self.max_count);
        self.set_cell_2bit_raw(slot, value);
        self.note_cell_transition(previous, value, slot);
    }

    fn set_cell_2bit_raw(&mut self, slot: usize, value: u64) {
        let word = slot >> 5;
        let offset = (slot & 31) << 1;
        let shifted_mask = 3u64 << offset;
        self.words[word] = (self.words[word] & !shifted_mask) | ((value & 3) << offset);
    }

    fn note_cell_transition(&mut self, previous: u64, value: u64, slot: usize) {
        match (previous == 0, value == 0) {
            (true, false) => {
                self.occupied_slots = self.occupied_slots.saturating_add(1);
                if let Some(slots) = &mut self.tracked_slots {
                    if slots.len() < PACKED_SKETCH_TRACKED_SLOT_LIMIT {
                        slots.push(slot);
                    } else {
                        self.tracked_slots = None;
                    }
                }
            }
            (false, true) => {
                self.occupied_slots = self.occupied_slots.saturating_sub(1);
                if let Some(slots) = &mut self.tracked_slots
                    && let Some(index) = slots.iter().position(|&tracked| tracked == slot)
                {
                    slots.swap_remove(index);
                }
            }
            _ => {}
        }
    }

    #[cfg(test)]
    fn depth_hist(&self, hist_len: usize) -> Vec<u64> {
        let Some(last_index) = hist_len.checked_sub(1) else {
            return Vec::new();
        };
        if let Some(slots) = &self.tracked_slots {
            let mut hist = slots
                .par_iter()
                .fold(Vec::new, |mut local, &slot| {
                    add_depth_to_dynamic_hist(&mut local, self.cell(slot), last_index);
                    local
                })
                .reduce(Vec::new, merge_dynamic_depth_hist);
            hist.resize(hist_len, 0);
            return hist;
        }
        let mut hist = (0..self.cells)
            .into_par_iter()
            .fold(Vec::new, |mut local, slot| {
                add_depth_to_dynamic_hist(&mut local, self.cell(slot), last_index);
                local
            })
            .reduce(Vec::new, merge_dynamic_depth_hist);
        hist.resize(hist_len, 0);
        hist
    }

    fn sparse_depth_hist(&self, hist_len: usize) -> SparseHist {
        let Some(last_index) = hist_len.checked_sub(1) else {
            return SparseHist::default();
        };
        if let Some(slots) = &self.tracked_slots {
            return slots
                .par_iter()
                .fold(SparseHist::default, |mut local, &slot| {
                    add_depth_to_sparse_hist(&mut local, self.cell(slot), last_index);
                    local
                })
                .reduce(SparseHist::default, merge_sparse_depth_hist);
        }
        (0..self.cells)
            .into_par_iter()
            .fold(SparseHist::default, |mut local, slot| {
                add_depth_to_sparse_hist(&mut local, self.cell(slot), last_index);
                local
            })
            .reduce(SparseHist::default, merge_sparse_depth_hist)
    }
}

impl PrefilterCountMinSketch {
    fn max_count(&self) -> u64 {
        match self {
            Self::Packed(sketch) => sketch.max_count,
            Self::AtomicPacked(sketch) => sketch.max_count,
        }
    }

    #[cfg(test)]
    fn bits(&self) -> u8 {
        match self {
            Self::Packed(sketch) => sketch.bits,
            Self::AtomicPacked(sketch) => sketch.bits,
        }
    }

    #[cfg(test)]
    fn update_mode(&self) -> CountMinUpdateMode {
        match self {
            Self::Packed(sketch) => sketch.update_mode,
            Self::AtomicPacked(sketch) => sketch.update_mode,
        }
    }

    fn layout_summary(
        &self,
        table: &'static str,
        prefilter_limit: Option<u64>,
    ) -> SketchLayoutSummary {
        match self {
            Self::Packed(sketch) => sketch.layout_summary(table, prefilter_limit),
            Self::AtomicPacked(sketch) => sketch.layout_summary(table, prefilter_limit),
        }
    }
}

impl CountLookup for PrefilterCountMinSketch {
    fn depth(&self, key: &KmerKey) -> u64 {
        match self {
            Self::Packed(sketch) => sketch.depth(key),
            Self::AtomicPacked(sketch) => sketch.depth(key),
        }
    }

    fn unique_kmers(&self) -> usize {
        match self {
            Self::Packed(sketch) => sketch.unique_kmers(),
            Self::AtomicPacked(sketch) => sketch.unique_kmers(),
        }
    }

    fn unique_kmers_at_least(&self, min_depth: u64) -> usize {
        match self {
            Self::Packed(sketch) => sketch.unique_kmers_at_least(min_depth),
            Self::AtomicPacked(sketch) => sketch.unique_kmers_at_least(min_depth),
        }
    }
}

impl AtomicCountMinSketch {
    #[cfg(test)]
    fn new(cells: usize, hashes: usize) -> Result<Self> {
        Self::new_with_min_arrays(cells, hashes, BBTOOLS_KCOUNT_ARRAY_MIN_ARRAYS)
    }

    #[cfg(test)]
    fn new_with_min_arrays(cells: usize, hashes: usize, min_arrays: usize) -> Result<Self> {
        Self::new_with_min_arrays_and_update_mode(
            cells,
            hashes,
            min_arrays,
            CountMinUpdateMode::Conservative,
            BBTOOLS_KCOUNT_ARRAY_FIRST_MASK_SEED,
        )
    }

    fn new_with_min_arrays_and_update_mode(
        cells: usize,
        hashes: usize,
        min_arrays: usize,
        update_mode: CountMinUpdateMode,
        mask_seed: u64,
    ) -> Result<Self> {
        let cells = cells.max(1);
        let hashes = hashes.clamp(1, BBTOOLS_KCOUNT_ARRAY_MAX_HASHES);
        let layout =
            KCountArrayLayout::new_with_min_arrays_and_mask_seed(cells, 32, min_arrays, mask_seed);
        let cells_by_hash =
            zeroed_atomic_u32_vec(cells).context("allocating atomic count-min sketch")?;
        let locks = atomic_count_min_locks(update_mode)?;
        Ok(Self {
            cells,
            hashes,
            max_count: i32::MAX as u32,
            layout,
            update_mode,
            parallel_replay: false,
            cells_by_hash,
            locks,
            increments: AtomicU64::new(0),
            occupied_slots: AtomicUsize::new(0),
        })
    }

    fn with_parallel_replay(mut self, parallel_replay: bool) -> Self {
        self.parallel_replay = parallel_replay;
        self
    }

    fn layout_summary(
        &self,
        table: &'static str,
        prefilter_limit: Option<u64>,
    ) -> SketchLayoutSummary {
        SketchLayoutSummary {
            table,
            kind: "atomic",
            cells: self.cells,
            hashes: self.hashes,
            bits: 32,
            arrays: self.layout.array_count(),
            cells_per_array: self.layout.cells_per_array,
            mask_seed: self.layout.mask_seed,
            update_mode: self.update_mode.as_str(),
            max_count: u64::from(self.max_count),
            memory_bytes: self
                .cells_by_hash
                .len()
                .saturating_mul(std::mem::size_of::<AtomicU32>())
                .saturating_add(
                    self.locks
                        .len()
                        .saturating_mul(std::mem::size_of::<Mutex<()>>()),
                )
                .saturating_add(std::mem::size_of::<AtomicUsize>()),
            prefilter_limit,
        }
    }

    fn increment_key(&self, key: &KmerKey) {
        self.add_key_count(key, 1);
    }

    fn add_key_count(&self, key: &KmerKey, count: u64) {
        let (_, newly_occupied) = self.increment_and_count_newly_occupied(key, count);
        self.add_occupied_slots(newly_occupied);
    }

    #[cfg(test)]
    fn increment_and_return_unincremented(&self, key: &KmerKey, count: u64) -> u64 {
        let (previous_min, newly_occupied) = self.increment_and_count_newly_occupied(key, count);
        self.add_occupied_slots(newly_occupied);
        previous_min
    }

    fn add_key_count_counting_newly_occupied(&self, key: &KmerKey, count: u64) -> usize {
        self.increment_and_count_newly_occupied(key, count).1
    }

    fn add_key_count_unlocked_counting_newly_occupied(&self, key: &KmerKey, count: u64) -> usize {
        if self.update_mode == CountMinUpdateMode::Independent {
            self.increment_independent_and_count_newly_occupied(key, count)
                .1
        } else {
            self.increment_conservative_unlocked_and_count_newly_occupied(key, count)
                .1
        }
    }

    fn increment_and_count_newly_occupied(&self, key: &KmerKey, count: u64) -> (u64, usize) {
        if count == 0 {
            return (self.depth(key), 0);
        }
        if self.update_mode == CountMinUpdateMode::Independent {
            return self.increment_independent_and_count_newly_occupied(key, count);
        }
        let _guard = self.lock_for_key(key);
        self.increment_conservative_unlocked_and_count_newly_occupied(key, count)
    }

    fn increment_conservative_unlocked_and_count_newly_occupied(
        &self,
        key: &KmerKey,
        count: u64,
    ) -> (u64, usize) {
        let target_increment = count.min(u64::from(self.max_count)) as u32;
        if self.hashes == 3 {
            return self.increment_conservative_three_unlocked_and_count_newly_occupied(
                key,
                target_increment,
            );
        }
        let mut slots = [0usize; 16];
        let mut min_depth = self.max_count;
        fill_count_min_buckets(key, self.hashes, self.layout, &mut slots);
        for slot in slots.iter().take(self.hashes) {
            min_depth = min_depth.min(self.cells_by_hash[*slot].load(Ordering::Relaxed));
        }
        if min_depth >= self.max_count {
            return (u64::from(min_depth), 0);
        }
        let target = min_depth
            .saturating_add(target_increment)
            .min(self.max_count);
        let mut previous_min = self.max_count;
        let mut newly_occupied = 0usize;
        for slot in slots.iter().take(self.hashes) {
            let (previous, cell_newly_occupied) =
                raise_atomic_cell_to_at_least(&self.cells_by_hash[*slot], target);
            previous_min = previous_min.min(previous);
            newly_occupied += usize::from(cell_newly_occupied);
        }
        (u64::from(previous_min), newly_occupied)
    }

    fn increment_conservative_three_unlocked_and_count_newly_occupied(
        &self,
        key: &KmerKey,
        target_increment: u32,
    ) -> (u64, usize) {
        let [first, second, third] = count_min_three_buckets(key, self.layout);
        let first_depth = self.cells_by_hash[first].load(Ordering::Relaxed);
        let second_depth = self.cells_by_hash[second].load(Ordering::Relaxed);
        let third_depth = self.cells_by_hash[third].load(Ordering::Relaxed);
        let min_depth = first_depth.min(second_depth).min(third_depth);
        if min_depth >= self.max_count {
            return (u64::from(min_depth), 0);
        }
        let target = min_depth
            .saturating_add(target_increment)
            .min(self.max_count);
        let (first_previous, first_new) =
            raise_atomic_cell_to_at_least(&self.cells_by_hash[first], target);
        let (second_previous, second_new) =
            raise_atomic_cell_to_at_least(&self.cells_by_hash[second], target);
        let (third_previous, third_new) =
            raise_atomic_cell_to_at_least(&self.cells_by_hash[third], target);
        (
            u64::from(first_previous.min(second_previous).min(third_previous)),
            usize::from(first_new) + usize::from(second_new) + usize::from(third_new),
        )
    }

    fn lock_for_key(&self, key: &KmerKey) -> std::sync::MutexGuard<'_, ()> {
        let lock_index = kcount_array_lock_index(key);
        self.locks[lock_index]
            .lock()
            .unwrap_or_else(|poisoned| poisoned.into_inner())
    }

    fn increment_independent_and_count_newly_occupied(
        &self,
        key: &KmerKey,
        count: u64,
    ) -> (u64, usize) {
        if count == 0 {
            return (self.depth(key), 0);
        }
        let increment = count.min(u64::from(self.max_count)) as u32;
        let mut previous_min = self.max_count;
        let mut newly_occupied = 0usize;
        let mut slots = [0usize; 16];
        fill_count_min_buckets(key, self.hashes, self.layout, &mut slots);
        for slot in slots.iter().take(self.hashes) {
            let (previous, cell_newly_occupied) = increment_atomic_cell_saturating(
                &self.cells_by_hash[*slot],
                increment,
                self.max_count,
            );
            previous_min = previous_min.min(previous);
            newly_occupied += usize::from(cell_newly_occupied);
        }
        (u64::from(previous_min), newly_occupied)
    }

    fn add_key_counts(&self, counts: &CountMap) {
        let newly_occupied =
            if self.parallel_replay && counts.len() >= ATOMIC_SKETCH_PAR_REPLAY_MIN_KEYS {
                counts
                    .par_iter()
                    .map(|(key, count)| self.add_key_count_counting_newly_occupied(key, *count))
                    .sum()
            } else {
                counts
                    .iter()
                    .map(|(key, count)| {
                        self.add_key_count_unlocked_counting_newly_occupied(key, *count)
                    })
                    .sum()
            };
        self.add_occupied_slots(newly_occupied);
    }

    fn add_key_increments(&self, key_increments: u64) {
        self.increments.fetch_add(key_increments, Ordering::Relaxed);
    }

    fn add_occupied_slots(&self, newly_occupied: usize) {
        if newly_occupied > 0 {
            self.occupied_slots
                .fetch_add(newly_occupied, Ordering::Relaxed);
        }
    }

    fn occupied_slots_at_least(&self, min_depth: u64) -> usize {
        if min_depth > u64::from(self.max_count) {
            return 0;
        }
        if min_depth <= 1 {
            return self.occupied_slots.load(Ordering::Relaxed);
        }
        let min_depth = min_depth.max(1) as u32;
        self.cells_by_hash
            .par_iter()
            .filter(|cell| cell.load(Ordering::Relaxed) >= min_depth)
            .count()
    }

    #[cfg(test)]
    fn depth_hist(&self, hist_len: usize) -> Vec<u64> {
        let Some(last_index) = hist_len.checked_sub(1) else {
            return Vec::new();
        };
        let mut hist = self
            .cells_by_hash
            .par_iter()
            .fold(Vec::new, |mut local, cell| {
                add_depth_to_dynamic_hist(
                    &mut local,
                    u64::from(cell.load(Ordering::Relaxed)),
                    last_index,
                );
                local
            })
            .reduce(Vec::new, merge_dynamic_depth_hist);
        hist.resize(hist_len, 0);
        hist
    }

    fn sparse_depth_hist(&self, hist_len: usize) -> SparseHist {
        let Some(last_index) = hist_len.checked_sub(1) else {
            return SparseHist::default();
        };
        self.cells_by_hash
            .par_iter()
            .fold(SparseHist::default, |mut local, cell| {
                add_depth_to_sparse_hist(
                    &mut local,
                    u64::from(cell.load(Ordering::Relaxed)),
                    last_index,
                );
                local
            })
            .reduce(SparseHist::default, merge_sparse_depth_hist)
    }
}

impl AtomicPackedCountMinSketch {
    fn new_with_min_arrays_and_update_mode(
        cells: usize,
        hashes: usize,
        bits: u8,
        min_arrays: usize,
        update_mode: CountMinUpdateMode,
        mask_seed: u64,
    ) -> Result<Self> {
        let cells = cells.max(1);
        let hashes = hashes.clamp(1, BBTOOLS_KCOUNT_ARRAY_MAX_HASHES);
        ensure!(
            bits.is_power_of_two() && bits <= 64,
            "atomic packed count-min sketches require power-of-two cell bits up to 64"
        );
        let layout = KCountArrayLayout::new_with_min_arrays_and_mask_seed(
            cells, bits, min_arrays, mask_seed,
        );
        let word_count = if bits == 64 {
            cells
        } else {
            let cells_per_word = 64 / bits as usize;
            cells.div_ceil(cells_per_word)
        };
        let words = zeroed_atomic_u64_vec(word_count)
            .context("allocating atomic packed count-min sketch")?;
        let locks = atomic_count_min_locks(update_mode)?;
        Ok(Self {
            cells,
            hashes,
            bits,
            max_count: count_min_max_count(bits),
            layout,
            update_mode,
            words,
            locks,
            increments: AtomicU64::new(0),
            occupied_slots: AtomicUsize::new(0),
        })
    }

    fn layout_summary(
        &self,
        table: &'static str,
        prefilter_limit: Option<u64>,
    ) -> SketchLayoutSummary {
        SketchLayoutSummary {
            table,
            kind: "atomic_packed",
            cells: self.cells,
            hashes: self.hashes,
            bits: self.bits,
            arrays: self.layout.array_count(),
            cells_per_array: self.layout.cells_per_array,
            mask_seed: self.layout.mask_seed,
            update_mode: self.update_mode.as_str(),
            max_count: self.max_count,
            memory_bytes: self
                .words
                .len()
                .saturating_mul(std::mem::size_of::<AtomicU64>())
                .saturating_add(
                    self.locks
                        .len()
                        .saturating_mul(std::mem::size_of::<Mutex<()>>()),
                )
                .saturating_add(std::mem::size_of::<AtomicUsize>()),
            prefilter_limit,
        }
    }

    #[cfg(test)]
    fn add_key_count(&self, key: &KmerKey, count: u64) {
        let (_, newly_occupied) = self.increment_and_count_newly_occupied(key, count);
        self.add_occupied_slots(newly_occupied);
    }

    fn add_key_count_counting_newly_occupied(&self, key: &KmerKey, count: u64) -> usize {
        self.increment_and_count_newly_occupied(key, count).1
    }

    fn increment_and_count_newly_occupied(&self, key: &KmerKey, count: u64) -> (u64, usize) {
        if count == 0 {
            return (self.depth(key), 0);
        }
        if self.update_mode == CountMinUpdateMode::Independent {
            return self.increment_independent_and_count_newly_occupied(key, count);
        }
        let _guard = self.lock_for_key(key);
        let target_increment = count.min(self.max_count);
        let mut slots = [0usize; 16];
        let mut min_depth = self.max_count;
        fill_count_min_buckets(key, self.hashes, self.layout, &mut slots);
        for slot in slots.iter().take(self.hashes) {
            min_depth = min_depth.min(self.cell(*slot));
        }
        if min_depth >= self.max_count {
            return (min_depth, 0);
        }
        let target = min_depth
            .saturating_add(target_increment)
            .min(self.max_count);
        let mut previous_min = self.max_count;
        let mut newly_occupied = 0usize;
        for slot in slots.iter().take(self.hashes) {
            let (previous, cell_newly_occupied) = self.raise_cell_to_at_least(*slot, target);
            previous_min = previous_min.min(previous);
            newly_occupied += usize::from(cell_newly_occupied);
        }
        (previous_min, newly_occupied)
    }

    fn increment_independent_and_count_newly_occupied(
        &self,
        key: &KmerKey,
        count: u64,
    ) -> (u64, usize) {
        if count == 0 {
            return (self.depth(key), 0);
        }
        let increment = count.min(self.max_count);
        let mut previous_min = self.max_count;
        let mut newly_occupied = 0usize;
        let mut slots = [0usize; 16];
        fill_count_min_buckets(key, self.hashes, self.layout, &mut slots);
        for slot in slots.iter().take(self.hashes) {
            let (previous, cell_newly_occupied) = self.increment_cell_saturating(*slot, increment);
            previous_min = previous_min.min(previous);
            newly_occupied += usize::from(cell_newly_occupied);
        }
        (previous_min, newly_occupied)
    }

    fn add_key_increments(&self, key_increments: u64) {
        self.increments.fetch_add(key_increments, Ordering::Relaxed);
    }

    fn add_occupied_slots(&self, newly_occupied: usize) {
        if newly_occupied > 0 {
            self.occupied_slots
                .fetch_add(newly_occupied, Ordering::Relaxed);
        }
    }

    fn lock_for_key(&self, key: &KmerKey) -> std::sync::MutexGuard<'_, ()> {
        let lock_index = kcount_array_lock_index(key);
        self.locks[lock_index]
            .lock()
            .unwrap_or_else(|poisoned| poisoned.into_inner())
    }

    fn cell(&self, slot: usize) -> u64 {
        let position = self.cell_position(slot);
        (self.words[position.word].load(Ordering::Relaxed) >> position.shift) & position.mask
    }

    fn raise_cell_to_at_least(&self, slot: usize, target: u64) -> (u64, bool) {
        let target = target.min(self.max_count);
        let position = self.cell_position(slot);
        let cell = &self.words[position.word];
        let mut current = cell.load(Ordering::Relaxed);
        loop {
            let previous = (current >> position.shift) & position.mask;
            if previous >= target {
                return (previous, false);
            }
            let next = replace_packed_cell(current, position, target);
            match cell.compare_exchange_weak(current, next, Ordering::Relaxed, Ordering::Relaxed) {
                Ok(_) => return (previous, previous == 0 && target > 0),
                Err(observed) => current = observed,
            }
        }
    }

    fn increment_cell_saturating(&self, slot: usize, increment: u64) -> (u64, bool) {
        let increment = increment.min(self.max_count);
        let position = self.cell_position(slot);
        let cell = &self.words[position.word];
        let mut current = cell.load(Ordering::Relaxed);
        loop {
            let previous = (current >> position.shift) & position.mask;
            if previous >= self.max_count {
                return (previous, false);
            }
            let next_value = previous.saturating_add(increment).min(self.max_count);
            let next = replace_packed_cell(current, position, next_value);
            match cell.compare_exchange_weak(current, next, Ordering::Relaxed, Ordering::Relaxed) {
                Ok(_) => return (previous, previous == 0 && next_value > 0),
                Err(observed) => current = observed,
            }
        }
    }

    fn cell_position(&self, slot: usize) -> PackedCellPosition {
        if self.bits == 64 {
            return PackedCellPosition {
                word: slot,
                shift: 0,
                mask: u64::MAX,
            };
        }
        let cells_per_word = 64 / self.bits as usize;
        let word = slot / cells_per_word;
        let shift = (slot % cells_per_word) * self.bits as usize;
        let mask = (1u64 << self.bits) - 1;
        PackedCellPosition { word, shift, mask }
    }

    fn occupied_slots_at_least(&self, min_depth: u64) -> usize {
        if min_depth > self.max_count {
            return 0;
        }
        if min_depth <= 1 {
            return self.occupied_slots.load(Ordering::Relaxed);
        }
        let min_depth = min_depth.max(1);
        (0..self.cells)
            .into_par_iter()
            .filter(|&slot| self.cell(slot) >= min_depth)
            .count()
    }
}

#[derive(Debug, Clone, Copy)]
struct PackedCellPosition {
    word: usize,
    shift: usize,
    mask: u64,
}

fn replace_packed_cell(word: u64, position: PackedCellPosition, value: u64) -> u64 {
    let shifted_mask = position.mask << position.shift;
    (word & !shifted_mask) | ((value & position.mask) << position.shift)
}

fn raise_atomic_cell_to_at_least(cell: &AtomicU32, target: u32) -> (u32, bool) {
    let mut current = cell.load(Ordering::Relaxed);
    loop {
        if current >= target {
            return (current, false);
        }
        match cell.compare_exchange_weak(current, target, Ordering::Relaxed, Ordering::Relaxed) {
            Ok(_) => return (current, current == 0 && target > 0),
            Err(observed) => current = observed,
        }
    }
}

fn increment_atomic_cell_saturating(
    cell: &AtomicU32,
    increment: u32,
    max_count: u32,
) -> (u32, bool) {
    let mut current = cell.load(Ordering::Relaxed);
    loop {
        if current >= max_count {
            return (current, false);
        }
        let next = current.saturating_add(increment).min(max_count);
        match cell.compare_exchange_weak(current, next, Ordering::Relaxed, Ordering::Relaxed) {
            Ok(_) => return (current, current == 0 && next > 0),
            Err(observed) => current = observed,
        }
    }
}

fn atomic_count_min_locks(update_mode: CountMinUpdateMode) -> Result<Vec<Mutex<()>>> {
    if update_mode == CountMinUpdateMode::Independent {
        return Ok(Vec::new());
    }
    let mut locks = Vec::new();
    locks
        .try_reserve_exact(BBTOOLS_KCOUNT_ARRAY_LOCKS)
        .context("allocating atomic count-min sketch locks")?;
    locks.resize_with(BBTOOLS_KCOUNT_ARRAY_LOCKS, || Mutex::new(()));
    Ok(locks)
}

impl CountLookup for PackedCountMinSketch {
    fn depth(&self, key: &KmerKey) -> u64 {
        if self.bits == 16 && self.hashes == 3 {
            return self.depth_16bit_3hash(key);
        }
        let mut slots = [0usize; 16];
        fill_count_min_buckets(key, self.hashes, self.layout, &mut slots);
        slots
            .iter()
            .take(self.hashes)
            .map(|&slot| self.cell(slot))
            .min()
            .unwrap_or(0)
    }

    fn unique_kmers(&self) -> usize {
        self.unique_kmers_at_least(1)
    }

    fn unique_kmers_at_least(&self, min_depth: u64) -> usize {
        let occupied = self.occupied_slots_at_least(min_depth);
        estimate_unique_kmers_from_occupied(self.cells, occupied, self.hashes, self.increments)
    }
}

impl CountLookup for AtomicCountMinSketch {
    fn depth(&self, key: &KmerKey) -> u64 {
        let mut slots = [0usize; 16];
        fill_count_min_buckets(key, self.hashes, self.layout, &mut slots);
        slots
            .iter()
            .take(self.hashes)
            .map(|&slot| u64::from(self.cells_by_hash[slot].load(Ordering::Relaxed)))
            .min()
            .unwrap_or(0)
    }

    fn unique_kmers(&self) -> usize {
        self.unique_kmers_at_least(1)
    }

    fn unique_kmers_at_least(&self, min_depth: u64) -> usize {
        let occupied = self.occupied_slots_at_least(min_depth);
        let increments = self.increments.load(Ordering::Relaxed);
        estimate_unique_kmers_from_occupied(self.cells, occupied, self.hashes, increments)
    }
}

impl CountLookup for AtomicPackedCountMinSketch {
    fn depth(&self, key: &KmerKey) -> u64 {
        let mut slots = [0usize; 16];
        fill_count_min_buckets(key, self.hashes, self.layout, &mut slots);
        slots
            .iter()
            .take(self.hashes)
            .map(|&slot| self.cell(slot))
            .min()
            .unwrap_or(0)
    }

    fn unique_kmers(&self) -> usize {
        self.unique_kmers_at_least(1)
    }

    fn unique_kmers_at_least(&self, min_depth: u64) -> usize {
        let occupied = self.occupied_slots_at_least(min_depth);
        let increments = self.increments.load(Ordering::Relaxed);
        estimate_unique_kmers_from_occupied(self.cells, occupied, self.hashes, increments)
    }
}

fn estimate_unique_kmers_from_occupied(
    total_slots: usize,
    occupied_slots: usize,
    hashes: usize,
    increments: u64,
) -> usize {
    // BBTools' KCountArray estimates cardinality from one shared used-cell
    // fraction, adjusted by hash count.
    if occupied_slots == 0 || total_slots == 0 {
        return 0;
    }
    let increment_cap = usize_from_u64_saturating(increments);
    if occupied_slots >= total_slots {
        return increment_cap;
    }
    let used_fraction = occupied_slots as f64 / total_slots as f64;
    let hash_count = hashes.max(1) as f64;
    let one_hash_fraction = 1.0 - (1.0 - used_fraction).powf(1.0 / hash_count);
    let estimate = (-(total_slots as f64) * (1.0 - one_hash_fraction).ln()).round();
    let estimate = estimate.max(1.0) as usize;
    estimate.min(increment_cap)
}

fn usize_from_u64_saturating(value: u64) -> usize {
    usize::try_from(value).unwrap_or(usize::MAX)
}

fn count_min_max_count(bits: u8) -> u64 {
    if bits >= 31 {
        i32::MAX as u64
    } else {
        (1u64 << bits.max(1)) - 1
    }
}

impl KCountArrayLayout {
    #[cfg(test)]
    fn new(cells: usize, bits: u8) -> Self {
        Self::new_with_min_arrays(cells, bits, BBTOOLS_KCOUNT_ARRAY_MIN_ARRAYS)
    }

    #[cfg(test)]
    fn new_with_min_arrays(cells: usize, bits: u8, min_arrays: usize) -> Self {
        Self::new_with_min_arrays_and_mask_seed(
            cells,
            bits,
            min_arrays,
            BBTOOLS_KCOUNT_ARRAY_FIRST_MASK_SEED,
        )
    }

    fn new_with_min_arrays_and_mask_seed(
        cells: usize,
        bits: u8,
        min_arrays: usize,
        mask_seed: u64,
    ) -> Self {
        let cells = cells.max(1);
        let arrays = kcount_array_count(cells, bits, min_arrays);
        let cells_per_array = (cells / arrays).max(1);
        Self {
            array_mask: arrays.saturating_sub(1) as u64,
            array_bits: arrays.trailing_zeros(),
            cells_per_array,
            mask_seed,
            masks: bbtools_hash_masks(mask_seed),
        }
    }

    fn array_count(self) -> usize {
        self.array_mask.saturating_add(1) as usize
    }

    fn bucket(self, hashed: u64) -> usize {
        if self.cells_per_array <= 1 && self.array_mask == 0 {
            return 0;
        }
        let array_num = (hashed & self.array_mask) as usize;
        let cell = ((hashed >> self.array_bits) % self.cells_per_array as u64) as usize;
        array_num * self.cells_per_array + cell
    }
}

#[cfg(test)]
fn count_min_bucket(key: &KmerKey, hash_index: usize, cells: usize) -> usize {
    count_min_bucket_with_layout(key, hash_index, KCountArrayLayout::new(cells, 32))
}

#[cfg(test)]
fn count_min_bucket_with_layout(
    key: &KmerKey,
    hash_index: usize,
    layout: KCountArrayLayout,
) -> usize {
    let hashed = bbtools_count_min_row_hash_with_masks(raw_kmer_key(key), hash_index, layout.masks);
    layout.bucket(hashed)
}

#[inline]
fn fill_count_min_buckets(
    key: &KmerKey,
    hashes: usize,
    layout: KCountArrayLayout,
    slots: &mut [usize; 16],
) {
    let hashes = hashes.min(slots.len());
    if hashes == 0 {
        return;
    }
    let mut hashed = bbtools_mask_hash_with_masks(raw_kmer_key(key), 0, layout.masks);
    slots[0] = layout.bucket(hashed);
    for (hash_index, slot) in slots.iter_mut().enumerate().take(hashes).skip(1) {
        hashed = hashed.rotate_right(BBTOOLS_HASH_BITS);
        hashed = bbtools_mask_hash_with_masks(hashed, hash_index, layout.masks);
        *slot = layout.bucket(hashed);
    }
}

#[inline]
fn count_min_three_buckets(key: &KmerKey, layout: KCountArrayLayout) -> [usize; 3] {
    count_min_three_buckets_raw(raw_kmer_key(key), layout)
}

#[inline]
fn count_min_three_buckets_raw(raw_key: u64, layout: KCountArrayLayout) -> [usize; 3] {
    let mut hashed = bbtools_mask_hash_with_masks(raw_key, 0, layout.masks);
    let first = layout.bucket(hashed);
    hashed = bbtools_mask_hash_with_masks(hashed.rotate_right(BBTOOLS_HASH_BITS), 1, layout.masks);
    let second = layout.bucket(hashed);
    hashed = bbtools_mask_hash_with_masks(hashed.rotate_right(BBTOOLS_HASH_BITS), 2, layout.masks);
    [first, second, layout.bucket(hashed)]
}

#[inline]
fn count_min_two_buckets(key: &KmerKey, layout: KCountArrayLayout) -> [usize; 2] {
    let mut hashed = bbtools_mask_hash_with_masks(raw_kmer_key(key), 0, layout.masks);
    let first = layout.bucket(hashed);
    hashed = bbtools_mask_hash_with_masks(hashed.rotate_right(BBTOOLS_HASH_BITS), 1, layout.masks);
    [first, layout.bucket(hashed)]
}

#[cfg(test)]
fn bbtools_count_min_row_hash_with_masks(
    raw_key: u64,
    hash_index: usize,
    masks: &BbtoolsHashMaskTable,
) -> u64 {
    let mut key = bbtools_mask_hash_with_masks(raw_key, 0, masks);
    for row in 1..=hash_index {
        key = key.rotate_right(BBTOOLS_HASH_BITS);
        key = bbtools_mask_hash_with_masks(key, row, masks);
    }
    key
}

#[cfg(test)]
fn bbtools_mask_hash(key: u64, row: usize, mask_seed: u64) -> u64 {
    let masks = bbtools_hash_masks(mask_seed);
    bbtools_mask_hash_with_masks(key, row, masks)
}

#[inline]
fn bbtools_mask_hash_with_masks(mut key: u64, row: usize, masks: &BbtoolsHashMaskTable) -> u64 {
    let row = row & 7;
    let mut cell =
        ((key & BBTOOLS_LONG_MAX_VALUE) % (BBTOOLS_HASH_ARRAY_LENGTH as u64 - 1)) as usize;

    if row == 0 {
        key ^= masks[(row + 4) & 7][cell];
        cell = ((key >> 5) & BBTOOLS_HASH_CELL_MASK) as usize;
    }

    key ^ masks[row][cell]
}

fn bbtools_hash_masks(mask_seed: u64) -> BbtoolsHashMaskRef {
    static SEED0_MASKS: OnceLock<BbtoolsHashMaskTable> = OnceLock::new();
    static SEED7_MASKS: OnceLock<BbtoolsHashMaskTable> = OnceLock::new();
    static SEED14_MASKS: OnceLock<BbtoolsHashMaskTable> = OnceLock::new();
    static OTHER_MASKS: OnceLock<Mutex<BbtoolsHashMaskCache>> = OnceLock::new();
    match mask_seed {
        BBTOOLS_KCOUNT_ARRAY_FIRST_MASK_SEED => {
            SEED0_MASKS.get_or_init(|| make_bbtools_hash_masks(mask_seed))
        }
        BBTOOLS_KCOUNT_ARRAY_SECOND_MASK_SEED => {
            SEED7_MASKS.get_or_init(|| make_bbtools_hash_masks(mask_seed))
        }
        BBTOOLS_KCOUNT_ARRAY_THIRD_MASK_SEED => {
            SEED14_MASKS.get_or_init(|| make_bbtools_hash_masks(mask_seed))
        }
        _ => {
            let cache = OTHER_MASKS.get_or_init(|| Mutex::new(FxHashMap::default()));
            let mut cache = cache
                .lock()
                .unwrap_or_else(|poisoned| poisoned.into_inner());
            if let Some(&masks) = cache.get(&mask_seed) {
                return masks;
            }
            let masks = Box::leak(Box::new(make_bbtools_hash_masks(mask_seed)));
            cache.insert(mask_seed, masks);
            masks
        }
    }
}

fn make_bbtools_hash_masks(mask_seed: u64) -> BbtoolsHashMaskTable {
    let mut masks = [[0u64; BBTOOLS_HASH_ARRAY_LENGTH]; 8];
    let mut rng = BbtoolsXoshiro::new(mask_seed);
    for row_masks in &mut masks {
        fill_bbtools_hash_mask_row(row_masks, &mut rng);
    }
    masks
}

fn fill_bbtools_hash_mask_row(
    row_masks: &mut [u64; BBTOOLS_HASH_ARRAY_LENGTH],
    rng: &mut BbtoolsXoshiro,
) {
    let mut low_cells = [0u8; BBTOOLS_HASH_ARRAY_LENGTH];
    let mut rotated_cells = [0u8; BBTOOLS_HASH_ARRAY_LENGTH];

    for mask in row_masks {
        let (value, low_cell, rotated_cell) = loop {
            let mut value = rng.next_long();
            while (value & 0xffff_ffff).count_ones() < 16 {
                value |= 1u64 << rng.next_power_of_two_int(32);
            }
            while (value & 0xffff_ffff).count_ones() > 16 {
                value &= !(1u64 << rng.next_power_of_two_int(32));
            }
            while (value & 0xffff_ffff_0000_0000).count_ones() < 16 {
                value |= 1u64 << (rng.next_power_of_two_int(32) + 32);
            }
            while (value & 0xffff_ffff_0000_0000).count_ones() > 16 {
                value &= !(1u64 << (rng.next_power_of_two_int(32) + 32));
            }

            let low_cell = (value & BBTOOLS_HASH_CELL_MASK) as usize;
            let rotated_cell =
                (((value as i64) >> BBTOOLS_HASH_BITS) as u64 & BBTOOLS_HASH_CELL_MASK) as usize;
            if low_cells[low_cell] == 0 && rotated_cells[rotated_cell] == 0 {
                break (value & BBTOOLS_LONG_MAX_VALUE, low_cell, rotated_cell);
            }
        };

        *mask = value;
        low_cells[low_cell] = low_cells[low_cell].saturating_add(1);
        rotated_cells[rotated_cell] = rotated_cells[rotated_cell].saturating_add(1);
    }
}

struct BbtoolsXoshiro {
    s0: u64,
    s1: u64,
    s2: u64,
    s3: u64,
}

impl BbtoolsXoshiro {
    fn new(seed: u64) -> Self {
        let mut rng = Self {
            s0: seed,
            s1: Self::mix_seed(seed),
            s2: 0,
            s3: 0,
        };
        rng.s2 = Self::mix_seed(rng.s1);
        rng.s3 = Self::mix_seed(rng.s2);
        if rng.s0 == 0 && rng.s1 == 0 && rng.s2 == 0 && rng.s3 == 0 {
            rng.s0 = 0x5DEECE66D;
            rng.s1 = 0xB;
            rng.s2 = 0xCCA;
            rng.s3 = 0xF00;
        }
        for _ in 0..4 {
            rng.next_long();
        }
        rng
    }

    fn mix_seed(mut value: u64) -> u64 {
        value = value.wrapping_add(0x9E37_79B9_7F4A_7C15);
        value = (value ^ (value >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
        value = (value ^ (value >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
        value ^ (value >> 31)
    }

    fn next_long(&mut self) -> u64 {
        let result = self.s0.wrapping_add(self.s3);
        let t = self.s1 << 17;

        self.s2 ^= self.s0;
        self.s3 ^= self.s1;
        self.s1 ^= self.s2;
        self.s0 ^= self.s3;

        self.s2 ^= t;
        self.s3 = self.s3.rotate_left(45);

        result
    }

    fn next_power_of_two_int(&mut self, bound: u32) -> u32 {
        debug_assert!(bound.is_power_of_two());
        (self.next_long() as u32) & (bound - 1)
    }
}

fn new_count_map(config: &Config) -> CountMap {
    let mut counts = CountMap::default();
    if let Some(capacity) = count_map_capacity_hint(config) {
        let _ = counts.try_reserve(capacity);
    }
    counts
}

fn count_map_with_capacity(capacity: usize) -> CountMap {
    let mut counts = CountMap::default();
    if capacity > 0 {
        let _ = counts.try_reserve(capacity);
    }
    counts
}

fn count_chunk_local_map(
    config: &Config,
    pairs: &[(SequenceRecord, Option<SequenceRecord>)],
) -> CountMap {
    count_map_with_capacity(count_chunk_local_map_capacity(config, pairs))
}

fn count_chunk_local_map_capacity(
    config: &Config,
    pairs: &[(SequenceRecord, Option<SequenceRecord>)],
) -> usize {
    let total_windows: usize = pairs
        .iter()
        .map(|(r1, r2)| pair_kmer_window_capacity(config, r1, r2.as_ref()))
        .sum();
    if total_windows == 0 {
        return 0;
    }
    total_windows
        .div_ceil(rayon::current_num_threads().max(1))
        .clamp(64, COUNT_CHUNK_LOCAL_MAP_MAX_CAPACITY)
}

fn count_map_capacity_hint(config: &Config) -> Option<usize> {
    let explicit = config.table_initial_size;
    let prealloc = preallocation_capacity_hint(config);
    explicit.max(prealloc)
}

fn preallocation_capacity_hint(config: &Config) -> Option<usize> {
    let fraction = config.table_prealloc_fraction?;
    let reads = config.table_reads.or(config.max_reads)?;
    let reads = usize::try_from(reads).ok()?;
    if reads == 0 || fraction <= 0.0 {
        return None;
    }
    let mates = if config.in2.is_some() || config.interleaved {
        2usize
    } else {
        1usize
    };
    let kmers_per_read_hint = 100usize.saturating_sub(config.k).saturating_add(1).max(1);
    let raw = reads
        .saturating_mul(mates)
        .saturating_mul(kmers_per_read_hint);
    Some(((raw as f64) * fraction).ceil().max(1.0) as usize)
}

fn count_primary(config: &Config, counts: &mut CountMap) -> Result<()> {
    if let Some(paths) = primary_input_lists(config) {
        if let Some(first) = paths.first.first() {
            if let Some(second) = paths.second.as_ref().and_then(|paths| paths.first()) {
                count_paired_files(config, first, second, counts, config.table_reads)?;
            } else {
                count_single_file(config, first, counts, config.table_reads)?;
            }
        }
        for path in paths.first.iter().skip(1) {
            count_single_file(config, path, counts, None)?;
        }
        if let Some(second) = &paths.second {
            for path in second.iter().skip(1) {
                count_single_file(config, path, counts, None)?;
            }
        }
        return Ok(());
    }

    let mut readers = PrimaryReaders::open(config, config.table_reads)?;
    let mut chunk = Vec::with_capacity(COUNT_PARALLEL_CHUNK_SIZE);
    while let Some((r1, r2)) = readers.next_pair()? {
        chunk.push((r1, r2));
        if chunk.len() >= COUNT_PARALLEL_CHUNK_SIZE {
            increment_counts_from_pair_chunk(config, counts, &chunk);
            chunk.clear();
        }
    }
    if !chunk.is_empty() {
        increment_counts_from_pair_chunk(config, counts, &chunk);
    }
    Ok(())
}

fn count_primary_sketch(
    config: &Config,
    sketch: &mut PackedCountMinSketch,
    prefilter: Option<PrefilterGate<'_>>,
) -> Result<()> {
    if let Some(paths) = primary_input_lists(config) {
        if let Some(first) = paths.first.first() {
            if let Some(second) = paths.second.as_ref().and_then(|paths| paths.first()) {
                count_paired_files_sketch(
                    config,
                    first,
                    second,
                    sketch,
                    config.table_reads,
                    prefilter,
                )?;
            } else {
                count_single_file_sketch(config, first, sketch, config.table_reads, prefilter)?;
            }
        }
        for path in paths.first.iter().skip(1) {
            count_single_file_sketch(config, path, sketch, None, prefilter)?;
        }
        if let Some(second) = &paths.second {
            for path in second.iter().skip(1) {
                count_single_file_sketch(config, path, sketch, None, prefilter)?;
            }
        }
        return Ok(());
    }

    let mut readers = PrimaryReaders::open(config, config.table_reads)?;
    let mut chunk = Vec::with_capacity(COUNT_PARALLEL_CHUNK_SIZE);
    while let Some((r1, r2)) = readers.next_pair()? {
        chunk.push((r1, r2));
        if chunk.len() >= COUNT_PARALLEL_CHUNK_SIZE {
            increment_sketch_from_pair_chunk(config, sketch, &chunk, prefilter);
            chunk.clear();
        }
    }
    if !chunk.is_empty() {
        increment_sketch_from_pair_chunk(config, sketch, &chunk, prefilter);
    }
    Ok(())
}

fn count_primary_prefilter_sketch(
    config: &Config,
    sketch: &mut PrefilterCountMinSketch,
) -> Result<()> {
    match sketch {
        PrefilterCountMinSketch::Packed(sketch) => count_primary_sketch(config, sketch, None),
        PrefilterCountMinSketch::AtomicPacked(sketch) => {
            count_primary_atomic_packed_sketch(config, sketch)
        }
    }
}

fn count_primary_atomic_packed_sketch(
    config: &Config,
    sketch: &AtomicPackedCountMinSketch,
) -> Result<()> {
    if let Some(paths) = primary_input_lists(config) {
        if let Some(first) = paths.first.first() {
            if let Some(second) = paths.second.as_ref().and_then(|paths| paths.first()) {
                count_paired_files_atomic_packed_sketch(
                    config,
                    first,
                    second,
                    sketch,
                    config.table_reads,
                )?;
            } else {
                count_single_file_atomic_packed_sketch(config, first, sketch, config.table_reads)?;
            }
        }
        for path in paths.first.iter().skip(1) {
            count_single_file_atomic_packed_sketch(config, path, sketch, None)?;
        }
        if let Some(second) = &paths.second {
            for path in second.iter().skip(1) {
                count_single_file_atomic_packed_sketch(config, path, sketch, None)?;
            }
        }
        return Ok(());
    }

    let mut readers = PrimaryReaders::open(config, config.table_reads)?;
    let mut chunk = Vec::with_capacity(COUNT_PARALLEL_CHUNK_SIZE);
    while let Some((r1, r2)) = readers.next_pair()? {
        chunk.push((r1, r2));
        if chunk.len() >= COUNT_PARALLEL_CHUNK_SIZE {
            increment_atomic_packed_sketch_from_pair_chunk(config, sketch, &chunk);
            chunk.clear();
        }
    }
    if !chunk.is_empty() {
        increment_atomic_packed_sketch_from_pair_chunk(config, sketch, &chunk);
    }
    Ok(())
}

fn count_primary_atomic_sketch(
    config: &Config,
    sketch: &AtomicCountMinSketch,
    prefilter: Option<PrefilterGate<'_>>,
) -> Result<()> {
    if let Some(paths) = primary_input_lists(config) {
        if let Some(first) = paths.first.first() {
            if let Some(second) = paths.second.as_ref().and_then(|paths| paths.first()) {
                count_paired_files_atomic_sketch(
                    config,
                    first,
                    second,
                    sketch,
                    config.table_reads,
                    prefilter,
                )?;
            } else {
                count_single_file_atomic_sketch(
                    config,
                    first,
                    sketch,
                    config.table_reads,
                    prefilter,
                )?;
            }
        }
        for path in paths.first.iter().skip(1) {
            count_single_file_atomic_sketch(config, path, sketch, None, prefilter)?;
        }
        if let Some(second) = &paths.second {
            for path in second.iter().skip(1) {
                count_single_file_atomic_sketch(config, path, sketch, None, prefilter)?;
            }
        }
        return Ok(());
    }

    let mut readers = PrimaryReaders::open(config, config.table_reads)?;
    let mut chunk = Vec::with_capacity(COUNT_PARALLEL_CHUNK_SIZE);
    while let Some((r1, r2)) = readers.next_pair()? {
        chunk.push((r1, r2));
        if chunk.len() >= COUNT_PARALLEL_CHUNK_SIZE {
            increment_atomic_sketch_from_pair_chunk(config, sketch, &chunk, prefilter);
            chunk.clear();
        }
    }
    if !chunk.is_empty() {
        increment_atomic_sketch_from_pair_chunk(config, sketch, &chunk, prefilter);
    }
    Ok(())
}

fn count_primary_gpu_reduced_runs_sketch(
    config: &Config,
    sketch: &mut PackedCountMinSketch,
) -> Result<()> {
    for_each_gpu_reduced_chunk_run(config, |key, count| {
        sketch.add_key_count(&key, count);
        sketch.add_key_increments(count);
    })
}

fn count_primary_gpu_reduced_runs_atomic_sketch(
    config: &Config,
    sketch: &AtomicCountMinSketch,
) -> Result<()> {
    for_each_gpu_reduced_chunk_run(config, |key, count| {
        sketch.add_key_count(&key, count);
        sketch.add_key_increments(count);
    })
}

fn for_each_gpu_reduced_chunk_run<F>(config: &Config, mut f: F) -> Result<()>
where
    F: FnMut(KmerKey, u64),
{
    let helper = config
        .gpu_helper
        .as_ref()
        .context("gpucounting=t requires gpuhelper=<cuda_kmer_reduce_runs binary>")?;
    if !helper.exists() {
        bail!("gpuhelper does not exist: {}", helper.display());
    }
    ensure!(
        config.k <= 31,
        "gpucounting=t currently supports short k-mers only (k<=31)"
    );
    ensure!(
        !use_prefilter_collision_estimates(config),
        "gpucounting=t currently supports the main bounded sketch without prefilter=t"
    );
    let temp_dir = config.temp_dir.clone().unwrap_or_else(std::env::temp_dir);
    fs::create_dir_all(&temp_dir)
        .with_context(|| format!("creating GPU counting temp dir {}", temp_dir.display()))?;
    let token = format!(
        "{}_{}",
        std::process::id(),
        SystemTime::now()
            .duration_since(UNIX_EPOCH)
            .unwrap_or_default()
            .as_nanos()
    );
    let kmers_path = temp_dir.join(format!("bbnorm-rs-gpu-kmers-{token}.u64"));
    let runs_path = temp_dir.join(format!("bbnorm-rs-gpu-runs-{token}.bin"));
    let result = (|| {
        let mut readers = PrimaryReaders::open(config, config.table_reads)?;
        let mut persistent = config
            .gpu_persistent
            .then(|| PersistentGpuReducer::start(helper))
            .transpose()?;
        let mut chunk = Vec::with_capacity(COUNT_PARALLEL_CHUNK_SIZE);
        while let Some((r1, r2)) = readers.next_pair()? {
            chunk.push((r1, r2));
            if chunk.len() >= COUNT_PARALLEL_CHUNK_SIZE {
                if let Some(reducer) = &mut persistent {
                    reduce_gpu_pair_chunk_persistent(config, reducer, &chunk, &mut f)?;
                } else {
                    reduce_gpu_pair_chunk(config, helper, &kmers_path, &runs_path, &chunk, &mut f)?;
                }
                chunk.clear();
            }
        }
        if !chunk.is_empty() {
            if let Some(reducer) = &mut persistent {
                reduce_gpu_pair_chunk_persistent(config, reducer, &chunk, &mut f)?;
            } else {
                reduce_gpu_pair_chunk(config, helper, &kmers_path, &runs_path, &chunk, &mut f)?;
            }
        }
        if let Some(reducer) = persistent {
            reducer.finish()?;
        }
        Ok(())
    })();
    let _ = fs::remove_file(&kmers_path);
    let _ = fs::remove_file(&runs_path);
    result
}

fn reduce_gpu_pair_chunk<F>(
    config: &Config,
    helper: &Path,
    kmers_path: &Path,
    runs_path: &Path,
    pairs: &[(SequenceRecord, Option<SequenceRecord>)],
    f: &mut F,
) -> Result<()>
where
    F: FnMut(KmerKey, u64),
{
    write_pair_chunk_short_kmers(config, pairs, kmers_path)?;
    if fs::metadata(kmers_path)?.len() == 0 {
        return Ok(());
    }
    let status = Command::new(helper)
        .arg(kmers_path)
        .arg(runs_path)
        .status()
        .with_context(|| format!("running GPU helper {}", helper.display()))?;
    if !status.success() {
        bail!("GPU helper failed with status {status}");
    }
    replay_reduced_runs_file(runs_path, f)?;
    let _ = fs::remove_file(kmers_path);
    let _ = fs::remove_file(runs_path);
    Ok(())
}

fn reduce_gpu_pair_chunk_persistent<F>(
    config: &Config,
    reducer: &mut PersistentGpuReducer,
    pairs: &[(SequenceRecord, Option<SequenceRecord>)],
    f: &mut F,
) -> Result<()>
where
    F: FnMut(KmerKey, u64),
{
    let mut keys = Vec::new();
    collect_pair_chunk_short_kmers(config, pairs, &mut keys)?;
    if keys.is_empty() {
        return Ok(());
    }
    reducer.reduce(&keys, f)
}

fn write_pair_chunk_short_kmers(
    config: &Config,
    pairs: &[(SequenceRecord, Option<SequenceRecord>)],
    path: &Path,
) -> Result<()> {
    let mut writer = BufWriter::new(
        fs::File::create(path).with_context(|| format!("create {}", path.display()))?,
    );
    let mut keys = Vec::new();
    collect_pair_chunk_short_kmers(config, pairs, &mut keys)?;
    for raw in keys {
        writer.write_all(&raw.to_le_bytes())?;
    }
    writer.flush()?;
    Ok(())
}

fn collect_pair_chunk_short_kmers(
    config: &Config,
    pairs: &[(SequenceRecord, Option<SequenceRecord>)],
    out: &mut Vec<u64>,
) -> Result<()> {
    out.clear();
    let mut keys = Vec::new();
    for (r1, r2) in pairs {
        if config.remove_duplicate_kmers {
            fill_unique_pair_kmers(config, r1, r2.as_ref(), &mut keys);
            for key in &keys {
                out.push(short_kmer_raw(key)?);
            }
        } else {
            let mut write_error = None;
            for_each_kmer_for_record(r1, config, |key| match short_kmer_raw(&key) {
                Ok(raw) => out.push(raw),
                Err(err) => {
                    write_error = Some(err);
                }
            });
            if let Some(mate) = r2 {
                for_each_kmer_for_record(mate, config, |key| match short_kmer_raw(&key) {
                    Ok(raw) => out.push(raw),
                    Err(err) => {
                        write_error = Some(err);
                    }
                });
            }
            if let Some(err) = write_error {
                return Err(err);
            }
        }
    }
    Ok(())
}

fn short_kmer_raw(key: &KmerKey) -> Result<u64> {
    let KmerKey::Short(raw) = key else {
        bail!("GPU counting helper only accepts short k-mer keys");
    };
    Ok(*raw)
}

struct PersistentGpuReducer {
    child: Child,
    stdin: BufWriter<ChildStdin>,
    stdout: BufReader<ChildStdout>,
}

impl PersistentGpuReducer {
    fn start(helper: &Path) -> Result<Self> {
        let mut child = Command::new(helper)
            .stdin(Stdio::piped())
            .stdout(Stdio::piped())
            .stderr(Stdio::inherit())
            .spawn()
            .with_context(|| format!("starting persistent GPU helper {}", helper.display()))?;
        let stdin = child
            .stdin
            .take()
            .context("persistent GPU helper stdin was not piped")?;
        let stdout = child
            .stdout
            .take()
            .context("persistent GPU helper stdout was not piped")?;
        Ok(Self {
            child,
            stdin: BufWriter::new(stdin),
            stdout: BufReader::new(stdout),
        })
    }

    fn reduce<F>(&mut self, keys: &[u64], f: &mut F) -> Result<()>
    where
        F: FnMut(KmerKey, u64),
    {
        let count = keys.len() as u64;
        self.stdin.write_all(&count.to_le_bytes())?;
        for key in keys {
            self.stdin.write_all(&key.to_le_bytes())?;
        }
        self.stdin.flush()?;

        let mut unique_buf = [0u8; 8];
        self.stdout
            .read_exact(&mut unique_buf)
            .context("reading persistent GPU helper unique count")?;
        let unique = u64::from_le_bytes(unique_buf);
        let mut record = [0u8; 12];
        for _ in 0..unique {
            self.stdout
                .read_exact(&mut record)
                .context("reading persistent GPU helper reduced run")?;
            let key = u64::from_le_bytes(record[0..8].try_into().unwrap());
            let count = u32::from_le_bytes(record[8..12].try_into().unwrap());
            f(KmerKey::Short(key), u64::from(count));
        }
        Ok(())
    }

    fn finish(mut self) -> Result<()> {
        self.stdin.write_all(&u64::MAX.to_le_bytes())?;
        self.stdin.flush()?;
        drop(self.stdin);
        let status = self
            .child
            .wait()
            .context("waiting for persistent GPU helper")?;
        if !status.success() {
            bail!("persistent GPU helper failed with status {status}");
        }
        Ok(())
    }
}

fn replay_reduced_runs_file<F>(path: &Path, f: &mut F) -> Result<()>
where
    F: FnMut(KmerKey, u64),
{
    let mut reader =
        BufReader::new(fs::File::open(path).with_context(|| format!("open {}", path.display()))?);
    let mut record = [0u8; 12];
    loop {
        match reader.read_exact(&mut record) {
            Ok(()) => {
                let key = u64::from_le_bytes(record[0..8].try_into().unwrap());
                let count = u32::from_le_bytes(record[8..12].try_into().unwrap());
                f(KmerKey::Short(key), u64::from(count));
            }
            Err(err) if err.kind() == ErrorKind::UnexpectedEof => break,
            Err(err) => return Err(err).context("reading GPU reduced runs"),
        }
    }
    Ok(())
}

fn count_single_file(
    config: &Config,
    path: &Path,
    counts: &mut CountMap,
    limit: Option<u64>,
) -> Result<()> {
    let mut reader = open_sequence_reader(config, path, sequence_settings(config))?;
    let mut reads_seen = 0u64;
    let mut chunk = Vec::with_capacity(COUNT_PARALLEL_CHUNK_SIZE);
    while let Some(record) = reader.next_record()? {
        if limit_reached(limit, reads_seen) {
            break;
        }
        chunk.push((record, None));
        if chunk.len() >= COUNT_PARALLEL_CHUNK_SIZE {
            increment_counts_from_pair_chunk(config, counts, &chunk);
            chunk.clear();
        }
        reads_seen += 1;
    }
    if !chunk.is_empty() {
        increment_counts_from_pair_chunk(config, counts, &chunk);
    }
    Ok(())
}

fn count_single_file_sketch(
    config: &Config,
    path: &Path,
    sketch: &mut PackedCountMinSketch,
    limit: Option<u64>,
    prefilter: Option<PrefilterGate<'_>>,
) -> Result<()> {
    let mut reader = open_sequence_reader(config, path, sequence_settings(config))?;
    let mut reads_seen = 0u64;
    let mut chunk = Vec::with_capacity(COUNT_PARALLEL_CHUNK_SIZE);
    while let Some(record) = reader.next_record()? {
        if limit_reached(limit, reads_seen) {
            break;
        }
        chunk.push((record, None));
        if chunk.len() >= COUNT_PARALLEL_CHUNK_SIZE {
            increment_sketch_from_pair_chunk(config, sketch, &chunk, prefilter);
            chunk.clear();
        }
        reads_seen += 1;
    }
    if !chunk.is_empty() {
        increment_sketch_from_pair_chunk(config, sketch, &chunk, prefilter);
    }
    Ok(())
}

fn count_single_file_prefilter_sketch(
    config: &Config,
    path: &Path,
    sketch: &mut PrefilterCountMinSketch,
    limit: Option<u64>,
) -> Result<()> {
    match sketch {
        PrefilterCountMinSketch::Packed(sketch) => {
            count_single_file_sketch(config, path, sketch, limit, None)
        }
        PrefilterCountMinSketch::AtomicPacked(sketch) => {
            count_single_file_atomic_packed_sketch(config, path, sketch, limit)
        }
    }
}

fn count_single_file_atomic_packed_sketch(
    config: &Config,
    path: &Path,
    sketch: &AtomicPackedCountMinSketch,
    limit: Option<u64>,
) -> Result<()> {
    let mut reader = open_sequence_reader(config, path, sequence_settings(config))?;
    let mut reads_seen = 0u64;
    let mut chunk = Vec::with_capacity(COUNT_PARALLEL_CHUNK_SIZE);
    while let Some(record) = reader.next_record()? {
        if limit_reached(limit, reads_seen) {
            break;
        }
        chunk.push((record, None));
        if chunk.len() >= COUNT_PARALLEL_CHUNK_SIZE {
            increment_atomic_packed_sketch_from_pair_chunk(config, sketch, &chunk);
            chunk.clear();
        }
        reads_seen += 1;
    }
    if !chunk.is_empty() {
        increment_atomic_packed_sketch_from_pair_chunk(config, sketch, &chunk);
    }
    Ok(())
}

fn count_single_file_atomic_sketch(
    config: &Config,
    path: &Path,
    sketch: &AtomicCountMinSketch,
    limit: Option<u64>,
    prefilter: Option<PrefilterGate<'_>>,
) -> Result<()> {
    let mut reader = open_sequence_reader(config, path, sequence_settings(config))?;
    let mut reads_seen = 0u64;
    let mut chunk = Vec::with_capacity(COUNT_PARALLEL_CHUNK_SIZE);
    while let Some(record) = reader.next_record()? {
        if limit_reached(limit, reads_seen) {
            break;
        }
        chunk.push((record, None));
        if chunk.len() >= COUNT_PARALLEL_CHUNK_SIZE {
            increment_atomic_sketch_from_pair_chunk(config, sketch, &chunk, prefilter);
            chunk.clear();
        }
        reads_seen += 1;
    }
    if !chunk.is_empty() {
        increment_atomic_sketch_from_pair_chunk(config, sketch, &chunk, prefilter);
    }
    Ok(())
}

fn count_paired_files(
    config: &Config,
    path1: &Path,
    path2: &Path,
    counts: &mut CountMap,
    limit: Option<u64>,
) -> Result<()> {
    let settings = sequence_settings(config);
    let (mut reader1, mut reader2) = open_paired_sequence_readers(config, path1, path2, settings)?;
    if reader1.format() != reader2.format() {
        bail!("paired inputs must use the same FASTA/FASTQ format");
    }

    let mut pairs_seen = 0u64;
    let mut chunk = Vec::with_capacity(COUNT_PARALLEL_CHUNK_SIZE);
    loop {
        if limit_reached(limit, pairs_seen) {
            break;
        }
        match (reader1.next_record()?, reader2.next_record()?) {
            (None, None) => break,
            (Some(read1), Some(read2)) => {
                chunk.push((read1, Some(read2)));
                if chunk.len() >= COUNT_PARALLEL_CHUNK_SIZE {
                    increment_counts_from_pair_chunk(config, counts, &chunk);
                    chunk.clear();
                }
                pairs_seen += 1;
            }
            (Some(_), None) => bail!(
                "{} has fewer records than {}",
                path2.display(),
                path1.display()
            ),
            (None, Some(_)) => bail!(
                "{} has fewer records than {}",
                path1.display(),
                path2.display()
            ),
        }
    }
    if !chunk.is_empty() {
        increment_counts_from_pair_chunk(config, counts, &chunk);
    }
    Ok(())
}

fn count_paired_files_sketch(
    config: &Config,
    path1: &Path,
    path2: &Path,
    sketch: &mut PackedCountMinSketch,
    limit: Option<u64>,
    prefilter: Option<PrefilterGate<'_>>,
) -> Result<()> {
    let settings = sequence_settings(config);
    let (mut reader1, mut reader2) = open_paired_sequence_readers(config, path1, path2, settings)?;
    if reader1.format() != reader2.format() {
        bail!("paired inputs must use the same FASTA/FASTQ format");
    }

    let mut pairs_seen = 0u64;
    let mut chunk = Vec::with_capacity(COUNT_PARALLEL_CHUNK_SIZE);
    loop {
        if limit_reached(limit, pairs_seen) {
            break;
        }
        match (reader1.next_record()?, reader2.next_record()?) {
            (None, None) => break,
            (Some(read1), Some(read2)) => {
                chunk.push((read1, Some(read2)));
                if chunk.len() >= COUNT_PARALLEL_CHUNK_SIZE {
                    increment_sketch_from_pair_chunk(config, sketch, &chunk, prefilter);
                    chunk.clear();
                }
                pairs_seen += 1;
            }
            (Some(_), None) => bail!(
                "{} has fewer records than {}",
                path2.display(),
                path1.display()
            ),
            (None, Some(_)) => bail!(
                "{} has fewer records than {}",
                path1.display(),
                path2.display()
            ),
        }
    }
    if !chunk.is_empty() {
        increment_sketch_from_pair_chunk(config, sketch, &chunk, prefilter);
    }
    Ok(())
}

fn count_paired_files_atomic_packed_sketch(
    config: &Config,
    path1: &Path,
    path2: &Path,
    sketch: &AtomicPackedCountMinSketch,
    limit: Option<u64>,
) -> Result<()> {
    let settings = sequence_settings(config);
    let (mut reader1, mut reader2) = open_paired_sequence_readers(config, path1, path2, settings)?;
    if reader1.format() != reader2.format() {
        bail!("paired inputs must use the same FASTA/FASTQ format");
    }

    let mut pairs_seen = 0u64;
    let mut chunk = Vec::with_capacity(COUNT_PARALLEL_CHUNK_SIZE);
    loop {
        if limit_reached(limit, pairs_seen) {
            break;
        }
        match (reader1.next_record()?, reader2.next_record()?) {
            (None, None) => break,
            (Some(read1), Some(read2)) => {
                chunk.push((read1, Some(read2)));
                if chunk.len() >= COUNT_PARALLEL_CHUNK_SIZE {
                    increment_atomic_packed_sketch_from_pair_chunk(config, sketch, &chunk);
                    chunk.clear();
                }
                pairs_seen += 1;
            }
            (Some(_), None) => bail!(
                "{} has fewer records than {}",
                path2.display(),
                path1.display()
            ),
            (None, Some(_)) => bail!(
                "{} has fewer records than {}",
                path1.display(),
                path2.display()
            ),
        }
    }
    if !chunk.is_empty() {
        increment_atomic_packed_sketch_from_pair_chunk(config, sketch, &chunk);
    }
    Ok(())
}

fn count_paired_files_atomic_sketch(
    config: &Config,
    path1: &Path,
    path2: &Path,
    sketch: &AtomicCountMinSketch,
    limit: Option<u64>,
    prefilter: Option<PrefilterGate<'_>>,
) -> Result<()> {
    let settings = sequence_settings(config);
    let (mut reader1, mut reader2) = open_paired_sequence_readers(config, path1, path2, settings)?;
    if reader1.format() != reader2.format() {
        bail!("paired inputs must use the same FASTA/FASTQ format");
    }

    let mut pairs_seen = 0u64;
    let mut chunk = Vec::with_capacity(COUNT_PARALLEL_CHUNK_SIZE);
    loop {
        if limit_reached(limit, pairs_seen) {
            break;
        }
        match (reader1.next_record()?, reader2.next_record()?) {
            (None, None) => break,
            (Some(read1), Some(read2)) => {
                chunk.push((read1, Some(read2)));
                if chunk.len() >= COUNT_PARALLEL_CHUNK_SIZE {
                    increment_atomic_sketch_from_pair_chunk(config, sketch, &chunk, prefilter);
                    chunk.clear();
                }
                pairs_seen += 1;
            }
            (Some(_), None) => bail!(
                "{} has fewer records than {}",
                path2.display(),
                path1.display()
            ),
            (None, Some(_)) => bail!(
                "{} has fewer records than {}",
                path1.display(),
                path2.display()
            ),
        }
    }
    if !chunk.is_empty() {
        increment_atomic_sketch_from_pair_chunk(config, sketch, &chunk, prefilter);
    }
    Ok(())
}

fn normalize_primary(
    config: &Config,
    input_counts: &dyn CountLookup,
    mut output_counts: Option<&mut OutputCounts>,
    mut output_cardinality: Option<&mut KmerCardinalityEstimator>,
    cardinality_config: &Config,
    random_seed: u64,
) -> Result<RunSummary> {
    let mut readers = PrimaryReaders::open(config, config.max_reads)?;
    let format1 = readers.format1();
    let format2 = readers.format2();
    let mut writers = OptionalWriters::open(config, format1, format2)?;
    let mut summary = RunSummary::default();
    let mut rng = JavaXoshiro::new(random_seed);
    let mut chunk = Vec::with_capacity(NORMALIZE_PARALLEL_CHUNK_SIZE);

    while let Some((r1, r2)) = readers.next_pair()? {
        chunk.push((readers.input_list_index(), r1, r2, rng.next_double()));
        if chunk.len() >= NORMALIZE_PARALLEL_CHUNK_SIZE {
            let pairs = normalize_pair_chunk(config, input_counts, &chunk);
            write_normalized_pairs(
                config,
                &mut writers,
                &mut output_counts,
                &mut output_cardinality,
                cardinality_config,
                &mut summary,
                &pairs,
            )?;
            chunk.clear();
        }
    }
    if !chunk.is_empty() {
        let pairs = normalize_pair_chunk(config, input_counts, &chunk);
        write_normalized_pairs(
            config,
            &mut writers,
            &mut output_counts,
            &mut output_cardinality,
            cardinality_config,
            &mut summary,
            &pairs,
        )?;
    }

    writers.flush()?;
    Ok(summary)
}

fn normalize_pair_chunk(
    config: &Config,
    input_counts: &dyn CountLookup,
    pairs: &[NormalizationInput],
) -> Vec<NormalizedPair> {
    pairs
        .par_iter()
        .map(|(input_list_index, r1, r2, rand)| {
            let mut r1 = r1.clone();
            let mut r2 = r2.clone();
            if !config.trim_after_marking {
                trim_pair(config, &mut r1, r2.as_mut());
            }
            let decision = decide_pair(config, input_counts, &r1, r2.as_ref(), Some(*rand));
            let mut correction = CorrectionResult::default();
            if config.error_correct && !decision.toss {
                correction =
                    correct_pair_errors_with_rollback(config, input_counts, &mut r1, r2.as_mut());
            }
            if config.trim_after_marking && config.error_correct {
                trim_pair(config, &mut r1, r2.as_mut());
            }
            let (out_r1, out_r2) = maybe_rename_pair(config, &r1, r2.as_ref(), &decision.analysis);
            let read_count = 1 + u64::from(r2.is_some());
            let base_count = r1.len() as u64 + r2.as_ref().map(|r| r.len() as u64).unwrap_or(0);
            NormalizedPair {
                input_list_index: *input_list_index,
                r1,
                r2,
                out_r1,
                out_r2,
                decision,
                uncorrectable: correction.uncorrectable,
                read_count,
                base_count,
            }
        })
        .collect()
}

fn write_normalized_pairs(
    config: &Config,
    writers: &mut OptionalWriters,
    output_counts: &mut Option<&mut OutputCounts>,
    output_cardinality: &mut Option<&mut KmerCardinalityEstimator>,
    cardinality_config: &Config,
    summary: &mut RunSummary,
    pairs: &[NormalizedPair],
) -> Result<()> {
    for pair in pairs {
        writers.sync_to_input_list_index(config, pair.input_list_index)?;
        summary.reads_in += pair.read_count;
        summary.bases_in += pair.base_count;

        if pair.decision.toss {
            summary.reads_tossed += pair.read_count;
            summary.bases_tossed += pair.base_count;
        } else {
            summary.reads_kept += pair.read_count;
            summary.bases_kept += pair.base_count;
        }

        writers.write_pair(pair.decision.toss, &pair.out_r1, pair.out_r2.as_ref())?;
        if pair.uncorrectable {
            writers.write_uncorrected(&pair.r1, pair.r2.as_ref())?;
        }
        if depth_bin_outputs_enabled(config) {
            writers.write_depth_bin(
                config,
                &pair.decision.analysis,
                &pair.out_r1,
                pair.out_r2.as_ref(),
            )?;
        }
    }
    if let Some(counts) = output_counts.as_mut() {
        increment_output_counts_from_normalized_chunk(config, counts, pairs);
    }
    if let Some(estimator) = output_cardinality.as_mut() {
        for pair in pairs.iter().filter(|pair| !pair.decision.toss) {
            estimator.observe_pair(cardinality_config, &pair.r1, pair.r2.as_ref());
        }
    }
    Ok(())
}

fn increment_output_counts_from_normalized_chunk(
    config: &Config,
    counts: &mut OutputCounts,
    pairs: &[NormalizedPair],
) {
    match counts {
        OutputCounts::Exact(counts) => {
            let chunk_counts = pairs
                .par_iter()
                .filter(|pair| !pair.decision.toss)
                .fold(CountMap::default, |mut local_counts, pair| {
                    increment_pair_counts(config, &mut local_counts, &pair.r1, pair.r2.as_ref());
                    local_counts
                })
                .reduce(CountMap::default, |mut left, right| {
                    merge_count_maps(&mut left, right);
                    left
                });
            merge_count_maps(counts, chunk_counts);
        }
        OutputCounts::Sketch(sketch) => {
            increment_sketch_from_normalized_chunk(config, sketch, pairs);
        }
        OutputCounts::AtomicSketch(sketch) => {
            increment_atomic_sketch_from_normalized_chunk(config, sketch, pairs);
        }
    }
}

fn increment_atomic_sketch_from_normalized_chunk(
    config: &Config,
    sketch: &AtomicCountMinSketch,
    pairs: &[NormalizedPair],
) {
    if !config.deterministic {
        let (key_increments, newly_occupied) = pairs
            .par_iter()
            .filter(|pair| !pair.decision.toss)
            .map(|pair| {
                increment_pair_atomic_sketch_direct(
                    config,
                    sketch,
                    &pair.r1,
                    pair.r2.as_ref(),
                    None,
                )
            })
            .reduce(
                || (0u64, 0usize),
                |left, right| {
                    (
                        left.0.saturating_add(right.0),
                        left.1.saturating_add(right.1),
                    )
                },
            );
        sketch.add_key_increments(key_increments);
        sketch.add_occupied_slots(newly_occupied);
        return;
    }

    let chunk_counts = pairs
        .par_iter()
        .filter(|pair| !pair.decision.toss)
        .fold(CountMap::default, |mut local_counts, pair| {
            increment_pair_counts(config, &mut local_counts, &pair.r1, pair.r2.as_ref());
            local_counts
        })
        .reduce(CountMap::default, |mut left, right| {
            merge_count_maps(&mut left, right);
            left
        });
    let key_increments = chunk_counts.values().copied().sum();
    sketch.add_key_counts(&chunk_counts);
    sketch.add_key_increments(key_increments);
}

fn increment_sketch_from_normalized_chunk(
    config: &Config,
    sketch: &mut PackedCountMinSketch,
    pairs: &[NormalizedPair],
) {
    let chunk_counts = pairs
        .par_iter()
        .filter(|pair| !pair.decision.toss)
        .fold(CountMap::default, |mut local_counts, pair| {
            increment_pair_counts(config, &mut local_counts, &pair.r1, pair.r2.as_ref());
            local_counts
        })
        .reduce(CountMap::default, |mut left, right| {
            merge_count_maps(&mut left, right);
            left
        });
    let key_increments = chunk_counts.values().copied().sum();
    sketch.add_key_counts(&chunk_counts);
    sketch.add_key_increments(key_increments);
}

#[cfg(test)]
fn collect_primary_hist(
    config: &Config,
    hist_counts: &dyn CountLookup,
    keep_filter_counts: Option<&dyn CountLookup>,
    random_seed: u64,
) -> Result<Vec<u64>> {
    let mut readers = PrimaryReaders::open(config, config.max_reads)?;
    let mut hist = vec![0u64; config.hist_len];
    let mut rng = JavaXoshiro::new(random_seed);
    let mut chunk = Vec::with_capacity(HIST_PARALLEL_CHUNK_SIZE);

    while let Some((mut r1, mut r2)) = readers.next_pair()? {
        trim_pair(config, &mut r1, r2.as_mut());
        let rand = keep_filter_counts.map(|_| rng.next_double());
        chunk.push((r1, r2, rand));
        if chunk.len() >= HIST_PARALLEL_CHUNK_SIZE {
            increment_hist_from_pair_chunk(
                config,
                hist_counts,
                keep_filter_counts,
                &mut hist,
                &chunk,
            );
            chunk.clear();
        }
    }
    if !chunk.is_empty() {
        increment_hist_from_pair_chunk(config, hist_counts, keep_filter_counts, &mut hist, &chunk);
    }

    Ok(hist)
}

fn collect_primary_sparse_hist(
    config: &Config,
    hist_counts: &dyn CountLookup,
    keep_filter_counts: Option<&dyn CountLookup>,
    random_seed: u64,
) -> Result<SparseHist> {
    let mut readers = PrimaryReaders::open(config, config.max_reads)?;
    let mut hist = SparseHist::default();
    let mut rng = JavaXoshiro::new(random_seed);
    let mut chunk = Vec::with_capacity(HIST_PARALLEL_CHUNK_SIZE);

    while let Some((mut r1, mut r2)) = readers.next_pair()? {
        trim_pair(config, &mut r1, r2.as_mut());
        let rand = keep_filter_counts.map(|_| rng.next_double());
        chunk.push((r1, r2, rand));
        if chunk.len() >= HIST_PARALLEL_CHUNK_SIZE {
            let chunk_hist =
                sparse_hist_from_pair_chunk(config, hist_counts, keep_filter_counts, &chunk);
            merge_sparse_hist(&mut hist, chunk_hist);
            chunk.clear();
        }
    }
    if !chunk.is_empty() {
        let chunk_hist =
            sparse_hist_from_pair_chunk(config, hist_counts, keep_filter_counts, &chunk);
        merge_sparse_hist(&mut hist, chunk_hist);
    }

    Ok(hist)
}

#[cfg(test)]
fn collect_primary_read_hist(
    config: &Config,
    hist_counts: &dyn CountLookup,
    keep_filter_counts: Option<&dyn CountLookup>,
    random_seed: u64,
) -> Result<ReadDepthHistogram> {
    let mut readers = PrimaryReaders::open(config, config.max_reads)?;
    let mut hist = ReadDepthHistogram::new(config.hist_len);
    let mut rng = JavaXoshiro::new(random_seed);
    let mut chunk = Vec::with_capacity(HIST_PARALLEL_CHUNK_SIZE);

    while let Some((mut r1, mut r2)) = readers.next_pair()? {
        trim_pair(config, &mut r1, r2.as_mut());
        let rand = keep_filter_counts.map(|_| rng.next_double());
        chunk.push((r1, r2, rand));
        if chunk.len() >= HIST_PARALLEL_CHUNK_SIZE {
            increment_read_hist_from_pair_chunk(
                config,
                hist_counts,
                keep_filter_counts,
                &mut hist,
                &chunk,
            );
            chunk.clear();
        }
    }
    if !chunk.is_empty() {
        increment_read_hist_from_pair_chunk(
            config,
            hist_counts,
            keep_filter_counts,
            &mut hist,
            &chunk,
        );
    }

    Ok(hist)
}

fn collect_primary_sparse_read_hist(
    config: &Config,
    hist_counts: &dyn CountLookup,
    keep_filter_counts: Option<&dyn CountLookup>,
    random_seed: u64,
) -> Result<SparseReadDepthHist> {
    let mut readers = PrimaryReaders::open(config, config.max_reads)?;
    let mut hist = SparseReadDepthHist::default();
    let mut rng = JavaXoshiro::new(random_seed);
    let mut chunk = Vec::with_capacity(HIST_PARALLEL_CHUNK_SIZE);

    while let Some((mut r1, mut r2)) = readers.next_pair()? {
        trim_pair(config, &mut r1, r2.as_mut());
        let rand = keep_filter_counts.map(|_| rng.next_double());
        chunk.push((r1, r2, rand));
        if chunk.len() >= HIST_PARALLEL_CHUNK_SIZE {
            let chunk_hist =
                sparse_read_hist_from_pair_chunk(config, hist_counts, keep_filter_counts, &chunk);
            merge_sparse_read_depth_hist(&mut hist, chunk_hist);
            chunk.clear();
        }
    }
    if !chunk.is_empty() {
        let chunk_hist =
            sparse_read_hist_from_pair_chunk(config, hist_counts, keep_filter_counts, &chunk);
        merge_sparse_read_depth_hist(&mut hist, chunk_hist);
    }

    Ok(hist)
}

#[cfg(test)]
fn collect_primary_hist_and_read_hist(
    config: &Config,
    hist_counts: &dyn CountLookup,
    keep_filter_counts: Option<&dyn CountLookup>,
    random_seed: u64,
) -> Result<(Vec<u64>, ReadDepthHistogram)> {
    let mut readers = PrimaryReaders::open(config, config.max_reads)?;
    let mut depth_hist = vec![0u64; config.hist_len];
    let mut read_hist = ReadDepthHistogram::new(config.hist_len);
    let mut rng = JavaXoshiro::new(random_seed);
    let mut chunk = Vec::with_capacity(HIST_PARALLEL_CHUNK_SIZE);

    while let Some((mut r1, mut r2)) = readers.next_pair()? {
        trim_pair(config, &mut r1, r2.as_mut());
        let rand = keep_filter_counts.map(|_| rng.next_double());
        chunk.push((r1, r2, rand));
        if chunk.len() >= HIST_PARALLEL_CHUNK_SIZE {
            increment_hist_and_read_hist_from_pair_chunk(
                config,
                hist_counts,
                keep_filter_counts,
                &mut depth_hist,
                &mut read_hist,
                &chunk,
            );
            chunk.clear();
        }
    }
    if !chunk.is_empty() {
        increment_hist_and_read_hist_from_pair_chunk(
            config,
            hist_counts,
            keep_filter_counts,
            &mut depth_hist,
            &mut read_hist,
            &chunk,
        );
    }

    Ok((depth_hist, read_hist))
}

fn collect_primary_sparse_hist_and_read_hist(
    config: &Config,
    hist_counts: &dyn CountLookup,
    keep_filter_counts: Option<&dyn CountLookup>,
    random_seed: u64,
) -> Result<(SparseHist, SparseReadDepthHist)> {
    let mut readers = PrimaryReaders::open(config, config.max_reads)?;
    let mut depth_hist = SparseHist::default();
    let mut read_hist = SparseReadDepthHist::default();
    let mut rng = JavaXoshiro::new(random_seed);
    let mut chunk = Vec::with_capacity(HIST_PARALLEL_CHUNK_SIZE);

    while let Some((mut r1, mut r2)) = readers.next_pair()? {
        trim_pair(config, &mut r1, r2.as_mut());
        let rand = keep_filter_counts.map(|_| rng.next_double());
        chunk.push((r1, r2, rand));
        if chunk.len() >= HIST_PARALLEL_CHUNK_SIZE {
            let (chunk_depth_hist, chunk_read_hist) = sparse_hist_and_read_hist_from_pair_chunk(
                config,
                hist_counts,
                keep_filter_counts,
                &chunk,
            );
            merge_sparse_hist(&mut depth_hist, chunk_depth_hist);
            merge_sparse_read_depth_hist(&mut read_hist, chunk_read_hist);
            chunk.clear();
        }
    }
    if !chunk.is_empty() {
        let (chunk_depth_hist, chunk_read_hist) = sparse_hist_and_read_hist_from_pair_chunk(
            config,
            hist_counts,
            keep_filter_counts,
            &chunk,
        );
        merge_sparse_hist(&mut depth_hist, chunk_depth_hist);
        merge_sparse_read_depth_hist(&mut read_hist, chunk_read_hist);
    }

    Ok((depth_hist, read_hist))
}

fn emit_read_local_side_outputs(config: &Config) -> Result<()> {
    if !read_local_side_outputs_enabled(config) {
        return Ok(());
    }

    let mut hist = collect_read_local_side_hists(config)?;
    if let Some(quality) = hist.quality.take() {
        emit_quality_side_outputs(config, &quality)?;
    }
    if let (Some(path), Some(length)) = (&config.length_hist_out, hist.length.as_ref()) {
        write_length_hist(path, length, config)?;
    }
    if let (Some(path), Some(gc)) = (&config.gc_hist_out, hist.gc.as_ref()) {
        write_gc_hist(path, gc, config)?;
    }
    if let (Some(path), Some(base)) = (&config.base_hist_out, hist.base.as_ref()) {
        write_base_content_hist(path, base, config)?;
    }
    if let (Some(path), Some(entropy)) = (&config.entropy_hist_out, hist.entropy.as_ref()) {
        write_entropy_hist(path, entropy, config)?;
    }
    if let (Some(path), Some(identity)) = (&config.identity_hist_out, hist.identity.as_ref()) {
        write_identity_hist(path, identity, config)?;
    }
    if let Some(alignment) = hist.alignment.as_ref() {
        emit_alignment_fallback_side_outputs(config, alignment)?;
    }
    if let (Some(path), Some(barcodes)) = (&config.barcode_stats_out, hist.barcodes.as_ref()) {
        write_barcode_stats(path, barcodes, config)?;
    }
    Ok(())
}

fn read_local_side_outputs_enabled(config: &Config) -> bool {
    config.quality_hist_out.is_some()
        || config.base_quality_hist_out.is_some()
        || config.quality_count_hist_out.is_some()
        || config.average_quality_hist_out.is_some()
        || config.overall_base_quality_hist_out.is_some()
        || config.length_hist_out.is_some()
        || config.gc_hist_out.is_some()
        || config.base_hist_out.is_some()
        || config.entropy_hist_out.is_some()
        || config.identity_hist_out.is_some()
        || config.barcode_stats_out.is_some()
        || alignment_fallback_side_outputs_enabled(config)
}

fn quality_side_outputs_enabled(config: &Config) -> bool {
    config.quality_hist_out.is_some()
        || config.base_quality_hist_out.is_some()
        || config.quality_count_hist_out.is_some()
        || config.average_quality_hist_out.is_some()
        || config.overall_base_quality_hist_out.is_some()
}

fn alignment_fallback_side_outputs_enabled(config: &Config) -> bool {
    config.match_hist_out.is_some()
        || config.insert_hist_out.is_some()
        || config.quality_accuracy_hist_out.is_some()
        || config.indel_hist_out.is_some()
        || config.error_hist_out.is_some()
}

fn emit_quality_side_outputs(config: &Config, hist: &QualitySideHistograms) -> Result<()> {
    if let Some(path) = &config.quality_hist_out {
        write_quality_hist(path, &hist.overall, config)?;
    }
    if let Some(path) = &config.quality_count_hist_out {
        write_quality_count_hist(
            path,
            &hist.first_counts,
            &hist.second_counts,
            hist.paired,
            config,
        )?;
    }
    if let Some(path) = &config.average_quality_hist_out {
        write_average_quality_hist(path, &hist.first_avg, &hist.second_avg, hist.paired, config)?;
    }
    if let Some(path) = &config.overall_base_quality_hist_out {
        write_overall_base_quality_hist(path, &hist.overall, config)?;
    }
    if let Some(path) = &config.base_quality_hist_out {
        write_base_quality_hist(path, hist, config)?;
    }
    Ok(())
}

fn collect_read_local_side_hists(config: &Config) -> Result<ReadLocalSideHistograms> {
    let mut readers = PrimaryReaders::open(config, config.max_reads)?;
    let quality_len = side_hist_len(config);
    let side_len = side_hist_len(config);
    let mut hist = ReadLocalSideHistograms {
        quality: quality_side_outputs_enabled(config).then(|| QualitySideHistograms {
            overall: vec![0; quality_len],
            first_counts: vec![0; quality_len],
            second_counts: vec![0; quality_len],
            first_avg: vec![0; quality_len],
            second_avg: vec![0; quality_len],
            first_by_pos: vec![vec![0; quality_len]; side_len],
            second_by_pos: vec![vec![0; quality_len]; side_len],
            paired: false,
        }),
        length: config
            .length_hist_out
            .is_some()
            .then(|| ReadDepthHistogram::new(side_len)),
        gc: config
            .gc_hist_out
            .is_some()
            .then(|| ReadDepthHistogram::new(gc_hist_len(config))),
        base: config
            .base_hist_out
            .is_some()
            .then(|| BaseContentHistogram {
                first: vec![BaseCounts::default(); side_len],
                second: vec![BaseCounts::default(); side_len],
            }),
        entropy: config
            .entropy_hist_out
            .is_some()
            .then(|| vec![0u64; config.entropy_bins.saturating_add(1).max(1)]),
        identity: config
            .identity_hist_out
            .is_some()
            .then(|| ReadDepthHistogram::new(config.identity_bins.saturating_add(1).max(1))),
        alignment: alignment_fallback_side_outputs_enabled(config).then(|| {
            AlignmentFallbackHistograms {
                first_match: vec![MatchCounts::default(); side_len],
                second_match: vec![MatchCounts::default(); side_len],
                quality_match: vec![0; quality_len],
                ..AlignmentFallbackHistograms::default()
            }
        }),
        barcodes: config.barcode_stats_out.is_some().then(BTreeMap::new),
    };

    while let Some((mut r1, mut r2)) = readers.next_pair()? {
        trim_pair(config, &mut r1, r2.as_mut());
        if let Some(barcodes) = hist.barcodes.as_mut() {
            increment_barcode_stats(barcodes, &r1, r2.is_some());
        }
        increment_read_local_side_hists(config, &mut hist, &r1, false);
        if let Some(mate) = r2.as_ref() {
            increment_read_local_side_hists(config, &mut hist, mate, true);
        }
    }

    Ok(hist)
}

fn side_hist_len(config: &Config) -> usize {
    config.side_hist_len.unwrap_or(config.hist_len).max(1)
}

fn gc_hist_len(config: &Config) -> usize {
    config.gc_bins.unwrap_or(101).max(1)
}

fn increment_length_hist(hist: &mut ReadDepthHistogram, read_len: usize) {
    let idx = read_len.min(hist.reads.len().saturating_sub(1));
    hist.reads[idx] += 1;
    hist.bases[idx] += read_len as u64;
}

fn increment_read_local_side_hists(
    config: &Config,
    hist: &mut ReadLocalSideHistograms,
    record: &SequenceRecord,
    second: bool,
) {
    if let Some(quality) = hist.quality.as_mut() {
        if second {
            quality.paired = true;
        }
        increment_quality_side_hists(config, quality, record, second);
    }
    if let Some(length) = hist.length.as_mut() {
        increment_length_hist(length, record.len());
    }
    if let Some(gc) = hist.gc.as_mut() {
        increment_gc_hist(gc, record);
    }
    if let Some(base) = hist.base.as_mut() {
        if second {
            increment_base_content_hist(&mut base.second, record);
        } else {
            increment_base_content_hist(&mut base.first, record);
        }
    }
    if let Some(entropy) = hist.entropy.as_mut() {
        increment_entropy_hist(config, entropy, record);
    }
    if let Some(identity) = hist.identity.as_mut() {
        increment_sequence_identity_hist(identity, record);
    }
    if let Some(alignment) = hist.alignment.as_mut() {
        increment_alignment_fallback_hists(config, alignment, record, second);
    }
}

fn increment_gc_hist(hist: &mut ReadDepthHistogram, record: &SequenceRecord) {
    let mut gc = 0usize;
    let mut acgt = 0usize;
    for base in &record.bases {
        match *base {
            b'G' | b'C' | b'g' | b'c' => {
                gc += 1;
                acgt += 1;
            }
            b'A' | b'T' | b'U' | b'a' | b't' | b'u' => acgt += 1,
            _ => {}
        }
    }
    let idx = if acgt == 0 {
        0
    } else {
        ((gc * hist.reads.len()) / acgt).min(hist.reads.len().saturating_sub(1))
    };
    hist.reads[idx] += 1;
    hist.bases[idx] += record.len() as u64;
}

fn increment_quality_side_hists(
    config: &Config,
    hist: &mut QualitySideHistograms,
    record: &SequenceRecord,
    second: bool,
) {
    if record.is_empty() {
        return;
    }

    let quality_len = hist.overall.len();
    let last_quality_idx = quality_len.saturating_sub(1);
    let (counts, avg_counts, by_pos) = if second {
        (
            &mut hist.second_counts,
            &mut hist.second_avg,
            &mut hist.second_by_pos,
        )
    } else {
        (
            &mut hist.first_counts,
            &mut hist.first_avg,
            &mut hist.first_by_pos,
        )
    };

    let mut sum = 0usize;
    for idx in 0..record.len() {
        let quality = record_quality_at(config, record, idx).min(last_quality_idx);
        hist.overall[quality] += 1;
        counts[quality] += 1;
        sum += quality;
        if idx < by_pos.len() {
            by_pos[idx][quality] += 1;
        }
    }

    let avg = ((sum as f64) / (record.len() as f64)).round() as usize;
    avg_counts[avg.min(last_quality_idx)] += 1;
}

fn record_quality_at(config: &Config, record: &SequenceRecord, idx: usize) -> usize {
    record
        .qualities
        .as_ref()
        .and_then(|qualities| qualities.get(idx))
        .map_or(config.fake_quality as usize, |quality| {
            quality.saturating_sub(33) as usize
        })
}

fn increment_base_content_hist(hist: &mut [BaseCounts], record: &SequenceRecord) {
    for (idx, base) in record.bases.iter().copied().enumerate().take(hist.len()) {
        let counts = &mut hist[idx];
        match base {
            b'A' | b'a' => counts.a += 1,
            b'C' | b'c' => counts.c += 1,
            b'G' | b'g' => counts.g += 1,
            b'T' | b't' | b'U' | b'u' => counts.t += 1,
            _ => counts.n += 1,
        }
    }
}

fn increment_entropy_hist(config: &Config, hist: &mut [u64], record: &SequenceRecord) {
    if record.is_empty() {
        return;
    }
    if let Some(entropy) = read_entropy(config, &record.bases) {
        let bins = hist.len().saturating_sub(1);
        let idx = ((entropy * hist.len() as f64) as usize).min(bins);
        hist[idx] += 1;
    }
}

fn increment_sequence_identity_hist(hist: &mut ReadDepthHistogram, record: &SequenceRecord) {
    let idx = hist.reads.len().saturating_sub(1);
    hist.reads[idx] += 1;
    hist.bases[idx] += record.len() as u64;
}

fn increment_barcode_stats(
    barcodes: &mut BTreeMap<String, u64>,
    record: &SequenceRecord,
    paired: bool,
) {
    let barcode = header_to_barcode(&record.id).unwrap_or("NONE");
    let count = if paired { 2 } else { 1 };
    *barcodes.entry(barcode.to_string()).or_insert(0) += count;
}

fn header_to_barcode(id: &str) -> Option<&str> {
    let loc = id.rfind(':')?;
    let loc2 = id
        .find(' ')
        .map(|idx| idx as isize)
        .unwrap_or(-1)
        .max(id.find('/').map(|idx| idx as isize).unwrap_or(-1));
    if (loc as isize) <= loc2 || loc >= id.len().saturating_sub(1) {
        return None;
    }
    let start = loc + 1;
    let stop = id[start..]
        .find([' ', '\t'])
        .map_or(id.len(), |offset| start + offset);
    Some(&id[start..stop])
}

fn increment_alignment_fallback_hists(
    config: &Config,
    hist: &mut AlignmentFallbackHistograms,
    record: &SequenceRecord,
    second: bool,
) {
    hist.read_count += 1;
    hist.base_count += record.len() as u64;
    if second {
        hist.paired = true;
        hist.pair_count += 1;
    }

    let match_hist = if second {
        &mut hist.second_match
    } else {
        &mut hist.first_match
    };
    for (idx, base) in record
        .bases
        .iter()
        .copied()
        .enumerate()
        .take(match_hist.len())
    {
        if is_acgt(base) {
            match_hist[idx].matches += 1;
        } else {
            match_hist[idx].n += 1;
        }
    }

    for idx in 0..record.len() {
        let quality =
            record_quality_at(config, record, idx).min(hist.quality_match.len().saturating_sub(1));
        hist.quality_match[quality] += 1;
    }
}

fn read_entropy(config: &Config, bases: &[u8]) -> Option<f64> {
    let k = config.entropy_k.clamp(1, 15);
    if bases.len() < k {
        return base_entropy(config, bases);
    }

    let window = config.entropy_window.max(k).min(bases.len());
    let mut sum = 0.0;
    let mut count = 0usize;
    for start in 0..=bases.len() - window {
        if let Some(entropy) = window_kmer_entropy(config, &bases[start..start + window], k) {
            sum += entropy;
            count += 1;
        }
    }

    if count == 0 {
        None
    } else {
        Some((sum / count as f64).clamp(0.0, 1.0))
    }
}

fn window_kmer_entropy(config: &Config, window: &[u8], k: usize) -> Option<f64> {
    if window.len() < k {
        return base_entropy(config, window);
    }

    let mut counts: FxHashMap<Vec<u8>, u64> = FxHashMap::default();
    let mut total = 0u64;
    for kmer in window.windows(k) {
        if !config.allow_entropy_ns && kmer.iter().any(|base| !is_acgt(*base)) {
            continue;
        }
        let key: Vec<u8> = kmer
            .iter()
            .copied()
            .map(|base| match base.to_ascii_uppercase() {
                b'A' | b'C' | b'G' | b'T' => base.to_ascii_uppercase(),
                _ => b'N',
            })
            .collect();
        *counts.entry(key).or_insert(0) += 1;
        total += 1;
    }

    if total == 0 {
        return None;
    }
    let entropy = shannon_entropy(counts.values().copied(), total);
    let max_entropy = (total as f64).ln();
    Some(if max_entropy > 0.0 {
        entropy / max_entropy
    } else {
        0.0
    })
}

fn base_entropy(config: &Config, bases: &[u8]) -> Option<f64> {
    let mut counts = [0u64; 5];
    let mut total = 0u64;
    for base in bases {
        let idx = match base.to_ascii_uppercase() {
            b'A' => Some(0),
            b'C' => Some(1),
            b'G' => Some(2),
            b'T' | b'U' => Some(3),
            _ if config.allow_entropy_ns => Some(4),
            _ => None,
        };
        if let Some(idx) = idx {
            counts[idx] += 1;
            total += 1;
        }
    }
    if total == 0 {
        return None;
    }
    let entropy = shannon_entropy(counts, total);
    let nonzero = counts.into_iter().filter(|count| *count > 0).count();
    let max_entropy = (nonzero.max(1) as f64).ln();
    Some(if max_entropy > 0.0 {
        entropy / max_entropy
    } else {
        0.0
    })
}

fn shannon_entropy(counts: impl IntoIterator<Item = u64>, total: u64) -> f64 {
    let total = total as f64;
    counts
        .into_iter()
        .filter(|count| *count > 0)
        .map(|count| {
            let p = count as f64 / total;
            -p * p.ln()
        })
        .sum()
}

fn is_acgt(base: u8) -> bool {
    matches!(base, b'A' | b'C' | b'G' | b'T' | b'a' | b'c' | b'g' | b't')
}

fn analyze_pair(
    config: &Config,
    counts: &dyn CountLookup,
    r1: &SequenceRecord,
    r2: Option<&SequenceRecord>,
) -> PairAnalysis {
    let (read1, read2) = match r2 {
        Some(record) if r1.len() + record.len() >= PAIRED_ANALYSIS_JOIN_MIN_BASES => {
            let (read1, read2) = rayon::join(
                || analyze_read(config, counts, r1),
                || analyze_read(config, counts, record),
            );
            (read1, Some(read2))
        }
        Some(record) => (
            analyze_read(config, counts, r1),
            Some(analyze_read(config, counts, record)),
        ),
        None => (analyze_read(config, counts, r1), None),
    };
    pair_analysis_from_reads(config, read1, read2)
}

fn analyze_pair_for_two_configs(
    config: &Config,
    other_config: &Config,
    counts: &dyn CountLookup,
    r1: &SequenceRecord,
    r2: Option<&SequenceRecord>,
) -> (PairAnalysis, PairAnalysis) {
    if !can_share_read_coverage(config, other_config) {
        return (
            analyze_pair(config, counts, r1, r2),
            analyze_pair(other_config, counts, r1, r2),
        );
    }

    let ((read1, other_read1), read2_pair) = match r2 {
        Some(record) if r1.len() + record.len() >= PAIRED_ANALYSIS_JOIN_MIN_BASES => {
            let (first, second) = rayon::join(
                || analyze_read_for_two_configs(config, other_config, counts, r1),
                || analyze_read_for_two_configs(config, other_config, counts, record),
            );
            (first, Some(second))
        }
        Some(record) => (
            analyze_read_for_two_configs(config, other_config, counts, r1),
            Some(analyze_read_for_two_configs(
                config,
                other_config,
                counts,
                record,
            )),
        ),
        None => (
            analyze_read_for_two_configs(config, other_config, counts, r1),
            None,
        ),
    };
    let (read2, other_read2) = read2_pair
        .map(|(read, other_read)| (Some(read), Some(other_read)))
        .unwrap_or((None, None));
    (
        pair_analysis_from_reads(config, read1, read2),
        pair_analysis_from_reads(other_config, other_read1, other_read2),
    )
}

fn can_share_read_coverage(config: &Config, other_config: &Config) -> bool {
    config.k == other_config.k
        && (config.canonical || config.k <= 31) == (other_config.canonical || other_config.k <= 31)
        && config.fix_spikes == other_config.fix_spikes
}

fn pair_analysis_from_reads(
    config: &Config,
    read1: ReadAnalysis,
    read2: Option<ReadAnalysis>,
) -> PairAnalysis {
    let depth_proxy_al = match (&read2, config.use_lower_depth) {
        (Some(read2), true) => min_option(read1.depth_al, read2.depth_al),
        (Some(read2), false) => max_option(read1.depth_al, read2.depth_al),
        (None, _) => read1.depth_al,
    };
    let max_true_depth = match &read2 {
        Some(read2) => max_option(read1.true_depth, read2.true_depth),
        None => read1.true_depth,
    };
    let low_kmer_count =
        read1.low_kmer_count + read2.as_ref().map(|read| read.low_kmer_count).unwrap_or(0);
    let total_kmer_count = read1.total_kmer_count
        + read2
            .as_ref()
            .map(|read| read.total_kmer_count)
            .unwrap_or(0);
    PairAnalysis {
        error1: read1.error,
        error2: read2.as_ref().is_some_and(|read| read.error),
        read1,
        read2,
        depth_proxy_al,
        max_true_depth,
        low_kmer_count,
        total_kmer_count,
    }
}

fn analyze_read(
    config: &Config,
    counts: &dyn CountLookup,
    record: &SequenceRecord,
) -> ReadAnalysis {
    let coverage = read_coverage_desc(config, counts, record);
    analyze_read_from_coverage(config, coverage.coverage_desc, coverage.had_kmer_windows)
}

fn analyze_read_for_two_configs(
    config: &Config,
    other_config: &Config,
    counts: &dyn CountLookup,
    record: &SequenceRecord,
) -> (ReadAnalysis, ReadAnalysis) {
    let coverage = read_coverage_desc(config, counts, record);
    let other_coverage = coverage.coverage_desc.clone();
    (
        analyze_read_from_coverage(config, coverage.coverage_desc, coverage.had_kmer_windows),
        analyze_read_from_coverage(other_config, other_coverage, coverage.had_kmer_windows),
    )
}

struct ReadCoverageDesc {
    coverage_desc: Vec<i64>,
    had_kmer_windows: bool,
}

fn read_coverage_desc(
    config: &Config,
    counts: &dyn CountLookup,
    record: &SequenceRecord,
) -> ReadCoverageDesc {
    let windows = unfiltered_kmer_windows_for_record(record, config);
    let mut coverage: Vec<i64> = windows
        .iter()
        .map(|window| match window {
            Some(kmer) => u64_to_i64_saturating(counts.depth(kmer)),
            None => -1,
        })
        .collect();
    if coverage.is_empty() {
        return ReadCoverageDesc {
            coverage_desc: coverage,
            had_kmer_windows: record.len() >= config.k,
        };
    }
    if config.fix_spikes {
        fix_spikes(&mut coverage, &windows, counts, config.k);
    }
    if coverage.len() >= COVERAGE_PAR_SORT_MIN_WINDOWS {
        coverage.par_sort_unstable_by(|a, b| b.cmp(a));
    } else {
        coverage.sort_unstable_by(|a, b| b.cmp(a));
    }
    ReadCoverageDesc {
        coverage_desc: coverage,
        had_kmer_windows: true,
    }
}

fn analyze_read_from_coverage(
    config: &Config,
    coverage: Vec<i64>,
    had_kmer_windows: bool,
) -> ReadAnalysis {
    if coverage.is_empty() {
        return ReadAnalysis {
            had_kmer_windows,
            ..ReadAnalysis::default()
        };
    }
    let cov_last = coverage.len() - 1;
    let high = coverage[percentile_index(cov_last, config.high_percentile)];
    let low = coverage[percentile_index(cov_last, config.low_percentile)];
    let true_depth = coverage[percentile_index(cov_last, config.depth_percentile)];
    let min_true_depth = low;
    let min_depth = u64_to_i64_saturating(config.min_depth)
        .max(high / u64_to_i64_saturating(config.error_detect_ratio));

    let mut above_limit = cov_last as isize;
    while above_limit >= 0 && coverage[above_limit as usize] < min_depth {
        above_limit -= 1;
    }

    let depth_al = if above_limit >= 0
        && ((above_limit as usize + 1) >= config.min_kmers_over_min_depth
            || config.min_kmers_over_min_depth > coverage.len())
    {
        let idx = ((above_limit as f64) * (1.0 - config.depth_percentile)) as usize;
        non_negative_depth(coverage[idx])
    } else {
        None
    };

    let low_thresh = u64_to_i64_saturating(config.low_thresh);
    let high_thresh = u64_to_i64_saturating(config.high_thresh);
    let error_detect_ratio = u64_to_i64_saturating(config.error_detect_ratio);
    let error = high <= low_thresh
        || (high >= high_thresh && low <= low_thresh)
        || high >= low.saturating_mul(error_detect_ratio);
    let low_kmer_count =
        low_kmer_count(&coverage, low_thresh, high_thresh, high, error_detect_ratio);

    ReadAnalysis {
        depth_al,
        true_depth: non_negative_depth(true_depth),
        min_true_depth: non_negative_depth(min_true_depth),
        low_kmer_count,
        total_kmer_count: coverage.len(),
        error,
        had_kmer_windows: true,
        coverage_desc: coverage,
    }
}

fn low_kmer_count(
    coverage_desc: &[i64],
    low_thresh: i64,
    high_thresh: i64,
    high_depth: i64,
    error_detect_ratio: i64,
) -> usize {
    if coverage_desc.is_empty() {
        return 0;
    }
    if coverage_desc[0] <= low_thresh {
        return coverage_desc.len();
    }
    if high_depth < high_thresh {
        return 0;
    }
    let limit = low_thresh.min(high_depth / error_detect_ratio.max(1));
    coverage_desc
        .iter()
        .rev()
        .take_while(|&&depth| depth <= limit)
        .count()
}

fn correct_pair_errors(
    config: &Config,
    counts: &dyn CountLookup,
    r1: &mut SequenceRecord,
    r2: Option<&mut SequenceRecord>,
) -> CorrectionResult {
    let mut result = CorrectionResult::default();
    let mut r2 = r2;
    if config.overlap_error_correct
        && !config.mark_errors_only
        && let Some(mate) = r2.as_deref_mut()
    {
        let overlap = correct_pair_by_overlap(config, r1, mate);
        result.corrected += overlap.corrected;
        result.marked += overlap.marked;
        result.uncorrectable |= overlap.uncorrectable;
    }

    let read_result = correct_read_errors(config, counts, r1);
    result.corrected += read_result.corrected;
    result.marked += read_result.marked;
    result.uncorrectable |= read_result.uncorrectable;
    if let Some(mate) = r2 {
        let mate_result = correct_read_errors(config, counts, mate);
        result.corrected += mate_result.corrected;
        result.marked += mate_result.marked;
        result.uncorrectable |= mate_result.uncorrectable;
    }
    result
}

fn correct_pair_errors_with_rollback(
    config: &Config,
    counts: &dyn CountLookup,
    r1: &mut SequenceRecord,
    mut r2: Option<&mut SequenceRecord>,
) -> CorrectionResult {
    let rollback =
        (!config.mark_uncorrectable_errors).then(|| (r1.clone(), r2.as_deref().cloned()));
    let correction = correct_pair_errors(config, counts, r1, r2.as_deref_mut());
    if correction.uncorrectable
        && let Some((original_r1, original_r2)) = rollback
    {
        *r1 = original_r1;
        if let (Some(mate), Some(original)) = (r2, original_r2) {
            *mate = original;
        }
    }
    correction
}

fn correct_pair_by_overlap(
    config: &Config,
    r1: &mut SequenceRecord,
    r2: &mut SequenceRecord,
) -> CorrectionResult {
    let Some(overlap) = best_pair_overlap(r1, r2) else {
        return CorrectionResult::default();
    };
    if overlap_expected_mismatch_rejects(r1, r2, &overlap) {
        return CorrectionResult::default();
    }
    if overlap_probability_rejects(r1, r2, &overlap) {
        return CorrectionResult::default();
    }
    let mut corrected = 0usize;

    for pair in overlap.pairs {
        let b1 = r1.bases[pair.r1_index].to_ascii_uppercase();
        let b2 = complement_base(r2.bases[pair.r2_index]).to_ascii_uppercase();
        let q1 = base_quality(r1, pair.r1_index);
        let q2 = base_quality(r2, pair.r2_index);
        let Some((merged_base, merged_quality)) =
            overlap_consensus_base_and_quality(config, b1, b2, q1, q2)
        else {
            continue;
        };

        let merged_r2_base = complement_base(merged_base);
        if r1.bases[pair.r1_index] != merged_base || r2.bases[pair.r2_index] != merged_r2_base {
            corrected += 1;
        }

        r1.bases[pair.r1_index] = merged_base;
        r2.bases[pair.r2_index] = merged_r2_base;

        if config.change_quality
            && let (Some(r1_qualities), Some(r2_qualities)) =
                (r1.qualities.as_mut(), r2.qualities.as_mut())
        {
            let merged_ascii = merged_quality.saturating_add(33);
            r1_qualities[pair.r1_index] = merged_ascii;
            r2_qualities[pair.r2_index] = merged_ascii;
        }
    }

    CorrectionResult {
        corrected,
        ..CorrectionResult::default()
    }
}

fn overlap_expected_mismatch_rejects(
    r1: &SequenceRecord,
    r2: &SequenceRecord,
    overlap: &PairOverlap,
) -> bool {
    let (Some(q1), Some(q2)) = (r1.qualities.as_ref(), r2.qualities.as_ref()) else {
        return false;
    };

    let mut expected = 0.0f64;
    for pair in &overlap.pairs {
        let b1 = r1.bases[pair.r1_index].to_ascii_uppercase();
        let b2 = complement_base(r2.bases[pair.r2_index]).to_ascii_uppercase();
        if !is_defined_base(b1) || !is_defined_base(b2) {
            continue;
        }
        let p1 = 1.0 - phred_error_probability(q1[pair.r1_index].saturating_sub(33));
        let p2 = 1.0 - phred_error_probability(q2[pair.r2_index].saturating_sub(33));
        expected += 1.0 - (p1 * p2);
    }

    (expected + 0.05) * 4.0 < overlap.mismatches as f64
}

fn overlap_probability_rejects(
    r1: &SequenceRecord,
    r2: &SequenceRecord,
    overlap: &PairOverlap,
) -> bool {
    const MIN_PROBABILITY: f64 = 0.0008;

    let (Some(q1), Some(q2)) = (r1.qualities.as_ref(), r2.qualities.as_ref()) else {
        return false;
    };

    let mut ln_actual = 0.0f64;
    let mut ln_common = 0.0f64;
    let mut measured = 0usize;

    for pair in &overlap.pairs {
        let b1 = r1.bases[pair.r1_index].to_ascii_uppercase();
        let b2 = complement_base(r2.bases[pair.r2_index]).to_ascii_uppercase();
        if !is_defined_base(b1) || !is_defined_base(b2) {
            continue;
        }

        let prob_correct = overlap_correctness_probability_v4(q1[pair.r1_index])
            * overlap_correctness_probability_v4(q2[pair.r2_index]);
        let prob_match = prob_correct + (1.0 - prob_correct) * 0.25;
        let prob_error = 1.0 - prob_match;

        ln_common += prob_match.max(prob_error).ln();
        ln_actual += if b1 == b2 { prob_match } else { prob_error }.ln();
        measured += 1;
    }

    if measured == 0 {
        return false;
    }

    0.5 * (ln_actual - ln_common) < MIN_PROBABILITY.ln()
}

fn overlap_consensus_base_and_quality(
    config: &Config,
    r1_base: u8,
    r2_base: u8,
    q1: u8,
    q2: u8,
) -> Option<(u8, u8)> {
    const MAX_MERGE_QUALITY: u8 = 50;

    if !is_defined_base(r1_base) && !is_defined_base(r2_base) {
        return None;
    }
    if !is_defined_base(r1_base) {
        return Some((r2_base, q2));
    }
    if !is_defined_base(r2_base) {
        return Some((r1_base, q1));
    }
    if r1_base == r2_base {
        let merged_quality = q1
            .max(q2)
            .saturating_add(q1.min(q2) / 4)
            .min(MAX_MERGE_QUALITY);
        return Some((r1_base, merged_quality));
    }
    if q1 == q2 {
        return Some((b'N', 0));
    }
    if q1 > q2 {
        if q2 > config.max_quality_to_correct {
            return None;
        }
        return Some((r1_base, q1.saturating_sub(q2)));
    }
    if q1 > config.max_quality_to_correct {
        return None;
    }
    Some((r2_base, q2.saturating_sub(q1)))
}

fn overlap_entropy_min_overlap(bases: &[u8]) -> usize {
    overlap_entropy_min_overlap_side(bases.iter().copied()).max(overlap_entropy_min_overlap_side(
        bases.iter().rev().copied(),
    ))
}

fn overlap_entropy_min_overlap_side(bases: impl IntoIterator<Item = u8>) -> usize {
    const K: usize = 3;
    const MASK: usize = (1 << (2 * K)) - 1;
    const MIN_SCORE: usize = 42;

    let mut counts = [0u16; 1 << (2 * K)];
    let mut kmer = 0usize;
    let mut len = 0usize;
    let mut ones = 0usize;
    let mut twos = 0usize;
    let mut seen = 0usize;

    for base in bases {
        let Some(bits) = base_to_two_bit(base) else {
            len = 0;
            kmer = 0;
            seen += 1;
            continue;
        };
        len += 1;
        kmer = ((kmer << 2) | bits) & MASK;
        if len >= K {
            counts[kmer] = counts[kmer].saturating_add(1);
            if counts[kmer] == 1 {
                ones += 1;
            } else if counts[kmer] == 2 {
                twos += 1;
            }
            if ones * 4 + twos >= MIN_SCORE {
                return seen;
            }
        }
        seen += 1;
    }

    seen + 1
}

fn base_to_two_bit(base: u8) -> Option<usize> {
    match base.to_ascii_uppercase() {
        b'A' => Some(0),
        b'C' => Some(1),
        b'G' => Some(2),
        b'T' => Some(3),
        _ => None,
    }
}

#[derive(Debug, Clone, Copy)]
struct OverlapBasePair {
    r1_index: usize,
    r2_index: usize,
}

#[derive(Debug, Clone)]
struct PairOverlap {
    pairs: Vec<OverlapBasePair>,
    mismatches: usize,
}

const OVERLAP_MAX_RATIO: f64 = 0.075;
const OVERLAP_MIN_SECOND_RATIO: f64 = 0.12;
const OVERLAP_RATIO_MARGIN: f64 = 7.5;
const OVERLAP_RATIO_OFFSET: f64 = 0.55;
const OVERLAP_PROB_CORRECT4: &[f64] = &[
    0.0000, 0.2501, 0.3690, 0.4988, 0.6019, 0.6838, 0.7488, 0.8005, 0.8415, 0.8741, 0.9000, 0.9206,
    0.9369, 0.9499, 0.9602, 0.9684, 0.9749, 0.9800, 0.9842, 0.9874, 0.9900, 0.9921, 0.9937, 0.9950,
    0.9960, 0.9968, 0.9975, 0.9980, 0.9984, 0.9987, 0.9990, 0.9992, 0.9994, 0.9995, 0.9996, 0.9997,
    0.9997, 0.9998, 0.9998, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,
    0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,
];

fn best_pair_overlap(r1: &SequenceRecord, r2: &SequenceRecord) -> Option<PairOverlap> {
    best_pair_overlap_without_qualities(&r1.bases, &r2.bases)
}

fn overlap_correctness_probability_v4(quality_ascii: u8) -> f64 {
    let phred = quality_ascii.saturating_sub(33) as usize;
    OVERLAP_PROB_CORRECT4[phred.min(OVERLAP_PROB_CORRECT4.len() - 1)]
}

fn best_pair_overlap_without_qualities(r1: &[u8], r2: &[u8]) -> Option<PairOverlap> {
    if r1.is_empty() || r2.is_empty() {
        return None;
    }
    if r1.len().min(r2.len()) < 35 {
        return None;
    }
    let min_overlap = 11usize
        .max(overlap_entropy_min_overlap(r1))
        .max(overlap_entropy_min_overlap(r2));
    let min_length = r1.len().min(r2.len());
    if min_overlap > min_length {
        return None;
    }

    let best_ratio_cap = find_best_overlap_ratio_without_qualities(r1, r2, min_overlap);
    if best_ratio_cap > OVERLAP_MAX_RATIO {
        return None;
    }
    let max_ratio = best_ratio_cap.min(OVERLAP_MAX_RATIO);
    let margin2 = (OVERLAP_RATIO_MARGIN + OVERLAP_RATIO_OFFSET) / min_length as f64;
    let mut best_insert: Option<usize> = None;
    let mut best_overlap = 0usize;
    let mut best_bad = min_length as f64;
    let mut best_good = 0.0f64;
    let mut best_ratio = 1.0f64;
    let mut best_mismatches = 0usize;
    let mut second_best_ratio = 1.0f64;
    let mut ambig = false;

    let largest_insert_to_test = r1.len() + r2.len() - 5;
    for insert in (25..=largest_insert_to_test).rev() {
        let start1 = if insert <= r2.len() {
            0
        } else {
            insert - r2.len()
        };
        let start2 = if insert >= r2.len() {
            0
        } else {
            r2.len() - insert
        };
        let overlap = (r1.len() - start1).min(r2.len() - start2).min(insert);
        if overlap < 5 {
            continue;
        }

        let bad_limit =
            1.2 * best_ratio.min(max_ratio) * OVERLAP_RATIO_MARGIN * overlap as f64 + 1.0;
        let mut good = 0.0f64;
        let mut bad = 0.0f64;
        let mut mismatches = 0usize;

        for step in 0..overlap {
            let r1_index = start1 + step;
            let r2_rc_index = start2 + step;
            let r2_index = r2.len() - 1 - r2_rc_index;
            let b1 = r1[r1_index].to_ascii_uppercase();
            let b2 = complement_base(r2[r2_index]).to_ascii_uppercase();
            if b1 == b2 {
                if b1 != b'N' {
                    good += 0.95;
                }
            } else {
                bad += 0.95;
                mismatches += 1;
                if bad > bad_limit {
                    break;
                }
            }
        }

        if bad > bad_limit {
            continue;
        }
        if bad == 0.0 && good > 5.0 && good < min_overlap as f64 {
            return None;
        }

        let ratio = (bad + OVERLAP_RATIO_OFFSET) / overlap as f64;
        if ratio < best_ratio * OVERLAP_RATIO_MARGIN {
            ambig = ratio * OVERLAP_RATIO_MARGIN >= best_ratio || good < min_overlap as f64;

            if ratio < best_ratio {
                second_best_ratio = best_ratio;
                best_insert = Some(insert);
                best_overlap = overlap;
                best_bad = bad;
                best_good = good;
                best_ratio = ratio;
                best_mismatches = mismatches;
            } else if ratio < second_best_ratio {
                second_best_ratio = ratio;
            }

            if (ambig && best_ratio < margin2) || second_best_ratio < OVERLAP_MIN_SECOND_RATIO {
                return None;
            }
        }
    }

    if second_best_ratio < OVERLAP_MIN_SECOND_RATIO {
        ambig = true;
    }
    if !ambig && best_ratio > max_ratio {
        return None;
    }

    let insert = best_insert?;
    let start1 = if insert <= r2.len() {
        0
    } else {
        insert - r2.len()
    };
    let start2 = if insert >= r2.len() {
        0
    } else {
        r2.len() - insert
    };
    let mut pairs = Vec::with_capacity(best_overlap);
    for step in 0..best_overlap {
        let r1_index = start1 + step;
        let r2_rc_index = start2 + step;
        let r2_index = r2.len() - 1 - r2_rc_index;
        pairs.push(OverlapBasePair { r1_index, r2_index });
    }

    let _ = (best_bad, best_good);
    Some(PairOverlap {
        pairs,
        mismatches: best_mismatches,
    })
}

fn find_best_overlap_ratio_without_qualities(r1: &[u8], r2: &[u8], min_overlap: usize) -> f64 {
    let mut best_ratio = OVERLAP_MAX_RATIO + 0.0001;
    let largest_insert_to_test = r1.len() + r2.len() - min_overlap;

    for insert in (35..=largest_insert_to_test).rev() {
        let start1 = if insert <= r2.len() {
            0
        } else {
            insert - r2.len()
        };
        let start2 = if insert >= r2.len() {
            0
        } else {
            r2.len() - insert
        };
        let overlap = (r1.len() - start1).min(r2.len() - start2).min(insert);
        if overlap < min_overlap {
            continue;
        }

        let mut good = 0.0f64;
        let mut bad = 0.0f64;
        let bad_limit = best_ratio * overlap as f64 + 1.0;

        for step in 0..overlap {
            let r1_index = start1 + step;
            let r2_rc_index = start2 + step;
            let r2_index = r2.len() - 1 - r2_rc_index;
            let b1 = r1[r1_index].to_ascii_uppercase();
            let b2 = complement_base(r2[r2_index]).to_ascii_uppercase();
            if b1 == b2 {
                if b1 != b'N' {
                    good += 0.95;
                }
            } else {
                bad += 0.95;
                if bad > bad_limit {
                    break;
                }
            }
        }

        if bad > bad_limit {
            continue;
        }
        if bad == 0.0 && good > 5.0 && good < min_overlap as f64 {
            return 100.0;
        }
        let ratio = (bad + OVERLAP_RATIO_OFFSET) / overlap as f64;
        if ratio < best_ratio {
            best_ratio = ratio;
            if good >= min_overlap as f64 && ratio < OVERLAP_MAX_RATIO * 0.5 {
                return best_ratio;
            }
        }
    }

    best_ratio
}

fn correct_read_errors(
    config: &Config,
    counts: &dyn CountLookup,
    record: &mut SequenceRecord,
) -> CorrectionResult {
    if config.max_errors_to_correct == 0 || record.len() < config.k || config.k > 31 {
        return CorrectionResult::default();
    }
    let mut coverage = coverage_windows_for_record(config, counts, record);
    if coverage.len() <= config.prefix_len.max(1) {
        return CorrectionResult::default();
    }
    if !has_error_discontinuity(config, &coverage) {
        return CorrectionResult::default();
    }

    if config.mark_errors_only {
        return mark_read_errors(config, record, &coverage);
    }

    let original_bases = record.bases.clone();
    let original_qualities = record.qualities.clone();
    let mut result = CorrectionResult::default();
    let mut remaining = config.max_errors_to_correct;

    if config.correct_from_left {
        let left = correct_errors_from_left(config, counts, record, &mut coverage, remaining);
        if left.uncorrectable {
            record.bases = original_bases;
            record.qualities = original_qualities;
            if config.mark_uncorrectable_errors {
                result.marked += mark_read_errors(config, record, &coverage).marked;
            }
            result.uncorrectable = true;
            return result;
        }
        remaining = remaining.saturating_sub(left.corrected);
        result.corrected += left.corrected;
    }

    if config.correct_from_right && remaining > 0 {
        let checkpoint_bases = record.bases.clone();
        let checkpoint_qualities = record.qualities.clone();
        let right = correct_errors_from_right(config, counts, record, &mut coverage, remaining);
        if right.uncorrectable {
            record.bases = checkpoint_bases;
            record.qualities = checkpoint_qualities;
            if config.mark_uncorrectable_errors {
                result.marked += mark_read_errors(config, record, &coverage).marked;
            }
            result.uncorrectable = true;
            return result;
        }
        result.corrected += right.corrected;
    }

    result
}

fn correct_errors_from_left(
    config: &Config,
    counts: &dyn CountLookup,
    record: &mut SequenceRecord,
    coverage: &mut Vec<i64>,
    max_to_correct: usize,
) -> CorrectionResult {
    let mut found = 0usize;
    let mut corrected = 0usize;
    let low = u64_to_i64_saturating(config.error_correct_low_thresh);
    let high = u64_to_i64_saturating(config.error_correct_high_thresh);
    let mult = u64_to_i64_saturating(config.error_correct_ratio);

    for i in config.prefix_len..coverage.len() {
        let a = min_coverage(&coverage[i - config.prefix_len..i]);
        let b = coverage[i];
        if !is_correction_discontinuity(a, b, low, high, mult) {
            continue;
        }
        found += 1;
        let loc = i + config.k - 1;
        if found > max_to_correct || base_quality(record, loc) > config.max_quality_to_correct {
            return CorrectionResult {
                corrected,
                uncorrectable: true,
                ..CorrectionResult::default()
            };
        }
        let target_lower = high.max(a / 2);
        let target_upper = a.saturating_mul(2);
        let target = CorrectionTarget {
            low,
            lower_bound: target_lower,
            upper_bound: target_upper,
            mult,
        };
        if try_correct_base(config, counts, record, loc, target) {
            corrected += 1;
            *coverage = coverage_windows_for_record(config, counts, record);
        } else {
            return CorrectionResult {
                corrected,
                uncorrectable: true,
                ..CorrectionResult::default()
            };
        }
    }

    CorrectionResult {
        corrected,
        ..CorrectionResult::default()
    }
}

fn correct_errors_from_right(
    config: &Config,
    counts: &dyn CountLookup,
    record: &mut SequenceRecord,
    coverage: &mut Vec<i64>,
    max_to_correct: usize,
) -> CorrectionResult {
    if coverage.len() <= config.prefix_len {
        return CorrectionResult::default();
    }
    let mut found = 0usize;
    let mut corrected = 0usize;
    let low = u64_to_i64_saturating(config.error_correct_low_thresh);
    let high = u64_to_i64_saturating(config.error_correct_high_thresh);
    let mult = u64_to_i64_saturating(config.error_correct_ratio);
    let start = coverage.len() - config.prefix_len - 1;

    for i in (0..=start).rev() {
        let a = min_coverage(&coverage[i + 1..=i + config.prefix_len]);
        let b = coverage[i];
        if !is_correction_discontinuity(a, b, low, high, mult) {
            continue;
        }
        found += 1;
        let loc = i;
        if found > max_to_correct || base_quality(record, loc) > config.max_quality_to_correct {
            return CorrectionResult {
                corrected,
                uncorrectable: true,
                ..CorrectionResult::default()
            };
        }
        let target_lower = high.max(a / 2);
        let target_upper = a.saturating_mul(2);
        let target = CorrectionTarget {
            low,
            lower_bound: target_lower,
            upper_bound: target_upper,
            mult,
        };
        if try_correct_base(config, counts, record, loc, target) {
            corrected += 1;
            *coverage = coverage_windows_for_record(config, counts, record);
        } else {
            return CorrectionResult {
                corrected,
                uncorrectable: true,
                ..CorrectionResult::default()
            };
        }
    }

    CorrectionResult {
        corrected,
        ..CorrectionResult::default()
    }
}

fn try_correct_base(
    config: &Config,
    counts: &dyn CountLookup,
    record: &mut SequenceRecord,
    loc: usize,
    target: CorrectionTarget,
) -> bool {
    let original = record.bases[loc];
    let mut candidates = [(b'A', 0i64), (b'C', 0), (b'G', 0), (b'T', 0)];
    for (base, support) in &mut candidates {
        *support = substitution_support(config, counts, record, loc, *base);
    }
    candidates.sort_by(|left, right| right.1.cmp(&left.1));
    let (best_base, best_support) = candidates[0];
    let second_best = candidates[1].1;
    if best_base == original.to_ascii_uppercase() {
        return false;
    }
    if best_support < target.lower_bound || best_support > target.upper_bound {
        return false;
    }
    if !(second_best <= target.low || second_best.saturating_mul(target.mult) <= best_support) {
        return false;
    }

    record.bases[loc] = best_base;
    if !is_defined_base(original)
        && let Some(qualities) = record.qualities.as_mut()
    {
        qualities[loc] = 20u8.saturating_add(33);
    }
    true
}

fn substitution_support(
    config: &Config,
    counts: &dyn CountLookup,
    record: &SequenceRecord,
    loc: usize,
    base: u8,
) -> i64 {
    let mut candidate = record.clone();
    candidate.bases[loc] = base;
    let windows = unfiltered_kmer_windows_for_record(&candidate, config);
    if windows.is_empty() {
        return 0;
    }
    let first = (loc + 1).saturating_sub(config.k);
    let last = loc.min(windows.len() - 1);
    let mut support = i64::MAX;
    for window in windows.iter().take(last + 1).skip(first) {
        let depth = window
            .as_ref()
            .map(|kmer| u64_to_i64_saturating(counts.depth(kmer)))
            .unwrap_or(0);
        support = support.min(depth);
    }
    if support == i64::MAX { 0 } else { support }
}

fn mark_read_errors(
    config: &Config,
    record: &mut SequenceRecord,
    coverage: &[i64],
) -> CorrectionResult {
    let low = u64_to_i64_saturating(config.error_correct_low_thresh);
    let high = u64_to_i64_saturating(config.error_correct_high_thresh);
    let mult = u64_to_i64_saturating(config.error_correct_ratio);
    let mut marked = 0usize;
    let mut marks = Vec::new();

    if config.correct_from_left {
        for i in config.prefix_len..coverage.len() {
            let a = min_coverage(&coverage[i - config.prefix_len..i]);
            let b = coverage[i];
            if is_correction_discontinuity(a, b, low, high, mult) {
                marks.push(i + config.k - 1);
            }
        }
    }
    if config.correct_from_right && coverage.len() > config.prefix_len {
        let start = coverage.len() - config.prefix_len - 1;
        for i in (0..=start).rev() {
            let a = min_coverage(&coverage[i + 1..=i + config.prefix_len]);
            let b = coverage[i];
            if is_correction_discontinuity(a, b, low, high, mult) {
                marks.push(i);
            }
        }
    }

    marks.sort_unstable();
    marks.dedup();
    for loc in marks {
        if let Some(qualities) = record.qualities.as_mut() {
            let phred = qualities[loc].saturating_sub(33);
            if phred == 0 {
                continue;
            }
            let new_phred = if config.mark_with_one {
                1
            } else {
                (phred / 2).saturating_sub(3).max(1)
            };
            qualities[loc] = new_phred.saturating_add(33);
        } else {
            record.bases[loc] = b'N';
        }
        marked += 1;
    }

    CorrectionResult {
        marked,
        ..CorrectionResult::default()
    }
}

fn coverage_windows_for_record(
    config: &Config,
    counts: &dyn CountLookup,
    record: &SequenceRecord,
) -> Vec<i64> {
    unfiltered_kmer_windows_for_record(record, config)
        .iter()
        .map(|window| {
            window
                .as_ref()
                .map(|kmer| u64_to_i64_saturating(counts.depth(kmer)))
                .unwrap_or(0)
        })
        .collect()
}

fn has_error_discontinuity(config: &Config, coverage: &[i64]) -> bool {
    let low = u64_to_i64_saturating(config.error_correct_low_thresh);
    let high = u64_to_i64_saturating(config.error_correct_high_thresh);
    let mult = u64_to_i64_saturating(config.error_correct_ratio);
    if coverage.len() <= config.prefix_len {
        return false;
    }
    for i in config.prefix_len..coverage.len() {
        if is_correction_discontinuity(
            min_coverage(&coverage[i - config.prefix_len..i]),
            coverage[i],
            low,
            high,
            mult,
        ) {
            return true;
        }
    }
    let start = coverage.len() - config.prefix_len - 1;
    for i in (0..=start).rev() {
        if is_correction_discontinuity(
            min_coverage(&coverage[i + 1..=i + config.prefix_len]),
            coverage[i],
            low,
            high,
            mult,
        ) {
            return true;
        }
    }
    false
}

fn is_correction_discontinuity(a: i64, b: i64, low: i64, high: i64, mult: i64) -> bool {
    a >= high && (b <= low || a >= b.saturating_mul(mult))
}

fn min_coverage(values: &[i64]) -> i64 {
    values.iter().copied().min().unwrap_or(0)
}

fn base_quality(record: &SequenceRecord, loc: usize) -> u8 {
    record
        .qualities
        .as_ref()
        .and_then(|qualities| qualities.get(loc))
        .copied()
        .map(|quality| quality.saturating_sub(33))
        .unwrap_or(10)
}

fn is_defined_base(base: u8) -> bool {
    matches!(base.to_ascii_uppercase(), b'A' | b'C' | b'G' | b'T')
}

fn complement_base(base: u8) -> u8 {
    match base.to_ascii_uppercase() {
        b'A' => b'T',
        b'C' => b'G',
        b'G' => b'C',
        b'T' => b'A',
        _ => b'N',
    }
}

fn fix_spikes(
    coverage: &mut [i64],
    windows: &[Option<KmerKey>],
    counts: &dyn CountLookup,
    k: usize,
) {
    if k == 0 || coverage.len() < 3 {
        return;
    }
    if coverage[1] - coverage[0] > 1 {
        coverage[0] = precise_kmer_count(windows[0].as_ref(), counts, k);
    }

    let last = coverage.len() - 1;
    if coverage[last] - coverage[last - 1] > 1 {
        coverage[last] = precise_kmer_count(windows[last].as_ref(), counts, k);
    }

    for i in 1..last {
        let b = coverage[i];
        if b <= 1 {
            continue;
        }
        let a = coverage[i - 1].max(1);
        let c = coverage[i + 1].max(1);
        if b > a && b > c && (b < 6 || b > a + 1 || b > c + 1) {
            coverage[i] = precise_min_kmer_count(windows[i].as_ref(), counts, k);
        }
    }
}

fn precise_kmer_count(window: Option<&KmerKey>, counts: &dyn CountLookup, k: usize) -> i64 {
    let Some(window) = window else {
        return 0;
    };
    let key = raw_kmer_key(window);
    let b = kmer_count(window, key, counts, k);
    if b < 1 {
        return b;
    }
    let a = left_kmer_count(window, key, counts, k);
    if a >= b {
        return b;
    }
    let c = right_kmer_count(window, key, counts, k);
    if c >= b {
        return b;
    }
    (a + c) / 2
}

fn precise_min_kmer_count(window: Option<&KmerKey>, counts: &dyn CountLookup, k: usize) -> i64 {
    let Some(window) = window else {
        return 0;
    };
    let key = raw_kmer_key(window);
    let b = kmer_count(window, key, counts, k);
    if b < 1 {
        return b;
    }
    let a = left_kmer_count(window, key, counts, k);
    if a < 1 {
        return a;
    }
    let c = right_kmer_count(window, key, counts, k);
    a.min(b).min(c)
}

fn raw_kmer_key(window: &KmerKey) -> u64 {
    match window {
        KmerKey::Short(key) | KmerKey::LongHash(key) => *key,
    }
}

fn kmer_count(template: &KmerKey, raw_key: u64, counts: &dyn CountLookup, k: usize) -> i64 {
    let key = match template {
        KmerKey::Short(_) => KmerKey::Short(canonical_short_code(raw_key, k)),
        KmerKey::LongHash(_) => KmerKey::LongHash(java_canonical_long_key(raw_key, k)),
    };
    u64_to_i64_saturating(counts.depth(&key))
}

fn left_kmer_count(template: &KmerKey, key: u64, counts: &dyn CountLookup, k: usize) -> i64 {
    let key2 = key >> 2;
    let shift = ((2 * (k - 1)) & 63) as u32;
    (0..4)
        .map(|base| kmer_count(template, key2 | (base << shift), counts, k))
        .fold(0i64, i64::saturating_add)
}

fn right_kmer_count(template: &KmerKey, key: u64, counts: &dyn CountLookup, k: usize) -> i64 {
    let mask = if k >= 32 {
        u64::MAX
    } else {
        (1u64 << (2 * k)) - 1
    };
    let key2 = (key << 2) & mask;
    (0..4)
        .map(|base| kmer_count(template, key2 | base, counts, k))
        .fold(0i64, i64::saturating_add)
}

fn java_canonical_long_key(key: u64, k: usize) -> u64 {
    let reverse = java_reverse_complement_binary_fast(key, k);
    key.max(reverse)
}

fn java_reverse_complement_binary_fast(key: u64, k: usize) -> u64 {
    let mut x = !key;
    x = ((x & 0x3333_3333_3333_3333) << 2) | ((x & 0xCCCC_CCCC_CCCC_CCCC) >> 2);
    x = ((x & 0x0F0F_0F0F_0F0F_0F0F) << 4) | ((x & 0xF0F0_F0F0_F0F0_F0F0) >> 4);
    x = ((x & 0x00FF_00FF_00FF_00FF) << 8) | ((x & 0xFF00_FF00_FF00_FF00) >> 8);
    x = ((x & 0x0000_FFFF_0000_FFFF) << 16) | ((x & 0xFFFF_0000_FFFF_0000) >> 16);
    x = x.rotate_right(32);

    let shift = (2usize.wrapping_mul(32usize.wrapping_sub(k)) & 63) as u32;
    x >> shift
}

fn decide_pair(
    config: &Config,
    input_counts: &dyn CountLookup,
    r1: &SequenceRecord,
    r2: Option<&SequenceRecord>,
    rand: Option<f64>,
) -> PairDecision {
    let analysis = analyze_pair(config, input_counts, r1, r2);
    decide_pair_from_analysis(config, r1, r2, analysis, rand)
}

fn decide_pair_from_analysis(
    config: &Config,
    r1: &SequenceRecord,
    r2: Option<&SequenceRecord>,
    analysis: PairAnalysis,
    rand: Option<f64>,
) -> PairDecision {
    let (target_depth, max_depth) = dynamic_depth_limits(config, &analysis);
    let mut toss = false;

    match analysis.depth_proxy_al {
        Some(depth) => {
            if depth > max_depth && (analysis.error1 || analysis.error2 || !config.discard_bad_only)
            {
                let coin = deterministic_coin(rand, depth);
                if coin > target_depth {
                    toss = true;
                }
            }
        }
        None => toss = true,
    }

    if r1.len() < config.min_length || r2.is_some_and(|mate| mate.len() < config.min_length) {
        toss = true;
    }

    if config.toss_error_reads && (analysis.error1 || analysis.error2) {
        let save_rare = config.save_rare_reads
            && analysis
                .depth_proxy_al
                .is_some_and(|depth| depth <= target_depth && depth >= config.high_thresh);
        if !save_rare
            && (!config.require_both_bad || r2.is_none() || (analysis.error1 && analysis.error2))
        {
            toss = true;
        }
    }

    if config.toss_by_low_true_depth && !config.save_rare_reads {
        let low_enough = analysis
            .max_true_depth
            .is_some_and(|depth| depth < config.min_depth);
        let required_bad = !config.require_both_bad
            || r2.is_none()
            || (depth_below_min(analysis.read1.min_true_depth, config.min_depth)
                && analysis
                    .read2
                    .as_ref()
                    .is_some_and(|read| depth_below_min(read.min_true_depth, config.min_depth)));
        if low_enough && required_bad {
            toss = true;
        }
    }

    if config.keep_all {
        toss = false;
    }

    PairDecision { toss, analysis }
}

fn dynamic_depth_limits(config: &Config, analysis: &PairAnalysis) -> (u64, u64) {
    let default_max_depth = config.max_depth.unwrap_or(config.target_depth);
    if analysis.low_kmer_count == 0 || analysis.total_kmer_count == 0 {
        return (config.target_depth, default_max_depth);
    }

    let low_target = ((config.target_depth as f64) * config.target_bad_percent_low)
        .round()
        .max(1.0);
    let high_target = ((config.target_depth as f64) * config.target_bad_percent_high)
        .round()
        .max(low_target)
        .min(config.target_depth as f64);
    let fraction_good = (analysis.total_kmer_count - analysis.low_kmer_count) as f64
        / analysis.total_kmer_count as f64;
    let adjusted = low_target + (high_target - low_target) * (fraction_good * fraction_good);
    let target = adjusted as u64;
    (target.max(1), target.max(1))
}

fn maybe_rename_pair(
    config: &Config,
    r1: &SequenceRecord,
    r2: Option<&SequenceRecord>,
    analysis: &PairAnalysis,
) -> (SequenceRecord, Option<SequenceRecord>) {
    if !config.rename_reads {
        return (r1.clone(), r2.cloned());
    }
    let d1 = depth_label(analysis.read1.depth_al);
    let out1 = match r2 {
        Some(_) => {
            let mut id = format!(
                "id={},d1={},d2={}",
                r1.numeric_id,
                d1,
                depth_label(analysis.read2.as_ref().and_then(|a| a.depth_al))
            );
            if config.error_correct {
                id.push_str(",e1=0,e2=0");
            }
            id.push_str(" /1");
            r1.renamed(id)
        }
        None => {
            let mut id = format!("id={},d1={}", r1.numeric_id, d1);
            if config.error_correct {
                id.push_str(",e1=0");
            }
            r1.renamed(id)
        }
    };
    let out2 = r2.map(|mate| {
        let mut id = format!(
            "id={},d1={},d2={}",
            r1.numeric_id,
            d1,
            depth_label(analysis.read2.as_ref().and_then(|a| a.depth_al))
        );
        if config.error_correct {
            id.push_str(",e1=0,e2=0");
        }
        id.push_str(" /2");
        mate.renamed(id)
    });
    (out1, out2)
}

fn depth_label(depth: Option<u64>) -> String {
    depth
        .map(|value| value.to_string())
        .unwrap_or_else(|| "-1".to_string())
}

fn increment_pair_counts(
    config: &Config,
    counts: &mut CountMap,
    r1: &SequenceRecord,
    r2: Option<&SequenceRecord>,
) {
    increment_pair_counts_with_prefilter(config, counts, r1, r2, None);
}

fn increment_pair_counts_with_prefilter(
    config: &Config,
    counts: &mut CountMap,
    r1: &SequenceRecord,
    r2: Option<&SequenceRecord>,
    prefilter: Option<PrefilterGate<'_>>,
) {
    if config.remove_duplicate_kmers && config.k <= 31 {
        for kmer in unique_pair_kmers(config, r1, r2) {
            if prefilter.is_none_or(|gate| gate.should_count_in_main(&kmer)) {
                *counts.entry(kmer).or_insert(0) += 1;
            }
        }
    } else {
        for_each_kmer_for_record(r1, config, |kmer| {
            if prefilter.is_none_or(|gate| gate.should_count_in_main(&kmer)) {
                *counts.entry(kmer).or_insert(0) += 1;
            }
        });
        if let Some(mate) = r2 {
            for_each_kmer_for_record(mate, config, |kmer| {
                if prefilter.is_none_or(|gate| gate.should_count_in_main(&kmer)) {
                    *counts.entry(kmer).or_insert(0) += 1;
                }
            });
        }
    }
}

fn increment_counts_from_pair_chunk(
    config: &Config,
    counts: &mut CountMap,
    pairs: &[(SequenceRecord, Option<SequenceRecord>)],
) {
    let chunk_counts = pairs
        .par_iter()
        .fold(
            || count_chunk_local_map(config, pairs),
            |mut local_counts, (r1, r2)| {
                increment_pair_counts(config, &mut local_counts, r1, r2.as_ref());
                local_counts
            },
        )
        .reduce(CountMap::default, |mut left, right| {
            merge_count_maps(&mut left, right);
            left
        });
    merge_count_maps(counts, chunk_counts);
}

fn increment_sketch_from_pair_chunk(
    config: &Config,
    sketch: &mut PackedCountMinSketch,
    pairs: &[(SequenceRecord, Option<SequenceRecord>)],
    prefilter: Option<PrefilterGate<'_>>,
) {
    if config.deterministic && sketch.update_mode == CountMinUpdateMode::Conservative {
        increment_sketch_from_pair_chunk_sorted_replay(config, sketch, pairs, prefilter);
        return;
    }
    let chunk_counts = pairs
        .par_iter()
        .fold(
            || count_chunk_local_map(config, pairs),
            |mut local_counts, (r1, r2)| {
                increment_pair_counts_with_prefilter(
                    config,
                    &mut local_counts,
                    r1,
                    r2.as_ref(),
                    prefilter,
                );
                local_counts
            },
        )
        .reduce(CountMap::default, |mut left, right| {
            merge_count_maps(&mut left, right);
            left
        });
    let key_increments = chunk_counts.values().copied().sum();
    sketch.add_key_counts(&chunk_counts);
    sketch.add_key_increments(key_increments);
}

fn increment_sketch_from_pair_chunk_sorted_replay(
    config: &Config,
    sketch: &mut PackedCountMinSketch,
    pairs: &[(SequenceRecord, Option<SequenceRecord>)],
    prefilter: Option<PrefilterGate<'_>>,
) {
    let mut entries = pairs
        .par_iter()
        .fold(
            || count_chunk_local_map(config, pairs),
            |mut local_counts, (r1, r2)| {
                increment_pair_counts_with_prefilter(
                    config,
                    &mut local_counts,
                    r1,
                    r2.as_ref(),
                    prefilter,
                );
                local_counts
            },
        )
        .map(|counts| counts.into_iter().collect::<Vec<_>>())
        .reduce(Vec::new, |mut left, mut right| {
            left.append(&mut right);
            left
        });
    entries.par_sort_unstable_by(|(left, _), (right, _)| left.cmp(right));

    let mut key_increments = 0u64;
    let mut iter = entries.into_iter();
    let Some((mut current_key, mut current_count)) = iter.next() else {
        return;
    };
    for (key, count) in iter {
        if key == current_key {
            current_count = current_count.saturating_add(count);
        } else {
            key_increments = key_increments.saturating_add(current_count);
            sketch.add_key_count(&current_key, current_count);
            current_key = key;
            current_count = count;
        }
    }
    key_increments = key_increments.saturating_add(current_count);
    sketch.add_key_count(&current_key, current_count);
    sketch.add_key_increments(key_increments);
}

fn increment_atomic_packed_sketch_from_pair_chunk(
    config: &Config,
    sketch: &AtomicPackedCountMinSketch,
    pairs: &[(SequenceRecord, Option<SequenceRecord>)],
) {
    let (key_increments, newly_occupied) = pairs
        .par_iter()
        .map(|(r1, r2)| increment_pair_atomic_packed_sketch(config, sketch, r1, r2.as_ref()))
        .reduce(
            || (0u64, 0usize),
            |left, right| {
                (
                    left.0.saturating_add(right.0),
                    left.1.saturating_add(right.1),
                )
            },
        );
    sketch.add_key_increments(key_increments);
    sketch.add_occupied_slots(newly_occupied);
}

fn increment_pair_atomic_packed_sketch(
    config: &Config,
    sketch: &AtomicPackedCountMinSketch,
    r1: &SequenceRecord,
    r2: Option<&SequenceRecord>,
) -> (u64, usize) {
    if config.remove_duplicate_kmers && config.k <= 31 {
        let keys = unique_pair_kmers(config, r1, r2);
        let mut newly_occupied = 0usize;
        for key in &keys {
            newly_occupied += sketch.add_key_count_counting_newly_occupied(key, 1);
        }
        return (keys.len() as u64, newly_occupied);
    }
    let mut key_increments = 0u64;
    let mut newly_occupied = 0usize;
    for_each_kmer_for_record(r1, config, |kmer| {
        newly_occupied += sketch.add_key_count_counting_newly_occupied(&kmer, 1);
        key_increments += 1;
    });
    if let Some(mate) = r2 {
        for_each_kmer_for_record(mate, config, |kmer| {
            newly_occupied += sketch.add_key_count_counting_newly_occupied(&kmer, 1);
            key_increments += 1;
        });
    }
    (key_increments, newly_occupied)
}

fn increment_atomic_sketch_from_pair_chunk(
    config: &Config,
    sketch: &AtomicCountMinSketch,
    pairs: &[(SequenceRecord, Option<SequenceRecord>)],
    prefilter: Option<PrefilterGate<'_>>,
) {
    if !config.deterministic {
        let (key_increments, newly_occupied) = pairs
            .par_iter()
            .map(|(r1, r2)| {
                increment_pair_atomic_sketch_direct(config, sketch, r1, r2.as_ref(), prefilter)
            })
            .reduce(
                || (0u64, 0usize),
                |left, right| {
                    (
                        left.0.saturating_add(right.0),
                        left.1.saturating_add(right.1),
                    )
                },
            );
        sketch.add_key_increments(key_increments);
        sketch.add_occupied_slots(newly_occupied);
        return;
    }

    let mut entries = pairs
        .par_iter()
        .fold(
            || count_chunk_local_map(config, pairs),
            |mut local_counts, (r1, r2)| {
                increment_pair_counts_with_prefilter(
                    config,
                    &mut local_counts,
                    r1,
                    r2.as_ref(),
                    prefilter,
                );
                local_counts
            },
        )
        .map(|counts| counts.into_iter().collect::<Vec<_>>())
        .reduce(Vec::new, |mut left, mut right| {
            left.append(&mut right);
            left
        });
    entries.par_sort_unstable_by(|(left, _), (right, _)| left.cmp(right));

    let mut key_increments = 0u64;
    let mut iter = entries.into_iter();
    let Some((mut current_key, mut current_count)) = iter.next() else {
        return;
    };
    for (key, count) in iter {
        if key == current_key {
            current_count = current_count.saturating_add(count);
        } else {
            key_increments = key_increments.saturating_add(current_count);
            sketch.add_key_count(&current_key, current_count);
            current_key = key;
            current_count = count;
        }
    }
    key_increments = key_increments.saturating_add(current_count);
    sketch.add_key_count(&current_key, current_count);
    sketch.add_key_increments(key_increments);
}

fn increment_pair_atomic_sketch_direct(
    config: &Config,
    sketch: &AtomicCountMinSketch,
    r1: &SequenceRecord,
    r2: Option<&SequenceRecord>,
    prefilter: Option<PrefilterGate<'_>>,
) -> (u64, usize) {
    if config.remove_duplicate_kmers && config.k <= 31 {
        let keys = unique_pair_kmers(config, r1, r2);
        let mut key_increments = 0u64;
        let mut newly_occupied = 0usize;
        for key in &keys {
            if prefilter.is_none_or(|gate| gate.should_count_in_main(key)) {
                newly_occupied += sketch.add_key_count_counting_newly_occupied(key, 1);
                key_increments += 1;
            }
        }
        return (key_increments, newly_occupied);
    }

    let mut key_increments = 0u64;
    let mut newly_occupied = 0usize;
    for_each_kmer_for_record(r1, config, |kmer| {
        if prefilter.is_none_or(|gate| gate.should_count_in_main(&kmer)) {
            newly_occupied += sketch.add_key_count_counting_newly_occupied(&kmer, 1);
            key_increments += 1;
        }
    });
    if let Some(mate) = r2 {
        for_each_kmer_for_record(mate, config, |kmer| {
            if prefilter.is_none_or(|gate| gate.should_count_in_main(&kmer)) {
                newly_occupied += sketch.add_key_count_counting_newly_occupied(&kmer, 1);
                key_increments += 1;
            }
        });
    }
    (key_increments, newly_occupied)
}

#[cfg(test)]
fn retain_prefilter_saturated_counts(counts: &mut CountMap, prefilter: Option<PrefilterGate<'_>>) {
    let Some(prefilter) = prefilter else {
        return;
    };
    counts.retain(|key, _| prefilter.should_count_in_main(key));
}

fn merge_count_maps(counts: &mut CountMap, source: CountMap) {
    for (kmer, count) in source {
        *counts.entry(kmer).or_insert(0) += count;
    }
}

fn trim_pair(config: &Config, r1: &mut SequenceRecord, r2: Option<&mut SequenceRecord>) {
    if !config.trim_left && !config.trim_right {
        return;
    }
    trim_record(config, r1);
    if let Some(mate) = r2 {
        trim_record(config, mate);
    }
}

fn trim_record(config: &Config, record: &mut SequenceRecord) {
    if record.is_empty() {
        return;
    }
    let (left0, right0) = if config.trim_optimal {
        optimal_trim_amounts(record, config)
    } else if config.trim_window {
        (0, window_trim_right_amount(record, config))
    } else {
        simple_trim_amounts(record, config)
    };
    let left = if config.trim_left { left0 } else { 0 };
    let right = if config.trim_right { right0 } else { 0 };
    trim_by_amount(record, left, right, 1);
}

fn optimal_trim_amounts(record: &SequenceRecord, config: &Config) -> (usize, usize) {
    let avg_error_rate = config
        .trim_optimal_bias
        .unwrap_or_else(|| phred_to_prob_error(config.trim_quality));
    if let Some(qualities) = record.qualities.as_deref() {
        let nprob = (avg_error_rate * 1.1).clamp(0.75, 1.0);
        let mut max_score = 0.0f64;
        let mut score = 0.0f64;
        let mut max_loc = 0usize;
        let mut max_count = 0usize;
        let mut count = 0usize;

        for (idx, (&base, &quality)) in record.bases.iter().zip(qualities).enumerate() {
            let phred = quality.saturating_sub(33);
            let prob_error = if base == b'N' || phred < 1 {
                nprob
            } else {
                phred_to_prob_error(f64::from(phred))
            };
            score += avg_error_rate - prob_error;
            if score > 0.0 {
                count += 1;
                if score > max_score || (score == max_score && count > max_count) {
                    max_score = score;
                    max_count = count;
                    max_loc = idx;
                }
            } else {
                score = 0.0;
                count = 0;
            }
        }

        if max_score > 0.0 {
            (max_loc + 1 - max_count, record.len() - max_loc - 1)
        } else {
            (0, record.len())
        }
    } else if avg_error_rate >= 1.0 {
        (0, 0)
    } else {
        (
            test_left_n(&record.bases, config.trim_min_good_interval),
            test_right_n(&record.bases, config.trim_min_good_interval),
        )
    }
}

fn simple_trim_amounts(record: &SequenceRecord, config: &Config) -> (usize, usize) {
    let trimq = config.trim_quality as u8;
    if let Some(qualities) = record.qualities.as_deref() {
        (
            test_left_quality(qualities, trimq, config.trim_min_good_interval),
            test_right_quality(qualities, trimq, config.trim_min_good_interval),
        )
    } else {
        (
            test_left_n(&record.bases, config.trim_min_good_interval),
            test_right_n(&record.bases, config.trim_min_good_interval),
        )
    }
}

fn window_trim_right_amount(record: &SequenceRecord, config: &Config) -> usize {
    let trimq = config.trim_quality as i32;
    let Some(qualities) = record.qualities.as_deref() else {
        return if trimq > 0 {
            0
        } else {
            test_right_n(&record.bases, config.trim_min_good_interval)
        };
    };
    if qualities.len() < config.trim_window_length {
        return if trimq > 0 {
            0
        } else {
            test_right_n(&record.bases, config.trim_min_good_interval)
        };
    }

    let Ok(window) = isize::try_from(config.trim_window_length) else {
        return 0;
    };
    let threshold = (config.trim_window_length as i32 * trimq).max(1);
    let mut sum = 0i32;
    for (idx, &quality) in qualities.iter().enumerate() {
        let Ok(idx) = isize::try_from(idx) else {
            return 0;
        };
        let j = idx - window;
        sum += i32::from(quality.saturating_sub(33));
        if j >= -1 {
            if j >= 0 {
                sum -= i32::from(qualities[j as usize].saturating_sub(33));
            }
            if sum < threshold {
                return qualities.len() - j as usize - 1;
            }
        }
    }
    0
}

fn test_left_quality(qualities: &[u8], trimq: u8, min_good_interval: usize) -> usize {
    let mut good = 0usize;
    let mut last_bad = None;
    for (idx, &quality) in qualities.iter().enumerate() {
        if good >= min_good_interval {
            break;
        }
        if quality.saturating_sub(33) > trimq {
            good += 1;
        } else {
            good = 0;
            last_bad = Some(idx);
        }
    }
    last_bad.map_or(0, |idx| idx + 1)
}

fn test_right_quality(qualities: &[u8], trimq: u8, min_good_interval: usize) -> usize {
    let mut good = 0usize;
    let mut last_bad = qualities.len();
    for (idx, &quality) in qualities.iter().enumerate().rev() {
        if good >= min_good_interval {
            break;
        }
        if quality.saturating_sub(33) > trimq {
            good += 1;
        } else {
            good = 0;
            last_bad = idx;
        }
    }
    qualities.len() - last_bad
}

fn test_left_n(bases: &[u8], min_good_interval: usize) -> usize {
    let mut good = 0usize;
    let mut last_bad = None;
    for (idx, &base) in bases.iter().enumerate() {
        if good >= min_good_interval {
            break;
        }
        if base != b'N' {
            good += 1;
        } else {
            good = 0;
            last_bad = Some(idx);
        }
    }
    last_bad.map_or(0, |idx| idx + 1)
}

fn test_right_n(bases: &[u8], min_good_interval: usize) -> usize {
    let mut good = 0usize;
    let mut last_bad = bases.len();
    for (idx, &base) in bases.iter().enumerate().rev() {
        if good >= min_good_interval {
            break;
        }
        if base != b'N' {
            good += 1;
        } else {
            good = 0;
            last_bad = idx;
        }
    }
    bases.len() - last_bad
}

fn trim_by_amount(
    record: &mut SequenceRecord,
    mut left_trim: usize,
    mut right_trim: usize,
    min_resulting_length: usize,
) -> usize {
    let len = record.len();
    if len == 0 {
        return 0;
    }
    let min_resulting_length = min_resulting_length.min(len);
    if left_trim + right_trim + min_resulting_length > len {
        right_trim = 1usize.max(len.saturating_sub(min_resulting_length));
        left_trim = 0;
    }
    let total = left_trim + right_trim;
    if total > 0 {
        record.bases = record.bases[left_trim..len - right_trim].to_vec();
        if let Some(qualities) = record.qualities.take() {
            let qlen = qualities.len();
            record.qualities = if total >= qlen {
                None
            } else {
                Some(qualities[left_trim..qlen - right_trim].to_vec())
            };
        }
    }
    total
}

fn phred_to_prob_error(q: f64) -> f64 {
    if q <= 0.0 {
        0.75
    } else if q <= 1.0 {
        0.75 - q * 0.05
    } else {
        0.7_f64.min(10_f64.powf(-0.1 * q))
    }
}

fn increment_sparse_hist_from_analysis(
    hist: &mut SparseHist,
    analysis: &ReadAnalysis,
    hist_len: usize,
) {
    for depth in &analysis.coverage_desc {
        if *depth < 0 {
            continue;
        }
        let idx = (*depth as usize).min(hist_len - 1);
        *hist.entry(idx).or_insert(0) += 1;
    }
}

#[cfg(test)]
fn increment_hist_from_pair_chunk(
    config: &Config,
    hist_counts: &dyn CountLookup,
    keep_filter_counts: Option<&dyn CountLookup>,
    hist: &mut [u64],
    pairs: &[AnalysisPair],
) {
    let chunk_hist = sparse_hist_from_pair_chunk(config, hist_counts, keep_filter_counts, pairs);
    merge_sparse_hist_into_dense(hist, chunk_hist);
}

fn sparse_hist_from_pair_chunk(
    config: &Config,
    hist_counts: &dyn CountLookup,
    keep_filter_counts: Option<&dyn CountLookup>,
    pairs: &[AnalysisPair],
) -> SparseHist {
    pairs
        .par_iter()
        .fold(SparseHist::default, |mut local_hist, (r1, r2, rand)| {
            if let Some(input_counts) = keep_filter_counts {
                let decision = decide_pair(config, input_counts, r1, r2.as_ref(), *rand);
                if decision.toss {
                    return local_hist;
                }
            }

            let analysis = analyze_pair(config, hist_counts, r1, r2.as_ref());
            increment_sparse_hist_from_analysis(&mut local_hist, &analysis.read1, config.hist_len);
            if let Some(read2) = &analysis.read2 {
                increment_sparse_hist_from_analysis(&mut local_hist, read2, config.hist_len);
            }
            local_hist
        })
        .reduce(SparseHist::default, |mut left, right| {
            merge_sparse_hist(&mut left, right);
            left
        })
}

fn merge_sparse_hist(target: &mut SparseHist, source: SparseHist) {
    for (idx, count) in source {
        *target.entry(idx).or_insert(0) += count;
    }
}

#[cfg(test)]
fn merge_sparse_hist_into_dense(target: &mut [u64], source: SparseHist) {
    for (idx, count) in source {
        target[idx] += count;
    }
}

fn increment_sparse_read_hist(
    hist: &mut SparseReadDepthHist,
    analysis: &ReadAnalysis,
    read_len: usize,
    hist_len: usize,
) {
    if !analysis.had_kmer_windows {
        return;
    }
    let depth = analysis.depth_al.or(analysis.true_depth).unwrap_or(0);
    let idx = (depth as usize).min(hist_len - 1);
    let entry = hist.entry(idx).or_insert((0, 0));
    entry.0 += 1;
    entry.1 += read_len as u64;
}

#[cfg(test)]
fn increment_read_hist_from_pair_chunk(
    config: &Config,
    hist_counts: &dyn CountLookup,
    keep_filter_counts: Option<&dyn CountLookup>,
    hist: &mut ReadDepthHistogram,
    pairs: &[AnalysisPair],
) {
    let chunk_hist =
        sparse_read_hist_from_pair_chunk(config, hist_counts, keep_filter_counts, pairs);
    merge_sparse_read_depth_hist_into_dense(hist, chunk_hist);
}

fn sparse_read_hist_from_pair_chunk(
    config: &Config,
    hist_counts: &dyn CountLookup,
    keep_filter_counts: Option<&dyn CountLookup>,
    pairs: &[AnalysisPair],
) -> SparseReadDepthHist {
    pairs
        .par_iter()
        .fold(
            SparseReadDepthHist::default,
            |mut local_hist, (r1, r2, rand)| {
                if let Some(input_counts) = keep_filter_counts {
                    let decision = decide_pair(config, input_counts, r1, r2.as_ref(), *rand);
                    if decision.toss {
                        return local_hist;
                    }
                }

                let analysis = analyze_pair(config, hist_counts, r1, r2.as_ref());
                increment_sparse_read_hist(
                    &mut local_hist,
                    &analysis.read1,
                    r1.len(),
                    config.hist_len,
                );
                if let (Some(read2_analysis), Some(read2)) = (&analysis.read2, r2.as_ref()) {
                    increment_sparse_read_hist(
                        &mut local_hist,
                        read2_analysis,
                        read2.len(),
                        config.hist_len,
                    );
                }
                local_hist
            },
        )
        .reduce(SparseReadDepthHist::default, |mut left, right| {
            merge_sparse_read_depth_hist(&mut left, right);
            left
        })
}

#[cfg(test)]
fn increment_hist_and_read_hist_from_pair_chunk(
    config: &Config,
    hist_counts: &dyn CountLookup,
    keep_filter_counts: Option<&dyn CountLookup>,
    depth_hist: &mut [u64],
    read_hist: &mut ReadDepthHistogram,
    pairs: &[AnalysisPair],
) {
    let (chunk_depth_hist, chunk_read_hist) =
        sparse_hist_and_read_hist_from_pair_chunk(config, hist_counts, keep_filter_counts, pairs);
    merge_sparse_hist_into_dense(depth_hist, chunk_depth_hist);
    merge_sparse_read_depth_hist_into_dense(read_hist, chunk_read_hist);
}

fn sparse_hist_and_read_hist_from_pair_chunk(
    config: &Config,
    hist_counts: &dyn CountLookup,
    keep_filter_counts: Option<&dyn CountLookup>,
    pairs: &[AnalysisPair],
) -> (SparseHist, SparseReadDepthHist) {
    pairs
        .par_iter()
        .fold(
            || (SparseHist::default(), SparseReadDepthHist::default()),
            |mut local, (r1, r2, rand)| {
                if let Some(input_counts) = keep_filter_counts {
                    let decision = decide_pair(config, input_counts, r1, r2.as_ref(), *rand);
                    if decision.toss {
                        return local;
                    }
                }

                let analysis = analyze_pair(config, hist_counts, r1, r2.as_ref());
                increment_sparse_hist_from_analysis(&mut local.0, &analysis.read1, config.hist_len);
                increment_sparse_read_hist(
                    &mut local.1,
                    &analysis.read1,
                    r1.len(),
                    config.hist_len,
                );
                if let Some(read2_analysis) = &analysis.read2 {
                    increment_sparse_hist_from_analysis(
                        &mut local.0,
                        read2_analysis,
                        config.hist_len,
                    );
                    if let Some(read2) = r2.as_ref() {
                        increment_sparse_read_hist(
                            &mut local.1,
                            read2_analysis,
                            read2.len(),
                            config.hist_len,
                        );
                    }
                }
                local
            },
        )
        .reduce(
            || (SparseHist::default(), SparseReadDepthHist::default()),
            |mut left, right| {
                merge_sparse_hist(&mut left.0, right.0);
                merge_sparse_read_depth_hist(&mut left.1, right.1);
                left
            },
        )
}

fn merge_sparse_read_depth_hist(target: &mut SparseReadDepthHist, source: SparseReadDepthHist) {
    for (idx, (reads, bases)) in source {
        let entry = target.entry(idx).or_insert((0, 0));
        entry.0 += reads;
        entry.1 += bases;
    }
}

#[cfg(test)]
fn merge_sparse_read_depth_hist_into_dense(
    target: &mut ReadDepthHistogram,
    source: SparseReadDepthHist,
) {
    for (idx, (reads, bases)) in source {
        target.reads[idx] += reads;
        target.bases[idx] += bases;
    }
}

#[cfg(test)]
fn write_depth_hist(path: &Path, raw_hist: &[u64], config: &Config) -> Result<()> {
    let mut writer = crate::seqio::create_output(path, config.overwrite || config.append)
        .with_context(|| format!("creating histogram {}", path.display()))?;
    match config.hist_columns {
        1 => writeln!(writer, "#tUnique_Kmers")?,
        2 => writeln!(writer, "#Depth\tUnique_Kmers")?,
        3 => writeln!(writer, "#Depth\tRaw_Count\tUnique_Kmers")?,
        _ => unreachable!("validated hist column count"),
    }

    let total_raw = raw_hist.iter().copied().fold(0u64, u64::saturating_add);
    let mut seen_raw = 0u64;
    let lim = raw_hist.len().saturating_sub(1);
    for depth in 0..lim {
        let raw = adjusted_depth_hist_raw(raw_hist, config.zero_bin, depth);
        seen_raw = seen_raw.saturating_add(raw);
        let unique = unique_from_raw(depth, raw);
        if config.print_zero_coverage || unique > 0 || config.hist_columns == 1 {
            write_hist_row(&mut writer, config.hist_columns, depth, raw, unique)?;
        }
        if seen_raw >= total_raw {
            break;
        }
    }

    let overflow_raw = (lim..raw_hist.len())
        .map(|depth| adjusted_depth_hist_raw(raw_hist, config.zero_bin, depth))
        .fold(0u64, u64::saturating_add);
    if overflow_raw > 0 {
        write_hist_row(
            &mut writer,
            config.hist_columns,
            lim,
            overflow_raw,
            unique_from_raw(lim, overflow_raw),
        )?;
    }
    writer.flush()?;
    Ok(())
}

fn write_sparse_depth_hist(
    path: &Path,
    raw_hist: &SparseHist,
    hist_len: usize,
    config: &Config,
) -> Result<()> {
    let mut writer = crate::seqio::create_output(path, config.overwrite || config.append)
        .with_context(|| format!("creating histogram {}", path.display()))?;
    match config.hist_columns {
        1 => writeln!(writer, "#tUnique_Kmers")?,
        2 => writeln!(writer, "#Depth\tUnique_Kmers")?,
        3 => writeln!(writer, "#Depth\tRaw_Count\tUnique_Kmers")?,
        _ => unreachable!("validated hist column count"),
    }

    let hist_len = hist_len.max(1);
    let lim = hist_len.saturating_sub(1);
    let total_raw = raw_hist.values().copied().fold(0u64, u64::saturating_add);
    let mut seen_raw = 0u64;

    if config.print_zero_coverage || config.hist_columns == 1 {
        for depth in 0..lim {
            let raw = adjusted_sparse_depth_hist_raw(raw_hist, hist_len, config.zero_bin, depth);
            seen_raw = seen_raw.saturating_add(raw);
            write_hist_row(
                &mut writer,
                config.hist_columns,
                depth,
                raw,
                unique_from_raw(depth, raw),
            )?;
            if seen_raw >= total_raw {
                break;
            }
        }
    } else {
        let mut depths: Vec<usize> = raw_hist
            .iter()
            .filter_map(|(&depth, &raw)| {
                let mapped_depth = if !config.zero_bin && hist_len > 1 && depth == 0 {
                    1
                } else {
                    depth
                };
                (mapped_depth < lim && raw > 0).then_some(mapped_depth)
            })
            .collect();
        depths.sort_unstable();
        depths.dedup();
        for depth in depths {
            let raw = adjusted_sparse_depth_hist_raw(raw_hist, hist_len, config.zero_bin, depth);
            seen_raw = seen_raw.saturating_add(raw);
            let unique = unique_from_raw(depth, raw);
            if unique > 0 {
                write_hist_row(&mut writer, config.hist_columns, depth, raw, unique)?;
            }
            if seen_raw >= total_raw {
                break;
            }
        }
    }

    let mut overflow_depths: Vec<usize> = raw_hist
        .keys()
        .copied()
        .filter_map(|depth| {
            let mapped_depth = if !config.zero_bin && hist_len > 1 && depth == 0 {
                1
            } else {
                depth
            };
            (mapped_depth >= lim).then_some(mapped_depth)
        })
        .collect();
    overflow_depths.sort_unstable();
    overflow_depths.dedup();
    let overflow_raw = overflow_depths.into_iter().fold(0u64, |sum, depth| {
        sum.saturating_add(adjusted_sparse_depth_hist_raw(
            raw_hist,
            hist_len,
            config.zero_bin,
            depth,
        ))
    });
    if overflow_raw > 0 {
        write_hist_row(
            &mut writer,
            config.hist_columns,
            lim,
            overflow_raw,
            unique_from_raw(lim, overflow_raw),
        )?;
    }
    writer.flush()?;
    Ok(())
}

#[cfg(test)]
fn adjusted_depth_hist_raw(raw_hist: &[u64], zero_bin: bool, depth: usize) -> u64 {
    let raw = raw_hist.get(depth).copied().unwrap_or(0);
    if zero_bin || raw_hist.len() <= 1 {
        return raw;
    }
    match depth {
        0 => 0,
        1 => raw.saturating_add(raw_hist[0]),
        _ => raw,
    }
}

fn adjusted_sparse_depth_hist_raw(
    raw_hist: &SparseHist,
    hist_len: usize,
    zero_bin: bool,
    depth: usize,
) -> u64 {
    let raw = raw_hist.get(&depth).copied().unwrap_or(0);
    if zero_bin || hist_len <= 1 {
        return raw;
    }
    match depth {
        0 => 0,
        1 => raw.saturating_add(raw_hist.get(&0).copied().unwrap_or(0)),
        _ => raw,
    }
}

#[cfg(test)]
fn sparse_hist_to_dense(raw_hist: &SparseHist, hist_len: usize) -> Vec<u64> {
    let mut dense = vec![0u64; hist_len.max(1)];
    for (&depth, &raw) in raw_hist {
        let idx = depth.min(dense.len() - 1);
        dense[idx] = dense[idx].saturating_add(raw);
    }
    dense
}

fn sparse_hist_to_peak_dense(raw_hist: &SparseHist, hist_len: usize) -> Vec<u64> {
    let hist_len = hist_len.max(1);
    let last_index = hist_len - 1;
    let last_nonzero = raw_hist
        .iter()
        .filter_map(|(&depth, &raw)| (raw > 0).then_some(depth.min(last_index)))
        .max()
        .unwrap_or(0);
    let dense_len = hist_len.min(
        last_nonzero
            .saturating_add(PEAK_COMPACT_ZERO_TAIL)
            .saturating_add(1),
    );
    let mut dense = vec![0u64; dense_len.max(1)];
    for (&depth, &raw) in raw_hist {
        if raw == 0 {
            continue;
        }
        let idx = depth.min(last_index);
        if idx < dense.len() {
            dense[idx] = dense[idx].saturating_add(raw);
        } else {
            dense.resize(idx + 1, 0);
            dense[idx] = dense[idx].saturating_add(raw);
        }
    }
    dense
}

fn write_hist_row(
    writer: &mut Box<dyn Write>,
    columns: u8,
    depth: usize,
    raw: u64,
    unique: u64,
) -> Result<()> {
    match columns {
        1 => writeln!(writer, "{unique}")?,
        2 => writeln!(writer, "{depth}\t{unique}")?,
        3 => writeln!(writer, "{depth}\t{raw}\t{unique}")?,
        _ => unreachable!("validated hist column count"),
    }
    Ok(())
}

#[cfg(test)]
fn write_read_depth_hist(path: &Path, hist: &ReadDepthHistogram, config: &Config) -> Result<()> {
    let mut writer = crate::seqio::create_output(path, config.overwrite || config.append)
        .with_context(|| format!("creating read histogram {}", path.display()))?;
    writeln!(writer, "#Depth\tReads\tBases")?;

    let total_reads: u64 = hist.reads.iter().sum();
    let mut seen_reads = 0u64;
    let lim = hist.reads.len().saturating_sub(1);

    for depth in 0..lim {
        let reads = hist.reads[depth];
        let bases = hist.bases[depth];
        seen_reads += reads;
        if config.print_zero_coverage || bases > 0 {
            writeln!(writer, "{depth}\t{reads}\t{bases}")?;
        }
        if seen_reads >= total_reads {
            break;
        }
    }

    let overflow_reads: u64 = hist.reads.iter().skip(lim).sum();
    let overflow_bases: u64 = hist.bases.iter().skip(lim).sum();
    if overflow_reads > 0 || overflow_bases > 0 {
        writeln!(writer, "{lim}\t{overflow_reads}\t{overflow_bases}")?;
    }
    writer.flush()?;
    Ok(())
}

fn write_sparse_read_depth_hist(
    path: &Path,
    hist: &SparseReadDepthHist,
    hist_len: usize,
    config: &Config,
) -> Result<()> {
    let mut writer = crate::seqio::create_output(path, config.overwrite || config.append)
        .with_context(|| format!("creating read histogram {}", path.display()))?;
    writeln!(writer, "#Depth\tReads\tBases")?;

    let hist_len = hist_len.max(1);
    let lim = hist_len.saturating_sub(1);
    let total_reads = hist
        .values()
        .map(|(reads, _)| *reads)
        .fold(0u64, u64::saturating_add);
    let mut seen_reads = 0u64;

    if config.print_zero_coverage {
        for depth in 0..lim {
            let (reads, bases) = hist.get(&depth).copied().unwrap_or_default();
            seen_reads = seen_reads.saturating_add(reads);
            writeln!(writer, "{depth}\t{reads}\t{bases}")?;
            if seen_reads >= total_reads {
                break;
            }
        }
    } else {
        let mut depths: Vec<usize> = hist.keys().copied().filter(|depth| *depth < lim).collect();
        depths.sort_unstable();
        for depth in depths {
            let (reads, bases) = hist.get(&depth).copied().unwrap_or_default();
            seen_reads = seen_reads.saturating_add(reads);
            if bases > 0 {
                writeln!(writer, "{depth}\t{reads}\t{bases}")?;
            }
            if seen_reads >= total_reads {
                break;
            }
        }
    }

    let (overflow_reads, overflow_bases) = hist.iter().filter(|(depth, _)| **depth >= lim).fold(
        (0u64, 0u64),
        |(read_sum, base_sum), (_, (reads, bases))| {
            (
                read_sum.saturating_add(*reads),
                base_sum.saturating_add(*bases),
            )
        },
    );
    if overflow_reads > 0 || overflow_bases > 0 {
        writeln!(writer, "{lim}\t{overflow_reads}\t{overflow_bases}")?;
    }
    writer.flush()?;
    Ok(())
}

fn write_quality_hist(path: &Path, hist: &[u64], config: &Config) -> Result<()> {
    let mut writer = crate::seqio::create_output(path, config.overwrite || config.append)
        .with_context(|| format!("creating quality histogram {}", path.display()))?;
    writeln!(writer, "#Quality\tBases")?;

    let total_bases: u64 = hist.iter().sum();
    let mut seen_bases = 0u64;
    let lim = hist.len().saturating_sub(1);

    for (quality, bases) in hist.iter().copied().enumerate().take(lim) {
        seen_bases += bases;
        if config.print_zero_coverage || bases > 0 {
            writeln!(writer, "{quality}\t{bases}")?;
        }
        if seen_bases >= total_bases {
            break;
        }
    }

    let overflow_bases: u64 = hist.iter().skip(lim).sum();
    if overflow_bases > 0 {
        writeln!(writer, "{lim}\t{overflow_bases}")?;
    }
    writer.flush()?;
    Ok(())
}

fn write_quality_count_hist(
    path: &Path,
    first: &[u64],
    second: &[u64],
    paired: bool,
    config: &Config,
) -> Result<()> {
    let mut writer = crate::seqio::create_output(path, config.overwrite || config.append)
        .with_context(|| format!("creating quality-count histogram {}", path.display()))?;
    writeln!(
        writer,
        "#Quality\tcount1\tfraction1{}",
        if paired { "\tcount2\tfraction2" } else { "" }
    )?;
    write_paired_quality_count_rows(&mut writer, first, second, paired, config)?;
    writer.flush()?;
    Ok(())
}

fn write_average_quality_hist(
    path: &Path,
    first: &[u64],
    second: &[u64],
    paired: bool,
    config: &Config,
) -> Result<()> {
    let mut writer = crate::seqio::create_output(path, config.overwrite || config.append)
        .with_context(|| format!("creating average-quality histogram {}", path.display()))?;
    writeln!(
        writer,
        "#Quality\tcount1\tfraction1{}",
        if paired { "\tcount2\tfraction2" } else { "" }
    )?;
    write_paired_quality_count_rows(&mut writer, first, second, paired, config)?;
    writer.flush()?;
    Ok(())
}

fn write_paired_quality_count_rows(
    writer: &mut Box<dyn Write>,
    first: &[u64],
    second: &[u64],
    paired: bool,
    config: &Config,
) -> Result<()> {
    let total1: u64 = first.iter().sum();
    let total2: u64 = second.iter().sum();
    let mut remaining = total1 + if paired { total2 } else { 0 };
    let denom1 = total1.max(1) as f64;
    let denom2 = total2.max(1) as f64;

    for (quality, count1) in first.iter().copied().enumerate() {
        let count2 = second.get(quality).copied().unwrap_or(0);
        if count1 > 0 || (paired && count2 > 0) || config.print_zero_coverage {
            write!(writer, "{quality}\t{count1}\t{:.5}", count1 as f64 / denom1)?;
            if paired {
                write!(writer, "\t{count2}\t{:.5}", count2 as f64 / denom2)?;
            }
            writeln!(writer)?;
        }
        remaining = remaining.saturating_sub(count1 + if paired { count2 } else { 0 });
        if remaining == 0 && !config.print_zero_coverage {
            break;
        }
    }
    Ok(())
}

fn write_overall_base_quality_hist(path: &Path, hist: &[u64], config: &Config) -> Result<()> {
    let mut writer = crate::seqio::create_output(path, config.overwrite || config.append)
        .with_context(|| format!("creating overall base-quality histogram {}", path.display()))?;
    let median = percentile_histogram(hist, 0.5);
    let mean = average_histogram(hist);
    let stdev = stdev_histogram(hist, mean, 0);
    let mean30 = average_histogram_min(hist, 30);
    let stdev30 = stdev_histogram(hist, mean30, 30);
    writeln!(writer, "#Median\t{median}")?;
    writeln!(writer, "#Mean\t{mean:.3}")?;
    writeln!(writer, "#STDev\t{stdev:.3}")?;
    writeln!(writer, "#Mean_30\t{mean30:.3}")?;
    writeln!(writer, "#STDev_30\t{stdev30:.3}")?;
    writeln!(writer, "#Quality\tbases\tfraction")?;

    let total: u64 = hist.iter().sum();
    let denom = total.max(1) as f64;
    let mut remaining = total;
    for (quality, bases) in hist.iter().copied().enumerate() {
        if bases > 0 || config.print_zero_coverage {
            writeln!(writer, "{quality}\t{bases}\t{:.5}", bases as f64 / denom)?;
        }
        remaining = remaining.saturating_sub(bases);
        if remaining == 0 && !config.print_zero_coverage {
            break;
        }
    }
    writer.flush()?;
    Ok(())
}

fn write_base_quality_hist(
    path: &Path,
    hist: &QualitySideHistograms,
    config: &Config,
) -> Result<()> {
    let mut writer = crate::seqio::create_output(path, config.overwrite || config.append)
        .with_context(|| format!("creating base-quality histogram {}", path.display()))?;
    write!(
        writer,
        "#BaseNum\tcount_1\tmin_1\tmax_1\tmean_1\tQ1_1\tmed_1\tQ3_1\tLW_1\tRW_1"
    )?;
    if hist.paired {
        write!(
            writer,
            "\tcount_2\tmin_2\tmax_2\tmean_2\tQ1_2\tmed_2\tQ3_2\tLW_2\tRW_2"
        )?;
    }
    writeln!(writer)?;

    for pos in 0..hist.first_by_pos.len() {
        let sum1: u64 = hist.first_by_pos[pos].iter().sum();
        let sum2: u64 = hist.second_by_pos[pos].iter().sum();
        if sum1 == 0 && sum2 == 0 && !config.print_zero_coverage {
            break;
        }
        write!(writer, "{pos}")?;
        write_base_quality_summary(&mut writer, &hist.first_by_pos[pos])?;
        if hist.paired {
            write_base_quality_summary(&mut writer, &hist.second_by_pos[pos])?;
        }
        writeln!(writer)?;
    }
    writer.flush()?;
    Ok(())
}

fn write_base_quality_summary(writer: &mut Box<dyn Write>, hist: &[u64]) -> Result<()> {
    let count: u64 = hist.iter().sum();
    let min = min_histogram(hist);
    let max = max_histogram(hist);
    let mean = average_histogram(hist);
    let q1 = percentile_histogram(hist, 0.25);
    let med = percentile_histogram(hist, 0.5);
    let q3 = percentile_histogram(hist, 0.75);
    let left_whisker = percentile_histogram(hist, 0.02);
    let right_whisker = percentile_histogram(hist, 0.98);
    write!(
        writer,
        "\t{count}\t{min}\t{max}\t{mean:.2}\t{q1}\t{med}\t{q3}\t{left_whisker}\t{right_whisker}"
    )?;
    Ok(())
}

fn min_histogram(hist: &[u64]) -> usize {
    hist.iter().position(|count| *count > 0).unwrap_or_default()
}

fn max_histogram(hist: &[u64]) -> usize {
    hist.iter()
        .rposition(|count| *count > 0)
        .unwrap_or_default()
}

fn mode_histogram(hist: &[u64]) -> usize {
    hist.iter()
        .copied()
        .enumerate()
        .max_by_key(|(_, count)| *count)
        .map_or(0, |(idx, _)| idx)
}

fn percentile_histogram(hist: &[u64], percentile: f64) -> usize {
    let total: u64 = hist.iter().sum();
    if total == 0 {
        return 0;
    }
    let threshold = ((total as f64) * percentile).ceil().max(1.0) as u64;
    let mut seen = 0u64;
    for (idx, count) in hist.iter().copied().enumerate() {
        seen += count;
        if seen >= threshold {
            return idx;
        }
    }
    hist.len().saturating_sub(1)
}

fn average_histogram(hist: &[u64]) -> f64 {
    average_histogram_min(hist, 0)
}

fn average_histogram_min(hist: &[u64], min_quality: usize) -> f64 {
    let mut count = 0u64;
    let mut sum = 0u64;
    for (quality, bases) in hist.iter().copied().enumerate().skip(min_quality) {
        count += bases;
        sum += quality as u64 * bases;
    }
    if count == 0 {
        0.0
    } else {
        sum as f64 / count as f64
    }
}

fn stdev_histogram(hist: &[u64], mean: f64, min_quality: usize) -> f64 {
    let mut count = 0u64;
    let mut sum = 0.0;
    for (quality, bases) in hist.iter().copied().enumerate().skip(min_quality) {
        count += bases;
        let delta = quality as f64 - mean;
        sum += delta * delta * bases as f64;
    }
    if count == 0 {
        0.0
    } else {
        (sum / count as f64).sqrt()
    }
}

fn write_length_hist(path: &Path, hist: &ReadDepthHistogram, config: &Config) -> Result<()> {
    let mut writer = crate::seqio::create_output(path, config.overwrite || config.append)
        .with_context(|| format!("creating length histogram {}", path.display()))?;
    writeln!(writer, "#Length\tReads\tBases")?;

    let total_reads: u64 = hist.reads.iter().sum();
    let mut seen_reads = 0u64;
    let lim = hist.reads.len().saturating_sub(1);

    for len in 0..lim {
        let reads = hist.reads[len];
        let bases = hist.bases[len];
        seen_reads += reads;
        if config.print_zero_coverage || reads > 0 {
            writeln!(writer, "{len}\t{reads}\t{bases}")?;
        }
        if seen_reads >= total_reads {
            break;
        }
    }

    let overflow_reads: u64 = hist.reads.iter().skip(lim).sum();
    let overflow_bases: u64 = hist.bases.iter().skip(lim).sum();
    if overflow_reads > 0 || overflow_bases > 0 {
        writeln!(writer, "{lim}\t{overflow_reads}\t{overflow_bases}")?;
    }
    writer.flush()?;
    Ok(())
}

fn write_gc_hist(path: &Path, hist: &ReadDepthHistogram, config: &Config) -> Result<()> {
    let mut writer = crate::seqio::create_output(path, config.overwrite || config.append)
        .with_context(|| format!("creating GC histogram {}", path.display()))?;
    writeln!(writer, "#GC_Bin\tReads\tBases")?;

    let total_reads: u64 = hist.reads.iter().sum();
    let mut seen_reads = 0u64;
    for (bin, reads) in hist.reads.iter().copied().enumerate() {
        let bases = hist.bases[bin];
        seen_reads += reads;
        if config.print_zero_coverage || reads > 0 {
            writeln!(writer, "{bin}\t{reads}\t{bases}")?;
        }
        if seen_reads >= total_reads {
            break;
        }
    }
    writer.flush()?;
    Ok(())
}

fn write_base_content_hist(
    path: &Path,
    hist: &BaseContentHistogram,
    config: &Config,
) -> Result<()> {
    let mut writer = crate::seqio::create_output(path, config.overwrite || config.append)
        .with_context(|| format!("creating base-content histogram {}", path.display()))?;
    writeln!(writer, "#Pos\tA\tC\tG\tT\tN")?;
    let first_rows = write_base_content_rows(&mut writer, &hist.first, 0, config)?;
    write_base_content_rows(&mut writer, &hist.second, first_rows, config)?;
    writer.flush()?;
    Ok(())
}

fn write_base_content_rows(
    writer: &mut Box<dyn Write>,
    hist: &[BaseCounts],
    offset: usize,
    config: &Config,
) -> Result<usize> {
    let rows = if config.print_zero_coverage {
        hist.len()
    } else {
        hist.iter()
            .rposition(|counts| counts.total() > 0)
            .map_or(0, |idx| idx + 1)
    };

    for (pos, counts) in hist.iter().copied().enumerate().take(rows) {
        let total = counts.total() as f64;
        let fraction = |value: u64| {
            if total == 0.0 {
                0.0
            } else {
                value as f64 / total
            }
        };
        writeln!(
            writer,
            "{}\t{:.5}\t{:.5}\t{:.5}\t{:.5}\t{:.5}",
            pos + offset,
            fraction(counts.a),
            fraction(counts.c),
            fraction(counts.g),
            fraction(counts.t),
            fraction(counts.n)
        )?;
    }
    Ok(rows)
}

fn write_entropy_hist(path: &Path, hist: &[u64], config: &Config) -> Result<()> {
    let mut writer = crate::seqio::create_output(path, config.overwrite || config.append)
        .with_context(|| format!("creating entropy histogram {}", path.display()))?;
    let bins = hist.len().saturating_sub(1).max(1);
    let mult = 1.0 / bins as f64;
    let mean = average_histogram(hist) * mult;
    let median = percentile_histogram(hist, 0.5) as f64 * mult;
    let mode = mode_histogram(hist) as f64 * mult;
    let stdev = stdev_histogram(hist, average_histogram(hist), 0) * mult;

    writeln!(writer, "#Mean\t{mean:.6}")?;
    writeln!(writer, "#Median\t{median:.6}")?;
    writeln!(writer, "#Mode\t{mode:.6}")?;
    writeln!(writer, "#STDev\t{stdev:.6}")?;
    writeln!(writer, "#Value\tCount")?;

    for (idx, count) in hist.iter().copied().enumerate() {
        if config.print_zero_coverage || count > 0 {
            writeln!(writer, "{:.4}\t{count}", idx as f64 * mult)?;
        }
    }
    writer.flush()?;
    Ok(())
}

fn write_identity_hist(path: &Path, hist: &ReadDepthHistogram, config: &Config) -> Result<()> {
    let mut writer = crate::seqio::create_output(path, config.overwrite || config.append)
        .with_context(|| format!("creating identity histogram {}", path.display()))?;
    let bins = hist.reads.len().saturating_sub(1).max(1);
    let mult = 100.0 / bins as f64;
    let mean_reads = average_histogram(&hist.reads) * mult;
    let mean_bases = average_histogram(&hist.bases) * mult;
    let median_reads = percentile_histogram(&hist.reads, 0.5) as f64 * mult;
    let median_bases = percentile_histogram(&hist.bases, 0.5) as f64 * mult;
    let mode_reads = mode_histogram(&hist.reads) as f64 * mult;
    let mode_bases = mode_histogram(&hist.bases) as f64 * mult;
    let stdev_reads = stdev_histogram(&hist.reads, average_histogram(&hist.reads), 0) * mult;
    let stdev_bases = stdev_histogram(&hist.bases, average_histogram(&hist.bases), 0) * mult;

    writeln!(writer, "#Mean_reads\t{mean_reads:.3}")?;
    writeln!(writer, "#Mean_bases\t{mean_bases:.3}")?;
    writeln!(writer, "#Median_reads\t{median_reads:.0}")?;
    writeln!(writer, "#Median_bases\t{median_bases:.0}")?;
    writeln!(writer, "#Mode_reads\t{mode_reads:.0}")?;
    writeln!(writer, "#Mode_bases\t{mode_bases:.0}")?;
    writeln!(writer, "#STDev_reads\t{stdev_reads:.3}")?;
    writeln!(writer, "#STDev_bases\t{stdev_bases:.3}")?;
    writeln!(writer, "#Identity\tReads\tBases")?;

    for (idx, reads) in hist.reads.iter().copied().enumerate() {
        let bases = hist.bases[idx];
        if config.print_zero_coverage || reads > 0 || bases > 0 {
            writeln!(writer, "{:.1}\t{reads}\t{bases}", idx as f64 * mult)?;
        }
    }
    writer.flush()?;
    Ok(())
}

fn emit_alignment_fallback_side_outputs(
    config: &Config,
    hist: &AlignmentFallbackHistograms,
) -> Result<()> {
    if let Some(path) = &config.match_hist_out {
        write_match_fallback_hist(path, hist, config)?;
    }
    if let Some(path) = &config.insert_hist_out {
        write_insert_fallback_hist(path, hist, config)?;
    }
    if let Some(path) = &config.quality_accuracy_hist_out {
        write_quality_accuracy_fallback_hist(path, hist, config)?;
    }
    if let Some(path) = &config.indel_hist_out {
        write_indel_fallback_hist(path, config)?;
    }
    if let Some(path) = &config.error_hist_out {
        write_error_fallback_hist(path, hist, config)?;
    }
    Ok(())
}

fn write_match_fallback_hist(
    path: &Path,
    hist: &AlignmentFallbackHistograms,
    config: &Config,
) -> Result<()> {
    let mut writer = crate::seqio::create_output(path, config.overwrite || config.append)
        .with_context(|| format!("creating match histogram {}", path.display()))?;
    if hist.paired {
        writeln!(
            writer,
            "#BaseNum\tMatch1\tSub1\tDel1\tIns1\tN1\tOther1\tMatch2\tSub2\tDel2\tIns2\tN2\tOther2"
        )?;
    } else {
        writeln!(writer, "#BaseNum\tMatch1\tSub1\tDel1\tIns1\tN1\tOther1")?;
    }

    for pos in 0..hist.first_match.len() {
        let first = hist.first_match[pos];
        let second = hist.second_match[pos];
        if first.matches + first.n + second.matches + second.n == 0 && !config.print_zero_coverage {
            break;
        }
        write!(writer, "{}", pos + 1)?;
        write_match_fallback_columns(&mut writer, first)?;
        if hist.paired {
            write_match_fallback_columns(&mut writer, second)?;
        }
        writeln!(writer)?;
    }
    writer.flush()?;
    Ok(())
}

fn write_match_fallback_columns(writer: &mut Box<dyn Write>, counts: MatchCounts) -> Result<()> {
    let total = (counts.matches + counts.n).max(1) as f64;
    write!(
        writer,
        "\t{:.5}\t0.00000\t0.00000\t0.00000\t{:.5}\t0.00000",
        counts.matches as f64 / total,
        counts.n as f64 / total
    )?;
    Ok(())
}

fn write_insert_fallback_hist(
    path: &Path,
    hist: &AlignmentFallbackHistograms,
    config: &Config,
) -> Result<()> {
    let mut writer = crate::seqio::create_output(path, config.overwrite || config.append)
        .with_context(|| format!("creating insert-size histogram {}", path.display()))?;
    let percent = if hist.read_count == 0 {
        0.0
    } else {
        (hist.pair_count * 2) as f64 * 100.0 / hist.read_count as f64
    };
    writeln!(writer, "#Mean\t0.000")?;
    writeln!(writer, "#Median\t0")?;
    writeln!(writer, "#Mode\t0")?;
    writeln!(writer, "#STDev\t0.000")?;
    writeln!(writer, "#PercentOfPairs\t{percent:.3}")?;
    writeln!(writer, "#InsertSize\tCount")?;
    writer.flush()?;
    Ok(())
}

fn write_quality_accuracy_fallback_hist(
    path: &Path,
    hist: &AlignmentFallbackHistograms,
    config: &Config,
) -> Result<()> {
    let mut writer = crate::seqio::create_output(path, config.overwrite || config.append)
        .with_context(|| format!("creating quality-accuracy histogram {}", path.display()))?;
    writeln!(writer, "#Deviation\t0.000")?;
    writeln!(writer, "#DeviationSub\t0.000")?;
    writeln!(writer, "#Avg_STDev\t0.000")?;
    writeln!(writer, "#Diversity\t0.000")?;
    writeln!(writer, "#Entropy\t0.000")?;
    writeln!(
        writer,
        "#Quality\tMatch\tSub\tIns\tDel\tTrueQuality\tTrueQualitySub"
    )?;

    let mut remaining: u64 = hist.quality_match.iter().sum();
    for (quality, matches) in hist.quality_match.iter().copied().enumerate() {
        if matches > 0 || config.print_zero_coverage {
            writeln!(writer, "{quality}\t{matches}\t0\t0\t0\t\t")?;
        }
        remaining = remaining.saturating_sub(matches);
        if remaining == 0 && !config.print_zero_coverage {
            break;
        }
    }
    writer.flush()?;
    Ok(())
}

fn write_indel_fallback_hist(path: &Path, config: &Config) -> Result<()> {
    let mut writer = crate::seqio::create_output(path, config.overwrite || config.append)
        .with_context(|| format!("creating indel histogram {}", path.display()))?;
    writeln!(writer, "#Length\tDeletions\tInsertions")?;
    if config.print_zero_coverage {
        writeln!(writer, "0\t0\t0")?;
    }
    writer.flush()?;
    Ok(())
}

fn write_error_fallback_hist(
    path: &Path,
    hist: &AlignmentFallbackHistograms,
    config: &Config,
) -> Result<()> {
    let mut writer = crate::seqio::create_output(path, config.overwrite || config.append)
        .with_context(|| format!("creating error histogram {}", path.display()))?;
    writeln!(writer, "#Errors\tCount")?;
    if hist.read_count > 0 || config.print_zero_coverage {
        writeln!(writer, "0\t{}", hist.read_count)?;
    }
    writer.flush()?;
    Ok(())
}

fn write_barcode_stats(
    path: &Path,
    barcodes: &BTreeMap<String, u64>,
    config: &Config,
) -> Result<()> {
    let mut writer = crate::seqio::create_output(path, config.overwrite || config.append)
        .with_context(|| format!("creating barcode stats {}", path.display()))?;
    let total: u64 = barcodes.values().copied().sum();
    writeln!(writer, "#Reads\t{total}")?;
    writeln!(writer, "#Barcodes\t{}", barcodes.len())?;

    let mut sorted: Vec<_> = barcodes.iter().collect();
    sorted.sort_by(|(left_name, left_count), (right_name, right_count)| {
        right_count
            .cmp(left_count)
            .then_with(|| left_name.cmp(right_name))
    });
    for (barcode, count) in sorted {
        writeln!(writer, "{barcode}\t{count}")?;
    }
    writer.flush()?;
    Ok(())
}

fn unique_from_raw(depth: usize, raw: u64) -> u64 {
    if depth < 1 {
        raw
    } else {
        (raw + (depth as u64 / 2)) / depth as u64
    }
}

fn percentile_index(cov_last: usize, percentile: f64) -> usize {
    ((cov_last as f64) * (1.0 - percentile)) as usize
}

fn deterministic_coin(rand: Option<f64>, depth: u64) -> u64 {
    debug_assert!(depth > 0);
    (((rand.unwrap_or(0.0) * depth as f64) as u64) + 1).min(depth)
}

fn non_negative_depth(depth: i64) -> Option<u64> {
    u64::try_from(depth).ok()
}

fn depth_below_min(depth: Option<u64>, min_depth: u64) -> bool {
    depth.is_none_or(|depth| depth < min_depth)
}

fn u64_to_i64_saturating(value: u64) -> i64 {
    i64::try_from(value).unwrap_or(i64::MAX)
}

fn min_option(a: Option<u64>, b: Option<u64>) -> Option<u64> {
    match (a, b) {
        (Some(a), Some(b)) => Some(a.min(b)),
        (Some(a), None) => Some(a),
        (None, Some(b)) => Some(b),
        (None, None) => None,
    }
}

fn max_option(a: Option<u64>, b: Option<u64>) -> Option<u64> {
    match (a, b) {
        (Some(a), Some(b)) => Some(a.max(b)),
        (Some(a), None) => Some(a),
        (None, Some(b)) => Some(b),
        (None, None) => None,
    }
}

fn limit_reached(limit: Option<u64>, reads_seen: u64) -> bool {
    limit.is_some_and(|limit| reads_seen >= limit)
}

fn primary_input_lists(config: &Config) -> Option<InputLists> {
    if config.interleaved {
        return None;
    }
    let input = config.in1.as_ref()?;
    if input.exists() {
        return None;
    }
    let text = input.to_string_lossy();
    if !text.contains(',') {
        return None;
    }
    let first = split_path_list(&text);
    if first.len() <= 1 {
        return None;
    }
    let second = config.in2.as_ref().map(|path| {
        let text = path.to_string_lossy();
        split_path_list(&text)
    });
    Some(InputLists { first, second })
}

fn split_path_list(value: &str) -> Vec<PathBuf> {
    value
        .split(',')
        .filter_map(|part| {
            let trimmed = part.trim();
            (!trimmed.is_empty()).then(|| PathBuf::from(trimmed))
        })
        .collect()
}

fn sequence_settings(config: &Config) -> SequenceSettings {
    SequenceSettings {
        bases: BaseSettings {
            u_to_t: config.u_to_t,
            to_upper_case: config.to_upper_case,
            lower_case_to_n: config.lower_case_to_n,
            dot_dash_x_to_n: config.dot_dash_x_to_n,
            iupac_to_n: config.iupac_to_n,
            fix_junk_and_iupac: config.fix_junk_and_iupac,
            junk_mode: config.junk_mode,
        },
        qualities: QualitySettings {
            input_offset: config.quality_in_offset,
            min_called: config.min_called_quality,
            max_called: config.max_called_quality,
            change_quality: config.change_quality,
        },
    }
}

fn open_sequence_writer(
    path: Option<&Path>,
    overwrite: bool,
    append: bool,
    quality_out_offset: u8,
    fake_quality: u8,
    fasta_wrap: usize,
    gzip_threads: Option<usize>,
) -> Result<Option<SequenceWriter>> {
    path.map(|path| {
        SequenceWriter::from_path_with_append_and_gzip_threads(
            path,
            overwrite,
            append,
            quality_out_offset,
            fake_quality,
            fasta_wrap,
            gzip_threads,
        )
    })
    .transpose()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::kmer::kmers_for_record;
    use crate::seqio::SequenceRecord;
    use std::fs;

    fn record(id: &str, bases: &[u8]) -> SequenceRecord {
        SequenceRecord {
            id: id.to_string(),
            numeric_id: 0,
            bases: bases.to_vec(),
            qualities: Some(vec![b'I'; bases.len()]),
        }
    }

    fn quality_record(id: &str, bases: &[u8], qualities: &[u8]) -> SequenceRecord {
        SequenceRecord {
            id: id.to_string(),
            numeric_id: 0,
            bases: bases.to_vec(),
            qualities: Some(qualities.to_vec()),
        }
    }

    #[test]
    fn gzip_threads_are_split_across_concurrent_gzip_streams() {
        assert_eq!(gzip_threads_for_streams(None, 2), None);
        assert_eq!(gzip_threads_for_streams(Some(1), 2), Some(1));
        assert_eq!(gzip_threads_for_streams(Some(8), 0), Some(8));
        assert_eq!(gzip_threads_for_streams(Some(8), 1), Some(8));
        assert_eq!(gzip_threads_for_streams(Some(8), 2), Some(4));
        assert_eq!(gzip_threads_for_streams(Some(8), 3), Some(2));
        assert_eq!(gzip_threads_for_streams(Some(2), 4), Some(1));

        assert_eq!(
            gzip_threads_for_paths(
                Some(8),
                [
                    Some(Path::new("reads_R1.fq.gz")),
                    Some(Path::new("reads_R2.fq.gz")),
                ],
            ),
            Some(4)
        );
        assert_eq!(
            gzip_threads_for_paths(
                Some(8),
                [
                    Some(Path::new("reads_R1.fq")),
                    Some(Path::new("reads_R2.fq.gz")),
                ],
            ),
            Some(8)
        );
    }

    #[test]
    fn write_depth_hist_folds_zero_bin_without_cloning_input_hist() {
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("hist.tsv");
        let hist = vec![5, 7, 4];
        let config = Config {
            overwrite: true,
            ..Config::default()
        };

        write_depth_hist(&path, &hist, &config).unwrap();

        assert_eq!(hist, vec![5, 7, 4]);
        assert_eq!(
            fs::read_to_string(path).unwrap(),
            "#Depth\tRaw_Count\tUnique_Kmers\n1\t12\t12\n2\t4\t2\n"
        );
    }

    #[test]
    fn write_depth_hist_preserves_zero_bin_when_requested() {
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("hist.tsv");
        let hist = vec![5, 7, 4];
        let config = Config {
            overwrite: true,
            zero_bin: true,
            ..Config::default()
        };

        write_depth_hist(&path, &hist, &config).unwrap();

        assert_eq!(
            fs::read_to_string(path).unwrap(),
            "#Depth\tRaw_Count\tUnique_Kmers\n0\t5\t5\n1\t7\t7\n2\t4\t2\n"
        );
    }

    #[test]
    fn write_sparse_depth_hist_matches_dense_output() {
        let dir = tempfile::tempdir().unwrap();
        let dense_path = dir.path().join("dense.hist.tsv");
        let sparse_path = dir.path().join("sparse.hist.tsv");
        let hist = vec![5, 7, 4];
        let sparse = SparseHist::from_iter([(0, 5), (1, 7), (2, 4)]);
        let config = Config {
            overwrite: true,
            ..Config::default()
        };

        write_depth_hist(&dense_path, &hist, &config).unwrap();
        write_sparse_depth_hist(&sparse_path, &sparse, hist.len(), &config).unwrap();

        assert_eq!(
            fs::read_to_string(sparse_path).unwrap(),
            fs::read_to_string(dense_path).unwrap()
        );
    }

    #[test]
    fn write_sparse_depth_hist_matches_dense_zero_coverage_columns_one() {
        let dir = tempfile::tempdir().unwrap();
        let dense_path = dir.path().join("dense.hist.tsv");
        let sparse_path = dir.path().join("sparse.hist.tsv");
        let hist = vec![0, 0, 6, 0, 4];
        let sparse = SparseHist::from_iter([(2, 6), (4, 4)]);
        let config = Config {
            overwrite: true,
            hist_columns: 1,
            print_zero_coverage: true,
            ..Config::default()
        };

        write_depth_hist(&dense_path, &hist, &config).unwrap();
        write_sparse_depth_hist(&sparse_path, &sparse, hist.len(), &config).unwrap();

        assert_eq!(
            fs::read_to_string(sparse_path).unwrap(),
            fs::read_to_string(dense_path).unwrap()
        );
    }

    #[test]
    fn output_counts_sparse_depth_hist_matches_dense_hist() {
        let hist_len = 5;
        let mut exact = CountMap::default();
        exact.insert(KmerKey::Short(1), 1);
        exact.insert(KmerKey::Short(2), 3);
        exact.insert(KmerKey::Short(3), 9);
        let exact = OutputCounts::Exact(exact);
        assert_eq!(
            sparse_hist_to_dense(&exact.sparse_depth_hist(hist_len), hist_len),
            exact.depth_hist(hist_len)
        );

        let mut packed = PackedCountMinSketch::new(8, 1, 4).unwrap();
        packed.set_cell(0, 1);
        packed.set_cell(1, 2);
        packed.set_cell(2, 9);
        let packed = OutputCounts::Sketch(packed);
        assert_eq!(
            sparse_hist_to_dense(&packed.sparse_depth_hist(hist_len), hist_len),
            packed.depth_hist(hist_len)
        );

        let atomic = AtomicCountMinSketch::new(64, 1).unwrap();
        atomic.add_key_count(&KmerKey::Short(7), 2);
        atomic.add_key_count(&KmerKey::Short(11), 4);
        atomic.add_key_count(&KmerKey::Short(13), 9);
        let atomic = OutputCounts::AtomicSketch(atomic);
        assert_eq!(
            sparse_hist_to_dense(&atomic.sparse_depth_hist(hist_len), hist_len),
            atomic.depth_hist(hist_len)
        );
    }

    #[test]
    fn sparse_peak_dense_trims_trailing_zero_histlen_without_changing_peaks() {
        let dir = tempfile::tempdir().unwrap();
        let dense_path = dir.path().join("dense.peaks.tsv");
        let compact_path = dir.path().join("compact.peaks.tsv");
        let hist_len = 10_000;
        let mut dense = vec![0u64; hist_len];
        dense[18] = 180;
        dense[19] = 380;
        dense[20] = 720;
        dense[21] = 380;
        dense[22] = 180;
        let sparse = SparseHist::from_iter(
            dense
                .iter()
                .copied()
                .enumerate()
                .filter_map(|(depth, raw)| (raw > 0).then_some((depth, raw))),
        );
        let compact = sparse_hist_to_peak_dense(&sparse, hist_len);
        let config = Config {
            overwrite: true,
            k: 5,
            peak_min_height: 1,
            peak_min_volume: 1,
            peak_min_width: 1,
            peak_min_peak: 1,
            peak_max_peak: 100,
            peak_max_count: 8,
            ..Config::default()
        };

        assert!(compact.len() < 128);
        write_peaks(&dense_path, &dense, &config).unwrap();
        write_peaks(&compact_path, &compact, &config).unwrap();

        assert_eq!(
            fs::read_to_string(compact_path).unwrap(),
            fs::read_to_string(dense_path).unwrap()
        );
    }

    #[test]
    fn write_sparse_read_depth_hist_matches_dense_output() {
        let dir = tempfile::tempdir().unwrap();
        let dense_path = dir.path().join("dense.rhist.tsv");
        let sparse_path = dir.path().join("sparse.rhist.tsv");
        let mut dense = ReadDepthHistogram::new(4);
        dense.reads[0] = 5;
        dense.bases[0] = 500;
        dense.reads[1] = 7;
        dense.bases[1] = 700;
        dense.reads[3] = 4;
        dense.bases[3] = 400;
        let mut sparse = SparseReadDepthHist::default();
        sparse.insert(0, (5, 500));
        sparse.insert(1, (7, 700));
        sparse.insert(3, (4, 400));
        let config = Config {
            overwrite: true,
            ..Config::default()
        };

        write_read_depth_hist(&dense_path, &dense, &config).unwrap();
        write_sparse_read_depth_hist(&sparse_path, &sparse, 4, &config).unwrap();

        assert_eq!(
            fs::read_to_string(sparse_path).unwrap(),
            fs::read_to_string(dense_path).unwrap()
        );
    }

    #[test]
    fn write_sparse_read_depth_hist_streams_zero_coverage_without_dense_histogram() {
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("sparse.rhist.tsv");
        let mut sparse = SparseReadDepthHist::default();
        sparse.insert(2, (1, 8));
        let config = Config {
            overwrite: true,
            print_zero_coverage: true,
            ..Config::default()
        };

        write_sparse_read_depth_hist(&path, &sparse, 8, &config).unwrap();

        assert_eq!(
            fs::read_to_string(path).unwrap(),
            "#Depth\tReads\tBases\n0\t0\t0\n1\t0\t0\n2\t1\t8\n"
        );
    }

    #[test]
    fn output_gzip_threads_are_split_across_all_active_output_streams() {
        fn plan(first: Option<&str>, second: Option<&str>) -> OutputPathPlan {
            OutputPathPlan {
                pairs: vec![OutputPathPair {
                    first: first.map(PathBuf::from),
                    second: second.map(PathBuf::from),
                }],
                fanout: false,
            }
        }

        let keep = plan(Some("keep1.fq.gz"), Some("keep2.fq.gz"));
        let toss = plan(Some("toss1.fq.gz"), Some("toss2.fq.gz"));
        let low = plan(Some("low.fq.gz"), None);
        let mid = plan(Some("mid.fq"), None);
        let high = plan(None, None);
        let uncorrected = plan(Some("uncorrected1.fq.gz"), Some("uncorrected2.fq.gz"));

        assert_eq!(
            output_gzip_threads_for_plans(
                Some(8),
                [&keep, &toss, &low, &mid, &high, &uncorrected],
                0
            )
            .unwrap(),
            Some(1)
        );

        assert_eq!(
            output_gzip_threads_for_plans(Some(8), [&keep, &toss], 0).unwrap(),
            Some(2)
        );
    }

    fn write_fastq(path: &Path, records: &[(&str, &[u8], &[u8])]) {
        let mut text = Vec::new();
        for (id, bases, qualities) in records {
            text.extend_from_slice(b"@");
            text.extend_from_slice(id.as_bytes());
            text.extend_from_slice(b"\n");
            text.extend_from_slice(bases);
            text.extend_from_slice(b"\n+\n");
            text.extend_from_slice(qualities);
            text.extend_from_slice(b"\n");
        }
        fs::write(path, text).unwrap();
    }

    fn write_repeated_fastq(
        path: &Path,
        prefix: &str,
        bases: &[u8],
        qualities: &[u8],
        count: usize,
    ) {
        let mut text = Vec::new();
        for index in 1..=count {
            text.extend_from_slice(b"@");
            text.extend_from_slice(format!("{prefix}{index}").as_bytes());
            text.extend_from_slice(b"\n");
            text.extend_from_slice(bases);
            text.extend_from_slice(b"\n+\n");
            text.extend_from_slice(qualities);
            text.extend_from_slice(b"\n");
        }
        fs::write(path, text).unwrap();
    }

    #[test]
    fn exact_counts_remove_duplicate_kmers_per_read() {
        let config = Config {
            k: 3,
            min_quality: 0,
            min_prob: 0.0,
            ..Config::default()
        };
        let mut counts = CountMap::default();
        increment_pair_counts(&config, &mut counts, &record("r1", b"AAAAAA"), None);
        assert_eq!(counts.values().copied().sum::<u64>(), 1);
    }

    #[test]
    fn exact_counts_keep_duplicate_long_kmers_like_java_bbnorm() {
        let config = Config {
            k: 40,
            min_quality: 0,
            min_prob: 0.0,
            ..Config::default()
        };
        let mut counts = CountMap::default();
        let record = record("r1", b"AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA");
        let kmers = kmers_for_record(&record, &config);
        assert!(kmers.len() > 1);
        assert!(kmers.windows(2).all(|pair| pair[0] == pair[1]));

        increment_pair_counts(&config, &mut counts, &record, None);

        assert_eq!(counts.len(), 1);
        assert_eq!(counts.values().copied().sum::<u64>(), kmers.len() as u64);
    }

    #[test]
    fn constrained_count_min_inflates_colliding_counts() {
        let config = Config {
            count_min: crate::cli::CountMinSettings {
                cells: Some(1),
                hashes: Some(2),
                bits: Some(8),
                memory_bytes: None,
            },
            ..Config::default()
        };
        let mut counts = CountMap::default();
        counts.insert(KmerKey::Short(7), 2);
        counts.insert(KmerKey::Short(11), 5);

        apply_count_min_collision_estimates(&config, &mut counts);

        assert_eq!(counts.get(&KmerKey::Short(7)), Some(&7));
        assert_eq!(counts.get(&KmerKey::Short(11)), Some(&7));
    }

    #[test]
    fn constrained_count_min_honors_cell_bit_saturation() {
        let config = Config {
            count_min: crate::cli::CountMinSettings {
                cells: Some(1),
                hashes: Some(1),
                bits: Some(2),
                memory_bytes: None,
            },
            ..Config::default()
        };
        let mut counts = CountMap::default();
        counts.insert(KmerKey::Short(7), 2);
        counts.insert(KmerKey::Short(11), 5);

        apply_count_min_collision_estimates(&config, &mut counts);

        assert_eq!(counts.get(&KmerKey::Short(7)), Some(&3));
        assert_eq!(counts.get(&KmerKey::Short(11)), Some(&3));
    }

    #[test]
    fn constrained_count_min_caps_wide_cells_like_kcountarray() {
        let config = Config {
            count_min: crate::cli::CountMinSettings {
                cells: Some(1),
                hashes: Some(1),
                bits: Some(32),
                memory_bytes: None,
            },
            ..Config::default()
        };
        let mut counts = CountMap::default();
        counts.insert(KmerKey::Short(7), i32::MAX as u64 + 10);
        counts.insert(KmerKey::Short(11), 1);

        apply_count_min_collision_estimates(&config, &mut counts);

        assert_eq!(counts.get(&KmerKey::Short(7)), Some(&(i32::MAX as u64)));
        assert_eq!(counts.get(&KmerKey::Short(11)), Some(&(i32::MAX as u64)));
        assert_eq!(count_min_max_count(31), i32::MAX as u64);
        assert_eq!(count_min_max_count(32), i32::MAX as u64);
        assert_eq!(count_min_max_count(64), i32::MAX as u64);
    }

    #[test]
    fn count_min_budget_guard_rejects_tables_above_safe_memory() {
        let available = 1_000_000usize;
        let safe_budget = safe_explicit_count_min_bytes(available);
        let fitting_cells = safe_budget / 4;
        assert!(
            ensure_count_min_budget_fits_ceiling("main", fitting_cells, 32, safe_budget).is_ok()
        );

        let oversized_cells = safe_budget.div_ceil(4) + 1;
        let err = ensure_count_min_budget_fits_ceiling("main", oversized_cells, 32, safe_budget)
            .unwrap_err()
            .to_string();
        assert!(
            err.contains("above safe memory budget"),
            "unexpected error: {err}"
        );
    }

    #[test]
    fn count_min_budget_guard_respects_configured_memory_below_available_ram() {
        let configured = 1_000_000usize;
        let available = 10_000_000usize;
        let safe_budget = count_min_safe_budget_bytes(Some(configured), Some(available)).unwrap();
        assert_eq!(safe_budget, configured);

        assert!(ensure_count_min_budget_fits_ceiling("main", 250_000, 32, safe_budget).is_ok());

        let cells_that_fit_available_but_not_configured = 250_001usize;
        let err = ensure_count_min_budget_fits_ceiling(
            "main",
            cells_that_fit_available_but_not_configured,
            32,
            safe_budget,
        )
        .unwrap_err()
        .to_string();
        assert!(
            err.contains("above safe memory budget"),
            "unexpected configured-budget error: {err}"
        );
    }

    #[test]
    fn count_min_budget_guard_rejects_size_overflow_before_prime_sizing() {
        let err = count_min_total_bytes(usize::MAX, 32)
            .unwrap_err()
            .to_string();
        assert!(
            err.contains("overflowed"),
            "unexpected overflow error: {err}"
        );
    }

    #[test]
    fn count_min_hash_uses_bbtools_row_rotation_masks() {
        let key = KmerKey::Short(0x1234_5678_9abc_def0);
        let first = count_min_bucket(&key, 0, 1024);
        let second = count_min_bucket(&key, 1, 1024);
        let third = count_min_bucket(&key, 2, 1024);

        assert!(first < 1024);
        assert!(second < 1024);
        assert!(third < 1024);
        assert_ne!(first, second);
        assert_ne!(second, third);

        let row0 = bbtools_mask_hash(raw_kmer_key(&key), 0, BBTOOLS_KCOUNT_ARRAY_FIRST_MASK_SEED);
        let row1 = bbtools_mask_hash(
            row0.rotate_right(BBTOOLS_HASH_BITS),
            1,
            BBTOOLS_KCOUNT_ARRAY_FIRST_MASK_SEED,
        );
        assert_eq!(
            count_min_bucket(&key, 1, 1024),
            KCountArrayLayout::new(1024, 32).bucket(row1)
        );

        let expected = [
            0x575a_4571_d954_c5e8,
            0x12bb_293c_ca33_0af3,
            0x0287_fcd8_b8b4_e1c9,
            0x2b62_7d06_2179_52bb,
            0x6bc1_463c_9db3_e422,
            0x710a_bca5_aeb9_5819,
            0x2487_597d_41ef_8ea1,
            0x653b_8694_aa03_bbf0,
        ];
        assert_eq!(
            &bbtools_hash_masks(BBTOOLS_KCOUNT_ARRAY_FIRST_MASK_SEED)[0][..8],
            expected.as_slice()
        );

        for row in bbtools_hash_masks(BBTOOLS_KCOUNT_ARRAY_FIRST_MASK_SEED) {
            for &mask in row {
                assert_eq!((mask & 0xffff_ffff).count_ones(), 16);
                assert!((15..=16).contains(&(mask >> 32).count_ones()));
                assert_eq!(mask >> 63, 0);
            }
        }
    }

    #[test]
    fn prefilter_and_main_sketches_use_independent_kcountarray_mask_seeds() {
        let config = Config {
            count_min: crate::cli::CountMinSettings {
                cells: Some(512),
                hashes: Some(2),
                bits: Some(32),
                memory_bytes: None,
            },
            prefilter: crate::cli::PrefilterSettings {
                enabled: true,
                force_disabled: false,
                ..Default::default()
            },
            ..Config::default()
        };

        let prefilter = new_prefilter_count_min_sketch(&config).unwrap();
        let main = new_atomic_count_min_sketch_with_mask_seed(
            &config,
            BBTOOLS_KCOUNT_ARRAY_SECOND_MASK_SEED,
        )
        .unwrap();
        let key = KmerKey::Short(0x1234_5678_9abc_def0);

        assert_eq!(
            prefilter.layout.mask_seed,
            BBTOOLS_KCOUNT_ARRAY_FIRST_MASK_SEED
        );
        assert_eq!(main.layout.mask_seed, BBTOOLS_KCOUNT_ARRAY_SECOND_MASK_SEED);
        assert_ne!(
            count_min_bucket_with_layout(&key, 0, prefilter.layout),
            count_min_bucket_with_layout(&key, 0, main.layout)
        );
    }

    #[test]
    fn nondeterministic_input_prefilter_uses_atomic_packed_sketch() {
        let config = Config {
            deterministic: false,
            count_min: crate::cli::CountMinSettings {
                cells: Some(512),
                hashes: Some(3),
                bits: Some(32),
                memory_bytes: None,
            },
            prefilter: crate::cli::PrefilterSettings {
                enabled: true,
                force_disabled: false,
                cells: Some(256),
                hashes: Some(2),
                bits: Some(2),
                memory_bytes: None,
                memory_fraction_micros: None,
            },
            ..Config::default()
        };

        let prefilter = new_input_prefilter_count_min_sketch(&config).unwrap();
        let layout = prefilter.layout_summary("input_prefilter", Some(prefilter.max_count()));

        assert!(matches!(
            prefilter,
            PrefilterCountMinSketch::AtomicPacked(_)
        ));
        assert_eq!(layout.kind, "atomic_packed");
        assert_eq!(layout.bits, 2);
        assert_eq!(layout.hashes, 2);
        assert_eq!(layout.update_mode, "conservative");
    }

    #[test]
    fn nondefault_kcountarray_mask_seeds_are_cached() {
        let seed = BBTOOLS_KCOUNT_ARRAY_FIRST_MASK_SEED + BBTOOLS_KCOUNT_ARRAY_MASK_SEED_STEP * 2;
        let first = bbtools_hash_masks(seed);
        let second = bbtools_hash_masks(seed);
        let third = bbtools_hash_masks(seed + BBTOOLS_KCOUNT_ARRAY_MASK_SEED_STEP);

        assert!(std::ptr::eq(first, second));
        assert!(!std::ptr::eq(first, third));
        assert_ne!(first[0][0], third[0][0]);
    }

    #[test]
    fn countup_prefilter_mask_seed_uses_dedicated_hot_cache() {
        let config = Config {
            count_up: true,
            prefilter: crate::cli::PrefilterSettings {
                enabled: true,
                force_disabled: false,
                ..Default::default()
            },
            count_min: crate::cli::CountMinSettings {
                cells: Some(10_000),
                bits: Some(32),
                ..Default::default()
            },
            ..Config::default()
        };

        let seed = countup_output_mask_seed(&config);
        assert_eq!(seed, BBTOOLS_KCOUNT_ARRAY_THIRD_MASK_SEED);
        assert!(std::ptr::eq(
            bbtools_hash_masks(seed),
            bbtools_hash_masks(BBTOOLS_KCOUNT_ARRAY_THIRD_MASK_SEED)
        ));
    }

    #[test]
    fn kcount_layout_carries_resolved_mask_table_for_bucket_fills() {
        let layout = KCountArrayLayout::new_with_min_arrays_and_mask_seed(
            4096,
            32,
            BBTOOLS_KCOUNT_ARRAY_MIN_ARRAYS,
            BBTOOLS_KCOUNT_ARRAY_THIRD_MASK_SEED,
        );

        assert!(std::ptr::eq(
            layout.masks,
            bbtools_hash_masks(BBTOOLS_KCOUNT_ARRAY_THIRD_MASK_SEED)
        ));
        assert_eq!(layout.mask_seed, BBTOOLS_KCOUNT_ARRAY_THIRD_MASK_SEED);
    }

    #[test]
    fn incremental_count_min_buckets_match_row_hash_replay() {
        let layout = KCountArrayLayout::new_with_min_arrays_and_mask_seed(
            4096,
            32,
            BBTOOLS_KCOUNT_ARRAY_MIN_ARRAYS,
            BBTOOLS_KCOUNT_ARRAY_SECOND_MASK_SEED,
        );
        for raw in [0, 1, 7, 31, 63, 255, 0x1234_5678_9abc_def0] {
            let key = KmerKey::Short(raw);
            let mut slots = [usize::MAX; 16];
            fill_count_min_buckets(&key, 8, layout, &mut slots);

            for (hash_index, slot) in slots.iter().enumerate().take(8) {
                assert_eq!(
                    *slot,
                    count_min_bucket_with_layout(&key, hash_index, layout)
                );
            }
        }
    }

    fn find_partial_row_collision(
        cells: usize,
        bits: u8,
    ) -> (KmerKey, KmerKey, usize, usize, usize) {
        let layout = KCountArrayLayout::new(cells, bits);
        let mut seen: Vec<Option<(KmerKey, usize)>> = vec![None; cells];
        for raw in 0..100_000u64 {
            let key = KmerKey::Short(raw);
            let row0 = count_min_bucket_with_layout(&key, 0, layout);
            let row1 = count_min_bucket_with_layout(&key, 1, layout);
            if let Some((previous, previous_row1)) = &seen[row0] {
                if *previous_row1 != row1 {
                    return (previous.clone(), key, row0, *previous_row1, row1);
                }
            } else {
                seen[row0] = Some((key, row1));
            }
        }
        panic!("expected to find a partial row collision for {cells} cells");
    }

    fn find_two_sided_partial_collisions(cells: usize, bits: u8) -> (KmerKey, KmerKey, KmerKey) {
        let layout = KCountArrayLayout::new(cells, bits);
        let base = KmerKey::Short(0);
        let base_row0 = count_min_bucket_with_layout(&base, 0, layout);
        let base_row1 = count_min_bucket_with_layout(&base, 1, layout);
        let mut row0_match = None;
        let mut row1_match = None;
        for raw in 1..200_000u64 {
            let key = KmerKey::Short(raw);
            let row0 = count_min_bucket_with_layout(&key, 0, layout);
            let row1 = count_min_bucket_with_layout(&key, 1, layout);
            if row0 == base_row0 && row1 != base_row1 && row0_match.is_none() {
                row0_match = Some(key.clone());
            }
            if row1 == base_row1 && row0 != base_row0 && row1_match.is_none() {
                row1_match = Some(key);
            }
            if let (Some(row0_match), Some(row1_match)) = (row0_match.clone(), row1_match.clone()) {
                return (base, row0_match, row1_match);
            }
        }
        panic!("expected to find two-sided partial row collisions for {cells} cells");
    }

    #[test]
    fn prefilter_sketch_defaults_to_kcountarray_locked_updates() {
        let config = Config {
            prefilter: crate::cli::PrefilterSettings {
                enabled: true,
                force_disabled: false,
                cells: Some(128),
                hashes: Some(2),
                bits: Some(2),
                memory_bytes: None,
                memory_fraction_micros: None,
            },
            threads: Some(2),
            ..Config::default()
        };
        let mut prefilter = new_prefilter_count_min_sketch(&config).unwrap();
        assert_eq!(prefilter.update_mode, CountMinUpdateMode::Conservative);
        let (left, right, row0, _, _) = find_partial_row_collision(prefilter.cells, prefilter.bits);

        prefilter.add_key_count(&left, 2);
        prefilter.add_key_count(&right, 1);

        assert_eq!(prefilter.cell(row0), 2);
    }

    #[test]
    fn lockedincrement_false_uses_independent_row_increments() {
        let config = Config {
            prefilter: crate::cli::PrefilterSettings {
                enabled: true,
                force_disabled: false,
                cells: Some(128),
                hashes: Some(2),
                bits: Some(2),
                memory_bytes: None,
                memory_fraction_micros: None,
            },
            locked_increment: Some(false),
            threads: Some(2),
            ..Config::default()
        };
        let mut unlocked = new_prefilter_count_min_sketch(&config).unwrap();
        assert_eq!(unlocked.update_mode, CountMinUpdateMode::Independent);
        let (left, right, row0, row1_left, row1_right) =
            find_partial_row_collision(unlocked.cells, unlocked.bits);

        let mut locked =
            PackedCountMinSketch::new(unlocked.cells, unlocked.hashes, unlocked.bits).unwrap();
        locked.add_key_count(&left, 2);
        locked.add_key_count(&right, 1);
        unlocked.add_key_count(&left, 2);
        unlocked.add_key_count(&right, 1);

        assert_eq!(locked.cell(row0), 2);
        assert_eq!(unlocked.cell(row0), 3);
        assert_eq!(unlocked.cell(row1_left), 2);
        assert_eq!(unlocked.cell(row1_right), 1);
    }

    #[test]
    fn atomic_count_min_honors_unlocked_independent_updates() {
        let config = Config {
            count_min: crate::cli::CountMinSettings {
                cells: Some(128),
                hashes: Some(2),
                bits: Some(32),
                memory_bytes: None,
            },
            locked_increment: Some(false),
            threads: Some(2),
            ..Config::default()
        };
        let unlocked = new_atomic_count_min_sketch(&config).unwrap();
        assert_eq!(unlocked.update_mode, CountMinUpdateMode::Independent);
        let (left, right, row0, row1_left, row1_right) =
            find_partial_row_collision(unlocked.cells, 32);

        let locked = AtomicCountMinSketch::new(unlocked.cells, unlocked.hashes).unwrap();
        locked.add_key_count(&left, 2);
        locked.add_key_count(&right, 1);
        unlocked.add_key_count(&left, 2);
        unlocked.add_key_count(&right, 1);

        assert_eq!(locked.cells_by_hash[row0].load(Ordering::Relaxed), 2);
        assert_eq!(unlocked.cells_by_hash[row0].load(Ordering::Relaxed), 3);
        assert_eq!(unlocked.cells_by_hash[row1_left].load(Ordering::Relaxed), 2);
        assert_eq!(
            unlocked.cells_by_hash[row1_right].load(Ordering::Relaxed),
            1
        );
    }

    #[test]
    fn atomic_count_min_allocates_locks_only_for_conservative_updates() {
        let conservative = new_atomic_count_min_sketch(&Config {
            count_min: crate::cli::CountMinSettings {
                cells: Some(128),
                hashes: Some(2),
                bits: Some(32),
                memory_bytes: None,
            },
            ..Config::default()
        })
        .unwrap();
        let independent = new_atomic_count_min_sketch(&Config {
            count_min: crate::cli::CountMinSettings {
                cells: Some(128),
                hashes: Some(2),
                bits: Some(32),
                memory_bytes: None,
            },
            locked_increment: Some(false),
            ..Config::default()
        })
        .unwrap();

        assert_eq!(conservative.locks.len(), BBTOOLS_KCOUNT_ARRAY_LOCKS);
        assert!(independent.locks.is_empty());
    }

    #[test]
    fn atomic_count_min_parallel_replay_requires_nondeterministic_mode() {
        let deterministic = new_atomic_count_min_sketch(&Config {
            count_min: crate::cli::CountMinSettings {
                cells: Some(128),
                hashes: Some(2),
                bits: Some(32),
                memory_bytes: None,
            },
            deterministic: true,
            ..Config::default()
        })
        .unwrap();
        let nondeterministic = new_atomic_count_min_sketch(&Config {
            count_min: crate::cli::CountMinSettings {
                cells: Some(128),
                hashes: Some(2),
                bits: Some(32),
                memory_bytes: None,
            },
            deterministic: false,
            ..Config::default()
        })
        .unwrap();

        assert!(!deterministic.parallel_replay);
        assert!(nondeterministic.parallel_replay);
    }

    #[test]
    fn packed_count_min_increment_returns_previous_min_like_kcountarray() {
        let key = KmerKey::Short(7);
        let mut sketch = PackedCountMinSketch::new(128, 2, 4).unwrap();

        assert_eq!(sketch.increment_and_return_unincremented(&key, 1), 0);
        assert_eq!(sketch.depth(&key), 1);
        assert_eq!(sketch.increment_and_return_unincremented(&key, 3), 1);
        assert_eq!(sketch.depth(&key), 4);
    }

    #[test]
    fn packed_count_min_increment_return_saturates_at_cell_max() {
        let key = KmerKey::Short(11);
        let mut sketch = PackedCountMinSketch::new(1, 2, 2).unwrap();

        assert_eq!(sketch.increment_and_return_unincremented(&key, 10), 0);
        assert_eq!(sketch.depth(&key), 3);
        assert_eq!(sketch.increment_and_return_unincremented(&key, 1), 3);
        assert_eq!(sketch.depth(&key), 3);
    }

    #[test]
    fn atomic_count_min_increment_returns_previous_min_like_kcountarray() {
        let key = KmerKey::Short(13);
        let sketch = AtomicCountMinSketch::new(128, 2).unwrap();

        assert_eq!(sketch.increment_and_return_unincremented(&key, 1), 0);
        assert_eq!(sketch.depth(&key), 1);
        assert_eq!(sketch.increment_and_return_unincremented(&key, 3), 1);
        assert_eq!(sketch.depth(&key), 4);
    }

    #[test]
    fn atomic_packed_count_min_matches_packed_sequential_updates() {
        let keys = [
            (KmerKey::Short(13), 1),
            (KmerKey::Short(29), 2),
            (KmerKey::Short(13), 1),
            (KmerKey::Short(47), 3),
        ];
        let mut packed = PackedCountMinSketch::new_with_min_arrays_and_mask_seed(
            4099,
            3,
            2,
            BBTOOLS_KCOUNT_ARRAY_MIN_ARRAYS,
            BBTOOLS_KCOUNT_ARRAY_FIRST_MASK_SEED,
        )
        .unwrap();
        let atomic = AtomicPackedCountMinSketch::new_with_min_arrays_and_update_mode(
            4099,
            3,
            2,
            BBTOOLS_KCOUNT_ARRAY_MIN_ARRAYS,
            CountMinUpdateMode::Conservative,
            BBTOOLS_KCOUNT_ARRAY_FIRST_MASK_SEED,
        )
        .unwrap();

        for (key, count) in &keys {
            packed.add_key_count(key, *count);
            atomic.add_key_count(key, *count);
        }
        let key_increments = keys.iter().map(|(_, count)| *count).sum();
        packed.add_key_increments(key_increments);
        atomic.add_key_increments(key_increments);

        for slot in 0..packed.cells {
            assert_eq!(atomic.cell(slot), packed.cell(slot));
        }
        let occupied = (0..packed.cells)
            .filter(|&slot| packed.cell(slot) > 0)
            .count();
        assert_eq!(atomic.occupied_slots_at_least(1), occupied);
        assert_eq!(atomic.unique_kmers(), packed.unique_kmers());
    }

    #[test]
    fn atomic_count_min_conservative_updates_are_key_locked_like_kcountarray() {
        let key = KmerKey::Short(13);
        let pool = rayon::ThreadPoolBuilder::new()
            .num_threads(4)
            .build()
            .unwrap();

        pool.install(|| {
            let sketch = AtomicCountMinSketch::new(128, 3).unwrap();

            (0..10_000u64)
                .into_par_iter()
                .for_each(|_| sketch.add_key_count(&key, 1));

            assert_eq!(sketch.depth(&key), 10_000);
        });
    }

    #[test]
    fn atomic_count_min_bulk_replay_matches_locked_sequential_updates() {
        let mut counts = CountMap::default();
        counts.insert(KmerKey::Short(13), 17);
        counts.insert(KmerKey::Short(29), 3);
        counts.insert(KmerKey::Short(31), 9);
        let locked = AtomicCountMinSketch::new(128, 3).unwrap();
        let bulk = AtomicCountMinSketch::new(128, 3).unwrap();

        for (key, count) in &counts {
            locked.add_key_count(key, *count);
        }
        bulk.add_key_counts(&counts);

        for slot in 0..locked.cells {
            assert_eq!(
                locked.cells_by_hash[slot].load(Ordering::Relaxed),
                bulk.cells_by_hash[slot].load(Ordering::Relaxed)
            );
        }
    }

    #[test]
    fn packed_count_min_reduced_sorted_replay_matches_individual_kmer_updates() {
        let keys = [
            KmerKey::Short(13),
            KmerKey::Short(29),
            KmerKey::Short(13),
            KmerKey::Short(31),
            KmerKey::Short(29),
            KmerKey::Short(29),
            KmerKey::Short(47),
        ];
        let mut individual = PackedCountMinSketch::new(4099, 3, 16).unwrap();
        let mut reduced = PackedCountMinSketch::new(4099, 3, 16).unwrap();

        for key in &keys {
            individual.increment(key);
        }
        for (key, count) in sorted_reduced_test_runs(keys) {
            reduced.add_key_count(&key, count);
            reduced.add_key_increments(count);
        }

        assert_eq!(reduced.increments, individual.increments);
        assert_eq!(reduced.occupied_slots, individual.occupied_slots);
        assert_eq!(reduced.words, individual.words);
    }

    #[test]
    fn atomic_count_min_reduced_sorted_replay_matches_individual_kmer_updates() {
        let keys = [
            KmerKey::Short(13),
            KmerKey::Short(29),
            KmerKey::Short(13),
            KmerKey::Short(31),
            KmerKey::Short(29),
            KmerKey::Short(29),
            KmerKey::Short(47),
        ];
        let individual = AtomicCountMinSketch::new(4099, 3).unwrap();
        let reduced = AtomicCountMinSketch::new(4099, 3).unwrap();

        for key in &keys {
            individual.increment_key(key);
            individual.add_key_increments(1);
        }
        for (key, count) in sorted_reduced_test_runs(keys) {
            reduced.add_key_count(&key, count);
            reduced.add_key_increments(count);
        }

        assert_eq!(
            reduced.increments.load(Ordering::Relaxed),
            individual.increments.load(Ordering::Relaxed)
        );
        assert_eq!(
            reduced.occupied_slots.load(Ordering::Relaxed),
            individual.occupied_slots.load(Ordering::Relaxed)
        );
        for slot in 0..individual.cells {
            assert_eq!(
                reduced.cells_by_hash[slot].load(Ordering::Relaxed),
                individual.cells_by_hash[slot].load(Ordering::Relaxed)
            );
        }
    }

    fn sorted_reduced_test_runs<const N: usize>(keys: [KmerKey; N]) -> Vec<(KmerKey, u64)> {
        let mut keys = keys;
        keys.sort_unstable();
        let mut runs = Vec::new();
        for key in keys {
            if let Some((last_key, count)) = runs.last_mut()
                && last_key == &key
            {
                *count += 1;
                continue;
            }
            runs.push((key, 1));
        }
        runs
    }

    #[test]
    fn exact_collision_estimates_follow_lockedincrement_mode() {
        let mut config = Config {
            count_min: crate::cli::CountMinSettings {
                cells: Some(128),
                hashes: Some(2),
                bits: Some(8),
                memory_bytes: None,
            },
            threads: Some(2),
            ..Config::default()
        };
        let cells = count_min_table_cells_from_total_bits(128, 8);
        let (left, right0, right1) = find_two_sided_partial_collisions(cells, 8);
        let mut locked = CountMap::default();
        locked.insert(left.clone(), 2);
        locked.insert(right0, 1);
        locked.insert(right1, 1);
        let mut unlocked = locked.clone();

        apply_count_min_collision_estimates(&config, &mut locked);
        config.locked_increment = Some(false);
        apply_count_min_collision_estimates(&config, &mut unlocked);

        assert_eq!(locked.get(&left), Some(&2));
        assert_eq!(unlocked.get(&left), Some(&3));
    }

    #[test]
    fn prefilter_exact_estimates_follow_lockedincrement_mode() {
        let mut config = Config {
            prefilter: crate::cli::PrefilterSettings {
                enabled: true,
                force_disabled: false,
                cells: Some(128),
                hashes: Some(2),
                bits: Some(8),
                memory_bytes: None,
                memory_fraction_micros: None,
            },
            threads: Some(2),
            ..Config::default()
        };
        let cells = count_min_table_cells_from_total_bits(128, 8);
        let (left, right0, right1) = find_two_sided_partial_collisions(cells, 8);
        let mut locked = CountMap::default();
        locked.insert(left.clone(), 2);
        locked.insert(right0, 1);
        locked.insert(right1, 1);
        let mut unlocked = locked.clone();

        apply_prefilter_collision_estimates(&config, &mut locked);
        config.locked_increment = Some(false);
        apply_prefilter_collision_estimates(&config, &mut unlocked);

        assert_eq!(locked.get(&left), Some(&2));
        assert_eq!(unlocked.get(&left), Some(&3));
    }

    #[test]
    fn prefilter_sketch_saturates_with_independent_row_increments_when_unlocked() {
        let config = Config {
            prefilter: crate::cli::PrefilterSettings {
                enabled: true,
                force_disabled: false,
                cells: Some(128),
                hashes: Some(2),
                bits: Some(2),
                memory_bytes: None,
                memory_fraction_micros: None,
            },
            locked_increment: Some(false),
            threads: Some(2),
            ..Config::default()
        };
        let mut prefilter = new_prefilter_count_min_sketch(&config).unwrap();
        let (left, right, row0, row1_left, row1_right) =
            find_partial_row_collision(prefilter.cells, prefilter.bits);

        let mut conservative =
            PackedCountMinSketch::new(prefilter.cells, prefilter.hashes, prefilter.bits).unwrap();
        conservative.add_key_count(&left, 2);
        conservative.add_key_count(&right, 1);
        prefilter.add_key_count(&left, 2);
        prefilter.add_key_count(&right, 1);

        assert_eq!(conservative.cell(row0), 2);
        assert_eq!(prefilter.cell(row0), 3);
        assert_eq!(prefilter.cell(row1_left), 2);
        assert_eq!(prefilter.cell(row1_right), 1);
    }

    #[test]
    fn packed_count_min_sketch_uses_fixed_saturating_cells() {
        let mut sketch = PackedCountMinSketch::new(1, 2, 3).unwrap();
        for _ in 0..10 {
            sketch.increment(&KmerKey::Short(7));
        }

        assert_eq!(sketch.words.len(), 1);
        assert_eq!(sketch.depth(&KmerKey::Short(7)), 7);
        assert_eq!(sketch.depth(&KmerKey::Short(11)), 7);
        assert_eq!(sketch.unique_kmers(), 10);
    }

    #[test]
    fn packed_count_min_depth_hist_uses_raw_depth_counts() {
        let mut sketch = PackedCountMinSketch::new(8, 2, 4).unwrap();
        sketch.set_cell(0, 1);
        sketch.set_cell(1, 2);
        sketch.set_cell(2, 2);
        sketch.set_cell(3, 5);

        assert_eq!(sketch.occupied_slots_at_least(1), 4);
        assert_eq!(sketch.tracked_slots.as_ref().unwrap().len(), 4);
        assert_eq!(sketch.depth_hist(4), vec![0, 1, 4, 5]);
    }

    #[test]
    fn packed_count_min_tracks_occupied_slots_without_duplicates() {
        let key = KmerKey::Short(17);
        let mut sketch = PackedCountMinSketch::new(128, 1, 4).unwrap();

        sketch.add_key_count(&key, 1);
        sketch.add_key_count(&key, 2);

        assert_eq!(sketch.occupied_slots_at_least(1), 1);
        assert_eq!(sketch.occupied_slots_at_least(3), 1);
        assert_eq!(sketch.tracked_slots.as_ref().unwrap().len(), 1);
        assert_eq!(sketch.depth_hist(5), vec![0, 0, 0, 3, 0]);
    }

    #[test]
    fn packed_count_min_disables_slot_tracking_for_large_tables() {
        let sketch = PackedCountMinSketch::new(PACKED_SKETCH_TRACKED_SLOT_LIMIT + 1, 1, 1).unwrap();

        assert!(sketch.tracked_slots.is_none());
        assert_eq!(sketch.tracked_slot_memory_bytes(), 0);
        assert_eq!(
            sketch.layout_summary("large", None).memory_bytes,
            sketch.words.len() * std::mem::size_of::<u64>()
        );
    }

    #[test]
    fn packed_count_min_layout_reports_tracked_slot_memory() {
        let key = KmerKey::Short(17);
        let mut sketch = PackedCountMinSketch::new(128, 1, 4).unwrap();

        sketch.add_key_count(&key, 1);

        let backing_bytes = sketch.words.len() * std::mem::size_of::<u64>();
        assert!(sketch.tracked_slot_memory_bytes() >= std::mem::size_of::<usize>());
        assert_eq!(
            sketch.layout_summary("small", None).memory_bytes,
            backing_bytes + sketch.tracked_slot_memory_bytes()
        );
    }

    #[test]
    fn packed_count_min_depth_hist_uses_compact_cell_bound_but_returns_requested_len() {
        let mut sketch = PackedCountMinSketch::new(16, 1, 4).unwrap();
        sketch.set_cell(0, 1);
        sketch.set_cell(1, 15);

        let hist = sketch.depth_hist(1024);

        assert_eq!(hist.len(), 1024);
        assert_eq!(hist[1], 1);
        assert_eq!(hist[15], 15);
        assert!(hist[16..].iter().all(|&value| value == 0));
    }

    #[test]
    fn packed_count_min_untracked_depth_hist_uses_compact_reducers() {
        let mut sketch = PackedCountMinSketch::new(16, 1, 4).unwrap();
        sketch.tracked_slots = None;
        sketch.set_cell(0, 1);
        sketch.set_cell(1, 15);

        let hist = sketch.depth_hist(1024);

        assert_eq!(hist.len(), 1024);
        assert_eq!(hist[1], 1);
        assert_eq!(hist[15], 15);
        assert!(hist[16..].iter().all(|&value| value == 0));
    }

    #[test]
    fn packed_count_min_depth_hist_uses_dynamic_reducers_for_wide_cells() {
        let mut sketch = PackedCountMinSketch::new(16, 1, 32).unwrap();
        sketch.set_cell(0, 1);
        sketch.set_cell(1, 4096);

        let hist = sketch.depth_hist(8192);

        assert_eq!(hist.len(), 8192);
        assert_eq!(hist[1], 1);
        assert_eq!(hist[4096], 4096);
        assert!(hist[4097..].iter().all(|&value| value == 0));
    }

    #[test]
    fn packed_count_min_untracked_depth_hist_uses_dynamic_reducers_for_wide_cells() {
        let mut sketch = PackedCountMinSketch::new(16, 1, 32).unwrap();
        sketch.tracked_slots = None;
        sketch.set_cell(0, 2);
        sketch.set_cell(1, 4096);

        let hist = sketch.depth_hist(8192);

        assert_eq!(hist.len(), 8192);
        assert_eq!(hist[2], 2);
        assert_eq!(hist[4096], 4096);
        assert!(hist[4097..].iter().all(|&value| value == 0));
    }

    #[test]
    fn atomic_count_min_depth_hist_uses_raw_depth_counts() {
        let sketch = AtomicCountMinSketch::new(8, 2).unwrap();
        sketch.cells_by_hash[0].store(1, Ordering::Relaxed);
        sketch.cells_by_hash[1].store(2, Ordering::Relaxed);
        sketch.cells_by_hash[2].store(2, Ordering::Relaxed);
        sketch.cells_by_hash[3].store(5, Ordering::Relaxed);

        assert_eq!(sketch.depth_hist(4), vec![0, 1, 4, 5]);
    }

    #[test]
    fn atomic_count_min_depth_hist_uses_compact_dynamic_reducers() {
        let sketch = AtomicCountMinSketch::new(16, 2).unwrap();
        sketch.cells_by_hash[0].store(1, Ordering::Relaxed);
        sketch.cells_by_hash[1].store(7, Ordering::Relaxed);

        let hist = sketch.depth_hist(8192);

        assert_eq!(hist.len(), 8192);
        assert_eq!(hist[1], 1);
        assert_eq!(hist[7], 7);
        assert!(hist[8..].iter().all(|&value| value == 0));
    }

    #[test]
    fn combined_primary_histograms_match_separate_collectors() {
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("reads.fq");
        write_fastq(
            &path,
            &[
                ("r1", b"ACGTACGT", b"IIIIIIII"),
                ("r2", b"ACGTTCGT", b"IIIIIIII"),
                ("r3", b"TTTTACGT", b"IIIIIIII"),
            ],
        );
        let config = Config {
            in1: Some(path.clone()),
            k: 3,
            min_quality: 0,
            min_prob: 0.0,
            ..Config::default()
        };
        let mut counts = CountMap::default();
        count_single_file(&config, &path, &mut counts, None).unwrap();

        let separate_hist = collect_primary_hist(&config, &counts, None, 0).unwrap();
        let sparse_hist = collect_primary_sparse_hist(&config, &counts, None, 0).unwrap();
        let separate_rhist = collect_primary_read_hist(&config, &counts, None, 0).unwrap();
        let sparse_rhist = collect_primary_sparse_read_hist(&config, &counts, None, 0).unwrap();
        let (sparse_combined_hist, sparse_combined_rhist) =
            collect_primary_sparse_hist_and_read_hist(&config, &counts, None, 0).unwrap();
        let (combined_hist, combined_rhist) =
            collect_primary_hist_and_read_hist(&config, &counts, None, 0).unwrap();

        assert_eq!(
            sparse_hist_to_dense(&sparse_hist, config.hist_len),
            separate_hist
        );
        assert_eq!(
            sparse_hist_to_dense(&sparse_combined_hist, config.hist_len),
            separate_hist
        );
        assert_eq!(combined_hist, separate_hist);
        assert_eq!(combined_rhist.reads, separate_rhist.reads);
        assert_eq!(combined_rhist.bases, separate_rhist.bases);
        let mut dense_sparse_rhist = ReadDepthHistogram::new(config.hist_len);
        merge_sparse_read_depth_hist_into_dense(&mut dense_sparse_rhist, sparse_rhist);
        assert_eq!(dense_sparse_rhist.reads, separate_rhist.reads);
        assert_eq!(dense_sparse_rhist.bases, separate_rhist.bases);
        let mut dense_sparse_combined_rhist = ReadDepthHistogram::new(config.hist_len);
        merge_sparse_read_depth_hist_into_dense(
            &mut dense_sparse_combined_rhist,
            sparse_combined_rhist,
        );
        assert_eq!(dense_sparse_combined_rhist.reads, separate_rhist.reads);
        assert_eq!(dense_sparse_combined_rhist.bases, separate_rhist.bases);
    }

    #[test]
    fn countup_work_source_collects_input_histograms_like_separate_collectors() {
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("reads.fq");
        write_fastq(
            &path,
            &[
                ("r1", b"ACGTACGT", b"IIIIIIII"),
                ("r2", b"ACGTTCGT", b"IIIIIIII"),
                ("r3", b"TTTTACGT", b"IIIIIIII"),
            ],
        );
        let config = Config {
            in1: Some(path.clone()),
            count_up: true,
            k: 3,
            min_quality: 0,
            min_prob: 0.0,
            hist_len: 64,
            ..Config::default()
        };
        let mut counts = CountMap::default();
        count_single_file(&config, &path, &mut counts, None).unwrap();

        let separate_hist = collect_primary_hist(&config, &counts, None, 0).unwrap();
        let separate_rhist = collect_primary_read_hist(&config, &counts, None, 0).unwrap();
        let build = collect_countup_work_source(&config, &counts, 0, true, true).unwrap();

        assert_eq!(build.format1, SeqFormat::Fastq);
        assert_eq!(build.format2, None);
        assert_eq!(
            sparse_hist_to_dense(&build.input_hist.unwrap(), config.hist_len),
            separate_hist
        );
        let mut combined_rhist = ReadDepthHistogram::new(config.hist_len);
        merge_sparse_read_depth_hist_into_dense(
            &mut combined_rhist,
            build.input_read_hist.unwrap(),
        );
        assert_eq!(combined_rhist.reads, separate_rhist.reads);
        assert_eq!(combined_rhist.bases, separate_rhist.bases);
    }

    #[test]
    fn combined_primary_histograms_with_keep_filter_match_separate_collectors() {
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("reads.fq");
        write_fastq(
            &path,
            &[
                ("r1", b"ACGTACGT", b"IIIIIIII"),
                ("r2", b"ACGTACGT", b"IIIIIIII"),
                ("r3", b"TTTTACGT", b"IIIIIIII"),
            ],
        );
        let config = Config {
            in1: Some(path.clone()),
            k: 3,
            min_quality: 0,
            min_prob: 0.0,
            ..Config::default()
        };
        let mut input_counts = CountMap::default();
        count_single_file(&config, &path, &mut input_counts, None).unwrap();
        let mut kept_counts = CountMap::default();
        increment_pair_counts(
            &config,
            &mut kept_counts,
            &record("kept", b"ACGTACGT"),
            None,
        );

        let separate_hist =
            collect_primary_hist(&config, &kept_counts, Some(&input_counts), 17).unwrap();
        let sparse_hist =
            collect_primary_sparse_hist(&config, &kept_counts, Some(&input_counts), 17).unwrap();
        let separate_rhist =
            collect_primary_read_hist(&config, &kept_counts, Some(&input_counts), 17).unwrap();
        let sparse_rhist =
            collect_primary_sparse_read_hist(&config, &kept_counts, Some(&input_counts), 17)
                .unwrap();
        let (sparse_combined_hist, sparse_combined_rhist) =
            collect_primary_sparse_hist_and_read_hist(
                &config,
                &kept_counts,
                Some(&input_counts),
                17,
            )
            .unwrap();
        let (combined_hist, combined_rhist) =
            collect_primary_hist_and_read_hist(&config, &kept_counts, Some(&input_counts), 17)
                .unwrap();

        assert_eq!(
            sparse_hist_to_dense(&sparse_hist, config.hist_len),
            separate_hist
        );
        assert_eq!(
            sparse_hist_to_dense(&sparse_combined_hist, config.hist_len),
            separate_hist
        );
        assert_eq!(combined_hist, separate_hist);
        assert_eq!(combined_rhist.reads, separate_rhist.reads);
        assert_eq!(combined_rhist.bases, separate_rhist.bases);
        let mut dense_sparse_rhist = ReadDepthHistogram::new(config.hist_len);
        merge_sparse_read_depth_hist_into_dense(&mut dense_sparse_rhist, sparse_rhist);
        assert_eq!(dense_sparse_rhist.reads, separate_rhist.reads);
        assert_eq!(dense_sparse_rhist.bases, separate_rhist.bases);
        let mut dense_sparse_combined_rhist = ReadDepthHistogram::new(config.hist_len);
        merge_sparse_read_depth_hist_into_dense(
            &mut dense_sparse_combined_rhist,
            sparse_combined_rhist,
        );
        assert_eq!(dense_sparse_combined_rhist.reads, separate_rhist.reads);
        assert_eq!(dense_sparse_combined_rhist.bases, separate_rhist.bases);
    }

    #[test]
    fn packed_count_min_unique_kmers_uses_bbtools_hash_adjusted_estimate() {
        let mut sketch = PackedCountMinSketch::new(1024, 4, 8).unwrap();
        for bucket in 0..256 {
            sketch.set_cell(bucket, 1);
        }
        sketch.increments = 1_000;

        let estimated = sketch.unique_kmers();
        assert!(
            (70..=80).contains(&estimated),
            "BBTools-style hash-adjusted estimate was {estimated}"
        );
    }

    #[test]
    fn packed_count_min_unique_kmers_honors_min_depth_threshold() {
        let mut sketch = PackedCountMinSketch::new(1024, 4, 8).unwrap();
        for bucket in 0..256 {
            let depth = if bucket < 128 { 3 } else { 1 };
            sketch.set_cell(bucket, depth);
        }
        sketch.increments = 1_000;

        let total_estimated = sketch.unique_kmers();
        let high_depth_estimated = sketch.unique_kmers_at_least(2);

        assert!(
            (70..=80).contains(&total_estimated),
            "all-depth estimate was {total_estimated}"
        );
        assert!(
            (30..=40).contains(&high_depth_estimated),
            "thresholded estimate was {high_depth_estimated}"
        );
        assert_eq!(sketch.unique_kmers_at_least(9), 0);
    }

    #[test]
    fn atomic_count_min_unique_kmers_honors_min_depth_threshold() {
        let sketch = AtomicCountMinSketch::new(1024, 4).unwrap();
        for bucket in 0..256 {
            let depth = if bucket < 128 { 3 } else { 1 };
            sketch.cells_by_hash[bucket].store(depth, Ordering::Relaxed);
        }
        sketch.occupied_slots.store(256, Ordering::Relaxed);
        sketch.add_key_increments(1_000);

        let total_estimated = sketch.unique_kmers();
        let high_depth_estimated = sketch.unique_kmers_at_least(2);

        assert!(
            (70..=80).contains(&total_estimated),
            "all-depth estimate was {total_estimated}"
        );
        assert!(
            (30..=40).contains(&high_depth_estimated),
            "thresholded estimate was {high_depth_estimated}"
        );
        assert_eq!(sketch.occupied_slots_at_least(1), 256);
    }

    #[test]
    fn cardinality_estimator_tracks_unique_keys_with_fixed_register_memory() {
        let config = Config {
            k: 31,
            cardinality: crate::cli::CardinalitySettings {
                input: true,
                buckets: 2048,
                seed: 42,
                ..Default::default()
            },
            ..Default::default()
        };
        let mut estimator = KmerCardinalityEstimator::from_config(&config);
        for key in 0..1_000 {
            estimator.observe_key(&KmerKey::Short(key));
            estimator.observe_key(&KmerKey::Short(key));
        }

        let estimate = estimator.estimate();
        assert_eq!(estimate.k, 31);
        assert_eq!(estimate.buckets, 2048);
        assert!(
            (900..=1_100).contains(&estimate.estimated_unique_kmers),
            "cardinality estimate was {}",
            estimate.estimated_unique_kmers
        );
        assert_eq!(estimator.registers.len(), 2048);
    }

    #[test]
    fn packed_count_min_sketch_packs_cells_across_word_boundaries() {
        let mut sketch = PackedCountMinSketch::new(17, 1, 5).unwrap();
        for slot in 0..17 {
            sketch.set_cell(slot, slot as u64);
        }

        for slot in 0..17 {
            assert_eq!(sketch.cell(slot), slot as u64);
        }
    }

    #[test]
    fn bounded_input_counts_builds_direct_sketch_when_cells_are_constrained() {
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("reads.fq");
        write_fastq(
            &path,
            &[
                ("r1", b"ACGTACGT", b"IIIIIIII"),
                ("r2", b"ACGTTCGT", b"IIIIIIII"),
            ],
        );
        let config = Config {
            in1: Some(path),
            k: 3,
            min_quality: 0,
            min_prob: 0.0,
            count_min: crate::cli::CountMinSettings {
                cells: Some(4),
                hashes: Some(2),
                bits: Some(4),
                memory_bytes: None,
            },
            ..Config::default()
        };
        let probe = kmers_for_record(&record("probe", b"ACGTACGT"), &config)
            .into_iter()
            .next()
            .unwrap();

        let counts = build_input_counts(&config).unwrap();

        let InputCounts::Sketch(sketch) = counts else {
            panic!("cells= should build a bounded packed count-min sketch");
        };
        assert_eq!(sketch.words.len(), 1);
        assert!(sketch.depth(&probe) > 0);
    }

    #[test]
    fn auto_count_min_uses_sketch_when_input_metadata_exceeds_threshold() {
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("reads.fq");
        write_fastq(
            &path,
            &[
                ("r1", b"ACGTACGT", b"IIIIIIII"),
                ("r2", b"ACGTTCGT", b"IIIIIIII"),
            ],
        );
        let config = Config {
            in1: Some(path),
            k: 3,
            min_quality: 0,
            min_prob: 0.0,
            auto_count_min_input_bytes: 1,
            auto_count_min_memory_bytes: Some(4096),
            ..Config::default()
        };

        let counts = build_input_counts(&config).unwrap();

        match counts {
            InputCounts::AtomicSketch(sketch) => {
                assert!(sketch.cells > 0);
                assert!(sketch.increments.load(Ordering::Relaxed) > 0);
            }
            InputCounts::Sketch(sketch) => {
                assert!(sketch.cells > 0);
                assert!(sketch.increments > 0);
            }
            InputCounts::PrefilteredSketch { .. } => {}
            InputCounts::Exact(_) => {
                panic!("large-input auto count-min should build a bounded sketch");
            }
        }
    }

    #[test]
    fn force_exact_counts_overrides_auto_and_explicit_sketch_settings() {
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("reads.fq");
        write_fastq(
            &path,
            &[
                ("r1", b"ACGTACGT", b"IIIIIIII"),
                ("r2", b"ACGTTCGT", b"IIIIIIII"),
            ],
        );
        let config = Config {
            in1: Some(path),
            k: 3,
            min_quality: 0,
            min_prob: 0.0,
            force_exact_counts: true,
            auto_count_min_input_bytes: 1,
            count_min: crate::cli::CountMinSettings {
                cells: Some(1),
                hashes: Some(2),
                bits: Some(4),
                memory_bytes: Some(1024),
            },
            ..Config::default()
        };

        let counts = build_input_counts(&config).unwrap();

        let InputCounts::Exact(counts) = counts else {
            panic!("force_exact_counts should override automatic and explicit sketch settings");
        };
        assert!(counts.len() > 1);
    }

    #[test]
    fn bounded_sketch_chunked_parallel_is_deterministic_and_conservative() {
        let config = Config {
            k: 3,
            min_quality: 0,
            min_prob: 0.0,
            count_min: crate::cli::CountMinSettings {
                cells: Some(32),
                hashes: Some(3),
                bits: Some(8),
                memory_bytes: None,
            },
            ..Config::default()
        };
        let pairs = vec![
            (
                record("r1/1", b"ACGTACGT"),
                Some(record("r1/2", b"TCGTACGA")),
            ),
            (record("r2/1", b"AAAAACCC"), None),
            (
                record("r3/1", b"GGGGTTTT"),
                Some(record("r3/2", b"CCCCAAAA")),
            ),
        ];
        let mut exact = CountMap::default();
        for (r1, r2) in &pairs {
            increment_pair_counts(&config, &mut exact, r1, r2.as_ref());
        }
        let mut chunked_a = new_bounded_count_min_sketch(&config).unwrap();
        let mut chunked_b = new_bounded_count_min_sketch(&config).unwrap();

        increment_sketch_from_pair_chunk(&config, &mut chunked_a, &pairs, None);
        increment_sketch_from_pair_chunk(&config, &mut chunked_b, &pairs, None);

        assert_eq!(chunked_a.words, chunked_b.words);
        assert_eq!(chunked_a.increments, exact.values().copied().sum::<u64>());
        for (key, exact_depth) in exact {
            assert!(chunked_a.depth(&key) >= exact_depth.min(chunked_a.max_count));
        }
    }

    #[test]
    fn atomic_count_min_chunked_parallel_matches_sequential_conservative_bits32() {
        let config = Config {
            k: 3,
            min_quality: 0,
            min_prob: 0.0,
            count_min: crate::cli::CountMinSettings {
                cells: Some(64),
                hashes: Some(3),
                bits: Some(32),
                memory_bytes: None,
            },
            ..Config::default()
        };
        let pairs = vec![
            (
                record("r1/1", b"ACGTACGT"),
                Some(record("r1/2", b"TCGTACGA")),
            ),
            (record("r2/1", b"AAAAACCC"), None),
            (
                record("r3/1", b"GGGGTTTT"),
                Some(record("r3/2", b"CCCCAAAA")),
            ),
        ];
        let sequential = new_atomic_count_min_sketch(&config).unwrap();
        let mut merged_counts = CountMap::default();
        for (r1, r2) in &pairs {
            let mut pair_counts = CountMap::default();
            increment_pair_counts(&config, &mut pair_counts, r1, r2.as_ref());
            merge_count_maps(&mut merged_counts, pair_counts);
        }
        let key_increments = merged_counts.values().copied().sum();
        sequential.add_key_counts(&merged_counts);
        sequential.add_key_increments(key_increments);
        let chunked = new_atomic_count_min_sketch(&config).unwrap();

        increment_atomic_sketch_from_pair_chunk(&config, &chunked, &pairs, None);

        assert_eq!(
            chunked.increments.load(Ordering::Relaxed),
            sequential.increments.load(Ordering::Relaxed)
        );
        assert_eq!(
            chunked.occupied_slots.load(Ordering::Relaxed),
            sequential.occupied_slots.load(Ordering::Relaxed)
        );
        for slot in 0..sequential.cells {
            assert_eq!(
                u64::from(chunked.cells_by_hash[slot].load(Ordering::Relaxed)),
                u64::from(sequential.cells_by_hash[slot].load(Ordering::Relaxed))
            );
        }
    }

    #[test]
    fn nondeterministic_atomic_count_min_direct_path_matches_sequential_without_collisions() {
        let config = Config {
            k: 5,
            min_quality: 0,
            min_prob: 0.0,
            deterministic: false,
            count_min: crate::cli::CountMinSettings {
                cells: Some(8192),
                hashes: Some(1),
                bits: Some(32),
                memory_bytes: None,
            },
            ..Config::default()
        };
        let pairs = vec![
            (
                record("r1/1", b"ACGTACGTAC"),
                Some(record("r1/2", b"TCGTACGAAA")),
            ),
            (record("r2/1", b"AAAAACCCCC"), None),
            (
                record("r3/1", b"GGGGTTTTAA"),
                Some(record("r3/2", b"CCCCAAAAGG")),
            ),
        ];
        let sequential = new_atomic_count_min_sketch(&Config {
            deterministic: true,
            ..config.clone()
        })
        .unwrap();
        let mut merged_counts = CountMap::default();
        for (r1, r2) in &pairs {
            increment_pair_counts(&config, &mut merged_counts, r1, r2.as_ref());
        }
        let key_increments = merged_counts.values().copied().sum();
        sequential.add_key_counts(&merged_counts);
        sequential.add_key_increments(key_increments);

        let direct = new_atomic_count_min_sketch(&config).unwrap();
        increment_atomic_sketch_from_pair_chunk(&config, &direct, &pairs, None);

        assert_eq!(
            direct.increments.load(Ordering::Relaxed),
            sequential.increments.load(Ordering::Relaxed)
        );
        assert_eq!(
            direct.occupied_slots.load(Ordering::Relaxed),
            sequential.occupied_slots.load(Ordering::Relaxed)
        );
        for slot in 0..sequential.cells {
            assert_eq!(
                u64::from(direct.cells_by_hash[slot].load(Ordering::Relaxed)),
                u64::from(sequential.cells_by_hash[slot].load(Ordering::Relaxed))
            );
        }
    }

    #[test]
    fn atomic_count_min_conservative_update_reduces_collision_inflation() {
        let config = Config {
            k: 3,
            min_quality: 0,
            min_prob: 0.0,
            count_min: crate::cli::CountMinSettings {
                cells: Some(1),
                hashes: Some(3),
                bits: Some(32),
                memory_bytes: None,
            },
            ..Config::default()
        };
        let key_a = KmerKey::Short(1);
        let key_b = KmerKey::Short(2);
        let sketch = new_atomic_count_min_sketch(&config).unwrap();

        sketch.add_key_count(&key_a, 5);
        sketch.add_key_count(&key_b, 1);

        assert_eq!(sketch.depth(&key_a), 6);
        assert_eq!(sketch.depth(&key_b), 6);
    }

    #[test]
    fn bounded_output_counts_uses_sketch_for_kept_kmers_when_cells_are_constrained() {
        let config = Config {
            k: 3,
            min_quality: 0,
            min_prob: 0.0,
            count_min: crate::cli::CountMinSettings {
                cells: Some(4),
                hashes: Some(2),
                bits: Some(4),
                memory_bytes: None,
            },
            ..Config::default()
        };
        let r1 = record("r1", b"ACGTACGT");
        let probe = kmers_for_record(&r1, &config).into_iter().next().unwrap();
        let pair = NormalizedPair {
            input_list_index: 0,
            r1: r1.clone(),
            r2: None,
            out_r1: r1,
            out_r2: None,
            decision: PairDecision::default(),
            uncorrectable: false,
            read_count: 1,
            base_count: 8,
        };
        let mut counts = new_output_counts(&config).unwrap();

        increment_output_counts_from_normalized_chunk(&config, &mut counts, &[pair]);

        let OutputCounts::Sketch(sketch) = counts else {
            panic!("cells= should use a bounded output sketch for kept-kmer side counts");
        };
        assert_eq!(sketch.words.len(), 1);
        assert!(sketch.depth(&probe) > 0);
    }

    #[test]
    fn nondeterministic_atomic_output_counts_direct_path_matches_sequential_without_collisions() {
        let config = Config {
            k: 5,
            min_quality: 0,
            min_prob: 0.0,
            deterministic: false,
            count_min: crate::cli::CountMinSettings {
                cells: Some(8192),
                hashes: Some(1),
                bits: Some(32),
                memory_bytes: None,
            },
            ..Config::default()
        };
        let kept_a = record("r1", b"ACGTACGTAC");
        let kept_b = record("r2", b"TTTTCCCCAA");
        let tossed = record("r3", b"GGGGAAAACC");
        let pairs = vec![
            NormalizedPair {
                input_list_index: 0,
                r1: kept_a.clone(),
                r2: None,
                out_r1: kept_a,
                out_r2: None,
                decision: PairDecision::default(),
                uncorrectable: false,
                read_count: 1,
                base_count: 10,
            },
            NormalizedPair {
                input_list_index: 0,
                r1: kept_b.clone(),
                r2: None,
                out_r1: kept_b,
                out_r2: None,
                decision: PairDecision::default(),
                uncorrectable: false,
                read_count: 1,
                base_count: 10,
            },
            NormalizedPair {
                input_list_index: 0,
                r1: tossed.clone(),
                r2: None,
                out_r1: tossed,
                out_r2: None,
                decision: PairDecision {
                    toss: true,
                    ..PairDecision::default()
                },
                uncorrectable: false,
                read_count: 1,
                base_count: 10,
            },
        ];
        let sequential_config = Config {
            deterministic: true,
            ..config.clone()
        };
        let mut sequential = new_output_counts(&sequential_config).unwrap();
        let mut direct = new_output_counts(&config).unwrap();

        increment_output_counts_from_normalized_chunk(&sequential_config, &mut sequential, &pairs);
        increment_output_counts_from_normalized_chunk(&config, &mut direct, &pairs);

        let (OutputCounts::AtomicSketch(sequential), OutputCounts::AtomicSketch(direct)) =
            (sequential, direct)
        else {
            panic!("bits=32 output counts should use atomic sketches");
        };
        assert_eq!(
            direct.increments.load(Ordering::Relaxed),
            sequential.increments.load(Ordering::Relaxed)
        );
        assert_eq!(
            direct.occupied_slots.load(Ordering::Relaxed),
            sequential.occupied_slots.load(Ordering::Relaxed)
        );
        for slot in 0..sequential.cells {
            assert_eq!(
                u64::from(direct.cells_by_hash[slot].load(Ordering::Relaxed)),
                u64::from(sequential.cells_by_hash[slot].load(Ordering::Relaxed))
            );
        }
    }

    #[test]
    fn bounded_sketch_memory_budget_derives_cell_count() {
        let config = Config {
            count_min: crate::cli::CountMinSettings {
                cells: None,
                hashes: Some(2),
                bits: Some(8),
                memory_bytes: Some(1000),
            },
            threads: Some(2),
            ..Config::default()
        };

        let sketch = new_bounded_count_min_sketch(&config).unwrap();

        assert_eq!(sketch.cells, 998);
        assert_eq!(sketch.words.len(), 125);
    }

    #[test]
    fn count_min_table_sizing_prime_adjusts_like_kcountarray() {
        assert_eq!(count_min_table_cells_from_total(1, 3), 1);
        assert_eq!(count_min_table_cells_from_total(9, 3), 7);
        assert_eq!(count_min_table_cells_from_total(64, 3), 62);
        assert_eq!(count_min_table_cells_from_total(1000, 2), 998);
    }

    #[test]
    fn non_prefiltered_short_kmer_sketch_caps_cells_to_kmer_space_like_bbnorm() {
        let config = Config {
            k: 3,
            count_min: crate::cli::CountMinSettings {
                cells: Some(10_000),
                hashes: Some(2),
                bits: Some(8),
                memory_bytes: None,
            },
            ..Config::default()
        };

        assert_eq!(short_kmer_space_cells(3), Some(64));
        assert_eq!(main_count_min_total_cells(&config, 8), 64);

        let sketch = new_bounded_count_min_sketch(&config).unwrap();
        assert!(sketch.cells <= 64);
    }

    #[test]
    fn prefiltered_short_kmer_sketch_preserves_requested_cells_like_bbnorm() {
        let config = Config {
            k: 3,
            count_min: crate::cli::CountMinSettings {
                cells: Some(10_000),
                hashes: Some(2),
                bits: Some(8),
                memory_bytes: None,
            },
            prefilter: crate::cli::PrefilterSettings {
                cells: Some(128),
                hashes: Some(2),
                bits: Some(2),
                ..Default::default()
            },
            ..Config::default()
        };

        assert_eq!(main_count_min_total_cells(&config, 8), 10_000);
    }

    #[test]
    fn kcount_array_min_arrays_rounds_threads_like_bbtools() {
        assert_eq!(kcount_array_min_arrays_for_threads(1), 2);
        assert_eq!(kcount_array_min_arrays_for_threads(2), 2);
        assert_eq!(kcount_array_min_arrays_for_threads(3), 4);
        assert_eq!(kcount_array_min_arrays_for_threads(8), 8);
        assert_eq!(kcount_array_min_arrays_for_threads(9), 16);
    }

    #[test]
    fn bounded_sketch_sizing_uses_configured_threads_for_kcount_arrays() {
        let config = Config {
            threads: Some(8),
            count_min: crate::cli::CountMinSettings {
                cells: Some(1000),
                hashes: Some(2),
                bits: Some(8),
                memory_bytes: None,
            },
            ..Config::default()
        };

        let sketch = new_bounded_count_min_sketch(&config).unwrap();

        assert_eq!(sketch.cells, 904);
        assert_eq!(sketch.words.len(), 113);
        assert_eq!(sketch.layout.array_mask, 7);
        assert_eq!(sketch.layout.array_bits, 3);
        assert_eq!(sketch.layout.cells_per_array, 113);
    }

    #[test]
    fn bounded_sketch_sizing_uses_active_rayon_threads_for_auto_threads() {
        let pool = rayon::ThreadPoolBuilder::new()
            .num_threads(3)
            .build()
            .unwrap();
        pool.install(|| {
            let config = Config {
                threads: None,
                count_min: crate::cli::CountMinSettings {
                    cells: Some(1000),
                    hashes: Some(2),
                    bits: Some(8),
                    memory_bytes: None,
                },
                ..Config::default()
            };

            let sketch = new_bounded_count_min_sketch(&config).unwrap();

            assert_eq!(kcount_array_min_arrays(&config), 4);
            assert_eq!(sketch.cells, 964);
            assert_eq!(sketch.words.len(), 121);
            assert_eq!(sketch.layout.array_mask, 3);
            assert_eq!(sketch.layout.array_bits, 2);
            assert_eq!(sketch.layout.cells_per_array, 241);
        });
    }

    #[test]
    fn explicit_count_min_cells_are_total_budget_like_bbtools() {
        let config = Config {
            count_min: crate::cli::CountMinSettings {
                cells: Some(9),
                hashes: Some(3),
                bits: Some(8),
                memory_bytes: None,
            },
            ..Config::default()
        };

        let packed = new_bounded_count_min_sketch(&config).unwrap();
        let atomic = new_atomic_count_min_sketch(&Config {
            count_min: crate::cli::CountMinSettings {
                bits: Some(32),
                ..config.count_min
            },
            ..Config::default()
        })
        .unwrap();

        assert_eq!(packed.cells, 7);
        assert_eq!(packed.words.len(), 1);
        assert_eq!(atomic.cells, 7);
        assert_eq!(atomic.cells_by_hash.len(), 7);
    }

    #[test]
    fn automatic_memory_budget_uses_bbtools_sizing_formula() {
        let config = Config {
            hist_in: Some(PathBuf::from("hist.tsv")),
            hist_len: 1000,
            threads: Some(3),
            build_passes: 2,
            ..Config::default()
        };

        let usable = bbtools_usable_table_memory_bytes(&config, 1_000_000_000);

        assert_eq!(usable, 329_944_000);
    }

    #[test]
    fn countup_auto_memory_budget_halves_filter_bytes_like_bbnorm() {
        let config = Config {
            auto_count_min_memory_bytes: Some(1_000_000_000),
            table_reads: Some(1_000_000),
            ..Config::default()
        };
        let countup_config = Config {
            count_up: true,
            ..config.clone()
        };

        assert_eq!(automatic_count_min_memory_bytes(&config), Some(659_920_000));
        assert_eq!(
            automatic_count_min_memory_bytes(&countup_config),
            Some(329_960_000)
        );
    }

    #[test]
    fn automatic_output_counts_use_side_budget_and_next_mask_seed() {
        let config = Config {
            auto_count_min_memory_bytes: Some(1_000_000_000),
            table_reads: Some(1_000_000),
            threads: Some(8),
            deterministic: false,
            ..Config::default()
        };

        assert_eq!(automatic_count_min_memory_bytes(&config), Some(659_920_000));
        assert_eq!(
            output_count_min_memory_bytes(&config, 32),
            Some(164_980_000)
        );

        let main = new_atomic_count_min_sketch(&config).unwrap();
        let output = new_output_counts(&config).unwrap();
        let OutputCounts::AtomicSketch(output) = output else {
            panic!("automatic bits=32 output counts should use atomic sketches");
        };
        let main_layout = main.layout_summary("input_main", None);
        let output_layout = output.layout_summary("output_kept", None);

        assert_eq!(
            output_layout.mask_seed,
            BBTOOLS_KCOUNT_ARRAY_SECOND_MASK_SEED
        );
        assert!(output_layout.memory_bytes < main_layout.memory_bytes / 2);
        assert!(output_layout.memory_bytes >= OUTPUT_COUNT_MIN_AUTO_MIN_MEMORY_BYTES);
    }

    #[test]
    fn explicit_output_count_memory_preserves_requested_budget() {
        let config = Config {
            count_min: crate::cli::CountMinSettings {
                cells: None,
                hashes: Some(3),
                bits: Some(32),
                memory_bytes: Some(128 * 1024 * 1024),
            },
            threads: Some(4),
            ..Config::default()
        };

        assert_eq!(
            output_count_min_memory_bytes(&config, 32),
            Some(128 * 1024 * 1024)
        );
        let main = new_atomic_count_min_sketch(&config).unwrap();
        let output = new_output_counts(&config).unwrap();
        let OutputCounts::AtomicSketch(output) = output else {
            panic!("explicit bits=32 output counts should use atomic sketches");
        };

        assert_eq!(output.cells, main.cells);
        assert_eq!(
            output.layout.mask_seed,
            BBTOOLS_KCOUNT_ARRAY_SECOND_MASK_SEED
        );
    }

    #[test]
    fn constrained_prefilter_inflates_unsaturated_colliding_counts() {
        let config = Config {
            prefilter: crate::cli::PrefilterSettings {
                enabled: false,
                force_disabled: false,
                cells: Some(1),
                hashes: Some(2),
                bits: Some(8),
                memory_bytes: None,
                memory_fraction_micros: None,
            },
            ..Config::default()
        };
        let mut counts = CountMap::default();
        counts.insert(KmerKey::Short(7), 2);
        counts.insert(KmerKey::Short(11), 5);

        apply_prefilter_collision_estimates(&config, &mut counts);

        assert_eq!(counts.get(&KmerKey::Short(7)), Some(&7));
        assert_eq!(counts.get(&KmerKey::Short(11)), Some(&7));
    }

    #[test]
    fn constrained_prefilter_keeps_exact_counts_after_saturation() {
        let config = Config {
            prefilter: crate::cli::PrefilterSettings {
                enabled: false,
                force_disabled: false,
                cells: Some(1),
                hashes: Some(1),
                bits: Some(2),
                memory_bytes: None,
                memory_fraction_micros: None,
            },
            ..Config::default()
        };
        let mut counts = CountMap::default();
        counts.insert(KmerKey::Short(7), 2);
        counts.insert(KmerKey::Short(11), 5);

        apply_prefilter_collision_estimates(&config, &mut counts);

        assert_eq!(counts.get(&KmerKey::Short(7)), Some(&2));
        assert_eq!(counts.get(&KmerKey::Short(11)), Some(&5));
    }

    #[test]
    fn prefilter_memory_budget_derives_prime_table_cells() {
        let config = Config {
            prefilter: crate::cli::PrefilterSettings {
                enabled: false,
                force_disabled: false,
                cells: None,
                hashes: Some(2),
                bits: Some(8),
                memory_bytes: Some(1000),
                memory_fraction_micros: None,
            },
            ..Config::default()
        };
        let mut counts = CountMap::default();
        counts.insert(KmerKey::Short(7), 2);
        counts.insert(KmerKey::Short(11), 5);

        let bits = config.prefilter.bits.unwrap();
        let total_cells = count_min_cells_from_memory(config.prefilter.memory_bytes, bits);
        let table_cells = count_min_table_cells_from_total_bits(total_cells, bits);

        assert_eq!(total_cells, 1000);
        assert_eq!(table_cells, 998);

        apply_prefilter_collision_estimates(&config, &mut counts);

        assert_eq!(counts.get(&KmerKey::Short(7)), Some(&2));
        assert_eq!(counts.get(&KmerKey::Short(11)), Some(&5));
    }

    #[test]
    fn prefilter_fraction_derives_memory_from_table_budget() {
        let config = Config {
            auto_count_min_memory_bytes: Some(10_000),
            prefilter: crate::cli::PrefilterSettings {
                enabled: true,
                force_disabled: false,
                cells: None,
                hashes: Some(2),
                bits: Some(8),
                memory_bytes: None,
                memory_fraction_micros: Some(350_000),
            },
            ..Config::default()
        };
        let mut counts = CountMap::default();
        counts.insert(KmerKey::Short(7), 2);
        counts.insert(KmerKey::Short(11), 5);

        let total_cells = prefilter_total_cells(&config, config.prefilter.bits.unwrap());
        let table_cells =
            count_min_table_cells_from_total_bits(total_cells, config.prefilter.bits.unwrap());

        assert_eq!(total_cells, 3500);
        assert_eq!(table_cells, 3494);

        apply_prefilter_collision_estimates(&config, &mut counts);

        assert_eq!(counts.get(&KmerKey::Short(7)), Some(&2));
        assert_eq!(counts.get(&KmerKey::Short(11)), Some(&5));
    }

    #[test]
    fn prefilter_fraction_partitions_main_cell_budget() {
        let config = Config {
            count_min: crate::cli::CountMinSettings {
                cells: Some(1000),
                hashes: Some(1),
                bits: Some(32),
                memory_bytes: None,
            },
            prefilter: crate::cli::PrefilterSettings {
                enabled: true,
                force_disabled: false,
                cells: None,
                hashes: Some(1),
                bits: Some(2),
                memory_bytes: None,
                memory_fraction_micros: Some(250_000),
            },
            threads: Some(2),
            ..Config::default()
        };

        assert_eq!(main_count_min_total_cells(&config, 32), 750);
        assert_eq!(prefilter_total_cells(&config, 2), 4000);

        let main = new_atomic_count_min_sketch(&config).unwrap();
        let prefilter = new_prefilter_count_min_sketch(&config).unwrap();

        assert_eq!(main.cells, count_min_table_cells_from_total_bits(750, 32));
        assert_eq!(
            prefilter.cells,
            count_min_table_cells_from_total_bits(4000, 2)
        );
        assert_eq!(prefilter.max_count, 3);
    }

    #[test]
    fn prefilter_flag_uses_bbtools_default_fraction_on_bounded_count_min_paths() {
        let config = Config {
            count_min: crate::cli::CountMinSettings {
                cells: Some(10_000),
                hashes: Some(2),
                bits: Some(32),
                memory_bytes: None,
            },
            prefilter: crate::cli::PrefilterSettings {
                enabled: true,
                force_disabled: false,
                cells: None,
                hashes: Some(2),
                bits: Some(2),
                memory_bytes: None,
                memory_fraction_micros: None,
            },
            ..Config::default()
        };

        assert!(use_prefilter_collision_estimates(&config));
        assert_eq!(main_count_min_total_cells(&config, 32), 6500);
        assert_eq!(prefilter_total_cells(&config, 2), 56_000);
    }

    #[test]
    fn zero_prefilter_fraction_does_not_force_prefilter_sketch() {
        let config = Config {
            count_min: crate::cli::CountMinSettings {
                cells: Some(10_000),
                bits: Some(32),
                ..Default::default()
            },
            prefilter: crate::cli::PrefilterSettings {
                enabled: false,
                force_disabled: false,
                memory_fraction_micros: Some(0),
                ..Default::default()
            },
            ..Config::default()
        };

        assert!(!use_prefilter_collision_estimates(&config));
        assert_eq!(main_count_min_total_cells(&config, 32), 10_000);
    }

    #[test]
    fn forced_off_prefilter_ignores_lingering_controls_like_bbnorm() {
        let config = Config {
            count_min: crate::cli::CountMinSettings {
                cells: Some(10_000),
                hashes: Some(3),
                bits: Some(32),
                ..Default::default()
            },
            prefilter: crate::cli::PrefilterSettings {
                enabled: false,
                force_disabled: true,
                cells: Some(1_000),
                hashes: Some(1),
                bits: Some(2),
                memory_bytes: None,
                memory_fraction_micros: Some(DEFAULT_PREFILTER_FRACTION_MICROS),
            },
            ..Config::default()
        };

        assert!(!use_prefilter_collision_estimates(&config));
        assert_eq!(prefilter_memory_fraction_micros(&config), None);
        assert_eq!(main_count_min_total_cells(&config, 32), 10_000);
    }

    #[test]
    fn prefilter_default_hashes_track_main_hashes_like_bbnorm() {
        let config = Config {
            count_min: crate::cli::CountMinSettings {
                cells: Some(10_000),
                hashes: Some(8),
                bits: Some(32),
                memory_bytes: None,
            },
            prefilter: crate::cli::PrefilterSettings {
                enabled: true,
                force_disabled: false,
                bits: Some(2),
                ..Default::default()
            },
            ..Config::default()
        };

        let prefilter = new_prefilter_count_min_sketch(&config).unwrap();

        assert_eq!(default_prefilter_hashes(&config), 4);
        assert_eq!(prefilter.hashes, 4);

        let explicit = Config {
            prefilter: crate::cli::PrefilterSettings {
                hashes: Some(1),
                ..config.prefilter
            },
            ..config
        };
        let prefilter = new_prefilter_count_min_sketch(&explicit).unwrap();
        assert_eq!(prefilter.hashes, 1);
    }

    #[test]
    fn explicit_prefilter_hashes_enable_default_partition_like_bbnorm() {
        let config = Config {
            count_min: crate::cli::CountMinSettings {
                cells: Some(10_000),
                hashes: Some(3),
                bits: Some(32),
                memory_bytes: None,
            },
            prefilter: crate::cli::PrefilterSettings {
                enabled: true,
                force_disabled: false,
                hashes: Some(1),
                bits: Some(2),
                ..Default::default()
            },
            ..Config::default()
        };

        assert_eq!(
            prefilter_memory_fraction_micros(&config),
            Some(DEFAULT_PREFILTER_FRACTION_MICROS)
        );
        assert_eq!(main_count_min_total_cells(&config, 32), 6500);
        assert_eq!(prefilter_total_cells(&config, 2), 56_000);
    }

    #[test]
    fn prefilter_flag_alone_keeps_small_exact_inputs_on_exact_path() {
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("reads.fq");
        write_fastq(&path, &[("r1", b"ACGTACGT", b"IIIIIIII")]);
        let config = Config {
            in1: Some(path),
            k: 3,
            min_quality: 0,
            min_prob: 0.0,
            prefilter: crate::cli::PrefilterSettings {
                enabled: true,
                force_disabled: false,
                ..Default::default()
            },
            ..Config::default()
        };

        let counts = build_input_counts(&config).unwrap();
        assert!(matches!(counts, InputCounts::Exact(_)));
    }

    #[test]
    fn prefilter_flag_builds_two_stage_sketch_when_count_min_is_bounded() {
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("reads.fq");
        write_fastq(
            &path,
            &[
                ("r1", b"ACGTACGT", b"IIIIIIII"),
                ("r2", b"ACGTACGT", b"IIIIIIII"),
                ("r3", b"ACGTACGT", b"IIIIIIII"),
            ],
        );
        let config = Config {
            in1: Some(path),
            k: 3,
            min_quality: 0,
            min_prob: 0.0,
            count_min: crate::cli::CountMinSettings {
                cells: Some(512),
                hashes: Some(2),
                bits: Some(32),
                memory_bytes: None,
            },
            prefilter: crate::cli::PrefilterSettings {
                enabled: true,
                force_disabled: false,
                ..Default::default()
            },
            ..Config::default()
        };

        let counts = build_input_counts(&config).unwrap();
        let InputCounts::PrefilteredSketch {
            prefilter,
            limit,
            main,
        } = counts
        else {
            panic!("prefilter=t plus bounded count-min should build a two-stage sketch");
        };
        assert_eq!(prefilter.bits(), DEFAULT_PREFILTER_BITS);
        assert_eq!(limit, prefilter.max_count());
        assert_eq!(prefilter_total_cells(&config, DEFAULT_PREFILTER_BITS), 2867);
        assert_eq!(main_count_min_total_cells(&config, 32), 332);
        assert!(matches!(*main, InputCounts::AtomicSketch(_)));
    }

    #[test]
    fn explicit_prefilter_memory_does_not_shrink_main_table_budget() {
        let config = Config {
            count_min: crate::cli::CountMinSettings {
                cells: Some(1000),
                hashes: Some(1),
                bits: Some(32),
                memory_bytes: None,
            },
            prefilter: crate::cli::PrefilterSettings {
                enabled: true,
                force_disabled: false,
                cells: None,
                hashes: Some(1),
                bits: Some(2),
                memory_bytes: Some(100),
                memory_fraction_micros: Some(250_000),
            },
            ..Config::default()
        };

        assert_eq!(main_count_min_total_cells(&config, 32), 1000);
        assert_eq!(prefilter_total_cells(&config, 2), 400);
    }

    #[test]
    fn prefiltered_input_counts_use_prefilter_until_saturation() {
        let low = KmerKey::Short(1);
        let high = KmerKey::Short(2);
        let mut prefilter = PackedCountMinSketch::new(4099, 2, 2).unwrap();
        prefilter.add_key_count(&low, 2);
        prefilter.add_key_count(&high, 3);

        let main = AtomicCountMinSketch::new(4099, 2).unwrap();
        main.add_key_count(&low, 99);
        main.add_key_count(&high, 5);

        let counts = InputCounts::PrefilteredSketch {
            limit: prefilter.max_count,
            prefilter: PrefilterCountMinSketch::Packed(prefilter),
            main: Box::new(InputCounts::AtomicSketch(main)),
        };

        assert_eq!(counts.depth(&low), 2);
        assert_eq!(counts.depth(&high), 5);
    }

    #[test]
    fn prefiltered_input_counts_honor_explicit_lower_prefilter_limit() {
        let key = KmerKey::Short(7);
        let mut prefilter = PackedCountMinSketch::new(4099, 2, 2).unwrap();
        prefilter.add_key_count(&key, 2);

        let main = AtomicCountMinSketch::new(4099, 2).unwrap();
        main.add_key_count(&key, 11);

        let counts = InputCounts::PrefilteredSketch {
            limit: 2,
            prefilter: PrefilterCountMinSketch::Packed(prefilter),
            main: Box::new(InputCounts::AtomicSketch(main)),
        };

        assert_eq!(counts.depth(&key), 11);
    }

    #[test]
    fn input_count_layout_summary_reports_prefilter_and_main_tables() {
        let prefilter =
            PackedCountMinSketch::new_with_min_arrays_and_mask_seed(4099, 2, 2, 4, 0).unwrap();
        let main = AtomicCountMinSketch::new_with_min_arrays_and_update_mode(
            8191,
            3,
            4,
            CountMinUpdateMode::Conservative,
            7,
        )
        .unwrap();
        let counts = InputCounts::PrefilteredSketch {
            limit: prefilter.max_count,
            prefilter: PrefilterCountMinSketch::Packed(prefilter),
            main: Box::new(InputCounts::AtomicSketch(main)),
        };

        let layouts = counts.sketch_layouts();

        assert_eq!(layouts.len(), 2);
        assert_eq!(layouts[0].table, "input_prefilter");
        assert_eq!(layouts[0].kind, "packed");
        assert_eq!(layouts[0].bits, 2);
        assert_eq!(layouts[0].hashes, 2);
        assert_eq!(layouts[0].mask_seed, 0);
        assert_eq!(layouts[0].update_mode, "conservative");
        assert_eq!(layouts[0].prefilter_limit, Some(3));
        assert!(layouts[0].memory_bytes > 0);
        assert_eq!(layouts[1].table, "input_main");
        assert_eq!(layouts[1].kind, "atomic");
        assert_eq!(layouts[1].bits, 32);
        assert_eq!(layouts[1].hashes, 3);
        assert_eq!(layouts[1].mask_seed, 7);
        assert_eq!(layouts[1].prefilter_limit, None);
        assert!(layouts[1].arrays >= 4);
        assert!(layouts[1].memory_bytes >= layouts[1].cells * std::mem::size_of::<AtomicU32>());
    }

    #[test]
    fn prefilter_gate_uses_explicit_limit_for_main_counts() {
        let below = KmerKey::Short(1);
        let at_limit = KmerKey::Short(2);
        let above = KmerKey::Short(3);
        let mut prefilter = PackedCountMinSketch::new(4099, 2, 2).unwrap();
        prefilter.add_key_count(&below, 1);
        prefilter.add_key_count(&at_limit, 2);
        prefilter.add_key_count(&above, 3);

        let mut counts = CountMap::default();
        counts.insert(below.clone(), 10);
        counts.insert(at_limit.clone(), 20);
        counts.insert(above.clone(), 30);

        let prefilter = PrefilterCountMinSketch::Packed(prefilter);
        retain_prefilter_saturated_counts(&mut counts, Some(PrefilterGate::new(&prefilter, 2)));

        assert!(!counts.contains_key(&below));
        assert_eq!(counts.get(&at_limit), Some(&20));
        assert_eq!(counts.get(&above), Some(&30));
    }

    #[test]
    fn prefilter_gate_during_collection_matches_post_retain() {
        let r1 = record("r1", b"ACGTACGTACGT");
        let r2 = record("r2", b"TGCATGCATGCA");

        for remove_duplicate_kmers in [false, true] {
            let config = Config {
                k: 3,
                min_quality: 0,
                min_prob: 0.0,
                remove_duplicate_kmers,
                ..Config::default()
            };
            let mut prefilter = PackedCountMinSketch::new(4099, 2, 2).unwrap();
            let keys = unique_pair_kmers(&config, &r1, Some(&r2));
            for key in keys.iter().step_by(2) {
                prefilter.add_key_count(key, prefilter.max_count);
            }
            let prefilter = PrefilterCountMinSketch::Packed(prefilter);
            let gate = PrefilterGate::new(&prefilter, prefilter.max_count());
            assert!(
                keys.iter().any(|key| !gate.should_count_in_main(key)),
                "fixture should include at least one prefilter-rejected k-mer"
            );

            let mut post_retain = CountMap::default();
            increment_pair_counts(&config, &mut post_retain, &r1, Some(&r2));
            retain_prefilter_saturated_counts(&mut post_retain, Some(gate));

            let mut during_collection = CountMap::default();
            increment_pair_counts_with_prefilter(
                &config,
                &mut during_collection,
                &r1,
                Some(&r2),
                Some(gate),
            );

            assert_eq!(during_collection, post_retain);
        }
    }

    #[test]
    fn prefiltered_input_counts_use_thresholded_main_unique_estimates_above_prefilter_max() {
        let mut prefilter = PackedCountMinSketch::new(1024, 4, 2).unwrap();
        for bucket in 0..256 {
            let depth = if bucket < 128 { prefilter.max_count } else { 1 };
            prefilter.set_cell(bucket, depth);
        }
        prefilter.increments = 1_000;

        let main = AtomicCountMinSketch::new(1024, 4).unwrap();
        for bucket in 0..128 {
            main.cells_by_hash[bucket].store(4, Ordering::Relaxed);
        }
        main.add_key_increments(1_000);

        let counts = InputCounts::PrefilteredSketch {
            limit: prefilter.max_count,
            prefilter: PrefilterCountMinSketch::Packed(prefilter),
            main: Box::new(InputCounts::AtomicSketch(main)),
        };

        let all_depth_estimated = counts.unique_kmers();
        let saturated_prefilter_estimated = counts.unique_kmers_at_least(2);
        let high_depth_estimated = counts.unique_kmers_at_least(4);
        let split = counts.unique_kmer_estimate_split().unwrap();

        assert!(
            (70..=80).contains(&all_depth_estimated),
            "prefilter all-depth estimate was {all_depth_estimated}"
        );
        assert!(
            (30..=40).contains(&saturated_prefilter_estimated),
            "prefilter threshold estimate was {saturated_prefilter_estimated}"
        );
        assert!(
            (30..=40).contains(&high_depth_estimated),
            "main high-depth estimate was {high_depth_estimated}"
        );
        assert_eq!(split.low_depth_max, 3);
        assert_eq!(split.high_depth_min, 4);
        assert_eq!(split.high_depth_kmers, high_depth_estimated);
        assert_eq!(
            split.low_depth_kmers,
            all_depth_estimated.saturating_sub(high_depth_estimated)
        );
        assert!(
            (30..=50).contains(&split.low_depth_kmers),
            "prefilter low-depth split estimate was {}",
            split.low_depth_kmers
        );
    }

    #[test]
    fn bounded_input_counts_builds_two_stage_prefiltered_sketch() {
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("reads.fq");
        write_fastq(
            &path,
            &[
                ("r1", b"ACGTACGT", b"IIIIIIII"),
                ("r2", b"ACGTACGT", b"IIIIIIII"),
                ("r3", b"ACGTACGT", b"IIIIIIII"),
            ],
        );
        let config = Config {
            in1: Some(path),
            k: 3,
            min_quality: 0,
            min_prob: 0.0,
            count_min: crate::cli::CountMinSettings {
                cells: Some(128),
                hashes: Some(2),
                bits: Some(32),
                memory_bytes: None,
            },
            prefilter: crate::cli::PrefilterSettings {
                enabled: false,
                force_disabled: false,
                cells: None,
                hashes: Some(2),
                bits: None,
                memory_bytes: Some(1024),
                memory_fraction_micros: None,
            },
            ..Config::default()
        };

        let counts = build_input_counts(&config).unwrap();
        let InputCounts::PrefilteredSketch {
            prefilter,
            limit,
            main,
        } = counts
        else {
            panic!("prefilter memory plus bounded count-min should build a two-stage sketch");
        };
        assert_eq!(prefilter.bits(), DEFAULT_PREFILTER_BITS);
        assert_eq!(prefilter.max_count(), 3);
        assert_eq!(limit, prefilter.max_count());
        assert_eq!(prefilter.update_mode(), CountMinUpdateMode::Conservative);
        assert!(matches!(*main, InputCounts::AtomicSketch(_)));
    }

    #[test]
    fn trusted_build_pass_filter_reduces_non_singleton_depths() {
        let config = Config {
            build_passes: 2,
            ..Config::default()
        };
        let mut counts = CountMap::default();
        counts.insert(KmerKey::Short(7), 1);
        counts.insert(KmerKey::Short(11), 2);
        counts.insert(KmerKey::Short(13), 3);

        apply_trusted_build_pass_filter(&config, &mut counts);

        assert_eq!(counts.get(&KmerKey::Short(7)), Some(&1));
        assert_eq!(counts.get(&KmerKey::Short(11)), Some(&1));
        assert_eq!(counts.get(&KmerKey::Short(13)), Some(&2));
    }

    #[test]
    fn ecco_auto_disables_overlap_repair_when_java_style_sample_is_empty() {
        let dir = tempfile::tempdir().unwrap();
        let r1_path = dir.path().join("r1.fq");
        let r2_path = dir.path().join("r2.fq");
        let r1 = b"TTAGTTGTGCCGCAGCGAAGTAGTGCTTGAAATATGCGAC";
        let r2 = b"GTCGCATATTTCAAGCACTAATTCGCTGCGGCACAACTAA";
        let q = b"IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII";
        write_fastq(
            &r1_path,
            &[
                ("overlap1/1", r1, q),
                ("overlap2/1", r1, q),
                ("overlap3/1", r1, q),
                ("overlap4/1", r1, q),
            ],
        );
        write_fastq(
            &r2_path,
            &[
                ("overlap1/2", r2, q),
                ("overlap2/2", r2, q),
                ("overlap3/2", r2, q),
                ("overlap4/2", r2, q),
            ],
        );
        let config = Config {
            in1: Some(r1_path),
            in2: Some(r2_path),
            error_correct: true,
            error_correct_first: true,
            error_correct_final: true,
            overlap_error_correct_auto: true,
            ..Config::default()
        };

        let resolved = resolve_overlap_error_correct_auto(&config).unwrap();

        assert!(!resolved.overlap_error_correct_auto);
        assert!(!resolved.overlap_error_correct);
    }

    #[test]
    fn ecco_auto_enables_overlap_repair_for_sampled_mergeable_pairs() {
        let dir = tempfile::tempdir().unwrap();
        let r1_path = dir.path().join("r1.fq");
        let r2_path = dir.path().join("r2.fq");
        let r1 = b"TTAGTTGTGCCGCAGCGAAGTAGTGCTTGAAATATGCGAC";
        let r2 = b"GTCGCATATTTCAAGCACTACTTCGCTGCGGCACAACTAA";
        let q = b"IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII";
        write_repeated_fastq(&r1_path, "overlap/1_", r1, q, 200);
        write_repeated_fastq(&r2_path, "overlap/2_", r2, q, 200);
        let config = Config {
            in1: Some(r1_path),
            in2: Some(r2_path),
            error_correct: true,
            error_correct_first: true,
            error_correct_final: true,
            overlap_error_correct_auto: true,
            ..Config::default()
        };

        let resolved = resolve_overlap_error_correct_auto(&config).unwrap();

        assert!(!resolved.overlap_error_correct_auto);
        assert!(resolved.overlap_error_correct);
    }

    #[test]
    fn countup_abrc_controls_tossed_read_table_updates() {
        let keys = vec![KmerKey::Short(7), KmerKey::Short(11)];
        let mut input_counts = CountMap::default();
        input_counts.insert(keys[0].clone(), 3);
        input_counts.insert(keys[1].clone(), 3);

        let base_config = Config {
            min_depth: 1,
            ..Config::default()
        };
        let mut kept_counts = OutputCounts::Exact(CountMap::default());
        update_countup_kept_counts_for_decision(
            &base_config,
            &mut kept_counts,
            &input_counts,
            &keys,
            true,
        );
        assert_eq!(kept_counts.unique_kmers(), 0);

        let add_bad_config = Config {
            add_bad_reads_countup: true,
            ..base_config.clone()
        };
        update_countup_kept_counts_for_decision(
            &add_bad_config,
            &mut kept_counts,
            &input_counts,
            &keys,
            true,
        );
        assert_eq!(kept_counts.depth(&keys[0]), 1);
        assert_eq!(kept_counts.depth(&keys[1]), 1);

        update_countup_kept_counts_for_decision(
            &base_config,
            &mut kept_counts,
            &input_counts,
            &keys,
            false,
        );
        assert_eq!(kept_counts.depth(&keys[0]), 2);
        assert_eq!(kept_counts.depth(&keys[1]), 2);
    }

    #[test]
    fn countup_decision_plan_reuses_input_depth_gate_for_kept_updates() {
        let keys = vec![KmerKey::Short(7), KmerKey::Short(11), KmerKey::Short(13)];
        let mut input_counts = CountMap::default();
        input_counts.insert(keys[0].clone(), 0);
        input_counts.insert(keys[1].clone(), 3);
        input_counts.insert(keys[2].clone(), 4);
        let kept_counts = CountMap::default();
        let config = Config {
            min_depth: 2,
            min_kmers_over_min_depth: 1,
            target_depth: 10,
            add_bad_reads_countup: true,
            ..Config::default()
        };

        let plan = countup_decision_plan(&config, &input_counts, &kept_counts, &keys, 10);
        assert_eq!(
            plan.toss,
            decide_countup_pair(&config, &input_counts, &kept_counts, &keys, 10)
        );
        assert_eq!(plan.eligible_key_indices, vec![1, 2]);

        let mut planned_counts = OutputCounts::Exact(CountMap::default());
        update_countup_kept_counts_for_plan(&config, &mut planned_counts, &keys, &plan);

        let mut replayed_counts = OutputCounts::Exact(CountMap::default());
        update_countup_kept_counts_for_decision(
            &config,
            &mut replayed_counts,
            &input_counts,
            &keys,
            plan.toss,
        );

        assert_eq!(
            planned_counts.unique_kmers(),
            replayed_counts.unique_kmers()
        );
        assert_eq!(planned_counts.depth(&keys[0]), 0);
        assert_eq!(planned_counts.depth(&keys[1]), 1);
        assert_eq!(planned_counts.depth(&keys[2]), 1);
    }

    #[test]
    fn countup_bounded_kept_counts_use_sketch_when_cells_are_constrained() {
        let keys = vec![KmerKey::Short(7), KmerKey::Short(11)];
        let mut input_counts = CountMap::default();
        input_counts.insert(keys[0].clone(), 3);
        input_counts.insert(keys[1].clone(), 3);
        let config = Config {
            min_depth: 1,
            count_min: crate::cli::CountMinSettings {
                cells: Some(1),
                hashes: Some(2),
                bits: Some(3),
                memory_bytes: None,
            },
            ..Config::default()
        };
        let mut kept_counts = new_output_counts(&config).unwrap();

        update_countup_kept_counts_for_decision(
            &config,
            &mut kept_counts,
            &input_counts,
            &keys,
            false,
        );

        let OutputCounts::Sketch(sketch) = kept_counts else {
            panic!("countup cells= should use a bounded kept-count sketch");
        };
        assert_eq!(sketch.words.len(), 1);
        assert_eq!(sketch.depth(&keys[0]), 2);
        assert_eq!(sketch.depth(&keys[1]), 2);
    }

    #[test]
    fn countup_kept_count_sketch_uses_java_target_sized_cells() {
        let config = Config {
            count_up: true,
            target_depth: 100,
            threads: Some(1),
            count_min: crate::cli::CountMinSettings {
                cells: Some(10_000),
                hashes: Some(8),
                bits: Some(32),
                memory_bytes: None,
            },
            ..Config::default()
        };

        let kept_counts = new_output_counts(&config).unwrap();

        let OutputCounts::Sketch(sketch) = kept_counts else {
            panic!("countup kept-count table should use a packed sketch");
        };
        assert_eq!(sketch.bits, 8);
        assert_eq!(sketch.hashes, 3);
        assert_eq!(sketch.cells, 9_998);
        assert_eq!(
            sketch.layout.mask_seed,
            BBTOOLS_KCOUNT_ARRAY_SECOND_MASK_SEED
        );
    }

    #[test]
    fn countup_kept_count_bits_use_adjusted_target_boundaries_like_bbnorm() {
        assert_eq!(
            countup_output_count_bits(&Config {
                count_up: true,
                target_depth: 16,
                ..Config::default()
            }),
            4
        );
        assert_eq!(
            countup_output_count_bits(&Config {
                count_up: true,
                target_depth: 17,
                ..Config::default()
            }),
            8
        );
        assert_eq!(
            countup_output_count_bits(&Config {
                count_up: true,
                target_depth: 268,
                ..Config::default()
            }),
            8
        );
        assert_eq!(
            countup_output_count_bits(&Config {
                count_up: true,
                target_depth: 269,
                ..Config::default()
            }),
            16
        );
    }

    #[test]
    fn output_pair_analysis_is_only_required_for_rename_or_depth_bins() {
        assert!(!needs_output_pair_analysis(&Config::default()));
        assert!(needs_output_pair_analysis(&Config {
            rename_reads: true,
            ..Config::default()
        }));
        assert!(needs_output_pair_analysis(&Config {
            out_low1: Some(PathBuf::from("low.fq")),
            ..Config::default()
        }));
        assert!(needs_output_pair_analysis(&Config {
            out_high2: Some(PathBuf::from("high2.fq")),
            ..Config::default()
        }));
    }

    #[test]
    fn countup_kept_count_sketch_uses_next_mask_seed_after_prefilter_and_main() {
        let config = Config {
            count_up: true,
            target_depth: 100,
            threads: Some(1),
            prefilter: crate::cli::PrefilterSettings {
                enabled: true,
                force_disabled: false,
                ..Default::default()
            },
            count_min: crate::cli::CountMinSettings {
                cells: Some(10_000),
                hashes: Some(3),
                bits: Some(32),
                memory_bytes: None,
            },
            ..Config::default()
        };

        let kept_counts = new_output_counts(&config).unwrap();

        let OutputCounts::Sketch(sketch) = kept_counts else {
            panic!("countup kept-count table should use a packed sketch");
        };
        assert_eq!(
            sketch.layout.mask_seed,
            BBTOOLS_KCOUNT_ARRAY_SECOND_MASK_SEED + BBTOOLS_KCOUNT_ARRAY_MASK_SEED_STEP
        );
    }

    #[test]
    fn multipass_caps_wide_count_min_bits_like_bbnorm() {
        let mut default_bits = Config {
            passes: 2,
            ..Config::default()
        };
        apply_bbtools_multipass_cell_bits_cap(&mut default_bits);
        assert_eq!(default_bits.count_min.bits, Some(16));

        let mut explicit_wide_bits = Config {
            passes: 2,
            count_min: crate::cli::CountMinSettings {
                bits: Some(32),
                ..Default::default()
            },
            ..Config::default()
        };
        apply_bbtools_multipass_cell_bits_cap(&mut explicit_wide_bits);
        assert_eq!(explicit_wide_bits.count_min.bits, Some(16));

        let mut explicit_narrow_bits = Config {
            passes: 2,
            count_min: crate::cli::CountMinSettings {
                bits: Some(8),
                ..Default::default()
            },
            ..Config::default()
        };
        apply_bbtools_multipass_cell_bits_cap(&mut explicit_narrow_bits);
        assert_eq!(explicit_narrow_bits.count_min.bits, Some(8));

        let mut single_pass = Config {
            passes: 1,
            ..Config::default()
        };
        apply_bbtools_multipass_cell_bits_cap(&mut single_pass);
        assert_eq!(single_pass.count_min.bits, None);
    }

    #[test]
    fn multipass_intermediate_pass_uses_bits1_like_bbnorm() {
        let config = Config {
            passes: 2,
            count_min_bits_first: Some(8),
            count_min: crate::cli::CountMinSettings {
                bits: Some(16),
                ..Default::default()
            },
            ..Config::default()
        };

        let pass_config = pass_config_for_intermediate(
            &config,
            1,
            Path::new("in1.fq"),
            None,
            false,
            PathBuf::from("out1.fq"),
            None,
            None,
            None,
        );

        assert_eq!(pass_config.count_min.bits, Some(8));
        assert_eq!(config.count_min.bits, Some(16));
    }

    #[test]
    fn count_map_capacity_hint_uses_initialsize_and_prealloc() {
        let explicit = Config {
            table_initial_size: Some(1234),
            ..Config::default()
        };
        assert_eq!(count_map_capacity_hint(&explicit), Some(1234));

        let paired_prealloc = Config {
            table_prealloc_fraction: Some(0.5),
            table_reads: Some(10),
            in2: Some(PathBuf::from("mate.fq")),
            k: 31,
            ..Config::default()
        };
        assert_eq!(preallocation_capacity_hint(&paired_prealloc), Some(700));

        let larger_prealloc = Config {
            table_initial_size: Some(100),
            table_prealloc_fraction: Some(1.0),
            table_reads: Some(10),
            in2: Some(PathBuf::from("mate.fq")),
            k: 31,
            ..Config::default()
        };
        assert_eq!(count_map_capacity_hint(&larger_prealloc), Some(1400));
    }

    #[test]
    fn countup_presort_prefers_low_error_reads_like_java() {
        let config = Config {
            k: 3,
            min_depth: 1,
            low_thresh: 1,
            high_thresh: 3,
            error_detect_ratio: 2,
            low_percentile: 0.20,
            ..Config::default()
        };
        let clean = SequenceRecord {
            id: "clean".to_string(),
            numeric_id: 2,
            bases: b"AAAAAAAAAA".to_vec(),
            qualities: Some(vec![b'I'; 10]),
        };
        let noisy = SequenceRecord {
            id: "noisy".to_string(),
            numeric_id: 1,
            bases: b"AAAAACCCCC".to_vec(),
            qualities: Some(vec![b'I'; 10]),
        };
        let mut input_counts = CountMap::default();
        for key in kmers_for_record(&clean, &config) {
            input_counts.insert(key, 10);
        }

        let mut pairs = [
            CountupWorkPair {
                input_list_index: 0,
                sort_key: countup_sort_key(&config, &input_counts, &noisy, None, 0),
                r1: noisy,
                r2: None,
            },
            CountupWorkPair {
                input_list_index: 0,
                sort_key: countup_sort_key(&config, &input_counts, &clean, None, 1),
                r1: clean,
                r2: None,
            },
        ];
        pairs.sort_by(compare_countup_work_pairs);

        assert_eq!(pairs[0].r1.id, "clean");
        assert_eq!(pairs[0].sort_key.errors, 0);
        assert!(pairs[1].sort_key.errors > pairs[0].sort_key.errors);
    }

    #[test]
    fn countup_presort_tie_breaks_by_record_id_without_duplicate_key_id() {
        fn tied_pair(id: &str, original_index: usize) -> CountupWorkPair {
            CountupWorkPair {
                input_list_index: 0,
                sort_key: CountupSortKey {
                    errors: 0,
                    total_len: 8,
                    expected_errors: 0.0,
                    numeric_id: 0,
                    original_index,
                },
                r1: record(id, b"ACGTACGT"),
                r2: None,
            }
        }

        let mut pairs = [tied_pair("read_b", 0), tied_pair("read_a", 1)];
        pairs.sort_by(compare_countup_work_pairs);

        assert_eq!(pairs[0].r1.id, "read_a");
        assert_eq!(pairs[1].r1.id, "read_b");
    }

    #[test]
    fn countup_spilled_runs_merge_like_in_memory_sort() {
        fn work_pair(
            id: &str,
            errors: usize,
            len: usize,
            original_index: usize,
        ) -> CountupWorkPair {
            CountupWorkPair {
                input_list_index: 0,
                sort_key: CountupSortKey {
                    errors,
                    total_len: len,
                    expected_errors: errors as f64,
                    numeric_id: original_index as u64,
                    original_index,
                },
                r1: record(id, b"ACGTACGT"),
                r2: None,
            }
        }

        let config = Config::default();
        let mut temp_dir = None;
        let mut run_paths = Vec::new();
        let mut spill_summary = CountupSpillSummary::default();
        let mut first_run = vec![work_pair("worse", 2, 8, 2), work_pair("best", 0, 8, 0)];
        let mut second_run = vec![work_pair("longer", 1, 12, 1), work_pair("shorter", 1, 8, 3)];
        let mut expected = first_run.clone();
        expected.extend(second_run.clone());
        expected.sort_by(compare_countup_work_pairs);

        spill_countup_run(
            &config,
            &mut temp_dir,
            &mut run_paths,
            &mut spill_summary,
            &mut first_run,
        )
        .unwrap();
        spill_countup_run(
            &config,
            &mut temp_dir,
            &mut run_paths,
            &mut spill_summary,
            &mut second_run,
        )
        .unwrap();
        spill_summary.final_runs = run_paths.len();
        let source = CountupWorkSource {
            temp_dir,
            inner: CountupWorkSourceInner::Spilled(run_paths),
        };
        let mut iter = source.into_iter().unwrap();
        let mut actual_ids = Vec::new();
        while let Some(pair) = iter.next_pair().unwrap() {
            actual_ids.push(pair.r1.id);
        }
        let expected_ids: Vec<_> = expected.into_iter().map(|pair| pair.r1.id).collect();

        assert_eq!(actual_ids, expected_ids);
        assert_eq!(actual_ids, ["best", "longer", "shorter", "worse"]);
        assert_eq!(spill_summary.initial_runs, 2);
        assert_eq!(spill_summary.merge_runs, 0);
        assert_eq!(spill_summary.final_runs, 2);
        assert!(spill_summary.bytes_written > 0);
        assert_eq!(
            spill_summary.peak_live_bytes,
            spill_summary.final_live_bytes
        );
    }

    #[test]
    fn countup_spill_live_limit_aborts_initial_run() {
        let config = Config {
            max_countup_spill_live_bytes: Some(0),
            ..Config::default()
        };
        let mut temp_dir = None;
        let mut run_paths = Vec::new();
        let mut spill_summary = CountupSpillSummary::default();
        let mut run = vec![CountupWorkPair {
            input_list_index: 0,
            sort_key: CountupSortKey {
                errors: 0,
                total_len: 8,
                expected_errors: 0.0,
                numeric_id: 0,
                original_index: 0,
            },
            r1: record("read", b"ACGTACGT"),
            r2: None,
        }];

        let err = spill_countup_run(
            &config,
            &mut temp_dir,
            &mut run_paths,
            &mut spill_summary,
            &mut run,
        )
        .unwrap_err()
        .to_string();

        assert!(err.contains("maxcountupspillbytes"), "{err}");
        assert_eq!(spill_summary.initial_runs, 1);
        assert!(spill_summary.peak_live_bytes > 0);
        assert_eq!(run_paths.len(), 1);
    }

    #[test]
    fn countup_spill_final_live_limit_aborts_initial_run() {
        let config = Config {
            max_countup_spill_final_live_bytes: Some(0),
            ..Config::default()
        };
        let mut temp_dir = None;
        let mut run_paths = Vec::new();
        let mut spill_summary = CountupSpillSummary::default();
        let mut run = vec![CountupWorkPair {
            input_list_index: 0,
            sort_key: CountupSortKey {
                errors: 0,
                total_len: 8,
                expected_errors: 0.0,
                numeric_id: 0,
                original_index: 0,
            },
            r1: record("read", b"ACGTACGT"),
            r2: None,
        }];

        let err = spill_countup_run(
            &config,
            &mut temp_dir,
            &mut run_paths,
            &mut spill_summary,
            &mut run,
        )
        .unwrap_err()
        .to_string();

        assert!(err.contains("maxcountupspillfinallivebytes"), "{err}");
        assert_eq!(spill_summary.initial_runs, 1);
        assert!(spill_summary.final_live_bytes > 0);
        assert_eq!(run_paths.len(), 1);
    }

    #[test]
    fn countup_spill_initial_run_limit_aborts_initial_run() {
        let config = Config {
            max_countup_spill_initial_runs: Some(0),
            ..Config::default()
        };
        let mut temp_dir = None;
        let mut run_paths = Vec::new();
        let mut spill_summary = CountupSpillSummary::default();
        let mut run = vec![CountupWorkPair {
            input_list_index: 0,
            sort_key: CountupSortKey {
                errors: 0,
                total_len: 8,
                expected_errors: 0.0,
                numeric_id: 0,
                original_index: 0,
            },
            r1: record("read", b"ACGTACGT"),
            r2: None,
        }];

        let err = spill_countup_run(
            &config,
            &mut temp_dir,
            &mut run_paths,
            &mut spill_summary,
            &mut run,
        )
        .unwrap_err()
        .to_string();

        assert!(err.contains("maxcountupspillinitialruns"), "{err}");
        assert_eq!(spill_summary.initial_runs, 1);
    }

    #[test]
    fn countup_compacted_run_group_preserves_sorted_order() {
        fn work_pair(
            id: &str,
            errors: usize,
            len: usize,
            original_index: usize,
        ) -> CountupWorkPair {
            CountupWorkPair {
                input_list_index: 0,
                sort_key: CountupSortKey {
                    errors,
                    total_len: len,
                    expected_errors: errors as f64,
                    numeric_id: original_index as u64,
                    original_index,
                },
                r1: record(id, b"ACGTACGT"),
                r2: None,
            }
        }

        let dir = tempfile::tempdir().unwrap();
        let mut all_pairs = Vec::new();
        let mut paths = Vec::new();
        for (run_index, mut run) in [
            vec![work_pair("c", 3, 8, 3), work_pair("a", 0, 8, 0)],
            vec![work_pair("d", 4, 8, 4), work_pair("b", 1, 8, 1)],
            vec![work_pair("e", 5, 8, 5), work_pair("aa", 1, 12, 2)],
        ]
        .into_iter()
        .enumerate()
        {
            all_pairs.extend(run.clone());
            run.sort_by(compare_countup_work_pairs);
            let path = dir.path().join(format!("run-{run_index}.bin"));
            write_countup_run(&path, &run).unwrap();
            paths.push(path);
        }
        all_pairs.sort_by(compare_countup_work_pairs);

        let merged = dir.path().join("merged.bin");
        let merged_bytes = merge_countup_run_group(&paths, &merged).unwrap();
        let mut reader = CountupRunReader::open(&merged).unwrap();
        let mut actual_ids = Vec::new();
        while let Some(pair) = reader.next_pair().unwrap() {
            actual_ids.push(pair.r1.id);
        }
        let expected_ids: Vec<_> = all_pairs.into_iter().map(|pair| pair.r1.id).collect();

        assert_eq!(actual_ids, expected_ids);
        assert_eq!(merged_bytes, merged.metadata().unwrap().len());
    }

    #[test]
    fn countup_compaction_tracks_peak_and_final_temp_bytes() {
        let dir = tempfile::tempdir().unwrap();
        let mut paths = Vec::new();
        let mut spill_summary = CountupSpillSummary::default();

        for run_index in 0..=COUNTUP_SORT_MERGE_FANIN {
            let path = dir.path().join(format!("run-{run_index}.bin"));
            let pair = CountupWorkPair {
                input_list_index: 0,
                sort_key: CountupSortKey {
                    errors: run_index,
                    total_len: 8,
                    expected_errors: run_index as f64,
                    numeric_id: run_index as u64,
                    original_index: run_index,
                },
                r1: record(&format!("read-{run_index}"), b"ACGTACGT"),
                r2: None,
            };
            let bytes = write_countup_run(&path, &[pair]).unwrap();
            spill_summary.note_initial_run(bytes);
            paths.push(path);
        }

        let initial_live_bytes = spill_summary.final_live_bytes;
        compact_countup_runs(&Config::default(), &mut paths, &mut spill_summary).unwrap();
        spill_summary.final_runs = paths.len();
        let final_live_from_files: u64 = paths
            .iter()
            .map(|path| path.metadata().unwrap().len())
            .sum();

        assert_eq!(spill_summary.initial_runs, COUNTUP_SORT_MERGE_FANIN + 1);
        assert_eq!(spill_summary.merge_runs, 2);
        assert_eq!(spill_summary.final_runs, 2);
        assert_eq!(spill_summary.final_live_bytes, final_live_from_files);
        assert!(spill_summary.bytes_written > initial_live_bytes);
        assert!(spill_summary.peak_live_bytes >= initial_live_bytes);
    }

    #[test]
    fn countup_spill_write_limit_aborts_compaction() {
        let dir = tempfile::tempdir().unwrap();
        let mut paths = Vec::new();
        let mut spill_summary = CountupSpillSummary::default();

        for run_index in 0..=COUNTUP_SORT_MERGE_FANIN {
            let path = dir.path().join(format!("run-{run_index}.bin"));
            let pair = CountupWorkPair {
                input_list_index: 0,
                sort_key: CountupSortKey {
                    errors: run_index,
                    total_len: 8,
                    expected_errors: run_index as f64,
                    numeric_id: run_index as u64,
                    original_index: run_index,
                },
                r1: record(&format!("read-{run_index}"), b"ACGTACGT"),
                r2: None,
            };
            let bytes = write_countup_run(&path, &[pair]).unwrap();
            spill_summary.note_initial_run(bytes);
            paths.push(path);
        }
        let config = Config {
            max_countup_spill_write_bytes: Some(spill_summary.bytes_written),
            ..Config::default()
        };

        let err = compact_countup_runs(&config, &mut paths, &mut spill_summary)
            .unwrap_err()
            .to_string();

        assert!(err.contains("maxcountupspillwritebytes"), "{err}");
        assert!(spill_summary.merge_runs > 0);
        assert!(spill_summary.bytes_written > config.max_countup_spill_write_bytes.unwrap());
    }

    #[test]
    fn countup_spill_run_limits_abort_compaction() {
        let dir = tempfile::tempdir().unwrap();
        let mut paths = Vec::new();
        let mut spill_summary = CountupSpillSummary::default();

        for run_index in 0..=COUNTUP_SORT_MERGE_FANIN {
            let path = dir.path().join(format!("run-{run_index}.bin"));
            let pair = CountupWorkPair {
                input_list_index: 0,
                sort_key: CountupSortKey {
                    errors: run_index,
                    total_len: 8,
                    expected_errors: run_index as f64,
                    numeric_id: run_index as u64,
                    original_index: run_index,
                },
                r1: record(&format!("read-{run_index}"), b"ACGTACGT"),
                r2: None,
            };
            let bytes = write_countup_run(&path, &[pair]).unwrap();
            spill_summary.note_initial_run(bytes);
            paths.push(path);
        }
        let merge_limited = Config {
            max_countup_spill_merge_runs: Some(0),
            ..Config::default()
        };
        let mut merge_limited_paths = paths.clone();
        let err =
            compact_countup_runs(&merge_limited, &mut merge_limited_paths, &mut spill_summary)
                .unwrap_err()
                .to_string();
        assert!(err.contains("maxcountupspillmergeruns"), "{err}");

        let mut spill_summary = CountupSpillSummary::default();
        let mut paths = Vec::new();
        for run_index in 0..=COUNTUP_SORT_MERGE_FANIN {
            let path = dir.path().join(format!("final-run-{run_index}.bin"));
            let pair = CountupWorkPair {
                input_list_index: 0,
                sort_key: CountupSortKey {
                    errors: run_index,
                    total_len: 8,
                    expected_errors: run_index as f64,
                    numeric_id: run_index as u64,
                    original_index: run_index,
                },
                r1: record(&format!("final-read-{run_index}"), b"ACGTACGT"),
                r2: None,
            };
            let bytes = write_countup_run(&path, &[pair]).unwrap();
            spill_summary.note_initial_run(bytes);
            paths.push(path);
        }
        let final_limited = Config {
            max_countup_spill_final_runs: Some(1),
            ..Config::default()
        };
        let err = compact_countup_runs(&final_limited, &mut paths, &mut spill_summary)
            .unwrap_err()
            .to_string();
        assert!(err.contains("maxcountupspillfinalruns"), "{err}");
    }

    #[test]
    fn countup_run_reader_uses_large_spill_buffer() {
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("run.bin");
        let pair = CountupWorkPair {
            input_list_index: 0,
            sort_key: CountupSortKey {
                errors: 0,
                total_len: 8,
                expected_errors: 0.0,
                numeric_id: 0,
                original_index: 0,
            },
            r1: record("read", b"ACGTACGT"),
            r2: None,
        };

        write_countup_run(&path, &[pair]).unwrap();
        let reader = CountupRunReader::open(&path).unwrap();

        assert_eq!(reader.reader.capacity(), COUNTUP_RUN_IO_BUFFER_CAPACITY);
    }

    #[test]
    fn countup_work_pair_memory_hint_tracks_payload_size() {
        let small = CountupWorkPair {
            input_list_index: 0,
            sort_key: CountupSortKey {
                errors: 0,
                total_len: 4,
                expected_errors: 0.0,
                numeric_id: 0,
                original_index: 0,
            },
            r1: record("small", b"ACGT"),
            r2: None,
        };
        let large = CountupWorkPair {
            input_list_index: 0,
            sort_key: CountupSortKey {
                errors: 0,
                total_len: 400,
                expected_errors: 0.0,
                numeric_id: 1,
                original_index: 1,
            },
            r1: record("large", &vec![b'A'; 400]),
            r2: Some(record("large/2", &vec![b'C'; 400])),
        };

        assert!(countup_work_pair_memory_hint(&large) > countup_work_pair_memory_hint(&small));
    }

    #[test]
    fn countup_work_candidate_memory_hint_tracks_payload_size() {
        let small = CountupWorkCandidate {
            input_list_index: 0,
            original_index: 0,
            rand: 0.0,
            r1: record("small", b"ACGT"),
            r2: None,
        };
        let large = CountupWorkCandidate {
            input_list_index: 0,
            original_index: 1,
            rand: 0.0,
            r1: record("large", &vec![b'A'; 400]),
            r2: Some(record("large/2", &vec![b'C'; 400])),
        };

        assert!(
            countup_work_candidate_memory_hint(&large) > countup_work_candidate_memory_hint(&small)
        );
    }

    #[test]
    fn countup_prepass_chunk_ready_respects_pair_and_byte_limits() {
        assert!(!countup_prepass_chunk_ready(
            COUNTUP_PREPASS_CHUNK_PAIR_LIMIT - 1,
            COUNTUP_PREPASS_CHUNK_BYTE_LIMIT - 1
        ));
        assert!(countup_prepass_chunk_ready(
            COUNTUP_PREPASS_CHUNK_PAIR_LIMIT,
            0
        ));
        assert!(countup_prepass_chunk_ready(
            1,
            COUNTUP_PREPASS_CHUNK_BYTE_LIMIT
        ));
    }

    #[test]
    fn countup_prepass_carries_tossed_reads_only_with_abrc() {
        let config = Config {
            k: 3,
            min_length: 11,
            target_depth: 2,
            min_depth: 1,
            min_kmers_over_min_depth: 3,
            ..Config::default()
        };
        let prepass = countup_prepass_config(&config);
        assert_eq!(prepass.target_depth, 8);
        assert_eq!(prepass.min_depth, 0);
        assert_eq!(prepass.min_kmers_over_min_depth, 1);

        let input_counts = CountMap::default();
        let mut filtered = record("short", b"AAAAAAAAAA");
        assert!(
            !countup_prepass_pair(&prepass, false, &input_counts, &mut filtered, None, 0.0,)
                .include
        );

        let mut carried = record("short", b"AAAAAAAAAA");
        assert!(
            countup_prepass_pair(&prepass, true, &input_counts, &mut carried, None, 0.0,).include
        );
    }

    #[test]
    fn countup_prepass_requires_both_mates_bad_like_java() {
        let config = Config {
            count_up: true,
            toss_error_reads: true,
            require_both_bad: false,
            k: 3,
            target_depth: 100,
            max_depth: Some(1000),
            min_depth: 1,
            min_kmers_over_min_depth: 1,
            error_detect_ratio: 2,
            high_thresh: 2,
            low_thresh: 1,
            ..Config::default()
        };
        let prepass = countup_prepass_config(&config);
        assert!(!config.require_both_bad);
        assert!(prepass.require_both_bad);

        let mut bad_mate = record("bad", b"AAACCC");
        let mut good_mate = record("good", b"GGGGGG");
        let mut input_counts = CountMap::default();
        let bad_keys = kmers_for_record(&bad_mate, &prepass);
        for key in &bad_keys {
            input_counts.insert(key.clone(), 10);
        }
        input_counts.insert(bad_keys[1].clone(), 1);
        input_counts.insert(bad_keys[2].clone(), 1);
        for key in kmers_for_record(&good_mate, &prepass) {
            input_counts.insert(key, 10);
        }

        assert!(analyze_pair(&prepass, &input_counts, &bad_mate, None).error1);
        assert!(!analyze_pair(&prepass, &input_counts, &good_mate, None).error1);
        assert!(
            countup_prepass_pair(
                &prepass,
                false,
                &input_counts,
                &mut bad_mate,
                Some(&mut good_mate),
                0.0,
            )
            .include
        );
    }

    #[test]
    fn countup_prepass_reuses_decision_analysis_for_sort_key_without_ecc() {
        let config = Config {
            count_up: true,
            k: 3,
            min_depth: 1,
            min_kmers_over_min_depth: 1,
            target_depth: 100,
            max_depth: Some(1000),
            ..Config::default()
        };
        let prepass = countup_prepass_config(&config);
        let mut read = record("read42", b"ACGTACGT");
        let mut input_counts = CountMap::default();
        for key in kmers_for_record(&read, &prepass) {
            input_counts.insert(key, 10);
        }

        let result = countup_prepass_pair(&prepass, false, &input_counts, &mut read, None, 0.0);
        let reused_key =
            countup_sort_key_from_analysis(&read, None, 42, result.sort_analysis.as_ref().unwrap());
        let replayed_key = countup_sort_key(&prepass, &input_counts, &read, None, 42);

        assert!(result.include);
        assert_eq!(reused_key.errors, replayed_key.errors);
        assert_eq!(reused_key.total_len, replayed_key.total_len);
        assert_eq!(reused_key.numeric_id, replayed_key.numeric_id);
        assert_eq!(reused_key.original_index, replayed_key.original_index);
        assert_eq!(reused_key.expected_errors, replayed_key.expected_errors);
    }

    #[test]
    fn countup_work_candidates_match_sequential_prepass_sort_keys() {
        let config = Config {
            count_up: true,
            k: 3,
            min_depth: 1,
            min_kmers_over_min_depth: 1,
            target_depth: 100,
            max_depth: Some(1000),
            ..Config::default()
        };
        let prepass = countup_prepass_config(&config);
        let clean = record("clean", b"ACGTACGT");
        let noisy = record("noisy", b"AAAACCCC");
        let mut input_counts = CountMap::default();
        for key in kmers_for_record(&clean, &prepass) {
            input_counts.insert(key, 10);
        }
        let candidates = vec![
            CountupWorkCandidate {
                input_list_index: 0,
                original_index: 0,
                rand: 0.0,
                r1: noisy.clone(),
                r2: None,
            },
            CountupWorkCandidate {
                input_list_index: 0,
                original_index: 1,
                rand: 0.0,
                r1: clean.clone(),
                r2: None,
            },
        ];
        let mut actual =
            process_countup_work_candidates(&config, &prepass, &input_counts, candidates);
        let mut expected = vec![
            CountupWorkPair {
                input_list_index: 0,
                sort_key: countup_sort_key(&prepass, &input_counts, &noisy, None, 0),
                r1: noisy,
                r2: None,
            },
            CountupWorkPair {
                input_list_index: 0,
                sort_key: countup_sort_key(&prepass, &input_counts, &clean, None, 1),
                r1: clean,
                r2: None,
            },
        ];
        actual.sort_by(compare_countup_work_pairs);
        expected.sort_by(compare_countup_work_pairs);

        let actual_ids: Vec<_> = actual.iter().map(|pair| pair.r1.id.as_str()).collect();
        let expected_ids: Vec<_> = expected.iter().map(|pair| pair.r1.id.as_str()).collect();
        assert_eq!(actual_ids, expected_ids);
        for (actual, expected) in actual.iter().zip(&expected) {
            assert_eq!(actual.sort_key.errors, expected.sort_key.errors);
            assert_eq!(actual.sort_key.total_len, expected.sort_key.total_len);
            assert_eq!(
                actual.sort_key.original_index,
                expected.sort_key.original_index
            );
        }
    }

    #[test]
    fn countup_length_filter_respects_keepall_override() {
        let read = record("short", b"ACGT");
        let filter_config = Config {
            min_length: 5,
            ..Config::default()
        };
        assert!(countup_length_toss(&filter_config, &read, None));

        let keepall_config = Config {
            keep_all: true,
            ..filter_config
        };
        assert!(!countup_length_toss(&keepall_config, &read, None));
    }

    #[test]
    fn countup_tossbadreads_applies_java_error_spike_rules() {
        let keys: Vec<_> = (0..20).map(KmerKey::Short).collect();
        let mut input_counts = CountMap::default();
        let mut kept_counts = CountMap::default();
        for (index, key) in keys.iter().enumerate() {
            let input_depth = if index < 8 { 1 } else { 10 };
            let kept_depth = if index < 8 { 0 } else { 10 };
            input_counts.insert(key.clone(), input_depth);
            kept_counts.insert(key.clone(), kept_depth);
        }
        let base_config = Config {
            min_depth: 1,
            min_kmers_over_min_depth: 1,
            target_depth: 10,
            low_thresh: 1,
            high_thresh: 10,
            error_detect_ratio: 2,
            ..Config::default()
        };

        assert!(!decide_countup_pair(
            &base_config,
            &input_counts,
            &kept_counts,
            &keys,
            10,
        ));

        let toss_errors_config = Config {
            toss_error_reads: true,
            ..base_config.clone()
        };
        assert!(decide_countup_pair(
            &toss_errors_config,
            &input_counts,
            &kept_counts,
            &keys,
            10,
        ));

        let keepall_config = Config {
            keep_all: true,
            ..toss_errors_config
        };
        assert!(!decide_countup_pair(
            &keepall_config,
            &input_counts,
            &kept_counts,
            &keys,
            10,
        ));
    }

    #[test]
    fn java_rng_matches_known_first_doubles() {
        let mut rng = JavaXoshiro::new(0);
        let values = [
            rng.next_double(),
            rng.next_double(),
            rng.next_double(),
            rng.next_double(),
        ];
        let expected = [
            0.02774461029305808,
            0.9419058303890074,
            0.3687890049137593,
            0.8390756877056451,
        ];
        for (actual, expected) in values.into_iter().zip(expected) {
            assert!((actual - expected).abs() < f64::EPSILON);
        }
    }

    #[test]
    fn nondeterministic_seed_varies_between_requests() {
        let first = nondeterministic_seed();
        let second = nondeterministic_seed();
        assert_ne!(first, second);
    }

    #[test]
    fn deterministic_coin_uses_java_read_rand_shape() {
        assert_eq!(deterministic_coin(Some(0.0), 7), 1);
        assert_eq!(deterministic_coin(Some(0.5), 7), 4);
        assert_eq!(deterministic_coin(Some(0.999_999), 7), 7);
    }

    #[test]
    fn qtrim_right_uses_java_optimal_quality_scoring() {
        let config = Config {
            trim_right: true,
            trim_quality: 10.0,
            ..Config::default()
        };
        let mut read = quality_record("r1", b"ACGTACGT", b"IIII!!!!");

        trim_record(&config, &mut read);

        assert_eq!(read.bases, b"ACGT");
        assert_eq!(read.qualities.as_deref(), Some(&b"IIII"[..]));
    }

    #[test]
    fn qtrim_left_uses_java_optimal_quality_scoring() {
        let config = Config {
            trim_left: true,
            trim_quality: 10.0,
            ..Config::default()
        };
        let mut read = quality_record("r1", b"ACGTACGT", b"!!!!IIII");

        trim_record(&config, &mut read);

        assert_eq!(read.bases, b"ACGT");
        assert_eq!(read.qualities.as_deref(), Some(&b"IIII"[..]));
    }

    #[test]
    fn ecc_corrects_single_substitution_from_exact_counts() {
        let clean = b"ACGTTGCATGTCAGTACCGTAACGTTGCA";
        let mut mutant = clean.to_vec();
        mutant[14] = b'A';
        assert_ne!(mutant, clean);

        let config = Config {
            k: 7,
            min_quality: 0,
            min_prob: 0.0,
            error_correct: true,
            passes: 1,
            ..Config::default()
        };
        let mut counts = CountMap::default();
        for i in 0..30 {
            increment_pair_counts(
                &config,
                &mut counts,
                &record(&format!("clean{i}"), clean),
                None,
            );
        }
        increment_pair_counts(&config, &mut counts, &record("mutant", &mutant), None);

        let mut read = record("mutant", &mutant);
        let result = correct_read_errors(&config, &counts, &mut read);

        assert_eq!(result.corrected, 1);
        assert!(!result.uncorrectable);
        assert_eq!(read.bases, clean);
    }

    #[test]
    fn ecc_flags_high_quality_suspect_error_as_uncorrectable() {
        let clean = b"ACGTTGCATGTCAGTACCGTAACGTTGCA";
        let mut mutant = clean.to_vec();
        mutant[14] = b'A';
        let config = Config {
            k: 7,
            min_quality: 0,
            min_prob: 0.0,
            error_correct: true,
            max_quality_to_correct: 0,
            passes: 1,
            ..Config::default()
        };
        let mut counts = CountMap::default();
        for i in 0..30 {
            increment_pair_counts(
                &config,
                &mut counts,
                &record(&format!("clean{i}"), clean),
                None,
            );
        }
        increment_pair_counts(&config, &mut counts, &record("mutant", &mutant), None);

        let mut read = record("mutant", &mutant);
        let result = correct_read_errors(&config, &counts, &mut read);

        assert_eq!(result.corrected, 0);
        assert!(result.uncorrectable);
        assert_eq!(read.bases, mutant);
    }

    #[test]
    fn ecc_pair_rollback_restores_corrected_mate_when_partner_is_uncorrectable() {
        let clean = b"ACGTTGCATGTCAGTACCGTAACGTTGCA";
        let mut mutant = clean.to_vec();
        mutant[14] = b'A';
        let config = Config {
            k: 7,
            min_quality: 0,
            min_prob: 0.0,
            error_correct: true,
            max_quality_to_correct: 20,
            passes: 1,
            ..Config::default()
        };
        let mut counts = CountMap::default();
        for i in 0..30 {
            increment_pair_counts(
                &config,
                &mut counts,
                &record(&format!("clean{i}"), clean),
                None,
            );
        }
        increment_pair_counts(&config, &mut counts, &record("mutant", &mutant), None);

        let low_quality = vec![b'!'; mutant.len()];
        let high_quality = vec![b'I'; mutant.len()];
        let mut correctable = quality_record("lowq", &mutant, &low_quality);
        let mut uncorrectable = quality_record("highq", &mutant, &high_quality);
        let original_correctable = correctable.clone();
        let original_uncorrectable = uncorrectable.clone();

        let result = correct_pair_errors_with_rollback(
            &config,
            &counts,
            &mut correctable,
            Some(&mut uncorrectable),
        );

        assert!(result.corrected > 0);
        assert!(result.uncorrectable);
        assert_eq!(correctable.bases, original_correctable.bases);
        assert_eq!(
            correctable.qualities.as_deref(),
            original_correctable.qualities.as_deref()
        );
        assert_eq!(uncorrectable.bases, original_uncorrectable.bases);
    }

    #[test]
    fn ecc_marks_uncorrectable_errors_when_requested() {
        let clean = b"ACGTTGCATGTCAGTACCGTAACGTTGCA";
        let mut mutant = clean.to_vec();
        mutant[14] = b'A';
        let config = Config {
            k: 7,
            min_quality: 0,
            min_prob: 0.0,
            error_correct: true,
            max_quality_to_correct: 0,
            mark_uncorrectable_errors: true,
            passes: 1,
            ..Config::default()
        };
        let mut counts = CountMap::default();
        for i in 0..30 {
            increment_pair_counts(
                &config,
                &mut counts,
                &record(&format!("clean{i}"), clean),
                None,
            );
        }
        increment_pair_counts(&config, &mut counts, &record("mutant", &mutant), None);

        let mut read = record("mutant", &mutant);
        let result = correct_read_errors(&config, &counts, &mut read);

        assert_eq!(result.corrected, 0);
        assert_eq!(result.marked, 1);
        assert!(result.uncorrectable);
        assert_eq!(read.bases, mutant);
        assert_eq!(read.qualities.as_ref().unwrap()[14], b'2');
    }

    #[test]
    fn ecc_mark_only_reduces_suspect_base_quality() {
        let clean = b"ACGTTGCATGTCAGTACCGTAACGTTGCA";
        let mut mutant = clean.to_vec();
        mutant[14] = b'A';
        let config = Config {
            k: 7,
            min_quality: 0,
            min_prob: 0.0,
            error_correct: true,
            mark_errors_only: true,
            passes: 1,
            ..Config::default()
        };
        let mut counts = CountMap::default();
        for i in 0..30 {
            increment_pair_counts(
                &config,
                &mut counts,
                &record(&format!("clean{i}"), clean),
                None,
            );
        }
        increment_pair_counts(&config, &mut counts, &record("mutant", &mutant), None);

        let mut read = record("mutant", &mutant);
        let result = correct_read_errors(&config, &counts, &mut read);

        assert_eq!(result.marked, 1);
        assert_eq!(read.bases, mutant);
        assert_eq!(read.qualities.as_ref().unwrap()[14], b'2');
    }

    #[test]
    fn ecc_mark_only_marks_all_detected_sites_even_when_ecclimit_is_low() {
        let config = Config {
            k: 7,
            prefix_len: 2,
            max_errors_to_correct: 1,
            correct_from_right: false,
            ..Config::default()
        };
        let mut read = quality_record("marked", b"ACGTTGCATGTC", b"IIIIIIIIIIII");
        let coverage = vec![30, 30, 0, 30, 30, 0];

        let result = mark_read_errors(&config, &mut read, &coverage);

        assert_eq!(result.marked, 2);
        let qualities = read.qualities.as_deref().unwrap();
        assert_eq!(qualities[8], b'2');
        assert_eq!(qualities[11], b'2');
    }

    #[test]
    fn overlap_ecc_repairs_lower_quality_mate_base() {
        let r1_bases = b"TTAGTTGTGCCGCAGCGAAGTAGTGCTTGAAATATGCGAC";
        let r2_clean = b"GTCGCATATTTCAAGCACTACTTCGCTGCGGCACAACTAA";
        let mut r2_bases = r2_clean.to_vec();
        r2_bases[20] = b'A';
        let mut r1 = quality_record("r1", r1_bases, b"IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII");
        let mut r2 = quality_record("r2", &r2_bases, b"IIIIIIIIIIIIIIIIIIII#IIIIIIIIIIIIIIIIIII");
        let config = Config {
            overlap_error_correct: true,
            max_quality_to_correct: 20,
            ..Config::default()
        };

        let result = correct_pair_by_overlap(&config, &mut r1, &mut r2);

        assert_eq!(result.corrected, 1);
        assert_eq!(r1.bases, r1_bases);
        assert_eq!(r2.bases, r2_clean);
        assert_eq!(
            r1.qualities.as_deref(),
            Some(&b"SSSSSSSSSSSSSSSSSSSGSSSSSSSSSSSSSSSSSSSS"[..])
        );
        assert_eq!(
            r2.qualities.as_deref(),
            Some(&b"SSSSSSSSSSSSSSSSSSSSGSSSSSSSSSSSSSSSSSSS"[..])
        );
    }

    #[test]
    fn overlap_ecc_skips_short_pairs_like_java_strict_mode() {
        let r1_bases = b"ACGTTGCATGTCAGTA";
        let r2_clean = b"TACTGACATGCAACGT";
        let mut r2_bases = r2_clean.to_vec();
        r2_bases[9] = b'T';
        let mut r1 = quality_record("r1", r1_bases, b"IIIIIIIIIIIIIIII");
        let mut r2 = quality_record("r2", &r2_bases, b"IIIIIIIII!IIIIII");
        let config = Config {
            overlap_error_correct: true,
            max_quality_to_correct: 20,
            ..Config::default()
        };

        let result = correct_pair_by_overlap(&config, &mut r1, &mut r2);

        assert_eq!(result.corrected, 0);
        assert_eq!(r1.bases, r1_bases);
        assert_eq!(r2.bases, r2_bases);
    }

    #[test]
    fn overlap_ecc_skips_ambiguous_repetitive_pairs_like_java_strict_mode() {
        let r1_bases = b"ACGTTGCATGTCAGTAACGTTGCATGTCAGTAACGTTGCA";
        let r2_clean = b"TGCAACGTTACTGACATGCAACGTTACTGACATGCAACGT";
        let mut r2_bases = r2_clean.to_vec();
        r2_bases[20] = b'C';
        let mut r1 = quality_record("r1", r1_bases, b"IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII");
        let mut r2 = quality_record("r2", &r2_bases, b"IIIIIIIIIIIIIIIIIIII!IIIIIIIIIIIIIIIIIII");
        let config = Config {
            overlap_error_correct: true,
            max_quality_to_correct: 20,
            ..Config::default()
        };

        let result = correct_pair_by_overlap(&config, &mut r1, &mut r2);

        assert_eq!(result.corrected, 0);
        assert_eq!(r1.bases, r1_bases);
        assert_eq!(r2.bases, r2_bases);
    }

    #[test]
    fn overlap_entropy_gate_keeps_java_strict_floor_for_high_entropy_fixture() {
        let bases = b"TTAGTTGTGCCGCAGCGAAGTAGTGCTTGAAATATGCGAC";
        assert_eq!(overlap_entropy_min_overlap(bases), 12);
    }

    #[test]
    fn overlap_entropy_gate_raises_min_overlap_for_low_complexity_reads() {
        let bases = b"AAAAAAAAAACCCCCCCCCCGGGGGGGGGGTTTTTTTTTT";
        assert_eq!(overlap_entropy_min_overlap(bases), 32);
    }

    #[test]
    fn overlap_ecc_rejects_high_confidence_mismatch_like_java_expected_filter() {
        let r1_bases = b"TTAGTTGTGCCGCAGCGAAGTAGTGCTTGAAATATGCGAC";
        let mut r2_bases = b"GTCGCATATTTCAAGCACTACTTCGCTGCGGCACAACTAA".to_vec();
        r2_bases[20] = b'A';
        let mut r1 = quality_record("r1", r1_bases, b"IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII");
        let mut r2 = quality_record("r2", &r2_bases, b"IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII");
        let config = Config {
            overlap_error_correct: true,
            max_quality_to_correct: 41,
            ..Config::default()
        };

        let result = correct_pair_by_overlap(&config, &mut r1, &mut r2);

        assert_eq!(result.corrected, 0);
        assert_eq!(r1.bases, r1_bases);
        assert_eq!(r2.bases, r2_bases);
        assert_eq!(
            r1.qualities.as_deref(),
            Some(&b"IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII"[..])
        );
        assert_eq!(
            r2.qualities.as_deref(),
            Some(&b"IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII"[..])
        );
    }

    #[test]
    fn overlap_ecc_rejects_low_confidence_tie_under_java_strict_mode() {
        let r1_bases = b"TTAGTTGTGCCGCAGCGAAGTAGTGCTTGAAATATGCGAC";
        let mut r2_bases = b"GTCGCATATTTCAAGCACTACTTCGCTGCGGCACAACTAA".to_vec();
        r2_bases[20] = b'A';
        let mut r1 = quality_record("r1", r1_bases, b"!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!");
        let mut r2 = quality_record("r2", &r2_bases, b"!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!");
        let config = Config {
            overlap_error_correct: true,
            max_quality_to_correct: 41,
            ..Config::default()
        };

        let result = correct_pair_by_overlap(&config, &mut r1, &mut r2);

        assert_eq!(result.corrected, 0);
        assert_eq!(r1.bases, r1_bases);
        assert_eq!(r2.bases, r2_bases);
        assert_eq!(
            r1.qualities.as_deref(),
            Some(&b"!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"[..])
        );
        assert_eq!(
            r2.qualities.as_deref(),
            Some(&b"!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"[..])
        );
    }

    #[test]
    fn overlap_ecc_rejects_quality_weighted_multimismatch_candidate_like_java() {
        let r1_bases = b"CAGTAACCAATGCCTGTTGAGATGCCAGACGCGTAACCAAAA";
        let r2_bases = b"TTTTGCTAACGCGTCTGGCATCTCAACAGGCATTGGTTAC";
        let mut r1 = quality_record(
            "r1",
            r1_bases,
            b"IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII",
        );
        let mut r2 = quality_record("r2", r2_bases, b"IIIII!I'IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII");
        let original_r1 = r1.clone();
        let original_r2 = r2.clone();
        let config = Config {
            overlap_error_correct: true,
            max_quality_to_correct: 41,
            ..Config::default()
        };

        let result = correct_pair_by_overlap(&config, &mut r1, &mut r2);

        assert_eq!(result.corrected, 0);
        assert_eq!(r1.bases, original_r1.bases);
        assert_eq!(r2.bases, original_r2.bases);
    }

    #[test]
    fn trim_after_marking_defers_qtrim_until_after_ecc_marking() {
        let clean = b"ACGTTGCATGTCAGTACCGTAACGTTGCA";
        let mut mutant = clean.to_vec();
        mutant[26] = b'A';
        let config = Config {
            k: 7,
            min_quality: 0,
            min_prob: 0.0,
            error_correct: true,
            mark_errors_only: true,
            trim_after_marking: true,
            trim_right: true,
            trim_optimal: false,
            trim_quality: 20.0,
            keep_all: true,
            passes: 1,
            ..Config::default()
        };
        let mut counts = CountMap::default();
        for i in 0..30 {
            increment_pair_counts(
                &config,
                &mut counts,
                &record(&format!("clean{i}"), clean),
                None,
            );
        }
        increment_pair_counts(&config, &mut counts, &record("mutant", &mutant), None);

        let input = vec![(0, record("mutant", &mutant), None, 0.0)];
        let pairs = normalize_pair_chunk(&config, &counts, &input);

        assert_eq!(pairs[0].out_r1.bases, b"ACGTTGCATGTCAGTACCGTAACGTTAC");
        assert_eq!(
            pairs[0].out_r1.qualities.as_deref(),
            Some(&b"IIIIIIIIIIIIIIIIIIIIIIIIIIII"[..])
        );
    }

    #[test]
    fn bad_kmer_fraction_lowers_dynamic_toss_target_like_bbnorm() {
        let config = Config {
            target_depth: 100,
            max_depth: Some(125),
            target_bad_percent_low: 0.2,
            target_bad_percent_high: 0.8,
            ..Config::default()
        };
        let clean = PairAnalysis::default();
        assert_eq!(dynamic_depth_limits(&config, &clean), (100, 125));

        let noisy = PairAnalysis {
            low_kmer_count: 5,
            total_kmer_count: 10,
            ..PairAnalysis::default()
        };
        assert_eq!(dynamic_depth_limits(&config, &noisy), (35, 35));
    }

    #[test]
    fn multipass_bad_depth_targets_match_java_pass_shape() {
        let config = Config {
            passes: 3,
            target_depth: 100,
            target_bad_percent_low: 0.2,
            target_bad_percent_high: 0.8,
            ..Config::default()
        };

        let first_target = intermediate_target_depth(&config, 1);
        assert_eq!(first_target, 400);
        assert_eq!(
            intermediate_bad_depth_targets(&config, 1, first_target),
            (30, 120)
        );

        let second_target = intermediate_target_depth(&config, 2);
        assert_eq!(second_target, 200);
        assert_eq!(
            intermediate_bad_depth_targets(&config, 2, second_target),
            (20, 80)
        );
    }

    #[test]
    fn qtrim_keeps_java_min_result_shape_for_all_bad_reads() {
        let config = Config {
            trim_right: true,
            trim_quality: 10.0,
            ..Config::default()
        };
        let mut read = quality_record("r1", b"ACGT", b"!!!!");

        trim_record(&config, &mut read);

        assert_eq!(read.bases, b"A");
        assert_eq!(read.qualities.as_deref(), Some(&b"!"[..]));
    }

    #[test]
    fn qtrim_window_uses_java_sliding_threshold() {
        let config = Config {
            trim_right: true,
            trim_quality: 10.0,
            trim_optimal: false,
            trim_window: true,
            trim_window_length: 4,
            ..Config::default()
        };
        let mut read = quality_record("r1", b"ACGTACGTACGT", b"IIIIIII!!!!!");

        trim_record(&config, &mut read);

        assert_eq!(read.bases, b"ACGTACG");
        assert_eq!(read.qualities.as_deref(), Some(&b"IIIIIII"[..]));
    }

    #[test]
    fn output_hash_patterns_match_bbnorm_pair_expansion() {
        let paths = prepare_output_paths(Some(Path::new("reads#.fq")), None, true);
        assert_eq!(paths.first, Some(PathBuf::from("reads1.fq")));
        assert_eq!(paths.second, Some(PathBuf::from("reads2.fq")));

        let paths = prepare_output_paths(
            Some(Path::new("reads#.fq")),
            Some(Path::new("mate.fq")),
            true,
        );
        assert_eq!(paths.first, Some(PathBuf::from("reads1.fq")));
        assert_eq!(paths.second, Some(PathBuf::from("mate.fq")));

        let paths = prepare_output_paths(Some(Path::new("single#.fq")), None, false);
        assert_eq!(paths.first, Some(PathBuf::from("single1.fq")));
        assert_eq!(paths.second, None);
    }
}

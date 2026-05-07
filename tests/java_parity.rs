use flate2::read::MultiGzDecoder;
use std::ffi::OsString;
use std::fs::{self, File};
use std::io::{BufRead, BufReader, Write};
use std::path::{Path, PathBuf};
use std::process::Command;
use tempfile::tempdir;

const SAMPLE1: &str = "vendor/BBTools-master/resources/sample1.fq.gz";
const SAMPLE2: &str = "vendor/BBTools-master/resources/sample2.fq.gz";
const SAMPLE_HASH: &str = "vendor/BBTools-master/resources/sample#.fq.gz";
const BBTOOLS_CP: &str = "vendor/BBTools-master/current";

#[test]
fn real_phi_x_pair_matches_java_bbnorm_for_single_pass_keepall_hist() {
    assert_keepall_hist_matches_java_with_extra(&[]);
}

#[test]
fn real_phi_x_pair_matches_java_bbnorm_for_histcol_one_header() {
    assert_keepall_hist_matches_java_with_extra(&["histcol=1"]);
}

#[test]
fn real_phi_x_pair_matches_java_bbnorm_for_histcol_two() {
    assert_keepall_hist_matches_java_with_extra(&["histcol=2"]);
}

#[test]
fn real_phi_x_pair_matches_java_bbnorm_for_zerobin_histogram() {
    assert_keepall_hist_matches_java_with_extra(&["zerobin=t"]);
}

#[test]
fn real_phi_x_pair_matches_java_bbnorm_for_print_zero_coverage_histogram() {
    assert_keepall_hist_matches_java_with_extra(&["printzerocoverage=t", "histlen=5"]);
}

#[test]
fn real_phi_x_pair_matches_java_bbnorm_for_reads_limit() {
    assert_keepall_hist_matches_java_with_extra(&["reads=10"]);
}

#[test]
fn real_phi_x_pair_matches_java_bbnorm_for_tablereads_limit() {
    assert_keepall_hist_matches_java_with_extra(&["tablereads=10"]);
}

#[test]
fn real_phi_x_pair_matches_java_bbnorm_for_kmg_reads_limit() {
    assert_keepall_hist_matches_java_with_extra(&["reads=0.01k"]);
}

#[test]
fn real_phi_x_pair_matches_java_bbnorm_for_kmg_tablereads_limit() {
    assert_keepall_hist_matches_java_with_extra(&["tablereads=0.01k"]);
}

#[test]
fn real_phi_x_pair_matches_java_bbnorm_for_buildstepsize_one() {
    assert_keepall_hist_matches_java_with_extra(&["buildstepsize=1"]);
}

#[test]
fn real_phi_x_pair_matches_java_bbnorm_for_stepsize_two_alias() {
    assert_keepall_hist_matches_java_with_extra(&["stepsize=2"]);
}

#[test]
fn real_phi_x_pair_matches_java_bbnorm_for_default_prefilter() {
    assert_keepall_hist_matches_java_with_extra(&["prefilter=t"]);
}

#[test]
fn real_phi_x_pair_matches_java_bbnorm_for_bare_prefilter() {
    assert_keepall_hist_matches_java_with_extra(&["prefilter"]);
}

#[test]
fn real_phi_x_pair_matches_java_bbnorm_for_runtime_noop_controls() {
    for extra in [["ordered=f"], ["verbose=t"], ["printcoverage=t"]] {
        assert_keepall_hist_matches_java_with_extra(&extra);
    }
}

#[test]
fn representative_header_trimming_controls_match_java_bbnorm() {
    require_file("vendor/BBTools-master/current/jgi/KmerNormalize.class");
    require_java();

    let dir = tempdir().expect("create tempdir");
    let input = dir.path().join("header_trim_input.fq");
    write_fastq_fixture(&input, "read1 comment here", "ACGTACGT");

    for (label, extra) in [
        ("trd_true", &["trd=t"][..]),
        ("trc_true", &["trc=t"][..]),
        ("trimreaddescription_true", &["trimreaddescription=t"][..]),
        ("trimreaddescriptions_true", &["trimreaddescriptions=t"][..]),
        ("trd_false", &["trd=f"][..]),
        ("trimrname_true", &["trimrname=t"][..]),
    ] {
        let java_keep = dir.path().join(format!("java.{label}.keep.fq"));
        let rust_keep = dir.path().join(format!("rust.{label}.keep.fq"));

        let java_args = quality_alias_keepall_args(&input, &java_keep, extra);
        let java_status = Command::new("java")
            .arg("-Xmx1g")
            .arg("-cp")
            .arg(BBTOOLS_CP)
            .arg("jgi.KmerNormalize")
            .args(&java_args)
            .output()
            .expect("run java bbnorm");
        assert_success("java bbnorm", &java_status);

        let rust_args = quality_alias_keepall_args(&input, &rust_keep, extra);
        let rust_status = Command::new(env!("CARGO_BIN_EXE_bbnorm-rs"))
            .args(&rust_args)
            .output()
            .expect("run bbnorm-rs");
        assert_success("bbnorm-rs", &rust_status);

        assert_same_file(&java_keep, &rust_keep);
    }
}

#[test]
fn real_phi_x_pair_matches_java_bbnorm_for_shared_io_runtime_noops() {
    assert_keepall_hist_matches_java_with_extra(&[
        "null",
        "monitor=f",
        "killswitch=600",
        "usejni=f",
        "bytefile1=t",
        "bytefile2=maybe",
        "bf1bufferlen=64k",
        "bfthreads=1",
        "readbufferlength=64k",
        "readbufferdata=1m",
        "readbuffers=1",
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
}

#[test]
fn real_phi_x_pair_matches_java_bbnorm_for_quality_recalibration_noops() {
    assert_keepall_hist_matches_java_with_extra(&[
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
}

#[test]
fn real_phi_x_pair_matches_java_bbnorm_for_shared_environment_runtime_noops() {
    assert_keepall_hist_matches_java_with_extra(&[
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
        "comment=bbnorm-rs",
        "taxpath=auto",
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
}

#[test]
fn real_phi_x_behavior_changing_sketch_controls_fall_back_to_exact_counting() {
    require_file(SAMPLE1);
    require_file(SAMPLE2);
    require_file("vendor/BBTools-master/current/jgi/KmerNormalize.class");
    require_java();

    let dir = tempdir().expect("create tempdir");
    let base_keep1 = dir.path().join("java.base.keep1.fq");
    let base_keep2 = dir.path().join("java.base.keep2.fq");
    let base_hist = dir.path().join("java.base.hist.tsv");
    let base_rhist = dir.path().join("java.base.rhist.tsv");
    let base_low1 = dir.path().join("java.base.low1.fq");
    let base_low2 = dir.path().join("java.base.low2.fq");
    let base_mid1 = dir.path().join("java.base.mid1.fq");
    let base_mid2 = dir.path().join("java.base.mid2.fq");
    let base_high1 = dir.path().join("java.base.high1.fq");
    let base_high2 = dir.path().join("java.base.high2.fq");
    let base_paths = OutputPaths {
        keep1: &base_keep1,
        keep2: &base_keep2,
        hist: &base_hist,
        rhist: &base_rhist,
        low1: &base_low1,
        low2: &base_low2,
        mid1: &base_mid1,
        mid2: &base_mid2,
        high1: &base_high1,
        high2: &base_high2,
    };

    let base_status = Command::new("java")
        .arg("-Xmx1g")
        .arg("-cp")
        .arg(BBTOOLS_CP)
        .arg("jgi.KmerNormalize")
        .args(shared_args(base_paths))
        .output()
        .expect("run java bbnorm baseline");
    assert_success("java bbnorm baseline", &base_status);
    let base_hist_bytes = fs::read(&base_hist).expect("read baseline hist");

    for option in ["prehashes=1", "prefilterhashes=1", "buildpasses=2"] {
        let label = option.replace('=', "_");
        let java_keep1 = dir.path().join(format!("java.{label}.keep1.fq"));
        let java_keep2 = dir.path().join(format!("java.{label}.keep2.fq"));
        let java_hist = dir.path().join(format!("java.{label}.hist.tsv"));
        let java_rhist = dir.path().join(format!("java.{label}.rhist.tsv"));
        let java_low1 = dir.path().join(format!("java.{label}.low1.fq"));
        let java_low2 = dir.path().join(format!("java.{label}.low2.fq"));
        let java_mid1 = dir.path().join(format!("java.{label}.mid1.fq"));
        let java_mid2 = dir.path().join(format!("java.{label}.mid2.fq"));
        let java_high1 = dir.path().join(format!("java.{label}.high1.fq"));
        let java_high2 = dir.path().join(format!("java.{label}.high2.fq"));
        let java_paths = OutputPaths {
            keep1: &java_keep1,
            keep2: &java_keep2,
            hist: &java_hist,
            rhist: &java_rhist,
            low1: &java_low1,
            low2: &java_low2,
            mid1: &java_mid1,
            mid2: &java_mid2,
            high1: &java_high1,
            high2: &java_high2,
        };
        let java_status = Command::new("java")
            .arg("-Xmx1g")
            .arg("-cp")
            .arg(BBTOOLS_CP)
            .arg("jgi.KmerNormalize")
            .args(shared_args_with_extra(java_paths, &[option]))
            .output()
            .expect("run java bbnorm sketch variant");
        assert_success("java bbnorm sketch variant", &java_status);
        assert_ne!(
            base_hist_bytes,
            fs::read(&java_hist).expect("read sketch hist"),
            "{option} should change the vendored Java histogram on the paired phiX fixture"
        );

        let rust_keep1 = dir.path().join(format!("rust.{label}.keep1.fq"));
        let rust_keep2 = dir.path().join(format!("rust.{label}.keep2.fq"));
        let rust_hist = dir.path().join(format!("rust.{label}.hist.tsv"));
        let rust_rhist = dir.path().join(format!("rust.{label}.rhist.tsv"));
        let rust_low1 = dir.path().join(format!("rust.{label}.low1.fq"));
        let rust_low2 = dir.path().join(format!("rust.{label}.low2.fq"));
        let rust_mid1 = dir.path().join(format!("rust.{label}.mid1.fq"));
        let rust_mid2 = dir.path().join(format!("rust.{label}.mid2.fq"));
        let rust_high1 = dir.path().join(format!("rust.{label}.high1.fq"));
        let rust_high2 = dir.path().join(format!("rust.{label}.high2.fq"));
        let rust_paths = OutputPaths {
            keep1: &rust_keep1,
            keep2: &rust_keep2,
            hist: &rust_hist,
            rhist: &rust_rhist,
            low1: &rust_low1,
            low2: &rust_low2,
            mid1: &rust_mid1,
            mid2: &rust_mid2,
            high1: &rust_high1,
            high2: &rust_high2,
        };
        let rust_status = Command::new(env!("CARGO_BIN_EXE_bbnorm-rs"))
            .args(shared_args_with_extra(rust_paths, &[option]))
            .output()
            .expect("run bbnorm-rs sketch fallback");
        assert_success("bbnorm-rs sketch fallback", &rust_status);
        let stderr = String::from_utf8_lossy(&rust_status.stderr);
        if option == "buildpasses=2" {
            assert!(
                stderr.contains("trusted-kmer filtering"),
                "unexpected stderr for {option}: {stderr}"
            );
            assert_ne!(
                base_hist_bytes,
                fs::read(&rust_hist).expect("read rust build-pass hist"),
                "{option} should now exercise real Rust trusted build-pass filtering"
            );
        } else {
            assert!(
                stderr.contains("prefilter collision estimates"),
                "unexpected stderr for {option}: {stderr}"
            );
            assert_ne!(
                base_hist_bytes,
                fs::read(&rust_hist).expect("read rust sketch hist"),
                "{option} should now exercise real Rust prefilter sketch behavior"
            );
        }
        for path in [
            &rust_keep1,
            &rust_keep2,
            &rust_rhist,
            &rust_low1,
            &rust_low2,
            &rust_mid1,
            &rust_mid2,
            &rust_high1,
            &rust_high2,
        ] {
            assert!(path.exists(), "expected output file for {option}: {path:?}");
        }
    }
}

#[test]
fn real_phi_x_pair_matches_java_bbnorm_for_threads_two() {
    assert_normalize_matches_java_with_extra(&["threads=2"]);
}

#[test]
fn real_phi_x_pair_matches_java_bbnorm_for_interleaved_true_with_in2() {
    assert_keepall_hist_matches_java_with_extra(&["interleaved=t"]);
}

#[test]
fn real_phi_x_pair_matches_java_bbnorm_for_canonical_false() {
    assert_keepall_hist_matches_java_with_extra(&["canonical=f"]);
}

#[test]
fn real_phi_x_pair_matches_java_bbnorm_for_keep_duplicate_kmers() {
    assert_keepall_hist_matches_java_with_extra(&["rdk=f"]);
}

#[test]
fn real_phi_x_pair_matches_java_bbnorm_for_canonical_false_keep_duplicate_kmers() {
    assert_keepall_hist_matches_java_with_extra(&["canonical=f", "rdk=f"]);
}

#[test]
fn real_phi_x_pair_matches_java_bbnorm_for_max_depth_clamped_to_target() {
    assert_keepall_hist_matches_java_with_extra(&["target=100", "max=50"]);
}

#[test]
fn real_phi_x_pair_matches_java_bbnorm_for_zero_minkmers_clamped_to_one() {
    assert_keepall_hist_matches_java_with_extra(&["minkmers=0"]);
}

#[test]
fn real_phi_x_pair_matches_java_bbnorm_for_fixspikes() {
    assert_keepall_hist_matches_java_with_extra(&["fixspikes=t"]);
}

#[test]
fn real_phi_x_pair_matches_java_bbnorm_for_inactive_trim_parser_options() {
    assert_keepall_hist_matches_java_with_extra(&[
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
}

#[test]
fn real_phi_x_pair_matches_java_bbnorm_for_inactive_read_filter_parser_options() {
    assert_keepall_hist_matches_java_with_extra(&[
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
}

#[test]
fn real_phi_x_pair_matches_java_bbnorm_for_hash_output_patterns() {
    require_file(SAMPLE1);
    require_file(SAMPLE2);
    require_file("vendor/BBTools-master/current/jgi/KmerNormalize.class");
    require_java();

    let dir = tempdir().expect("create tempdir");
    let java_keep_hash = dir.path().join("java.keep#.fq");
    let java_keep1 = dir.path().join("java.keep1.fq");
    let java_keep2 = dir.path().join("java.keep2.fq");
    let java_low_hash = dir.path().join("java.low#.fq");
    let java_low1 = dir.path().join("java.low1.fq");
    let java_low2 = dir.path().join("java.low2.fq");
    let rust_keep_hash = dir.path().join("rust.keep#.fq");
    let rust_keep1 = dir.path().join("rust.keep1.fq");
    let rust_keep2 = dir.path().join("rust.keep2.fq");
    let rust_low_hash = dir.path().join("rust.low#.fq");
    let rust_low1 = dir.path().join("rust.low1.fq");
    let rust_low2 = dir.path().join("rust.low2.fq");

    let shared = hash_output_pattern_args(&java_keep_hash, &java_low_hash);
    let java_status = Command::new("java")
        .arg("-Xmx1g")
        .arg("-cp")
        .arg(BBTOOLS_CP)
        .arg("jgi.KmerNormalize")
        .args(&shared)
        .output()
        .expect("run java bbnorm");
    assert_success("java bbnorm", &java_status);

    let rust_args = hash_output_pattern_args(&rust_keep_hash, &rust_low_hash);
    let rust_status = Command::new(env!("CARGO_BIN_EXE_bbnorm-rs"))
        .args(&rust_args)
        .output()
        .expect("run bbnorm-rs");
    assert_success("bbnorm-rs", &rust_status);

    assert_same_file(&java_keep1, &rust_keep1);
    assert_same_file(&java_keep2, &rust_keep2);
    assert_same_file(&java_low1, &rust_low1);
    assert_same_file(&java_low2, &rust_low2);
}

#[test]
fn real_phi_x_pair_matches_java_bbnorm_for_hash_input_pattern() {
    require_file(SAMPLE1);
    require_file(SAMPLE2);
    require_file("vendor/BBTools-master/current/jgi/KmerNormalize.class");
    require_java();

    let dir = tempdir().expect("create tempdir");
    let java_keep1 = dir.path().join("java.keep1.fq");
    let java_keep2 = dir.path().join("java.keep2.fq");
    let java_toss1 = dir.path().join("java.toss1.fq");
    let java_toss2 = dir.path().join("java.toss2.fq");
    let rust_keep1 = dir.path().join("rust.keep1.fq");
    let rust_keep2 = dir.path().join("rust.keep2.fq");
    let rust_toss1 = dir.path().join("rust.toss1.fq");
    let rust_toss2 = dir.path().join("rust.toss2.fq");

    let java_paths = NormalizePaths {
        keep1: &java_keep1,
        keep2: &java_keep2,
        toss1: &java_toss1,
        toss2: &java_toss2,
    };
    let shared = normalize_hash_input_args(java_paths);
    let java_status = Command::new("java")
        .arg("-Xmx1g")
        .arg("-cp")
        .arg(BBTOOLS_CP)
        .arg("jgi.KmerNormalize")
        .args(&shared)
        .output()
        .expect("run java bbnorm");
    assert_success("java bbnorm", &java_status);

    let rust_paths = NormalizePaths {
        keep1: &rust_keep1,
        keep2: &rust_keep2,
        toss1: &rust_toss1,
        toss2: &rust_toss2,
    };
    let rust_args = normalize_hash_input_args(rust_paths);
    let rust_status = Command::new(env!("CARGO_BIN_EXE_bbnorm-rs"))
        .args(&rust_args)
        .output()
        .expect("run bbnorm-rs");
    assert_success("bbnorm-rs", &rust_status);

    assert_same_file(&java_keep1, &rust_keep1);
    assert_same_file(&java_keep2, &rust_keep2);
    assert_same_file(&java_toss1, &rust_toss1);
    assert_same_file(&java_toss2, &rust_toss2);
}

#[test]
fn real_phi_x_pair_matches_java_bbnorm_for_implicit_interleaved_output() {
    require_file(SAMPLE1);
    require_file(SAMPLE2);
    require_file("vendor/BBTools-master/current/jgi/KmerNormalize.class");
    require_java();

    let dir = tempdir().expect("create tempdir");
    let java_keep = dir.path().join("java.keep.fq");
    let rust_keep = dir.path().join("rust.keep.fq");

    let shared = implicit_interleaved_paired_output_args(&java_keep);
    let java_status = Command::new("java")
        .arg("-Xmx1g")
        .arg("-cp")
        .arg(BBTOOLS_CP)
        .arg("jgi.KmerNormalize")
        .args(&shared)
        .output()
        .expect("run java bbnorm");
    assert_success("java bbnorm", &java_status);

    let rust_args = implicit_interleaved_paired_output_args(&rust_keep);
    let rust_status = Command::new(env!("CARGO_BIN_EXE_bbnorm-rs"))
        .args(&rust_args)
        .output()
        .expect("run bbnorm-rs");
    assert_success("bbnorm-rs", &rust_status);

    assert_same_file(&java_keep, &rust_keep);
}

#[test]
fn real_phi_x_pair_matches_java_bbnorm_for_outuncorrected_without_ecc() {
    require_file(SAMPLE1);
    require_file(SAMPLE2);
    require_file("vendor/BBTools-master/current/jgi/KmerNormalize.class");
    require_java();

    let dir = tempdir().expect("create tempdir");
    let java_keep1 = dir.path().join("java.keep1.fq");
    let java_keep2 = dir.path().join("java.keep2.fq");
    let java_unc1 = dir.path().join("java.unc1.fq");
    let java_unc2 = dir.path().join("java.unc2.fq");
    let rust_keep1 = dir.path().join("rust.keep1.fq");
    let rust_keep2 = dir.path().join("rust.keep2.fq");
    let rust_unc1 = dir.path().join("rust.unc1.fq");
    let rust_unc2 = dir.path().join("rust.unc2.fq");

    let shared = outuncorrected_noecc_args(&java_keep1, &java_keep2, &java_unc1, &java_unc2);
    let java_status = Command::new("java")
        .arg("-Xmx1g")
        .arg("-cp")
        .arg(BBTOOLS_CP)
        .arg("jgi.KmerNormalize")
        .args(&shared)
        .output()
        .expect("run java bbnorm");
    assert_success("java bbnorm", &java_status);

    let rust_args = outuncorrected_noecc_args(&rust_keep1, &rust_keep2, &rust_unc1, &rust_unc2);
    let rust_status = Command::new(env!("CARGO_BIN_EXE_bbnorm-rs"))
        .args(&rust_args)
        .output()
        .expect("run bbnorm-rs");
    assert_success("bbnorm-rs", &rust_status);

    assert_same_file(&java_keep1, &rust_keep1);
    assert_same_file(&java_keep2, &rust_keep2);
    assert_same_file(&java_unc1, &rust_unc1);
    assert_same_file(&java_unc2, &rust_unc2);
    assert_eq!(fs::metadata(&rust_unc1).expect("unc1 metadata").len(), 0);
    assert_eq!(fs::metadata(&rust_unc2).expect("unc2 metadata").len(), 0);
}

#[test]
fn representative_overlap_only_ecco_repetitive_fixture_matches_java_bbnorm() {
    require_file("vendor/BBTools-master/current/jgi/KmerNormalize.class");
    require_java();

    let dir = tempdir().expect("create tempdir");
    let input1 = dir.path().join("overlap.r1.fq");
    let input2 = dir.path().join("overlap.r2.fq");
    write_overlap_only_ecco_fixture(&input1, &input2);

    let java_keep1 = dir.path().join("java.keep1.fq");
    let java_keep2 = dir.path().join("java.keep2.fq");
    let rust_keep1 = dir.path().join("rust.keep1.fq");
    let rust_keep2 = dir.path().join("rust.keep2.fq");

    let java_status = Command::new("java")
        .arg("-Xmx1g")
        .arg("-cp")
        .arg(BBTOOLS_CP)
        .arg("jgi.KmerNormalize")
        .args(overlap_only_ecco_args(
            &input1,
            &input2,
            &java_keep1,
            &java_keep2,
            true,
        ))
        .output()
        .expect("run java bbnorm");
    assert_success("java bbnorm", &java_status);

    let rust_status = Command::new(env!("CARGO_BIN_EXE_bbnorm-rs"))
        .args(overlap_only_ecco_args(
            &input1,
            &input2,
            &rust_keep1,
            &rust_keep2,
            true,
        ))
        .output()
        .expect("run bbnorm-rs");
    assert_success("bbnorm-rs", &rust_status);

    assert_same_file(&java_keep1, &rust_keep1);
    assert_same_file(&java_keep2, &rust_keep2);
}

#[test]
fn representative_overlap_only_ecco_high_entropy_sequence_matches_java_bbnorm() {
    require_file("vendor/BBTools-master/current/jgi/KmerNormalize.class");
    require_java();

    let dir = tempdir().expect("create tempdir");
    let input1 = dir.path().join("overlap_hi.r1.fq");
    let input2 = dir.path().join("overlap_hi.r2.fq");
    write_overlap_only_ecco_high_entropy_fixture(&input1, &input2);

    let java_keep1 = dir.path().join("java.keep1.fq");
    let java_keep2 = dir.path().join("java.keep2.fq");
    let rust_keep1 = dir.path().join("rust.keep1.fq");
    let rust_keep2 = dir.path().join("rust.keep2.fq");

    for ecco in [false, true] {
        let java_status = Command::new("java")
            .arg("-Xmx1g")
            .arg("-cp")
            .arg(BBTOOLS_CP)
            .arg("jgi.KmerNormalize")
            .args(overlap_only_ecco_args(
                &input1,
                &input2,
                &java_keep1,
                &java_keep2,
                ecco,
            ))
            .output()
            .expect("run java bbnorm");
        assert_success("java bbnorm", &java_status);

        let rust_status = Command::new(env!("CARGO_BIN_EXE_bbnorm-rs"))
            .args(overlap_only_ecco_args(
                &input1,
                &input2,
                &rust_keep1,
                &rust_keep2,
                ecco,
            ))
            .output()
            .expect("run bbnorm-rs");
        assert_success("bbnorm-rs", &rust_status);

        let java_mate1 = fs::read_to_string(&java_keep1)
            .expect("read java overlap output")
            .lines()
            .nth(1)
            .expect("java overlap mate-1 sequence")
            .to_string();
        let rust_mate1 = fs::read_to_string(&rust_keep1)
            .expect("read rust overlap output")
            .lines()
            .nth(1)
            .expect("rust overlap mate-1 sequence")
            .to_string();
        let java_mate2 = fs::read_to_string(&java_keep2)
            .expect("read java overlap output")
            .lines()
            .nth(1)
            .expect("java overlap mate-2 sequence")
            .to_string();
        let rust_mate2 = fs::read_to_string(&rust_keep2)
            .expect("read rust overlap output")
            .lines()
            .nth(1)
            .expect("rust overlap mate-2 sequence")
            .to_string();
        assert_eq!(java_mate1, "TTAGTTGTGCCGCAGCGAAGTAGTGCTTGAAATATGCGAC");
        assert_eq!(rust_mate1, "TTAGTTGTGCCGCAGCGAAGTAGTGCTTGAAATATGCGAC");
        if ecco {
            assert_eq!(java_mate2, "GTCGCATATTTCAAGCACTACTTCGCTGCGGCACAACTAA");
            assert_eq!(rust_mate2, "GTCGCATATTTCAAGCACTACTTCGCTGCGGCACAACTAA");
        } else {
            assert_eq!(java_mate2, "GTCGCATATTTCAAGCACTAATTCGCTGCGGCACAACTAA");
            assert_eq!(rust_mate2, "GTCGCATATTTCAAGCACTAATTCGCTGCGGCACAACTAA");
        }
    }
}

#[test]
fn representative_overlap_only_ecco_auto_matches_java_bbnorm() {
    require_file("vendor/BBTools-master/current/jgi/KmerNormalize.class");
    require_java();

    let dir = tempdir().expect("create tempdir");
    let input1 = dir.path().join("overlap_auto.r1.fq");
    let input2 = dir.path().join("overlap_auto.r2.fq");
    write_overlap_only_ecco_high_entropy_fixture(&input1, &input2);

    let java_keep1 = dir.path().join("java.keep1.fq");
    let java_keep2 = dir.path().join("java.keep2.fq");
    let rust_keep1 = dir.path().join("rust.keep1.fq");
    let rust_keep2 = dir.path().join("rust.keep2.fq");

    let java_status = Command::new("java")
        .arg("-Xmx1g")
        .arg("-cp")
        .arg(BBTOOLS_CP)
        .arg("jgi.KmerNormalize")
        .args(overlap_only_ecco_value_args(
            &input1,
            &input2,
            &java_keep1,
            &java_keep2,
            "auto",
        ))
        .output()
        .expect("run java bbnorm");
    assert_success("java bbnorm", &java_status);

    let rust_status = Command::new(env!("CARGO_BIN_EXE_bbnorm-rs"))
        .args(overlap_only_ecco_value_args(
            &input1,
            &input2,
            &rust_keep1,
            &rust_keep2,
            "auto",
        ))
        .output()
        .expect("run bbnorm-rs");
    assert_success("bbnorm-rs", &rust_status);

    assert_same_file(&java_keep1, &rust_keep1);
    assert_same_file(&java_keep2, &rust_keep2);
}

#[test]
fn representative_overlap_only_ecco_high_entropy_quality_matches_java_bbnorm() {
    require_file("vendor/BBTools-master/current/jgi/KmerNormalize.class");
    require_java();

    let dir = tempdir().expect("create tempdir");
    let input1 = dir.path().join("overlap_hi.r1.fq");
    let input2 = dir.path().join("overlap_hi.r2.fq");
    write_overlap_only_ecco_high_entropy_fixture(&input1, &input2);

    let java_keep1 = dir.path().join("java.keep1.fq");
    let java_keep2 = dir.path().join("java.keep2.fq");
    let rust_keep1 = dir.path().join("rust.keep1.fq");
    let rust_keep2 = dir.path().join("rust.keep2.fq");

    let java_status = Command::new("java")
        .arg("-Xmx1g")
        .arg("-cp")
        .arg(BBTOOLS_CP)
        .arg("jgi.KmerNormalize")
        .args(overlap_only_ecco_args(
            &input1,
            &input2,
            &java_keep1,
            &java_keep2,
            true,
        ))
        .output()
        .expect("run java bbnorm");
    assert_success("java bbnorm", &java_status);

    let rust_status = Command::new(env!("CARGO_BIN_EXE_bbnorm-rs"))
        .args(overlap_only_ecco_args(
            &input1,
            &input2,
            &rust_keep1,
            &rust_keep2,
            true,
        ))
        .output()
        .expect("run bbnorm-rs");
    assert_success("bbnorm-rs", &rust_status);

    let java_quality = fs::read_to_string(&java_keep2)
        .expect("read java overlap output")
        .lines()
        .nth(3)
        .expect("java overlap mate-2 quality")
        .to_string();
    let rust_quality = fs::read_to_string(&rust_keep2)
        .expect("read rust overlap output")
        .lines()
        .nth(3)
        .expect("rust overlap mate-2 quality")
        .to_string();
    assert_eq!(java_quality, "SSSSSSSSSSSSSSSSSSSSGSSSSSSSSSSSSSSSSSSS");
    assert_eq!(rust_quality, java_quality);
}

#[test]
fn representative_overlap_only_ecco_high_confidence_mismatch_matches_java_bbnorm() {
    require_file("vendor/BBTools-master/current/jgi/KmerNormalize.class");
    require_java();

    let dir = tempdir().expect("create tempdir");
    let input1 = dir.path().join("overlap_hi_conf.r1.fq");
    let input2 = dir.path().join("overlap_hi_conf.r2.fq");
    write_overlap_only_ecco_high_confidence_fixture(&input1, &input2);

    let java_keep1 = dir.path().join("java.keep1.fq");
    let java_keep2 = dir.path().join("java.keep2.fq");
    let rust_keep1 = dir.path().join("rust.keep1.fq");
    let rust_keep2 = dir.path().join("rust.keep2.fq");

    let java_status = Command::new("java")
        .arg("-Xmx1g")
        .arg("-cp")
        .arg(BBTOOLS_CP)
        .arg("jgi.KmerNormalize")
        .args(overlap_only_ecco_args(
            &input1,
            &input2,
            &java_keep1,
            &java_keep2,
            true,
        ))
        .output()
        .expect("run java bbnorm");
    assert_success("java bbnorm", &java_status);

    let rust_status = Command::new(env!("CARGO_BIN_EXE_bbnorm-rs"))
        .args(overlap_only_ecco_args(
            &input1,
            &input2,
            &rust_keep1,
            &rust_keep2,
            true,
        ))
        .output()
        .expect("run bbnorm-rs");
    assert_success("bbnorm-rs", &rust_status);

    assert_same_file(&java_keep1, &rust_keep1);
    assert_same_file(&java_keep2, &rust_keep2);
}

#[test]
fn representative_overlap_only_ecco_quality_threshold_boundary_matches_java_bbnorm() {
    require_file("vendor/BBTools-master/current/jgi/KmerNormalize.class");
    require_java();

    let dir = tempdir().expect("create tempdir");

    for (label, mismatch_quality, expect_corrected) in
        [("repair_q7", 7u8, true), ("reject_q8", 8u8, false)]
    {
        let input1 = dir.path().join(format!("{label}.r1.fq"));
        let input2 = dir.path().join(format!("{label}.r2.fq"));
        write_overlap_only_ecco_quality_fixture(&input1, &input2, mismatch_quality);

        let java_keep1 = dir.path().join(format!("{label}.java.keep1.fq"));
        let java_keep2 = dir.path().join(format!("{label}.java.keep2.fq"));
        let rust_keep1 = dir.path().join(format!("{label}.rust.keep1.fq"));
        let rust_keep2 = dir.path().join(format!("{label}.rust.keep2.fq"));

        let java_status = Command::new("java")
            .arg("-Xmx1g")
            .arg("-cp")
            .arg(BBTOOLS_CP)
            .arg("jgi.KmerNormalize")
            .args(overlap_only_ecco_args(
                &input1,
                &input2,
                &java_keep1,
                &java_keep2,
                true,
            ))
            .output()
            .expect("run java bbnorm");
        assert_success("java bbnorm", &java_status);

        let rust_status = Command::new(env!("CARGO_BIN_EXE_bbnorm-rs"))
            .args(overlap_only_ecco_args(
                &input1,
                &input2,
                &rust_keep1,
                &rust_keep2,
                true,
            ))
            .output()
            .expect("run bbnorm-rs");
        assert_success("bbnorm-rs", &rust_status);

        assert_same_file(&java_keep1, &rust_keep1);
        assert_same_file(&java_keep2, &rust_keep2);

        let mate2 = fs::read_to_string(&rust_keep2)
            .expect("read rust overlap output")
            .lines()
            .nth(1)
            .expect("rust overlap mate-2 sequence")
            .to_string();
        let corrected = mate2 == "GTCGCATATTTCAAGCACTACTTCGCTGCGGCACAACTAA";
        assert_eq!(corrected, expect_corrected, "{label} corrected state");
    }
}

#[test]
fn representative_overlap_only_ecco_competing_short_overlap_matches_java_bbnorm() {
    require_file("vendor/BBTools-master/current/jgi/KmerNormalize.class");
    require_java();

    let dir = tempdir().expect("create tempdir");
    let input1 = dir.path().join("overlap_weighted_multi.r1.fq");
    let input2 = dir.path().join("overlap_weighted_multi.r2.fq");
    write_overlap_only_ecco_quality_weighted_multimismatch_fixture(&input1, &input2);

    let java_keep1 = dir.path().join("java.keep1.fq");
    let java_keep2 = dir.path().join("java.keep2.fq");
    let rust_keep1 = dir.path().join("rust.keep1.fq");
    let rust_keep2 = dir.path().join("rust.keep2.fq");

    let java_status = Command::new("java")
        .arg("-Xmx1g")
        .arg("-cp")
        .arg(BBTOOLS_CP)
        .arg("jgi.KmerNormalize")
        .args(overlap_only_ecco_args(
            &input1,
            &input2,
            &java_keep1,
            &java_keep2,
            true,
        ))
        .output()
        .expect("run java bbnorm");
    assert_success("java bbnorm", &java_status);

    let rust_status = Command::new(env!("CARGO_BIN_EXE_bbnorm-rs"))
        .args(overlap_only_ecco_args(
            &input1,
            &input2,
            &rust_keep1,
            &rust_keep2,
            true,
        ))
        .output()
        .expect("run bbnorm-rs");
    assert_success("bbnorm-rs", &rust_status);

    assert_same_file(&java_keep1, &rust_keep1);
    assert_same_file(&java_keep2, &rust_keep2);
}

#[test]
fn representative_multipass_first_stage_markuncorrectable_matches_java_outuncorrected_routing() {
    assert_multipass_markuncorrectable_matches_java_with_stage(&["ecc1=t", "eccf=f"], 0);
}

#[test]
fn representative_multipass_final_stage_markuncorrectable_matches_java_outuncorrected_routing() {
    assert_multipass_markuncorrectable_matches_java_with_stage(&["ecc1=f", "eccf=t"], 1);
}

#[test]
fn representative_multipass_both_stage_markuncorrectable_matches_java_outuncorrected_routing() {
    assert_multipass_markuncorrectable_matches_java_with_stage(&["ecc=t"], 1);
}

fn assert_multipass_markuncorrectable_matches_java_with_stage(
    stage_args: &[&str],
    expected_uncorrected_records: usize,
) {
    require_file("vendor/BBTools-master/current/jgi/KmerNormalize.class");
    require_java();

    let dir = tempdir().expect("create tempdir");
    let input1 = dir.path().join("multipass_uncorrectable.r1.fq");
    let input2 = dir.path().join("multipass_uncorrectable.r2.fq");
    write_multipass_uncorrectable_paired_fixture(&input1, &input2);

    let java_keep1 = dir.path().join("java.keep1.fq");
    let java_keep2 = dir.path().join("java.keep2.fq");
    let java_unc1 = dir.path().join("java.unc1.fq");
    let java_unc2 = dir.path().join("java.unc2.fq");
    let rust_keep1 = dir.path().join("rust.keep1.fq");
    let rust_keep2 = dir.path().join("rust.keep2.fq");
    let rust_unc1 = dir.path().join("rust.unc1.fq");
    let rust_unc2 = dir.path().join("rust.unc2.fq");

    let java_args = multipass_markuncorrectable_args(
        &input1,
        &input2,
        &java_keep1,
        &java_keep2,
        &java_unc1,
        &java_unc2,
        stage_args,
    );
    let java_status = Command::new("java")
        .arg("-Xmx1g")
        .arg("-cp")
        .arg(BBTOOLS_CP)
        .arg("jgi.KmerNormalize")
        .args(&java_args)
        .output()
        .expect("run java bbnorm");
    assert_success("java bbnorm", &java_status);

    let rust_args = multipass_markuncorrectable_args(
        &input1,
        &input2,
        &rust_keep1,
        &rust_keep2,
        &rust_unc1,
        &rust_unc2,
        stage_args,
    );
    let rust_status = Command::new(env!("CARGO_BIN_EXE_bbnorm-rs"))
        .args(&rust_args)
        .output()
        .expect("run bbnorm-rs");
    assert_success("bbnorm-rs", &rust_status);

    assert_same_file(&java_keep1, &rust_keep1);
    assert_same_file(&java_keep2, &rust_keep2);
    assert_same_file(&java_unc1, &rust_unc1);
    assert_same_file(&java_unc2, &rust_unc2);
    assert_eq!(
        fastq_record_count(&java_unc1),
        expected_uncorrected_records,
        "unexpected java mate-1 outuncorrected record count for staged multipass ECC"
    );
    assert_eq!(
        fastq_record_count(&java_unc2),
        expected_uncorrected_records,
        "unexpected java mate-2 outuncorrected record count for staged multipass ECC"
    );
    assert_eq!(
        fastq_record_count(&rust_unc1),
        expected_uncorrected_records,
        "unexpected rust mate-1 outuncorrected record count for staged multipass ECC"
    );
    assert_eq!(
        fastq_record_count(&rust_unc2),
        expected_uncorrected_records,
        "unexpected rust mate-2 outuncorrected record count for staged multipass ECC"
    );
}

fn assert_keepall_hist_matches_java_with_extra(extra: &[&str]) {
    require_file(SAMPLE1);
    require_file(SAMPLE2);
    require_file("vendor/BBTools-master/current/jgi/KmerNormalize.class");
    require_java();

    let dir = tempdir().expect("create tempdir");
    let java_keep1 = dir.path().join("java.keep1.fq");
    let java_keep2 = dir.path().join("java.keep2.fq");
    let java_hist = dir.path().join("java.hist.tsv");
    let java_rhist = dir.path().join("java.rhist.tsv");
    let java_low1 = dir.path().join("java.low1.fq");
    let java_low2 = dir.path().join("java.low2.fq");
    let java_mid1 = dir.path().join("java.mid1.fq");
    let java_mid2 = dir.path().join("java.mid2.fq");
    let java_high1 = dir.path().join("java.high1.fq");
    let java_high2 = dir.path().join("java.high2.fq");
    let rust_keep1 = dir.path().join("rust.keep1.fq");
    let rust_keep2 = dir.path().join("rust.keep2.fq");
    let rust_hist = dir.path().join("rust.hist.tsv");
    let rust_rhist = dir.path().join("rust.rhist.tsv");
    let rust_low1 = dir.path().join("rust.low1.fq");
    let rust_low2 = dir.path().join("rust.low2.fq");
    let rust_mid1 = dir.path().join("rust.mid1.fq");
    let rust_mid2 = dir.path().join("rust.mid2.fq");
    let rust_high1 = dir.path().join("rust.high1.fq");
    let rust_high2 = dir.path().join("rust.high2.fq");

    let java_paths = OutputPaths {
        keep1: &java_keep1,
        keep2: &java_keep2,
        hist: &java_hist,
        rhist: &java_rhist,
        low1: &java_low1,
        low2: &java_low2,
        mid1: &java_mid1,
        mid2: &java_mid2,
        high1: &java_high1,
        high2: &java_high2,
    };
    let shared = shared_args_with_extra(java_paths, extra);
    let java_status = Command::new("java")
        .arg("-Xmx1g")
        .arg("-cp")
        .arg(BBTOOLS_CP)
        .arg("jgi.KmerNormalize")
        .args(&shared)
        .output()
        .expect("run java bbnorm");
    assert_success("java bbnorm", &java_status);

    let rust_paths = OutputPaths {
        keep1: &rust_keep1,
        keep2: &rust_keep2,
        hist: &rust_hist,
        rhist: &rust_rhist,
        low1: &rust_low1,
        low2: &rust_low2,
        mid1: &rust_mid1,
        mid2: &rust_mid2,
        high1: &rust_high1,
        high2: &rust_high2,
    };
    let rust_args = shared_args_with_extra(rust_paths, extra);
    let rust_status = Command::new(env!("CARGO_BIN_EXE_bbnorm-rs"))
        .args(&rust_args)
        .output()
        .expect("run bbnorm-rs");
    assert_success("bbnorm-rs", &rust_status);

    assert_same_file(&java_keep1, &rust_keep1);
    assert_same_file(&java_keep2, &rust_keep2);
    assert_same_file(&java_hist, &rust_hist);
    assert_same_file(&java_rhist, &rust_rhist);
    assert_same_file(&java_low1, &rust_low1);
    assert_same_file(&java_low2, &rust_low2);
    assert_same_file(&java_mid1, &rust_mid1);
    assert_same_file(&java_mid2, &rust_mid2);
    assert_same_file(&java_high1, &rust_high1);
    assert_same_file(&java_high2, &rust_high2);
}

#[test]
fn real_phi_x_pair_matches_java_bbnorm_for_renamed_keepall_output() {
    require_file(SAMPLE1);
    require_file(SAMPLE2);
    require_file("vendor/BBTools-master/current/jgi/KmerNormalize.class");
    require_java();

    let dir = tempdir().expect("create tempdir");
    let java_keep1 = dir.path().join("java.keep1.fq");
    let java_keep2 = dir.path().join("java.keep2.fq");
    let rust_keep1 = dir.path().join("rust.keep1.fq");
    let rust_keep2 = dir.path().join("rust.keep2.fq");

    let shared = rename_paired_keepall_args(&java_keep1, &java_keep2);
    let java_status = Command::new("java")
        .arg("-Xmx1g")
        .arg("-cp")
        .arg(BBTOOLS_CP)
        .arg("jgi.KmerNormalize")
        .args(&shared)
        .output()
        .expect("run java bbnorm");
    assert_success("java bbnorm", &java_status);

    let rust_args = rename_paired_keepall_args(&rust_keep1, &rust_keep2);
    let rust_status = Command::new(env!("CARGO_BIN_EXE_bbnorm-rs"))
        .args(&rust_args)
        .output()
        .expect("run bbnorm-rs");
    assert_success("bbnorm-rs", &rust_status);

    assert_same_file(&java_keep1, &rust_keep1);
    assert_same_file(&java_keep2, &rust_keep2);
}

#[test]
fn real_phi_x_pair_matches_java_bbnorm_for_single_pass_error_toss_normalization() {
    require_file(SAMPLE1);
    require_file(SAMPLE2);
    require_file("vendor/BBTools-master/current/jgi/KmerNormalize.class");
    require_java();

    let dir = tempdir().expect("create tempdir");
    let java_keep1 = dir.path().join("java.keep1.fq");
    let java_keep2 = dir.path().join("java.keep2.fq");
    let java_toss1 = dir.path().join("java.toss1.fq");
    let java_toss2 = dir.path().join("java.toss2.fq");
    let rust_keep1 = dir.path().join("rust.keep1.fq");
    let rust_keep2 = dir.path().join("rust.keep2.fq");
    let rust_toss1 = dir.path().join("rust.toss1.fq");
    let rust_toss2 = dir.path().join("rust.toss2.fq");

    let java_paths = NormalizePaths {
        keep1: &java_keep1,
        keep2: &java_keep2,
        toss1: &java_toss1,
        toss2: &java_toss2,
    };
    let shared = normalize_args(java_paths);
    let java_status = Command::new("java")
        .arg("-Xmx1g")
        .arg("-cp")
        .arg(BBTOOLS_CP)
        .arg("jgi.KmerNormalize")
        .args(&shared)
        .output()
        .expect("run java bbnorm");
    assert_success("java bbnorm", &java_status);

    let rust_paths = NormalizePaths {
        keep1: &rust_keep1,
        keep2: &rust_keep2,
        toss1: &rust_toss1,
        toss2: &rust_toss2,
    };
    let rust_args = normalize_args(rust_paths);
    let rust_status = Command::new(env!("CARGO_BIN_EXE_bbnorm-rs"))
        .args(&rust_args)
        .output()
        .expect("run bbnorm-rs");
    assert_success("bbnorm-rs", &rust_status);

    assert_same_file(&java_keep1, &rust_keep1);
    assert_same_file(&java_keep2, &rust_keep2);
    assert_same_file(&java_toss1, &rust_toss1);
    assert_same_file(&java_toss2, &rust_toss2);
}

#[test]
fn real_phi_x_pair_matches_java_bbnorm_for_lowthresh_zero_error_toss() {
    assert_normalize_matches_java_with_extra(&["lowthresh=0"]);
}

#[test]
fn real_phi_x_pair_matches_java_bbnorm_for_highthresh_one_error_toss() {
    assert_normalize_matches_java_with_extra(&["highthresh=1"]);
}

#[test]
fn real_phi_x_pair_matches_java_bbnorm_for_errordetectratio_two_error_toss() {
    assert_normalize_matches_java_with_extra(&["errordetectratio=2"]);
}

#[test]
fn real_phi_x_pair_matches_java_bbnorm_for_mixed_error_thresholds() {
    assert_normalize_matches_java_with_extra(&[
        "lowthresh=1",
        "highthresh=1",
        "errordetectratio=2",
    ]);
}

#[test]
fn real_phi_x_pair_matches_java_bbnorm_for_output_histograms_after_error_toss() {
    require_file(SAMPLE1);
    require_file(SAMPLE2);
    require_file("vendor/BBTools-master/current/jgi/KmerNormalize.class");
    require_java();

    let dir = tempdir().expect("create tempdir");
    let java_keep1 = dir.path().join("java.keep1.fq");
    let java_keep2 = dir.path().join("java.keep2.fq");
    let java_toss1 = dir.path().join("java.toss1.fq");
    let java_toss2 = dir.path().join("java.toss2.fq");
    let java_hist = dir.path().join("java.hist.tsv");
    let java_rhist = dir.path().join("java.rhist.tsv");
    let java_histout = dir.path().join("java.histout.tsv");
    let java_rhistout = dir.path().join("java.rhistout.tsv");
    let rust_keep1 = dir.path().join("rust.keep1.fq");
    let rust_keep2 = dir.path().join("rust.keep2.fq");
    let rust_toss1 = dir.path().join("rust.toss1.fq");
    let rust_toss2 = dir.path().join("rust.toss2.fq");
    let rust_hist = dir.path().join("rust.hist.tsv");
    let rust_rhist = dir.path().join("rust.rhist.tsv");
    let rust_histout = dir.path().join("rust.histout.tsv");
    let rust_rhistout = dir.path().join("rust.rhistout.tsv");

    let shared = normalize_with_histograms_args(
        NormalizePaths {
            keep1: &java_keep1,
            keep2: &java_keep2,
            toss1: &java_toss1,
            toss2: &java_toss2,
        },
        &java_hist,
        &java_rhist,
        &java_histout,
        &java_rhistout,
        &["highthresh=1"],
    );
    let java_status = Command::new("java")
        .arg("-Xmx1g")
        .arg("-cp")
        .arg(BBTOOLS_CP)
        .arg("jgi.KmerNormalize")
        .args(&shared)
        .output()
        .expect("run java bbnorm");
    assert_success("java bbnorm", &java_status);

    let rust_args = normalize_with_histograms_args(
        NormalizePaths {
            keep1: &rust_keep1,
            keep2: &rust_keep2,
            toss1: &rust_toss1,
            toss2: &rust_toss2,
        },
        &rust_hist,
        &rust_rhist,
        &rust_histout,
        &rust_rhistout,
        &["highthresh=1"],
    );
    let rust_status = Command::new(env!("CARGO_BIN_EXE_bbnorm-rs"))
        .args(&rust_args)
        .output()
        .expect("run bbnorm-rs");
    assert_success("bbnorm-rs", &rust_status);

    assert_same_file(&java_keep1, &rust_keep1);
    assert_same_file(&java_keep2, &rust_keep2);
    assert_same_file(&java_toss1, &rust_toss1);
    assert_same_file(&java_toss2, &rust_toss2);
    assert_same_file(&java_hist, &rust_hist);
    assert_same_file(&java_rhist, &rust_rhist);
    assert_same_file(&java_histout, &rust_histout);
    assert_same_file(&java_rhistout, &rust_rhistout);
}

#[test]
fn real_phi_x_pair_matches_java_bbnorm_for_requirebothbad_error_toss() {
    assert_normalize_matches_java_with_extra(&["requirebothbad=t"]);
}

#[test]
fn real_phi_x_pair_matches_java_bbnorm_for_removeifeitherbad_alias() {
    assert_normalize_matches_java_with_extra(&["requirebothbad=t", "removeifeitherbad=t"]);
}

#[test]
fn real_phi_x_pair_matches_java_bbnorm_for_save_rare_reads_error_toss() {
    assert_normalize_matches_java_with_extra(&["saverarereads=t", "highthresh=1"]);
}

#[test]
fn real_phi_x_pair_matches_java_bbnorm_for_save_rare_requireboth_modes() {
    let cases: [&[&str]; 2] = [
        &["saverarereads=t", "highthresh=1", "requirebothbad=t"],
        &[
            "saverarereads=t",
            "highthresh=1",
            "requirebothbad=t",
            "removeifeitherbad=t",
        ],
    ];

    for extra in cases {
        assert_normalize_matches_java_with_extra(extra);
    }
}

#[test]
fn real_phi_x_pair_matches_java_bbnorm_for_final_stage_toss_error_alias() {
    assert_normalize_without_default_toss_matches_java_with_extra(&["tossbadreadsf=t"]);
}

#[test]
fn real_phi_x_pair_matches_java_bbnorm_for_final_stage_toss_error_alias_variants() {
    for extra in [&["tossbadreads2=t"][..], &["ter2=t"][..], &["tbr2=t"][..]] {
        assert_normalize_without_default_toss_matches_java_with_extra(extra);
    }
}

#[test]
fn real_phi_x_pair_matches_java_bbnorm_for_first_stage_toss_error_alias_noop() {
    assert_normalize_without_default_toss_matches_java_with_extra(&["tossbadreads1=t"]);
}

#[test]
fn real_phi_x_pair_matches_java_bbnorm_for_single_pass_multipass_only_controls() {
    assert_keepall_hist_matches_java_with_extra(&[
        "target1=7",
        "targetbadpercentilelow=20",
        "targetbadpercentilehigh=80",
        "abrc=t",
        "discardbadonly1=t",
    ]);
}

#[test]
fn real_phi_x_pair_matches_java_bbnorm_for_inactive_ecc_tuning_controls() {
    assert_keepall_hist_matches_java_with_extra(&[
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
}

#[test]
fn real_phi_x_duplicated_high_depth_matches_java_bbnorm_for_discardbadonly() {
    require_file(SAMPLE1);
    require_file("vendor/BBTools-master/current/jgi/KmerNormalize.class");
    require_java();

    let dir = tempdir().expect("create tempdir");
    let input = dir.path().join("real_duplicated_high_depth.fq");
    write_duplicated_n_free_single_end_fixture(SAMPLE1, &input, 8);

    let java_keep = dir.path().join("java.keep.fq");
    let java_toss = dir.path().join("java.toss.fq");
    let rust_keep = dir.path().join("rust.keep.fq");
    let rust_toss = dir.path().join("rust.toss.fq");

    let shared = discard_bad_only_single_end_args(&input, &java_keep, &java_toss);
    let java_status = Command::new("java")
        .arg("-Xmx1g")
        .arg("-cp")
        .arg(BBTOOLS_CP)
        .arg("jgi.KmerNormalize")
        .args(&shared)
        .output()
        .expect("run java bbnorm");
    assert_success("java bbnorm", &java_status);

    let rust_args = discard_bad_only_single_end_args(&input, &rust_keep, &rust_toss);
    let rust_status = Command::new(env!("CARGO_BIN_EXE_bbnorm-rs"))
        .args(&rust_args)
        .output()
        .expect("run bbnorm-rs");
    assert_success("bbnorm-rs", &rust_status);

    assert_same_file(&java_keep, &rust_keep);
    assert_same_file(&java_toss, &rust_toss);
}

#[test]
fn real_phi_x_duplicated_high_depth_matches_java_bbnorm_for_minlen_toss() {
    assert_min_length_single_end_matches_java("101");
}

#[test]
fn real_phi_x_duplicated_high_depth_matches_java_bbnorm_for_kmg_minlen_toss() {
    assert_min_length_single_end_matches_java("0.101k");
}

#[test]
fn real_phi_x_quality_tail_matches_java_bbnorm_for_qtrim_right() {
    require_file(SAMPLE1);
    require_file("vendor/BBTools-master/current/jgi/KmerNormalize.class");
    require_java();

    let dir = tempdir().expect("create tempdir");
    let input = dir.path().join("real_qtrim_tail.fq");
    write_qtrim_right_fixture(SAMPLE1, &input, 4, 12);

    let java_keep = dir.path().join("java.keep.fq");
    let java_hist = dir.path().join("java.hist.tsv");
    let java_rhist = dir.path().join("java.rhist.tsv");
    let rust_keep = dir.path().join("rust.keep.fq");
    let rust_hist = dir.path().join("rust.hist.tsv");
    let rust_rhist = dir.path().join("rust.rhist.tsv");

    let shared =
        qtrim_keepall_single_end_args(&input, &java_keep, &java_hist, &java_rhist, "qtrim=r");
    let java_status = Command::new("java")
        .arg("-Xmx1g")
        .arg("-cp")
        .arg(BBTOOLS_CP)
        .arg("jgi.KmerNormalize")
        .args(&shared)
        .output()
        .expect("run java bbnorm");
    assert_success("java bbnorm", &java_status);

    let rust_args =
        qtrim_keepall_single_end_args(&input, &rust_keep, &rust_hist, &rust_rhist, "qtrim=r");
    let rust_status = Command::new(env!("CARGO_BIN_EXE_bbnorm-rs"))
        .args(&rust_args)
        .output()
        .expect("run bbnorm-rs");
    assert_success("bbnorm-rs", &rust_status);

    assert_same_file(&java_keep, &rust_keep);
    assert_same_file(&java_hist, &rust_hist);
    assert_same_file(&java_rhist, &rust_rhist);
}

#[test]
fn real_phi_x_quality_tail_matches_java_bbnorm_for_qtrim_window() {
    require_file(SAMPLE1);
    require_file("vendor/BBTools-master/current/jgi/KmerNormalize.class");
    require_java();

    let dir = tempdir().expect("create tempdir");
    let input = dir.path().join("real_qtrim_window_tail.fq");
    write_qtrim_right_fixture(SAMPLE1, &input, 4, 12);

    let java_keep = dir.path().join("java.keep.fq");
    let java_hist = dir.path().join("java.hist.tsv");
    let java_rhist = dir.path().join("java.rhist.tsv");
    let rust_keep = dir.path().join("rust.keep.fq");
    let rust_hist = dir.path().join("rust.hist.tsv");
    let rust_rhist = dir.path().join("rust.rhist.tsv");

    let shared =
        qtrim_keepall_single_end_args(&input, &java_keep, &java_hist, &java_rhist, "qtrim=w,4");
    let java_status = Command::new("java")
        .arg("-Xmx1g")
        .arg("-cp")
        .arg(BBTOOLS_CP)
        .arg("jgi.KmerNormalize")
        .args(&shared)
        .output()
        .expect("run java bbnorm");
    assert_success("java bbnorm", &java_status);

    let rust_args =
        qtrim_keepall_single_end_args(&input, &rust_keep, &rust_hist, &rust_rhist, "qtrim=w,4");
    let rust_status = Command::new(env!("CARGO_BIN_EXE_bbnorm-rs"))
        .args(&rust_args)
        .output()
        .expect("run bbnorm-rs");
    assert_success("bbnorm-rs", &rust_status);

    assert_same_file(&java_keep, &rust_keep);
    assert_same_file(&java_hist, &rust_hist);
    assert_same_file(&java_rhist, &rust_rhist);
}

#[test]
fn representative_qout64_matches_java_bbnorm() {
    require_file("vendor/BBTools-master/current/jgi/KmerNormalize.class");
    require_java();

    let dir = tempdir().expect("create tempdir");
    let input = dir.path().join("representative_qout64.fq");
    write_qout64_fixture(&input);

    let java_keep = dir.path().join("java.keep.fq");
    let rust_keep = dir.path().join("rust.keep.fq");

    let shared = qout64_keepall_args(&input, &java_keep);
    let java_status = Command::new("java")
        .arg("-Xmx1g")
        .arg("-cp")
        .arg(BBTOOLS_CP)
        .arg("jgi.KmerNormalize")
        .args(&shared)
        .output()
        .expect("run java bbnorm");
    assert_success("java bbnorm", &java_status);

    let rust_args = qout64_keepall_args(&input, &rust_keep);
    let rust_status = Command::new(env!("CARGO_BIN_EXE_bbnorm-rs"))
        .args(&rust_args)
        .output()
        .expect("run bbnorm-rs");
    assert_success("bbnorm-rs", &rust_status);

    assert_same_file(&java_keep, &rust_keep);
}

#[test]
fn representative_qin64_matches_java_bbnorm() {
    require_file("vendor/BBTools-master/current/jgi/KmerNormalize.class");
    require_java();

    let dir = tempdir().expect("create tempdir");
    let input = dir.path().join("representative_qin64.fq");
    write_qin64_fixture(&input);

    let java_keep = dir.path().join("java.keep.fq");
    let rust_keep = dir.path().join("rust.keep.fq");

    let shared = qin64_keepall_args(&input, &java_keep);
    let java_status = Command::new("java")
        .arg("-Xmx1g")
        .arg("-cp")
        .arg(BBTOOLS_CP)
        .arg("jgi.KmerNormalize")
        .args(&shared)
        .output()
        .expect("run java bbnorm");
    assert_success("java bbnorm", &java_status);

    let rust_args = qin64_keepall_args(&input, &rust_keep);
    let rust_status = Command::new(env!("CARGO_BIN_EXE_bbnorm-rs"))
        .args(&rust_args)
        .output()
        .expect("run bbnorm-rs");
    assert_success("bbnorm-rs", &rust_status);

    assert_same_file(&java_keep, &rust_keep);
}

#[test]
fn representative_quality_aliases_match_java_bbnorm() {
    require_file("vendor/BBTools-master/current/jgi/KmerNormalize.class");
    require_java();

    let dir = tempdir().expect("create tempdir");
    let input33 = dir.path().join("representative_q33.fq");
    let input64 = dir.path().join("representative_q64.fq");
    write_qout64_fixture(&input33);
    write_qin64_fixture(&input64);

    let cases: [(&str, &Path, &[&str]); 5] = [
        ("ascii64", &input64, &["ascii=64"]),
        ("asciiin64", &input64, &["asciiin=64"]),
        ("asciiout64", &input33, &["asciiout=64", "qin=33"]),
        ("qinauto", &input64, &["qin=auto"]),
        ("qauto", &input64, &["qauto=t"]),
    ];

    for (label, input, extra) in cases {
        let java_keep = dir.path().join(format!("java.{label}.keep.fq"));
        let rust_keep = dir.path().join(format!("rust.{label}.keep.fq"));

        let shared = quality_alias_keepall_args(input, &java_keep, extra);
        let java_status = Command::new("java")
            .arg("-Xmx1g")
            .arg("-cp")
            .arg(BBTOOLS_CP)
            .arg("jgi.KmerNormalize")
            .args(&shared)
            .output()
            .expect("run java bbnorm");
        assert_success("java bbnorm", &java_status);

        let rust_args = quality_alias_keepall_args(input, &rust_keep, extra);
        let rust_status = Command::new(env!("CARGO_BIN_EXE_bbnorm-rs"))
            .args(&rust_args)
            .output()
            .expect("run bbnorm-rs");
        assert_success("bbnorm-rs", &rust_status);

        assert_same_file(&java_keep, &rust_keep);
    }
}

#[test]
fn representative_quality_auto_input_lists_match_java_bbnorm() {
    require_file("vendor/BBTools-master/current/jgi/KmerNormalize.class");
    require_java();

    let dir = tempdir().expect("create tempdir");
    let input33 = dir.path().join("representative_q33_list.fq");
    let input64 = dir.path().join("representative_q64_list.fq");
    write_qout64_fixture(&input33);
    write_qin64_fixture(&input64);

    let cases: [(&str, [&Path; 2], &[&str]); 2] = [
        ("qinauto_q33_then_q64", [&input33, &input64], &["qin=auto"]),
        ("qauto_q64_then_q33", [&input64, &input33], &["qauto=t"]),
    ];

    for (label, inputs, extra) in cases {
        let java_keep_a = dir.path().join(format!("java.{label}.a.fq"));
        let java_keep_b = dir.path().join(format!("java.{label}.b.fq"));
        let java_hist = dir.path().join(format!("java.{label}.hist.tsv"));
        let rust_keep_a = dir.path().join(format!("rust.{label}.a.fq"));
        let rust_keep_b = dir.path().join(format!("rust.{label}.b.fq"));
        let rust_hist = dir.path().join(format!("rust.{label}.hist.tsv"));

        let shared = single_end_input_output_list_args_with_options(
            &inputs,
            &[&java_keep_a, &java_keep_b],
            &java_hist,
            extra,
        );
        let java_status = Command::new("java")
            .arg("-Xmx1g")
            .arg("-cp")
            .arg(BBTOOLS_CP)
            .arg("jgi.KmerNormalize")
            .args(&shared)
            .output()
            .expect("run java bbnorm");
        assert_success("java bbnorm", &java_status);

        let rust_args = single_end_input_output_list_args_with_options(
            &inputs,
            &[&rust_keep_a, &rust_keep_b],
            &rust_hist,
            extra,
        );
        let rust_status = Command::new(env!("CARGO_BIN_EXE_bbnorm-rs"))
            .args(&rust_args)
            .output()
            .expect("run bbnorm-rs");
        assert_success("bbnorm-rs", &rust_status);

        assert_same_file(&java_keep_a, &rust_keep_a);
        assert_same_file(&java_keep_b, &rust_keep_b);
        assert_same_file(&java_hist, &rust_hist);
    }
}

#[test]
fn representative_quality_auto_paired_and_interleaved_match_java_bbnorm() {
    require_file("vendor/BBTools-master/current/jgi/KmerNormalize.class");
    require_java();

    let dir = tempdir().expect("create tempdir");
    let r1 = dir.path().join("representative_qauto_r1.fq");
    let r2 = dir.path().join("representative_qauto_r2.fq");
    let interleaved = dir.path().join("representative_qauto_interleaved.fq");
    write_quality_auto_pair_fixtures(&r1, &r2, &interleaved);

    let paired_cases: [(&str, &[&str]); 1] = [("paired_qinauto", &["qin=auto"])];
    for (label, extra) in paired_cases {
        let java_keep1 = dir.path().join(format!("java.{label}.1.fq"));
        let java_keep2 = dir.path().join(format!("java.{label}.2.fq"));
        let java_hist = dir.path().join(format!("java.{label}.hist.tsv"));
        let rust_keep1 = dir.path().join(format!("rust.{label}.1.fq"));
        let rust_keep2 = dir.path().join(format!("rust.{label}.2.fq"));
        let rust_hist = dir.path().join(format!("rust.{label}.hist.tsv"));

        let shared = paired_input_list_args_with_options(
            &[&r1],
            &[&r2],
            &java_keep1,
            &java_keep2,
            &java_hist,
            extra,
        );
        let java_status = Command::new("java")
            .arg("-Xmx1g")
            .arg("-cp")
            .arg(BBTOOLS_CP)
            .arg("jgi.KmerNormalize")
            .args(&shared)
            .output()
            .expect("run java bbnorm");
        assert_success("java bbnorm", &java_status);

        let rust_args = paired_input_list_args_with_options(
            &[&r1],
            &[&r2],
            &rust_keep1,
            &rust_keep2,
            &rust_hist,
            extra,
        );
        let rust_status = Command::new(env!("CARGO_BIN_EXE_bbnorm-rs"))
            .args(&rust_args)
            .output()
            .expect("run bbnorm-rs");
        assert_success("bbnorm-rs", &rust_status);

        assert_same_file(&java_keep1, &rust_keep1);
        assert_same_file(&java_keep2, &rust_keep2);
        assert_same_file(&java_hist, &rust_hist);
    }

    let interleaved_cases: [(&str, &str, &[&str]); 3] = [
        ("interleaved_explicit_qauto", "interleaved=t", &["qauto=t"]),
        ("interleaved_int_alias_qauto", "int=t", &["qauto=t"]),
        (
            "interleaved_auto_qinauto",
            "interleaved=auto",
            &["qin=auto"],
        ),
    ];
    for (label, interleaved_arg, extra) in interleaved_cases {
        let java_keep = dir.path().join(format!("java.{label}.fq"));
        let java_hist = dir.path().join(format!("java.{label}.hist.tsv"));
        let rust_keep = dir.path().join(format!("rust.{label}.fq"));
        let rust_hist = dir.path().join(format!("rust.{label}.hist.tsv"));

        let shared = interleaved_quality_keepall_args(
            &interleaved,
            &java_keep,
            &java_hist,
            interleaved_arg,
            extra,
        );
        let java_status = Command::new("java")
            .arg("-Xmx1g")
            .arg("-cp")
            .arg(BBTOOLS_CP)
            .arg("jgi.KmerNormalize")
            .args(&shared)
            .output()
            .expect("run java bbnorm");
        assert_success("java bbnorm", &java_status);

        let rust_args = interleaved_quality_keepall_args(
            &interleaved,
            &rust_keep,
            &rust_hist,
            interleaved_arg,
            extra,
        );
        let rust_status = Command::new(env!("CARGO_BIN_EXE_bbnorm-rs"))
            .args(&rust_args)
            .output()
            .expect("run bbnorm-rs");
        assert_success("bbnorm-rs", &rust_status);

        assert_same_file(&java_keep, &rust_keep);
        assert_same_file(&java_hist, &rust_hist);
    }
}

#[test]
fn representative_quality_change_controls_match_java_bbnorm() {
    require_file("vendor/BBTools-master/current/jgi/KmerNormalize.class");
    require_java();

    let dir = tempdir().expect("create tempdir");
    let input = dir.path().join("representative_quality_controls.fq");
    write_qout64_fixture(&input);

    let cases: [(&str, &[&str]); 5] = [
        ("changequality_false", &["changequality=f"]),
        ("cq_false_alias", &["cq=f"]),
        ("ignorebadquality_true", &["ignorebadquality=t"]),
        ("ibq_true", &["ibq=t"]),
        ("min5_max30", &["mincalledquality=5", "maxcalledquality=30"]),
    ];

    for (label, extra) in cases {
        let java_keep = dir.path().join(format!("java.{label}.keep.fq"));
        let rust_keep = dir.path().join(format!("rust.{label}.keep.fq"));

        let shared = quality_alias_keepall_args(&input, &java_keep, extra);
        let java_status = Command::new("java")
            .arg("-Xmx1g")
            .arg("-cp")
            .arg(BBTOOLS_CP)
            .arg("jgi.KmerNormalize")
            .args(&shared)
            .output()
            .expect("run java bbnorm");
        assert_success("java bbnorm", &java_status);

        let rust_args = quality_alias_keepall_args(&input, &rust_keep, extra);
        let rust_status = Command::new(env!("CARGO_BIN_EXE_bbnorm-rs"))
            .args(&rust_args)
            .output()
            .expect("run bbnorm-rs");
        assert_success("bbnorm-rs", &rust_status);

        assert_same_file(&java_keep, &rust_keep);
    }
}

#[test]
fn representative_base_cleanup_matches_java_bbnorm() {
    require_file("vendor/BBTools-master/current/jgi/KmerNormalize.class");
    require_java();

    let dir = tempdir().expect("create tempdir");
    let iupac = dir.path().join("representative_base_iupac.fq");
    let ddx = dir.path().join("representative_base_dotdashx.fq");
    let junk = dir.path().join("representative_base_junk.fq");
    let combo = dir.path().join("representative_base_combo.fq");
    let lower = dir.path().join("representative_base_lower.fq");
    write_fastq_fixture(&iupac, "iupac", "acgtuURYSWKMBDHVNn");
    write_fastq_fixture(&ddx, "ddx", "ACGT.-XxNn");
    write_fastq_fixture(&junk, "junk", "ACGT?ZN");
    write_fastq_fixture(&combo, "combo", "acgtuUnN.-XxRrYy");
    write_fastq_fixture(&lower, "lower", "acgtnu");

    let cases: [(&str, &Path, &[&str]); 19] = [
        ("utot", &iupac, &["utot=t"]),
        ("touppercase", &iupac, &["touppercase=t"]),
        ("tuc_alias", &iupac, &["tuc=t"]),
        ("lowercaseton", &iupac, &["lowercaseton=t"]),
        ("lctn_alias", &iupac, &["lctn=t"]),
        ("iupacton", &iupac, &["iupacton=t"]),
        ("undefinedton_alias", &iupac, &["undefinedton=t"]),
        ("dotdashxton", &ddx, &["dotdashxton=t"]),
        ("fixjunk", &junk, &["fixjunk=t"]),
        ("ignorejunk", &junk, &["ignorejunk=t"]),
        ("flagjunk", &junk, &["flagjunk=t"]),
        ("tossjunk", &junk, &["tossjunk=t"]),
        ("junk_flag", &junk, &["junk=flag"]),
        ("junk_discard", &junk, &["junk=discard"]),
        ("crashjunk_false", &junk, &["crashjunk=f"]),
        ("failjunk_false", &junk, &["failjunk=f"]),
        ("junk_crash_valid", &lower, &["junk=crash"]),
        ("junk_fail_valid", &lower, &["junk=fail"]),
        (
            "lowercaseton_changequality_false",
            &lower,
            &["lowercaseton=t", "changequality=f"],
        ),
    ];

    for (label, input, extra) in cases {
        let java_keep = dir.path().join(format!("java.{label}.keep.fq"));
        let rust_keep = dir.path().join(format!("rust.{label}.keep.fq"));

        let shared = quality_alias_keepall_args(input, &java_keep, extra);
        let java_status = Command::new("java")
            .arg("-Xmx1g")
            .arg("-cp")
            .arg(BBTOOLS_CP)
            .arg("jgi.KmerNormalize")
            .args(&shared)
            .output()
            .expect("run java bbnorm");
        assert_success("java bbnorm", &java_status);

        let rust_args = quality_alias_keepall_args(input, &rust_keep, extra);
        let rust_status = Command::new(env!("CARGO_BIN_EXE_bbnorm-rs"))
            .args(&rust_args)
            .output()
            .expect("run bbnorm-rs");
        assert_success("bbnorm-rs", &rust_status);

        assert_same_file(&java_keep, &rust_keep);
    }

    let combo_cases: [(&str, &[&str]); 3] = [
        ("dotdashxton_iupacton", &["dotdashxton=t", "iupacton=t"]),
        (
            "dotdashxton_iupacton_lctn",
            &["dotdashxton=t", "iupacton=t", "lowercaseton=t"],
        ),
        ("junk_iupacton", &["junk=iupacton"]),
    ];

    for (label, extra) in combo_cases {
        let java_keep = dir.path().join(format!("java.{label}.keep.fq"));
        let rust_keep = dir.path().join(format!("rust.{label}.keep.fq"));

        let shared = quality_alias_keepall_args(&combo, &java_keep, extra);
        let java_status = Command::new("java")
            .arg("-Xmx1g")
            .arg("-cp")
            .arg(BBTOOLS_CP)
            .arg("jgi.KmerNormalize")
            .args(&shared)
            .output()
            .expect("run java bbnorm");
        assert_success("java bbnorm", &java_status);

        let rust_args = quality_alias_keepall_args(&combo, &rust_keep, extra);
        let rust_status = Command::new(env!("CARGO_BIN_EXE_bbnorm-rs"))
            .args(&rust_args)
            .output()
            .expect("run bbnorm-rs");
        assert_success("bbnorm-rs", &rust_status);

        assert_same_file(&java_keep, &rust_keep);
    }
}

#[test]
fn representative_output_format_and_fake_quality_match_java_bbnorm() {
    require_file("vendor/BBTools-master/current/jgi/KmerNormalize.class");
    require_java();

    let dir = tempdir().expect("create tempdir");
    let fasta = dir.path().join("representative_output.fa");
    let fastq = dir.path().join("representative_output.fq");
    let long_fastq = dir.path().join("representative_output_wrap.fq");
    write_fasta_fixture(&fasta, &[("seq1", "ACGTNN"), ("seq2", "TTTTAC")]);
    write_fastq_fixture(&fastq, "seq1", "ACGTNN");
    let long_bases = "ACGT".repeat(25);
    write_fastq_fixture(&long_fastq, "long", &long_bases);

    let cases: [(&str, &Path, &str, &[&str]); 18] = [
        ("fasta_to_fastq_default", &fasta, "keep.fq", &[]),
        ("fasta_to_fastq_qout64", &fasta, "keep.fq", &["qout=64"]),
        (
            "fasta_to_fastq_fakefasta20",
            &fasta,
            "keep.fq",
            &["fakefastaquality=20"],
        ),
        ("fasta_to_fastq_qfake15", &fasta, "keep.fq", &["qfake=15"]),
        ("fasta_to_fasta", &fasta, "keep.fa", &[]),
        (
            "fasta_to_fasta_fastareadlen",
            &fasta,
            "keep.fa",
            &["fastareadlen=4"],
        ),
        (
            "fasta_to_fasta_fastareadlength",
            &fasta,
            "keep.fa",
            &["fastareadlength=4"],
        ),
        (
            "fasta_to_fasta_fastaminread",
            &fasta,
            "keep.fa",
            &["fastaminread=1"],
        ),
        (
            "fasta_to_fasta_fastaminlen",
            &fasta,
            "keep.fa",
            &["fastaminlen=1"],
        ),
        (
            "fasta_to_fasta_fastaminlength",
            &fasta,
            "keep.fa",
            &["fastaminlength=1"],
        ),
        (
            "fasta_to_fasta_forcesectionname",
            &fasta,
            "keep.fa",
            &["forcesectionname=t"],
        ),
        (
            "fasta_to_fasta_fastadump_false",
            &fasta,
            "keep.fa",
            &["fastadump=f"],
        ),
        ("fasta_to_unknown_txt", &fasta, "keep.txt", &[]),
        ("fastq_to_fasta", &fastq, "keep.fa", &[]),
        ("fastq_to_fasta_default_wrap", &long_fastq, "keep.fa", &[]),
        (
            "fastq_to_fasta_fastawrap20",
            &long_fastq,
            "keep.fa",
            &["fastawrap=20"],
        ),
        (
            "fastq_to_fasta_wrap_alias20",
            &long_fastq,
            "keep.fa",
            &["wrap=20"],
        ),
        (
            "fastq_to_fasta_fastawrap0",
            &long_fastq,
            "keep.fa",
            &["fastawrap=0"],
        ),
    ];

    for (label, input, output_name, extra) in cases {
        let java_keep = dir.path().join(format!("java.{label}.{output_name}"));
        let rust_keep = dir.path().join(format!("rust.{label}.{output_name}"));

        let shared = quality_alias_keepall_args(input, &java_keep, extra);
        let java_status = Command::new("java")
            .arg("-Xmx1g")
            .arg("-cp")
            .arg(BBTOOLS_CP)
            .arg("jgi.KmerNormalize")
            .args(&shared)
            .output()
            .expect("run java bbnorm");
        assert_success("java bbnorm", &java_status);

        let rust_args = quality_alias_keepall_args(input, &rust_keep, extra);
        let rust_status = Command::new(env!("CARGO_BIN_EXE_bbnorm-rs"))
            .args(&rust_args)
            .output()
            .expect("run bbnorm-rs");
        assert_success("bbnorm-rs", &rust_status);

        assert_same_file(&java_keep, &rust_keep);
    }
}

#[test]
fn representative_fasta_parser_numeric_controls_reject_malformed_like_java_bbnorm() {
    require_file("vendor/BBTools-master/current/jgi/KmerNormalize.class");
    require_java();

    let dir = tempdir().expect("create tempdir");
    let input = dir.path().join("fasta_parser_malformed.fa");
    write_fasta_fixture(&input, &[("seq1", "ACGTNN"), ("seq2", "TTTTAC")]);

    for option in ["fastaminread=abc", "fastaminlen=abc", "fastaminlength=abc"] {
        let java_keep = dir
            .path()
            .join(format!("java.{}.keep.fa", option.replace('=', "_")));
        let rust_keep = dir
            .path()
            .join(format!("rust.{}.keep.fa", option.replace('=', "_")));

        let java_args = quality_alias_keepall_args(&input, &java_keep, &[option]);
        let java_status = Command::new("java")
            .arg("-Xmx1g")
            .arg("-cp")
            .arg(BBTOOLS_CP)
            .arg("jgi.KmerNormalize")
            .args(&java_args)
            .output()
            .expect("run java bbnorm");
        assert!(
            !java_status.status.success(),
            "java should reject malformed {option}; stdout:\n{}\nstderr:\n{}",
            String::from_utf8_lossy(&java_status.stdout),
            String::from_utf8_lossy(&java_status.stderr)
        );

        let rust_args = quality_alias_keepall_args(&input, &rust_keep, &[option]);
        let rust_status = Command::new(env!("CARGO_BIN_EXE_bbnorm-rs"))
            .args(&rust_args)
            .output()
            .expect("run bbnorm-rs");
        assert!(
            !rust_status.status.success(),
            "bbnorm-rs should reject malformed {option}; stdout:\n{}\nstderr:\n{}",
            String::from_utf8_lossy(&rust_status.stdout),
            String::from_utf8_lossy(&rust_status.stderr)
        );
        assert!(!rust_keep.exists());
    }
}

#[test]
fn representative_append_outputs_match_java_bbnorm() {
    require_file("vendor/BBTools-master/current/jgi/KmerNormalize.class");
    require_java();

    let dir = tempdir().expect("create tempdir");
    let input = dir.path().join("append_input.fq");
    write_fastq_fixture(&input, "append_probe", "ACGTACGT");

    let java_keep = dir.path().join("java.append.keep.fq");
    let rust_keep = dir.path().join("rust.append.keep.fq");
    let java_hist = dir.path().join("java.append.hist.tsv");
    let rust_hist = dir.path().join("rust.append.hist.tsv");
    fs::write(&java_keep, b"PREFASTQ\n").expect("seed java keep");
    fs::write(&rust_keep, b"PREFASTQ\n").expect("seed rust keep");
    fs::write(&java_hist, b"PREHIST\n").expect("seed java hist");
    fs::write(&rust_hist, b"PREHIST\n").expect("seed rust hist");

    let shared = append_keepall_args(&input, &java_keep, &java_hist);
    let java_status = Command::new("java")
        .arg("-Xmx1g")
        .arg("-cp")
        .arg(BBTOOLS_CP)
        .arg("jgi.KmerNormalize")
        .args(&shared)
        .output()
        .expect("run java bbnorm");
    assert_success("java bbnorm", &java_status);

    let rust_args = append_keepall_args(&input, &rust_keep, &rust_hist);
    let rust_status = Command::new(env!("CARGO_BIN_EXE_bbnorm-rs"))
        .args(&rust_args)
        .output()
        .expect("run bbnorm-rs");
    assert_success("bbnorm-rs", &rust_status);

    assert_same_file(&java_keep, &rust_keep);
    assert_same_file(&java_hist, &rust_hist);
    assert!(
        fs::read_to_string(&rust_keep)
            .expect("read appended keep")
            .starts_with("PREFASTQ\n@append_probe\n")
    );
    assert!(
        !fs::read_to_string(&rust_hist)
            .expect("read overwritten hist")
            .starts_with("PREHIST\n")
    );
}

#[test]
fn representative_default_multipass_matches_java_bbnorm() {
    require_file("vendor/BBTools-master/current/jgi/KmerNormalize.class");
    require_java();

    let dir = tempdir().expect("create tempdir");
    let input = dir.path().join("default_multipass_input.fq");
    write_fastq_fixture(&input, "default_pass_probe", "ACGTACGT");

    let java_keep = dir.path().join("java.default.keep.fq");
    let java_hist = dir.path().join("java.default.hist.tsv");
    let rust_keep = dir.path().join("rust.default.keep.fq");
    let rust_hist = dir.path().join("rust.default.hist.tsv");

    let java_args = default_multipass_fallback_args(&input, &java_keep, &java_hist);
    let java_status = Command::new("java")
        .arg("-Xmx1g")
        .arg("-cp")
        .arg(BBTOOLS_CP)
        .arg("jgi.KmerNormalize")
        .args(&java_args)
        .output()
        .expect("run java bbnorm");
    assert_success("java bbnorm", &java_status);

    let rust_args = default_multipass_fallback_args(&input, &rust_keep, &rust_hist);
    let rust_status = Command::new(env!("CARGO_BIN_EXE_bbnorm-rs"))
        .args(&rust_args)
        .output()
        .expect("run bbnorm-rs");
    assert_success("bbnorm-rs default multipass", &rust_status);
    assert_same_file(&java_keep, &rust_keep);
    assert_same_file(&java_hist, &rust_hist);
}

#[test]
fn representative_wrapper_sampling_options_fall_back_to_supported_normalization() {
    require_file("vendor/BBTools-master/current/jgi/KmerNormalize.class");
    require_java();

    let dir = tempdir().expect("create tempdir");
    let input = dir.path().join("sampling_guard_input.fq");
    write_fastq_fixture(&input, "sampling_guard", "ACGTACGT");
    let java_base_keep = dir.path().join("java.base.keep.fq");
    let java_base_args = sampling_guard_args(&input, &java_base_keep, "");
    let java_base_status = Command::new("java")
        .arg("-Xmx1g")
        .arg("-cp")
        .arg(BBTOOLS_CP)
        .arg("jgi.KmerNormalize")
        .args(java_base_args.into_iter().filter(|arg| !arg.is_empty()))
        .output()
        .expect("run java bbnorm baseline");
    assert_success("java bbnorm sampling baseline", &java_base_status);

    for option in [
        "sampleoutput=1",
        "readsample=1",
        "kmersample=1",
        "samplerate=0.5",
        "sample=0.5",
        "sampleseed=1",
        "seed=1",
    ] {
        let java_keep = dir
            .path()
            .join(format!("java.{}.keep.fq", option.replace('=', "_")));
        let rust_keep = dir
            .path()
            .join(format!("rust.{}.keep.fq", option.replace('=', "_")));

        let java_args = sampling_guard_args(&input, &java_keep, option);
        let java_status = Command::new("java")
            .arg("-Xmx1g")
            .arg("-cp")
            .arg(BBTOOLS_CP)
            .arg("jgi.KmerNormalize")
            .args(&java_args)
            .output()
            .expect("run java bbnorm");
        assert!(
            !java_status.status.success(),
            "java should reject {option}; stdout:\n{}\nstderr:\n{}",
            String::from_utf8_lossy(&java_status.stdout),
            String::from_utf8_lossy(&java_status.stderr)
        );
        assert!(
            String::from_utf8_lossy(&java_status.stderr).contains("Unknown parameter"),
            "unexpected java stderr for {option}: {}",
            String::from_utf8_lossy(&java_status.stderr)
        );

        let rust_args = sampling_guard_args(&input, &rust_keep, option);
        let rust_status = Command::new(env!("CARGO_BIN_EXE_bbnorm-rs"))
            .args(&rust_args)
            .output()
            .expect("run bbnorm-rs");
        assert_success("bbnorm-rs sampling fallback", &rust_status);
        assert!(
            String::from_utf8_lossy(&rust_status.stderr).contains("Rust ignores it"),
            "unexpected rust stderr for {option}: {}",
            String::from_utf8_lossy(&rust_status.stderr)
        );
        assert!(!java_keep.exists());
        assert_same_file(&java_base_keep, &rust_keep);
    }
}

#[test]
fn representative_inactive_numeric_options_reject_malformed_like_java_bbnorm() {
    require_file("vendor/BBTools-master/current/jgi/KmerNormalize.class");
    require_java();

    let dir = tempdir().expect("create tempdir");
    let input = dir.path().join("multipass_only_malformed_input.fq");
    write_fastq_fixture(&input, "multipass_only_malformed", "ACGTACGT");

    for option in [
        "target1=abc",
        "targetbadpercentilelow=abc",
        "tbph=abc",
        "bits1=abc",
        "cbits1=abc",
        "cellbits1=abc",
        "stepsize=abc",
        "buildstepsize=abc",
    ] {
        let java_keep = dir
            .path()
            .join(format!("java.{}.keep.fq", option.replace('=', "_")));
        let rust_keep = dir
            .path()
            .join(format!("rust.{}.keep.fq", option.replace('=', "_")));

        let java_args = quality_alias_keepall_args(&input, &java_keep, &[option]);
        let java_status = Command::new("java")
            .arg("-Xmx1g")
            .arg("-cp")
            .arg(BBTOOLS_CP)
            .arg("jgi.KmerNormalize")
            .args(&java_args)
            .output()
            .expect("run java bbnorm");
        assert!(
            !java_status.status.success(),
            "java should reject malformed {option}; stdout:\n{}\nstderr:\n{}",
            String::from_utf8_lossy(&java_status.stdout),
            String::from_utf8_lossy(&java_status.stderr)
        );

        let rust_args = quality_alias_keepall_args(&input, &rust_keep, &[option]);
        let rust_status = Command::new(env!("CARGO_BIN_EXE_bbnorm-rs"))
            .args(&rust_args)
            .output()
            .expect("run bbnorm-rs");
        assert!(
            !rust_status.status.success(),
            "bbnorm-rs should reject malformed {option}; stdout:\n{}\nstderr:\n{}",
            String::from_utf8_lossy(&rust_status.stdout),
            String::from_utf8_lossy(&rust_status.stderr)
        );
        assert!(!rust_keep.exists());
    }
}

#[test]
fn representative_inactive_ecc_numeric_options_reject_malformed_like_java_bbnorm() {
    require_file("vendor/BBTools-master/current/jgi/KmerNormalize.class");
    require_java();

    let dir = tempdir().expect("create tempdir");
    let input = dir.path().join("ecc_tuning_malformed_input.fq");
    write_fastq_fixture(&input, "ecc_tuning_malformed", "ACGTACGT");

    for option in [
        "ecclimit=abc",
        "eccmaxqual=abc",
        "ecr=abc",
        "echthresh=abc",
        "eclt=abc",
        "suflen=abc",
        "prelen=abc",
    ] {
        let java_keep = dir
            .path()
            .join(format!("java.{}.keep.fq", option.replace('=', "_")));
        let rust_keep = dir
            .path()
            .join(format!("rust.{}.keep.fq", option.replace('=', "_")));

        let java_args = quality_alias_keepall_args(&input, &java_keep, &[option]);
        let java_status = Command::new("java")
            .arg("-Xmx1g")
            .arg("-cp")
            .arg(BBTOOLS_CP)
            .arg("jgi.KmerNormalize")
            .args(&java_args)
            .output()
            .expect("run java bbnorm");
        assert!(
            !java_status.status.success(),
            "java should reject malformed {option}; stdout:\n{}\nstderr:\n{}",
            String::from_utf8_lossy(&java_status.stdout),
            String::from_utf8_lossy(&java_status.stderr)
        );

        let rust_args = quality_alias_keepall_args(&input, &rust_keep, &[option]);
        let rust_status = Command::new(env!("CARGO_BIN_EXE_bbnorm-rs"))
            .args(&rust_args)
            .output()
            .expect("run bbnorm-rs");
        assert!(
            !rust_status.status.success(),
            "bbnorm-rs should reject malformed {option}; stdout:\n{}\nstderr:\n{}",
            String::from_utf8_lossy(&rust_status.stdout),
            String::from_utf8_lossy(&rust_status.stderr)
        );
        assert!(!rust_keep.exists());
    }
}

#[test]
fn representative_shared_io_runtime_controls_reject_malformed_like_java_bbnorm() {
    require_file("vendor/BBTools-master/current/jgi/KmerNormalize.class");
    require_java();

    let dir = tempdir().expect("create tempdir");
    let input = dir.path().join("shared_io_runtime_malformed_input.fq");
    write_fastq_fixture(&input, "shared_io_runtime_malformed", "ACGTACGT");

    for option in [
        "bf1bufferlen=abc",
        "bfthreads=abc",
        "readbufferlength=abc",
        "readbuffers=abc",
    ] {
        let java_keep = dir
            .path()
            .join(format!("java.{}.keep.fq", option.replace('=', "_")));
        let rust_keep = dir
            .path()
            .join(format!("rust.{}.keep.fq", option.replace('=', "_")));

        let java_args = quality_alias_keepall_args(&input, &java_keep, &[option]);
        let java_status = Command::new("java")
            .arg("-Xmx1g")
            .arg("-cp")
            .arg(BBTOOLS_CP)
            .arg("jgi.KmerNormalize")
            .args(&java_args)
            .output()
            .expect("run java bbnorm");
        assert!(
            !java_status.status.success(),
            "java should reject malformed {option}; stdout:\n{}\nstderr:\n{}",
            String::from_utf8_lossy(&java_status.stdout),
            String::from_utf8_lossy(&java_status.stderr)
        );

        let rust_args = quality_alias_keepall_args(&input, &rust_keep, &[option]);
        let rust_status = Command::new(env!("CARGO_BIN_EXE_bbnorm-rs"))
            .args(&rust_args)
            .output()
            .expect("run bbnorm-rs");
        assert!(
            !rust_status.status.success(),
            "bbnorm-rs should reject malformed {option}; stdout:\n{}\nstderr:\n{}",
            String::from_utf8_lossy(&rust_status.stdout),
            String::from_utf8_lossy(&rust_status.stderr)
        );
        assert!(!rust_keep.exists());
    }
}

#[test]
fn representative_shared_environment_runtime_controls_reject_malformed_like_java_bbnorm() {
    require_file("vendor/BBTools-master/current/jgi/KmerNormalize.class");
    require_java();

    let dir = tempdir().expect("create tempdir");
    let input = dir.path().join("shared_environment_malformed_input.fq");
    write_fastq_fixture(&input, "shared_environment_malformed", "ACGTACGT");

    for option in ["entropyk=abc", "entropywindow=abc"] {
        let java_keep = dir
            .path()
            .join(format!("java.{}.keep.fq", option.replace('=', "_")));
        let rust_keep = dir
            .path()
            .join(format!("rust.{}.keep.fq", option.replace('=', "_")));

        let java_args = quality_alias_keepall_args(&input, &java_keep, &[option]);
        let java_status = Command::new("java")
            .arg("-Xmx1g")
            .arg("-cp")
            .arg(BBTOOLS_CP)
            .arg("jgi.KmerNormalize")
            .args(&java_args)
            .output()
            .expect("run java bbnorm");
        assert!(
            !java_status.status.success(),
            "java should reject malformed {option}; stdout:\n{}\nstderr:\n{}",
            String::from_utf8_lossy(&java_status.stdout),
            String::from_utf8_lossy(&java_status.stderr)
        );

        let rust_args = quality_alias_keepall_args(&input, &rust_keep, &[option]);
        let rust_status = Command::new(env!("CARGO_BIN_EXE_bbnorm-rs"))
            .args(&rust_args)
            .output()
            .expect("run bbnorm-rs");
        assert!(
            !rust_status.status.success(),
            "bbnorm-rs should reject malformed {option}; stdout:\n{}\nstderr:\n{}",
            String::from_utf8_lossy(&rust_status.stdout),
            String::from_utf8_lossy(&rust_status.stderr)
        );
        assert!(!rust_keep.exists());
    }
}

#[test]
fn representative_extra_inputs_match_java_bbnorm() {
    require_file("vendor/BBTools-master/current/jgi/KmerNormalize.class");
    require_java();

    let dir = tempdir().expect("create tempdir");
    let input = dir.path().join("representative_extra_main.fq");
    let extra1 = dir.path().join("representative_extra_1.fq");
    let extra2 = dir.path().join("representative_extra_2.fq");
    write_fastq_records(&input, &[("main", "ACGTACGT")]);
    write_fastq_records(&extra1, &[("extra1", "ACGTACGT"), ("extra2", "ACGTACGT")]);
    write_fastq_records(&extra2, &[("extra3", "TTTTTTTT")]);

    let java_keep = dir.path().join("java.keep.fq");
    let java_hist = dir.path().join("java.hist.tsv");
    let rust_keep = dir.path().join("rust.keep.fq");
    let rust_hist = dir.path().join("rust.hist.tsv");

    let shared = extra_input_keepall_args(&input, &[&extra1, &extra2], &java_keep, &java_hist);
    let java_status = Command::new("java")
        .arg("-Xmx1g")
        .arg("-cp")
        .arg(BBTOOLS_CP)
        .arg("jgi.KmerNormalize")
        .args(&shared)
        .output()
        .expect("run java bbnorm");
    assert_success("java bbnorm", &java_status);

    let rust_args = extra_input_keepall_args(&input, &[&extra1, &extra2], &rust_keep, &rust_hist);
    let rust_status = Command::new(env!("CARGO_BIN_EXE_bbnorm-rs"))
        .args(&rust_args)
        .output()
        .expect("run bbnorm-rs");
    assert_success("bbnorm-rs", &rust_status);

    assert_same_file(&java_keep, &rust_keep);
    assert_same_file(&java_hist, &rust_hist);
}

#[test]
fn representative_literal_comma_extra_input_matches_java_bbnorm() {
    require_file("vendor/BBTools-master/current/jgi/KmerNormalize.class");
    require_java();

    let dir = tempdir().expect("create tempdir");
    let input = dir.path().join("representative_extra_main.fq");
    let extra = dir.path().join("representative,extra,literal.fq");
    write_fastq_records(&input, &[("main", "ACGTACGT")]);
    write_fastq_records(&extra, &[("extra1", "ACGTACGT"), ("extra2", "ACGTACGT")]);

    let java_keep = dir.path().join("java.keep.fq");
    let java_hist = dir.path().join("java.hist.tsv");
    let rust_keep = dir.path().join("rust.keep.fq");
    let rust_hist = dir.path().join("rust.hist.tsv");

    let shared = extra_input_keepall_args(&input, &[&extra], &java_keep, &java_hist);
    let java_status = Command::new("java")
        .arg("-Xmx1g")
        .arg("-cp")
        .arg(BBTOOLS_CP)
        .arg("jgi.KmerNormalize")
        .args(&shared)
        .output()
        .expect("run java bbnorm");
    assert_success("java bbnorm", &java_status);

    let rust_args = extra_input_keepall_args(&input, &[&extra], &rust_keep, &rust_hist);
    let rust_status = Command::new(env!("CARGO_BIN_EXE_bbnorm-rs"))
        .args(&rust_args)
        .output()
        .expect("run bbnorm-rs");
    assert_success("bbnorm-rs", &rust_status);

    assert_same_file(&java_keep, &rust_keep);
    assert_same_file(&java_hist, &rust_hist);
}

#[test]
fn representative_extra_inputs_ignore_tablereads_limit_like_java_bbnorm() {
    require_file("vendor/BBTools-master/current/jgi/KmerNormalize.class");
    require_java();

    let dir = tempdir().expect("create tempdir");
    let input = dir.path().join("representative_extra_tablereads_main.fq");
    let extra = dir.path().join("representative_extra_tablereads_extra.fq");
    write_fastq_records(&input, &[("main1", "ACGTACGT"), ("main2", "ACGTACGT")]);
    write_fastq_records(&extra, &[("extra1", "ACGTACGT"), ("extra2", "ACGTACGT")]);

    let java_keep = dir.path().join("java.keep.fq");
    let java_hist = dir.path().join("java.hist.tsv");
    let rust_keep = dir.path().join("rust.keep.fq");
    let rust_hist = dir.path().join("rust.hist.tsv");

    let shared = extra_input_keepall_args_with_options(
        &input,
        &[&extra],
        &java_keep,
        &java_hist,
        &["tablereads=1"],
    );
    let java_status = Command::new("java")
        .arg("-Xmx1g")
        .arg("-cp")
        .arg(BBTOOLS_CP)
        .arg("jgi.KmerNormalize")
        .args(&shared)
        .output()
        .expect("run java bbnorm");
    assert_success("java bbnorm", &java_status);

    let rust_args = extra_input_keepall_args_with_options(
        &input,
        &[&extra],
        &rust_keep,
        &rust_hist,
        &["tablereads=1"],
    );
    let rust_status = Command::new(env!("CARGO_BIN_EXE_bbnorm-rs"))
        .args(&rust_args)
        .output()
        .expect("run bbnorm-rs");
    assert_success("bbnorm-rs", &rust_status);

    assert_same_file(&java_keep, &rust_keep);
    assert_same_file(&java_hist, &rust_hist);
}

#[test]
fn representative_missing_hash_extra_input_is_rejected_like_java_bbnorm() {
    require_file("vendor/BBTools-master/current/jgi/KmerNormalize.class");
    require_java();

    let dir = tempdir().expect("create tempdir");
    let input = dir.path().join("representative_extra_hash_main.fq");
    let extra_hash = dir.path().join("representative_extra_hash#.fq");
    let extra1 = dir.path().join("representative_extra_hash1.fq");
    let extra2 = dir.path().join("representative_extra_hash2.fq");
    write_fastq_records(&input, &[("main", "ACGTACGT")]);
    write_fastq_records(&extra1, &[("extra1", "ACGTACGT")]);
    write_fastq_records(&extra2, &[("extra2", "ACGTACGT")]);

    let java_keep = dir.path().join("java.keep.fq");
    let java_hist = dir.path().join("java.hist.tsv");
    let rust_keep = dir.path().join("rust.keep.fq");
    let rust_hist = dir.path().join("rust.hist.tsv");

    let shared = extra_input_keepall_args(&input, &[&extra_hash], &java_keep, &java_hist);
    let java_status = Command::new("java")
        .arg("-Xmx1g")
        .arg("-cp")
        .arg(BBTOOLS_CP)
        .arg("jgi.KmerNormalize")
        .args(&shared)
        .output()
        .expect("run java bbnorm");
    assert!(
        !java_status.status.success(),
        "java bbnorm unexpectedly accepted missing hash extra input"
    );

    let rust_args = extra_input_keepall_args(&input, &[&extra_hash], &rust_keep, &rust_hist);
    let rust_status = Command::new(env!("CARGO_BIN_EXE_bbnorm-rs"))
        .args(&rust_args)
        .output()
        .expect("run bbnorm-rs");
    assert!(
        !rust_status.status.success(),
        "bbnorm-rs unexpectedly accepted missing hash extra input"
    );
}

#[test]
fn representative_single_end_input_list_matches_java_bbnorm() {
    require_file("vendor/BBTools-master/current/jgi/KmerNormalize.class");
    require_java();

    let dir = tempdir().expect("create tempdir");
    let input1 = dir.path().join("representative_list_a.fq");
    let input2 = dir.path().join("representative_list_b.fq");
    write_fastq_records(&input1, &[("a1", "ACGTACGT"), ("a2", "CCCCCCCC")]);
    write_fastq_records(&input2, &[("b1", "TTTTTTTT"), ("b2", "GGGGGGGG")]);

    let java_keep = dir.path().join("java.keep.fq");
    let java_hist = dir.path().join("java.hist.tsv");
    let rust_keep = dir.path().join("rust.keep.fq");
    let rust_hist = dir.path().join("rust.hist.tsv");

    let shared = single_end_input_list_args(&[&input1, &input2], &java_keep, &java_hist);
    let java_status = Command::new("java")
        .arg("-Xmx1g")
        .arg("-cp")
        .arg(BBTOOLS_CP)
        .arg("jgi.KmerNormalize")
        .args(&shared)
        .output()
        .expect("run java bbnorm");
    assert_success("java bbnorm", &java_status);

    let rust_args = single_end_input_list_args(&[&input1, &input2], &rust_keep, &rust_hist);
    let rust_status = Command::new(env!("CARGO_BIN_EXE_bbnorm-rs"))
        .args(&rust_args)
        .output()
        .expect("run bbnorm-rs");
    assert_success("bbnorm-rs", &rust_status);

    assert_same_file(&java_keep, &rust_keep);
    assert_same_file(&java_hist, &rust_hist);
}

#[test]
fn representative_paired_input_list_matches_java_bbnorm() {
    require_file("vendor/BBTools-master/current/jgi/KmerNormalize.class");
    require_java();

    let dir = tempdir().expect("create tempdir");
    let r1a = dir.path().join("representative_pair_list_a_1.fq");
    let r2a = dir.path().join("representative_pair_list_a_2.fq");
    let r1b = dir.path().join("representative_pair_list_b_1.fq");
    let r2b = dir.path().join("representative_pair_list_b_2.fq");
    write_fastq_records(&r1a, &[("a1/1", "AAA"), ("a2/1", "CCC")]);
    write_fastq_records(&r2a, &[("a1/2", "GGG"), ("a2/2", "TTT")]);
    write_fastq_records(&r1b, &[("b1/1", "ACG"), ("b2/1", "CGT")]);
    write_fastq_records(&r2b, &[("b1/2", "GTA"), ("b2/2", "TAC")]);

    let java_keep1 = dir.path().join("java.keep1.fq");
    let java_keep2 = dir.path().join("java.keep2.fq");
    let java_hist = dir.path().join("java.hist.tsv");
    let rust_keep1 = dir.path().join("rust.keep1.fq");
    let rust_keep2 = dir.path().join("rust.keep2.fq");
    let rust_hist = dir.path().join("rust.hist.tsv");

    let shared = paired_input_list_args(
        &[&r1a, &r1b],
        &[&r2a, &r2b],
        &java_keep1,
        &java_keep2,
        &java_hist,
    );
    let java_status = Command::new("java")
        .arg("-Xmx1g")
        .arg("-cp")
        .arg(BBTOOLS_CP)
        .arg("jgi.KmerNormalize")
        .args(&shared)
        .output()
        .expect("run java bbnorm");
    assert_success("java bbnorm", &java_status);

    let rust_args = paired_input_list_args(
        &[&r1a, &r1b],
        &[&r2a, &r2b],
        &rust_keep1,
        &rust_keep2,
        &rust_hist,
    );
    let rust_status = Command::new(env!("CARGO_BIN_EXE_bbnorm-rs"))
        .args(&rust_args)
        .output()
        .expect("run bbnorm-rs");
    assert_success("bbnorm-rs", &rust_status);

    assert_same_file(&java_keep1, &rust_keep1);
    assert_same_file(&java_keep2, &rust_keep2);
    assert_same_file(&java_hist, &rust_hist);
}

#[test]
fn representative_single_end_input_list_output_fanout_matches_java_bbnorm() {
    require_file("vendor/BBTools-master/current/jgi/KmerNormalize.class");
    require_java();

    let dir = tempdir().expect("create tempdir");
    let input1 = dir.path().join("representative_list_a.fq");
    let input2 = dir.path().join("representative_list_b.fq");
    write_fastq_records(&input1, &[("a1", "ACGTACGT"), ("a2", "CCCCCCCC")]);
    write_fastq_records(&input2, &[("b1", "TTTTTTTT"), ("b2", "GGGGGGGG")]);

    let java_keep_a = dir.path().join("java.keep.a.fq");
    let java_keep_b = dir.path().join("java.keep.b.fq");
    let java_hist = dir.path().join("java.output_list.hist.tsv");
    let rust_keep_a = dir.path().join("rust.keep.a.fq");
    let rust_keep_b = dir.path().join("rust.keep.b.fq");
    let rust_hist = dir.path().join("rust.output_list.hist.tsv");

    let shared = single_end_input_output_list_args(
        &[&input1, &input2],
        &[&java_keep_a, &java_keep_b],
        &java_hist,
    );
    let java_status = Command::new("java")
        .arg("-Xmx1g")
        .arg("-cp")
        .arg(BBTOOLS_CP)
        .arg("jgi.KmerNormalize")
        .args(&shared)
        .output()
        .expect("run java bbnorm");
    assert_success("java bbnorm", &java_status);

    let rust_args = single_end_input_output_list_args(
        &[&input1, &input2],
        &[&rust_keep_a, &rust_keep_b],
        &rust_hist,
    );
    let rust_status = Command::new(env!("CARGO_BIN_EXE_bbnorm-rs"))
        .args(&rust_args)
        .output()
        .expect("run bbnorm-rs");
    assert_success("bbnorm-rs", &rust_status);

    assert_same_file(&java_keep_a, &rust_keep_a);
    assert_same_file(&java_keep_b, &rust_keep_b);
    assert_same_file(&java_hist, &rust_hist);
}

#[test]
fn representative_single_end_input_list_low_bin_fanout_matches_java_bbnorm() {
    require_file("vendor/BBTools-master/current/jgi/KmerNormalize.class");
    require_java();

    let dir = tempdir().expect("create tempdir");
    let input1 = dir.path().join("representative_list_a.fq");
    let input2 = dir.path().join("representative_list_b.fq");
    write_fastq_records(&input1, &[("a1", "ACGTACGT"), ("a2", "CCCCCCCC")]);
    write_fastq_records(&input2, &[("b1", "TTTTTTTT"), ("b2", "GGGGGGGG")]);

    let java_keep = dir.path().join("java.keep.fq");
    let java_low_a = dir.path().join("java.low.a.fq");
    let java_low_b = dir.path().join("java.low.b.fq");
    let java_hist = dir.path().join("java.low_output_list.hist.tsv");
    let rust_keep = dir.path().join("rust.keep.fq");
    let rust_low_a = dir.path().join("rust.low.a.fq");
    let rust_low_b = dir.path().join("rust.low.b.fq");
    let rust_hist = dir.path().join("rust.low_output_list.hist.tsv");

    let shared = single_end_input_low_bin_output_list_args(
        &[&input1, &input2],
        &java_keep,
        &[&java_low_a, &java_low_b],
        &java_hist,
    );
    let java_status = Command::new("java")
        .arg("-Xmx1g")
        .arg("-cp")
        .arg(BBTOOLS_CP)
        .arg("jgi.KmerNormalize")
        .args(&shared)
        .output()
        .expect("run java bbnorm");
    assert_success("java bbnorm", &java_status);

    let rust_args = single_end_input_low_bin_output_list_args(
        &[&input1, &input2],
        &rust_keep,
        &[&rust_low_a, &rust_low_b],
        &rust_hist,
    );
    let rust_status = Command::new(env!("CARGO_BIN_EXE_bbnorm-rs"))
        .args(&rust_args)
        .output()
        .expect("run bbnorm-rs");
    assert_success("bbnorm-rs", &rust_status);

    assert_same_file(&java_keep, &rust_keep);
    assert_same_file(&java_low_a, &rust_low_a);
    assert_same_file(&java_low_b, &rust_low_b);
    assert_same_file(&java_hist, &rust_hist);
}

#[test]
fn representative_single_end_input_list_mid_and_high_bin_fanout_match_java_bbnorm() {
    require_file("vendor/BBTools-master/current/jgi/KmerNormalize.class");
    require_java();

    let dir = tempdir().expect("create tempdir");
    let input1 = dir.path().join("representative_list_a.fq");
    let input2 = dir.path().join("representative_list_b.fq");
    write_fastq_records(&input1, &[("a1", "ACGTACGT"), ("a2", "CCCCCCCC")]);
    write_fastq_records(&input2, &[("b1", "TTTTTTTT"), ("b2", "GGGGGGGG")]);

    for (label, output_key, high_bin_depth) in [("mid", "outmid", "80"), ("high", "outhigh", "0")] {
        let java_keep = dir.path().join(format!("java.{label}.keep.fq"));
        let java_bin_a = dir.path().join(format!("java.{label}.a.fq"));
        let java_bin_b = dir.path().join(format!("java.{label}.b.fq"));
        let java_hist = dir
            .path()
            .join(format!("java.{label}_output_list.hist.tsv"));
        let rust_keep = dir.path().join(format!("rust.{label}.keep.fq"));
        let rust_bin_a = dir.path().join(format!("rust.{label}.a.fq"));
        let rust_bin_b = dir.path().join(format!("rust.{label}.b.fq"));
        let rust_hist = dir
            .path()
            .join(format!("rust.{label}_output_list.hist.tsv"));

        let shared = single_end_input_depth_bin_output_list_args(
            &[&input1, &input2],
            &java_keep,
            output_key,
            &[&java_bin_a, &java_bin_b],
            &java_hist,
            high_bin_depth,
        );
        let java_status = Command::new("java")
            .arg("-Xmx1g")
            .arg("-cp")
            .arg(BBTOOLS_CP)
            .arg("jgi.KmerNormalize")
            .args(&shared)
            .output()
            .expect("run java bbnorm");
        assert_success("java bbnorm", &java_status);

        let rust_args = single_end_input_depth_bin_output_list_args(
            &[&input1, &input2],
            &rust_keep,
            output_key,
            &[&rust_bin_a, &rust_bin_b],
            &rust_hist,
            high_bin_depth,
        );
        let rust_status = Command::new(env!("CARGO_BIN_EXE_bbnorm-rs"))
            .args(&rust_args)
            .output()
            .expect("run bbnorm-rs");
        assert_success("bbnorm-rs", &rust_status);

        assert_same_file(&java_keep, &rust_keep);
        assert_same_file(&java_bin_a, &rust_bin_a);
        assert_same_file(&java_bin_b, &rust_bin_b);
        assert_same_file(&java_hist, &rust_hist);
    }
}

#[test]
fn representative_single_end_input_list_toss_fanout_matches_java_bbnorm() {
    require_file("vendor/BBTools-master/current/jgi/KmerNormalize.class");
    require_java();

    let dir = tempdir().expect("create tempdir");
    let input1 = dir.path().join("representative_list_a.fq");
    let input2 = dir.path().join("representative_list_b.fq");
    write_fastq_records(&input1, &[("a1", "ACGTACGT"), ("a2", "CCCCCCCC")]);
    write_fastq_records(&input2, &[("b1", "TTTTTTTT"), ("b2", "GGGGGGGG")]);

    let java_keep = dir.path().join("java.keep.fq");
    let java_toss_a = dir.path().join("java.toss.a.fq");
    let java_toss_b = dir.path().join("java.toss.b.fq");
    let java_hist = dir.path().join("java.toss_output_list.hist.tsv");
    let rust_keep = dir.path().join("rust.keep.fq");
    let rust_toss_a = dir.path().join("rust.toss.a.fq");
    let rust_toss_b = dir.path().join("rust.toss.b.fq");
    let rust_hist = dir.path().join("rust.toss_output_list.hist.tsv");

    let shared = single_end_input_toss_output_list_args(
        &[&input1, &input2],
        &java_keep,
        &[&java_toss_a, &java_toss_b],
        &java_hist,
    );
    let java_status = Command::new("java")
        .arg("-Xmx1g")
        .arg("-cp")
        .arg(BBTOOLS_CP)
        .arg("jgi.KmerNormalize")
        .args(&shared)
        .output()
        .expect("run java bbnorm");
    assert_success("java bbnorm", &java_status);

    let rust_args = single_end_input_toss_output_list_args(
        &[&input1, &input2],
        &rust_keep,
        &[&rust_toss_a, &rust_toss_b],
        &rust_hist,
    );
    let rust_status = Command::new(env!("CARGO_BIN_EXE_bbnorm-rs"))
        .args(&rust_args)
        .output()
        .expect("run bbnorm-rs");
    assert_success("bbnorm-rs", &rust_status);

    assert_same_file(&java_keep, &rust_keep);
    assert_same_file(&java_toss_a, &rust_toss_a);
    assert_same_file(&java_toss_b, &rust_toss_b);
    assert_same_file(&java_hist, &rust_hist);
}

#[test]
fn representative_paired_input_list_output_fanout_matches_java_bbnorm() {
    require_file("vendor/BBTools-master/current/jgi/KmerNormalize.class");
    require_java();

    let dir = tempdir().expect("create tempdir");
    let r1a = dir.path().join("representative_pair_list_a_1.fq");
    let r2a = dir.path().join("representative_pair_list_a_2.fq");
    let r1b = dir.path().join("representative_pair_list_b_1.fq");
    let r2b = dir.path().join("representative_pair_list_b_2.fq");
    write_fastq_records(&r1a, &[("a1/1", "AAA"), ("a2/1", "CCC")]);
    write_fastq_records(&r2a, &[("a1/2", "GGG"), ("a2/2", "TTT")]);
    write_fastq_records(&r1b, &[("b1/1", "ACG"), ("b2/1", "CGT")]);
    write_fastq_records(&r2b, &[("b1/2", "GTA"), ("b2/2", "TAC")]);

    let java_keep1_a = dir.path().join("java.keep1.a.fq");
    let java_keep1_b = dir.path().join("java.keep1.b.fq");
    let java_keep2_a = dir.path().join("java.keep2.a.fq");
    let java_keep2_b = dir.path().join("java.keep2.b.fq");
    let java_hist = dir.path().join("java.paired_output_list.hist.tsv");
    let rust_keep1_a = dir.path().join("rust.keep1.a.fq");
    let rust_keep1_b = dir.path().join("rust.keep1.b.fq");
    let rust_keep2_a = dir.path().join("rust.keep2.a.fq");
    let rust_keep2_b = dir.path().join("rust.keep2.b.fq");
    let rust_hist = dir.path().join("rust.paired_output_list.hist.tsv");

    let shared = paired_input_output_list_args(
        &[&r1a, &r1b],
        &[&r2a, &r2b],
        &[&java_keep1_a, &java_keep1_b],
        &[&java_keep2_a, &java_keep2_b],
        &java_hist,
    );
    let java_status = Command::new("java")
        .arg("-Xmx1g")
        .arg("-cp")
        .arg(BBTOOLS_CP)
        .arg("jgi.KmerNormalize")
        .args(&shared)
        .output()
        .expect("run java bbnorm");
    assert_success("java bbnorm", &java_status);

    let rust_args = paired_input_output_list_args(
        &[&r1a, &r1b],
        &[&r2a, &r2b],
        &[&rust_keep1_a, &rust_keep1_b],
        &[&rust_keep2_a, &rust_keep2_b],
        &rust_hist,
    );
    let rust_status = Command::new(env!("CARGO_BIN_EXE_bbnorm-rs"))
        .args(&rust_args)
        .output()
        .expect("run bbnorm-rs");
    assert_success("bbnorm-rs", &rust_status);

    assert_same_file(&java_keep1_a, &rust_keep1_a);
    assert_same_file(&java_keep1_b, &rust_keep1_b);
    assert_same_file(&java_keep2_a, &rust_keep2_a);
    assert_same_file(&java_keep2_b, &rust_keep2_b);
    assert_same_file(&java_hist, &rust_hist);
}

#[test]
fn representative_paired_input_list_low_bin_fanout_matches_java_bbnorm() {
    require_file("vendor/BBTools-master/current/jgi/KmerNormalize.class");
    require_java();

    let dir = tempdir().expect("create tempdir");
    let r1a = dir.path().join("representative_pair_list_a_1.fq");
    let r2a = dir.path().join("representative_pair_list_a_2.fq");
    let r1b = dir.path().join("representative_pair_list_b_1.fq");
    let r2b = dir.path().join("representative_pair_list_b_2.fq");
    write_fastq_records(&r1a, &[("a1/1", "AAA"), ("a2/1", "CCC")]);
    write_fastq_records(&r2a, &[("a1/2", "GGG"), ("a2/2", "TTT")]);
    write_fastq_records(&r1b, &[("b1/1", "ACG"), ("b2/1", "CGT")]);
    write_fastq_records(&r2b, &[("b1/2", "GTA"), ("b2/2", "TAC")]);

    let java_keep1 = dir.path().join("java.keep1.fq");
    let java_keep2 = dir.path().join("java.keep2.fq");
    let java_low1_a = dir.path().join("java.low1.a.fq");
    let java_low1_b = dir.path().join("java.low1.b.fq");
    let java_low2_a = dir.path().join("java.low2.a.fq");
    let java_low2_b = dir.path().join("java.low2.b.fq");
    let java_hist = dir.path().join("java.paired_low_output_list.hist.tsv");
    let rust_keep1 = dir.path().join("rust.keep1.fq");
    let rust_keep2 = dir.path().join("rust.keep2.fq");
    let rust_low1_a = dir.path().join("rust.low1.a.fq");
    let rust_low1_b = dir.path().join("rust.low1.b.fq");
    let rust_low2_a = dir.path().join("rust.low2.a.fq");
    let rust_low2_b = dir.path().join("rust.low2.b.fq");
    let rust_hist = dir.path().join("rust.paired_low_output_list.hist.tsv");

    let shared = paired_input_low_bin_output_list_args(
        &[&r1a, &r1b],
        &[&r2a, &r2b],
        &java_keep1,
        &java_keep2,
        &[&java_low1_a, &java_low1_b],
        &[&java_low2_a, &java_low2_b],
        &java_hist,
    );
    let java_status = Command::new("java")
        .arg("-Xmx1g")
        .arg("-cp")
        .arg(BBTOOLS_CP)
        .arg("jgi.KmerNormalize")
        .args(&shared)
        .output()
        .expect("run java bbnorm");
    assert_success("java bbnorm", &java_status);

    let rust_args = paired_input_low_bin_output_list_args(
        &[&r1a, &r1b],
        &[&r2a, &r2b],
        &rust_keep1,
        &rust_keep2,
        &[&rust_low1_a, &rust_low1_b],
        &[&rust_low2_a, &rust_low2_b],
        &rust_hist,
    );
    let rust_status = Command::new(env!("CARGO_BIN_EXE_bbnorm-rs"))
        .args(&rust_args)
        .output()
        .expect("run bbnorm-rs");
    assert_success("bbnorm-rs", &rust_status);

    assert_same_file(&java_keep1, &rust_keep1);
    assert_same_file(&java_keep2, &rust_keep2);
    assert_same_file(&java_low1_a, &rust_low1_a);
    assert_same_file(&java_low1_b, &rust_low1_b);
    assert_same_file(&java_low2_a, &rust_low2_a);
    assert_same_file(&java_low2_b, &rust_low2_b);
    assert_same_file(&java_hist, &rust_hist);
}

#[test]
fn representative_paired_input_list_mid_and_high_bin_fanout_match_java_bbnorm() {
    require_file("vendor/BBTools-master/current/jgi/KmerNormalize.class");
    require_java();

    let dir = tempdir().expect("create tempdir");
    let r1a = dir.path().join("representative_pair_list_a_1.fq");
    let r2a = dir.path().join("representative_pair_list_a_2.fq");
    let r1b = dir.path().join("representative_pair_list_b_1.fq");
    let r2b = dir.path().join("representative_pair_list_b_2.fq");
    write_fastq_records(&r1a, &[("a1/1", "AAA"), ("a2/1", "CCC")]);
    write_fastq_records(&r2a, &[("a1/2", "GGG"), ("a2/2", "TTT")]);
    write_fastq_records(&r1b, &[("b1/1", "ACG"), ("b2/1", "CGT")]);
    write_fastq_records(&r2b, &[("b1/2", "GTA"), ("b2/2", "TAC")]);

    for (label, output_key1, output_key2, high_bin_depth) in [
        ("mid", "outmid", "outmid2", "80"),
        ("high", "outhigh", "outhigh2", "0"),
    ] {
        let java_keep1 = dir.path().join(format!("java.{label}.keep1.fq"));
        let java_keep2 = dir.path().join(format!("java.{label}.keep2.fq"));
        let java_bin1_a = dir.path().join(format!("java.{label}1.a.fq"));
        let java_bin1_b = dir.path().join(format!("java.{label}1.b.fq"));
        let java_bin2_a = dir.path().join(format!("java.{label}2.a.fq"));
        let java_bin2_b = dir.path().join(format!("java.{label}2.b.fq"));
        let java_hist = dir
            .path()
            .join(format!("java.paired_{label}_output_list.hist.tsv"));
        let rust_keep1 = dir.path().join(format!("rust.{label}.keep1.fq"));
        let rust_keep2 = dir.path().join(format!("rust.{label}.keep2.fq"));
        let rust_bin1_a = dir.path().join(format!("rust.{label}1.a.fq"));
        let rust_bin1_b = dir.path().join(format!("rust.{label}1.b.fq"));
        let rust_bin2_a = dir.path().join(format!("rust.{label}2.a.fq"));
        let rust_bin2_b = dir.path().join(format!("rust.{label}2.b.fq"));
        let rust_hist = dir
            .path()
            .join(format!("rust.paired_{label}_output_list.hist.tsv"));

        let shared = paired_input_depth_bin_output_list_args(
            &[&r1a, &r1b],
            &[&r2a, &r2b],
            &java_keep1,
            &java_keep2,
            output_key1,
            output_key2,
            &[&java_bin1_a, &java_bin1_b],
            &[&java_bin2_a, &java_bin2_b],
            &java_hist,
            high_bin_depth,
        );
        let java_status = Command::new("java")
            .arg("-Xmx1g")
            .arg("-cp")
            .arg(BBTOOLS_CP)
            .arg("jgi.KmerNormalize")
            .args(&shared)
            .output()
            .expect("run java bbnorm");
        assert_success("java bbnorm", &java_status);

        let rust_args = paired_input_depth_bin_output_list_args(
            &[&r1a, &r1b],
            &[&r2a, &r2b],
            &rust_keep1,
            &rust_keep2,
            output_key1,
            output_key2,
            &[&rust_bin1_a, &rust_bin1_b],
            &[&rust_bin2_a, &rust_bin2_b],
            &rust_hist,
            high_bin_depth,
        );
        let rust_status = Command::new(env!("CARGO_BIN_EXE_bbnorm-rs"))
            .args(&rust_args)
            .output()
            .expect("run bbnorm-rs");
        assert_success("bbnorm-rs", &rust_status);

        assert_same_file(&java_keep1, &rust_keep1);
        assert_same_file(&java_keep2, &rust_keep2);
        assert_same_file(&java_bin1_a, &rust_bin1_a);
        assert_same_file(&java_bin1_b, &rust_bin1_b);
        assert_same_file(&java_bin2_a, &rust_bin2_a);
        assert_same_file(&java_bin2_b, &rust_bin2_b);
        assert_same_file(&java_hist, &rust_hist);
    }
}

#[test]
fn representative_paired_input_list_toss_fanout_matches_java_bbnorm() {
    require_file("vendor/BBTools-master/current/jgi/KmerNormalize.class");
    require_java();

    let dir = tempdir().expect("create tempdir");
    let r1a = dir.path().join("representative_pair_list_a_1.fq");
    let r2a = dir.path().join("representative_pair_list_a_2.fq");
    let r1b = dir.path().join("representative_pair_list_b_1.fq");
    let r2b = dir.path().join("representative_pair_list_b_2.fq");
    write_fastq_records(&r1a, &[("a1/1", "AAA"), ("a2/1", "CCC")]);
    write_fastq_records(&r2a, &[("a1/2", "GGG"), ("a2/2", "TTT")]);
    write_fastq_records(&r1b, &[("b1/1", "ACG"), ("b2/1", "CGT")]);
    write_fastq_records(&r2b, &[("b1/2", "GTA"), ("b2/2", "TAC")]);

    let java_keep1 = dir.path().join("java.keep1.fq");
    let java_keep2 = dir.path().join("java.keep2.fq");
    let java_toss1_a = dir.path().join("java.toss1.a.fq");
    let java_toss1_b = dir.path().join("java.toss1.b.fq");
    let java_toss2_a = dir.path().join("java.toss2.a.fq");
    let java_toss2_b = dir.path().join("java.toss2.b.fq");
    let java_hist = dir.path().join("java.paired_toss_output_list.hist.tsv");
    let rust_keep1 = dir.path().join("rust.keep1.fq");
    let rust_keep2 = dir.path().join("rust.keep2.fq");
    let rust_toss1_a = dir.path().join("rust.toss1.a.fq");
    let rust_toss1_b = dir.path().join("rust.toss1.b.fq");
    let rust_toss2_a = dir.path().join("rust.toss2.a.fq");
    let rust_toss2_b = dir.path().join("rust.toss2.b.fq");
    let rust_hist = dir.path().join("rust.paired_toss_output_list.hist.tsv");

    let shared = paired_input_toss_output_list_args(
        &[&r1a, &r1b],
        &[&r2a, &r2b],
        &java_keep1,
        &java_keep2,
        &[&java_toss1_a, &java_toss1_b],
        &[&java_toss2_a, &java_toss2_b],
        &java_hist,
    );
    let java_status = Command::new("java")
        .arg("-Xmx1g")
        .arg("-cp")
        .arg(BBTOOLS_CP)
        .arg("jgi.KmerNormalize")
        .args(&shared)
        .output()
        .expect("run java bbnorm");
    assert_success("java bbnorm", &java_status);

    let rust_args = paired_input_toss_output_list_args(
        &[&r1a, &r1b],
        &[&r2a, &r2b],
        &rust_keep1,
        &rust_keep2,
        &[&rust_toss1_a, &rust_toss1_b],
        &[&rust_toss2_a, &rust_toss2_b],
        &rust_hist,
    );
    let rust_status = Command::new(env!("CARGO_BIN_EXE_bbnorm-rs"))
        .args(&rust_args)
        .output()
        .expect("run bbnorm-rs");
    assert_success("bbnorm-rs", &rust_status);

    assert_same_file(&java_keep1, &rust_keep1);
    assert_same_file(&java_keep2, &rust_keep2);
    assert_same_file(&java_toss1_a, &rust_toss1_a);
    assert_same_file(&java_toss1_b, &rust_toss1_b);
    assert_same_file(&java_toss2_a, &rust_toss2_a);
    assert_same_file(&java_toss2_b, &rust_toss2_b);
    assert_same_file(&java_hist, &rust_hist);
}

#[test]
fn representative_paired_input_list_second_output_list_without_first_fanout_matches_java_bbnorm() {
    require_file("vendor/BBTools-master/current/jgi/KmerNormalize.class");
    require_java();

    let dir = tempdir().expect("create tempdir");
    let r1a = dir.path().join("representative_pair_list_a_1.fq");
    let r2a = dir.path().join("representative_pair_list_a_2.fq");
    let r1b = dir.path().join("representative_pair_list_b_1.fq");
    let r2b = dir.path().join("representative_pair_list_b_2.fq");
    write_fastq_records(&r1a, &[("a1/1", "AAA"), ("a2/1", "CCC")]);
    write_fastq_records(&r2a, &[("a1/2", "GGG"), ("a2/2", "TTT")]);
    write_fastq_records(&r1b, &[("b1/1", "ACG"), ("b2/1", "CGT")]);
    write_fastq_records(&r2b, &[("b1/2", "GTA"), ("b2/2", "TAC")]);

    let java_keep1 = dir.path().join("java.keep1.fq");
    let java_keep2_a = dir.path().join("java.keep2.a.fq");
    let java_keep2_b = dir.path().join("java.keep2.b.fq");
    let java_hist = dir.path().join("java.second_output_list.hist.tsv");
    let rust_keep1 = dir.path().join("rust.keep1.fq");
    let rust_keep2_a = dir.path().join("rust.keep2.a.fq");
    let rust_keep2_b = dir.path().join("rust.keep2.b.fq");
    let rust_hist = dir.path().join("rust.second_output_list.hist.tsv");

    let shared = paired_input_second_output_list_args(
        &[&r1a, &r1b],
        &[&r2a, &r2b],
        &java_keep1,
        &[&java_keep2_a, &java_keep2_b],
        &java_hist,
    );
    let java_status = Command::new("java")
        .arg("-Xmx1g")
        .arg("-cp")
        .arg(BBTOOLS_CP)
        .arg("jgi.KmerNormalize")
        .args(&shared)
        .output()
        .expect("run java bbnorm");
    assert_success("java bbnorm", &java_status);

    let rust_args = paired_input_second_output_list_args(
        &[&r1a, &r1b],
        &[&r2a, &r2b],
        &rust_keep1,
        &[&rust_keep2_a, &rust_keep2_b],
        &rust_hist,
    );
    let rust_status = Command::new(env!("CARGO_BIN_EXE_bbnorm-rs"))
        .args(&rust_args)
        .output()
        .expect("run bbnorm-rs");
    assert_success("bbnorm-rs", &rust_status);

    assert_same_file(&java_keep1, &rust_keep1);
    assert_same_file(&java_keep2_a, &rust_keep2_a);
    assert_eq!(
        java_keep2_b.exists(),
        rust_keep2_b.exists(),
        "second out2 list tail should be ignored when out is not a list"
    );
    assert_same_file(&java_hist, &rust_hist);
}

#[test]
fn representative_input_list_short_output_list_fails_lazily_like_java_bbnorm() {
    require_file("vendor/BBTools-master/current/jgi/KmerNormalize.class");
    require_java();

    let dir = tempdir().expect("create tempdir");
    let input1 = dir.path().join("representative_list_a.fq");
    let input2 = dir.path().join("representative_list_b.fq");
    let input3 = dir.path().join("representative_list_c.fq");
    write_fastq_records(&input1, &[("a1", "ACGT")]);
    write_fastq_records(&input2, &[("b1", "ACGT")]);
    write_fastq_records(&input3, &[("c1", "ACGT")]);

    let java_keep_a = dir.path().join("java.keep.a.fq");
    let java_keep_b = dir.path().join("java.keep.b.fq");
    let java_hist = dir.path().join("java.short_output_list.hist.tsv");
    let rust_keep_a = dir.path().join("rust.keep.a.fq");
    let rust_keep_b = dir.path().join("rust.keep.b.fq");
    let rust_hist = dir.path().join("rust.short_output_list.hist.tsv");

    let shared = single_end_input_output_list_args(
        &[&input1, &input2, &input3],
        &[&java_keep_a, &java_keep_b],
        &java_hist,
    );
    let java_status = Command::new("java")
        .arg("-Xmx1g")
        .arg("-cp")
        .arg(BBTOOLS_CP)
        .arg("jgi.KmerNormalize")
        .args(&shared)
        .output()
        .expect("run java bbnorm");
    assert!(
        !java_status.status.success(),
        "java bbnorm unexpectedly accepted too-short output list"
    );

    let rust_args = single_end_input_output_list_args(
        &[&input1, &input2, &input3],
        &[&rust_keep_a, &rust_keep_b],
        &rust_hist,
    );
    let rust_status = Command::new(env!("CARGO_BIN_EXE_bbnorm-rs"))
        .args(&rust_args)
        .output()
        .expect("run bbnorm-rs");
    assert!(
        !rust_status.status.success(),
        "bbnorm-rs unexpectedly accepted too-short output list"
    );

    assert_same_file(&java_keep_a, &rust_keep_a);
    assert_same_file(&java_keep_b, &rust_keep_b);
}

#[test]
fn representative_paired_input_list_short_second_output_list_fails_lazily_like_java_bbnorm() {
    require_file("vendor/BBTools-master/current/jgi/KmerNormalize.class");
    require_java();

    let dir = tempdir().expect("create tempdir");
    let r1a = dir.path().join("representative_pair_list_a_1.fq");
    let r2a = dir.path().join("representative_pair_list_a_2.fq");
    let r1b = dir.path().join("representative_pair_list_b_1.fq");
    let r2b = dir.path().join("representative_pair_list_b_2.fq");
    let r1c = dir.path().join("representative_pair_list_c_1.fq");
    let r2c = dir.path().join("representative_pair_list_c_2.fq");
    write_fastq_records(&r1a, &[("a1/1", "AAA")]);
    write_fastq_records(&r2a, &[("a1/2", "GGG")]);
    write_fastq_records(&r1b, &[("b1/1", "AAA")]);
    write_fastq_records(&r2b, &[("b1/2", "GGG")]);
    write_fastq_records(&r1c, &[("c1/1", "AAA")]);
    write_fastq_records(&r2c, &[("c1/2", "GGG")]);

    let java_keep1_a = dir.path().join("java.keep1.a.fq");
    let java_keep1_b = dir.path().join("java.keep1.b.fq");
    let java_keep1_c = dir.path().join("java.keep1.c.fq");
    let java_keep2_a = dir.path().join("java.keep2.a.fq");
    let java_keep2_b = dir.path().join("java.keep2.b.fq");
    let java_hist = dir.path().join("java.short_second_output_list.hist.tsv");
    let rust_keep1_a = dir.path().join("rust.keep1.a.fq");
    let rust_keep1_b = dir.path().join("rust.keep1.b.fq");
    let rust_keep1_c = dir.path().join("rust.keep1.c.fq");
    let rust_keep2_a = dir.path().join("rust.keep2.a.fq");
    let rust_keep2_b = dir.path().join("rust.keep2.b.fq");
    let rust_hist = dir.path().join("rust.short_second_output_list.hist.tsv");

    let shared = paired_input_output_list_args(
        &[&r1a, &r1b, &r1c],
        &[&r2a, &r2b, &r2c],
        &[&java_keep1_a, &java_keep1_b, &java_keep1_c],
        &[&java_keep2_a, &java_keep2_b],
        &java_hist,
    );
    let java_status = Command::new("java")
        .arg("-Xmx1g")
        .arg("-cp")
        .arg(BBTOOLS_CP)
        .arg("jgi.KmerNormalize")
        .args(&shared)
        .output()
        .expect("run java bbnorm");
    assert!(
        !java_status.status.success(),
        "java bbnorm unexpectedly accepted too-short out2 list"
    );

    let rust_args = paired_input_output_list_args(
        &[&r1a, &r1b, &r1c],
        &[&r2a, &r2b, &r2c],
        &[&rust_keep1_a, &rust_keep1_b, &rust_keep1_c],
        &[&rust_keep2_a, &rust_keep2_b],
        &rust_hist,
    );
    let rust_status = Command::new(env!("CARGO_BIN_EXE_bbnorm-rs"))
        .args(&rust_args)
        .output()
        .expect("run bbnorm-rs");
    assert!(
        !rust_status.status.success(),
        "bbnorm-rs unexpectedly accepted too-short out2 list"
    );

    assert_same_file(&java_keep1_a, &rust_keep1_a);
    assert_same_file(&java_keep1_b, &rust_keep1_b);
    assert_eq!(
        java_keep1_c.exists(),
        rust_keep1_c.exists(),
        "third first-output file should not be opened after out2 indexing fails"
    );
    assert_same_file(&java_keep2_a, &rust_keep2_a);
    assert_same_file(&java_keep2_b, &rust_keep2_b);
}

fn assert_min_length_single_end_matches_java(minlen: &str) {
    require_file(SAMPLE1);
    require_file("vendor/BBTools-master/current/jgi/KmerNormalize.class");
    require_java();

    let dir = tempdir().expect("create tempdir");
    let input = dir.path().join("real_duplicated_short_for_minlen.fq");
    write_duplicated_n_free_single_end_fixture(SAMPLE1, &input, 4);

    let java_keep = dir.path().join("java.keep.fq");
    let java_toss = dir.path().join("java.toss.fq");
    let rust_keep = dir.path().join("rust.keep.fq");
    let rust_toss = dir.path().join("rust.toss.fq");

    let shared = min_length_single_end_args(&input, &java_keep, &java_toss, minlen);
    let java_status = Command::new("java")
        .arg("-Xmx1g")
        .arg("-cp")
        .arg(BBTOOLS_CP)
        .arg("jgi.KmerNormalize")
        .args(&shared)
        .output()
        .expect("run java bbnorm");
    assert_success("java bbnorm", &java_status);

    let rust_args = min_length_single_end_args(&input, &rust_keep, &rust_toss, minlen);
    let rust_status = Command::new(env!("CARGO_BIN_EXE_bbnorm-rs"))
        .args(&rust_args)
        .output()
        .expect("run bbnorm-rs");
    assert_success("bbnorm-rs", &rust_status);

    assert_same_file(&java_keep, &rust_keep);
    assert_same_file(&java_toss, &rust_toss);
}

#[test]
fn real_phi_x_duplicated_high_depth_matches_java_bbnorm_for_long_kmer_hist() {
    require_file(SAMPLE1);
    require_file("vendor/BBTools-master/current/jgi/KmerNormalize.class");
    require_java();

    let dir = tempdir().expect("create tempdir");
    let input = dir.path().join("real_duplicated_high_depth.fq");
    write_duplicated_n_free_single_end_fixture(SAMPLE1, &input, 8);

    let java_keep = dir.path().join("java.keep.fq");
    let java_hist = dir.path().join("java.hist.tsv");
    let java_rhist = dir.path().join("java.rhist.tsv");
    let rust_keep = dir.path().join("rust.keep.fq");
    let rust_hist = dir.path().join("rust.hist.tsv");
    let rust_rhist = dir.path().join("rust.rhist.tsv");

    let shared = long_kmer_keepall_single_end_args(&input, &java_keep, &java_hist, &java_rhist);
    let java_status = Command::new("java")
        .arg("-Xmx1g")
        .arg("-cp")
        .arg(BBTOOLS_CP)
        .arg("jgi.KmerNormalize")
        .args(&shared)
        .output()
        .expect("run java bbnorm");
    assert_success("java bbnorm", &java_status);

    let rust_args = long_kmer_keepall_single_end_args(&input, &rust_keep, &rust_hist, &rust_rhist);
    let rust_status = Command::new(env!("CARGO_BIN_EXE_bbnorm-rs"))
        .args(&rust_args)
        .output()
        .expect("run bbnorm-rs");
    assert_success("bbnorm-rs", &rust_status);

    assert_same_file(&java_keep, &rust_keep);
    assert_same_file(&java_hist, &rust_hist);
    assert_same_file(&java_rhist, &rust_rhist);
}

#[test]
fn representative_long_kmer_spike_matches_java_bbnorm_for_fixspikes() {
    require_file("vendor/BBTools-master/current/jgi/KmerNormalize.class");
    require_java();

    let dir = tempdir().expect("create tempdir");
    let input = dir.path().join("long_kmer_spike.fq");
    write_long_kmer_spike_fixture(&input);

    let java_keep = dir.path().join("java.keep.fq");
    let java_hist = dir.path().join("java.hist.tsv");
    let java_rhist = dir.path().join("java.rhist.tsv");
    let rust_keep = dir.path().join("rust.keep.fq");
    let rust_hist = dir.path().join("rust.hist.tsv");
    let rust_rhist = dir.path().join("rust.rhist.tsv");

    let shared = long_kmer_keepall_single_end_args_with_extra(
        &input,
        &java_keep,
        &java_hist,
        &java_rhist,
        &["fixspikes=t"],
    );
    let java_status = Command::new("java")
        .arg("-Xmx1g")
        .arg("-cp")
        .arg(BBTOOLS_CP)
        .arg("jgi.KmerNormalize")
        .args(&shared)
        .output()
        .expect("run java bbnorm");
    assert_success("java bbnorm", &java_status);

    let rust_args = long_kmer_keepall_single_end_args_with_extra(
        &input,
        &rust_keep,
        &rust_hist,
        &rust_rhist,
        &["fixspikes=t"],
    );
    let rust_status = Command::new(env!("CARGO_BIN_EXE_bbnorm-rs"))
        .args(&rust_args)
        .output()
        .expect("run bbnorm-rs");
    assert_success("bbnorm-rs", &rust_status);

    assert_same_file(&java_keep, &rust_keep);
    assert_same_file(&java_hist, &rust_hist);
    assert_same_file(&java_rhist, &rust_rhist);
}

#[test]
fn representative_peak_output_matches_java_bbnorm() {
    require_file("vendor/BBTools-master/current/jgi/KmerNormalize.class");
    require_java();

    let dir = tempdir().expect("create tempdir");
    let input = dir.path().join("representative_peak.fq");
    write_representative_peak_fixture(&input);

    let java_keep = dir.path().join("java.keep.fq");
    let java_hist = dir.path().join("java.hist.tsv");
    let java_peaks = dir.path().join("java.peaks.tsv");
    let java_peaksout = dir.path().join("java.peaksout.tsv");
    let rust_keep = dir.path().join("rust.keep.fq");
    let rust_hist = dir.path().join("rust.hist.tsv");
    let rust_peaks = dir.path().join("rust.peaks.tsv");
    let rust_peaksout = dir.path().join("rust.peaksout.tsv");

    let shared =
        representative_peak_args(&input, &java_keep, &java_hist, &java_peaks, &java_peaksout);
    let java_status = Command::new("java")
        .arg("-Xmx1g")
        .arg("-cp")
        .arg(BBTOOLS_CP)
        .arg("jgi.KmerNormalize")
        .args(&shared)
        .output()
        .expect("run java bbnorm");
    assert_success("java bbnorm", &java_status);

    let rust_args =
        representative_peak_args(&input, &rust_keep, &rust_hist, &rust_peaks, &rust_peaksout);
    let rust_status = Command::new(env!("CARGO_BIN_EXE_bbnorm-rs"))
        .args(&rust_args)
        .output()
        .expect("run bbnorm-rs");
    assert_success("bbnorm-rs", &rust_status);

    assert_same_file(&java_keep, &rust_keep);
    assert_same_file(&java_hist, &rust_hist);
    assert_same_file(&java_peaks, &rust_peaks);
    assert_same_file(&java_peaksout, &rust_peaksout);
}

#[test]
fn representative_multi_peak_output_matches_java_bbnorm() {
    require_file("vendor/BBTools-master/current/jgi/KmerNormalize.class");
    require_java();

    let dir = tempdir().expect("create tempdir");
    let input = dir.path().join("representative_multi_peak.fq");
    write_representative_multi_peak_fixture(&input);

    let java_keep = dir.path().join("java.keep.fq");
    let java_hist = dir.path().join("java.hist.tsv");
    let java_peaks = dir.path().join("java.peaks.tsv");
    let java_peaksout = dir.path().join("java.peaksout.tsv");
    let rust_keep = dir.path().join("rust.keep.fq");
    let rust_hist = dir.path().join("rust.hist.tsv");
    let rust_peaks = dir.path().join("rust.peaks.tsv");
    let rust_peaksout = dir.path().join("rust.peaksout.tsv");

    let shared = multi_peak_args(&input, &java_keep, &java_hist, &java_peaks, &java_peaksout);
    let java_status = Command::new("java")
        .arg("-Xmx1g")
        .arg("-cp")
        .arg(BBTOOLS_CP)
        .arg("jgi.KmerNormalize")
        .args(&shared)
        .output()
        .expect("run java bbnorm");
    assert_success("java bbnorm", &java_status);

    let rust_args = multi_peak_args(&input, &rust_keep, &rust_hist, &rust_peaks, &rust_peaksout);
    let rust_status = Command::new(env!("CARGO_BIN_EXE_bbnorm-rs"))
        .args(&rust_args)
        .output()
        .expect("run bbnorm-rs");
    assert_success("bbnorm-rs", &rust_status);

    assert_same_file(&java_keep, &rust_keep);
    assert_same_file(&java_hist, &rust_hist);
    assert_same_file(&java_peaks, &rust_peaks);
    assert_same_file(&java_peaksout, &rust_peaksout);
}

#[test]
fn real_phi_x_mixed_depth_pair_matches_java_bbnorm_for_use_lower_depth_save_rare() {
    require_file(SAMPLE1);
    require_file(SAMPLE2);
    require_file("vendor/BBTools-master/current/jgi/KmerNormalize.class");
    require_java();

    let dir = tempdir().expect("create tempdir");
    let input1 = dir.path().join("mixed_depth.r1.fq");
    let input2 = dir.path().join("mixed_depth.r2.fq");
    write_mixed_depth_paired_fixture(SAMPLE1, SAMPLE2, &input1, &input2, 8);

    for use_lower_depth in [true, false] {
        let label = if use_lower_depth {
            "lower_depth"
        } else {
            "higher_depth"
        };
        let java_keep1 = dir.path().join(format!("java.{label}.keep1.fq"));
        let java_keep2 = dir.path().join(format!("java.{label}.keep2.fq"));
        let java_toss1 = dir.path().join(format!("java.{label}.toss1.fq"));
        let java_toss2 = dir.path().join(format!("java.{label}.toss2.fq"));
        let rust_keep1 = dir.path().join(format!("rust.{label}.keep1.fq"));
        let rust_keep2 = dir.path().join(format!("rust.{label}.keep2.fq"));
        let rust_toss1 = dir.path().join(format!("rust.{label}.toss1.fq"));
        let rust_toss2 = dir.path().join(format!("rust.{label}.toss2.fq"));

        let java_paths = NormalizePaths {
            keep1: &java_keep1,
            keep2: &java_keep2,
            toss1: &java_toss1,
            toss2: &java_toss2,
        };
        let shared = mixed_depth_save_rare_args(&input1, &input2, java_paths, use_lower_depth);
        let java_status = Command::new("java")
            .arg("-Xmx1g")
            .arg("-cp")
            .arg(BBTOOLS_CP)
            .arg("jgi.KmerNormalize")
            .args(&shared)
            .output()
            .expect("run java bbnorm");
        assert_success("java bbnorm", &java_status);

        let rust_paths = NormalizePaths {
            keep1: &rust_keep1,
            keep2: &rust_keep2,
            toss1: &rust_toss1,
            toss2: &rust_toss2,
        };
        let rust_args = mixed_depth_save_rare_args(&input1, &input2, rust_paths, use_lower_depth);
        let rust_status = Command::new(env!("CARGO_BIN_EXE_bbnorm-rs"))
            .args(&rust_args)
            .output()
            .expect("run bbnorm-rs");
        assert_success("bbnorm-rs", &rust_status);

        assert_same_file(&java_keep1, &rust_keep1);
        assert_same_file(&java_keep2, &rust_keep2);
        assert_same_file(&java_toss1, &rust_toss1);
        assert_same_file(&java_toss2, &rust_toss2);
    }
}

#[test]
fn real_phi_x_mixed_depth_pair_matches_java_bbnorm_for_percentile_save_rare() {
    require_file(SAMPLE1);
    require_file(SAMPLE2);
    require_file("vendor/BBTools-master/current/jgi/KmerNormalize.class");
    require_java();

    let dir = tempdir().expect("create tempdir");
    let input1 = dir.path().join("mixed_depth.r1.fq");
    let input2 = dir.path().join("mixed_depth.r2.fq");
    write_mixed_depth_paired_fixture(SAMPLE1, SAMPLE2, &input1, &input2, 8);

    let cases: [(&str, bool, &[&str]); 4] = [
        ("dp20", true, &["percentile=20"]),
        ("dp80", true, &["percentile=80"]),
        (
            "hdp50_ldp10",
            true,
            &["highdepthpercentile=50", "lowdepthpercentile=10"],
        ),
        (
            "mixed_percentiles_higher_depth",
            false,
            &[
                "percentile=80",
                "highdepthpercentile=50",
                "lowdepthpercentile=10",
            ],
        ),
    ];

    for (label, use_lower_depth, extra) in cases {
        let java_keep1 = dir.path().join(format!("java.{label}.keep1.fq"));
        let java_keep2 = dir.path().join(format!("java.{label}.keep2.fq"));
        let java_toss1 = dir.path().join(format!("java.{label}.toss1.fq"));
        let java_toss2 = dir.path().join(format!("java.{label}.toss2.fq"));
        let rust_keep1 = dir.path().join(format!("rust.{label}.keep1.fq"));
        let rust_keep2 = dir.path().join(format!("rust.{label}.keep2.fq"));
        let rust_toss1 = dir.path().join(format!("rust.{label}.toss1.fq"));
        let rust_toss2 = dir.path().join(format!("rust.{label}.toss2.fq"));

        let java_paths = NormalizePaths {
            keep1: &java_keep1,
            keep2: &java_keep2,
            toss1: &java_toss1,
            toss2: &java_toss2,
        };
        let shared = mixed_depth_save_rare_args_with_extra(
            &input1,
            &input2,
            java_paths,
            use_lower_depth,
            extra,
        );
        let java_status = Command::new("java")
            .arg("-Xmx1g")
            .arg("-cp")
            .arg(BBTOOLS_CP)
            .arg("jgi.KmerNormalize")
            .args(&shared)
            .output()
            .expect("run java bbnorm");
        assert_success("java bbnorm", &java_status);

        let rust_paths = NormalizePaths {
            keep1: &rust_keep1,
            keep2: &rust_keep2,
            toss1: &rust_toss1,
            toss2: &rust_toss2,
        };
        let rust_args = mixed_depth_save_rare_args_with_extra(
            &input1,
            &input2,
            rust_paths,
            use_lower_depth,
            extra,
        );
        let rust_status = Command::new(env!("CARGO_BIN_EXE_bbnorm-rs"))
            .args(&rust_args)
            .output()
            .expect("run bbnorm-rs");
        assert_success("bbnorm-rs", &rust_status);

        assert_same_file(&java_keep1, &rust_keep1);
        assert_same_file(&java_keep2, &rust_keep2);
        assert_same_file(&java_toss1, &rust_toss1);
        assert_same_file(&java_toss2, &rust_toss2);
    }
}

#[test]
fn real_phi_x_interleaved_pair_matches_java_bbnorm_for_single_pass_normalization() {
    require_file(SAMPLE1);
    require_file(SAMPLE2);
    require_file("vendor/BBTools-master/current/jgi/KmerNormalize.class");
    require_java();

    let dir = tempdir().expect("create tempdir");
    let interleaved = dir.path().join("interleaved.fq");
    write_interleaved_fixture(SAMPLE1, SAMPLE2, &interleaved);

    let java_keep = dir.path().join("java.keep.fq");
    let java_toss = dir.path().join("java.toss.fq");
    let rust_keep = dir.path().join("rust.keep.fq");
    let rust_toss = dir.path().join("rust.toss.fq");

    let shared = normalize_interleaved_args(&interleaved, &java_keep, &java_toss);
    let java_status = Command::new("java")
        .arg("-Xmx1g")
        .arg("-cp")
        .arg(BBTOOLS_CP)
        .arg("jgi.KmerNormalize")
        .args(&shared)
        .output()
        .expect("run java bbnorm");
    assert_success("java bbnorm", &java_status);

    let rust_args = normalize_interleaved_args(&interleaved, &rust_keep, &rust_toss);
    let rust_status = Command::new(env!("CARGO_BIN_EXE_bbnorm-rs"))
        .args(&rust_args)
        .output()
        .expect("run bbnorm-rs");
    assert_success("bbnorm-rs", &rust_status);

    assert_same_file(&java_keep, &rust_keep);
    assert_same_file(&java_toss, &rust_toss);
}

#[test]
fn real_phi_x_interleaved_pair_matches_java_bbnorm_for_renamed_keepall_output() {
    require_file(SAMPLE1);
    require_file(SAMPLE2);
    require_file("vendor/BBTools-master/current/jgi/KmerNormalize.class");
    require_java();

    let dir = tempdir().expect("create tempdir");
    let interleaved = dir.path().join("interleaved.fq");
    write_interleaved_fixture(SAMPLE1, SAMPLE2, &interleaved);

    let java_keep = dir.path().join("java.keep.fq");
    let rust_keep = dir.path().join("rust.keep.fq");

    let shared = rename_interleaved_keepall_args(&interleaved, &java_keep);
    let java_status = Command::new("java")
        .arg("-Xmx1g")
        .arg("-cp")
        .arg(BBTOOLS_CP)
        .arg("jgi.KmerNormalize")
        .args(&shared)
        .output()
        .expect("run java bbnorm");
    assert_success("java bbnorm", &java_status);

    let rust_args = rename_interleaved_keepall_args(&interleaved, &rust_keep);
    let rust_status = Command::new(env!("CARGO_BIN_EXE_bbnorm-rs"))
        .args(&rust_args)
        .output()
        .expect("run bbnorm-rs");
    assert_success("bbnorm-rs", &rust_status);

    assert_same_file(&java_keep, &rust_keep);
}

#[test]
fn real_phi_x_auto_interleaved_pair_matches_java_bbnorm_for_single_pass_normalization() {
    require_file(SAMPLE1);
    require_file(SAMPLE2);
    require_file("vendor/BBTools-master/current/jgi/KmerNormalize.class");
    require_java();

    let dir = tempdir().expect("create tempdir");
    let interleaved = dir.path().join("interleaved.fq");
    write_interleaved_fixture_with_pair_names(SAMPLE1, SAMPLE2, &interleaved);

    let java_keep = dir.path().join("java.keep.fq");
    let java_toss = dir.path().join("java.toss.fq");
    let rust_keep = dir.path().join("rust.keep.fq");
    let rust_toss = dir.path().join("rust.toss.fq");

    let shared = normalize_interleaved_auto_args(&interleaved, &java_keep, &java_toss);
    let java_status = Command::new("java")
        .arg("-Xmx1g")
        .arg("-cp")
        .arg(BBTOOLS_CP)
        .arg("jgi.KmerNormalize")
        .args(&shared)
        .output()
        .expect("run java bbnorm");
    assert_success("java bbnorm", &java_status);

    let rust_args = normalize_interleaved_auto_args(&interleaved, &rust_keep, &rust_toss);
    let rust_status = Command::new(env!("CARGO_BIN_EXE_bbnorm-rs"))
        .args(&rust_args)
        .output()
        .expect("run bbnorm-rs");
    assert_success("bbnorm-rs", &rust_status);

    assert_same_file(&java_keep, &rust_keep);
    assert_same_file(&java_toss, &rust_toss);
}

#[derive(Clone, Copy)]
struct OutputPaths<'a> {
    keep1: &'a Path,
    keep2: &'a Path,
    hist: &'a Path,
    rhist: &'a Path,
    low1: &'a Path,
    low2: &'a Path,
    mid1: &'a Path,
    mid2: &'a Path,
    high1: &'a Path,
    high2: &'a Path,
}

#[derive(Clone, Copy)]
struct NormalizePaths<'a> {
    keep1: &'a Path,
    keep2: &'a Path,
    toss1: &'a Path,
    toss2: &'a Path,
}

fn shared_args(paths: OutputPaths<'_>) -> Vec<OsString> {
    [
        format!("in={SAMPLE1}"),
        format!("in2={SAMPLE2}"),
        format!("out={}", paths.keep1.display()),
        format!("out2={}", paths.keep2.display()),
        format!("hist={}", paths.hist.display()),
        format!("rhist={}", paths.rhist.display()),
        format!("outlow={}", paths.low1.display()),
        format!("outlow2={}", paths.low2.display()),
        format!("outmid={}", paths.mid1.display()),
        format!("outmid2={}", paths.mid2.display()),
        format!("outhigh={}", paths.high1.display()),
        format!("outhigh2={}", paths.high2.display()),
        "passes=1".to_string(),
        "keepall=t".to_string(),
        "lowbindepth=1".to_string(),
        "highbindepth=1".to_string(),
        "k=31".to_string(),
        "minq=0".to_string(),
        "minprob=0".to_string(),
        "min=0".to_string(),
        "minkmers=1".to_string(),
        "target=999999999".to_string(),
        "max=999999999".to_string(),
        "threads=1".to_string(),
        "overwrite=t".to_string(),
        "bits=32".to_string(),
    ]
    .into_iter()
    .map(OsString::from)
    .collect()
}

fn shared_args_with_extra(paths: OutputPaths<'_>, extra: &[&str]) -> Vec<OsString> {
    let mut args = shared_args(paths);
    args.extend(extra.iter().map(OsString::from));
    args
}

fn hash_output_pattern_args(keep_hash: &Path, low_hash: &Path) -> Vec<OsString> {
    [
        format!("in={SAMPLE1}"),
        format!("in2={SAMPLE2}"),
        format!("out={}", keep_hash.display()),
        format!("outlow={}", low_hash.display()),
        "passes=1".to_string(),
        "keepall=t".to_string(),
        "lowbindepth=1".to_string(),
        "highbindepth=1".to_string(),
        "k=31".to_string(),
        "minq=0".to_string(),
        "minprob=0".to_string(),
        "min=0".to_string(),
        "minkmers=1".to_string(),
        "target=999999999".to_string(),
        "max=999999999".to_string(),
        "threads=1".to_string(),
        "overwrite=t".to_string(),
        "bits=32".to_string(),
    ]
    .into_iter()
    .map(OsString::from)
    .collect()
}

fn implicit_interleaved_paired_output_args(keep: &Path) -> Vec<OsString> {
    [
        format!("in={SAMPLE1}"),
        format!("in2={SAMPLE2}"),
        format!("out={}", keep.display()),
        "passes=1".to_string(),
        "keepall=t".to_string(),
        "k=31".to_string(),
        "minq=0".to_string(),
        "minprob=0".to_string(),
        "min=0".to_string(),
        "minkmers=1".to_string(),
        "target=999999999".to_string(),
        "max=999999999".to_string(),
        "threads=1".to_string(),
        "overwrite=t".to_string(),
        "bits=32".to_string(),
    ]
    .into_iter()
    .map(OsString::from)
    .collect()
}

fn outuncorrected_noecc_args(
    keep1: &Path,
    keep2: &Path,
    unc1: &Path,
    unc2: &Path,
) -> Vec<OsString> {
    [
        format!("in={SAMPLE1}"),
        format!("in2={SAMPLE2}"),
        format!("out={}", keep1.display()),
        format!("out2={}", keep2.display()),
        format!("outuncorrected={}", unc1.display()),
        format!("outuncorrected2={}", unc2.display()),
        "passes=1".to_string(),
        "keepall=t".to_string(),
        "k=31".to_string(),
        "minq=0".to_string(),
        "minprob=0".to_string(),
        "min=0".to_string(),
        "minkmers=1".to_string(),
        "target=999999999".to_string(),
        "max=999999999".to_string(),
        "threads=1".to_string(),
        "overwrite=t".to_string(),
        "bits=32".to_string(),
    ]
    .into_iter()
    .map(OsString::from)
    .collect()
}

fn multipass_markuncorrectable_args(
    input1: &Path,
    input2: &Path,
    keep1: &Path,
    keep2: &Path,
    unc1: &Path,
    unc2: &Path,
    stage_args: &[&str],
) -> Vec<OsString> {
    let mut args = vec![
        format!("in={}", input1.display()),
        format!("in2={}", input2.display()),
        format!("out={}", keep1.display()),
        format!("out2={}", keep2.display()),
        format!("outuncorrected={}", unc1.display()),
        format!("outuncorrected2={}", unc2.display()),
        "passes=2".to_string(),
        "keepall=t".to_string(),
        "ecco=f".to_string(),
        "eccmaxqual=0".to_string(),
        "markuncorrectableerrors=t".to_string(),
        "k=7".to_string(),
        "minq=0".to_string(),
        "minprob=0".to_string(),
        "min=0".to_string(),
        "minkmers=1".to_string(),
        "target=999999999".to_string(),
        "max=999999999".to_string(),
        "threads=1".to_string(),
        "overwrite=t".to_string(),
        "bits=32".to_string(),
    ];
    args.extend(stage_args.iter().map(|arg| (*arg).to_string()));
    args.into_iter().map(OsString::from).collect()
}

fn overlap_only_ecco_args(
    input1: &Path,
    input2: &Path,
    keep1: &Path,
    keep2: &Path,
    ecco: bool,
) -> Vec<OsString> {
    overlap_only_ecco_value_args(input1, input2, keep1, keep2, if ecco { "t" } else { "f" })
}

fn overlap_only_ecco_value_args(
    input1: &Path,
    input2: &Path,
    keep1: &Path,
    keep2: &Path,
    ecco: &str,
) -> Vec<OsString> {
    [
        format!("in={}", input1.display()),
        format!("in2={}", input2.display()),
        format!("out={}", keep1.display()),
        format!("out2={}", keep2.display()),
        "passes=1".to_string(),
        "keepall=t".to_string(),
        "ecc=t".to_string(),
        format!("ecco={ecco}"),
        "k=62".to_string(),
        "minq=0".to_string(),
        "minprob=0".to_string(),
        "min=0".to_string(),
        "minkmers=1".to_string(),
        "target=999999999".to_string(),
        "max=999999999".to_string(),
        "threads=1".to_string(),
        "overwrite=t".to_string(),
        "bits=32".to_string(),
    ]
    .into_iter()
    .map(OsString::from)
    .collect()
}

fn normalize_args(paths: NormalizePaths<'_>) -> Vec<OsString> {
    [
        format!("in={SAMPLE1}"),
        format!("in2={SAMPLE2}"),
        format!("out={}", paths.keep1.display()),
        format!("out2={}", paths.keep2.display()),
        format!("outt={}", paths.toss1.display()),
        format!("outt2={}", paths.toss2.display()),
        "passes=1".to_string(),
        "k=31".to_string(),
        "minq=0".to_string(),
        "minprob=0".to_string(),
        "min=0".to_string(),
        "minkmers=1".to_string(),
        "target=999999999".to_string(),
        "max=999999999".to_string(),
        "tossbadreads=t".to_string(),
        "threads=1".to_string(),
        "overwrite=t".to_string(),
        "bits=32".to_string(),
    ]
    .into_iter()
    .map(OsString::from)
    .collect()
}

fn normalize_args_without_default_toss(paths: NormalizePaths<'_>) -> Vec<OsString> {
    [
        format!("in={SAMPLE1}"),
        format!("in2={SAMPLE2}"),
        format!("out={}", paths.keep1.display()),
        format!("out2={}", paths.keep2.display()),
        format!("outt={}", paths.toss1.display()),
        format!("outt2={}", paths.toss2.display()),
        "passes=1".to_string(),
        "k=31".to_string(),
        "minq=0".to_string(),
        "minprob=0".to_string(),
        "min=0".to_string(),
        "minkmers=1".to_string(),
        "target=999999999".to_string(),
        "max=999999999".to_string(),
        "threads=1".to_string(),
        "overwrite=t".to_string(),
        "bits=32".to_string(),
    ]
    .into_iter()
    .map(OsString::from)
    .collect()
}

fn normalize_hash_input_args(paths: NormalizePaths<'_>) -> Vec<OsString> {
    [
        format!("in={SAMPLE_HASH}"),
        format!("out={}", paths.keep1.display()),
        format!("out2={}", paths.keep2.display()),
        format!("outt={}", paths.toss1.display()),
        format!("outt2={}", paths.toss2.display()),
        "passes=1".to_string(),
        "k=31".to_string(),
        "minq=0".to_string(),
        "minprob=0".to_string(),
        "min=0".to_string(),
        "minkmers=1".to_string(),
        "target=999999999".to_string(),
        "max=999999999".to_string(),
        "threads=1".to_string(),
        "overwrite=t".to_string(),
        "bits=32".to_string(),
    ]
    .into_iter()
    .map(OsString::from)
    .collect()
}

fn rename_paired_keepall_args(keep1: &Path, keep2: &Path) -> Vec<OsString> {
    [
        format!("in={SAMPLE1}"),
        format!("in2={SAMPLE2}"),
        format!("out={}", keep1.display()),
        format!("out2={}", keep2.display()),
        "passes=1".to_string(),
        "keepall=t".to_string(),
        "rename=t".to_string(),
        "k=31".to_string(),
        "minq=0".to_string(),
        "minprob=0".to_string(),
        "min=0".to_string(),
        "minkmers=1".to_string(),
        "target=999999999".to_string(),
        "max=999999999".to_string(),
        "threads=1".to_string(),
        "overwrite=t".to_string(),
        "bits=32".to_string(),
    ]
    .into_iter()
    .map(OsString::from)
    .collect()
}

fn normalize_args_with_extra(paths: NormalizePaths<'_>, extra: &[&str]) -> Vec<OsString> {
    let mut args = normalize_args(paths);
    args.extend(extra.iter().map(OsString::from));
    args
}

fn normalize_with_histograms_args(
    paths: NormalizePaths<'_>,
    hist: &Path,
    rhist: &Path,
    histout: &Path,
    rhistout: &Path,
    extra: &[&str],
) -> Vec<OsString> {
    let mut args = normalize_args(paths);
    args.extend([
        OsString::from(format!("hist={}", hist.display())),
        OsString::from(format!("rhist={}", rhist.display())),
        OsString::from(format!("histout={}", histout.display())),
        OsString::from(format!("rhistout={}", rhistout.display())),
    ]);
    args.extend(extra.iter().map(OsString::from));
    args
}

fn normalize_args_without_default_toss_with_extra(
    paths: NormalizePaths<'_>,
    extra: &[&str],
) -> Vec<OsString> {
    let mut args = normalize_args_without_default_toss(paths);
    args.extend(extra.iter().map(OsString::from));
    args
}

fn normalize_interleaved_args(input: &Path, keep: &Path, toss: &Path) -> Vec<OsString> {
    [
        format!("in={}", input.display()),
        format!("out={}", keep.display()),
        format!("outt={}", toss.display()),
        "interleaved=t".to_string(),
        "passes=1".to_string(),
        "k=31".to_string(),
        "minq=0".to_string(),
        "minprob=0".to_string(),
        "min=0".to_string(),
        "minkmers=1".to_string(),
        "target=1".to_string(),
        "max=1".to_string(),
        "threads=1".to_string(),
        "overwrite=t".to_string(),
        "bits=32".to_string(),
    ]
    .into_iter()
    .map(OsString::from)
    .collect()
}

fn normalize_interleaved_auto_args(input: &Path, keep: &Path, toss: &Path) -> Vec<OsString> {
    let mut args = normalize_interleaved_args(input, keep, toss);
    if let Some(interleaved) = args
        .iter_mut()
        .find(|arg| arg.to_string_lossy() == "interleaved=t")
    {
        *interleaved = OsString::from("interleaved=auto");
    }
    args
}

fn rename_interleaved_keepall_args(input: &Path, keep: &Path) -> Vec<OsString> {
    [
        format!("in={}", input.display()),
        format!("out={}", keep.display()),
        "interleaved=t".to_string(),
        "passes=1".to_string(),
        "keepall=t".to_string(),
        "rename=t".to_string(),
        "k=31".to_string(),
        "minq=0".to_string(),
        "minprob=0".to_string(),
        "min=0".to_string(),
        "minkmers=1".to_string(),
        "target=999999999".to_string(),
        "max=999999999".to_string(),
        "threads=1".to_string(),
        "overwrite=t".to_string(),
        "bits=32".to_string(),
    ]
    .into_iter()
    .map(OsString::from)
    .collect()
}

fn interleaved_quality_keepall_args(
    input: &Path,
    keep: &Path,
    hist: &Path,
    interleaved_arg: &str,
    extra: &[&str],
) -> Vec<OsString> {
    [
        format!("in={}", input.display()),
        format!("out={}", keep.display()),
        format!("hist={}", hist.display()),
        interleaved_arg.to_string(),
        "passes=1".to_string(),
        "keepall=t".to_string(),
        "k=3".to_string(),
        "minq=0".to_string(),
        "minprob=0".to_string(),
        "min=0".to_string(),
        "minkmers=1".to_string(),
        "target=999999999".to_string(),
        "max=999999999".to_string(),
        "threads=1".to_string(),
        "overwrite=t".to_string(),
        "bits=32".to_string(),
    ]
    .into_iter()
    .chain(extra.iter().map(|arg| arg.to_string()))
    .map(OsString::from)
    .collect()
}

fn discard_bad_only_single_end_args(input: &Path, keep: &Path, toss: &Path) -> Vec<OsString> {
    [
        format!("in={}", input.display()),
        format!("out={}", keep.display()),
        format!("outt={}", toss.display()),
        "passes=1".to_string(),
        "k=31".to_string(),
        "minq=0".to_string(),
        "minprob=0".to_string(),
        "min=0".to_string(),
        "minkmers=1".to_string(),
        "target=1".to_string(),
        "max=1".to_string(),
        "discardbadonly=t".to_string(),
        "threads=1".to_string(),
        "overwrite=t".to_string(),
        "bits=32".to_string(),
    ]
    .into_iter()
    .map(OsString::from)
    .collect()
}

fn min_length_single_end_args(
    input: &Path,
    keep: &Path,
    toss: &Path,
    minlen: &str,
) -> Vec<OsString> {
    [
        format!("in={}", input.display()),
        format!("out={}", keep.display()),
        format!("outt={}", toss.display()),
        "passes=1".to_string(),
        "k=31".to_string(),
        "minq=0".to_string(),
        "minprob=0".to_string(),
        "min=0".to_string(),
        "minkmers=1".to_string(),
        "target=999999999".to_string(),
        "max=999999999".to_string(),
        format!("minlen={minlen}"),
        "threads=1".to_string(),
        "overwrite=t".to_string(),
        "bits=32".to_string(),
    ]
    .into_iter()
    .map(OsString::from)
    .collect()
}

fn long_kmer_keepall_single_end_args(
    input: &Path,
    keep: &Path,
    hist: &Path,
    rhist: &Path,
) -> Vec<OsString> {
    long_kmer_keepall_single_end_args_with_extra(input, keep, hist, rhist, &[])
}

fn long_kmer_keepall_single_end_args_with_extra(
    input: &Path,
    keep: &Path,
    hist: &Path,
    rhist: &Path,
    extra: &[&str],
) -> Vec<OsString> {
    [
        format!("in={}", input.display()),
        format!("out={}", keep.display()),
        format!("hist={}", hist.display()),
        format!("rhist={}", rhist.display()),
        "passes=1".to_string(),
        "keepall=t".to_string(),
        "k=40".to_string(),
        "minq=0".to_string(),
        "minprob=0".to_string(),
        "min=0".to_string(),
        "minkmers=1".to_string(),
        "target=999999999".to_string(),
        "max=999999999".to_string(),
        "threads=1".to_string(),
        "overwrite=t".to_string(),
        "bits=32".to_string(),
    ]
    .into_iter()
    .chain(extra.iter().map(|arg| arg.to_string()))
    .map(OsString::from)
    .collect()
}

fn representative_peak_args(
    input: &Path,
    keep: &Path,
    hist: &Path,
    peaks: &Path,
    peaksout: &Path,
) -> Vec<OsString> {
    [
        format!("in={}", input.display()),
        format!("out={}", keep.display()),
        format!("hist={}", hist.display()),
        format!("peaks={}", peaks.display()),
        format!("peaksout={}", peaksout.display()),
        "passes=1".to_string(),
        "keepall=t".to_string(),
        "k=5".to_string(),
        "minq=0".to_string(),
        "minprob=0".to_string(),
        "min=0".to_string(),
        "minkmers=1".to_string(),
        "target=999999999".to_string(),
        "max=999999999".to_string(),
        "threads=1".to_string(),
        "overwrite=t".to_string(),
        "bits=32".to_string(),
        "minheight=1".to_string(),
        "minvolume=1".to_string(),
        "minwidth=1".to_string(),
        "minpeak=1".to_string(),
        "maxpeak=100".to_string(),
        "maxpeakcount=8".to_string(),
    ]
    .into_iter()
    .map(OsString::from)
    .collect()
}

fn multi_peak_args(
    input: &Path,
    keep: &Path,
    hist: &Path,
    peaks: &Path,
    peaksout: &Path,
) -> Vec<OsString> {
    [
        format!("in={}", input.display()),
        format!("out={}", keep.display()),
        format!("hist={}", hist.display()),
        format!("peaks={}", peaks.display()),
        format!("peaksout={}", peaksout.display()),
        "passes=1".to_string(),
        "keepall=t".to_string(),
        "k=11".to_string(),
        "minq=0".to_string(),
        "minprob=0".to_string(),
        "min=0".to_string(),
        "minkmers=1".to_string(),
        "target=999999999".to_string(),
        "max=999999999".to_string(),
        "threads=1".to_string(),
        "overwrite=t".to_string(),
        "bits=32".to_string(),
        "minheight=1".to_string(),
        "minvolume=1".to_string(),
        "minwidth=1".to_string(),
        "minpeak=1".to_string(),
        "maxpeak=100".to_string(),
        "maxpeakcount=8".to_string(),
    ]
    .into_iter()
    .map(OsString::from)
    .collect()
}

fn qtrim_keepall_single_end_args(
    input: &Path,
    keep: &Path,
    hist: &Path,
    rhist: &Path,
    qtrim: &str,
) -> Vec<OsString> {
    [
        format!("in={}", input.display()),
        format!("out={}", keep.display()),
        format!("hist={}", hist.display()),
        format!("rhist={}", rhist.display()),
        "passes=1".to_string(),
        "keepall=t".to_string(),
        qtrim.to_string(),
        "trimq=10".to_string(),
        "k=31".to_string(),
        "minq=0".to_string(),
        "minprob=0".to_string(),
        "min=0".to_string(),
        "minkmers=1".to_string(),
        "target=999999999".to_string(),
        "max=999999999".to_string(),
        "threads=1".to_string(),
        "overwrite=t".to_string(),
        "bits=32".to_string(),
    ]
    .into_iter()
    .map(OsString::from)
    .collect()
}

fn qout64_keepall_args(input: &Path, keep: &Path) -> Vec<OsString> {
    [
        format!("in={}", input.display()),
        format!("out={}", keep.display()),
        "passes=1".to_string(),
        "keepall=t".to_string(),
        "k=3".to_string(),
        "minq=0".to_string(),
        "minprob=0".to_string(),
        "min=0".to_string(),
        "minkmers=1".to_string(),
        "target=999999999".to_string(),
        "max=999999999".to_string(),
        "threads=1".to_string(),
        "overwrite=t".to_string(),
        "bits=32".to_string(),
        "qin=33".to_string(),
        "qout=64".to_string(),
    ]
    .into_iter()
    .map(OsString::from)
    .collect()
}

fn qin64_keepall_args(input: &Path, keep: &Path) -> Vec<OsString> {
    [
        format!("in={}", input.display()),
        format!("out={}", keep.display()),
        "passes=1".to_string(),
        "keepall=t".to_string(),
        "k=3".to_string(),
        "minq=0".to_string(),
        "minprob=0".to_string(),
        "min=0".to_string(),
        "minkmers=1".to_string(),
        "target=999999999".to_string(),
        "max=999999999".to_string(),
        "threads=1".to_string(),
        "overwrite=t".to_string(),
        "bits=32".to_string(),
        "qin=64".to_string(),
    ]
    .into_iter()
    .map(OsString::from)
    .collect()
}

fn quality_alias_keepall_args(input: &Path, keep: &Path, extra: &[&str]) -> Vec<OsString> {
    [
        format!("in={}", input.display()),
        format!("out={}", keep.display()),
        "passes=1".to_string(),
        "keepall=t".to_string(),
        "k=3".to_string(),
        "minq=0".to_string(),
        "minprob=0".to_string(),
        "min=0".to_string(),
        "minkmers=1".to_string(),
        "target=999999999".to_string(),
        "max=999999999".to_string(),
        "threads=1".to_string(),
        "overwrite=t".to_string(),
        "bits=32".to_string(),
    ]
    .into_iter()
    .chain(extra.iter().map(|arg| arg.to_string()))
    .map(OsString::from)
    .collect()
}

fn append_keepall_args(input: &Path, keep: &Path, hist: &Path) -> Vec<OsString> {
    [
        format!("in={}", input.display()),
        format!("out={}", keep.display()),
        format!("hist={}", hist.display()),
        "passes=1".to_string(),
        "keepall=t".to_string(),
        "k=3".to_string(),
        "minq=0".to_string(),
        "minprob=0".to_string(),
        "min=0".to_string(),
        "minkmers=1".to_string(),
        "target=999999999".to_string(),
        "max=999999999".to_string(),
        "threads=1".to_string(),
        "append=t".to_string(),
        "bits=32".to_string(),
    ]
    .into_iter()
    .map(OsString::from)
    .collect()
}

fn default_multipass_fallback_args(input: &Path, keep: &Path, hist: &Path) -> Vec<OsString> {
    [
        format!("in={}", input.display()),
        format!("out={}", keep.display()),
        format!("hist={}", hist.display()),
        "keepall=t".to_string(),
        "k=3".to_string(),
        "minq=0".to_string(),
        "minprob=0".to_string(),
        "min=0".to_string(),
        "minkmers=1".to_string(),
        "target=999999999".to_string(),
        "max=999999999".to_string(),
        "threads=1".to_string(),
        "bits=32".to_string(),
    ]
    .into_iter()
    .map(OsString::from)
    .collect()
}

fn sampling_guard_args(input: &Path, keep: &Path, option: &str) -> Vec<OsString> {
    [
        format!("in={}", input.display()),
        format!("out={}", keep.display()),
        "passes=1".to_string(),
        "keepall=t".to_string(),
        "k=3".to_string(),
        "minq=0".to_string(),
        "minprob=0".to_string(),
        "min=0".to_string(),
        "minkmers=1".to_string(),
        "target=999999999".to_string(),
        "max=999999999".to_string(),
        "threads=1".to_string(),
        "overwrite=t".to_string(),
        "bits=32".to_string(),
        option.to_string(),
    ]
    .into_iter()
    .map(OsString::from)
    .collect()
}

fn extra_input_keepall_args(
    input: &Path,
    extra: &[&Path],
    keep: &Path,
    hist: &Path,
) -> Vec<OsString> {
    extra_input_keepall_args_with_options(input, extra, keep, hist, &[])
}

fn extra_input_keepall_args_with_options(
    input: &Path,
    extra: &[&Path],
    keep: &Path,
    hist: &Path,
    options: &[&str],
) -> Vec<OsString> {
    let extra_arg = extra
        .iter()
        .map(|path| path.display().to_string())
        .collect::<Vec<_>>()
        .join(",");
    [
        format!("in={}", input.display()),
        format!("extra={extra_arg}"),
        format!("out={}", keep.display()),
        format!("hist={}", hist.display()),
        "passes=1".to_string(),
        "keepall=t".to_string(),
        "k=3".to_string(),
        "minq=0".to_string(),
        "minprob=0".to_string(),
        "min=0".to_string(),
        "minkmers=1".to_string(),
        "target=999999999".to_string(),
        "max=999999999".to_string(),
        "threads=1".to_string(),
        "overwrite=t".to_string(),
        "bits=32".to_string(),
    ]
    .into_iter()
    .chain(options.iter().map(|arg| arg.to_string()))
    .map(OsString::from)
    .collect()
}

fn single_end_input_list_args(inputs: &[&Path], keep: &Path, hist: &Path) -> Vec<OsString> {
    let input_arg = join_path_args(inputs);
    [
        format!("in={input_arg}"),
        format!("out={}", keep.display()),
        format!("hist={}", hist.display()),
        "passes=1".to_string(),
        "keepall=t".to_string(),
        "k=3".to_string(),
        "minq=0".to_string(),
        "minprob=0".to_string(),
        "min=0".to_string(),
        "minkmers=1".to_string(),
        "target=999999999".to_string(),
        "max=999999999".to_string(),
        "threads=1".to_string(),
        "overwrite=t".to_string(),
        "bits=32".to_string(),
        "reads=1".to_string(),
    ]
    .into_iter()
    .map(OsString::from)
    .collect()
}

fn single_end_input_output_list_args(
    inputs: &[&Path],
    keep_outputs: &[&Path],
    hist: &Path,
) -> Vec<OsString> {
    single_end_input_output_list_args_with_options(inputs, keep_outputs, hist, &[])
}

fn single_end_input_output_list_args_with_options(
    inputs: &[&Path],
    keep_outputs: &[&Path],
    hist: &Path,
    options: &[&str],
) -> Vec<OsString> {
    let input_arg = join_path_args(inputs);
    let output_arg = join_path_args(keep_outputs);
    [
        format!("in={input_arg}"),
        format!("out={output_arg}"),
        format!("hist={}", hist.display()),
        "passes=1".to_string(),
        "keepall=t".to_string(),
        "k=3".to_string(),
        "minq=0".to_string(),
        "minprob=0".to_string(),
        "min=0".to_string(),
        "minkmers=1".to_string(),
        "target=999999999".to_string(),
        "max=999999999".to_string(),
        "threads=1".to_string(),
        "overwrite=t".to_string(),
        "bits=32".to_string(),
        "reads=1".to_string(),
    ]
    .into_iter()
    .chain(options.iter().map(|arg| arg.to_string()))
    .map(OsString::from)
    .collect()
}

fn single_end_input_low_bin_output_list_args(
    inputs: &[&Path],
    keep: &Path,
    low_outputs: &[&Path],
    hist: &Path,
) -> Vec<OsString> {
    let input_arg = join_path_args(inputs);
    let low_arg = join_path_args(low_outputs);
    [
        format!("in={input_arg}"),
        format!("out={}", keep.display()),
        format!("outlow={low_arg}"),
        format!("hist={}", hist.display()),
        "passes=1".to_string(),
        "keepall=t".to_string(),
        "k=3".to_string(),
        "minq=0".to_string(),
        "minprob=0".to_string(),
        "min=0".to_string(),
        "minkmers=1".to_string(),
        "target=999999999".to_string(),
        "max=999999999".to_string(),
        "threads=1".to_string(),
        "overwrite=t".to_string(),
        "bits=32".to_string(),
        "reads=1".to_string(),
    ]
    .into_iter()
    .map(OsString::from)
    .collect()
}

fn single_end_input_depth_bin_output_list_args(
    inputs: &[&Path],
    keep: &Path,
    output_key: &str,
    bin_outputs: &[&Path],
    hist: &Path,
    high_bin_depth: &str,
) -> Vec<OsString> {
    let input_arg = join_path_args(inputs);
    let bin_arg = join_path_args(bin_outputs);
    [
        format!("in={input_arg}"),
        format!("out={}", keep.display()),
        format!("{output_key}={bin_arg}"),
        format!("hist={}", hist.display()),
        "passes=1".to_string(),
        "keepall=t".to_string(),
        "k=3".to_string(),
        "minq=0".to_string(),
        "minprob=0".to_string(),
        "min=0".to_string(),
        "minkmers=1".to_string(),
        "target=999999999".to_string(),
        "max=999999999".to_string(),
        "threads=1".to_string(),
        "overwrite=t".to_string(),
        "bits=32".to_string(),
        "reads=1".to_string(),
        "lowbindepth=0".to_string(),
        format!("highbindepth={high_bin_depth}"),
    ]
    .into_iter()
    .map(OsString::from)
    .collect()
}

fn single_end_input_toss_output_list_args(
    inputs: &[&Path],
    keep: &Path,
    toss_outputs: &[&Path],
    hist: &Path,
) -> Vec<OsString> {
    let input_arg = join_path_args(inputs);
    let toss_arg = join_path_args(toss_outputs);
    [
        format!("in={input_arg}"),
        format!("out={}", keep.display()),
        format!("outt={toss_arg}"),
        format!("hist={}", hist.display()),
        "passes=1".to_string(),
        "k=3".to_string(),
        "minq=0".to_string(),
        "minprob=0".to_string(),
        "min=0".to_string(),
        "minkmers=1".to_string(),
        "target=999999999".to_string(),
        "max=999999999".to_string(),
        "threads=1".to_string(),
        "overwrite=t".to_string(),
        "bits=32".to_string(),
        "reads=1".to_string(),
        "minlen=9".to_string(),
    ]
    .into_iter()
    .map(OsString::from)
    .collect()
}

fn paired_input_list_args(
    inputs1: &[&Path],
    inputs2: &[&Path],
    keep1: &Path,
    keep2: &Path,
    hist: &Path,
) -> Vec<OsString> {
    paired_input_list_args_with_options(inputs1, inputs2, keep1, keep2, hist, &[])
}

fn paired_input_list_args_with_options(
    inputs1: &[&Path],
    inputs2: &[&Path],
    keep1: &Path,
    keep2: &Path,
    hist: &Path,
    options: &[&str],
) -> Vec<OsString> {
    let input1_arg = join_path_args(inputs1);
    let input2_arg = join_path_args(inputs2);
    [
        format!("in={input1_arg}"),
        format!("in2={input2_arg}"),
        format!("out={}", keep1.display()),
        format!("out2={}", keep2.display()),
        format!("hist={}", hist.display()),
        "passes=1".to_string(),
        "keepall=t".to_string(),
        "k=3".to_string(),
        "minq=0".to_string(),
        "minprob=0".to_string(),
        "min=0".to_string(),
        "minkmers=1".to_string(),
        "target=999999999".to_string(),
        "max=999999999".to_string(),
        "threads=1".to_string(),
        "overwrite=t".to_string(),
        "bits=32".to_string(),
        "reads=1".to_string(),
    ]
    .into_iter()
    .chain(options.iter().map(|arg| arg.to_string()))
    .map(OsString::from)
    .collect()
}

fn paired_input_low_bin_output_list_args(
    inputs1: &[&Path],
    inputs2: &[&Path],
    keep1: &Path,
    keep2: &Path,
    low1_outputs: &[&Path],
    low2_outputs: &[&Path],
    hist: &Path,
) -> Vec<OsString> {
    let input1_arg = join_path_args(inputs1);
    let input2_arg = join_path_args(inputs2);
    let low1_arg = join_path_args(low1_outputs);
    let low2_arg = join_path_args(low2_outputs);
    [
        format!("in={input1_arg}"),
        format!("in2={input2_arg}"),
        format!("out={}", keep1.display()),
        format!("out2={}", keep2.display()),
        format!("outlow={low1_arg}"),
        format!("outlow2={low2_arg}"),
        format!("hist={}", hist.display()),
        "passes=1".to_string(),
        "keepall=t".to_string(),
        "k=3".to_string(),
        "minq=0".to_string(),
        "minprob=0".to_string(),
        "min=0".to_string(),
        "minkmers=1".to_string(),
        "target=999999999".to_string(),
        "max=999999999".to_string(),
        "threads=1".to_string(),
        "overwrite=t".to_string(),
        "bits=32".to_string(),
        "reads=1".to_string(),
    ]
    .into_iter()
    .map(OsString::from)
    .collect()
}

#[expect(
    clippy::too_many_arguments,
    reason = "test helper mirrors paired CLI streams"
)]
fn paired_input_depth_bin_output_list_args(
    inputs1: &[&Path],
    inputs2: &[&Path],
    keep1: &Path,
    keep2: &Path,
    output_key1: &str,
    output_key2: &str,
    bin1_outputs: &[&Path],
    bin2_outputs: &[&Path],
    hist: &Path,
    high_bin_depth: &str,
) -> Vec<OsString> {
    let input1_arg = join_path_args(inputs1);
    let input2_arg = join_path_args(inputs2);
    let bin1_arg = join_path_args(bin1_outputs);
    let bin2_arg = join_path_args(bin2_outputs);
    [
        format!("in={input1_arg}"),
        format!("in2={input2_arg}"),
        format!("out={}", keep1.display()),
        format!("out2={}", keep2.display()),
        format!("{output_key1}={bin1_arg}"),
        format!("{output_key2}={bin2_arg}"),
        format!("hist={}", hist.display()),
        "passes=1".to_string(),
        "keepall=t".to_string(),
        "k=3".to_string(),
        "minq=0".to_string(),
        "minprob=0".to_string(),
        "min=0".to_string(),
        "minkmers=1".to_string(),
        "target=999999999".to_string(),
        "max=999999999".to_string(),
        "threads=1".to_string(),
        "overwrite=t".to_string(),
        "bits=32".to_string(),
        "reads=1".to_string(),
        "lowbindepth=0".to_string(),
        format!("highbindepth={high_bin_depth}"),
    ]
    .into_iter()
    .map(OsString::from)
    .collect()
}

fn paired_input_toss_output_list_args(
    inputs1: &[&Path],
    inputs2: &[&Path],
    keep1: &Path,
    keep2: &Path,
    toss1_outputs: &[&Path],
    toss2_outputs: &[&Path],
    hist: &Path,
) -> Vec<OsString> {
    let input1_arg = join_path_args(inputs1);
    let input2_arg = join_path_args(inputs2);
    let toss1_arg = join_path_args(toss1_outputs);
    let toss2_arg = join_path_args(toss2_outputs);
    [
        format!("in={input1_arg}"),
        format!("in2={input2_arg}"),
        format!("out={}", keep1.display()),
        format!("out2={}", keep2.display()),
        format!("outt={toss1_arg}"),
        format!("outt2={toss2_arg}"),
        format!("hist={}", hist.display()),
        "passes=1".to_string(),
        "k=3".to_string(),
        "minq=0".to_string(),
        "minprob=0".to_string(),
        "min=0".to_string(),
        "minkmers=1".to_string(),
        "target=999999999".to_string(),
        "max=999999999".to_string(),
        "threads=1".to_string(),
        "overwrite=t".to_string(),
        "bits=32".to_string(),
        "reads=1".to_string(),
        "minlen=4".to_string(),
    ]
    .into_iter()
    .map(OsString::from)
    .collect()
}

fn paired_input_output_list_args(
    inputs1: &[&Path],
    inputs2: &[&Path],
    keep1_outputs: &[&Path],
    keep2_outputs: &[&Path],
    hist: &Path,
) -> Vec<OsString> {
    let input1_arg = join_path_args(inputs1);
    let input2_arg = join_path_args(inputs2);
    let output1_arg = join_path_args(keep1_outputs);
    let output2_arg = join_path_args(keep2_outputs);
    [
        format!("in={input1_arg}"),
        format!("in2={input2_arg}"),
        format!("out={output1_arg}"),
        format!("out2={output2_arg}"),
        format!("hist={}", hist.display()),
        "passes=1".to_string(),
        "keepall=t".to_string(),
        "k=3".to_string(),
        "minq=0".to_string(),
        "minprob=0".to_string(),
        "min=0".to_string(),
        "minkmers=1".to_string(),
        "target=999999999".to_string(),
        "max=999999999".to_string(),
        "threads=1".to_string(),
        "overwrite=t".to_string(),
        "bits=32".to_string(),
        "reads=1".to_string(),
    ]
    .into_iter()
    .map(OsString::from)
    .collect()
}

fn paired_input_second_output_list_args(
    inputs1: &[&Path],
    inputs2: &[&Path],
    keep1: &Path,
    keep2_outputs: &[&Path],
    hist: &Path,
) -> Vec<OsString> {
    let input1_arg = join_path_args(inputs1);
    let input2_arg = join_path_args(inputs2);
    let output2_arg = join_path_args(keep2_outputs);
    [
        format!("in={input1_arg}"),
        format!("in2={input2_arg}"),
        format!("out={}", keep1.display()),
        format!("out2={output2_arg}"),
        format!("hist={}", hist.display()),
        "passes=1".to_string(),
        "keepall=t".to_string(),
        "k=3".to_string(),
        "minq=0".to_string(),
        "minprob=0".to_string(),
        "min=0".to_string(),
        "minkmers=1".to_string(),
        "target=999999999".to_string(),
        "max=999999999".to_string(),
        "threads=1".to_string(),
        "overwrite=t".to_string(),
        "bits=32".to_string(),
        "reads=1".to_string(),
    ]
    .into_iter()
    .map(OsString::from)
    .collect()
}

fn join_path_args(paths: &[&Path]) -> String {
    paths
        .iter()
        .map(|path| path.display().to_string())
        .collect::<Vec<_>>()
        .join(",")
}

fn mixed_depth_save_rare_args(
    input1: &Path,
    input2: &Path,
    paths: NormalizePaths<'_>,
    use_lower_depth: bool,
) -> Vec<OsString> {
    mixed_depth_save_rare_args_with_extra(input1, input2, paths, use_lower_depth, &[])
}

fn mixed_depth_save_rare_args_with_extra(
    input1: &Path,
    input2: &Path,
    paths: NormalizePaths<'_>,
    use_lower_depth: bool,
    extra: &[&str],
) -> Vec<OsString> {
    [
        format!("in={}", input1.display()),
        format!("in2={}", input2.display()),
        format!("out={}", paths.keep1.display()),
        format!("out2={}", paths.keep2.display()),
        format!("outt={}", paths.toss1.display()),
        format!("outt2={}", paths.toss2.display()),
        "passes=1".to_string(),
        "k=31".to_string(),
        "minq=0".to_string(),
        "minprob=0".to_string(),
        "min=0".to_string(),
        "minkmers=1".to_string(),
        "target=1".to_string(),
        "max=999999999".to_string(),
        "tossbadreads=t".to_string(),
        "saverarereads=t".to_string(),
        "lowthresh=1".to_string(),
        "highthresh=1".to_string(),
        "errordetectratio=2".to_string(),
        format!("uselowerdepth={}", if use_lower_depth { "t" } else { "f" }),
        "threads=1".to_string(),
        "overwrite=t".to_string(),
        "bits=32".to_string(),
    ]
    .into_iter()
    .chain(extra.iter().map(|arg| arg.to_string()))
    .map(OsString::from)
    .collect()
}

fn assert_normalize_matches_java_with_extra(extra: &[&str]) {
    require_file(SAMPLE1);
    require_file(SAMPLE2);
    require_file("vendor/BBTools-master/current/jgi/KmerNormalize.class");
    require_java();

    let dir = tempdir().expect("create tempdir");
    let java_keep1 = dir.path().join("java.keep1.fq");
    let java_keep2 = dir.path().join("java.keep2.fq");
    let java_toss1 = dir.path().join("java.toss1.fq");
    let java_toss2 = dir.path().join("java.toss2.fq");
    let rust_keep1 = dir.path().join("rust.keep1.fq");
    let rust_keep2 = dir.path().join("rust.keep2.fq");
    let rust_toss1 = dir.path().join("rust.toss1.fq");
    let rust_toss2 = dir.path().join("rust.toss2.fq");

    let java_paths = NormalizePaths {
        keep1: &java_keep1,
        keep2: &java_keep2,
        toss1: &java_toss1,
        toss2: &java_toss2,
    };
    let shared = normalize_args_with_extra(java_paths, extra);
    let java_status = Command::new("java")
        .arg("-Xmx1g")
        .arg("-cp")
        .arg(BBTOOLS_CP)
        .arg("jgi.KmerNormalize")
        .args(&shared)
        .output()
        .expect("run java bbnorm");
    assert_success("java bbnorm", &java_status);

    let rust_paths = NormalizePaths {
        keep1: &rust_keep1,
        keep2: &rust_keep2,
        toss1: &rust_toss1,
        toss2: &rust_toss2,
    };
    let rust_args = normalize_args_with_extra(rust_paths, extra);
    let rust_status = Command::new(env!("CARGO_BIN_EXE_bbnorm-rs"))
        .args(&rust_args)
        .output()
        .expect("run bbnorm-rs");
    assert_success("bbnorm-rs", &rust_status);

    assert_same_file(&java_keep1, &rust_keep1);
    assert_same_file(&java_keep2, &rust_keep2);
    assert_same_file(&java_toss1, &rust_toss1);
    assert_same_file(&java_toss2, &rust_toss2);
}

fn assert_normalize_without_default_toss_matches_java_with_extra(extra: &[&str]) {
    require_file(SAMPLE1);
    require_file(SAMPLE2);
    require_file("vendor/BBTools-master/current/jgi/KmerNormalize.class");
    require_java();

    let dir = tempdir().expect("create tempdir");
    let java_keep1 = dir.path().join("java.keep1.fq");
    let java_keep2 = dir.path().join("java.keep2.fq");
    let java_toss1 = dir.path().join("java.toss1.fq");
    let java_toss2 = dir.path().join("java.toss2.fq");
    let rust_keep1 = dir.path().join("rust.keep1.fq");
    let rust_keep2 = dir.path().join("rust.keep2.fq");
    let rust_toss1 = dir.path().join("rust.toss1.fq");
    let rust_toss2 = dir.path().join("rust.toss2.fq");

    let java_paths = NormalizePaths {
        keep1: &java_keep1,
        keep2: &java_keep2,
        toss1: &java_toss1,
        toss2: &java_toss2,
    };
    let shared = normalize_args_without_default_toss_with_extra(java_paths, extra);
    let java_status = Command::new("java")
        .arg("-Xmx1g")
        .arg("-cp")
        .arg(BBTOOLS_CP)
        .arg("jgi.KmerNormalize")
        .args(&shared)
        .output()
        .expect("run java bbnorm");
    assert_success("java bbnorm", &java_status);

    let rust_paths = NormalizePaths {
        keep1: &rust_keep1,
        keep2: &rust_keep2,
        toss1: &rust_toss1,
        toss2: &rust_toss2,
    };
    let rust_args = normalize_args_without_default_toss_with_extra(rust_paths, extra);
    let rust_status = Command::new(env!("CARGO_BIN_EXE_bbnorm-rs"))
        .args(&rust_args)
        .output()
        .expect("run bbnorm-rs");
    assert_success("bbnorm-rs", &rust_status);

    assert_same_file(&java_keep1, &rust_keep1);
    assert_same_file(&java_keep2, &rust_keep2);
    assert_same_file(&java_toss1, &rust_toss1);
    assert_same_file(&java_toss2, &rust_toss2);
}

fn write_interleaved_fixture(r1_path: &str, r2_path: &str, output: &Path) {
    let file1 = File::open(r1_path).expect("open r1 fixture");
    let file2 = File::open(r2_path).expect("open r2 fixture");
    let mut reader1 = BufReader::new(MultiGzDecoder::new(file1));
    let mut reader2 = BufReader::new(MultiGzDecoder::new(file2));
    let mut writer = File::create(output).expect("create interleaved fixture");

    loop {
        let r1 = read_fastq_record(&mut reader1);
        let r2 = read_fastq_record(&mut reader2);
        match (r1, r2) {
            (Some(record1), Some(record2)) => {
                for line in record1.iter().chain(record2.iter()) {
                    writer
                        .write_all(line.as_bytes())
                        .expect("write interleaved record");
                }
            }
            (None, None) => break,
            _ => panic!("paired fixtures have different record counts"),
        }
    }
}

fn write_interleaved_fixture_with_pair_names(r1_path: &str, r2_path: &str, output: &Path) {
    let file1 = File::open(r1_path).expect("open r1 fixture");
    let file2 = File::open(r2_path).expect("open r2 fixture");
    let mut reader1 = BufReader::new(MultiGzDecoder::new(file1));
    let mut reader2 = BufReader::new(MultiGzDecoder::new(file2));
    let mut writer = File::create(output).expect("create interleaved fixture");
    let mut pair_index = 0usize;

    loop {
        let mut r1 = read_fastq_record(&mut reader1);
        let mut r2 = read_fastq_record(&mut reader2);
        match (&mut r1, &mut r2) {
            (Some(record1), Some(record2)) => {
                record1[0] = format!("@pair{pair_index}/1\n");
                record2[0] = format!("@pair{pair_index}/2\n");
                for line in record1.iter().chain(record2.iter()) {
                    writer
                        .write_all(line.as_bytes())
                        .expect("write interleaved record");
                }
                pair_index += 1;
            }
            (None, None) => break,
            _ => panic!("paired fixtures have different record counts"),
        }
    }
}

fn write_duplicated_n_free_single_end_fixture(input_path: &str, output: &Path, copies: usize) {
    let file = File::open(input_path).expect("open single-end fixture");
    let mut reader = BufReader::new(MultiGzDecoder::new(file));
    let mut writer = File::create(output).expect("create duplicated fixture");

    while let Some(record) = read_fastq_record(&mut reader) {
        if record[1].bytes().any(|base| base == b'N' || base == b'n') {
            continue;
        }
        for copy in 0..copies {
            writeln!(writer, "@real_dup_{copy}").expect("write duplicate header");
            writer
                .write_all(record[1].as_bytes())
                .expect("write duplicate bases");
            writer
                .write_all(b"+\n")
                .expect("write duplicate quality header");
            writer
                .write_all(record[3].as_bytes())
                .expect("write duplicate qualities");
        }
        return;
    }

    panic!("no N-free read found in fixture: {input_path}");
}

fn write_long_kmer_spike_fixture(output: &Path) {
    let center = "ACGTTGCAAGTCGATCGTAGCTAGGATCCGATGCTAGTCA";
    assert_eq!(center.len(), 40);

    let mut writer = File::create(output).expect("create long-kmer spike fixture");
    let target = format!("A{center}C");
    writeln!(writer, "@target").expect("write target header");
    writeln!(writer, "{target}").expect("write target bases");
    writeln!(writer, "+").expect("write target quality header");
    writeln!(writer, "{}", "I".repeat(target.len())).expect("write target qualities");

    for copy in 0..3 {
        writeln!(writer, "@center_dup_{copy}").expect("write duplicate header");
        writeln!(writer, "{center}").expect("write duplicate bases");
        writeln!(writer, "+").expect("write duplicate quality header");
        writeln!(writer, "{}", "I".repeat(center.len())).expect("write duplicate qualities");
    }
}

fn write_representative_peak_fixture(output: &Path) {
    let mut writer = File::create(output).expect("create representative peak fixture");
    let seq = "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT";
    let qual = "I".repeat(seq.len());
    for idx in 0..20 {
        writeln!(writer, "@rep{idx}\n{seq}\n+\n{qual}").expect("write representative peak read");
    }
}

fn write_representative_multi_peak_fixture(output: &Path) {
    let mut writer = File::create(output).expect("create representative multi-peak fixture");
    let fixtures = [
        (
            "low",
            "GCTAAAGACAATTACATAACATACACGTCAGCACGAAACTTGTTGGCCCAGTGTGAATCGCTTAAGGGTTAAGTAAGTGT",
            10,
        ),
        (
            "high",
            "GATGCATACGCCTTTACTTGCTGTGTCCACCCCATCGGACTGGCATTTTTATTACACTCAGAAACAGAACTCGGGTAATT",
            25,
        ),
    ];
    for (label, seq, copies) in fixtures {
        let qual = "I".repeat(seq.len());
        for idx in 0..copies {
            writeln!(writer, "@{label}_{idx}\n{seq}\n+\n{qual}")
                .expect("write representative multi-peak read");
        }
    }
}

fn write_qtrim_right_fixture(input_path: &str, output: &Path, copies: usize, low_tail: usize) {
    let file = File::open(input_path).expect("open single-end fixture");
    let mut reader = BufReader::new(MultiGzDecoder::new(file));
    let mut writer = File::create(output).expect("create qtrim fixture");

    while let Some(record) = read_fastq_record(&mut reader) {
        let bases = record[1].trim_end();
        if bases.bytes().any(|base| base == b'N' || base == b'n') || bases.len() <= low_tail + 31 {
            continue;
        }
        let qualities = format!(
            "{}{}\n",
            "I".repeat(bases.len() - low_tail),
            "!".repeat(low_tail)
        );
        for copy in 0..copies {
            writeln!(writer, "@real_qtrim_{copy}").expect("write qtrim header");
            writer
                .write_all(record[1].as_bytes())
                .expect("write qtrim bases");
            writer.write_all(b"+\n").expect("write qtrim plus");
            writer
                .write_all(qualities.as_bytes())
                .expect("write qtrim qualities");
        }
        return;
    }

    panic!("no long N-free read found in fixture: {input_path}");
}

fn write_qout64_fixture(output: &Path) {
    let mut writer = File::create(output).expect("create qout64 fixture");
    writeln!(writer, "@qout64\nACGTNNACGT\n+\n!#I~IIIIII").expect("write qout64 read");
}

fn write_qin64_fixture(output: &Path) {
    let mut writer = File::create(output).expect("create qin64 fixture");
    writeln!(writer, "@qin64\nACGTNNACGT\n+\n@Bh|hhhhhh").expect("write qin64 read");
}

fn write_fastq_fixture(output: &Path, id: &str, bases: &str) {
    let mut writer = File::create(output).expect("create FASTQ fixture");
    writeln!(writer, "@{id}").expect("write FASTQ header");
    writeln!(writer, "{bases}").expect("write FASTQ bases");
    writeln!(writer, "+").expect("write FASTQ plus");
    writeln!(writer, "{}", "I".repeat(bases.len())).expect("write FASTQ qualities");
}

fn write_fastq_records(output: &Path, records: &[(&str, &str)]) {
    let mut writer = File::create(output).expect("create FASTQ records fixture");
    for (id, bases) in records {
        writeln!(writer, "@{id}").expect("write FASTQ header");
        writeln!(writer, "{bases}").expect("write FASTQ bases");
        writeln!(writer, "+").expect("write FASTQ plus");
        writeln!(writer, "{}", "I".repeat(bases.len())).expect("write FASTQ qualities");
    }
}

fn write_quality_auto_pair_fixtures(r1: &Path, r2: &Path, interleaved: &Path) {
    let r1_records = [
        ("pair0/1", "ACGTNNACGT", "!#I~IIIIII"),
        ("pair1/1", "TGCATGCATG", "IIIIIIIIII"),
    ];
    let r2_records = [
        ("pair0/2", "TGCANNACGT", "@Bh|hhhhhh"),
        ("pair1/2", "ACGTACGTAC", "hhhhhhhhhh"),
    ];
    write_fastq_records_with_qualities(r1, &r1_records);
    write_fastq_records_with_qualities(r2, &r2_records);

    let mut writer = File::create(interleaved).expect("create qauto interleaved fixture");
    for (&record1, &record2) in r1_records.iter().zip(r2_records.iter()) {
        write_fastq_record_with_quality(&mut writer, record1);
        write_fastq_record_with_quality(&mut writer, record2);
    }
}

fn write_multipass_uncorrectable_paired_fixture(r1: &Path, r2: &Path) {
    let clean1 = "ACGTTGCATGTCAGTACCGTAACGTTGCA";
    let clean2 = "TGCAACGTTACGGTACTGACATGCAACGT";
    let mutant1 = "ACGTTGCATGTCAGAACCGTAACGTTGCA";
    let quality1 = "I".repeat(clean1.len());
    let quality2 = "I".repeat(clean2.len());

    let mut r1_records = Vec::new();
    let mut r2_records = Vec::new();
    for i in 0..30 {
        r1_records.push((
            format!("clean{:03}/1", i + 1),
            clean1.to_string(),
            quality1.clone(),
        ));
        r2_records.push((
            format!("clean{:03}/2", i + 1),
            clean2.to_string(),
            quality2.clone(),
        ));
    }
    r1_records.push(("uncorrectable/1".to_string(), mutant1.to_string(), quality1));
    r2_records.push(("uncorrectable/2".to_string(), clean2.to_string(), quality2));

    let r1_records_ref: Vec<(&str, &str, &str)> = r1_records
        .iter()
        .map(|(id, bases, qualities)| (id.as_str(), bases.as_str(), qualities.as_str()))
        .collect();
    let r2_records_ref: Vec<(&str, &str, &str)> = r2_records
        .iter()
        .map(|(id, bases, qualities)| (id.as_str(), bases.as_str(), qualities.as_str()))
        .collect();
    write_fastq_records_with_qualities(r1, &r1_records_ref);
    write_fastq_records_with_qualities(r2, &r2_records_ref);
}

fn write_overlap_only_ecco_fixture(r1: &Path, r2: &Path) {
    let read1 = "ACGTTGCATGTCAGTAACGTTGCATGTCAGTAACGTTGCA";
    let read2_mutant = "TGCAACGTTACTGACATGCACCGTTACTGACATGCAACGT";
    let qual1 = "IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII";
    let qual2 = "IIIIIIIIIIIIIIIIIIII!IIIIIIIIIIIIIIIIIII";
    let r1_records = [
        ("overlap1/1", read1, qual1),
        ("overlap2/1", read1, qual1),
        ("overlap3/1", read1, qual1),
        ("overlap4/1", read1, qual1),
        ("overlap5/1", read1, qual1),
    ];
    let r2_records = [
        ("overlap1/2", read2_mutant, qual2),
        ("overlap2/2", read2_mutant, qual2),
        ("overlap3/2", read2_mutant, qual2),
        ("overlap4/2", read2_mutant, qual2),
        ("overlap5/2", read2_mutant, qual2),
    ];
    write_fastq_records_with_qualities(r1, &r1_records);
    write_fastq_records_with_qualities(r2, &r2_records);
}

fn write_overlap_only_ecco_high_entropy_fixture(r1: &Path, r2: &Path) {
    let read1 = "TTAGTTGTGCCGCAGCGAAGTAGTGCTTGAAATATGCGAC";
    let read2_mutant = "GTCGCATATTTCAAGCACTAATTCGCTGCGGCACAACTAA";
    let qual1 = "IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII";
    let qual2 = "IIIIIIIIIIIIIIIIIIII!IIIIIIIIIIIIIIIIIII";
    let r1_records = [
        ("overlap1/1", read1, qual1),
        ("overlap2/1", read1, qual1),
        ("overlap3/1", read1, qual1),
        ("overlap4/1", read1, qual1),
        ("overlap5/1", read1, qual1),
    ];
    let r2_records = [
        ("overlap1/2", read2_mutant, qual2),
        ("overlap2/2", read2_mutant, qual2),
        ("overlap3/2", read2_mutant, qual2),
        ("overlap4/2", read2_mutant, qual2),
        ("overlap5/2", read2_mutant, qual2),
    ];
    write_fastq_records_with_qualities(r1, &r1_records);
    write_fastq_records_with_qualities(r2, &r2_records);
}

fn write_overlap_only_ecco_high_confidence_fixture(r1: &Path, r2: &Path) {
    let read1 = "TTAGTTGTGCCGCAGCGAAGTAGTGCTTGAAATATGCGAC";
    let read2_mutant = "GTCGCATATTTCAAGCACTAATTCGCTGCGGCACAACTAA";
    let qual1 = "IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII";
    let qual2 = "IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII";
    let r1_records = [
        ("overlap1/1", read1, qual1),
        ("overlap2/1", read1, qual1),
        ("overlap3/1", read1, qual1),
        ("overlap4/1", read1, qual1),
        ("overlap5/1", read1, qual1),
    ];
    let r2_records = [
        ("overlap1/2", read2_mutant, qual2),
        ("overlap2/2", read2_mutant, qual2),
        ("overlap3/2", read2_mutant, qual2),
        ("overlap4/2", read2_mutant, qual2),
        ("overlap5/2", read2_mutant, qual2),
    ];
    write_fastq_records_with_qualities(r1, &r1_records);
    write_fastq_records_with_qualities(r2, &r2_records);
}

fn write_overlap_only_ecco_quality_fixture(r1: &Path, r2: &Path, mismatch_quality: u8) {
    let read1 = "TTAGTTGTGCCGCAGCGAAGTAGTGCTTGAAATATGCGAC";
    let read2_mutant = "GTCGCATATTTCAAGCACTAATTCGCTGCGGCACAACTAA";
    let qual1 = "IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII";
    let mut qual2 = vec![b'I'; 40];
    qual2[20] = mismatch_quality.saturating_add(33);
    let qual2 = String::from_utf8(qual2).expect("quality fixture utf8");
    let r1_records = [
        ("overlap1/1", read1, qual1),
        ("overlap2/1", read1, qual1),
        ("overlap3/1", read1, qual1),
        ("overlap4/1", read1, qual1),
        ("overlap5/1", read1, qual1),
    ];
    let r2_records = [
        ("overlap1/2", read2_mutant, qual2.as_str()),
        ("overlap2/2", read2_mutant, qual2.as_str()),
        ("overlap3/2", read2_mutant, qual2.as_str()),
        ("overlap4/2", read2_mutant, qual2.as_str()),
        ("overlap5/2", read2_mutant, qual2.as_str()),
    ];
    write_fastq_records_with_qualities(r1, &r1_records);
    write_fastq_records_with_qualities(r2, &r2_records);
}

fn write_overlap_only_ecco_quality_weighted_multimismatch_fixture(r1: &Path, r2: &Path) {
    let read1 = "CAGTAACCAATGCCTGTTGAGATGCCAGACGCGTAACCAAAA";
    let read2_mutant = "TTTTGCTAACGCGTCTGGCATCTCAACAGGCATTGGTTAC";
    let qual1 = "IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII";
    let qual2 = "IIIII!I'IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII";
    let r1_records = [
        ("overlap1/1", read1, qual1),
        ("overlap2/1", read1, qual1),
        ("overlap3/1", read1, qual1),
        ("overlap4/1", read1, qual1),
        ("overlap5/1", read1, qual1),
    ];
    let r2_records = [
        ("overlap1/2", read2_mutant, qual2),
        ("overlap2/2", read2_mutant, qual2),
        ("overlap3/2", read2_mutant, qual2),
        ("overlap4/2", read2_mutant, qual2),
        ("overlap5/2", read2_mutant, qual2),
    ];
    write_fastq_records_with_qualities(r1, &r1_records);
    write_fastq_records_with_qualities(r2, &r2_records);
}

fn write_fastq_records_with_qualities(output: &Path, records: &[(&str, &str, &str)]) {
    let mut writer = File::create(output).expect("create FASTQ quality fixture");
    for &record in records {
        write_fastq_record_with_quality(&mut writer, record);
    }
}

fn write_fastq_record_with_quality(writer: &mut File, record: (&str, &str, &str)) {
    let (id, bases, qualities) = record;
    writeln!(writer, "@{id}").expect("write FASTQ header");
    writeln!(writer, "{bases}").expect("write FASTQ bases");
    writeln!(writer, "+").expect("write FASTQ plus");
    writeln!(writer, "{qualities}").expect("write FASTQ qualities");
}

fn write_fasta_fixture(output: &Path, records: &[(&str, &str)]) {
    let mut writer = File::create(output).expect("create FASTA fixture");
    for (id, bases) in records {
        writeln!(writer, ">{id}").expect("write FASTA header");
        writeln!(writer, "{bases}").expect("write FASTA bases");
    }
}

fn write_mixed_depth_paired_fixture(
    r1_path: &str,
    r2_path: &str,
    output1: &Path,
    output2: &Path,
    pairs: usize,
) {
    let r1_template = first_n_free_fastq_records(r1_path, 1)
        .into_iter()
        .next()
        .expect("one N-free r1 record");
    let r2_records = first_n_free_fastq_records(r2_path, pairs);
    let mut writer1 = File::create(output1).expect("create mixed-depth r1 fixture");
    let mut writer2 = File::create(output2).expect("create mixed-depth r2 fixture");
    let mutated_r1_bases = mutate_middle_base(&r1_template[1]);

    for (pair, r2) in r2_records.iter().enumerate().take(pairs) {
        writeln!(writer1, "@mixed_r1_{pair}").expect("write r1 header");
        writer1
            .write_all(if pair == 0 {
                mutated_r1_bases.as_bytes()
            } else {
                r1_template[1].as_bytes()
            })
            .expect("write r1 bases");
        writer1.write_all(b"+\n").expect("write r1 plus");
        writer1
            .write_all(r1_template[3].as_bytes())
            .expect("write r1 qualities");

        writeln!(writer2, "@mixed_r2_{pair}").expect("write r2 header");
        writer2.write_all(r2[1].as_bytes()).expect("write r2 bases");
        writer2.write_all(b"+\n").expect("write r2 plus");
        writer2
            .write_all(r2[3].as_bytes())
            .expect("write r2 qualities");
    }
}

fn first_n_free_fastq_records(path: &str, count: usize) -> Vec<Vec<String>> {
    let file = File::open(path).expect("open FASTQ fixture");
    let mut reader = BufReader::new(MultiGzDecoder::new(file));
    let mut records = Vec::with_capacity(count);
    while records.len() < count {
        let Some(record) = read_fastq_record(&mut reader) else {
            break;
        };
        if !record[1].bytes().any(|base| base == b'N' || base == b'n')
            && record[1].trim_end().len() >= 80
        {
            records.push(record);
        }
    }
    assert_eq!(
        records.len(),
        count,
        "not enough N-free fixture reads in {path}"
    );
    records
}

fn mutate_middle_base(bases_line: &str) -> String {
    let mut bases = bases_line.trim_end().as_bytes().to_vec();
    let middle = bases.len() / 2;
    bases[middle] = match bases[middle] {
        b'A' | b'a' => b'C',
        b'C' | b'c' => b'G',
        b'G' | b'g' => b'T',
        b'T' | b't' | b'U' | b'u' => b'A',
        _ => b'A',
    };
    let mut mutated = String::from_utf8(bases).expect("fixture bases are UTF-8");
    mutated.push('\n');
    mutated
}

fn read_fastq_record<R: BufRead>(reader: &mut R) -> Option<Vec<String>> {
    let mut record = Vec::with_capacity(4);
    for _ in 0..4 {
        let mut line = String::new();
        let bytes = reader.read_line(&mut line).expect("read FASTQ line");
        if bytes == 0 {
            assert!(record.is_empty(), "truncated FASTQ record in fixture");
            return None;
        }
        record.push(line);
    }
    Some(record)
}

fn require_file(path: &str) {
    assert!(Path::new(path).exists(), "missing required fixture: {path}");
}

fn require_java() {
    let status = Command::new("java")
        .arg("-version")
        .output()
        .expect("java must be installed for real parity test");
    assert_success("java -version", &status);
}

fn assert_success(label: &str, output: &std::process::Output) {
    assert!(
        output.status.success(),
        "{label} failed\nstatus: {}\nstdout:\n{}\nstderr:\n{}",
        output.status,
        String::from_utf8_lossy(&output.stdout),
        String::from_utf8_lossy(&output.stderr)
    );
}

fn assert_same_file(left: &PathBuf, right: &PathBuf) {
    let left_bytes = fs::read(left).unwrap_or_else(|err| panic!("read {}: {err}", left.display()));
    let right_bytes =
        fs::read(right).unwrap_or_else(|err| panic!("read {}: {err}", right.display()));
    if left_bytes != right_bytes {
        let first_diff = left_bytes
            .iter()
            .zip(&right_bytes)
            .position(|(left, right)| left != right)
            .unwrap_or_else(|| left_bytes.len().min(right_bytes.len()));
        panic!(
            "files differ: {} vs {}; lengths {} vs {}; first differing byte offset {}",
            left.display(),
            right.display(),
            left_bytes.len(),
            right_bytes.len(),
            first_diff
        );
    }
}

fn fastq_record_count(path: &Path) -> usize {
    let bytes = fs::read(path).unwrap_or_else(|err| panic!("read {}: {err}", path.display()));
    let line_count = String::from_utf8(bytes)
        .unwrap_or_else(|err| panic!("decode {} as UTF-8: {err}", path.display()))
        .lines()
        .count();
    assert_eq!(
        line_count % 4,
        0,
        "FASTQ line count must be divisible by 4 in {}",
        path.display()
    );
    line_count / 4
}

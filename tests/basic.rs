use bbnorm_rs::{parse_args, run};
use std::ffi::OsString;
use std::fs;
use std::path::PathBuf;
use tempfile::tempdir;

#[test]
fn writes_histograms_and_keeps_all_fastq_records() {
    let dir = tempdir().unwrap();
    let input = dir.path().join("reads.fq");
    let keep = dir.path().join("keep.fq");
    let hist = dir.path().join("hist.tsv");
    let histout = dir.path().join("histout.tsv");
    let rhist = dir.path().join("rhist.tsv");

    fs::write(
        &input,
        b"@r1\nAAAAAA\n+\nIIIIII\n@r2\nAAAAAA\n+\nIIIIII\n@r3\nCCCCCC\n+\nIIIIII\n",
    )
    .unwrap();

    let args = [
        format!("in={}", input.display()),
        format!("out={}", keep.display()),
        format!("hist={}", hist.display()),
        format!("histout={}", histout.display()),
        format!("rhist={}", rhist.display()),
        "k=3".to_string(),
        "minq=0".to_string(),
        "minprob=0".to_string(),
        "minkmers=1".to_string(),
        "min=1".to_string(),
        "target=1".to_string(),
        "max=1".to_string(),
        "keepall=t".to_string(),
        "histlen=10".to_string(),
        "overwrite=t".to_string(),
        "passes=1".to_string(),
    ];
    let config = parse_args(args.into_iter().map(OsString::from)).unwrap();
    let summary = run(&config).unwrap();

    assert_eq!(summary.reads_in, 3);
    assert_eq!(summary.reads_kept, 3);
    assert_eq!(summary.reads_tossed, 0);
    assert_eq!(summary.unique_kmers_in, 2);
    assert_eq!(summary.unique_kmers_out, Some(2));

    let keep_text = fs::read_to_string(&keep).unwrap();
    assert_eq!(keep_text.matches('@').count(), 3);

    let hist_text = fs::read_to_string(&hist).unwrap();
    assert!(hist_text.contains("#Depth\tRaw_Count\tUnique_Kmers"));
    assert!(hist_text.contains("1\t4\t4"));
    assert!(hist_text.contains("2\t8\t4"));
    assert_eq!(hist_text, fs::read_to_string(&histout).unwrap());

    let rhist_text = fs::read_to_string(&rhist).unwrap();
    assert_eq!(rhist_text, "#Depth\tReads\tBases\n1\t1\t6\n2\t2\t12\n");
}

#[test]
fn multipass_toss_output_includes_intermediate_tosses() {
    let dir = tempdir().unwrap();
    let input = dir.path().join("reads.fq");
    let keep = dir.path().join("keep.fq");
    let toss = dir.path().join("toss.fq");
    let mut records = String::new();
    for i in 0..20 {
        records.push_str(&format!(
            "@r{i}\nACGTACGTACGTACGTACGT\n+\nIIIIIIIIIIIIIIIIIIII\n"
        ));
    }
    fs::write(&input, records).unwrap();

    let args = [
        format!("in={}", input.display()),
        format!("out={}", keep.display()),
        format!("outt={}", toss.display()),
        "passes=2".to_string(),
        "k=7".to_string(),
        "minq=0".to_string(),
        "minprob=0".to_string(),
        "min=0".to_string(),
        "minkmers=1".to_string(),
        "target=1".to_string(),
        "max=1".to_string(),
        "overwrite=t".to_string(),
    ];
    let config = parse_args(args.into_iter().map(OsString::from)).unwrap();
    run(&config).unwrap();

    let keep_records = fs::read_to_string(&keep).unwrap().lines().count() / 4;
    let toss_records = fs::read_to_string(&toss).unwrap().lines().count() / 4;
    assert!(toss_records > 0);
    assert_eq!(keep_records + toss_records, 20);
}

#[test]
fn parses_peak_calling_options() {
    let config = parse_args(
        [
            "in=reads.fq",
            "passes=1",
            "peaks=peaks.tsv",
            "peaksout=peaksout.tsv",
            "minheight=1",
            "minvolume=1",
            "minwidth=1",
            "minpeak=1",
            "maxpeak=100",
            "ploidy=1",
            "maxpeakcount=8",
        ]
        .into_iter()
        .map(OsString::from),
    )
    .unwrap();

    assert_eq!(config.peaks_in.unwrap(), PathBuf::from("peaks.tsv"));
    assert_eq!(config.peaks_out.unwrap(), PathBuf::from("peaksout.tsv"));
    assert_eq!(config.peak_min_height, 1);
    assert_eq!(config.peak_min_volume, 1);
    assert_eq!(config.peak_min_width, 1);
    assert_eq!(config.peak_min_peak, 1);
    assert_eq!(config.peak_max_peak, 100);
    assert_eq!(config.peak_ploidy, 1);
    assert_eq!(config.peak_max_count, 8);
}

#[test]
fn parses_outuncorrected_output_paths() {
    let config = parse_args(
        [
            "in=reads.fq",
            "passes=1",
            "in2=reads2.fq",
            "outuncorrected=unc1.fq",
            "outuncorrected2=unc2.fq",
        ]
        .into_iter()
        .map(OsString::from),
    )
    .unwrap();

    assert_eq!(config.out_uncorrected1.unwrap(), PathBuf::from("unc1.fq"));
    assert_eq!(config.out_uncorrected2.unwrap(), PathBuf::from("unc2.fq"));
}

#[test]
fn countup_allows_side_output_streams() {
    let config = parse_args(
        [
            "in=reads.fq",
            "in2=reads2.fq",
            "passes=1",
            "countup=t",
            "outlow=low1.fq",
            "outlow2=low2.fq",
            "outmid=mid1.fq",
            "outmid2=mid2.fq",
            "outhigh=high1.fq",
            "outhigh2=high2.fq",
            "outuncorrected=unc1.fq",
            "outuncorrected2=unc2.fq",
        ]
        .into_iter()
        .map(OsString::from),
    )
    .unwrap();

    assert!(config.count_up);
    assert_eq!(config.out_low1.unwrap(), PathBuf::from("low1.fq"));
    assert_eq!(config.out_uncorrected2.unwrap(), PathBuf::from("unc2.fq"));
}

#[test]
fn countup_writes_histout_and_peaksout() {
    let dir = tempdir().unwrap();
    let input = dir.path().join("reads.fq");
    let histout = dir.path().join("histout.tsv");
    let peaksout = dir.path().join("peaksout.tsv");
    fs::write(
        &input,
        b"@r1\nACGTACGTAC\n+\nIIIIIIIIII\n@r2\nACGTACGTAC\n+\nIIIIIIIIII\n@r3\nTTTTACGTAC\n+\nIIIIIIIIII\n",
    )
    .unwrap();

    let args = [
        format!("in={}", input.display()),
        format!("histout={}", histout.display()),
        format!("peaksout={}", peaksout.display()),
        "countup=t".to_string(),
        "k=3".to_string(),
        "minq=0".to_string(),
        "minprob=0".to_string(),
        "minkmers=1".to_string(),
        "min=1".to_string(),
        "target=2".to_string(),
        "max=3".to_string(),
        "histlen=20".to_string(),
        "overwrite=t".to_string(),
    ];
    let config = parse_args(args.into_iter().map(OsString::from)).unwrap();
    let summary = run(&config).unwrap();

    assert!(summary.reads_in > 0);
    assert!(summary.reads_in <= 3);
    assert!(summary.unique_kmers_out.is_some());
    assert!(
        fs::read_to_string(&histout)
            .unwrap()
            .contains("#Depth\tRaw_Count\tUnique_Kmers")
    );
    assert!(
        fs::read_to_string(&peaksout)
            .unwrap()
            .contains("#start\tcenter\tstop\tmax\tvolume")
    );
}

#[test]
fn multipass_allows_depth_bin_side_outputs() {
    let config = parse_args(
        [
            "in=reads.fq",
            "in2=reads2.fq",
            "passes=2",
            "outlow=low1.fq",
            "outlow2=low2.fq",
            "outmid=mid1.fq",
            "outmid2=mid2.fq",
            "outhigh=high1.fq",
            "outhigh2=high2.fq",
        ]
        .into_iter()
        .map(OsString::from),
    )
    .unwrap();

    assert_eq!(config.passes, 2);
    assert_eq!(config.out_low1.unwrap(), PathBuf::from("low1.fq"));
    assert_eq!(config.out_mid2.unwrap(), PathBuf::from("mid2.fq"));
    assert_eq!(config.out_high2.unwrap(), PathBuf::from("high2.fq"));
}

#[test]
fn parses_callpeaks_short_peak_aliases() {
    let config = parse_args(
        [
            "in=reads.fq",
            "passes=1",
            "h=1",
            "v=2",
            "w=3",
            "minp=4",
            "maxp=100",
        ]
        .into_iter()
        .map(OsString::from),
    )
    .unwrap();

    assert_eq!(config.peak_min_height, 1);
    assert_eq!(config.peak_min_volume, 2);
    assert_eq!(config.peak_min_width, 3);
    assert_eq!(config.peak_min_peak, 4);
    assert_eq!(config.peak_max_peak, 100);
}

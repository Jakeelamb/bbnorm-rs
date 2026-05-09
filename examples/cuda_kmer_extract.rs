use anyhow::{Context, Result, bail};
use bbnorm_rs::Config;
use bbnorm_rs::kmer::{KmerKey, for_each_kmer_for_record};
use bbnorm_rs::seqio::{QualitySettings, SequenceReader, SequenceSettings};
use std::env;
use std::fs::File;
use std::io::{BufWriter, Write, stderr};
use std::path::PathBuf;
use std::time::Instant;

#[derive(Debug)]
struct Args {
    r1: PathBuf,
    r2: Option<PathBuf>,
    out: PathBuf,
    reads: u64,
    k: usize,
    gzip_threads: Option<usize>,
}

fn parse_args() -> Result<Args> {
    let mut r1 = None;
    let mut r2 = None;
    let mut out = None;
    let mut reads = None;
    let mut k = 31usize;
    let mut gzip_threads = Some(4usize);

    let mut args = env::args().skip(1);
    while let Some(arg) = args.next() {
        match arg.as_str() {
            "--r1" => r1 = Some(PathBuf::from(value(&mut args, "--r1")?)),
            "--r2" => r2 = Some(PathBuf::from(value(&mut args, "--r2")?)),
            "--out" => out = Some(PathBuf::from(value(&mut args, "--out")?)),
            "--reads" => reads = Some(value(&mut args, "--reads")?.parse()?),
            "--k" => k = value(&mut args, "--k")?.parse()?,
            "--gzip-threads" => {
                let parsed: usize = value(&mut args, "--gzip-threads")?.parse()?;
                gzip_threads = (parsed > 0).then_some(parsed);
            }
            "--help" | "-h" => {
                println!(
                    "usage: cargo run --release --example cuda_kmer_extract -- --r1 R1 [--r2 R2] --out kmers.u64 [--reads N] [--k 31] [--gzip-threads 4]"
                );
                std::process::exit(0);
            }
            other => bail!("unknown argument: {other}"),
        }
    }

    let r1 = r1.context("--r1 is required")?;
    let out = out.context("--out is required")?;
    Ok(Args {
        r1,
        r2,
        out,
        reads: reads.unwrap_or(u64::MAX),
        k,
        gzip_threads,
    })
}

fn value(args: &mut impl Iterator<Item = String>, name: &str) -> Result<String> {
    args.next()
        .with_context(|| format!("{name} requires a value"))
}

fn main() -> Result<()> {
    let args = parse_args()?;
    if args.k == 0 || args.k > 31 {
        bail!("--k must be in 1..31 for the short-kmer CUDA probe");
    }

    let settings = SequenceSettings {
        qualities: QualitySettings {
            input_offset: 33,
            min_called: 2,
            max_called: 50,
            change_quality: true,
        },
        ..SequenceSettings::default()
    };
    let config = Config {
        k: args.k,
        ..Config::default()
    };

    let started = Instant::now();
    let mut reader1 =
        SequenceReader::from_path_with_gzip_threads(&args.r1, settings, args.gzip_threads)
            .with_context(|| format!("open R1 {}", args.r1.display()))?;
    let mut reader2 = args
        .r2
        .as_ref()
        .map(|path| {
            SequenceReader::from_path_with_gzip_threads(path, settings, args.gzip_threads)
                .with_context(|| format!("open R2 {}", path.display()))
        })
        .transpose()?;
    let out_is_stdout = args.out == PathBuf::from("-");
    let mut writer: Box<dyn Write> = if out_is_stdout {
        Box::new(BufWriter::new(std::io::stdout()))
    } else {
        let file =
            File::create(&args.out).with_context(|| format!("create {}", args.out.display()))?;
        Box::new(BufWriter::new(file))
    };

    let mut reads = 0u64;
    let mut kmers = 0u64;
    while reads < args.reads {
        let Some(record1) = reader1.next_record()? else {
            break;
        };
        kmers += write_record_kmers(&record1, &config, &mut writer)?;
        if let Some(reader2) = &mut reader2 {
            let record2 = reader2
                .next_record()?
                .with_context(|| "R2 has fewer records than R1")?;
            kmers += write_record_kmers(&record2, &config, &mut writer)?;
        }
        reads += 1;
    }
    writer.flush()?;

    let extract_seconds = started.elapsed().as_secs_f64();
    if out_is_stdout {
        let mut stats = stderr().lock();
        writeln!(stats, "extractor\trust")?;
        writeln!(stats, "reads\t{reads}")?;
        writeln!(stats, "extracted_kmers\t{kmers}")?;
        writeln!(stats, "extract_seconds\t{extract_seconds:.6}")?;
    } else {
        println!("extractor\trust");
        println!("reads\t{reads}");
        println!("extracted_kmers\t{kmers}");
        println!("extract_seconds\t{extract_seconds:.6}");
    }
    Ok(())
}

fn write_record_kmers(
    record: &bbnorm_rs::seqio::SequenceRecord,
    config: &Config,
    writer: &mut impl Write,
) -> Result<u64> {
    let mut written = 0u64;
    let mut error = None;
    for_each_kmer_for_record(record, config, |kmer| match kmer {
        KmerKey::Short(raw) => {
            if let Err(err) = writer.write_all(&raw.to_le_bytes()) {
                error = Some(err);
            } else {
                written += 1;
            }
        }
        KmerKey::LongHash(_) => {
            error = Some(std::io::Error::other("long k-mer emitted for k <= 31"));
        }
    });
    if let Some(err) = error {
        return Err(err.into());
    }
    Ok(written)
}

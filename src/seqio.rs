use anyhow::{Context, Result, bail};
use flate2::Compression;
use flate2::read::MultiGzDecoder;
use flate2::write::GzEncoder;
use gzp::deflate::Gzip;
use gzp::par::compress::ParCompressBuilder;
use std::env;
use std::ffi::OsString;
use std::fs::{File, OpenOptions};
use std::io::{self, BufRead, BufReader, BufWriter, Read, Write};
use std::path::{Path, PathBuf};
use std::process::{Child, ChildStdout, Command, Stdio};

const PHRED33_OFFSET: u8 = 33;
const HARD_MIN_CALLED_QUALITY: u8 = 2;
const QUALITY_CONVERT_BUFFER: usize = 8192;
const IO_BUFFER_CAPACITY: usize = 1024 * 1024;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum JunkMode {
    Crash,
    Fix,
    Flag,
    Ignore,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct BaseSettings {
    pub u_to_t: bool,
    pub to_upper_case: bool,
    pub lower_case_to_n: bool,
    pub dot_dash_x_to_n: bool,
    pub iupac_to_n: bool,
    pub fix_junk_and_iupac: bool,
    pub junk_mode: JunkMode,
}

impl Default for BaseSettings {
    fn default() -> Self {
        Self {
            u_to_t: false,
            to_upper_case: false,
            lower_case_to_n: false,
            dot_dash_x_to_n: false,
            iupac_to_n: false,
            fix_junk_and_iupac: false,
            junk_mode: JunkMode::Crash,
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct QualitySettings {
    pub input_offset: u8,
    pub min_called: u8,
    pub max_called: u8,
    pub change_quality: bool,
}

impl Default for QualitySettings {
    fn default() -> Self {
        Self {
            input_offset: PHRED33_OFFSET,
            min_called: 2,
            max_called: 50,
            change_quality: true,
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub struct SequenceSettings {
    pub bases: BaseSettings,
    pub qualities: QualitySettings,
}

impl From<QualitySettings> for SequenceSettings {
    fn from(qualities: QualitySettings) -> Self {
        Self {
            qualities,
            ..Self::default()
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SeqFormat {
    Fasta,
    Fastq,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct SequenceRecord {
    pub id: String,
    pub numeric_id: u64,
    pub bases: Vec<u8>,
    pub qualities: Option<Vec<u8>>,
}

impl SequenceRecord {
    pub fn len(&self) -> usize {
        self.bases.len()
    }

    pub fn is_empty(&self) -> bool {
        self.bases.is_empty()
    }

    pub fn renamed(&self, id: String) -> Self {
        let mut clone = self.clone();
        clone.id = id;
        clone
    }
}

pub struct SequenceReader {
    reader: Box<dyn BufRead>,
    format: SeqFormat,
    settings: SequenceSettings,
    pending_header: Option<String>,
    path: PathBuf,
    next_numeric_id: u64,
}

impl SequenceReader {
    pub fn from_path(path: &Path, settings: impl Into<SequenceSettings>) -> Result<Self> {
        Self::from_path_with_gzip_threads(path, settings, None)
    }

    pub fn from_path_with_gzip_threads(
        path: &Path,
        settings: impl Into<SequenceSettings>,
        gzip_threads: Option<usize>,
    ) -> Result<Self> {
        let mut reader = open_input_with_gzip_threads(path, gzip_threads)?;
        let format = detect_format(&mut reader, path)?;
        Ok(Self {
            reader,
            format,
            settings: settings.into(),
            pending_header: None,
            path: path.to_path_buf(),
            next_numeric_id: 0,
        })
    }

    pub fn format(&self) -> SeqFormat {
        self.format
    }

    pub fn next_record(&mut self) -> Result<Option<SequenceRecord>> {
        let mut record = match self.format {
            SeqFormat::Fastq => self.next_fastq(),
            SeqFormat::Fasta => self.next_fasta(),
        }?;
        if let Some(record) = &mut record {
            record.numeric_id = self.next_numeric_id;
            self.next_numeric_id += 1;
        }
        Ok(record)
    }

    fn next_fastq(&mut self) -> Result<Option<SequenceRecord>> {
        let mut header = String::new();
        loop {
            header.clear();
            let read = self.reader.read_line(&mut header)?;
            if read == 0 {
                return Ok(None);
            }
            trim_line(&mut header);
            if !header.is_empty() {
                break;
            }
        }
        if !header.starts_with('@') {
            bail!(
                "{}: expected FASTQ header starting with @, got {header:?}",
                self.path.display()
            );
        }

        let mut seq = String::new();
        let mut plus = String::new();
        let mut qual = String::new();
        if self.reader.read_line(&mut seq)? == 0
            || self.reader.read_line(&mut plus)? == 0
            || self.reader.read_line(&mut qual)? == 0
        {
            bail!("{}: truncated FASTQ record {header:?}", self.path.display());
        }
        trim_line(&mut seq);
        trim_line(&mut plus);
        trim_line(&mut qual);
        if !plus.starts_with('+') {
            bail!(
                "{}: expected FASTQ plus line for {header:?}",
                self.path.display()
            );
        }
        if seq.len() != qual.len() {
            bail!(
                "{}: sequence and quality lengths differ for {header:?} ({} vs {})",
                self.path.display(),
                seq.len(),
                qual.len()
            );
        }

        let mut bases = seq.into_bytes();
        let mut qualities = qual.into_bytes();
        normalize_fastq_qualities(&bases, &mut qualities, self.settings.qualities);
        normalize_bases(&mut bases, self.settings.bases)
            .with_context(|| format!("validating bases for FASTQ record {header:?}"))?;
        if self.settings.qualities.change_quality {
            cap_fastq_qualities(&bases, &mut qualities, self.settings.qualities);
        }

        Ok(Some(SequenceRecord {
            id: header[1..].to_string(),
            numeric_id: 0,
            bases,
            qualities: Some(qualities),
        }))
    }

    fn next_fasta(&mut self) -> Result<Option<SequenceRecord>> {
        let header = match self.pending_header.take() {
            Some(header) => header,
            None => loop {
                let mut line = String::new();
                let read = self.reader.read_line(&mut line)?;
                if read == 0 {
                    return Ok(None);
                }
                trim_line(&mut line);
                if line.is_empty() {
                    continue;
                }
                if !line.starts_with('>') {
                    bail!(
                        "{}: expected FASTA header starting with >, got {line:?}",
                        self.path.display()
                    );
                }
                break line;
            },
        };

        let mut bases = Vec::new();
        loop {
            let mut line = String::new();
            let read = self.reader.read_line(&mut line)?;
            if read == 0 {
                break;
            }
            trim_line(&mut line);
            if line.starts_with('>') {
                self.pending_header = Some(line);
                break;
            }
            if !line.is_empty() {
                bases.extend_from_slice(line.as_bytes());
            }
        }

        normalize_bases(&mut bases, self.settings.bases)
            .with_context(|| format!("validating bases for FASTA record {header:?}"))?;

        Ok(Some(SequenceRecord {
            id: header[1..].to_string(),
            numeric_id: 0,
            bases,
            qualities: None,
        }))
    }
}

pub fn detect_interleaved_input(
    path: &Path,
    settings: impl Into<SequenceSettings>,
) -> Result<bool> {
    detect_interleaved_input_with_gzip_threads(path, settings, None)
}

pub fn detect_interleaved_input_with_gzip_threads(
    path: &Path,
    settings: impl Into<SequenceSettings>,
    gzip_threads: Option<usize>,
) -> Result<bool> {
    if path.to_string_lossy().contains("_interleaved.") {
        return Ok(true);
    }

    let mut reader = SequenceReader::from_path_with_gzip_threads(path, settings, gzip_threads)?;
    let Some(first) = reader.next_record()? else {
        return Ok(false);
    };
    let Some(second) = reader.next_record()? else {
        return Ok(false);
    };
    Ok(record_ids_look_paired(&first.id, &second.id, false))
}

pub fn record_ids_look_paired(id1: &str, id2: &str, allow_identical: bool) -> bool {
    if id1.len() != id2.len() {
        return false;
    }

    let bytes1 = id1.as_bytes();
    let bytes2 = id2.as_bytes();
    let slash1 = id1.rfind('/');
    let slash2 = id2.rfind('/');
    let space1 = id1.find(' ');
    let space2 = id2.find(' ');

    if let (Some(idx1), Some(idx2)) = (space1, space2)
        && idx1 == idx2
        && idx1 > 0
        && id1.len() >= idx1 + 3
        && bytes1[idx1 + 1] == b'1'
        && bytes1[idx1 + 2] == b':'
        && bytes2[idx2 + 1] == b'2'
        && bytes2[idx2 + 2] == b':'
        && bytes1[..idx1] == bytes2[..idx2]
    {
        return true;
    }

    if let (Some(idx1), Some(idx2)) = (slash1, slash2)
        && idx1 == idx2
        && idx1 > 0
        && id1.len() >= idx1 + 2
        && bytes1[idx1 + 1] == b'1'
        && bytes2[idx2 + 1] == b'2'
        && bytes1[..idx1] == bytes2[..idx2]
        && bytes1[idx1 + 2..] == bytes2[idx2 + 2..]
    {
        return true;
    }

    allow_identical && id1 == id2
}

pub struct SequenceWriter {
    writer: Box<dyn Write>,
    format: SeqFormat,
    quality_out_offset: u8,
    fake_quality: u8,
    fasta_wrap: usize,
}

impl SequenceWriter {
    pub fn from_path(
        path: &Path,
        overwrite: bool,
        quality_out_offset: u8,
        fake_quality: u8,
        fasta_wrap: usize,
    ) -> Result<Self> {
        Self::from_path_with_append(
            path,
            overwrite,
            false,
            quality_out_offset,
            fake_quality,
            fasta_wrap,
        )
    }

    pub fn from_path_with_append(
        path: &Path,
        overwrite: bool,
        append: bool,
        quality_out_offset: u8,
        fake_quality: u8,
        fasta_wrap: usize,
    ) -> Result<Self> {
        Self::from_path_with_append_and_gzip_threads(
            path,
            overwrite,
            append,
            quality_out_offset,
            fake_quality,
            fasta_wrap,
            None,
        )
    }

    pub fn from_path_with_append_and_gzip_threads(
        path: &Path,
        overwrite: bool,
        append: bool,
        quality_out_offset: u8,
        fake_quality: u8,
        fasta_wrap: usize,
        gzip_threads: Option<usize>,
    ) -> Result<Self> {
        Ok(Self {
            writer: create_output_with_append_and_gzip_threads(
                path,
                overwrite,
                append,
                gzip_threads,
            )?,
            format: output_format(path),
            quality_out_offset,
            fake_quality,
            fasta_wrap,
        })
    }

    pub fn write_record(&mut self, record: &SequenceRecord) -> Result<()> {
        match self.format {
            SeqFormat::Fastq => {
                writeln!(self.writer, "@{}", record.id)?;
                self.writer.write_all(&record.bases)?;
                writeln!(self.writer)?;
                writeln!(self.writer, "+")?;
                if let Some(qualities) = record.qualities.as_ref() {
                    write_qualities(&mut self.writer, qualities, self.quality_out_offset)?;
                } else {
                    write_fake_qualities(
                        &mut self.writer,
                        &record.bases,
                        self.quality_out_offset,
                        self.fake_quality,
                    )?;
                }
                writeln!(self.writer)?;
            }
            SeqFormat::Fasta => {
                writeln!(self.writer, ">{}", record.id)?;
                write_wrapped(&mut self.writer, &record.bases, self.fasta_wrap)?;
            }
        }
        Ok(())
    }

    pub fn flush(&mut self) -> Result<()> {
        self.writer.flush()?;
        Ok(())
    }
}

fn normalize_fastq_qualities(bases: &[u8], qualities: &mut [u8], settings: QualitySettings) {
    for (&base, quality) in bases.iter().zip(qualities) {
        let normalized = if is_fully_defined_base(base) {
            quality
                .saturating_sub(settings.input_offset)
                .max(HARD_MIN_CALLED_QUALITY)
        } else {
            0
        };
        *quality = normalized.saturating_add(PHRED33_OFFSET);
    }
}

fn cap_fastq_qualities(bases: &[u8], qualities: &mut [u8], settings: QualitySettings) {
    for (&base, quality) in bases.iter().zip(qualities) {
        let normalized = quality.saturating_sub(PHRED33_OFFSET);
        *quality = if is_fully_defined_base(base) {
            normalized
                .clamp(settings.min_called, settings.max_called)
                .saturating_add(PHRED33_OFFSET)
        } else {
            PHRED33_OFFSET
        };
    }
}

fn normalize_bases(bases: &mut [u8], settings: BaseSettings) -> Result<()> {
    if settings.u_to_t {
        for base in bases.iter_mut() {
            *base = match *base {
                b'U' => b'T',
                b'u' => b't',
                other => other,
            };
        }
    }

    if settings.lower_case_to_n {
        for base in bases.iter_mut() {
            if base.is_ascii_lowercase() {
                *base = b'N';
            }
        }
    } else if settings.to_upper_case {
        for base in bases.iter_mut() {
            *base = base.to_ascii_uppercase();
        }
    }

    if settings.fix_junk_and_iupac || (settings.dot_dash_x_to_n && settings.iupac_to_n) {
        for base in bases.iter_mut() {
            *base = base_to_acgtn(*base);
        }
    } else {
        if settings.iupac_to_n {
            for base in bases.iter_mut() {
                if is_iupac_nocall(*base) {
                    *base = b'N';
                }
            }
        }
        if settings.dot_dash_x_to_n {
            for base in bases.iter_mut() {
                *base = match *base {
                    b'.' | b'-' | b'X' | b'x' | b'n' => b'N',
                    other => other,
                };
            }
        }
    }

    match settings.junk_mode {
        JunkMode::Flag | JunkMode::Ignore => Ok(()),
        JunkMode::Fix => {
            for base in bases.iter_mut() {
                if !is_extended_nucleotide(*base) {
                    *base = b'N';
                }
            }
            Ok(())
        }
        JunkMode::Crash => {
            if let Some(&base) = bases.iter().find(|&&base| !is_extended_nucleotide(base)) {
                bail!(
                    "invalid nucleotide byte {} ({:?}); use fixjunk=t or ignorejunk=t to match BBTools junk handling",
                    base,
                    base as char
                );
            }
            Ok(())
        }
    }
}

fn is_fully_defined_base(base: u8) -> bool {
    matches!(
        base,
        b'A' | b'C' | b'G' | b'T' | b'U' | b'a' | b'c' | b'g' | b't' | b'u'
    )
}

fn is_iupac_nocall(base: u8) -> bool {
    is_extended_nucleotide(base) && !is_fully_defined_base(base)
}

fn is_extended_nucleotide(base: u8) -> bool {
    matches!(
        base,
        b'A' | b'C'
            | b'G'
            | b'T'
            | b'U'
            | b'a'
            | b'c'
            | b'g'
            | b't'
            | b'u'
            | b'R'
            | b'Y'
            | b'S'
            | b'W'
            | b'K'
            | b'M'
            | b'B'
            | b'D'
            | b'H'
            | b'V'
            | b'N'
            | b'r'
            | b'y'
            | b's'
            | b'w'
            | b'k'
            | b'm'
            | b'b'
            | b'd'
            | b'h'
            | b'v'
            | b'n'
    )
}

fn base_to_acgtn(base: u8) -> u8 {
    match base {
        b'A' | b'a' => b'A',
        b'C' | b'c' => b'C',
        b'G' | b'g' => b'G',
        b'T' | b't' | b'U' | b'u' => b'T',
        _ => b'N',
    }
}

fn write_qualities(writer: &mut Box<dyn Write>, qualities: &[u8], output_offset: u8) -> Result<()> {
    if output_offset == PHRED33_OFFSET {
        writer.write_all(qualities)?;
        return Ok(());
    }

    let mut converted = [0u8; QUALITY_CONVERT_BUFFER];
    for chunk in qualities.chunks(QUALITY_CONVERT_BUFFER) {
        for (slot, &quality) in converted.iter_mut().zip(chunk) {
            *slot = quality
                .saturating_sub(PHRED33_OFFSET)
                .saturating_add(output_offset);
        }
        writer.write_all(&converted[..chunk.len()])?;
    }
    Ok(())
}

fn write_fake_qualities(
    writer: &mut Box<dyn Write>,
    bases: &[u8],
    output_offset: u8,
    fake_quality: u8,
) -> Result<()> {
    let mut converted = [0u8; QUALITY_CONVERT_BUFFER];
    for chunk in bases.chunks(QUALITY_CONVERT_BUFFER) {
        for (slot, &base) in converted.iter_mut().zip(chunk) {
            let quality = if is_fully_defined_base(base) {
                fake_quality
            } else {
                0
            };
            *slot = quality.saturating_add(output_offset);
        }
        writer.write_all(&converted[..chunk.len()])?;
    }
    Ok(())
}

fn output_format(path: &Path) -> SeqFormat {
    let stem = if has_gz_extension(path) {
        path.file_stem()
            .and_then(|stem| Path::new(stem).extension())
            .and_then(|ext| ext.to_str())
    } else {
        path.extension().and_then(|ext| ext.to_str())
    };

    match stem.map(str::to_ascii_lowercase).as_deref() {
        Some("fa" | "fasta" | "fas" | "fna" | "ffn" | "frn" | "fsa" | "seq") => SeqFormat::Fasta,
        _ => SeqFormat::Fastq,
    }
}

pub fn create_output(path: &Path, overwrite: bool) -> Result<Box<dyn Write>> {
    create_output_with_append(path, overwrite, false)
}

pub fn create_output_with_append(
    path: &Path,
    overwrite: bool,
    append: bool,
) -> Result<Box<dyn Write>> {
    create_output_with_append_and_gzip_threads(path, overwrite, append, None)
}

pub fn create_output_with_append_and_gzip_threads(
    path: &Path,
    overwrite: bool,
    append: bool,
    gzip_threads: Option<usize>,
) -> Result<Box<dyn Write>> {
    if is_null_output(path) {
        return Ok(Box::new(io::sink()));
    }
    if is_stdout(path) {
        return Ok(Box::new(io::stdout()));
    }

    let mut options = OpenOptions::new();
    options.write(true).create(true);
    if append {
        options.append(true);
    } else {
        options.truncate(true);
    }
    if !overwrite && !append {
        options.create_new(true);
    }
    let file = options
        .open(path)
        .with_context(|| format!("opening output {}", path.display()))?;
    let buffered = BufWriter::with_capacity(IO_BUFFER_CAPACITY, file);
    if has_gz_extension(path) {
        if let Some(threads) = gzip_threads.filter(|threads| *threads > 1) {
            let writer = ParCompressBuilder::<Gzip>::new()
                .num_threads(threads)
                .context("configuring parallel gzip output")?
                .from_writer(buffered);
            return Ok(Box::new(writer));
        }
        Ok(Box::new(GzEncoder::new(buffered, Compression::default())))
    } else {
        Ok(Box::new(buffered))
    }
}

fn open_input_with_gzip_threads(
    path: &Path,
    gzip_threads: Option<usize>,
) -> Result<Box<dyn BufRead>> {
    let file = File::open(path).with_context(|| format!("opening input {}", path.display()))?;
    if has_gz_extension(path) {
        if let Some(threads) = gzip_threads.filter(|threads| *threads > 1)
            && let Some(reader) = spawn_external_gzip_decoder(path, threads)
        {
            return Ok(Box::new(reader));
        }
        let decoder = MultiGzDecoder::new(file);
        Ok(Box::new(BufReader::with_capacity(
            IO_BUFFER_CAPACITY,
            decoder,
        )))
    } else {
        Ok(Box::new(BufReader::with_capacity(IO_BUFFER_CAPACITY, file)))
    }
}

fn spawn_external_gzip_decoder(path: &Path, threads: usize) -> Option<ChildGzipReader> {
    for candidate in external_gzip_decoder_candidates() {
        let mut command = Command::new(&candidate.executable);
        if candidate.program == "pigz" {
            command.arg("-dc");
        } else {
            command.arg("-c");
        }
        command
            .arg("-p")
            .arg(threads.to_string())
            .arg(path)
            .stdin(Stdio::null())
            .stdout(Stdio::piped())
            .stderr(Stdio::null());

        if let Ok(mut child) = command.spawn()
            && let Some(stdout) = child.stdout.take()
        {
            return Some(ChildGzipReader::new(
                candidate.executable.display().to_string(),
                child,
                stdout,
            ));
        }
    }
    None
}

#[derive(Debug, Clone, PartialEq, Eq)]
struct ExternalGzipDecoderCandidate {
    program: &'static str,
    executable: PathBuf,
}

fn external_gzip_decoder_candidates() -> Vec<ExternalGzipDecoderCandidate> {
    external_gzip_decoder_candidates_from_env(env::var_os("PATH"), env::var_os("HOME"))
}

fn external_gzip_decoder_candidates_from_env(
    path: Option<OsString>,
    home: Option<OsString>,
) -> Vec<ExternalGzipDecoderCandidate> {
    let mut candidates = Vec::new();
    for program in ["pigz", "unpigz"] {
        add_path_decoder_candidates(&mut candidates, program, path.as_ref());
    }
    if let Some(home) = home {
        let local_bin = PathBuf::from(home).join(".local").join("bin");
        for program in ["pigz", "unpigz"] {
            add_decoder_candidate(&mut candidates, program, local_bin.join(program));
        }
    }
    candidates
}

fn add_path_decoder_candidates(
    candidates: &mut Vec<ExternalGzipDecoderCandidate>,
    program: &'static str,
    path: Option<&OsString>,
) {
    let Some(path) = path else {
        return;
    };
    for dir in env::split_paths(path) {
        add_decoder_candidate(candidates, program, dir.join(program));
    }
}

fn add_decoder_candidate(
    candidates: &mut Vec<ExternalGzipDecoderCandidate>,
    program: &'static str,
    executable: PathBuf,
) {
    if !executable.is_file() {
        return;
    }
    if candidates
        .iter()
        .any(|candidate| candidate.executable == executable)
    {
        return;
    }
    candidates.push(ExternalGzipDecoderCandidate {
        program,
        executable,
    });
}

struct ChildGzipReader {
    program: String,
    child: Child,
    stdout: Option<BufReader<ChildStdout>>,
    status_checked: bool,
}

impl ChildGzipReader {
    fn new(program: String, child: Child, stdout: ChildStdout) -> Self {
        Self {
            program,
            child,
            stdout: Some(BufReader::with_capacity(IO_BUFFER_CAPACITY, stdout)),
            status_checked: false,
        }
    }

    fn stdout_mut(&mut self) -> io::Result<&mut BufReader<ChildStdout>> {
        self.stdout.as_mut().ok_or_else(|| {
            io::Error::new(
                io::ErrorKind::BrokenPipe,
                "parallel gzip decoder stdout was already closed",
            )
        })
    }

    fn check_child_status(&mut self) -> io::Result<()> {
        if self.status_checked {
            return Ok(());
        }
        self.status_checked = true;
        let status = self.child.wait()?;
        if status.success() {
            Ok(())
        } else {
            Err(io::Error::new(
                io::ErrorKind::InvalidData,
                format!("{} exited with status {status}", self.program),
            ))
        }
    }
}

impl Read for ChildGzipReader {
    fn read(&mut self, buf: &mut [u8]) -> io::Result<usize> {
        let read = self.stdout_mut()?.read(buf)?;
        if read == 0 {
            self.check_child_status()?;
        }
        Ok(read)
    }
}

impl BufRead for ChildGzipReader {
    fn fill_buf(&mut self) -> io::Result<&[u8]> {
        let empty = self.stdout_mut()?.fill_buf()?.is_empty();
        if empty {
            self.check_child_status()?;
        }
        self.stdout_mut()?.fill_buf()
    }

    fn consume(&mut self, amt: usize) {
        if let Some(stdout) = self.stdout.as_mut() {
            stdout.consume(amt);
        }
    }
}

impl Drop for ChildGzipReader {
    fn drop(&mut self) {
        self.stdout.take();
        if !self.status_checked {
            let _ = self.child.kill();
            let _ = self.child.wait();
        }
    }
}

fn detect_format(reader: &mut Box<dyn BufRead>, path: &Path) -> Result<SeqFormat> {
    loop {
        let buf = reader.fill_buf()?;
        if buf.is_empty() {
            bail!("{}: empty sequence file", path.display());
        }
        match buf[0] {
            b'>' => return Ok(SeqFormat::Fasta),
            b'@' => return Ok(SeqFormat::Fastq),
            b'\n' | b'\r' | b'\t' | b' ' => {
                reader.consume(1);
            }
            other => bail!(
                "{}: could not detect FASTA/FASTQ format from leading byte {:?}",
                path.display(),
                other as char
            ),
        }
    }
}

fn write_wrapped(writer: &mut Box<dyn Write>, bytes: &[u8], width: usize) -> Result<()> {
    if bytes.is_empty() {
        writeln!(writer)?;
        return Ok(());
    }
    if width == 0 {
        writer.write_all(bytes)?;
        writeln!(writer)?;
        return Ok(());
    }
    for chunk in bytes.chunks(width) {
        writer.write_all(chunk)?;
        writeln!(writer)?;
    }
    Ok(())
}

fn trim_line(line: &mut String) {
    while matches!(line.as_bytes().last(), Some(b'\n' | b'\r')) {
        line.pop();
    }
}

fn has_gz_extension(path: &Path) -> bool {
    path.extension()
        .and_then(|ext| ext.to_str())
        .is_some_and(|ext| ext.eq_ignore_ascii_case("gz"))
}

fn is_stdout(path: &Path) -> bool {
    path == Path::new("-") || path == Path::new("stdout") || path == Path::new("stdout.fq")
}

fn is_null_output(path: &Path) -> bool {
    path.to_str()
        .is_some_and(|text| text.eq_ignore_ascii_case("null"))
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Write;
    use tempfile::{Builder, NamedTempFile};

    #[test]
    fn reads_fastq_records() {
        let mut file = NamedTempFile::new().unwrap();
        writeln!(file, "@r1\nACGT\n+\nIIII").unwrap();
        let mut reader =
            SequenceReader::from_path(file.path(), QualitySettings::default()).unwrap();
        let record = reader.next_record().unwrap().unwrap();
        assert_eq!(record.id, "r1");
        assert_eq!(record.bases, b"ACGT");
        assert_eq!(record.qualities.as_deref(), Some(&b"IIII"[..]));
        assert!(reader.next_record().unwrap().is_none());
    }

    #[test]
    fn normalizes_fastq_qualities_like_bbnorm() {
        let mut file = NamedTempFile::new().unwrap();
        writeln!(file, "@r1\nACGTN\n+\n!#I~I").unwrap();
        let mut reader =
            SequenceReader::from_path(file.path(), QualitySettings::default()).unwrap();
        let record = reader.next_record().unwrap().unwrap();
        assert_eq!(record.qualities.as_deref(), Some(&b"##IS!"[..]));

        let mut file = NamedTempFile::new().unwrap();
        writeln!(file, "@r1\nACGTN\n+\n@Bh|h").unwrap();
        let settings = QualitySettings {
            input_offset: 64,
            ..QualitySettings::default()
        };
        let mut reader = SequenceReader::from_path(file.path(), settings).unwrap();
        let record = reader.next_record().unwrap().unwrap();
        assert_eq!(record.qualities.as_deref(), Some(&b"##IS!"[..]));
    }

    #[test]
    fn applies_quality_change_controls_like_bbnorm() {
        let mut file = NamedTempFile::new().unwrap();
        writeln!(file, "@r1\nACGTN\n+\n!#I~I").unwrap();
        let settings = QualitySettings {
            min_called: 5,
            max_called: 30,
            ..QualitySettings::default()
        };
        let mut reader = SequenceReader::from_path(file.path(), settings).unwrap();
        let record = reader.next_record().unwrap().unwrap();
        assert_eq!(record.qualities.as_deref(), Some(&b"&&??!"[..]));

        let mut file = NamedTempFile::new().unwrap();
        writeln!(file, "@r1\nACGTN\n+\n!#I~I").unwrap();
        let settings = QualitySettings {
            change_quality: false,
            ..QualitySettings::default()
        };
        let mut reader = SequenceReader::from_path(file.path(), settings).unwrap();
        let record = reader.next_record().unwrap().unwrap();
        assert_eq!(record.qualities.as_deref(), Some(&b"##I~!"[..]));
    }

    #[test]
    fn applies_base_cleanup_like_bbnorm() {
        let cases: [(&str, BaseSettings, &[u8], &[u8]); 7] = [
            (
                "default",
                BaseSettings::default(),
                b"acgtuURYSWKMBDHVNn",
                b"IIIIII!!!!!!!!!!!!",
            ),
            (
                "utot",
                BaseSettings {
                    u_to_t: true,
                    ..BaseSettings::default()
                },
                b"acgttTRYSWKMBDHVNn",
                b"IIIIII!!!!!!!!!!!!",
            ),
            (
                "touppercase",
                BaseSettings {
                    to_upper_case: true,
                    ..BaseSettings::default()
                },
                b"ACGTUURYSWKMBDHVNN",
                b"IIIIII!!!!!!!!!!!!",
            ),
            (
                "lowercaseton",
                BaseSettings {
                    lower_case_to_n: true,
                    ..BaseSettings::default()
                },
                b"NNNNNURYSWKMBDHVNN",
                b"!!!!!I!!!!!!!!!!!!",
            ),
            (
                "iupacton",
                BaseSettings {
                    iupac_to_n: true,
                    ..BaseSettings::default()
                },
                b"acgtuUNNNNNNNNNNNN",
                b"IIIIII!!!!!!!!!!!!",
            ),
            (
                "iupacton_utot_touppercase",
                BaseSettings {
                    u_to_t: true,
                    to_upper_case: true,
                    iupac_to_n: true,
                    ..BaseSettings::default()
                },
                b"ACGTTTNNNNNNNNNNNN",
                b"IIIIII!!!!!!!!!!!!",
            ),
            (
                "dotdashxton",
                BaseSettings {
                    dot_dash_x_to_n: true,
                    ..BaseSettings::default()
                },
                b"ACGTNNNNNN",
                b"IIII!!!!!!",
            ),
        ];

        for (label, base_settings, expected_bases, expected_qualities) in cases {
            let input_bases = if label == "dotdashxton" {
                "ACGT.-XxNn"
            } else {
                "acgtuURYSWKMBDHVNn"
            };
            let mut file = NamedTempFile::new().unwrap();
            writeln!(
                file,
                "@r1\n{input_bases}\n+\n{}",
                "I".repeat(input_bases.len())
            )
            .unwrap();
            let settings = SequenceSettings {
                bases: base_settings,
                ..SequenceSettings::default()
            };
            let mut reader = SequenceReader::from_path(file.path(), settings).unwrap();
            let record = reader.next_record().unwrap().unwrap();
            assert_eq!(record.bases, expected_bases, "{label} bases");
            assert_eq!(
                record.qualities.as_deref(),
                Some(expected_qualities),
                "{label} qualities"
            );
        }
    }

    #[test]
    fn lowercaseton_respects_changequality_false_like_bbnorm() {
        let mut file = NamedTempFile::new().unwrap();
        writeln!(file, "@r1\nacgtnu\n+\nIIIIII").unwrap();
        let settings = SequenceSettings {
            bases: BaseSettings {
                lower_case_to_n: true,
                ..BaseSettings::default()
            },
            qualities: QualitySettings {
                change_quality: false,
                ..QualitySettings::default()
            },
        };
        let mut reader = SequenceReader::from_path(file.path(), settings).unwrap();
        let record = reader.next_record().unwrap().unwrap();
        assert_eq!(record.bases, b"NNNNNN");
        assert_eq!(record.qualities.as_deref(), Some(&b"IIII!I"[..]));
    }

    #[test]
    fn handles_junk_modes_like_bbnorm() {
        let mut file = NamedTempFile::new().unwrap();
        writeln!(file, "@r1\nACGT?ZN\n+\nIIIIIII").unwrap();
        let err = SequenceReader::from_path(file.path(), SequenceSettings::default())
            .and_then(|mut reader| reader.next_record().map(|_| ()));
        assert!(err.is_err());

        let mut file = NamedTempFile::new().unwrap();
        writeln!(file, "@r1\nACGT?ZN\n+\nIIIIIII").unwrap();
        let settings = SequenceSettings {
            bases: BaseSettings {
                junk_mode: JunkMode::Fix,
                ..BaseSettings::default()
            },
            ..SequenceSettings::default()
        };
        let mut reader = SequenceReader::from_path(file.path(), settings).unwrap();
        let record = reader.next_record().unwrap().unwrap();
        assert_eq!(record.bases, b"ACGTNNN");
        assert_eq!(record.qualities.as_deref(), Some(&b"IIII!!!"[..]));

        let mut file = NamedTempFile::new().unwrap();
        writeln!(file, "@r1\nACGT?ZN\n+\nIIIIIII").unwrap();
        let settings = SequenceSettings {
            bases: BaseSettings {
                junk_mode: JunkMode::Flag,
                ..BaseSettings::default()
            },
            ..SequenceSettings::default()
        };
        let mut reader = SequenceReader::from_path(file.path(), settings).unwrap();
        let record = reader.next_record().unwrap().unwrap();
        assert_eq!(record.bases, b"ACGT?ZN");
        assert_eq!(record.qualities.as_deref(), Some(&b"IIII!!!"[..]));

        let mut file = NamedTempFile::new().unwrap();
        writeln!(file, "@r1\nACGT?ZN\n+\nIIIIIII").unwrap();
        let settings = SequenceSettings {
            bases: BaseSettings {
                junk_mode: JunkMode::Ignore,
                ..BaseSettings::default()
            },
            ..SequenceSettings::default()
        };
        let mut reader = SequenceReader::from_path(file.path(), settings).unwrap();
        let record = reader.next_record().unwrap().unwrap();
        assert_eq!(record.bases, b"ACGT?ZN");
        assert_eq!(record.qualities.as_deref(), Some(&b"IIII!!!"[..]));
    }

    #[test]
    fn writer_infers_output_format_and_fakes_fasta_qualities_like_bbnorm() {
        let record = SequenceRecord {
            id: "r1".to_string(),
            numeric_id: 0,
            bases: b"ACGTNN".to_vec(),
            qualities: None,
        };

        let file = Builder::new().suffix(".fq").tempfile().unwrap();
        let mut writer = SequenceWriter::from_path(file.path(), true, 33, 30, 70).unwrap();
        writer.write_record(&record).unwrap();
        writer.flush().unwrap();
        assert_eq!(
            std::fs::read_to_string(file.path()).unwrap(),
            "@r1\nACGTNN\n+\n????!!\n"
        );

        let file = Builder::new().suffix(".fq").tempfile().unwrap();
        let mut writer = SequenceWriter::from_path(file.path(), true, 64, 30, 70).unwrap();
        writer.write_record(&record).unwrap();
        writer.flush().unwrap();
        assert_eq!(
            std::fs::read_to_string(file.path()).unwrap(),
            "@r1\nACGTNN\n+\n^^^^@@\n"
        );

        let record = SequenceRecord {
            id: "r2".to_string(),
            numeric_id: 0,
            bases: b"ACGT".to_vec(),
            qualities: Some(b"IIII".to_vec()),
        };
        let file = Builder::new().suffix(".fa").tempfile().unwrap();
        let mut writer = SequenceWriter::from_path(file.path(), true, 33, 30, 70).unwrap();
        writer.write_record(&record).unwrap();
        writer.flush().unwrap();
        assert_eq!(std::fs::read_to_string(file.path()).unwrap(), ">r2\nACGT\n");
    }

    #[test]
    fn writer_appends_sequence_outputs_like_bbnorm() {
        let record = SequenceRecord {
            id: "r1".to_string(),
            numeric_id: 0,
            bases: b"ACGT".to_vec(),
            qualities: Some(b"IIII".to_vec()),
        };
        let file = Builder::new().suffix(".fq").tempfile().unwrap();
        std::fs::write(file.path(), b"PREEXISTING\n").unwrap();

        let mut writer =
            SequenceWriter::from_path_with_append(file.path(), false, true, 33, 30, 70).unwrap();
        writer.write_record(&record).unwrap();
        writer.flush().unwrap();

        assert_eq!(
            std::fs::read_to_string(file.path()).unwrap(),
            "PREEXISTING\n@r1\nACGT\n+\nIIII\n"
        );
    }

    #[test]
    fn writer_parallel_gzip_output_round_trips_fastq() {
        let record = SequenceRecord {
            id: "r1".to_string(),
            numeric_id: 0,
            bases: b"ACGTNNACGT".to_vec(),
            qualities: Some(b"IIII!!IIII".to_vec()),
        };
        let file = Builder::new().suffix(".fq.gz").tempfile().unwrap();

        let mut writer = SequenceWriter::from_path_with_append_and_gzip_threads(
            file.path(),
            true,
            false,
            33,
            30,
            70,
            Some(2),
        )
        .unwrap();
        writer.write_record(&record).unwrap();
        writer.flush().unwrap();
        drop(writer);

        let mut reader =
            SequenceReader::from_path(file.path(), SequenceSettings::default()).unwrap();
        let round_tripped = reader.next_record().unwrap().unwrap();
        assert_eq!(round_tripped.id, record.id);
        assert_eq!(round_tripped.bases, record.bases);
        assert_eq!(round_tripped.qualities, record.qualities);
        assert!(reader.next_record().unwrap().is_none());
    }

    #[test]
    fn reader_threaded_gzip_input_round_trips_fastq() {
        let record = SequenceRecord {
            id: "r1".to_string(),
            numeric_id: 0,
            bases: b"ACGTACGT".to_vec(),
            qualities: Some(b"IIIIIIII".to_vec()),
        };
        let file = Builder::new().suffix(".fq.gz").tempfile().unwrap();
        {
            let mut writer = SequenceWriter::from_path(file.path(), true, 33, 30, 70).unwrap();
            writer.write_record(&record).unwrap();
            writer.flush().unwrap();
        }

        let mut reader = SequenceReader::from_path_with_gzip_threads(
            file.path(),
            SequenceSettings::default(),
            Some(2),
        )
        .unwrap();
        let round_tripped = reader.next_record().unwrap().unwrap();
        assert_eq!(round_tripped.id, record.id);
        assert_eq!(round_tripped.bases, record.bases);
        assert_eq!(round_tripped.qualities, record.qualities);
        assert!(reader.next_record().unwrap().is_none());
    }

    #[test]
    fn external_gzip_decoder_candidates_include_home_local_bin() {
        let path_dir = Builder::new().prefix("bbnorm-rs-path-").tempdir().unwrap();
        let home_dir = Builder::new().prefix("bbnorm-rs-home-").tempdir().unwrap();
        let path_pigz = path_dir.path().join("pigz");
        std::fs::write(&path_pigz, b"fake pigz").unwrap();
        let home_bin = home_dir.path().join(".local").join("bin");
        std::fs::create_dir_all(&home_bin).unwrap();
        let home_unpigz = home_bin.join("unpigz");
        std::fs::write(&home_unpigz, b"fake unpigz").unwrap();

        let path = std::env::join_paths([path_dir.path()]).unwrap();
        let candidates = external_gzip_decoder_candidates_from_env(
            Some(path),
            Some(home_dir.path().as_os_str().to_os_string()),
        );

        assert_eq!(
            candidates.first(),
            Some(&ExternalGzipDecoderCandidate {
                program: "pigz",
                executable: path_pigz,
            })
        );
        assert!(
            candidates
                .iter()
                .any(|candidate| candidate.program == "unpigz"
                    && candidate.executable == home_unpigz),
            "expected ~/.local/bin/unpigz fallback in {candidates:?}"
        );
    }

    #[test]
    fn writer_treats_null_outputs_as_sinks_like_bbnorm() {
        let record = SequenceRecord {
            id: "r1".to_string(),
            numeric_id: 0,
            bases: b"ACGT".to_vec(),
            qualities: Some(b"IIII".to_vec()),
        };
        for path in [Path::new("null"), Path::new("NULL"), Path::new("Null")] {
            let _ = std::fs::remove_file(path);
            let mut writer = SequenceWriter::from_path(path, false, 33, 30, 70).unwrap();
            writer.write_record(&record).unwrap();
            writer.flush().unwrap();
            assert!(
                !path.exists(),
                "{} should be treated as a null output sink, not a real file",
                path.display()
            );
        }

        let none_path = Path::new("none");
        let _ = std::fs::remove_file(none_path);
        let mut writer = SequenceWriter::from_path(none_path, false, 33, 30, 70).unwrap();
        writer.write_record(&record).unwrap();
        writer.flush().unwrap();
        assert!(
            none_path.exists(),
            "BBTools treats none as a literal output path, not a null sink"
        );
        std::fs::remove_file(none_path).unwrap();
    }

    #[test]
    fn writer_wraps_fasta_output_like_bbnorm() {
        let record = SequenceRecord {
            id: "long".to_string(),
            numeric_id: 0,
            bases: b"A".repeat(100),
            qualities: Some(b"I".repeat(100)),
        };

        let seventy = "A".repeat(70);
        let thirty = "A".repeat(30);
        let file = Builder::new().suffix(".fa").tempfile().unwrap();
        let mut writer = SequenceWriter::from_path(file.path(), true, 33, 30, 70).unwrap();
        writer.write_record(&record).unwrap();
        writer.flush().unwrap();
        assert_eq!(
            std::fs::read_to_string(file.path()).unwrap(),
            format!(">long\n{seventy}\n{thirty}\n")
        );

        let twenty = "A".repeat(20);
        let file = Builder::new().suffix(".fa").tempfile().unwrap();
        let mut writer = SequenceWriter::from_path(file.path(), true, 33, 30, 20).unwrap();
        writer.write_record(&record).unwrap();
        writer.flush().unwrap();
        assert_eq!(
            std::fs::read_to_string(file.path()).unwrap(),
            format!(">long\n{twenty}\n{twenty}\n{twenty}\n{twenty}\n{twenty}\n")
        );

        let hundred = "A".repeat(100);
        let file = Builder::new().suffix(".fa").tempfile().unwrap();
        let mut writer = SequenceWriter::from_path(file.path(), true, 33, 30, 0).unwrap();
        writer.write_record(&record).unwrap();
        writer.flush().unwrap();
        assert_eq!(
            std::fs::read_to_string(file.path()).unwrap(),
            format!(">long\n{hundred}\n")
        );
    }

    #[test]
    fn reads_multiline_fasta() {
        let mut file = NamedTempFile::new().unwrap();
        writeln!(file, ">r1\nAC\nGT\n>r2\nTT").unwrap();
        let mut reader =
            SequenceReader::from_path(file.path(), QualitySettings::default()).unwrap();
        assert_eq!(reader.next_record().unwrap().unwrap().bases, b"ACGT");
        assert_eq!(reader.next_record().unwrap().unwrap().bases, b"TT");
        assert!(reader.next_record().unwrap().is_none());
    }

    #[test]
    fn detects_common_pair_name_shapes() {
        assert!(record_ids_look_paired("read/1", "read/2", false));
        assert!(record_ids_look_paired(
            "read 1:N:0:ATCACG",
            "read 2:N:0:ATCACG",
            false
        ));
        assert!(!record_ids_look_paired("read/1", "other/2", false));
        assert!(!record_ids_look_paired("same", "same", false));
        assert!(record_ids_look_paired("same", "same", true));
    }
}

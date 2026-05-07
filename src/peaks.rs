use crate::cli::Config;
use crate::seqio::create_output;
use anyhow::{Context, Result};
use std::io::Write;
use std::path::Path;

const MAX_WIDTH_MULT: f64 = 2.5;
const PROGRESSIVE_MULT: f64 = 2.0;
const MAX_RADIUS: i32 = 10;

#[derive(Clone, Debug)]
struct Peak {
    start: usize,
    center: usize,
    stop: usize,
    max_pos: usize,
    max_height: u64,
    start_height: u64,
    stop_height: u64,
    left_min: u64,
    right_min: u64,
    volume: u64,
    volume2: u64,
}

impl Peak {
    #[allow(clippy::too_many_arguments)]
    fn new(
        center: usize,
        start: usize,
        stop: usize,
        max_pos: usize,
        max_height: u64,
        start_height: u64,
        stop_height: u64,
        left_min: u64,
        right_min: u64,
        volume: u64,
        volume2: u64,
    ) -> Self {
        Self {
            start,
            center,
            stop,
            max_pos,
            max_height,
            start_height,
            stop_height,
            left_min,
            right_min,
            volume,
            volume2,
        }
    }

    fn recalculate(&mut self, array: &[u64]) {
        self.max_height = array[self.center];
        self.start_height = array[self.start];
        self.stop_height = array[self.stop];
        self.left_min = self.start_height;
        self.right_min = self.stop_height;
        self.max_pos = self.center;
        self.volume = 0;
        self.volume2 = 0;

        for (i, &x) in array.iter().enumerate().take(self.stop).skip(self.start) {
            if x > self.max_height {
                self.max_pos = i;
                self.max_height = x;
            }
            if i < self.center {
                self.left_min = self.left_min.min(x);
            } else if i > self.center {
                self.right_min = self.right_min.min(x);
            }
            self.volume = self.volume.saturating_add(x);
            self.volume2 = self.volume2.saturating_add(x.saturating_mul(i as u64));
        }
    }

    fn compatible_with(&self, other: &Self) -> bool {
        let min = self.center.min(other.stop) as f64;
        let max = self.stop.max(other.center) as f64;
        min * MAX_WIDTH_MULT >= max
    }

    fn absorb(&mut self, other: &Self) {
        if self.center > other.center {
            if self.start > other.start {
                self.start = other.start;
                self.start_height = other.start_height;
            }
            self.left_min = self.left_min.min(other.left_min);
        } else {
            if self.stop < other.stop {
                self.stop = other.stop;
                self.stop_height = other.stop_height;
            }
            self.right_min = self.right_min.min(other.right_min);
        }

        if self.max_height < other.max_height {
            self.max_height = other.max_height;
            self.max_pos = other.max_pos;
        }
        self.volume = self.volume.saturating_add(other.volume);
        self.volume2 = self.volume2.saturating_add(other.volume2);
    }
}

pub(crate) fn write_peaks(path: &Path, raw_hist: &[u64], config: &Config) -> Result<()> {
    let unique_hist: Vec<u64> = raw_hist
        .iter()
        .copied()
        .enumerate()
        .map(|(depth, raw)| unique_from_raw(depth, raw))
        .collect();
    let peaks = call_peaks(&unique_hist, config);
    let unique_kmers: u64 = unique_hist.iter().sum();

    let mut writer = create_output(path, config.overwrite || config.append)
        .with_context(|| format!("creating peaks {}", path.display()))?;
    print_peaks(&mut writer, &peaks, config, unique_kmers, &unique_hist)?;
    writer.flush()?;
    Ok(())
}

fn call_peaks(original: &[u64], config: &Config) -> Vec<Peak> {
    let array = smooth_progressive(original, 1);
    let length = original.len();
    let Some(dip0) = (1..length).find_map(|i| (array[i - 1] < array[i]).then_some(i - 1)) else {
        return Vec::new();
    };

    let mut peaks = Vec::new();
    let mut mode_up = true;
    let mut start = dip0;
    let mut center = 0usize;
    let mut prev = array[dip0];
    let mut sum = prev;
    let mut sum2 = prev.saturating_mul(dip0 as u64);
    let mut i = dip0 + 1;

    while i < length {
        let x = array[i];
        if mode_up {
            if x < prev {
                mode_up = false;
                center = i - 1;
            }
        } else if x > prev {
            let mut stop = i - 1;
            let max = array[center];
            if peak_passes_filters(center, start, stop, max, sum, config) {
                center = middle_of_mesa(&array, center, max);
                stop = middle_of_valley(&array, stop);
                let height1 = array[start];
                let height2 = array[stop];
                peaks.push(Peak::new(
                    center, start, stop, center, max, height1, height2, height1, height2, sum, sum2,
                ));
            }
            start = stop;
            sum = 0;
            sum2 = 0;
            if i > config.peak_max_peak {
                break;
            }
            while i < array.len() && array[i] == 0 {
                i += 1;
            }
            mode_up = true;
        }

        sum = sum.saturating_add(x);
        sum2 = sum2.saturating_add(x.saturating_mul(i as u64));
        prev = x;
        i += 1;
    }

    if !mode_up {
        let mut stop = length;
        let max = array[center];
        center = middle_of_mesa(&array, center, max);
        let valley = array[stop - 1];
        for j in (0..stop).rev() {
            if array[j] != valley {
                stop = if valley == 0 {
                    j + 1
                } else {
                    (stop - 1 + j + 2) / 2
                };
                break;
            }
        }
        let clamped_stop = stop.min(length - 1);
        if peak_passes_filters(center, start, stop, max, sum, config) {
            let height1 = array[start];
            let height2 = array[clamped_stop];
            peaks.push(Peak::new(
                center,
                start,
                clamped_stop,
                center,
                max,
                height1,
                height2,
                height1,
                height2,
                sum,
                sum2,
            ));
        }
    }

    cap_width(&mut peaks, &array);
    if config.peak_max_count < peaks.len() {
        peaks = condense(&peaks, config.peak_max_count);
    }
    cap_width(&mut peaks, &array);

    if peaks.len() > 1 {
        let biggest_volume = peaks[biggest_peak(&peaks)].volume as f64;
        while peaks.len() > 1 && (peaks[0].volume as f64) < 0.0001 * biggest_volume {
            peaks.remove(0);
        }
    }

    for peak in &mut peaks {
        peak.recalculate(original);
    }
    peaks
        .into_iter()
        .filter(|peak| peak.volume >= config.peak_min_volume)
        .collect()
}

fn peak_passes_filters(
    center: usize,
    start: usize,
    stop: usize,
    max: u64,
    sum: u64,
    config: &Config,
) -> bool {
    center >= config.peak_min_peak
        && center <= config.peak_max_peak
        && max >= config.peak_min_height
        && stop.saturating_sub(start) >= config.peak_min_width
        && sum >= config.peak_min_volume
}

fn middle_of_mesa(array: &[u64], center: usize, max: u64) -> usize {
    for j in (0..center).rev() {
        if array[j] != max {
            return (center + j + 2) / 2;
        }
    }
    center
}

fn middle_of_valley(array: &[u64], stop: usize) -> usize {
    let valley = array[stop];
    for j in (0..=stop).rev() {
        if array[j] != valley {
            return if valley == 0 {
                j + 1
            } else {
                (stop + j + 2) / 2
            };
        }
    }
    stop
}

fn cap_width(peaks: &mut [Peak], counts: &[u64]) {
    let mult = 1.0 / MAX_WIDTH_MULT;
    for peak in peaks {
        peak.start = java_round_f64((peak.start as f64).max(peak.center as f64 * mult)) as usize;
        peak.stop =
            java_round_f64((peak.stop as f64).min(peak.center as f64 * MAX_WIDTH_MULT)) as usize;
        peak.stop = peak.stop.min(counts.len().saturating_sub(1));
        peak.start = peak.start.min(peak.center);
        peak.recalculate(counts);
    }
}

fn condense(input: &[Peak], max_count: usize) -> Vec<Peak> {
    if input.is_empty() {
        return Vec::new();
    }
    let max_count = max_count.min(input.len()).max(1);
    let mut heights: Vec<u64> = input.iter().map(|peak| peak.max_height).collect();
    heights.sort_unstable();
    let hlimit = heights[heights.len() - max_count];

    let mc2 = max_count.div_ceil(2);
    let mut volumes: Vec<u64> = input.iter().map(|peak| peak.volume).collect();
    volumes.sort_unstable();
    let vlimit = volumes[volumes.len() - mc2];

    let mut output = Vec::with_capacity(max_count.min(input.len()));
    for peak in input {
        if peak.volume >= vlimit || peak.max_height >= hlimit {
            output.push(peak.clone());
        }
    }

    for peak in input {
        if peak.volume < vlimit && peak.max_height < hlimit {
            let Some((closest_idx, _)) = output
                .iter()
                .enumerate()
                .min_by_key(|(_, kept)| peak.center.abs_diff(kept.center))
            else {
                continue;
            };
            if output[closest_idx].compatible_with(peak) {
                output[closest_idx].absorb(peak);
            }
        }
    }

    output
}

fn print_peaks(
    writer: &mut Box<dyn Write>,
    peaks: &[Peak],
    config: &Config,
    unique_kmers: u64,
    hist: &[u64],
) -> Result<()> {
    if !peaks.is_empty() {
        let min_volume_fraction = single_copy_kmer_fraction(0.0003, config.k, 2).min(1.0);
        let ploidy_estimate = calc_ploidy(peaks, min_volume_fraction);
        let ploidy = if config.peak_ploidy > 0 {
            config.peak_ploidy
        } else {
            ploidy_estimate
        };
        let haploid_peak_center = haploid_peak_center(peaks, ploidy);
        let error_kmers = error_kmers(peaks, hist, min_volume_fraction);
        let genome_size_in_peaks = genome_size_in_peaks(peaks, haploid_peak_center);
        let total_genome_size = genome_size2(peaks, haploid_peak_center, hist);
        let repeat_size = repeat_size(peaks, ploidy, haploid_peak_center);
        let repeat_size2 = repeat_size2(peaks, ploidy, haploid_peak_center, hist);
        let haploid_size = total_genome_size / i64::from(ploidy);
        let het_locs = calc_het_locations(peaks, ploidy, haploid_peak_center, config.k);
        let het_rate = (het_locs as f64 / haploid_size as f64) / 2.0;
        let repeat_rate = repeat_size as f64 / genome_size_in_peaks as f64;
        let repeat_rate2 = repeat_size2 as f64 / total_genome_size as f64;

        let mut main_peak = &peaks[0];
        let mut ploidy_peak = &peaks[0];
        let target = haploid_peak_center * f64::from(ploidy);
        for peak in peaks {
            if peak.volume > main_peak.volume {
                main_peak = peak;
            }
            if abs_diff(peak.center as f64, target) < abs_diff(ploidy_peak.center as f64, target) {
                ploidy_peak = peak;
            }
        }
        let haploid_cov = if target.max(ploidy_peak.center as f64)
            / target.min(ploidy_peak.center as f64)
            < 1.3
        {
            ploidy_peak.center as i64
        } else {
            target as i64
        };

        writeln!(writer, "#k\t{}", config.k)?;
        writeln!(writer, "#unique_kmers\t{unique_kmers}")?;
        writeln!(writer, "#error_kmers\t{error_kmers}")?;
        writeln!(
            writer,
            "#genomic_kmers\t{}",
            unique_kmers as i64 - error_kmers
        )?;
        writeln!(writer, "#main_peak\t{}", main_peak.center)?;
        writeln!(writer, "#genome_size_in_peaks\t{genome_size_in_peaks}")?;
        writeln!(writer, "#genome_size\t{total_genome_size}")?;
        writeln!(writer, "#haploid_genome_size\t{haploid_size}")?;
        writeln!(
            writer,
            "#fold_coverage\t{}",
            java_round_f64(haploid_peak_center)
        )?;
        writeln!(writer, "#haploid_fold_coverage\t{haploid_cov}")?;
        writeln!(writer, "#ploidy\t{ploidy}")?;
        if ploidy != ploidy_estimate {
            writeln!(writer, "#ploidy_detected\t{ploidy_estimate}")?;
        }
        if ploidy > 1 {
            writeln!(writer, "#het_rate\t{het_rate:.5}")?;
        }
        writeln!(
            writer,
            "#percent_repeat_in_peaks\t{:.3}",
            100.0 * repeat_rate
        )?;
        if repeat_size2 >= 0 {
            writeln!(writer, "#percent_repeat\t{:.3}", 100.0 * repeat_rate2)?;
        }
    }

    writeln!(writer, "#start\tcenter\tstop\tmax\tvolume")?;
    for peak in peaks {
        if peak.volume >= config.peak_min_volume {
            writeln!(
                writer,
                "{}\t{}\t{}\t{}\t{}",
                peak.start, peak.center, peak.stop, peak.max_height, peak.volume
            )?;
        }
    }
    Ok(())
}

fn single_copy_kmer_fraction(het_rate: f64, k: usize, ploidy: i32) -> f64 {
    if ploidy < 2 {
        return 1.0;
    }
    let single_copy_kmers = het_rate * k as f64;
    let asymptote = single_copy_kmers / (1.0 + single_copy_kmers);
    asymptote * 2.0
}

fn error_kmers(peaks: &[Peak], hist: &[u64], min_volume_fraction: f64) -> i64 {
    let Some(first) = first_genomic_peak(peaks, min_volume_fraction) else {
        return 0;
    };
    hist.iter().take(first.start).sum::<u64>() as i64
}

fn first_genomic_peak(peaks: &[Peak], min_fraction: f64) -> Option<&Peak> {
    let biggest = &peaks[biggest_peak(peaks)];
    let min_volume = (biggest.volume as f64 * min_fraction) as u64;
    peaks.iter().find(|peak| peak.volume >= min_volume)
}

fn genome_size_in_peaks(peaks: &[Peak], haploid_peak_center: f64) -> i64 {
    let mult = 1.0 / haploid_peak_center.max(1.0);
    peaks
        .iter()
        .map(|peak| peak.volume as i64 * java_round_f64(peak.center as f64 * mult))
        .sum()
}

fn genome_size2(peaks: &[Peak], haploid_peak_center: f64, hist: &[u64]) -> i64 {
    let mult = 1.0 / haploid_peak_center.max(1.0);
    hist.iter()
        .enumerate()
        .skip(peaks[0].start)
        .map(|(i, &count)| count as i64 * java_round_f64(i as f64 * mult).max(1))
        .sum()
}

fn repeat_size(peaks: &[Peak], ploidy: i32, haploid_peak_center: f64) -> i64 {
    if peaks.len() < 2 {
        return 0;
    }
    let mult = 1.0 / haploid_peak_center.max(1.0);
    let homozygous_loc = homozygous_peak(peaks, ploidy, haploid_peak_center);
    peaks
        .iter()
        .skip(homozygous_loc + 1)
        .map(|peak| peak.volume as i64 * (java_round_f64(peak.center as f64 * mult) - 1))
        .sum()
}

fn repeat_size2(peaks: &[Peak], ploidy: i32, haploid_peak_center: f64, hist: &[u64]) -> i64 {
    let mult = 1.0 / haploid_peak_center.max(1.0);
    let mut valley =
        (haploid_peak_center * f64::from(ploidy) * (1.2 + 1.0 / 2.0_f64.max(f64::from(ploidy))))
            .ceil() as usize;
    let homozygous_loc = homozygous_peak(peaks, ploidy, haploid_peak_center);
    if ploidy > 1 {
        valley = peaks[homozygous_loc].stop + 1;
    }
    hist.iter()
        .enumerate()
        .skip(valley)
        .map(|(i, &count)| count as i64 * (java_round_f64(i as f64 * mult) - 1))
        .sum()
}

fn biggest_peak(peaks: &[Peak]) -> usize {
    if peaks.len() < 2 {
        return peaks.len().saturating_sub(1);
    }
    let mut loc = 0;
    let mut biggest = &peaks[0];
    for (i, peak) in peaks.iter().enumerate().skip(1) {
        if peak.volume > biggest.volume {
            loc = i;
            biggest = peak;
        }
    }
    loc
}

fn second_biggest_peak(peaks: &[Peak]) -> usize {
    if peaks.len() < 2 {
        return peaks.len().saturating_sub(1);
    }
    let mut biggest = &peaks[0];
    let mut second = &peaks[1];
    let mut bloc = 0;
    let mut sloc = 1;
    if second.volume > biggest.volume {
        std::mem::swap(&mut second, &mut biggest);
        bloc = 1;
        sloc = 0;
    }
    for (i, peak) in peaks.iter().enumerate().skip(2) {
        if peak.volume > second.volume {
            sloc = i;
            second = peak;
            if second.volume > biggest.volume {
                std::mem::swap(&mut second, &mut biggest);
                sloc = bloc;
                bloc = i;
            }
        }
    }
    sloc
}

fn homozygous_peak(peaks: &[Peak], ploidy: i32, haploid_peak_center: f64) -> usize {
    if peaks.len() < 2 {
        return peaks.len().saturating_sub(1);
    }
    let target = haploid_peak_center * f64::from(ploidy);
    let mut best_diff = f64::from(i32::MAX);
    let mut loc = 0;
    for (i, peak) in peaks.iter().enumerate() {
        let diff = abs_diff(target, peak.center as f64);
        if diff < best_diff {
            best_diff = diff;
            loc = i;
        }
    }
    loc
}

fn haploid_peak_center(peaks: &[Peak], ploidy: i32) -> f64 {
    let biggest = &peaks[biggest_peak(peaks)];
    let second = &peaks[second_biggest_peak(peaks)];
    if second.volume.saturating_mul(4) >= biggest.volume {
        biggest.center.min(second.center) as f64
    } else {
        biggest.center as f64 / f64::from(ploidy)
    }
}

fn calc_ploidy(peaks: &[Peak], min_volume_fraction: f64) -> i32 {
    if peaks.len() < 2 {
        return 1;
    }
    let biggest = &peaks[biggest_peak(peaks)];
    let second = &peaks[second_biggest_peak(peaks)];

    if std::ptr::eq(second, biggest) {
        return 1;
    }
    if second.center < biggest.center {
        if (second.volume as f64) < (biggest.volume as f64) * min_volume_fraction {
            return 1;
        }
    } else if second.volume.saturating_mul(4) < biggest.volume {
        return 1;
    }

    let max = biggest.center.max(second.center);
    let min = biggest.center.min(second.center);
    java_round_f64(max as f64 / min as f64).max(1) as i32
}

fn calc_het_locations(peaks: &[Peak], ploidy: i32, haploid_peak_center: f64, k: usize) -> i64 {
    if peaks.len() < 2 {
        return 0;
    }
    let homozygous_loc = homozygous_peak(peaks, ploidy, haploid_peak_center);
    let homo_peak = &peaks[homozygous_loc];
    let mut sum = 0i64;
    let lim = ploidy / 2;
    for peak in peaks.iter().take(homozygous_loc) {
        let copy_count =
            java_round_f64((peak.center as f64 * f64::from(ploidy)) / homo_peak.center as f64);
        if copy_count > i64::from(lim) {
            break;
        }
        sum += peak.volume as i64;
    }
    sum / k as i64
}

fn smooth_progressive(data: &[u64], radius0: i32) -> Vec<u64> {
    let mut radius = radius0;
    let mut div = i64::from(radius * radius);
    let mut mult = 1.0 / div as f64;
    let mut smoothed = Vec::with_capacity(data.len());
    let mut next = 5i32;
    for i in 0..data.len() {
        let sum = sum_point(data, i, radius);
        smoothed.push(java_round_f64(sum as f64 * mult) as u64);
        if i as i32 > next {
            next = (1.0 + f64::from(next) * PROGRESSIVE_MULT).ceil() as i32;
            radius += 1;
            div = i64::from(radius * radius);
            mult = 1.0 / div as f64;
            if radius > MAX_RADIUS {
                next = i32::MAX;
            }
        }
    }
    smoothed
}

fn sum_point(data: &[u64], loc: usize, radius: i32) -> u64 {
    let mut sum = 0u64;
    let start = loc as i32 - radius + 1;
    let stop = loc as i32 + radius - 1;
    for (x, i) in (start..loc as i32).enumerate() {
        let i2 = i.max(0) as usize;
        sum = sum.saturating_add(data[i2].saturating_mul((x + 1) as u64));
    }
    let max = data.len().saturating_sub(1) as i32;
    let mut x = radius;
    for i in loc as i32..=stop {
        let i2 = i.min(max) as usize;
        sum = sum.saturating_add(data[i2].saturating_mul(x as u64));
        x -= 1;
    }
    sum
}

fn unique_from_raw(depth: usize, raw: u64) -> u64 {
    if depth < 1 {
        raw
    } else {
        (raw + (depth as u64 / 2)) / depth as u64
    }
}

fn java_round_f64(value: f64) -> i64 {
    (value + 0.5).floor() as i64
}

fn abs_diff(a: f64, b: f64) -> f64 {
    (a - b).abs()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn representative_histogram_calls_java_shaped_peak() {
        let config = Config {
            k: 5,
            peak_min_height: 1,
            peak_min_volume: 1,
            peak_min_width: 1,
            peak_min_peak: 1,
            peak_max_peak: 100,
            peak_max_count: 8,
            ..Default::default()
        };
        let mut raw = vec![0; 128];
        raw[20] = 720;
        let unique: Vec<u64> = raw
            .iter()
            .copied()
            .enumerate()
            .map(|(depth, raw)| unique_from_raw(depth, raw))
            .collect();
        let peaks = call_peaks(&unique, &config);
        assert_eq!(peaks.len(), 1);
        assert_eq!(peaks[0].start, 17);
        assert_eq!(peaks[0].center, 20);
        assert_eq!(peaks[0].stop, 23);
        assert_eq!(peaks[0].max_height, 36);
        assert_eq!(peaks[0].volume, 36);
    }
}

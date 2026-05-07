use crate::cli::Config;
use crate::seqio::SequenceRecord;
use std::sync::OnceLock;

const INLINE_LONG_WORDS: usize = 8;

#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum KmerKey {
    Short(u64),
    LongHash(u64),
}

pub fn kmers_for_record(record: &SequenceRecord, config: &Config) -> Vec<KmerKey> {
    let mut out = Vec::new();
    for_each_kmer_for_record(record, config, |kmer| out.push(kmer));
    out
}

pub fn for_each_kmer_for_record<F>(record: &SequenceRecord, config: &Config, f: F)
where
    F: FnMut(KmerKey),
{
    for_each_kmer(
        &record.bases,
        record.qualities.as_deref(),
        config.k,
        config.min_quality,
        config.min_prob,
        bbnorm_count_canonical(config),
        f,
    );
}

pub fn unfiltered_kmer_windows_for_record(
    record: &SequenceRecord,
    config: &Config,
) -> Vec<Option<KmerKey>> {
    unfiltered_kmer_windows(&record.bases, config.k, bbnorm_count_canonical(config))
}

fn bbnorm_count_canonical(config: &Config) -> bool {
    config.canonical || config.k <= 31
}

pub fn kmers(
    bases: &[u8],
    qualities: Option<&[u8]>,
    k: usize,
    min_quality: u8,
    min_prob: f64,
    canonical: bool,
) -> Vec<KmerKey> {
    let mut out = Vec::new();
    for_each_kmer(
        bases,
        qualities,
        k,
        min_quality,
        min_prob,
        canonical,
        |kmer| out.push(kmer),
    );
    out
}

pub fn for_each_kmer<F>(
    bases: &[u8],
    qualities: Option<&[u8]>,
    k: usize,
    min_quality: u8,
    min_prob: f64,
    canonical: bool,
    mut f: F,
) where
    F: FnMut(KmerKey),
{
    if k == 0 || bases.len() < k {
        return;
    }
    if k > 31 {
        for_each_long_kmer(bases, qualities, k, min_quality, min_prob, canonical, f);
        return;
    }

    let mask = short_kmer_mask(k);
    let shift2 = 2 * (k - 1);
    let mut forward = 0u64;
    let mut reverse = 0u64;
    let mut len = 0usize;
    let mut prob = 1.0f64;
    let mut zero_prob_count = 0usize;
    let track_probability = qualities.is_some() && min_prob > 0.0;

    for (index, &base) in bases.iter().enumerate() {
        let Some(bits) = base_bits(base) else {
            len = 0;
            forward = 0;
            reverse = 0;
            prob = 1.0;
            zero_prob_count = 0;
            continue;
        };

        let q = qualities
            .and_then(|qual| qual.get(index).copied())
            .map(|quality| quality.saturating_sub(33))
            .unwrap_or(50);
        if q < min_quality {
            len = 0;
            forward = 0;
            reverse = 0;
            prob = 1.0;
            zero_prob_count = 0;
            continue;
        }

        forward = ((forward << 2) | u64::from(bits)) & mask;
        reverse = ((reverse >> 2) | (u64::from(3 - bits) << shift2)) & mask;

        if track_probability {
            add_quality_probability(q, &mut prob, &mut zero_prob_count);
            if len >= k {
                let old_q = qualities
                    .and_then(|qual| qual.get(index - k).copied())
                    .map(|quality| quality.saturating_sub(33))
                    .unwrap_or(50);
                remove_quality_probability(old_q, &mut prob, &mut zero_prob_count);
            }
        }

        len += 1;
        if len >= k && (!track_probability || window_probability(prob, zero_prob_count) >= min_prob)
        {
            let code = if canonical {
                forward.max(reverse)
            } else {
                forward
            };
            f(KmerKey::Short(code));
        }
    }
}

pub fn unfiltered_kmer_windows(bases: &[u8], k: usize, canonical: bool) -> Vec<Option<KmerKey>> {
    if k == 0 || bases.len() < k {
        return Vec::new();
    }
    if k > 31 {
        return unfiltered_long_kmer_windows(bases, k, canonical);
    }

    let mut out = Vec::with_capacity(bases.len() - k + 1);
    let mask = short_kmer_mask(k);
    let shift2 = 2 * (k - 1);
    let mut forward = 0u64;
    let mut reverse = 0u64;
    let mut len = 0usize;
    for (index, &base) in bases.iter().enumerate() {
        if let Some(bits) = base_bits(base) {
            forward = ((forward << 2) | u64::from(bits)) & mask;
            reverse = ((reverse >> 2) | (u64::from(3 - bits) << shift2)) & mask;
            len += 1;
        } else {
            forward = 0;
            reverse = 0;
            len = 0;
        }
        if index + 1 < k {
            continue;
        }
        out.push((len >= k).then(|| {
            let code = if canonical {
                forward.max(reverse)
            } else {
                forward
            };
            KmerKey::Short(code)
        }));
    }
    out
}

fn for_each_long_kmer<F>(
    bases: &[u8],
    qualities: Option<&[u8]>,
    k: usize,
    min_quality: u8,
    min_prob: f64,
    _canonical: bool,
    mut f: F,
) where
    F: FnMut(KmerKey),
{
    let mut roller = LongKmerRoller::new(k);
    let mut len = 0usize;
    let mut prob = 1.0f64;
    let mut zero_prob_count = 0usize;
    let track_quality = qualities.is_some();
    let track_probability = track_quality && min_prob > 0.0;

    for (index, &base) in bases.iter().enumerate() {
        let Some(bits) = base_bits(base) else {
            roller.clear();
            len = 0;
            prob = 1.0;
            zero_prob_count = 0;
            continue;
        };

        let q = qualities
            .and_then(|qual| qual.get(index).copied())
            .map(|quality| quality.saturating_sub(33))
            .unwrap_or(50);
        if track_quality && q < min_quality {
            roller.clear();
            len = 0;
            prob = 1.0;
            zero_prob_count = 0;
            continue;
        }

        roller.add_right_bits(bits);

        if track_probability {
            add_quality_probability(q, &mut prob, &mut zero_prob_count);
            if len >= k {
                let old_q = qualities
                    .and_then(|qual| qual.get(index - k).copied())
                    .map(|quality| quality.saturating_sub(33))
                    .unwrap_or(50);
                remove_quality_probability(old_q, &mut prob, &mut zero_prob_count);
            }
        }

        len += 1;
        if len >= k && (!track_probability || window_probability(prob, zero_prob_count) >= min_prob)
        {
            f(KmerKey::LongHash(roller.xor_key()));
        }
    }
}

fn unfiltered_long_kmer_windows(bases: &[u8], k: usize, _canonical: bool) -> Vec<Option<KmerKey>> {
    let mut out = Vec::with_capacity(bases.len() - k + 1);
    let mut roller = LongKmerRoller::new(k);
    let mut len = 0usize;
    for (index, &base) in bases.iter().enumerate() {
        if let Some(bits) = base_bits(base) {
            roller.add_right_bits(bits);
            len += 1;
        } else {
            roller.clear();
            len = 0;
        }
        if index + 1 < k {
            continue;
        }
        out.push((len >= k).then(|| KmerKey::LongHash(roller.xor_key())));
    }
    out
}

struct LongKmerRoller {
    mult: usize,
    shift2: u32,
    mask: u64,
    forward_inline: [u64; INLINE_LONG_WORDS],
    reverse_inline: [u64; INLINE_LONG_WORDS],
    forward_heap: Vec<u64>,
    reverse_heap: Vec<u64>,
}

impl LongKmerRoller {
    fn new(requested_k: usize) -> Self {
        let (word_len, mult, _) = java_long_layout(requested_k);
        let (forward_heap, reverse_heap) = if mult <= INLINE_LONG_WORDS {
            (Vec::new(), Vec::new())
        } else {
            (vec![0; mult], vec![0; mult])
        };
        Self {
            mult,
            shift2: (2 * (word_len - 1)) as u32,
            mask: short_kmer_mask(word_len),
            forward_inline: [0; INLINE_LONG_WORDS],
            reverse_inline: [0; INLINE_LONG_WORDS],
            forward_heap,
            reverse_heap,
        }
    }

    fn clear(&mut self) {
        if self.mult <= INLINE_LONG_WORDS {
            self.forward_inline[..self.mult].fill(0);
            self.reverse_inline[..self.mult].fill(0);
        } else {
            self.forward_heap.fill(0);
            self.reverse_heap.fill(0);
        }
    }

    fn add_right_bits(&mut self, bits: u8) {
        let mut x = u64::from(bits);
        let mut x2 = u64::from(3 - bits);
        let mult = self.mult;
        let shift2 = self.shift2;
        let mask = self.mask;
        let (forward, reverse) = self.words_mut();

        for (j, reverse_word) in reverse.iter_mut().enumerate() {
            let i = mult - 1 - j;
            let y = (forward[i] >> shift2) & 3;
            let y2 = *reverse_word & 3;

            forward[i] = ((forward[i] << 2) | x) & mask;
            *reverse_word = ((*reverse_word >> 2) | (x2 << shift2)) & mask;

            x = y;
            x2 = y2;
        }
    }

    fn xor_key(&self) -> u64 {
        let (forward, reverse) = self.words();
        if java_words_less(forward, reverse) {
            java_words_xor(reverse)
        } else {
            java_words_xor(forward)
        }
    }

    fn words(&self) -> (&[u64], &[u64]) {
        if self.mult <= INLINE_LONG_WORDS {
            (
                &self.forward_inline[..self.mult],
                &self.reverse_inline[..self.mult],
            )
        } else {
            (&self.forward_heap, &self.reverse_heap)
        }
    }

    fn words_mut(&mut self) -> (&mut [u64], &mut [u64]) {
        if self.mult <= INLINE_LONG_WORDS {
            (
                &mut self.forward_inline[..self.mult],
                &mut self.reverse_inline[..self.mult],
            )
        } else {
            (&mut self.forward_heap, &mut self.reverse_heap)
        }
    }
}

#[cfg(test)]
fn java_long_kmer_hash_window(bases: &[u8], start: usize, k: usize) -> u64 {
    let (word_len, mult, effective_len) = java_long_layout(k);
    let effective_start = start + k - effective_len;

    if mult <= INLINE_LONG_WORDS {
        let mut forward = [0u64; INLINE_LONG_WORDS];
        let mut reverse = [0u64; INLINE_LONG_WORDS];
        let forward = &mut forward[..mult];
        let reverse = &mut reverse[..mult];
        fill_java_long_words_forward(bases, effective_start, word_len, forward);
        fill_java_long_words_reverse_complement(bases, effective_start, word_len, reverse);
        if java_words_less(forward, reverse) {
            java_words_xor(reverse)
        } else {
            java_words_xor(forward)
        }
    } else {
        let mut forward = vec![0u64; mult];
        let mut reverse = vec![0u64; mult];
        fill_java_long_words_forward(bases, effective_start, word_len, &mut forward);
        fill_java_long_words_reverse_complement(bases, effective_start, word_len, &mut reverse);
        if java_words_less(&forward, &reverse) {
            java_words_xor(&reverse)
        } else {
            java_words_xor(&forward)
        }
    }
}

fn java_long_layout(requested_k: usize) -> (usize, usize, usize) {
    let mult = java_long_mult(requested_k);
    let word_len = requested_k / mult;
    (word_len, mult, word_len * mult)
}

fn java_long_mult(kbig: usize) -> usize {
    let word = 31usize;
    let mult1 = kbig.div_ceil(word);
    let mult2 = (kbig / word).max(1);
    if mult1 == mult2 {
        return mult1;
    }
    let k1 = word.min(kbig / mult1);
    let k2 = word.min(kbig / mult2);
    let kbig1 = k1 * mult1;
    let kbig2 = k2 * mult2;
    if kbig2 >= kbig1 { mult2 } else { mult1 }
}

#[cfg(test)]
fn fill_java_long_words_forward(bases: &[u8], start: usize, word_len: usize, words: &mut [u64]) {
    for (word_index, word_ref) in words.iter_mut().enumerate() {
        let mut word = 0u64;
        let offset = start + word_index * word_len;
        for &base in bases.iter().skip(offset).take(word_len) {
            let base = base_bits(base).unwrap_or(0);
            word = (word << 2) | u64::from(base);
        }
        *word_ref = word;
    }
}

#[cfg(test)]
fn fill_java_long_words_reverse_complement(
    bases: &[u8],
    start: usize,
    word_len: usize,
    words: &mut [u64],
) {
    let mult = words.len();
    let effective_len = word_len * mult;
    for (word_index, word_ref) in words.iter_mut().enumerate() {
        let mut word = 0u64;
        let offset = word_index * word_len;
        for idx in offset..offset + word_len {
            let base = base_bits(bases[start + effective_len - 1 - idx]).unwrap_or(0);
            word = (word << 2) | u64::from(3 - base);
        }
        *word_ref = word;
    }
}

fn java_words_less(left: &[u64], right: &[u64]) -> bool {
    for (&left_word, &right_word) in left.iter().zip(right) {
        if left_word < right_word {
            return true;
        }
        if left_word > right_word {
            return false;
        }
    }
    false
}

fn java_words_xor(words: &[u64]) -> u64 {
    let mut xor = words[0];
    for &word in &words[1..] {
        xor = xor.rotate_left(25) ^ word;
    }
    xor & i64::MAX as u64
}

fn base_bits(base: u8) -> Option<u8> {
    match base {
        b'A' | b'a' => Some(0),
        b'C' | b'c' => Some(1),
        b'G' | b'g' => Some(2),
        b'T' | b't' | b'U' | b'u' => Some(3),
        _ => None,
    }
}

fn base_correct_probability(q: u8) -> f64 {
    1.0 - 10f64.powf(-(f64::from(q)) / 10.0)
}

fn quality_probability_table() -> &'static [f64; 256] {
    static TABLE: OnceLock<[f64; 256]> = OnceLock::new();
    TABLE.get_or_init(|| {
        let mut table = [0.0; 256];
        for (q, probability) in table.iter_mut().enumerate() {
            *probability = base_correct_probability(q as u8);
        }
        table
    })
}

fn add_quality_probability(q: u8, prob: &mut f64, zero_prob_count: &mut usize) {
    let probability = quality_probability_table()[q as usize];
    if probability == 0.0 {
        *zero_prob_count += 1;
    } else {
        *prob *= probability;
    }
}

fn remove_quality_probability(q: u8, prob: &mut f64, zero_prob_count: &mut usize) {
    let probability = quality_probability_table()[q as usize];
    if probability == 0.0 {
        *zero_prob_count = zero_prob_count.saturating_sub(1);
    } else {
        *prob /= probability;
    }
}

fn window_probability(prob: f64, zero_prob_count: usize) -> f64 {
    if zero_prob_count == 0 { prob } else { 0.0 }
}

fn short_kmer_mask(k: usize) -> u64 {
    if k >= 32 {
        u64::MAX
    } else {
        (1u64 << (2 * k)) - 1
    }
}

pub(crate) fn canonical_short_code(code: u64, k: usize) -> u64 {
    let rc = reverse_complement_code(code, k);
    code.max(rc)
}

fn reverse_complement_code(mut code: u64, k: usize) -> u64 {
    let mut rc = 0u64;
    for _ in 0..k {
        let base = (!code) & 0b11;
        rc = (rc << 2) | base;
        code >>= 2;
    }
    rc
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::cli::Config;
    use crate::seqio::SequenceRecord;

    #[test]
    fn canonicalizes_reverse_complements() {
        let cfg = Config {
            k: 4,
            min_quality: 0,
            min_prob: 0.0,
            canonical: true,
            ..Config::default()
        };
        let a = SequenceRecord {
            id: "a".to_string(),
            numeric_id: 0,
            bases: b"ACGT".to_vec(),
            qualities: None,
        };
        let b = SequenceRecord {
            id: "b".to_string(),
            numeric_id: 1,
            bases: b"ACGT".to_vec(),
            qualities: None,
        };
        assert_eq!(kmers_for_record(&a, &cfg), kmers_for_record(&b, &cfg));
    }

    #[test]
    fn bbnorm_short_kmers_remain_canonical_when_canonical_flag_is_false() {
        let cfg = Config {
            k: 4,
            min_quality: 0,
            min_prob: 0.0,
            canonical: false,
            ..Config::default()
        };
        let a = SequenceRecord {
            id: "a".to_string(),
            numeric_id: 0,
            bases: b"ACGA".to_vec(),
            qualities: Some(vec![b'I'; 4]),
        };
        let b = SequenceRecord {
            id: "b".to_string(),
            numeric_id: 0,
            bases: b"TCGT".to_vec(),
            qualities: Some(vec![b'I'; 4]),
        };
        assert_eq!(kmers_for_record(&a, &cfg), kmers_for_record(&b, &cfg));
    }

    #[test]
    fn canonicalizes_long_reverse_complements() {
        let cfg = Config {
            k: 32,
            min_quality: 0,
            min_prob: 0.0,
            canonical: true,
            ..Config::default()
        };
        let a = SequenceRecord {
            id: "a".to_string(),
            numeric_id: 0,
            bases: b"ACGTACGTACGTACGTACGTACGTACGTACGA".to_vec(),
            qualities: None,
        };
        let b = SequenceRecord {
            id: "b".to_string(),
            numeric_id: 1,
            bases: b"TCGTACGTACGTACGTACGTACGTACGTACGT".to_vec(),
            qualities: None,
        };
        assert_eq!(kmers_for_record(&a, &cfg), kmers_for_record(&b, &cfg));
    }

    #[test]
    fn long_kmer_hash_matches_allocating_java_word_reference() {
        let bases = b"ACGTGCAATTCGCGATATCGGATCCGATTAACCGGTACGTTAGCATCGATCGGCTAGCTAGCTT".repeat(5);
        for &k in &[32usize, 40, 63, 124, 248, 280] {
            for &start in &[0usize, 3, 17] {
                assert!(start + k <= bases.len());
                assert_eq!(
                    java_long_kmer_hash_window(&bases, start, k),
                    slow_java_long_kmer_hash_window(&bases, start, k),
                    "k={k}, start={start}"
                );
            }
        }
    }

    #[test]
    fn rolling_long_kmers_match_window_scan_with_quality_resets() {
        let mut bases =
            b"ACGTGCAATTCGCGATATCGGATCCGATTAACCGGTACGTTAGCATCGATCGGCTAGCTAGCTT".repeat(3);
        bases[47] = b'N';
        let mut qualities = vec![b'I'; bases.len()];
        qualities[21] = b'"';
        qualities[91] = b'#';

        let rolling = kmers(&bases, Some(&qualities), 40, 2, 0.90, true);
        let scanned = slow_long_kmers(&bases, Some(&qualities), 40, 2, 0.90);

        assert_eq!(rolling, scanned);
    }

    #[test]
    fn rolling_unfiltered_long_windows_match_window_scan() {
        let mut bases =
            b"ACGTGCAATTCGCGATATCGGATCCGATTAACCGGTACGTTAGCATCGATCGGCTAGCTAGCTT".repeat(2);
        bases[7] = b'N';
        bases[72] = b'N';

        assert_eq!(
            unfiltered_kmer_windows(&bases, 40, true),
            slow_unfiltered_long_windows(&bases, 40)
        );
    }

    #[test]
    fn filters_low_quality_windows() {
        let cfg = Config {
            k: 3,
            min_quality: 30,
            min_prob: 0.0,
            ..Config::default()
        };
        let record = SequenceRecord {
            id: "r".to_string(),
            numeric_id: 0,
            bases: b"ACGT".to_vec(),
            qualities: Some(b"I!II".to_vec()),
        };
        assert!(kmers_for_record(&record, &cfg).is_empty());
    }

    #[test]
    fn unfiltered_windows_preserve_invalid_positions() {
        let windows = unfiltered_kmer_windows(b"ACNT", 2, false);
        assert_eq!(windows.len(), 3);
        assert!(windows[0].is_some());
        assert!(windows[1].is_none());
        assert!(windows[2].is_none());
    }

    #[test]
    fn unfiltered_windows_ignore_quality_filters() {
        let cfg = Config {
            k: 3,
            min_quality: 30,
            min_prob: 0.99,
            ..Config::default()
        };
        let record = SequenceRecord {
            id: "r".to_string(),
            numeric_id: 0,
            bases: b"ACGT".to_vec(),
            qualities: Some(b"!!!!".to_vec()),
        };
        assert!(kmers_for_record(&record, &cfg).is_empty());
        assert_eq!(unfiltered_kmer_windows_for_record(&record, &cfg).len(), 2);
        assert!(
            unfiltered_kmer_windows_for_record(&record, &cfg)
                .into_iter()
                .all(|window| window.is_some())
        );
    }

    #[test]
    fn rolling_short_kmers_match_window_scan_with_quality_resets() {
        let bases = b"ACGTNACGTACGTT";
        let mut qualities = vec![b'I'; bases.len()];
        qualities[6] = b'"';
        qualities[10] = b'#';

        let rolling = kmers(bases, Some(&qualities), 4, 2, 0.90, true);
        let scanned = slow_short_kmers(bases, Some(&qualities), 4, 2, 0.90, true);

        assert_eq!(rolling, scanned);
    }

    #[test]
    fn rolling_unfiltered_short_windows_match_window_scan() {
        let bases = b"NACGTNACGTA";
        assert_eq!(
            unfiltered_kmer_windows(bases, 3, true),
            slow_unfiltered_short_windows(bases, 3, true)
        );
    }

    fn slow_short_kmers(
        bases: &[u8],
        qualities: Option<&[u8]>,
        k: usize,
        min_quality: u8,
        min_prob: f64,
        canonical: bool,
    ) -> Vec<KmerKey> {
        let mut out = Vec::new();
        for start in 0..=bases.len() - k {
            let mut code = 0u64;
            let mut prob = 1.0f64;
            let mut valid = true;
            for offset in 0..k {
                let idx = start + offset;
                let Some(bits) = base_bits(bases[idx]) else {
                    valid = false;
                    break;
                };
                let q = qualities
                    .and_then(|qual| qual.get(idx).copied())
                    .map(|quality| quality.saturating_sub(33))
                    .unwrap_or(50);
                if q < min_quality {
                    valid = false;
                    break;
                }
                if qualities.is_some() && min_prob > 0.0 {
                    prob *= base_correct_probability(q);
                    if prob < min_prob {
                        valid = false;
                        break;
                    }
                }
                code = (code << 2) | u64::from(bits);
            }
            if valid {
                out.push(KmerKey::Short(if canonical {
                    canonical_short_code(code, k)
                } else {
                    code
                }));
            }
        }
        out
    }

    fn slow_unfiltered_short_windows(
        bases: &[u8],
        k: usize,
        canonical: bool,
    ) -> Vec<Option<KmerKey>> {
        let mut out = Vec::new();
        for start in 0..=bases.len() - k {
            let mut code = 0u64;
            let mut valid = true;
            for offset in 0..k {
                let Some(bits) = base_bits(bases[start + offset]) else {
                    valid = false;
                    break;
                };
                code = (code << 2) | u64::from(bits);
            }
            out.push(valid.then(|| {
                KmerKey::Short(if canonical {
                    canonical_short_code(code, k)
                } else {
                    code
                })
            }));
        }
        out
    }

    fn slow_java_long_kmer_hash_window(bases: &[u8], start: usize, k: usize) -> u64 {
        let (word_len, mult, effective_len) = java_long_layout(k);
        let effective_start = start + k - effective_len;
        let forward = slow_encode_java_long_words_forward(bases, effective_start, word_len, mult);
        let reverse =
            slow_encode_java_long_words_reverse_complement(bases, effective_start, word_len, mult);
        if java_words_less(&forward, &reverse) {
            java_words_xor(&reverse)
        } else {
            java_words_xor(&forward)
        }
    }

    fn slow_long_kmers(
        bases: &[u8],
        qualities: Option<&[u8]>,
        k: usize,
        min_quality: u8,
        min_prob: f64,
    ) -> Vec<KmerKey> {
        let mut out = Vec::new();
        for start in 0..=bases.len() - k {
            let mut prob = 1.0f64;
            let mut valid = true;
            for offset in 0..k {
                let idx = start + offset;
                if base_bits(bases[idx]).is_none() {
                    valid = false;
                    break;
                }
                if let Some(qual) = qualities {
                    let q = qual
                        .get(idx)
                        .copied()
                        .map(|quality| quality.saturating_sub(33))
                        .unwrap_or(50);
                    if q < min_quality {
                        valid = false;
                        break;
                    }
                    if min_prob > 0.0 {
                        prob *= base_correct_probability(q);
                        if prob < min_prob {
                            valid = false;
                            break;
                        }
                    }
                }
            }
            if valid {
                out.push(KmerKey::LongHash(java_long_kmer_hash_window(
                    bases, start, k,
                )));
            }
        }
        out
    }

    fn slow_unfiltered_long_windows(bases: &[u8], k: usize) -> Vec<Option<KmerKey>> {
        let mut out = Vec::new();
        for start in 0..=bases.len() - k {
            let mut valid = true;
            for offset in 0..k {
                if base_bits(bases[start + offset]).is_none() {
                    valid = false;
                    break;
                }
            }
            out.push(valid.then(|| KmerKey::LongHash(java_long_kmer_hash_window(bases, start, k))));
        }
        out
    }

    fn slow_encode_java_long_words_forward(
        bases: &[u8],
        start: usize,
        word_len: usize,
        mult: usize,
    ) -> Vec<u64> {
        let mut words = Vec::with_capacity(mult);
        for word_index in 0..mult {
            let mut word = 0u64;
            let offset = start + word_index * word_len;
            for &base in bases.iter().skip(offset).take(word_len) {
                let base = base_bits(base).unwrap_or(0);
                word = (word << 2) | u64::from(base);
            }
            words.push(word);
        }
        words
    }

    fn slow_encode_java_long_words_reverse_complement(
        bases: &[u8],
        start: usize,
        word_len: usize,
        mult: usize,
    ) -> Vec<u64> {
        let mut words = Vec::with_capacity(mult);
        let effective_len = word_len * mult;
        for word_index in 0..mult {
            let mut word = 0u64;
            let offset = word_index * word_len;
            for idx in offset..offset + word_len {
                let base = base_bits(bases[start + effective_len - 1 - idx]).unwrap_or(0);
                word = (word << 2) | u64::from(3 - base);
            }
            words.push(word);
        }
        words
    }
}

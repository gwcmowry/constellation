use serde::{Deserialize, Serialize};

#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
#[repr(u8)]
pub enum Base {
    A = 0,
    C = 1,
    G = 2,
    T = 3,
}

impl Base {
    pub fn from_ascii(b: u8) -> Option<Self> {
        match b.to_ascii_uppercase() {
            b'A' => Some(Self::A),
            b'C' => Some(Self::C),
            b'G' => Some(Self::G),
            b'T' | b'U' => Some(Self::T),
            _ => None,
        }
    }

    pub fn to_ascii(self) -> u8 {
        match self {
            Self::A => b'A',
            Self::C => b'C',
            Self::G => b'G',
            Self::T => b'T',
        }
    }

    pub fn complement(self) -> Self {
        match self {
            Self::A => Self::T,
            Self::C => Self::G,
            Self::G => Self::C,
            Self::T => Self::A,
        }
    }
}

#[derive(Debug, Clone, Default, PartialEq, Eq, Serialize, Deserialize)]
pub struct EncodedSeq {
    pub words: Vec<u64>,
    pub len: u32,
    pub n_mask_words: Vec<u64>,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Kmer {
    pub code: u64,
    pub pos: u32,
    pub k: u8,
}

pub fn encode_acgt(seq: &[u8]) -> EncodedSeq {
    let mut encoded = EncodedSeq::default();
    encode_acgt_into(seq, &mut encoded);
    encoded
}

/// Encode a sequence while retaining the buffers across reads of any length.
pub fn encode_acgt_into(seq: &[u8], encoded: &mut EncodedSeq) {
    encoded.len = seq.len() as u32;
    let word_count = seq.len().div_ceil(32);
    let mask_word_count = seq.len().div_ceil(64);
    encoded.words.resize(word_count, 0);
    encoded.words.fill(0);
    encoded.n_mask_words.resize(mask_word_count, 0);
    encoded.n_mask_words.fill(0);
    let words = &mut encoded.words;
    let n_mask_words = &mut encoded.n_mask_words;
    let mut has_n = false;

    for (i, &base) in seq.iter().enumerate() {
        let word_idx = i / 32;
        let shift = (i % 32) * 2;
        match Base::from_ascii(base) {
            Some(b) => words[word_idx] |= (b as u64) << shift,
            None => {
                has_n = true;
                n_mask_words[i / 64] |= 1_u64 << (i % 64);
            }
        }
    }

    if !has_n {
        n_mask_words.clear();
    }
}

const BASE_CODES: [u8; 256] = {
    let mut codes = [4; 256];
    codes[b'A' as usize] = 0;
    codes[b'a' as usize] = 0;
    codes[b'C' as usize] = 1;
    codes[b'c' as usize] = 1;
    codes[b'G' as usize] = 2;
    codes[b'g' as usize] = 2;
    codes[b'T' as usize] = 3;
    codes[b't' as usize] = 3;
    codes[b'U' as usize] = 3;
    codes[b'u' as usize] = 3;
    codes
};

/// A/C/G/T/U code, or 4 for an invalid base.
#[inline]
pub fn acgt_code(base: u8) -> u8 {
    BASE_CODES[base as usize]
}

/// Encode and compute the mapper's exact low-complexity predicate in one pass.
/// U encodes as T but is excluded from A/C/G/T composition, as in the sketcher.
pub fn encode_and_check_complexity(seq: &[u8], k: u8, encoded: &mut EncodedSeq) -> bool {
    encoded.len = seq.len() as u32;
    encoded.words.resize(seq.len().div_ceil(32), 0);
    encoded.words.fill(0);
    encoded.n_mask_words.resize(seq.len().div_ceil(64), 0);
    encoded.n_mask_words.fill(0);
    let mut counts = [0_usize; 5];
    let mut valid_run = 0_usize;
    let mut kmer_count = 0_usize;
    let mut has_n = false;
    for (i, &base) in seq.iter().enumerate() {
        let code = acgt_code(base);
        let composition_code = if base.eq_ignore_ascii_case(&b'U') {
            4
        } else {
            code
        };
        counts[composition_code as usize] += 1;
        encoded.words[i / 32] |= ((code & 3) as u64) << ((i % 32) * 2);
        if code < 4 {
            valid_run += 1;
            kmer_count += usize::from(k > 0 && valid_run >= k as usize);
        } else {
            valid_run = 0;
            has_n = true;
            encoded.n_mask_words[i / 64] |= 1_u64 << (i % 64);
        }
    }
    if !has_n {
        encoded.n_mask_words.clear();
    }
    let valid_bases = counts[..4].iter().sum::<usize>();
    kmer_count <= 1
        || valid_bases == 0
        || counts[..4].iter().copied().max().unwrap_or(0) * 5 >= valid_bases * 4
        || counts[..4].iter().filter(|&&n| n > 0).count() <= 1
}

pub fn base_at(encoded: &EncodedSeq, pos: usize) -> Option<Base> {
    if pos >= encoded.len as usize {
        return None;
    }
    if !encoded.n_mask_words.is_empty()
        && (encoded.n_mask_words[pos / 64] & (1_u64 << (pos % 64))) != 0
    {
        return None;
    }
    let code = (encoded.words[pos / 32] >> ((pos % 32) * 2)) & 0b11;
    match code {
        0 => Some(Base::A),
        1 => Some(Base::C),
        2 => Some(Base::G),
        3 => Some(Base::T),
        _ => unreachable!(),
    }
}

pub fn decode_to_vec(encoded: &EncodedSeq) -> Vec<u8> {
    (0..encoded.len as usize)
        .map(|pos| base_at(encoded, pos).map_or(b'N', Base::to_ascii))
        .collect()
}

pub fn reverse_complement(encoded: &EncodedSeq) -> EncodedSeq {
    let mut out = Vec::with_capacity(encoded.len as usize);
    for pos in (0..encoded.len as usize).rev() {
        let b = base_at(encoded, pos)
            .map(|base| base.complement().to_ascii())
            .unwrap_or(b'N');
        out.push(b);
    }
    encode_acgt(&out)
}

pub fn iter_kmers_2bit(encoded: &EncodedSeq, k: u8) -> KmerIter<'_> {
    KmerIter {
        encoded,
        k,
        pos: 0,
        code: 0,
        valid_run: 0,
    }
}

/// Iterates closed syncmers in one pass without allocating per sequence.
///
/// A k-mer is selected when the smallest mixed canonical s-mer order value
/// within it occurs at either end of the k-mer. Ties at an endpoint are
/// selected. The rule is sequence- and strand-intrinsic, so the same k-mer is
/// selected in a reference and query regardless of its absolute position.
pub fn iter_closed_syncmers_2bit(encoded: &EncodedSeq, k: u8, s: u8) -> ClosedSyncmerIter<'_> {
    ClosedSyncmerIter {
        encoded,
        k,
        s,
        pos: 0,
        kmer_code: 0,
        smer_code: 0,
        smer_rc_code: 0,
        valid_run: 0,
        deque_positions: [0; 33],
        deque_orders: [0; 33],
        deque_head: 0,
        deque_len: 0,
    }
}

pub fn hamming_2bit(a: &EncodedSeq, b: &EncodedSeq) -> u32 {
    let shared = a.len.min(b.len) as usize;
    let mut mismatches = a.len.abs_diff(b.len);
    for pos in 0..shared {
        if base_at(a, pos) != base_at(b, pos) {
            mismatches += 1;
        }
    }
    mismatches
}

pub struct KmerIter<'a> {
    encoded: &'a EncodedSeq,
    k: u8,
    pos: usize,
    code: u64,
    valid_run: u8,
}

impl Iterator for KmerIter<'_> {
    type Item = Kmer;

    fn next(&mut self) -> Option<Self::Item> {
        if self.k == 0 || self.k > 32 {
            return None;
        }
        let mask = if self.k == 32 {
            u64::MAX
        } else {
            (1_u64 << (self.k * 2)) - 1
        };

        while self.pos < self.encoded.len as usize {
            let pos = self.pos;
            self.pos += 1;
            match base_at(self.encoded, pos) {
                Some(base) => {
                    self.code = ((self.code << 2) | base as u64) & mask;
                    self.valid_run = self.valid_run.saturating_add(1).min(self.k);
                    if self.valid_run >= self.k {
                        return Some(Kmer {
                            code: self.code,
                            pos: (pos + 1 - self.k as usize) as u32,
                            k: self.k,
                        });
                    }
                }
                None => {
                    self.code = 0;
                    self.valid_run = 0;
                }
            }
        }
        None
    }
}

pub struct ClosedSyncmerIter<'a> {
    encoded: &'a EncodedSeq,
    k: u8,
    s: u8,
    pos: usize,
    kmer_code: u64,
    smer_code: u64,
    smer_rc_code: u64,
    valid_run: u8,
    deque_positions: [u32; 33],
    deque_orders: [u64; 33],
    deque_head: usize,
    deque_len: usize,
}

impl Iterator for ClosedSyncmerIter<'_> {
    type Item = Kmer;

    fn next(&mut self) -> Option<Self::Item> {
        if self.k == 0 || self.k > 32 || self.s == 0 || self.s > self.k {
            return None;
        }
        let kmer_mask = low_2bit_mask(self.k);
        let smer_mask = low_2bit_mask(self.s);

        while self.pos < self.encoded.len as usize {
            let pos = self.pos;
            self.pos += 1;
            let Some(base) = base_at(self.encoded, pos) else {
                self.kmer_code = 0;
                self.smer_code = 0;
                self.smer_rc_code = 0;
                self.valid_run = 0;
                self.deque_head = 0;
                self.deque_len = 0;
                continue;
            };

            self.kmer_code = ((self.kmer_code << 2) | base as u64) & kmer_mask;
            self.smer_code = ((self.smer_code << 2) | base as u64) & smer_mask;
            self.smer_rc_code =
                (self.smer_rc_code >> 2) | ((3 - base as u64) << (2 * (self.s as usize - 1)));
            self.valid_run = self.valid_run.saturating_add(1).min(self.k);

            let mut current_smer_order = 0;
            if self.valid_run >= self.s {
                let smer_pos = (pos + 1 - self.s as usize) as u32;
                current_smer_order = syncmer_order(self.smer_code, self.smer_rc_code);
                while self.deque_len > 0 {
                    let back = (self.deque_head + self.deque_len - 1) % 33;
                    if self.deque_orders[back] <= current_smer_order {
                        break;
                    }
                    self.deque_len -= 1;
                }
                let tail = (self.deque_head + self.deque_len) % 33;
                self.deque_positions[tail] = smer_pos;
                self.deque_orders[tail] = current_smer_order;
                self.deque_len += 1;
            }

            if self.valid_run < self.k {
                continue;
            }

            let kmer_pos = (pos + 1 - self.k as usize) as u32;
            while self.deque_len > 0 && self.deque_positions[self.deque_head] < kmer_pos {
                self.deque_head = (self.deque_head + 1) % 33;
                self.deque_len -= 1;
            }
            debug_assert!(self.deque_len > 0);
            let min_order = self.deque_orders[self.deque_head];
            let left_is_min = self.deque_positions[self.deque_head] == kmer_pos;
            let right_is_min = current_smer_order == min_order;
            if left_is_min || right_is_min {
                return Some(Kmer {
                    code: self.kmer_code,
                    pos: kmer_pos,
                    k: self.k,
                });
            }
        }
        None
    }
}

#[inline]
fn low_2bit_mask(len: u8) -> u64 {
    if len == 32 {
        u64::MAX
    } else {
        (1_u64 << (len * 2)) - 1
    }
}

#[inline]
fn syncmer_order(smer_code: u64, smer_rc_code: u64) -> u64 {
    let canonical = smer_code.min(smer_rc_code);
    let mut value = canonical.wrapping_add(0x9e37_79b9_7f4a_7c15);
    value = (value ^ (value >> 30)).wrapping_mul(0xbf58_476d_1ce4_e5b9);
    value = (value ^ (value >> 27)).wrapping_mul(0x94d0_49bb_1331_11eb);
    value ^ (value >> 31)
}

pub fn kmer_code_from_ascii(seq: &[u8]) -> Option<u64> {
    if seq.len() > 32 {
        return None;
    }
    let mut code = 0_u64;
    for &b in seq {
        code = (code << 2) | Base::from_ascii(b)? as u64;
    }
    Some(code)
}

pub fn reverse_complement_kmer_code(mut code: u64, k: u8) -> u64 {
    let mut out = 0_u64;
    for _ in 0..k {
        let base = code & 0b11;
        out = (out << 2) | (base ^ 0b11);
        code >>= 2;
    }
    out
}

pub fn decode_kmer(mut code: u64, k: u8) -> Vec<u8> {
    let mut out = vec![b'A'; k as usize];
    for idx in (0..k as usize).rev() {
        out[idx] = match code & 0b11 {
            0 => b'A',
            1 => b'C',
            2 => b'G',
            3 => b'T',
            _ => unreachable!(),
        };
        code >>= 2;
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn fused_encoding_matches_independent_encoding_and_complexity() {
        let mut encoded = EncodedSeq::default();
        let mut state = 9173_u64;
        let alphabet = b"ACGTUacgtuNN-?";
        for iteration in 0..500 {
            let len = (iteration * 137) % 701;
            let seq: Vec<_> = (0..len)
                .map(|_| {
                    state ^= state << 13;
                    state ^= state >> 7;
                    state ^= state << 17;
                    if iteration % 11 == 0 {
                        state as u8
                    } else {
                        alphabet[state as usize % alphabet.len()]
                    }
                })
                .collect();
            for k in [1, 7, 31, 32] {
                let low_complexity = encode_and_check_complexity(&seq, k, &mut encoded);
                assert_eq!(encoded, encode_acgt(&seq));
                assert_eq!(
                    low_complexity,
                    crate::sketch::is_low_complexity_read(&seq, k)
                );
            }
        }
    }

    use proptest::prelude::*;

    #[test]
    fn encode_decode_roundtrip() {
        let seq = b"ACGTNacgtn";
        let encoded = encode_acgt(seq);
        assert_eq!(decode_to_vec(&encoded), b"ACGTNACGTN");
        assert_eq!(base_at(&encoded, 4), None);
    }

    #[test]
    fn reverse_complement_involution() {
        let seq = encode_acgt(b"ACGTNNGATTACA");
        let twice = reverse_complement(&reverse_complement(&seq));
        assert_eq!(decode_to_vec(&twice), decode_to_vec(&seq));
    }

    #[test]
    fn reverse_complement_kmer_code_matches_sequence() {
        let acg = kmer_code_from_ascii(b"ACG").unwrap();
        assert_eq!(decode_kmer(reverse_complement_kmer_code(acg, 3), 3), b"CGT");
    }

    #[test]
    fn kmer_iterator_matches_string_reference() {
        let encoded = encode_acgt(b"ACGTNACGT");
        let kmers: Vec<_> = iter_kmers_2bit(&encoded, 3)
            .map(|k| (k.pos, decode_kmer(k.code, k.k)))
            .collect();
        assert_eq!(
            kmers,
            vec![
                (0, b"ACG".to_vec()),
                (1, b"CGT".to_vec()),
                (5, b"ACG".to_vec()),
                (6, b"CGT".to_vec()),
            ]
        );
    }

    #[test]
    fn closed_syncmer_iterator_matches_naive_selection() {
        let encoded = encode_acgt(b"ACGTTGCATGTCGCATGATGCATGAGAGCTNACGTACGTACGT");
        for (k, s) in [(7, 3), (11, 7), (16, 9), (31, 21), (32, 32)] {
            let expected = iter_kmers_2bit(&encoded, k)
                .filter(|kmer| naive_closed_syncmer(kmer.code, k, s))
                .collect::<Vec<_>>();
            let actual = iter_closed_syncmers_2bit(&encoded, k, s).collect::<Vec<_>>();
            assert_eq!(actual, expected, "k={k}, s={s}");
        }
    }

    fn naive_closed_syncmer(kmer_code: u64, k: u8, s: u8) -> bool {
        let smer_mask = low_2bit_mask(s);
        let num_smers = k - s + 1;
        let mut min_smer = u64::MAX;
        let mut left = 0;
        let mut right = 0;
        for offset in 0..num_smers {
            let shift = 2 * (num_smers - 1 - offset);
            let smer = (kmer_code >> shift) & smer_mask;
            let order = syncmer_order(smer, reverse_complement_kmer_code(smer, s));
            min_smer = min_smer.min(order);
            if offset == 0 {
                left = order;
            }
            if offset + 1 == num_smers {
                right = order;
            }
        }
        left == min_smer || right == min_smer
    }

    #[test]
    fn hamming_counts_ns_and_length_difference() {
        let a = encode_acgt(b"ACGTN");
        let b = encode_acgt(b"ACGTAAC");
        assert_eq!(hamming_2bit(&a, &b), 3);
    }

    proptest! {
        #[test]
        fn prop_encode_decode_acgt(seq in proptest::collection::vec("[ACGT]".prop_map(|s| s.into_bytes()[0]), 0..200)) {
            let encoded = encode_acgt(&seq);
            prop_assert_eq!(decode_to_vec(&encoded), seq);
        }

        #[test]
        fn prop_reverse_complement_twice(seq in proptest::collection::vec("[ACGTN]".prop_map(|s| s.into_bytes()[0]), 0..200)) {
            let encoded = encode_acgt(&seq);
            let twice = reverse_complement(&reverse_complement(&encoded));
            prop_assert_eq!(decode_to_vec(&twice), decode_to_vec(&encoded));
        }
    }
}

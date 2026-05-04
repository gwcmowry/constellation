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

#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
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
    let len = seq.len() as u32;
    let word_count = seq.len().div_ceil(32);
    let mask_word_count = seq.len().div_ceil(64);
    let mut words = vec![0_u64; word_count];
    let mut n_mask_words = vec![0_u64; mask_word_count];
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

    EncodedSeq {
        words,
        len,
        n_mask_words,
    }
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

pub fn normalize_base(base: u8) -> Option<u8> {
    match base {
        b'A' | b'a' => Some(b'A'),
        b'C' | b'c' => Some(b'C'),
        b'G' | b'g' => Some(b'G'),
        b'T' | b't' | b'U' | b'u' => Some(b'T'),
        _ => None,
    }
}

pub fn normalize_seq(seq: &[u8]) -> Vec<u8> {
    seq.iter()
        .filter_map(|&base| normalize_base(base))
        .collect()
}

pub fn encode_kmer(kmer: &[u8]) -> Option<u64> {
    if kmer.len() > 31 {
        return None;
    }
    let mut value = 0u64;
    for &base in kmer {
        let bits = match base {
            b'A' | b'a' => 0,
            b'C' | b'c' => 1,
            b'G' | b'g' => 2,
            b'T' | b't' | b'U' | b'u' => 3,
            _ => return None,
        };
        value = (value << 2) | bits;
    }
    Some(value)
}

pub fn revcomp(seq: &[u8]) -> Vec<u8> {
    seq.iter()
        .rev()
        .map(|base| match base {
            b'A' | b'a' => b'T',
            b'C' | b'c' => b'G',
            b'G' | b'g' => b'C',
            b'T' | b't' | b'U' | b'u' => b'A',
            _ => b'N',
        })
        .collect()
}

pub fn kmer_keys(seq: &[u8], k: usize) -> Vec<u64> {
    if k == 0 || seq.len() < k {
        return Vec::new();
    }
    seq.windows(k).filter_map(encode_kmer).collect()
}

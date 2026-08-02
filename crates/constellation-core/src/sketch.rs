use crate::dna::{encode_acgt, iter_kmers_2bit, EncodedSeq};
use crate::ReadId;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct SketchRecord {
    pub read_id: ReadId,
    pub key_primary: u64,
    pub key_secondary: u64,
    pub rare_anchor: u64,
    pub gc_bin: u8,
    pub len: u16,
    pub flags: u16,
}

pub fn sketch_read(read_id: ReadId, seq: &[u8], k: u8) -> SketchRecord {
    let encoded = encode_acgt(seq);
    sketch_encoded_read(read_id, seq, &encoded, k)
}

pub fn sketch_encoded_read(
    read_id: ReadId,
    seq: &[u8],
    encoded: &EncodedSeq,
    k: u8,
) -> SketchRecord {
    let mut min_hash = u64::MAX;
    let mut sim = 0_u64;
    let mut count = 0_u64;
    for kmer in iter_kmers_2bit(encoded, k) {
        let h = stable_hash(kmer.code);
        min_hash = min_hash.min(h);
        sim ^= h.rotate_left((count % 63) as u32);
        count += 1;
    }

    let gc = seq
        .iter()
        .filter(|&&b| matches!(b.to_ascii_uppercase(), b'G' | b'C'))
        .count();
    let gc_bin = if seq.is_empty() {
        0
    } else {
        ((gc * 10) / seq.len()).min(9) as u8
    };
    let low_complexity = count <= 1 || min_hash == u64::MAX || is_low_complexity(seq);

    SketchRecord {
        read_id,
        key_primary: min_hash,
        key_secondary: sim,
        rare_anchor: min_hash,
        gc_bin,
        len: seq.len().min(u16::MAX as usize) as u16,
        flags: u16::from(low_complexity),
    }
}

/// Returns the low-complexity decision used by [`sketch_encoded_read`] without
/// calculating sketch hashes that the read-at-a-time gene-EC mapper does not
/// consume.
pub fn is_low_complexity_read(seq: &[u8], k: u8) -> bool {
    let mut counts = [0_usize; 4];
    let mut valid_bases = 0_usize;
    let mut valid_run = 0_usize;
    let mut kmer_count = 0_usize;

    for &base in seq {
        let upper = base.to_ascii_uppercase();
        match upper {
            b'A' => counts[0] += 1,
            b'C' => counts[1] += 1,
            b'G' => counts[2] += 1,
            b'T' => counts[3] += 1,
            _ => {}
        }
        if matches!(upper, b'A' | b'C' | b'G' | b'T') {
            valid_bases += 1;
        }

        if matches!(upper, b'A' | b'C' | b'G' | b'T' | b'U') {
            valid_run += 1;
            if k > 0 && valid_run >= k as usize {
                kmer_count = (kmer_count + 1).min(2);
            }
        } else {
            valid_run = 0;
        }
    }

    if kmer_count <= 1 || valid_bases == 0 {
        return true;
    }
    let max_base = counts.into_iter().max().unwrap_or(0);
    let distinct = counts.into_iter().filter(|&count| count > 0).count();
    max_base * 5 >= valid_bases * 4 || distinct <= 1
}

pub fn sort_sketch_records(records: &mut [SketchRecord]) {
    records.sort_by_key(|record| {
        (
            record.key_primary,
            record.key_secondary,
            record.rare_anchor,
            record.read_id,
        )
    });
}

fn stable_hash(mut x: u64) -> u64 {
    x = x.wrapping_add(0x9E3779B97F4A7C15);
    x = (x ^ (x >> 30)).wrapping_mul(0xBF58476D1CE4E5B9);
    x = (x ^ (x >> 27)).wrapping_mul(0x94D049BB133111EB);
    x ^ (x >> 31)
}

fn is_low_complexity(seq: &[u8]) -> bool {
    if seq.is_empty() {
        return true;
    }
    let mut counts = [0_usize; 4];
    let mut valid = 0_usize;
    for &base in seq {
        match base.to_ascii_uppercase() {
            b'A' => counts[0] += 1,
            b'C' => counts[1] += 1,
            b'G' => counts[2] += 1,
            b'T' => counts[3] += 1,
            _ => continue,
        }
        valid += 1;
    }
    if valid == 0 {
        return true;
    }
    let max_base = counts.into_iter().max().unwrap_or(0);
    let distinct = counts.into_iter().filter(|&count| count > 0).count();
    max_base * 5 >= valid * 4 || distinct <= 1
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn deterministic_keys_for_fixed_read() {
        assert_eq!(
            sketch_read(7, b"ACGTACGT", 3),
            sketch_read(7, b"ACGTACGT", 3)
        );
    }

    #[test]
    fn encoded_sketch_matches_ascii_sketch() {
        let seq = b"ACGTACGTACGT";
        let encoded = encode_acgt(seq);
        assert_eq!(
            sketch_encoded_read(9, seq, &encoded, 4),
            sketch_read(9, seq, 4)
        );
    }

    #[test]
    fn flags_low_complexity_homopolymer() {
        assert_ne!(sketch_read(0, b"AAAAAAAAAAAAAAAAAAAAAAAA", 3).flags & 1, 0);
    }

    #[test]
    fn lightweight_low_complexity_check_matches_full_sketch() {
        for seq in [
            b"".as_slice(),
            b"AAAAAAAAAAAAAAAAAAAAAAAA",
            b"ACGTACGTACGTACGTACGTACGT",
            b"ACGTNNNNACGTACGTACGTACGT",
            b"UUUUACGTACGTACGTACGTACGT",
            b"acgtacgtacgtacgtacgtacgt",
        ] {
            for k in [1, 3, 15, 21, 31] {
                assert_eq!(
                    is_low_complexity_read(seq, k),
                    sketch_read(0, seq, k).flags & 1 != 0,
                    "sequence={:?}, k={k}",
                    String::from_utf8_lossy(seq)
                );
            }
        }
    }
}

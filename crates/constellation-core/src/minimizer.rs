use crate::dna::{encode_acgt, iter_kmers_2bit};

pub fn minimizer_code(seq: &[u8], k: u8) -> Option<u64> {
    let encoded = encode_acgt(seq);
    iter_kmers_2bit(&encoded, k).map(|kmer| kmer.code).min()
}

use constellation_core::dna::kmer_code_from_ascii;
use constellation_core::index::{KmerEntry, Posting, TranscriptIndex, TranscriptMeta};
use criterion::{criterion_group, criterion_main, Criterion};
use std::hint::black_box;

fn bench_lookup(c: &mut Criterion) {
    let code = kmer_code_from_ascii(b"ACGTACGTACGTACG").unwrap();
    let index = TranscriptIndex {
        format_version: 1,
        k: 15,
        max_kmer_frequency: 256,
        genes: vec!["GENE".to_owned()],
        transcripts: vec![TranscriptMeta {
            transcript_id: 0,
            gene_id: 0,
            name: "tx".to_owned(),
            len: 64,
        }],
        transcript_sequences: vec!["ACGT".repeat(16)],
        kmers: vec![KmerEntry {
            kmer_code: code,
            postings_start: 0,
            postings_len: 16,
            freq_class: 0,
        }],
        postings: (0..16)
            .map(|pos| Posting {
                transcript_id: 0,
                gene_id: 0,
                pos,
                strand: 0,
            })
            .collect(),
    };
    c.bench_function("lookup_single_kmer", |b| {
        b.iter(|| black_box(index.lookup(black_box(code))).len());
    });
}

criterion_group!(benches, bench_lookup);
criterion_main!(benches);

use constellation_core::candidate::{make_candidate_locus_buckets, CandidateHit};
use criterion::{criterion_group, criterion_main, Criterion};
use std::hint::black_box;

fn bench_bucket(c: &mut Criterion) {
    let hits: Vec<_> = (0..4096)
        .rev()
        .map(|i| CandidateHit {
            read_id: i as u64,
            transcript_id: (i % 32) as u32,
            gene_id: (i % 16) as u32,
            pos: (i * 7 % 2048) as u32,
            strand: 0,
            seed_count: 1,
            seed_score: 1,
        })
        .collect();
    c.bench_function("candidate_locus_bucket_4096", |b| {
        b.iter(|| make_candidate_locus_buckets(black_box(hits.clone()), black_box(64)));
    });
}

criterion_group!(benches, bench_bucket);
criterion_main!(benches);

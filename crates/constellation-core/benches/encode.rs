use constellation_core::dna::encode_acgt;
use criterion::{criterion_group, criterion_main, Criterion};
use std::hint::black_box;

fn bench_encode(c: &mut Criterion) {
    let seq = b"ACGT".repeat(4096);
    c.bench_function("encode_16kb", |b| {
        b.iter(|| encode_acgt(black_box(&seq)));
    });
}

criterion_group!(benches, bench_encode);
criterion_main!(benches);

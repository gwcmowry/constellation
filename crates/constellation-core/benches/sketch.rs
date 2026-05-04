use constellation_core::sketch::sketch_read;
use criterion::{criterion_group, criterion_main, Criterion};
use std::hint::black_box;

fn bench_sketch(c: &mut Criterion) {
    let seq = b"ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT".to_vec();
    c.bench_function("sketch_read_k15", |b| {
        b.iter(|| sketch_read(black_box(0), black_box(&seq), black_box(15)));
    });
}

criterion_group!(benches, bench_sketch);
criterion_main!(benches);

use constellation_core::candidate::{generate_candidate_hits, make_candidate_locus_buckets};
use constellation_core::index_build::build_transcript_index;
use constellation_core::score::CandidateScorer;
use constellation_core::score_scalar::ScalarScorer;
use criterion::{criterion_group, criterion_main, Criterion};
use std::hint::black_box;
use std::io::Write;

fn bench_end_to_end(c: &mut Criterion) {
    let mut fasta = tempfile::NamedTempFile::new().unwrap();
    writeln!(
        fasta,
        ">tx1|gene=GENE_A\n{}\n>tx2|gene=GENE_B\n{}",
        "ACGT".repeat(256),
        "TGCA".repeat(256)
    )
    .unwrap();
    let index = build_transcript_index(fasta.path(), 15, 256).unwrap();
    let reads: Vec<_> = (0..512)
        .map(|read_id| (read_id, b"ACGTACGTACGTACGTACGTACGTACGTACGT".to_vec()))
        .collect();
    let quals: Vec<_> = reads
        .iter()
        .map(|(read_id, seq)| (*read_id, vec![b'I'; seq.len()]))
        .collect();

    c.bench_function("end_to_end_synthetic_512", |b| {
        b.iter(|| {
            let mut hits = Vec::new();
            for (read_id, seq) in &reads {
                hits.extend(generate_candidate_hits(&index, *read_id, seq, 8, 256));
            }
            let buckets = make_candidate_locus_buckets(hits, 64);
            let mut scored = Vec::new();
            for bucket in &buckets {
                ScalarScorer::default().score_bucket(
                    black_box(&reads),
                    black_box(&quals),
                    black_box(&index),
                    black_box(bucket),
                    black_box(&mut scored),
                    None,
                );
            }
            scored
        });
    });
}

criterion_group!(benches, bench_end_to_end);
criterion_main!(benches);

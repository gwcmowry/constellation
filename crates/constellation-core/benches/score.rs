use constellation_core::candidate::{CandidateBucketKey, CandidateHit, CandidateLocusBucket};
use constellation_core::index::{TranscriptIndex, TranscriptMeta};
use constellation_core::score::CandidateScorer;
use constellation_core::score_scalar::ScalarScorer;
use constellation_core::score_simd::PulpScorer;
use criterion::{criterion_group, criterion_main, Criterion};
use std::hint::black_box;

fn bench_score(c: &mut Criterion) {
    let index = TranscriptIndex {
        format_version: 1,
        k: 15,
        max_kmer_frequency: 256,
        genes: vec!["GENE".to_owned()],
        transcripts: vec![TranscriptMeta {
            transcript_id: 0,
            gene_id: 0,
            name: "tx".to_owned(),
            len: 4096,
        }],
        transcript_sequences: vec!["ACGT".repeat(1024)],
        kmers: vec![],
        postings: vec![],
    };
    let reads: Vec<_> = (0..1024)
        .map(|read_id| (read_id, b"ACGTACGTACGTACGTACGTACGTACGTACGT".to_vec()))
        .collect();
    let quals: Vec<_> = reads
        .iter()
        .map(|(read_id, seq)| (*read_id, vec![b'I'; seq.len()]))
        .collect();
    let bucket = CandidateLocusBucket {
        key: CandidateBucketKey {
            transcript_id: 0,
            pos_bin: 0,
            strand: 0,
        },
        hits: reads
            .iter()
            .map(|(read_id, _)| CandidateHit {
                read_id: *read_id,
                transcript_id: 0,
                gene_id: 0,
                pos: 0,
                strand: 0,
                seed_count: 1,
                seed_score: 1,
            })
            .collect(),
    };

    c.bench_function("score_scalar_1024", |b| {
        b.iter(|| {
            let mut out = Vec::new();
            ScalarScorer::default().score_bucket(
                black_box(&reads),
                black_box(&quals),
                black_box(&index),
                black_box(&bucket),
                black_box(&mut out),
                None,
            );
            out
        });
    });

    c.bench_function("score_pulp_trait_1024", |b| {
        b.iter(|| {
            let mut out = Vec::new();
            PulpScorer::default().score_bucket(
                black_box(&reads),
                black_box(&quals),
                black_box(&index),
                black_box(&bucket),
                black_box(&mut out),
                None,
            );
            out
        });
    });
}

criterion_group!(benches, bench_score);
criterion_main!(benches);

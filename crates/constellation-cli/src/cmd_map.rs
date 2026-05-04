use crate::{MapArgs, MapMode, ScoreMode};
use anyhow::{anyhow, Result};
use constellation_core::assign::{assign_read, flagged_assignment, Assignment, AssignmentType};
use constellation_core::candidate::{
    generate_candidate_hits_with_quality_stats, make_candidate_locus_buckets,
    make_single_hit_buckets, CandidateEarlyStopConfig, CandidateGenerationStats,
    CandidateLocusBucket,
};
use constellation_core::chemistry::{parse_tenx_3p_v3_r1, BarcodeUmi, Chemistry};
use constellation_core::fastq::read_fastq;
use constellation_core::index::{IndexAccess, LoadedIndex};
use constellation_core::metrics::MapMetrics;
use constellation_core::score::{CandidateScorer, ScoredCandidate};
use constellation_core::score_scalar::ScalarScorer;
use constellation_core::score_simd::PulpScorer;
use constellation_core::sketch::{sketch_read, sort_sketch_records};
use rayon::prelude::*;
use std::collections::BTreeMap;
use std::fs;
use std::time::{Duration, Instant};

pub fn run_map(args: MapArgs) -> Result<()> {
    let _chemistry = Chemistry::parse(&args.chemistry)?;
    let mut function_timings = FunctionTimings::default();
    let total_started = Instant::now();
    let index_started = Instant::now();
    let index = LoadedIndex::load(&args.index)?;
    let index_load_seconds = index_started.elapsed().as_secs_f64();
    function_timings.add_elapsed("LoadedIndex::load", index_started);
    let r1_started = Instant::now();
    let r1 = read_fastq(&args.r1)?;
    function_timings.add_elapsed("read_fastq_r1", r1_started);
    let r2_started = Instant::now();
    let r2 = read_fastq(&args.r2)?;
    function_timings.add_elapsed("read_fastq_r2", r2_started);
    let fastq_load_seconds = r1_started.elapsed().as_secs_f64();
    if r1.len() != r2.len() {
        return Err(anyhow!(
            "R1/R2 record count mismatch: {} != {}",
            r1.len(),
            r2.len()
        ));
    }

    let started = Instant::now();
    let mut reads = Vec::with_capacity(r2.len());
    let mut quals = Vec::with_capacity(r2.len());
    let mut barcode_umis = Vec::with_capacity(r1.len());
    let mut low_quality = Vec::with_capacity(r2.len());
    let mut invalid_barcode_umi = Vec::with_capacity(r2.len());
    let mut sketches = Vec::with_capacity(r2.len());
    let mut sketch_flags_by_read = Vec::with_capacity(r2.len());
    let mut total_bases = 0_u64;

    let prep_started = Instant::now();
    let prepared_reads: Vec<_> = r1
        .par_iter()
        .zip(r2.par_iter())
        .enumerate()
        .map(|(idx, (r1_record, r2_record))| {
            let read_id = idx as u64;
            let (barcode_umi, invalid_bc_umi) = parse_tenx_3p_v3_r1(&r1_record.seq)
                .map(|parsed| (parsed, false))
                .unwrap_or_else(|_| (lossy_barcode_umi(&r1_record.seq), true));
            let low_quality = mean_phred_quality(&r2_record.qual) < args.min_mean_quality;
            let sketch = sketch_read(read_id, &r2_record.seq, index.k());
            PreparedRead {
                read_id,
                seq: r2_record.seq.clone(),
                qual: r2_record.qual.clone(),
                cell_barcode: barcode_umi.cell_barcode_seq,
                umi: barcode_umi.umi_seq,
                invalid_barcode_umi: invalid_bc_umi,
                low_quality,
                sketch,
            }
        })
        .collect();
    function_timings.add_elapsed("prepare_reads_parallel", prep_started);

    for prepared in prepared_reads {
        total_bases += prepared.seq.len() as u64;
        barcode_umis.push((prepared.cell_barcode, prepared.umi));
        invalid_barcode_umi.push(prepared.invalid_barcode_umi);
        low_quality.push(prepared.low_quality);
        reads.push((prepared.read_id, prepared.seq));
        quals.push((prepared.read_id, prepared.qual));
        sketch_flags_by_read.push(prepared.sketch.flags);
        sketches.push(prepared.sketch);
    }
    let preprocess_seconds = started.elapsed().as_secs_f64();

    let scorer = make_scorer(args.score_mode);
    let early_stop = early_stop_config(&args);
    let mut stats = CandidateGenerationStats::default();
    let mut stage_times = MappingStageTimes::default();
    let mut bucket_sizes = Vec::new();
    let mapping_started = Instant::now();
    let assignments = match args.mode {
        MapMode::ReadAtATime => map_read_at_a_time(
            &args,
            &reads,
            &quals,
            &barcode_umis,
            &low_quality,
            &invalid_barcode_umi,
            &sketch_flags_by_read,
            &index,
            scorer.as_ref(),
            early_stop,
            &mut stats,
            &mut stage_times,
            &mut function_timings,
            &mut bucket_sizes,
        ),
        MapMode::SketchBucket => {
            let sort_started = Instant::now();
            sort_sketch_records(&mut sketches);
            function_timings.add_elapsed("sort_sketch_records", sort_started);
            map_bucketed(
                &args,
                &reads,
                &quals,
                &barcode_umis,
                &low_quality,
                &invalid_barcode_umi,
                &sketch_flags_by_read,
                &sketches,
                &index,
                scorer.as_ref(),
                early_stop,
                false,
                &mut stats,
                &mut stage_times,
                &mut function_timings,
                &mut bucket_sizes,
            )
        }
        MapMode::CandidateLocus => {
            let sort_started = Instant::now();
            sort_sketch_records(&mut sketches);
            function_timings.add_elapsed("sort_sketch_records", sort_started);
            map_bucketed(
                &args,
                &reads,
                &quals,
                &barcode_umis,
                &low_quality,
                &invalid_barcode_umi,
                &sketch_flags_by_read,
                &sketches,
                &index,
                scorer.as_ref(),
                early_stop,
                true,
                &mut stats,
                &mut stage_times,
                &mut function_timings,
                &mut bucket_sizes,
            )
        }
    };
    let mapping_seconds = mapping_started.elapsed().as_secs_f64();

    let write_started = Instant::now();
    let mut assignment_tsv = String::from(Assignment::tsv_header());
    for assignment in &assignments {
        assignment_tsv.push_str(&assignment.to_tsv_row());
    }
    fs::write(&args.out, assignment_tsv)?;
    let write_seconds = write_started.elapsed().as_secs_f64();

    if let Some(metrics_path) = args.emit_metrics {
        let mut metrics = MapMetrics {
            mode: args.mode.as_str().to_owned(),
            score_mode: args.score_mode.as_str().to_owned(),
            num_reads: reads.len() as u64,
            index_load_seconds,
            fastq_load_seconds,
            preprocess_seconds,
            mapping_seconds,
            candidate_generation_seconds: stage_times.candidate_generation_seconds,
            bucket_build_seconds: stage_times.bucket_build_seconds,
            scoring_seconds: stage_times.scoring_seconds,
            assignment_seconds: stage_times.assignment_seconds,
            write_seconds,
            seed_lookups_per_read: stats.seed_lookups as f64 / reads.len().max(1) as f64,
            seed_candidates_considered_per_read: stats.seed_candidates_considered as f64
                / reads.len().max(1) as f64,
            candidate_hits_per_read: per_read(stats.candidate_hits as usize, reads.len()),
            candidate_buckets_per_read: per_read(bucket_sizes.len(), reads.len()),
            postings_skipped_due_to_frequency: stats.postings_skipped_due_to_frequency,
            seeds_skipped_due_to_quality: stats.seeds_skipped_due_to_quality,
            early_stopped_reads: stats.early_stopped_reads,
            seed_lookups_saved_by_early_stop: stats.seed_lookups_saved_by_early_stop,
            scored_candidates_per_read: assignments
                .iter()
                .map(|assignment| assignment.candidate_count as u64)
                .sum::<u64>() as f64
                / reads.len().max(1) as f64,
            mean_candidate_bucket_size: mean(&bucket_sizes),
            median_candidate_bucket_size: median(bucket_sizes.clone()),
            max_candidate_bucket_size: bucket_sizes.iter().copied().max().unwrap_or(0) as u64,
            unique_gene_rate: rate(&assignments, AssignmentType::UniqueGene),
            ambiguous_gene_rate: rate(&assignments, AssignmentType::AmbiguousGene)
                + rate(&assignments, AssignmentType::AmbiguousTranscriptSameGene),
            low_complexity_rate: rate(&assignments, AssignmentType::LowComplexity),
            low_quality_rate: rate(&assignments, AssignmentType::LowQuality),
            unmapped_rate: rate(&assignments, AssignmentType::Unmapped),
            function_seconds: function_timings.into_seconds(),
            ..MapMetrics::default()
        };
        metrics.finish_rates(total_started.elapsed(), total_bases);
        fs::write(metrics_path, serde_json::to_vec_pretty(&metrics)?)?;
    }

    Ok(())
}

#[derive(Debug, Default)]
struct MappingStageTimes {
    candidate_generation_seconds: f64,
    bucket_build_seconds: f64,
    scoring_seconds: f64,
    assignment_seconds: f64,
}

struct PreparedRead {
    read_id: u64,
    seq: Vec<u8>,
    qual: Vec<u8>,
    cell_barcode: String,
    umi: String,
    invalid_barcode_umi: bool,
    low_quality: bool,
    sketch: constellation_core::sketch::SketchRecord,
}

#[derive(Debug, Default)]
struct FunctionTimings {
    totals: BTreeMap<&'static str, Duration>,
}

impl FunctionTimings {
    fn add_elapsed(&mut self, name: &'static str, started: Instant) {
        *self.totals.entry(name).or_default() += started.elapsed();
    }

    fn into_seconds(self) -> BTreeMap<String, f64> {
        self.totals
            .into_iter()
            .map(|(name, duration)| (name.to_owned(), duration.as_secs_f64()))
            .collect()
    }
}

fn make_scorer(score_mode: ScoreMode) -> Box<dyn CandidateScorer> {
    match score_mode {
        ScoreMode::Scalar => Box::new(ScalarScorer::default()),
        ScoreMode::Pulp => Box::new(PulpScorer::default()),
    }
}

fn early_stop_config(args: &MapArgs) -> Option<CandidateEarlyStopConfig> {
    let config = CandidateEarlyStopConfig {
        posterior_threshold: args.early_stop_posterior,
        prior_alpha: args.early_stop_prior_alpha,
        min_seed_lookups: args.early_stop_min_seed_lookups,
        min_top_seed_count: args.early_stop_min_top_seed_count,
    };
    config.is_enabled().then_some(config)
}

fn map_read_at_a_time(
    args: &MapArgs,
    reads: &[(u64, Vec<u8>)],
    quals: &[(u64, Vec<u8>)],
    barcode_umis: &[(String, String)],
    low_quality: &[bool],
    invalid_barcode_umi: &[bool],
    sketch_flags_by_read: &[u16],
    index: &dyn IndexAccess,
    scorer: &dyn CandidateScorer,
    early_stop: Option<CandidateEarlyStopConfig>,
    stats: &mut CandidateGenerationStats,
    stage_times: &mut MappingStageTimes,
    function_timings: &mut FunctionTimings,
    bucket_sizes: &mut Vec<usize>,
) -> Vec<Assignment> {
    let mut assignments = Vec::with_capacity(reads.len());
    for (read_id, seq) in reads {
        let (cell_barcode, umi) = barcode_umis[*read_id as usize].clone();
        if low_quality[*read_id as usize] || invalid_barcode_umi[*read_id as usize] {
            let assignment_started = Instant::now();
            assignments.push(flagged_assignment(
                *read_id,
                cell_barcode,
                umi,
                AssignmentType::LowQuality,
                if invalid_barcode_umi[*read_id as usize] {
                    4
                } else {
                    2
                },
            ));
            stage_times.assignment_seconds += assignment_started.elapsed().as_secs_f64();
            function_timings.add_elapsed("flagged_assignment", assignment_started);
            continue;
        }
        if sketch_flags_by_read[*read_id as usize] & 1 != 0 {
            let assignment_started = Instant::now();
            assignments.push(flagged_assignment(
                *read_id,
                cell_barcode,
                umi,
                AssignmentType::LowComplexity,
                sketch_flags_by_read[*read_id as usize],
            ));
            stage_times.assignment_seconds += assignment_started.elapsed().as_secs_f64();
            function_timings.add_elapsed("flagged_assignment", assignment_started);
            continue;
        }

        let candidate_started = Instant::now();
        let (hits, read_stats) = generate_candidate_hits_with_quality_stats(
            index,
            *read_id,
            seq,
            Some(&quals[*read_id as usize].1),
            args.min_seed_quality,
            args.max_seeds_per_read,
            args.max_postings_per_seed,
            args.search_reverse_complement,
            early_stop,
        );
        stage_times.candidate_generation_seconds += candidate_started.elapsed().as_secs_f64();
        function_timings.add_elapsed(
            "generate_candidate_hits_with_quality_stats",
            candidate_started,
        );
        add_stats(stats, read_stats);
        let bucket_started = Instant::now();
        let buckets = make_single_hit_buckets(hits, args.candidate_bin_size);
        stage_times.bucket_build_seconds += bucket_started.elapsed().as_secs_f64();
        function_timings.add_elapsed("make_single_hit_buckets", bucket_started);
        bucket_sizes.extend(buckets.iter().map(|bucket| bucket.hits.len()));
        let mut scored = Vec::new();
        let scoring_started = Instant::now();
        for bucket in &buckets {
            let score_bucket_started = Instant::now();
            scorer.score_bucket(reads, index, bucket, &mut scored);
            function_timings.add_elapsed("CandidateScorer::score_bucket", score_bucket_started);
        }
        stage_times.scoring_seconds += scoring_started.elapsed().as_secs_f64();
        let sort_started = Instant::now();
        scored.sort();
        function_timings.add_elapsed("ScoredCandidate::sort", sort_started);
        let assignment_started = Instant::now();
        assignments.push(assign_read(*read_id, cell_barcode, umi, &scored));
        stage_times.assignment_seconds += assignment_started.elapsed().as_secs_f64();
        function_timings.add_elapsed("assign_read", assignment_started);
    }
    assignments
}

#[allow(clippy::too_many_arguments)]
fn map_bucketed(
    args: &MapArgs,
    reads: &[(u64, Vec<u8>)],
    quals: &[(u64, Vec<u8>)],
    barcode_umis: &[(String, String)],
    low_quality: &[bool],
    invalid_barcode_umi: &[bool],
    sketch_flags_by_read: &[u16],
    sketches: &[constellation_core::sketch::SketchRecord],
    index: &dyn IndexAccess,
    scorer: &dyn CandidateScorer,
    early_stop: Option<CandidateEarlyStopConfig>,
    regroup_by_locus: bool,
    stats: &mut CandidateGenerationStats,
    stage_times: &mut MappingStageTimes,
    function_timings: &mut FunctionTimings,
    bucket_sizes: &mut Vec<usize>,
) -> Vec<Assignment> {
    let candidate_started = Instant::now();
    let candidate_work: Vec<_> = sketches
        .par_iter()
        .map(|sketch| {
            if low_quality[sketch.read_id as usize] || invalid_barcode_umi[sketch.read_id as usize]
            {
                return CandidateWork::LowQuality(sketch.read_id);
            }
            if sketch.flags & 1 != 0 {
                return CandidateWork::LowComplexity(sketch.read_id);
            }
            let (_, seq) = &reads[sketch.read_id as usize];
            let (hits, read_stats) = generate_candidate_hits_with_quality_stats(
                index,
                sketch.read_id,
                seq,
                Some(&quals[sketch.read_id as usize].1),
                args.min_seed_quality,
                args.max_seeds_per_read,
                args.max_postings_per_seed,
                args.search_reverse_complement,
                early_stop,
            );
            CandidateWork::Hits { hits, read_stats }
        })
        .collect();
    stage_times.candidate_generation_seconds += candidate_started.elapsed().as_secs_f64();
    function_timings.add_elapsed(
        "generate_candidate_hits_with_quality_stats_parallel",
        candidate_started,
    );

    let mut all_hits = Vec::new();
    let mut low_complexity = Vec::new();
    let mut low_quality_reads = Vec::new();
    for work in candidate_work {
        match work {
            CandidateWork::Hits { hits, read_stats } => {
                add_stats(stats, read_stats);
                all_hits.extend(hits);
            }
            CandidateWork::LowQuality(read_id) => low_quality_reads.push(read_id),
            CandidateWork::LowComplexity(read_id) => low_complexity.push(read_id),
        }
    }
    low_complexity.sort_unstable();
    low_quality_reads.sort_unstable();

    let bucket_started = Instant::now();
    let buckets = if regroup_by_locus {
        make_candidate_locus_buckets(all_hits, args.candidate_bin_size)
    } else {
        make_single_hit_buckets(all_hits, args.candidate_bin_size)
    };
    stage_times.bucket_build_seconds += bucket_started.elapsed().as_secs_f64();
    if regroup_by_locus {
        function_timings.add_elapsed("make_candidate_locus_buckets", bucket_started);
    } else {
        function_timings.add_elapsed("make_single_hit_buckets", bucket_started);
    }
    bucket_sizes.extend(buckets.iter().map(|bucket| bucket.hits.len()));
    let scoring_started = Instant::now();
    let by_read = score_buckets(reads, index, &buckets, scorer, function_timings);
    stage_times.scoring_seconds += scoring_started.elapsed().as_secs_f64();

    let mut assignments = Vec::with_capacity(reads.len());
    let assignment_started = Instant::now();
    for (read_id, _) in reads {
        let (cell_barcode, umi) = barcode_umis[*read_id as usize].clone();
        if low_quality_reads.binary_search(read_id).is_ok() {
            assignments.push(flagged_assignment(
                *read_id,
                cell_barcode,
                umi,
                AssignmentType::LowQuality,
                if invalid_barcode_umi[*read_id as usize] {
                    4
                } else {
                    2
                },
            ));
            continue;
        }
        if low_complexity.binary_search(read_id).is_ok() {
            assignments.push(flagged_assignment(
                *read_id,
                cell_barcode,
                umi,
                AssignmentType::LowComplexity,
                sketch_flags_by_read[*read_id as usize],
            ));
            continue;
        }
        let candidates = by_read
            .get(*read_id as usize)
            .map(Vec::as_slice)
            .unwrap_or(&[]);
        assignments.push(assign_read(*read_id, cell_barcode, umi, &candidates));
    }
    stage_times.assignment_seconds += assignment_started.elapsed().as_secs_f64();
    function_timings.add_elapsed("assignments_by_read", assignment_started);

    assignments
}

enum CandidateWork {
    Hits {
        hits: Vec<constellation_core::candidate::CandidateHit>,
        read_stats: CandidateGenerationStats,
    },
    LowQuality(u64),
    LowComplexity(u64),
}

fn score_buckets(
    reads: &[(u64, Vec<u8>)],
    index: &dyn IndexAccess,
    buckets: &[CandidateLocusBucket],
    scorer: &dyn CandidateScorer,
    function_timings: &mut FunctionTimings,
) -> Vec<Vec<ScoredCandidate>> {
    let mut scored = Vec::new();
    for bucket in buckets {
        let score_bucket_started = Instant::now();
        scorer.score_bucket(reads, index, bucket, &mut scored);
        function_timings.add_elapsed("CandidateScorer::score_bucket", score_bucket_started);
    }
    let sort_started = Instant::now();
    scored.sort();
    function_timings.add_elapsed("ScoredCandidate::sort", sort_started);

    let group_started = Instant::now();
    let mut by_read = vec![Vec::new(); reads.len()];
    for candidate in scored {
        if let Some(candidates) = by_read.get_mut(candidate.read_id as usize) {
            candidates.push(candidate);
        }
    }
    function_timings.add_elapsed("group_scored_candidates_by_read", group_started);
    by_read
}

fn add_stats(total: &mut CandidateGenerationStats, next: CandidateGenerationStats) {
    total.seed_lookups += next.seed_lookups;
    total.seed_candidates_considered += next.seed_candidates_considered;
    total.candidate_hits += next.candidate_hits;
    total.postings_skipped_due_to_frequency += next.postings_skipped_due_to_frequency;
    total.seeds_skipped_due_to_quality += next.seeds_skipped_due_to_quality;
    total.early_stopped_reads += next.early_stopped_reads;
    total.seed_lookups_saved_by_early_stop += next.seed_lookups_saved_by_early_stop;
}

fn per_read(count: usize, reads: usize) -> f64 {
    count as f64 / reads.max(1) as f64
}

fn mean_phred_quality(qual: &[u8]) -> f64 {
    if qual.is_empty() {
        return 0.0;
    }
    qual.iter()
        .map(|&byte| byte.saturating_sub(33) as u64)
        .sum::<u64>() as f64
        / qual.len() as f64
}

fn lossy_barcode_umi(seq: &[u8]) -> BarcodeUmi {
    let mut padded = vec![b'N'; 28];
    for (dst, src) in padded.iter_mut().zip(seq.iter().copied()) {
        *dst = src.to_ascii_uppercase();
    }
    let cell_barcode_seq = String::from_utf8_lossy(&padded[..16]).into_owned();
    let umi_seq = String::from_utf8_lossy(&padded[16..28]).into_owned();
    BarcodeUmi {
        cell_barcode: 0,
        umi: 0,
        cell_barcode_seq,
        umi_seq,
    }
}

fn mean(values: &[usize]) -> f64 {
    if values.is_empty() {
        0.0
    } else {
        values.iter().sum::<usize>() as f64 / values.len() as f64
    }
}

fn median(mut values: Vec<usize>) -> f64 {
    if values.is_empty() {
        return 0.0;
    }
    values.sort_unstable();
    let mid = values.len() / 2;
    if values.len() % 2 == 0 {
        (values[mid - 1] + values[mid]) as f64 / 2.0
    } else {
        values[mid] as f64
    }
}

fn rate(assignments: &[constellation_core::assign::Assignment], ty: AssignmentType) -> f64 {
    assignments
        .iter()
        .filter(|assignment| assignment.assignment_type == ty)
        .count() as f64
        / assignments.len().max(1) as f64
}

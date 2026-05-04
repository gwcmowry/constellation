use crate::{
    CandidatePruningArg, CandidateSearchArg, LibraryStrandArg, MapArgs, MapMode, ReadChunkingArg,
    RetrievalModeArg, ScoreMode, SeedPlannerArg,
};
use anyhow::{anyhow, Result};
use constellation_core::assign::{assign_read, flagged_assignment, Assignment, AssignmentType};
use constellation_core::candidate::{
    best_sparse_probe_seed_key, generate_candidate_hits_seed_batched,
    generate_candidate_hits_with_quality_stats, make_candidate_locus_buckets,
    make_single_hit_buckets, CandidateEarlyStopConfig, CandidateGenerationStats,
    CandidateLocusBucket, CandidatePruningMode, CandidateSearchMode, SeedPlanner,
    SparseProbeConfig, WandLiteConfig,
};
use constellation_core::chemistry::{parse_tenx_3p_v3_r1, BarcodeUmi, Chemistry};
use constellation_core::fastq::read_fastq;
use constellation_core::index::{IndexAccess, LoadedIndex};
use constellation_core::metrics::MapMetrics;
use constellation_core::score::{
    tso_prefix_trim_len, CandidateScorer, LibraryStrand, ScoreConfig, ScoreFailureStats,
    ScoredCandidate,
};
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
            let tso_trim = if args.trim_tso {
                tso_prefix_trim_len(
                    &r2_record.seq,
                    args.min_tso_match_len,
                    args.max_tso_mismatches,
                )
            } else {
                0
            };
            let seq = r2_record.seq[tso_trim..].to_vec();
            let qual = r2_record.qual[tso_trim..].to_vec();
            let low_quality = mean_phred_quality(&qual) < args.min_mean_quality;
            let sketch = sketch_read(read_id, &seq, index.k());
            PreparedRead {
                read_id,
                seq,
                qual,
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

    let score_config = score_config(&args);
    let scorer = make_scorer(args.score_mode, score_config);
    let early_stop = early_stop_config(&args);
    let mut stats = CandidateGenerationStats::default();
    let mut score_stats = ScoreFailureStats::default();
    let mut stage_times = MappingStageTimes::default();
    let mut bucket_sizes = Vec::new();
    let mut read_candidate_stats = vec![CandidateGenerationStats::default(); reads.len()];
    let mut read_score_stats = vec![ScoreFailureStats::default(); reads.len()];
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
            &mut score_stats,
            &mut stage_times,
            &mut function_timings,
            &mut bucket_sizes,
            &mut read_candidate_stats,
            &mut read_score_stats,
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
                &mut score_stats,
                &mut stage_times,
                &mut function_timings,
                &mut bucket_sizes,
                &mut read_candidate_stats,
                &mut read_score_stats,
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
                &mut score_stats,
                &mut stage_times,
                &mut function_timings,
                &mut bucket_sizes,
                &mut read_candidate_stats,
                &mut read_score_stats,
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

    if let Some(path) = args.emit_unmapped_diagnostics.as_ref() {
        fs::write(
            path,
            unmapped_diagnostics_tsv(
                &assignments,
                &reads,
                &read_candidate_stats,
                &read_score_stats,
            ),
        )?;
    }

    if let Some(metrics_path) = args.emit_metrics {
        let unique_gene_rate = rate(&assignments, AssignmentType::UniqueGene);
        let same_gene_multitranscript_rate =
            rate(&assignments, AssignmentType::AmbiguousTranscriptSameGene);
        let multi_gene_ambiguous_rate = rate(&assignments, AssignmentType::AmbiguousGene);
        let antisense_gene_rate = rate(&assignments, AssignmentType::AntisenseGene);
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
            reads_with_no_valid_kmers: stats.reads_with_no_valid_kmers,
            reads_with_no_selected_seeds: stats.reads_with_no_selected_seeds,
            reads_with_no_seed_postings: stats.reads_with_no_seed_postings,
            reads_all_selected_seeds_over_frequency_cap: stats
                .reads_all_selected_seeds_over_frequency_cap,
            reads_with_candidate_hits: stats.reads_with_candidate_hits,
            selected_seeds_absent_from_index: stats.selected_seeds_absent_from_index,
            selected_seeds_over_frequency_cap: stats.selected_seeds_over_frequency_cap,
            selected_seeds_with_postings: stats.selected_seeds_with_postings,
            query_seed_occurrences_per_read: stats.query_seed_occurrences as f64
                / reads.len().max(1) as f64,
            distinct_seed_lists_loaded: stats.distinct_seed_lists_loaded,
            posting_list_reuse_factor: stats.query_seed_occurrences as f64
                / stats.distinct_seed_lists_loaded.max(1) as f64,
            candidate_votes_per_read: stats.candidate_votes as f64 / reads.len().max(1) as f64,
            seed_groups_skipped_due_to_frequency: stats.seed_groups_skipped_due_to_frequency,
            sparse_probe_attempted_reads: stats.sparse_probe_attempted_reads,
            sparse_probe_accepted_reads: stats.sparse_probe_accepted_reads,
            sparse_probe_fallback_reads: stats.sparse_probe_fallback_reads,
            candidate_hits_pruned: stats.candidate_hits_pruned,
            score_candidates_seen: score_stats.candidates_seen,
            score_candidates_out_of_bounds: score_stats.candidates_out_of_bounds,
            score_candidates_failed_mismatch: score_stats.candidates_failed_mismatch,
            score_candidates_failed_min_scored_len: score_stats.candidates_failed_min_scored_len,
            score_candidates_failed_trimmed_out_of_bounds: score_stats
                .candidates_failed_trimmed_out_of_bounds,
            score_candidates_full_length_mismatch_only: score_stats
                .candidates_full_length_mismatch_only,
            score_candidates_passed_full_length: score_stats.candidates_passed_full_length,
            score_candidates_passed_trimmed_or_softclipped: score_stats
                .candidates_passed_trimmed_or_softclipped,
            unique_gene_rate,
            ambiguous_gene_rate: multi_gene_ambiguous_rate,
            same_gene_multitranscript_rate,
            multi_gene_ambiguous_rate,
            antisense_gene_rate,
            gene_countable_rate: unique_gene_rate + same_gene_multitranscript_rate,
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

    fn add_all(&mut self, next: FunctionTimings) {
        for (name, duration) in next.totals {
            *self.totals.entry(name).or_default() += duration;
        }
    }

    fn into_seconds(self) -> BTreeMap<String, f64> {
        self.totals
            .into_iter()
            .map(|(name, duration)| (name.to_owned(), duration.as_secs_f64()))
            .collect()
    }
}

fn make_scorer(score_mode: ScoreMode, config: ScoreConfig) -> Box<dyn CandidateScorer> {
    match score_mode {
        ScoreMode::Scalar => Box::new(ScalarScorer { config }),
        ScoreMode::Pulp => Box::new(PulpScorer { config }),
    }
}

fn score_config(args: &MapArgs) -> ScoreConfig {
    ScoreConfig {
        max_mismatches: args.max_mismatches,
        max_right_softclip: args.max_right_softclip,
        max_left_softclip: args.max_left_softclip,
        trim_poly_a: args.trim_poly_a,
        trim_poly_t: args.trim_poly_t,
        trim_low_quality_tail: args.trim_low_quality_tail,
        trim_tso: args.trim_tso,
        max_tso_mismatches: args.max_tso_mismatches,
        min_tso_match_len: args.min_tso_match_len,
        library_strand: match args.library_strand {
            LibraryStrandArg::Unstranded => LibraryStrand::Unstranded,
            LibraryStrandArg::Forward => LibraryStrand::Forward,
            LibraryStrandArg::Reverse => LibraryStrand::Reverse,
        },
        min_scored_len: args.min_scored_length,
        min_tail_phred: args.min_tail_phred,
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

fn seed_planner(arg: SeedPlannerArg) -> SeedPlanner {
    match arg {
        SeedPlannerArg::RawFrequency => SeedPlanner::RawFrequency,
        SeedPlannerArg::GeneIdf => SeedPlanner::GeneIdf,
    }
}

fn candidate_search(args: &MapArgs) -> CandidateSearchMode {
    match args.candidate_search {
        CandidateSearchArg::Full => CandidateSearchMode::Full,
        CandidateSearchArg::SparseProbe => CandidateSearchMode::SparseProbe(SparseProbeConfig {
            stride: args.sparse_probe_stride,
            max_seeds: args.sparse_probe_max_seeds,
            min_seed_hits: args.sparse_probe_min_seed_hits,
        }),
    }
}

fn candidate_pruning(args: &MapArgs) -> CandidatePruningMode {
    match args.candidate_pruning {
        CandidatePruningArg::None => CandidatePruningMode::None,
        CandidatePruningArg::WandLite => CandidatePruningMode::WandLite(WandLiteConfig {
            min_top_seed_count: args.wand_min_top_seed_count,
            score_ratio_percent: args.wand_score_ratio_percent,
            max_loci_per_gene: args.wand_max_loci_per_gene,
        }),
        CandidatePruningArg::GeneWandLite => CandidatePruningMode::GeneWandLite(WandLiteConfig {
            min_top_seed_count: args.wand_min_top_seed_count,
            score_ratio_percent: args.wand_score_ratio_percent,
            max_loci_per_gene: args.wand_max_loci_per_gene,
        }),
    }
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
    score_stats: &mut ScoreFailureStats,
    stage_times: &mut MappingStageTimes,
    function_timings: &mut FunctionTimings,
    bucket_sizes: &mut Vec<usize>,
    read_candidate_stats: &mut [CandidateGenerationStats],
    read_score_stats: &mut [ScoreFailureStats],
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
            seed_planner(args.seed_planner),
        );
        stage_times.candidate_generation_seconds += candidate_started.elapsed().as_secs_f64();
        function_timings.add_elapsed(
            "generate_candidate_hits_with_quality_stats",
            candidate_started,
        );
        add_stats(stats, read_stats);
        if let Some(slot) = read_candidate_stats.get_mut(*read_id as usize) {
            *slot = read_stats;
        }
        let bucket_started = Instant::now();
        let buckets = make_single_hit_buckets(hits, args.candidate_bin_size);
        stage_times.bucket_build_seconds += bucket_started.elapsed().as_secs_f64();
        function_timings.add_elapsed("make_single_hit_buckets", bucket_started);
        bucket_sizes.extend(buckets.iter().map(|bucket| bucket.hits.len()));
        let mut scored = Vec::new();
        let scoring_started = Instant::now();
        for bucket in &buckets {
            let score_bucket_started = Instant::now();
            let bucket_stats = scorer.score_bucket(
                reads,
                quals,
                index,
                bucket,
                &mut scored,
                Some(read_score_stats),
            );
            add_score_stats(score_stats, bucket_stats);
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
    score_stats: &mut ScoreFailureStats,
    stage_times: &mut MappingStageTimes,
    function_timings: &mut FunctionTimings,
    bucket_sizes: &mut Vec<usize>,
    read_candidate_stats: &mut [CandidateGenerationStats],
    read_score_stats: &mut [ScoreFailureStats],
) -> Vec<Assignment> {
    if args.retrieval_mode == RetrievalModeArg::SeedBatched && early_stop.is_none() {
        return map_seed_batched_streamed(
            args,
            reads,
            quals,
            barcode_umis,
            low_quality,
            invalid_barcode_umi,
            sketch_flags_by_read,
            sketches,
            index,
            scorer,
            regroup_by_locus,
            stats,
            score_stats,
            stage_times,
            function_timings,
            bucket_sizes,
            read_candidate_stats,
            read_score_stats,
        );
    }

    let candidate_started = Instant::now();
    let mut all_hits = Vec::new();
    let mut low_complexity = Vec::new();
    let mut low_quality_reads = Vec::new();
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
                seed_planner(args.seed_planner),
            );
            CandidateWork::Hits {
                read_id: sketch.read_id,
                hits,
                read_stats,
            }
        })
        .collect();

    for work in candidate_work {
        match work {
            CandidateWork::Hits {
                read_id,
                hits,
                read_stats,
            } => {
                add_stats(stats, read_stats);
                if let Some(slot) = read_candidate_stats.get_mut(read_id as usize) {
                    *slot = read_stats;
                }
                all_hits.extend(hits);
            }
            CandidateWork::LowQuality(read_id) => low_quality_reads.push(read_id),
            CandidateWork::LowComplexity(read_id) => low_complexity.push(read_id),
        }
    }
    stage_times.candidate_generation_seconds += candidate_started.elapsed().as_secs_f64();
    function_timings.add_elapsed(
        "generate_candidate_hits_with_quality_stats_parallel",
        candidate_started,
    );
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
    let (by_read, bucket_score_stats) = score_buckets(
        reads,
        quals,
        index,
        &buckets,
        scorer,
        function_timings,
        if args.emit_unmapped_diagnostics.is_some() {
            Some(read_score_stats)
        } else {
            None
        },
    );
    add_score_stats(score_stats, bucket_score_stats);
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

#[allow(clippy::too_many_arguments)]
fn map_seed_batched_streamed(
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
    regroup_by_locus: bool,
    stats: &mut CandidateGenerationStats,
    score_stats: &mut ScoreFailureStats,
    stage_times: &mut MappingStageTimes,
    function_timings: &mut FunctionTimings,
    bucket_sizes: &mut Vec<usize>,
    read_candidate_stats: &mut [CandidateGenerationStats],
    read_score_stats: &mut [ScoreFailureStats],
) -> Vec<Assignment> {
    let mut active_read_ids = Vec::new();
    let mut assignments: Vec<Option<Assignment>> = vec![None; reads.len()];

    let flag_assignment_started = Instant::now();
    for sketch in sketches {
        let read_id = sketch.read_id;
        let (cell_barcode, umi) = barcode_umis[read_id as usize].clone();
        if low_quality[read_id as usize] || invalid_barcode_umi[read_id as usize] {
            assignments[read_id as usize] = Some(flagged_assignment(
                read_id,
                cell_barcode,
                umi,
                AssignmentType::LowQuality,
                if invalid_barcode_umi[read_id as usize] {
                    4
                } else {
                    2
                },
            ));
        } else if sketch.flags & 1 != 0 {
            assignments[read_id as usize] = Some(flagged_assignment(
                read_id,
                cell_barcode,
                umi,
                AssignmentType::LowComplexity,
                sketch_flags_by_read[read_id as usize],
            ));
        } else {
            active_read_ids.push(read_id);
        }
    }
    stage_times.assignment_seconds += flag_assignment_started.elapsed().as_secs_f64();
    function_timings.add_elapsed("flagged_assignments_streamed", flag_assignment_started);

    if args.read_chunking == ReadChunkingArg::BestSeed {
        let planner = seed_planner(args.seed_planner);
        let best_seed_sort_started = Instant::now();
        let mut keyed_read_ids: Vec<_> = active_read_ids
            .par_iter()
            .map(|&read_id| {
                let (_, seq) = &reads[read_id as usize];
                (
                    best_sparse_probe_seed_key(
                        index,
                        seq,
                        Some(&quals[read_id as usize].1),
                        args.min_seed_quality,
                        planner,
                        args.sparse_probe_stride,
                    ),
                    read_id,
                )
            })
            .collect();
        keyed_read_ids.par_sort_unstable();
        active_read_ids = keyed_read_ids
            .into_iter()
            .map(|(_, read_id)| read_id)
            .collect();
        function_timings.add_elapsed("sort_reads_by_best_sparse_seed", best_seed_sort_started);
    }

    let chunk_size = args.batch_size.max(1);
    if args.emit_unmapped_diagnostics.is_none() {
        let chunk_results: Vec<_> = active_read_ids
            .par_chunks(chunk_size)
            .map(|chunk_read_ids| {
                map_seed_batched_chunk(
                    args,
                    reads,
                    quals,
                    barcode_umis,
                    chunk_read_ids,
                    index,
                    scorer,
                    regroup_by_locus,
                    None,
                )
            })
            .collect();
        for chunk_result in chunk_results {
            merge_streamed_chunk_result(
                chunk_result,
                stats,
                score_stats,
                stage_times,
                function_timings,
                bucket_sizes,
                read_candidate_stats,
                &mut assignments,
            );
        }
    } else {
        for chunk_read_ids in active_read_ids.chunks(chunk_size) {
            let chunk_result = map_seed_batched_chunk(
                args,
                reads,
                quals,
                barcode_umis,
                chunk_read_ids,
                index,
                scorer,
                regroup_by_locus,
                Some(&mut *read_score_stats),
            );
            merge_streamed_chunk_result(
                chunk_result,
                stats,
                score_stats,
                stage_times,
                function_timings,
                bucket_sizes,
                read_candidate_stats,
                &mut assignments,
            );
        }
    }

    assignments
        .into_iter()
        .enumerate()
        .map(|(read_id, assignment)| {
            assignment.unwrap_or_else(|| {
                let (cell_barcode, umi) = barcode_umis[read_id].clone();
                assign_read(read_id as u64, cell_barcode, umi, &[])
            })
        })
        .collect()
}

struct StreamedChunkResult {
    assignments: Vec<(u64, Assignment)>,
    stats: CandidateGenerationStats,
    per_read_stats: Vec<(u64, CandidateGenerationStats)>,
    score_stats: ScoreFailureStats,
    stage_times: MappingStageTimes,
    function_timings: FunctionTimings,
    bucket_sizes: Vec<usize>,
}

#[allow(clippy::too_many_arguments)]
fn map_seed_batched_chunk(
    args: &MapArgs,
    reads: &[(u64, Vec<u8>)],
    quals: &[(u64, Vec<u8>)],
    barcode_umis: &[(String, String)],
    chunk_read_ids: &[u64],
    index: &dyn IndexAccess,
    scorer: &dyn CandidateScorer,
    regroup_by_locus: bool,
    read_score_stats: Option<&mut [ScoreFailureStats]>,
) -> StreamedChunkResult {
    let mut function_timings = FunctionTimings::default();
    let mut stage_times = MappingStageTimes::default();

    let candidate_started = Instant::now();
    let (hits, stats, per_read_stats) = generate_candidate_hits_seed_batched(
        index,
        reads,
        quals,
        chunk_read_ids,
        args.min_seed_quality,
        args.max_seeds_per_read,
        args.max_postings_per_seed,
        args.search_reverse_complement,
        seed_planner(args.seed_planner),
        candidate_search(args),
        candidate_pruning(args),
    );
    stage_times.candidate_generation_seconds += candidate_started.elapsed().as_secs_f64();
    function_timings.add_elapsed("generate_candidate_hits_seed_batched", candidate_started);

    let bucket_started = Instant::now();
    let buckets = if regroup_by_locus {
        make_candidate_locus_buckets(hits, args.candidate_bin_size)
    } else {
        make_single_hit_buckets(hits, args.candidate_bin_size)
    };
    stage_times.bucket_build_seconds += bucket_started.elapsed().as_secs_f64();
    if regroup_by_locus {
        function_timings.add_elapsed("make_candidate_locus_buckets", bucket_started);
    } else {
        function_timings.add_elapsed("make_single_hit_buckets", bucket_started);
    }
    let bucket_sizes = buckets.iter().map(|bucket| bucket.hits.len()).collect();

    let scoring_started = Instant::now();
    let (scored, score_stats) = score_buckets_flat(
        reads,
        quals,
        index,
        &buckets,
        scorer,
        &mut function_timings,
        read_score_stats,
    );
    stage_times.scoring_seconds += scoring_started.elapsed().as_secs_f64();

    let assignment_started = Instant::now();
    let assignments = make_streamed_chunk_assignments(chunk_read_ids, &scored, barcode_umis);
    stage_times.assignment_seconds += assignment_started.elapsed().as_secs_f64();
    function_timings.add_elapsed("assign_streamed_chunk", assignment_started);

    StreamedChunkResult {
        assignments,
        stats,
        per_read_stats,
        score_stats,
        stage_times,
        function_timings,
        bucket_sizes,
    }
}

fn make_streamed_chunk_assignments(
    chunk_read_ids: &[u64],
    scored: &[ScoredCandidate],
    barcode_umis: &[(String, String)],
) -> Vec<(u64, Assignment)> {
    let mut assignments = Vec::with_capacity(chunk_read_ids.len());
    for &read_id in chunk_read_ids {
        let start = scored.partition_point(|candidate| candidate.read_id < read_id);
        let end = start + scored[start..].partition_point(|candidate| candidate.read_id == read_id);
        let (cell_barcode, umi) = barcode_umis[read_id as usize].clone();
        assignments.push((
            read_id,
            assign_read(read_id, cell_barcode, umi, &scored[start..end]),
        ));
    }
    assignments
}

#[allow(clippy::too_many_arguments)]
fn merge_streamed_chunk_result(
    chunk_result: StreamedChunkResult,
    stats: &mut CandidateGenerationStats,
    score_stats: &mut ScoreFailureStats,
    stage_times: &mut MappingStageTimes,
    function_timings: &mut FunctionTimings,
    bucket_sizes: &mut Vec<usize>,
    read_candidate_stats: &mut [CandidateGenerationStats],
    assignments: &mut [Option<Assignment>],
) {
    add_stats(stats, chunk_result.stats);
    add_score_stats(score_stats, chunk_result.score_stats);
    stage_times.candidate_generation_seconds +=
        chunk_result.stage_times.candidate_generation_seconds;
    stage_times.bucket_build_seconds += chunk_result.stage_times.bucket_build_seconds;
    stage_times.scoring_seconds += chunk_result.stage_times.scoring_seconds;
    stage_times.assignment_seconds += chunk_result.stage_times.assignment_seconds;
    function_timings.add_all(chunk_result.function_timings);
    bucket_sizes.extend(chunk_result.bucket_sizes);
    for (read_id, read_stats) in chunk_result.per_read_stats {
        if let Some(slot) = read_candidate_stats.get_mut(read_id as usize) {
            *slot = read_stats;
        }
    }
    for (read_id, assignment) in chunk_result.assignments {
        assignments[read_id as usize] = Some(assignment);
    }
}

enum CandidateWork {
    Hits {
        read_id: u64,
        hits: Vec<constellation_core::candidate::CandidateHit>,
        read_stats: CandidateGenerationStats,
    },
    LowQuality(u64),
    LowComplexity(u64),
}

fn score_buckets(
    reads: &[(u64, Vec<u8>)],
    quals: &[(u64, Vec<u8>)],
    index: &dyn IndexAccess,
    buckets: &[CandidateLocusBucket],
    scorer: &dyn CandidateScorer,
    function_timings: &mut FunctionTimings,
    read_score_stats: Option<&mut [ScoreFailureStats]>,
) -> (Vec<Vec<ScoredCandidate>>, ScoreFailureStats) {
    let (scored, score_stats) = score_buckets_flat(
        reads,
        quals,
        index,
        buckets,
        scorer,
        function_timings,
        read_score_stats,
    );

    let group_started = Instant::now();
    let mut by_read = vec![Vec::new(); reads.len()];
    for candidate in scored {
        if let Some(candidates) = by_read.get_mut(candidate.read_id as usize) {
            candidates.push(candidate);
        }
    }
    function_timings.add_elapsed("group_scored_candidates_by_read", group_started);
    (by_read, score_stats)
}

fn score_buckets_flat(
    reads: &[(u64, Vec<u8>)],
    quals: &[(u64, Vec<u8>)],
    index: &dyn IndexAccess,
    buckets: &[CandidateLocusBucket],
    scorer: &dyn CandidateScorer,
    function_timings: &mut FunctionTimings,
    read_score_stats: Option<&mut [ScoreFailureStats]>,
) -> (Vec<ScoredCandidate>, ScoreFailureStats) {
    if read_score_stats.is_none() {
        let score_bucket_started = Instant::now();
        let scored_work: Vec<_> = buckets
            .par_iter()
            .map(|bucket| {
                let mut scored = Vec::new();
                let stats = scorer.score_bucket(reads, quals, index, bucket, &mut scored, None);
                (scored, stats)
            })
            .collect();
        function_timings.add_elapsed(
            "CandidateScorer::score_bucket_parallel",
            score_bucket_started,
        );

        let mut scored = Vec::new();
        let mut score_stats = ScoreFailureStats::default();
        for (mut bucket_scored, bucket_stats) in scored_work {
            scored.append(&mut bucket_scored);
            add_score_stats(&mut score_stats, bucket_stats);
        }
        let sort_started = Instant::now();
        scored.sort();
        function_timings.add_elapsed("ScoredCandidate::sort", sort_started);
        return (scored, score_stats);
    }

    let read_score_stats = read_score_stats.expect("checked above");
    let mut scored = Vec::new();
    let mut score_stats = ScoreFailureStats::default();
    for bucket in buckets {
        let score_bucket_started = Instant::now();
        let bucket_stats = scorer.score_bucket(
            reads,
            quals,
            index,
            bucket,
            &mut scored,
            Some(read_score_stats),
        );
        add_score_stats(&mut score_stats, bucket_stats);
        function_timings.add_elapsed("CandidateScorer::score_bucket", score_bucket_started);
    }
    let sort_started = Instant::now();
    scored.sort();
    function_timings.add_elapsed("ScoredCandidate::sort", sort_started);
    (scored, score_stats)
}

fn add_stats(total: &mut CandidateGenerationStats, next: CandidateGenerationStats) {
    total.seed_lookups += next.seed_lookups;
    total.seed_candidates_considered += next.seed_candidates_considered;
    total.candidate_hits += next.candidate_hits;
    total.postings_skipped_due_to_frequency += next.postings_skipped_due_to_frequency;
    total.seeds_skipped_due_to_quality += next.seeds_skipped_due_to_quality;
    total.early_stopped_reads += next.early_stopped_reads;
    total.seed_lookups_saved_by_early_stop += next.seed_lookups_saved_by_early_stop;
    total.selected_seeds += next.selected_seeds;
    total.reads_with_no_valid_kmers += next.reads_with_no_valid_kmers;
    total.reads_with_no_selected_seeds += next.reads_with_no_selected_seeds;
    total.reads_with_no_seed_postings += next.reads_with_no_seed_postings;
    total.reads_all_selected_seeds_over_frequency_cap +=
        next.reads_all_selected_seeds_over_frequency_cap;
    total.reads_with_candidate_hits += next.reads_with_candidate_hits;
    total.selected_seeds_absent_from_index += next.selected_seeds_absent_from_index;
    total.selected_seeds_over_frequency_cap += next.selected_seeds_over_frequency_cap;
    total.selected_seeds_with_postings += next.selected_seeds_with_postings;
    total.query_seed_occurrences += next.query_seed_occurrences;
    total.distinct_seed_lists_loaded += next.distinct_seed_lists_loaded;
    total.candidate_votes += next.candidate_votes;
    total.seed_groups_skipped_due_to_frequency += next.seed_groups_skipped_due_to_frequency;
    total.sparse_probe_attempted_reads += next.sparse_probe_attempted_reads;
    total.sparse_probe_accepted_reads += next.sparse_probe_accepted_reads;
    total.sparse_probe_fallback_reads += next.sparse_probe_fallback_reads;
    total.candidate_hits_pruned += next.candidate_hits_pruned;
}

fn add_score_stats(total: &mut ScoreFailureStats, next: ScoreFailureStats) {
    total.candidates_seen += next.candidates_seen;
    total.candidates_out_of_bounds += next.candidates_out_of_bounds;
    total.candidates_failed_mismatch += next.candidates_failed_mismatch;
    total.candidates_failed_min_scored_len += next.candidates_failed_min_scored_len;
    total.candidates_failed_trimmed_out_of_bounds += next.candidates_failed_trimmed_out_of_bounds;
    total.candidates_full_length_mismatch_only += next.candidates_full_length_mismatch_only;
    total.candidates_passed_full_length += next.candidates_passed_full_length;
    total.candidates_passed_trimmed_or_softclipped += next.candidates_passed_trimmed_or_softclipped;
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

fn unmapped_diagnostics_tsv(
    assignments: &[Assignment],
    reads: &[(u64, Vec<u8>)],
    read_candidate_stats: &[CandidateGenerationStats],
    read_score_stats: &[ScoreFailureStats],
) -> String {
    let mut out = String::from(
        "read_id\treason\tseq_len\tnum_valid_kmers\tselected_seeds\tseed_hits\tcandidate_hits\tscored_candidates\tscore_seen\tscore_oob\tscore_mismatch\tscore_min_len\tscore_trimmed_oob\tscore_full_length_mismatch_only\tflags\n",
    );
    for assignment in assignments {
        if matches!(
            assignment.assignment_type,
            AssignmentType::UniqueGene | AssignmentType::AmbiguousTranscriptSameGene
        ) {
            continue;
        }
        let read_id = assignment.read_id as usize;
        let seq_len = reads.get(read_id).map_or(0, |(_, seq)| seq.len());
        let stats = read_candidate_stats
            .get(read_id)
            .copied()
            .unwrap_or_default();
        let score_stats = read_score_stats.get(read_id).copied().unwrap_or_default();
        out.push_str(&format!(
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n",
            assignment.read_id,
            diagnostic_reason(assignment, stats, score_stats),
            seq_len,
            stats.seed_candidates_considered,
            stats.selected_seeds,
            stats.selected_seeds_with_postings,
            stats.candidate_hits,
            assignment.candidate_count,
            score_stats.candidates_seen,
            score_stats.candidates_out_of_bounds,
            score_stats.candidates_failed_mismatch,
            score_stats.candidates_failed_min_scored_len,
            score_stats.candidates_failed_trimmed_out_of_bounds,
            score_stats.candidates_full_length_mismatch_only,
            assignment.flags
        ));
    }
    out
}

fn diagnostic_reason(
    assignment: &Assignment,
    stats: CandidateGenerationStats,
    score_stats: ScoreFailureStats,
) -> &'static str {
    match assignment.assignment_type {
        AssignmentType::LowQuality if assignment.flags & 4 != 0 => "invalid_barcode_umi",
        AssignmentType::LowQuality => "low_quality",
        AssignmentType::LowComplexity => "low_complexity",
        AssignmentType::Unmapped if stats.reads_with_no_valid_kmers > 0 => "no_valid_kmers",
        AssignmentType::Unmapped if stats.reads_with_no_selected_seeds > 0 => "no_selected_seeds",
        AssignmentType::Unmapped if stats.reads_all_selected_seeds_over_frequency_cap > 0 => {
            "all_selected_seeds_over_frequency_cap"
        }
        AssignmentType::Unmapped if stats.selected_seeds_with_postings == 0 => {
            "no_seed_with_postings"
        }
        AssignmentType::Unmapped
            if stats.candidate_hits > 0
                && assignment.candidate_count == 0
                && score_stats.candidates_failed_min_scored_len > 0 =>
        {
            "score_failed_min_scored_len"
        }
        AssignmentType::Unmapped
            if stats.candidate_hits > 0
                && assignment.candidate_count == 0
                && score_stats.candidates_failed_trimmed_out_of_bounds > 0 =>
        {
            "score_failed_trimmed_out_of_bounds"
        }
        AssignmentType::Unmapped
            if stats.candidate_hits > 0
                && assignment.candidate_count == 0
                && score_stats.candidates_out_of_bounds > 0 =>
        {
            "score_failed_full_length_out_of_bounds"
        }
        AssignmentType::Unmapped
            if stats.candidate_hits > 0
                && assignment.candidate_count == 0
                && score_stats.candidates_failed_mismatch > 0 =>
        {
            "score_failed_mismatch"
        }
        AssignmentType::Unmapped
            if stats.candidate_hits > 0
                && assignment.candidate_count == 0
                && score_stats.candidates_full_length_mismatch_only > 0 =>
        {
            "score_failed_full_length_mismatch"
        }
        AssignmentType::Unmapped if stats.candidate_hits > 0 && assignment.candidate_count == 0 => {
            "candidate_hits_but_score_failed"
        }
        AssignmentType::Unmapped => "unmapped_unknown",
        AssignmentType::AmbiguousGene => "multi_gene_ambiguous",
        AssignmentType::AntisenseGene => "antisense_gene",
        AssignmentType::UniqueGene | AssignmentType::AmbiguousTranscriptSameGene => {
            "gene_countable"
        }
    }
}

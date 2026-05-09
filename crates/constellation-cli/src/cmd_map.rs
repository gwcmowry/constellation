use crate::{
    CandidatePruningArg, CandidateSearchArg, LibraryStrandArg, MapArgs, MapMode,
    MapOutputFormatArg, OutputCompressionArg, ReadChunkingArg, RetrievalModeArg, ScoreMode,
    SeedPlannerArg,
};
use anyhow::{anyhow, Context, Result};
use constellation_core::assign::{assign_read, flagged_assignment, Assignment, AssignmentType};
use constellation_core::candidate::{
    best_sparse_probe_seed_key, generate_candidate_hits_seed_batched,
    generate_candidate_hits_with_quality_stats, make_candidate_locus_buckets,
    make_single_hit_buckets, CandidateEarlyStopConfig, CandidateGenerationStats,
    CandidateLocusBucket, CandidatePruningMode, CandidateSearchMode, SeedPlanner,
    SparseProbeConfig, WandLiteConfig,
};
use constellation_core::chemistry::{parse_tenx_3p_v3_r1, BarcodeUmi, Chemistry};
use constellation_core::dna::{encode_acgt, iter_kmers_2bit, reverse_complement};
use constellation_core::ec_index::{EcGeneScore, EcMappingStats, LoadedEcIndex};
use constellation_core::fastq::read_fastq;
use constellation_core::index::{IndexAccess, LoadedIndex};
use constellation_core::metrics::MapMetrics;
use constellation_core::score::{
    tso_prefix_trim_len, CandidateScorer, LibraryStrand, ScoreConfig, ScoreFailureStats,
    ScoredCandidate, SCORE_FLAG_ANTISENSE, SCORE_FLAG_RESCUED, SCORE_FLAG_TARGET_EXON,
    SCORE_FLAG_TARGET_GENE_BODY, SCORE_FLAG_TARGET_INTRON,
};
use constellation_core::score_scalar::{target_class_flag, ScalarScorer};
use constellation_core::score_simd::PulpScorer;
use constellation_core::sketch::{sketch_read, sort_sketch_records};
use flate2::read::MultiGzDecoder;
use rayon::prelude::*;
use std::collections::{BTreeMap, HashMap};
use std::fs::{self, File};
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::path::{Path, PathBuf};
use std::process::{Child, ChildStdin, Command, Stdio};
use std::sync::mpsc::{self, Receiver};
use std::thread::{self, JoinHandle};
use std::time::{Duration, Instant};

pub fn run_map(args: MapArgs) -> Result<()> {
    let _chemistry = Chemistry::parse(&args.chemistry)?;
    if args.mode != MapMode::GeneEc {
        return Err(anyhow!(
            "legacy positional map mode '{}' is no longer supported by the default CLI; use --mode gene-ec with an EC index",
            args.mode.as_str()
        ));
    }
    if args.mode == MapMode::GeneEc {
        return run_map_gene_ec(args);
    }
    let mut function_timings = FunctionTimings::default();
    let total_started = Instant::now();
    let index_started = Instant::now();
    let index = LoadedIndex::load(&args.index)?;
    function_timings.add_elapsed("LoadedIndex::load", index_started);
    let hot_index = if let Some(path) = args.hot_index.as_ref() {
        let hot_index_started = Instant::now();
        let loaded = LoadedIndex::load(path)?;
        function_timings.add_elapsed("LoadedIndex::load_hot", hot_index_started);
        let _ = build_two_tier_gene_map(&loaded, &index)?;
        Some(loaded)
    } else {
        None
    };
    let index_load_seconds = index_started.elapsed().as_secs_f64();
    let primary_index: &dyn IndexAccess = hot_index
        .as_ref()
        .map(|loaded| loaded as &dyn IndexAccess)
        .unwrap_or(&index);
    let cold_to_primary_gene = hot_index
        .as_ref()
        .map(|hot| build_two_tier_gene_map(hot, &index))
        .transpose()?;
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
            let sketch = sketch_read(read_id, &seq, primary_index.k());
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
    let scorer = make_scorer(args.score_mode, score_config, primary_index);
    let cold_scorer = hot_index
        .is_some()
        .then(|| make_scorer(args.score_mode, score_config, &index));
    let mut cold_preload_seconds = 0.0;
    if hot_index.is_some() && args.preload_cold_index {
        let preload_started = Instant::now();
        index.preload()?;
        cold_preload_seconds = preload_started.elapsed().as_secs_f64();
        function_timings.add_elapsed("LoadedIndex::preload_cold", preload_started);
    }
    let early_stop = early_stop_config(&args);
    let mut stats = CandidateGenerationStats::default();
    let mut score_stats = ScoreFailureStats::default();
    let mut stage_times = MappingStageTimes::default();
    let mut bucket_sizes = Vec::new();
    let mut read_candidate_stats = vec![CandidateGenerationStats::default(); reads.len()];
    let mut read_score_stats = vec![ScoreFailureStats::default(); reads.len()];
    let mapping_started = Instant::now();
    let primary_mapping_started = Instant::now();
    let mut assignments = match args.mode {
        MapMode::ReadAtATime => map_read_at_a_time(
            &args,
            &reads,
            &quals,
            &barcode_umis,
            &low_quality,
            &invalid_barcode_umi,
            &sketch_flags_by_read,
            primary_index,
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
                primary_index,
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
                primary_index,
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
        MapMode::GeneEc => unreachable!("gene-ec mode returns before positional index mapping"),
    };
    let hot_mapping_seconds = primary_mapping_started.elapsed().as_secs_f64();
    let mut cold_fallback_seconds = 0.0;
    let mut cold_fallback_reads = 0_u64;
    let mut cold_fallback_replaced_reads = 0_u64;
    let mut cold_fallback_shard_groups = 0_u64;
    if let (Some(_), Some(cold_scorer)) = (hot_index.as_ref(), cold_scorer.as_ref()) {
        let cold_started = Instant::now();
        let cold_summary = run_cold_fallback(
            &args,
            &reads,
            &quals,
            &barcode_umis,
            &index,
            cold_scorer.as_ref(),
            cold_to_primary_gene.as_deref(),
            true,
            &mut stats,
            &mut score_stats,
            &mut stage_times,
            &mut function_timings,
            &mut bucket_sizes,
            &mut read_candidate_stats,
            if args.emit_unmapped_diagnostics.is_some() {
                Some(&mut read_score_stats)
            } else {
                None
            },
            &mut assignments,
        );
        cold_fallback_seconds = cold_started.elapsed().as_secs_f64();
        cold_fallback_reads = cold_summary.fallback_reads;
        cold_fallback_replaced_reads = cold_summary.replaced_reads;
        cold_fallback_shard_groups = cold_summary.shard_groups;
    }
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
            hot_index_enabled: hot_index.is_some(),
            hot_mapping_seconds,
            cold_fallback_seconds,
            cold_fallback_reads,
            cold_fallback_replaced_reads,
            cold_fallback_shard_groups,
            cold_shard_prefix_bits: args.cold_shard_prefix_bits,
            cold_evict_interval: args.cold_evict_interval,
            cold_shard_scheduling: args.cold_shard_scheduling,
            cold_preload_seconds,
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
            exon_target_assignment_rate: target_assignment_rate(
                &assignments,
                SCORE_FLAG_TARGET_EXON,
            ),
            intron_target_assignment_rate: target_assignment_rate(
                &assignments,
                SCORE_FLAG_TARGET_INTRON,
            ),
            gene_body_target_assignment_rate: target_assignment_rate(
                &assignments,
                SCORE_FLAG_TARGET_GENE_BODY,
            ),
            mixed_target_assignment_rate: mixed_target_assignment_rate(&assignments),
            unknown_target_assignment_rate: unknown_target_assignment_rate(&assignments),
            rescued_assignment_rate: target_assignment_rate(&assignments, SCORE_FLAG_RESCUED),
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

fn run_map_gene_ec(args: MapArgs) -> Result<()> {
    let mut function_timings = FunctionTimings::default();
    let total_started = Instant::now();

    let index_started = Instant::now();
    let index = LoadedEcIndex::load(&args.index)?;
    function_timings.add_elapsed("LoadedEcIndex::load", index_started);
    let index_load_seconds = index_started.elapsed().as_secs_f64();

    let stream_batch_size = args.batch_size.max(65_536);
    let mut fastq_pipeline =
        FastqPairBatchPipeline::new(args.r1.clone(), args.r2.clone(), stream_batch_size);
    let mut assignment_writer = if args.skip_assignments {
        None
    } else {
        Some(GeneEcOutputWriter::create(
            args.output_format,
            args.output_compression,
            args.zstd_level,
            &args.out,
        )?)
    };
    let collect_diagnostics = args.emit_unmapped_diagnostics.is_some();
    let mut diagnostic_assignments = Vec::new();
    let mut diagnostic_reads = Vec::new();
    let mut diagnostic_candidate_stats = Vec::new();
    let mut diagnostic_score_stats = Vec::new();
    let mut ec_stats = EcMappingStats::default();
    let mut assignment_counts = AssignmentCounts::default();
    let mut total_gene_candidates = 0_u64;
    let mut num_reads = 0_u64;
    let mut total_bases = 0_u64;
    let mut fastq_load_seconds = 0.0;
    let preprocess_seconds = 0.0;
    let mut mapping_seconds = 0.0;
    let mut assignment_seconds = 0.0;
    let mut write_seconds = 0.0;
    let sort_by_locality = matches!(args.read_chunking, ReadChunkingArg::BestSeed);

    loop {
        let Some((batch, read_elapsed)) = fastq_pipeline.next_batch()? else {
            break;
        };
        fastq_load_seconds += read_elapsed.as_secs_f64();
        function_timings.add_duration("stream_fastq_pair_batches", read_elapsed);

        let map_started = Instant::now();
        let mapped: Vec<_> = if sort_by_locality {
            let mut prepared_batch: Vec<_> = batch
                .par_iter()
                .map(|raw| {
                    prepare_gene_ec_raw(
                        &index,
                        raw,
                        args.trim_tso,
                        args.min_tso_match_len,
                        args.max_tso_mismatches,
                        args.min_mean_quality,
                        args.sparse_probe_stride,
                        true,
                    )
                })
                .collect();
            prepared_batch.sort_unstable_by_key(|read| (read.locality_key, read.prepared.read_id));
            let mut mapped: Vec<_> = prepared_batch
                .par_iter()
                .map(|prepared| {
                    map_prepared_gene_ec(
                        &index,
                        prepared,
                        args.search_reverse_complement,
                        args.min_mean_quality,
                        args.sparse_probe_stride,
                        args.max_postings_per_seed as u32,
                        collect_diagnostics,
                        matches!(args.output_format, MapOutputFormatArg::GeneEcRad),
                    )
                })
                .collect();
            mapped.sort_unstable_by_key(|mapped_read| mapped_read.assignment.read_id);
            mapped
        } else {
            batch
                .par_iter()
                .map(|raw| {
                    prepare_and_map_gene_ec_raw(
                        &index,
                        raw,
                        args.trim_tso,
                        args.min_tso_match_len,
                        args.max_tso_mismatches,
                        args.min_mean_quality,
                        args.search_reverse_complement,
                        args.sparse_probe_stride,
                        args.max_postings_per_seed as u32,
                        collect_diagnostics,
                        matches!(args.output_format, MapOutputFormatArg::GeneEcRad),
                    )
                })
                .collect()
        };
        let map_elapsed = map_started.elapsed();
        mapping_seconds += map_elapsed.as_secs_f64();
        function_timings.add_duration("prepare_and_map_gene_ec_batches", map_elapsed);

        let assignment_started = Instant::now();
        let mut batch_tsv = if assignment_writer.is_some() {
            Some(Vec::with_capacity(output_batch_capacity(
                args.output_format,
                mapped.len(),
            )))
        } else {
            None
        };
        for mapped_read in mapped {
            add_ec_stats(&mut ec_stats, mapped_read.stats);
            total_gene_candidates += mapped_read.gene_candidate_count as u64;
            total_bases += mapped_read.seq_len as u64;
            num_reads += 1;
            assignment_counts.add(&mapped_read.assignment);
            write_gene_ec_record_to_buffer(&mut batch_tsv, args.output_format, &mapped_read)?;
            if collect_diagnostics {
                diagnostic_reads.push((mapped_read.assignment.read_id, mapped_read.seq));
                diagnostic_candidate_stats.push(mapped_read.candidate_stats);
                diagnostic_score_stats.push(ScoreFailureStats::default());
                diagnostic_assignments.push(mapped_read.assignment);
            }
        }
        let assignment_elapsed = assignment_started.elapsed();
        assignment_seconds += assignment_elapsed.as_secs_f64();
        function_timings.add_duration("assign_gene_ec_reads", assignment_elapsed);
        if let (Some(writer), Some(tsv)) = (assignment_writer.as_mut(), batch_tsv) {
            let write_started = Instant::now();
            writer.write_all(&tsv)?;
            write_seconds += write_started.elapsed().as_secs_f64();
        }
    }

    fastq_pipeline.finish()?;
    if let Some(writer) = assignment_writer {
        let flush_started = Instant::now();
        writer.finish()?;
        write_seconds += flush_started.elapsed().as_secs_f64();
    }

    if let Some(path) = args.emit_unmapped_diagnostics.as_ref() {
        fs::write(
            path,
            unmapped_diagnostics_tsv(
                &diagnostic_assignments,
                &diagnostic_reads,
                &diagnostic_candidate_stats,
                &diagnostic_score_stats,
            ),
        )?;
    }

    if let Some(metrics_path) = args.emit_metrics {
        let unique_gene_rate = assignment_counts.rate(num_reads, AssignmentType::UniqueGene);
        let same_gene_multitranscript_rate =
            assignment_counts.rate(num_reads, AssignmentType::AmbiguousTranscriptSameGene);
        let multi_gene_ambiguous_rate =
            assignment_counts.rate(num_reads, AssignmentType::AmbiguousGene);
        let mut metrics = MapMetrics {
            mode: args.mode.as_str().to_owned(),
            score_mode: "gene-ec".to_owned(),
            num_reads,
            index_load_seconds,
            fastq_load_seconds,
            preprocess_seconds,
            mapping_seconds,
            candidate_generation_seconds: mapping_seconds,
            assignment_seconds,
            write_seconds,
            seed_lookups_per_read: per_read_u64(ec_stats.seed_lookups, num_reads),
            seed_candidates_considered_per_read: ec_stats.seed_lookups as f64
                / num_reads.max(1) as f64,
            candidate_hits_per_read: per_read_u64(total_gene_candidates, num_reads),
            reads_with_no_seed_postings: num_reads.saturating_sub(ec_stats.reads_with_gene_support),
            reads_with_candidate_hits: ec_stats.reads_with_gene_support,
            selected_seeds_absent_from_index: ec_stats
                .seed_lookups
                .saturating_sub(ec_stats.seeds_found),
            selected_seeds_over_frequency_cap: ec_stats.seeds_skipped_high_gene_df,
            selected_seeds_with_postings: ec_stats.seeds_found,
            query_seed_occurrences_per_read: ec_stats.seed_lookups as f64 / num_reads.max(1) as f64,
            distinct_seed_lists_loaded: ec_stats.seeds_found,
            posting_list_reuse_factor: 1.0,
            candidate_votes_per_read: per_read_u64(ec_stats.gene_votes, num_reads),
            scored_candidates_per_read: per_read_u64(total_gene_candidates, num_reads),
            unique_gene_rate,
            ambiguous_gene_rate: multi_gene_ambiguous_rate,
            same_gene_multitranscript_rate,
            multi_gene_ambiguous_rate,
            gene_countable_rate: unique_gene_rate + same_gene_multitranscript_rate,
            exon_target_assignment_rate: assignment_counts
                .flag_rate(num_reads, SCORE_FLAG_TARGET_EXON),
            intron_target_assignment_rate: assignment_counts
                .flag_rate(num_reads, SCORE_FLAG_TARGET_INTRON),
            mixed_target_assignment_rate: assignment_counts.mixed_target_rate(num_reads),
            unknown_target_assignment_rate: assignment_counts.unknown_target_rate(num_reads),
            low_complexity_rate: assignment_counts.rate(num_reads, AssignmentType::LowComplexity),
            low_quality_rate: assignment_counts.rate(num_reads, AssignmentType::LowQuality),
            unmapped_rate: assignment_counts.rate(num_reads, AssignmentType::Unmapped),
            function_seconds: function_timings.into_seconds(),
            ..MapMetrics::default()
        };
        metrics.finish_rates(total_started.elapsed(), total_bases);
        fs::write(metrics_path, serde_json::to_vec_pretty(&metrics)?)?;
    }

    Ok(())
}

fn map_gene_ec_read(
    index: &LoadedEcIndex,
    prepared: &PreparedRead,
    search_reverse_complement: bool,
    min_mean_quality: f64,
    sparse_probe_stride: u32,
    max_gene_df: u32,
    emit_top_gene_ids: bool,
) -> (Assignment, EcMappingStats, u32, Vec<u32>) {
    let cell_barcode = prepared.cell_barcode.clone();
    let umi = prepared.umi.clone();
    if prepared.low_quality || prepared.invalid_barcode_umi {
        return (
            flagged_assignment(
                prepared.read_id,
                cell_barcode,
                umi,
                AssignmentType::LowQuality,
                if prepared.invalid_barcode_umi { 4 } else { 2 },
            ),
            EcMappingStats::default(),
            0,
            Vec::new(),
        );
    }
    if prepared.sketch.flags & 1 != 0 {
        return (
            flagged_assignment(
                prepared.read_id,
                cell_barcode,
                umi,
                AssignmentType::LowComplexity,
                prepared.sketch.flags,
            ),
            EcMappingStats::default(),
            0,
            Vec::new(),
        );
    }
    if mean_phred_quality(&prepared.qual) < min_mean_quality {
        return (
            flagged_assignment(
                prepared.read_id,
                cell_barcode,
                umi,
                AssignmentType::LowQuality,
                2,
            ),
            EcMappingStats::default(),
            0,
            Vec::new(),
        );
    }

    let (mapping, stats) = index.map_read_sparse(
        &prepared.seq,
        search_reverse_complement,
        sparse_probe_stride,
        max_gene_df,
    );
    let gene_candidate_count = mapping.gene_scores.len().min(u32::MAX as usize) as u32;
    if mapping.best_score == 0 {
        return (
            Assignment {
                read_id: prepared.read_id,
                cell_barcode,
                umi,
                assignment_type: AssignmentType::Unmapped,
                gene_id: None,
                transcript_id: None,
                candidate_count: 0,
                score: 0,
                flags: 0,
            },
            stats,
            gene_candidate_count,
            Vec::new(),
        );
    }
    let top_gene_ids = if emit_top_gene_ids {
        top_gene_ids_for_ec(&mapping.gene_scores, mapping.best_score)
    } else {
        Vec::new()
    };
    let resolved = resolve_gene_ec_assignment(&mapping.gene_scores, mapping.best_score);
    (
        Assignment {
            read_id: prepared.read_id,
            cell_barcode,
            umi,
            assignment_type: resolved.assignment_type,
            gene_id: resolved.gene_id,
            transcript_id: None,
            candidate_count: gene_candidate_count,
            score: resolved.score,
            flags: resolved.flags,
        },
        stats,
        gene_candidate_count,
        top_gene_ids,
    )
}

fn top_gene_ids_for_ec(scores: &[EcGeneScore], best_score: u16) -> Vec<u32> {
    if best_score == 0 {
        return Vec::new();
    }
    scores
        .iter()
        .filter(|score| score.total_score == best_score)
        .map(|score| score.gene_id)
        .collect()
}

#[derive(Debug, Clone, Copy)]
struct GeneEcResolution {
    assignment_type: AssignmentType,
    gene_id: Option<u32>,
    score: u16,
    flags: u16,
}

fn resolve_gene_ec_assignment(
    scores: &[EcGeneScore],
    fallback_best_score: u16,
) -> GeneEcResolution {
    if let Some(score) = unique_top_by(scores, fallback_best_score, |score| score.exon_score) {
        return GeneEcResolution {
            assignment_type: AssignmentType::UniqueGene,
            gene_id: Some(score.gene_id),
            score: score.exon_score,
            flags: SCORE_FLAG_TARGET_EXON,
        };
    }
    if let Some(score) = unique_top_by(scores, fallback_best_score, |score| score.intron_score) {
        return GeneEcResolution {
            assignment_type: AssignmentType::UniqueGene,
            gene_id: Some(score.gene_id),
            score: score.intron_score,
            flags: SCORE_FLAG_TARGET_INTRON,
        };
    }
    if let Some(score) = unique_top_by(scores, fallback_best_score, |score| score.unknown_score) {
        return GeneEcResolution {
            assignment_type: AssignmentType::UniqueGene,
            gene_id: Some(score.gene_id),
            score: score.unknown_score,
            flags: 0,
        };
    }
    if let Some(score) = unique_top_by(scores, fallback_best_score, |score| score.total_score) {
        return GeneEcResolution {
            assignment_type: AssignmentType::UniqueGene,
            gene_id: Some(score.gene_id),
            score: score.total_score,
            flags: target_flags_for_gene_score(score),
        };
    }
    GeneEcResolution {
        assignment_type: AssignmentType::AmbiguousGene,
        gene_id: None,
        score: fallback_best_score,
        flags: aggregate_gene_ec_flags(scores),
    }
}

fn unique_top_by(
    scores: &[EcGeneScore],
    required_total_score: u16,
    value: impl Fn(&EcGeneScore) -> u16,
) -> Option<&EcGeneScore> {
    let max_score = scores
        .iter()
        .filter(|score| score.total_score == required_total_score)
        .map(&value)
        .max()
        .unwrap_or(0);
    if max_score == 0 {
        return None;
    }
    let mut top = scores
        .iter()
        .filter(|score| score.total_score == required_total_score && value(score) == max_score);
    let first = top.next()?;
    top.next().is_none().then_some(first)
}

fn target_flags_for_gene_score(score: &EcGeneScore) -> u16 {
    let mut flags = 0;
    if score.exon_score > 0 {
        flags |= SCORE_FLAG_TARGET_EXON;
    }
    if score.intron_score > 0 {
        flags |= SCORE_FLAG_TARGET_INTRON;
    }
    flags
}

fn aggregate_gene_ec_flags(scores: &[EcGeneScore]) -> u16 {
    scores.iter().fold(0, |flags, score| {
        flags
            | target_flags_for_gene_score(score)
            | if score.antisense_score > 0 {
                SCORE_FLAG_ANTISENSE
            } else {
                0
            }
    })
}

fn add_ec_stats(total: &mut EcMappingStats, next: EcMappingStats) {
    total.seed_lookups += next.seed_lookups;
    total.seeds_found += next.seeds_found;
    total.seeds_skipped_high_gene_df += next.seeds_skipped_high_gene_df;
    total.reads_with_gene_support += next.reads_with_gene_support;
    total.ec_hits += next.ec_hits;
    total.gene_votes += next.gene_votes;
}

#[derive(Debug, Default)]
struct MappingStageTimes {
    candidate_generation_seconds: f64,
    bucket_build_seconds: f64,
    scoring_seconds: f64,
    assignment_seconds: f64,
}

struct GeneEcOutputWriter {
    sink: OutputSink,
}

enum OutputSink {
    Plain(BufWriter<File>),
    Zstd {
        writer: BufWriter<ChildStdin>,
        child: Child,
    },
}

impl GeneEcOutputWriter {
    fn create(
        format: MapOutputFormatArg,
        compression: OutputCompressionArg,
        zstd_level: i32,
        path: &Path,
    ) -> Result<Self> {
        let mut sink = OutputSink::create(compression, zstd_level, path)?;
        match format {
            MapOutputFormatArg::Tsv => {
                sink.write_all(Assignment::tsv_header().as_bytes())?;
            }
            MapOutputFormatArg::GeneEcRad => {
                write_gene_ec_rad_header(&mut sink)?;
            }
        }
        Ok(Self { sink })
    }

    fn write_all(&mut self, bytes: &[u8]) -> Result<()> {
        self.sink.write_all(bytes)
    }

    fn finish(self) -> Result<()> {
        self.sink.finish()
    }
}

impl OutputSink {
    fn create(compression: OutputCompressionArg, zstd_level: i32, path: &Path) -> Result<Self> {
        match compression {
            OutputCompressionArg::None => Ok(Self::Plain(BufWriter::with_capacity(
                1 << 20,
                File::create(path)?,
            ))),
            OutputCompressionArg::Zstd => Self::create_zstd(zstd_level, path),
        }
    }

    fn create_zstd(zstd_level: i32, path: &Path) -> Result<Self> {
        if !(1..=22).contains(&zstd_level) {
            return Err(anyhow!(
                "--zstd-level must be between 1 and 22, got {zstd_level}"
            ));
        }
        let output = File::create(path)?;
        let mut child = Command::new("zstd")
            .arg("-q")
            .arg("-T0")
            .arg(format!("-{zstd_level}"))
            .arg("-c")
            .stdin(Stdio::piped())
            .stdout(Stdio::from(output))
            .spawn()
            .with_context(|| "failed to spawn zstd for compressed assignment output")?;
        let stdin = child
            .stdin
            .take()
            .ok_or_else(|| anyhow!("failed to open zstd stdin"))?;
        Ok(Self::Zstd {
            writer: BufWriter::with_capacity(1 << 20, stdin),
            child,
        })
    }

    fn write_all(&mut self, bytes: &[u8]) -> Result<()> {
        match self {
            Self::Plain(writer) => writer.write_all(bytes)?,
            Self::Zstd { writer, .. } => writer.write_all(bytes)?,
        }
        Ok(())
    }
}

impl Write for OutputSink {
    fn write(&mut self, buf: &[u8]) -> std::io::Result<usize> {
        match self {
            Self::Plain(writer) => writer.write(buf),
            Self::Zstd { writer, .. } => writer.write(buf),
        }
    }

    fn flush(&mut self) -> std::io::Result<()> {
        match self {
            Self::Plain(writer) => writer.flush(),
            Self::Zstd { writer, .. } => writer.flush(),
        }
    }
}

impl OutputSink {
    fn finish(mut self) -> Result<()> {
        match &mut self {
            Self::Plain(writer) => writer.flush()?,
            Self::Zstd { writer, .. } => writer.flush()?,
        }
        if let Self::Zstd { writer, mut child } = self {
            drop(writer);
            let status = child
                .wait()
                .with_context(|| "failed waiting for zstd assignment output")?;
            if !status.success() {
                return Err(anyhow!(
                    "zstd assignment output failed with status {status}"
                ));
            }
        }
        Ok(())
    }
}

fn output_batch_capacity(format: MapOutputFormatArg, reads: usize) -> usize {
    match format {
        MapOutputFormatArg::Tsv => reads * 64,
        MapOutputFormatArg::GeneEcRad => reads * 32,
    }
}

fn write_gene_ec_record_to_buffer(
    buffer: &mut Option<Vec<u8>>,
    format: MapOutputFormatArg,
    mapped_read: &MappedGeneEcRead,
) -> Result<()> {
    let Some(buffer) = buffer.as_mut() else {
        return Ok(());
    };
    match format {
        MapOutputFormatArg::Tsv => {
            buffer.extend_from_slice(mapped_read.assignment.to_tsv_row().as_bytes())
        }
        MapOutputFormatArg::GeneEcRad => write_gene_ec_rad_record(buffer, mapped_read)?,
    }
    Ok(())
}

fn write_gene_ec_rad_header(writer: &mut impl Write) -> Result<()> {
    writer.write_all(b"CSTRAD1\0")?;
    write_le_u32(writer, 1)?; // version
    write_le_u32(writer, 28)?; // fixed bytes per record before variable gene ids
    write_le_u32(writer, 1)?; // read IDs are implicit by record order
    write_le_u32(writer, 0)?; // reserved
    Ok(())
}

fn write_gene_ec_rad_record(buffer: &mut Vec<u8>, mapped_read: &MappedGeneEcRead) -> Result<()> {
    let assignment = &mapped_read.assignment;
    let (barcode_umi_code, barcode_umi_valid_mask) =
        pack_barcode_umi(&assignment.cell_barcode, &assignment.umi);
    let primary_gene_id = assignment.gene_id.unwrap_or(u32::MAX);
    let top_gene_count = u16::try_from(mapped_read.top_gene_ids.len())
        .map_err(|_| anyhow!("gene-EC RAD top gene set exceeds u16::MAX"))?;

    buffer.extend_from_slice(&barcode_umi_code.to_le_bytes());
    buffer.extend_from_slice(&barcode_umi_valid_mask.to_le_bytes());
    buffer.extend_from_slice(&primary_gene_id.to_le_bytes());
    buffer.extend_from_slice(&assignment.candidate_count.to_le_bytes());
    buffer.extend_from_slice(&assignment.score.to_le_bytes());
    buffer.extend_from_slice(&assignment.flags.to_le_bytes());
    buffer.push(assignment_type_code(assignment.assignment_type));
    buffer.push(0);
    buffer.extend_from_slice(&top_gene_count.to_le_bytes());
    for &gene_id in &mapped_read.top_gene_ids {
        buffer.extend_from_slice(&gene_id.to_le_bytes());
    }
    Ok(())
}

fn pack_barcode_umi(cell_barcode: &str, umi: &str) -> (u64, u32) {
    let mut code = 0_u64;
    let mut valid_mask = 0_u32;
    for (idx, base) in cell_barcode.bytes().chain(umi.bytes()).take(28).enumerate() {
        if let Some(bits) = base_2bit(base) {
            code |= (bits as u64) << (idx * 2);
            valid_mask |= 1_u32 << idx;
        }
    }
    (code, valid_mask)
}

fn base_2bit(base: u8) -> Option<u8> {
    match base.to_ascii_uppercase() {
        b'A' => Some(0),
        b'C' => Some(1),
        b'G' => Some(2),
        b'T' => Some(3),
        _ => None,
    }
}

fn assignment_type_code(assignment_type: AssignmentType) -> u8 {
    match assignment_type {
        AssignmentType::UniqueGene => 1,
        AssignmentType::AmbiguousGene => 2,
        AssignmentType::AmbiguousTranscriptSameGene => 3,
        AssignmentType::AntisenseGene => 4,
        AssignmentType::Unmapped => 5,
        AssignmentType::LowComplexity => 6,
        AssignmentType::LowQuality => 7,
    }
}

fn write_le_u32(writer: &mut impl Write, value: u32) -> Result<()> {
    writer.write_all(&value.to_le_bytes())?;
    Ok(())
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

struct RawReadPair {
    read_id: u64,
    r1_seq: Vec<u8>,
    r2_seq: Vec<u8>,
    r2_qual: Vec<u8>,
}

struct MappedGeneEcRead {
    assignment: Assignment,
    stats: EcMappingStats,
    gene_candidate_count: u32,
    top_gene_ids: Vec<u32>,
    candidate_stats: CandidateGenerationStats,
    seq: Vec<u8>,
    seq_len: usize,
}

struct PreparedGeneEcRead {
    prepared: PreparedRead,
    seq_len: usize,
    locality_key: u32,
}

struct StreamFastqRecord {
    seq: Vec<u8>,
    qual: Vec<u8>,
}

struct FastqRecordBatch {
    records: Vec<StreamFastqRecord>,
    records_read: u64,
}

type FastqBatchResult = std::result::Result<FastqRecordBatch, String>;

struct FastqPairBatchPipeline {
    r1_rx: Receiver<FastqBatchResult>,
    r2_rx: Receiver<FastqBatchResult>,
    r1_handle: Option<JoinHandle<()>>,
    r2_handle: Option<JoinHandle<()>>,
    next_read_id: u64,
    finished: bool,
}

struct FastqStreamReader {
    reader: Box<dyn BufRead>,
    records_read: u64,
    id_buf: String,
    seq_buf: String,
    plus_buf: String,
    qual_buf: String,
}

impl FastqStreamReader {
    fn open(path: &Path) -> Result<Self> {
        const BUFFER_SIZE: usize = 1 << 20;
        let file = File::open(path)?;
        let reader: Box<dyn BufRead> = if path.display().to_string().ends_with(".gz") {
            let file_reader = BufReader::with_capacity(BUFFER_SIZE, file);
            let decoder = MultiGzDecoder::new(file_reader);
            Box::new(BufReader::with_capacity(BUFFER_SIZE, decoder))
        } else {
            Box::new(BufReader::with_capacity(BUFFER_SIZE, file))
        };
        Ok(Self {
            reader,
            records_read: 0,
            id_buf: String::new(),
            seq_buf: String::new(),
            plus_buf: String::new(),
            qual_buf: String::new(),
        })
    }

    fn next_record(&mut self) -> Result<Option<StreamFastqRecord>> {
        self.id_buf.clear();
        self.seq_buf.clear();
        self.plus_buf.clear();
        self.qual_buf.clear();
        if self.reader.read_line(&mut self.id_buf)? == 0 {
            return Ok(None);
        }
        self.records_read += 1;
        if self.reader.read_line(&mut self.seq_buf)? == 0 {
            return Err(anyhow!(
                "malformed FASTQ near record {}: missing sequence line",
                self.records_read
            ));
        }
        if self.reader.read_line(&mut self.plus_buf)? == 0 {
            return Err(anyhow!(
                "malformed FASTQ near record {}: missing plus line",
                self.records_read
            ));
        }
        if self.reader.read_line(&mut self.qual_buf)? == 0 {
            return Err(anyhow!(
                "malformed FASTQ near record {}: missing quality line",
                self.records_read
            ));
        }
        let id = trim_fastq_line(self.id_buf.as_bytes());
        let seq = trim_fastq_line(self.seq_buf.as_bytes());
        let plus = trim_fastq_line(self.plus_buf.as_bytes());
        let qual = trim_fastq_line(self.qual_buf.as_bytes());
        if !id.starts_with(b"@") {
            return Err(anyhow!(
                "malformed FASTQ near record {}: identifier line must start with @",
                self.records_read
            ));
        }
        if !plus.starts_with(b"+") {
            return Err(anyhow!(
                "malformed FASTQ near record {}: separator line must start with +",
                self.records_read
            ));
        }
        if seq.len() != qual.len() {
            return Err(anyhow!(
                "malformed FASTQ near record {}: sequence length {} != quality length {}",
                self.records_read,
                seq.len(),
                qual.len()
            ));
        }
        Ok(Some(StreamFastqRecord {
            seq: seq.to_ascii_uppercase(),
            qual: qual.to_vec(),
        }))
    }
}

impl FastqPairBatchPipeline {
    fn new(r1: PathBuf, r2: PathBuf, batch_size: usize) -> Self {
        let r1_batch_size = batch_size.max(1);
        let r2_batch_size = batch_size.max(1);
        let (r1_tx, r1_rx) = mpsc::sync_channel(2);
        let (r2_tx, r2_rx) = mpsc::sync_channel(2);
        let r1_handle = thread::spawn(move || {
            stream_fastq_batches("R1", r1, r1_batch_size, r1_tx);
        });
        let r2_handle = thread::spawn(move || {
            stream_fastq_batches("R2", r2, r2_batch_size, r2_tx);
        });
        Self {
            r1_rx,
            r2_rx,
            r1_handle: Some(r1_handle),
            r2_handle: Some(r2_handle),
            next_read_id: 0,
            finished: false,
        }
    }

    fn next_batch(&mut self) -> Result<Option<(Vec<RawReadPair>, Duration)>> {
        if self.finished {
            return Ok(None);
        }
        let wait_started = Instant::now();
        let r1_batch = recv_fastq_batch(&self.r1_rx, "R1")?;
        let r2_batch = recv_fastq_batch(&self.r2_rx, "R2")?;
        let read_elapsed = wait_started.elapsed();

        if r1_batch.records.len() != r2_batch.records.len() {
            return Err(anyhow!(
                "R1/R2 record count mismatch: {} != {}",
                r1_batch.records_read,
                r2_batch.records_read
            ));
        }
        if r1_batch.records.is_empty() {
            self.finished = true;
            return Ok(None);
        }

        let mut batch = Vec::with_capacity(r1_batch.records.len());
        for (r1_record, r2_record) in r1_batch.records.into_iter().zip(r2_batch.records) {
            let read_id = self.next_read_id;
            self.next_read_id += 1;
            batch.push(RawReadPair {
                read_id,
                r1_seq: r1_record.seq,
                r2_seq: r2_record.seq,
                r2_qual: r2_record.qual,
            });
        }
        Ok(Some((batch, read_elapsed)))
    }

    fn finish(mut self) -> Result<()> {
        if let Some(handle) = self.r1_handle.take() {
            handle
                .join()
                .map_err(|_| anyhow!("R1 FASTQ reader thread panicked"))?;
        }
        if let Some(handle) = self.r2_handle.take() {
            handle
                .join()
                .map_err(|_| anyhow!("R2 FASTQ reader thread panicked"))?;
        }
        Ok(())
    }
}

fn recv_fastq_batch(rx: &Receiver<FastqBatchResult>, label: &str) -> Result<FastqRecordBatch> {
    rx.recv()
        .map_err(|_| anyhow!("{label} FASTQ reader thread stopped"))?
        .map_err(|err| anyhow!("{err}"))
}

fn stream_fastq_batches(
    label: &'static str,
    path: PathBuf,
    batch_size: usize,
    tx: mpsc::SyncSender<FastqBatchResult>,
) {
    let result = (|| -> Result<()> {
        let mut reader = FastqStreamReader::open(&path)?;
        loop {
            let records = read_fastq_record_batch(&mut reader, batch_size)?;
            let is_empty = records.is_empty();
            let records_read = reader.records_read;
            if tx
                .send(Ok(FastqRecordBatch {
                    records,
                    records_read,
                }))
                .is_err()
            {
                break;
            }
            if is_empty {
                break;
            }
        }
        Ok(())
    })();
    if let Err(err) = result {
        let _ = tx.send(Err(format!("{label} FASTQ reader failed: {err:#}")));
    }
}

fn read_fastq_record_batch(
    reader: &mut FastqStreamReader,
    batch_size: usize,
) -> Result<Vec<StreamFastqRecord>> {
    let mut records = Vec::with_capacity(batch_size.max(1));
    for _ in 0..batch_size.max(1) {
        let Some(record) = reader.next_record()? else {
            break;
        };
        records.push(record);
    }
    Ok(records)
}

fn trim_fastq_line(mut bytes: &[u8]) -> &[u8] {
    if bytes.ends_with(b"\n") {
        bytes = &bytes[..bytes.len() - 1];
    }
    if bytes.ends_with(b"\r") {
        bytes = &bytes[..bytes.len() - 1];
    }
    bytes
}

#[allow(clippy::too_many_arguments)]
fn prepare_and_map_gene_ec_raw(
    index: &LoadedEcIndex,
    raw: &RawReadPair,
    trim_tso: bool,
    min_tso_match_len: u8,
    max_tso_mismatches: u8,
    min_mean_quality: f64,
    search_reverse_complement: bool,
    sparse_probe_stride: u32,
    max_gene_df: u32,
    collect_diagnostics: bool,
    emit_top_gene_ids: bool,
) -> MappedGeneEcRead {
    let (barcode_umi, invalid_bc_umi) = parse_tenx_3p_v3_r1(&raw.r1_seq)
        .map(|parsed| (parsed, false))
        .unwrap_or_else(|_| (lossy_barcode_umi(&raw.r1_seq), true));
    let tso_trim = if trim_tso {
        tso_prefix_trim_len(&raw.r2_seq, min_tso_match_len, max_tso_mismatches)
    } else {
        0
    };
    let seq = raw.r2_seq[tso_trim..].to_vec();
    let qual = raw.r2_qual[tso_trim..].to_vec();
    let low_quality = mean_phred_quality(&qual) < min_mean_quality;
    let sketch = sketch_read(raw.read_id, &seq, index.k());
    let seq_len = seq.len();
    let seq_for_diagnostics = if collect_diagnostics {
        seq.clone()
    } else {
        Vec::new()
    };
    let prepared = PreparedRead {
        read_id: raw.read_id,
        seq,
        qual,
        cell_barcode: barcode_umi.cell_barcode_seq,
        umi: barcode_umi.umi_seq,
        invalid_barcode_umi: invalid_bc_umi,
        low_quality,
        sketch,
    };
    let (assignment, stats, gene_candidate_count, top_gene_ids) = map_gene_ec_read(
        index,
        &prepared,
        search_reverse_complement,
        min_mean_quality,
        sparse_probe_stride,
        max_gene_df,
        emit_top_gene_ids,
    );
    let candidate_stats = ec_candidate_stats(stats, gene_candidate_count);
    MappedGeneEcRead {
        assignment,
        stats,
        gene_candidate_count,
        top_gene_ids,
        candidate_stats,
        seq: seq_for_diagnostics,
        seq_len,
    }
}

#[allow(clippy::too_many_arguments)]
fn prepare_gene_ec_raw(
    index: &LoadedEcIndex,
    raw: &RawReadPair,
    trim_tso: bool,
    min_tso_match_len: u8,
    max_tso_mismatches: u8,
    min_mean_quality: f64,
    sparse_probe_stride: u32,
    sort_by_locality: bool,
) -> PreparedGeneEcRead {
    let (barcode_umi, invalid_bc_umi) = parse_tenx_3p_v3_r1(&raw.r1_seq)
        .map(|parsed| (parsed, false))
        .unwrap_or_else(|_| (lossy_barcode_umi(&raw.r1_seq), true));
    let tso_trim = if trim_tso {
        tso_prefix_trim_len(&raw.r2_seq, min_tso_match_len, max_tso_mismatches)
    } else {
        0
    };
    let seq = raw.r2_seq[tso_trim..].to_vec();
    let qual = raw.r2_qual[tso_trim..].to_vec();
    let low_quality = mean_phred_quality(&qual) < min_mean_quality;
    let sketch = sketch_read(raw.read_id, &seq, index.k());
    let locality_key = if sort_by_locality {
        gene_ec_locality_key(index, &seq, sparse_probe_stride)
    } else {
        raw.read_id.min(u32::MAX as u64) as u32
    };
    let seq_len = seq.len();
    let prepared = PreparedRead {
        read_id: raw.read_id,
        seq,
        qual,
        cell_barcode: barcode_umi.cell_barcode_seq,
        umi: barcode_umi.umi_seq,
        invalid_barcode_umi: invalid_bc_umi,
        low_quality,
        sketch,
    };
    PreparedGeneEcRead {
        prepared,
        seq_len,
        locality_key,
    }
}

fn map_prepared_gene_ec(
    index: &LoadedEcIndex,
    prepared_read: &PreparedGeneEcRead,
    search_reverse_complement: bool,
    min_mean_quality: f64,
    sparse_probe_stride: u32,
    max_gene_df: u32,
    collect_diagnostics: bool,
    emit_top_gene_ids: bool,
) -> MappedGeneEcRead {
    let (assignment, stats, gene_candidate_count, top_gene_ids) = map_gene_ec_read(
        index,
        &prepared_read.prepared,
        search_reverse_complement,
        min_mean_quality,
        sparse_probe_stride,
        max_gene_df,
        emit_top_gene_ids,
    );
    let candidate_stats = ec_candidate_stats(stats, gene_candidate_count);
    MappedGeneEcRead {
        assignment,
        stats,
        gene_candidate_count,
        top_gene_ids,
        candidate_stats,
        seq: if collect_diagnostics {
            prepared_read.prepared.seq.clone()
        } else {
            Vec::new()
        },
        seq_len: prepared_read.seq_len,
    }
}

fn gene_ec_locality_key(index: &LoadedEcIndex, seq: &[u8], sparse_probe_stride: u32) -> u32 {
    let encoded = encode_acgt(seq);
    let stride = sparse_probe_stride.max(1);
    iter_kmers_2bit(&encoded, index.k())
        .find(|kmer| kmer.pos % stride == 0)
        .map(|kmer| index.kmer_prefix(kmer.code))
        .unwrap_or(u32::MAX)
}

fn ec_candidate_stats(
    stats: EcMappingStats,
    gene_candidate_count: u32,
) -> CandidateGenerationStats {
    let mut out = CandidateGenerationStats {
        seed_lookups: stats.seed_lookups,
        seed_candidates_considered: stats.seed_lookups,
        query_seed_occurrences: stats.seed_lookups,
        selected_seeds: stats.seeds_found,
        selected_seeds_with_postings: stats.seeds_found,
        selected_seeds_absent_from_index: stats.seed_lookups.saturating_sub(stats.seeds_found),
        selected_seeds_over_frequency_cap: stats.seeds_skipped_high_gene_df,
        candidate_hits: gene_candidate_count as u64,
        candidate_votes: stats.gene_votes,
        ..CandidateGenerationStats::default()
    };
    if gene_candidate_count > 0 {
        out.reads_with_candidate_hits = 1;
    } else {
        out.reads_with_no_seed_postings = 1;
    }
    out
}

#[derive(Debug, Default)]
struct AssignmentCounts {
    unique_gene: u64,
    ambiguous_gene: u64,
    same_gene_multitranscript: u64,
    antisense_gene: u64,
    unmapped: u64,
    low_complexity: u64,
    low_quality: u64,
    exon_target: u64,
    intron_target: u64,
    gene_body_target: u64,
    mixed_target: u64,
    unknown_target: u64,
    rescued: u64,
}

impl AssignmentCounts {
    fn add(&mut self, assignment: &Assignment) {
        match assignment.assignment_type {
            AssignmentType::UniqueGene => self.unique_gene += 1,
            AssignmentType::AmbiguousGene => self.ambiguous_gene += 1,
            AssignmentType::AmbiguousTranscriptSameGene => self.same_gene_multitranscript += 1,
            AssignmentType::AntisenseGene => self.antisense_gene += 1,
            AssignmentType::Unmapped => self.unmapped += 1,
            AssignmentType::LowComplexity => self.low_complexity += 1,
            AssignmentType::LowQuality => self.low_quality += 1,
        }
        if assignment.flags & SCORE_FLAG_TARGET_EXON != 0 {
            self.exon_target += 1;
        }
        if assignment.flags & SCORE_FLAG_TARGET_INTRON != 0 {
            self.intron_target += 1;
        }
        if assignment.flags & SCORE_FLAG_TARGET_GENE_BODY != 0 {
            self.gene_body_target += 1;
        }
        if assignment.flags & SCORE_FLAG_RESCUED != 0 {
            self.rescued += 1;
        }
        let target_count = [
            SCORE_FLAG_TARGET_EXON,
            SCORE_FLAG_TARGET_INTRON,
            SCORE_FLAG_TARGET_GENE_BODY,
        ]
        .iter()
        .filter(|&&flag| assignment.flags & flag != 0)
        .count();
        if target_count > 1 {
            self.mixed_target += 1;
        }
        if is_target_assignment_type(assignment.assignment_type) && target_count == 0 {
            self.unknown_target += 1;
        }
    }

    fn rate(&self, total: u64, assignment_type: AssignmentType) -> f64 {
        let count = match assignment_type {
            AssignmentType::UniqueGene => self.unique_gene,
            AssignmentType::AmbiguousGene => self.ambiguous_gene,
            AssignmentType::AmbiguousTranscriptSameGene => self.same_gene_multitranscript,
            AssignmentType::AntisenseGene => self.antisense_gene,
            AssignmentType::Unmapped => self.unmapped,
            AssignmentType::LowComplexity => self.low_complexity,
            AssignmentType::LowQuality => self.low_quality,
        };
        count as f64 / total.max(1) as f64
    }

    fn flag_rate(&self, total: u64, flag: u16) -> f64 {
        let count = match flag {
            SCORE_FLAG_TARGET_EXON => self.exon_target,
            SCORE_FLAG_TARGET_INTRON => self.intron_target,
            SCORE_FLAG_TARGET_GENE_BODY => self.gene_body_target,
            SCORE_FLAG_RESCUED => self.rescued,
            _ => 0,
        };
        count as f64 / total.max(1) as f64
    }

    fn mixed_target_rate(&self, total: u64) -> f64 {
        self.mixed_target as f64 / total.max(1) as f64
    }

    fn unknown_target_rate(&self, total: u64) -> f64 {
        self.unknown_target as f64 / total.max(1) as f64
    }
}

fn is_target_assignment_type(assignment_type: AssignmentType) -> bool {
    matches!(
        assignment_type,
        AssignmentType::UniqueGene
            | AssignmentType::AmbiguousTranscriptSameGene
            | AssignmentType::AmbiguousGene
            | AssignmentType::AntisenseGene
    )
}

#[derive(Debug, Default)]
struct FunctionTimings {
    totals: BTreeMap<&'static str, Duration>,
}

impl FunctionTimings {
    fn add_elapsed(&mut self, name: &'static str, started: Instant) {
        *self.totals.entry(name).or_default() += started.elapsed();
    }

    fn add_duration(&mut self, name: &'static str, duration: Duration) {
        *self.totals.entry(name).or_default() += duration;
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

fn build_two_tier_gene_map(hot: &dyn IndexAccess, cold: &dyn IndexAccess) -> Result<Vec<u32>> {
    if hot.k() != cold.k() {
        return Err(anyhow!(
            "hot and cold indexes must use the same k: hot={}, cold={}",
            hot.k(),
            cold.k()
        ));
    }
    let hot_genes: HashMap<_, _> = (0..hot.num_genes())
        .filter_map(|gene_id| {
            hot.gene_name(gene_id as u32)
                .map(|name| (name.to_owned(), gene_id as u32))
        })
        .collect();
    let mut cold_to_hot = Vec::with_capacity(cold.num_genes());
    for cold_gene_id in 0..cold.num_genes() {
        let Some(name) = cold.gene_name(cold_gene_id as u32) else {
            return Err(anyhow!("cold index gene_id {} has no name", cold_gene_id));
        };
        let Some(&hot_gene_id) = hot_genes.get(name) else {
            return Err(anyhow!(
                "cold index gene {} at gene_id {} is absent from hot index",
                name,
                cold_gene_id
            ));
        };
        cold_to_hot.push(hot_gene_id);
    }
    Ok(cold_to_hot)
}

fn make_scorer(
    score_mode: ScoreMode,
    config: ScoreConfig,
    index: &dyn IndexAccess,
) -> Box<dyn CandidateScorer> {
    let target_class_flags: &'static [u16] = Box::leak(
        (0..index.num_transcripts())
            .map(|idx| target_class_flag(index.transcript_target_class(idx as u32)))
            .collect::<Vec<_>>()
            .into_boxed_slice(),
    );
    match score_mode {
        ScoreMode::Scalar => Box::new(ScalarScorer {
            config,
            target_class_flags,
        }),
        ScoreMode::Pulp => Box::new(PulpScorer {
            config,
            target_class_flags,
        }),
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
        CandidateSearchArg::ExonFirst => CandidateSearchMode::ExonFirst(SparseProbeConfig {
            stride: args.sparse_probe_stride,
            max_seeds: args.sparse_probe_max_seeds,
            min_seed_hits: args.sparse_probe_min_seed_hits,
        }),
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

#[derive(Debug, Default)]
struct ColdFallbackSummary {
    fallback_reads: u64,
    replaced_reads: u64,
    shard_groups: u64,
}

#[allow(clippy::too_many_arguments)]
fn run_cold_fallback(
    args: &MapArgs,
    reads: &[(u64, Vec<u8>)],
    quals: &[(u64, Vec<u8>)],
    barcode_umis: &[(String, String)],
    cold_index: &LoadedIndex,
    cold_scorer: &dyn CandidateScorer,
    cold_to_primary_gene: Option<&[u32]>,
    regroup_by_locus: bool,
    stats: &mut CandidateGenerationStats,
    score_stats: &mut ScoreFailureStats,
    stage_times: &mut MappingStageTimes,
    function_timings: &mut FunctionTimings,
    bucket_sizes: &mut Vec<usize>,
    read_candidate_stats: &mut [CandidateGenerationStats],
    mut read_score_stats: Option<&mut [ScoreFailureStats]>,
    assignments: &mut [Assignment],
) -> ColdFallbackSummary {
    let sort_started = Instant::now();
    let prefix_bits = args
        .cold_shard_prefix_bits
        .min(cold_index.k().saturating_mul(2));
    let mut keyed_read_ids: Vec<(u64, u64, u64)> = assignments
        .par_iter()
        .filter(|assignment| should_try_cold_fallback(assignment))
        .map(|assignment| {
            let read_id = assignment.read_id;
            if args.cold_shard_scheduling {
                let (_, seq) = &reads[read_id as usize];
                let best_seed = first_sparse_seed_key(
                    seq,
                    cold_index.k(),
                    args.sparse_probe_stride,
                    args.search_reverse_complement,
                );
                (
                    cold_shard_key(best_seed, cold_index.k(), prefix_bits),
                    best_seed,
                    read_id,
                )
            } else {
                (0, read_id, read_id)
            }
        })
        .collect();
    if args.cold_shard_scheduling {
        keyed_read_ids.par_sort_unstable();
        function_timings.add_elapsed("cold_sort_reads_by_seed_shard", sort_started);
    }

    let mut summary = ColdFallbackSummary {
        fallback_reads: keyed_read_ids.len() as u64,
        ..ColdFallbackSummary::default()
    };
    if keyed_read_ids.is_empty() {
        return summary;
    }

    let chunk_size = args.cold_batch_size.unwrap_or(args.batch_size).max(1);
    let mut group_start = 0;
    while group_start < keyed_read_ids.len() {
        let group_end = if args.cold_shard_scheduling {
            let shard = keyed_read_ids[group_start].0;
            let mut group_end = group_start + 1;
            while group_end < keyed_read_ids.len() && keyed_read_ids[group_end].0 == shard {
                group_end += 1;
            }
            group_end
        } else {
            keyed_read_ids.len()
        };
        summary.shard_groups += 1;
        for chunk in keyed_read_ids[group_start..group_end].chunks(chunk_size) {
            let chunk_read_ids: Vec<_> = chunk.iter().map(|(_, _, read_id)| *read_id).collect();
            let chunk_result = map_cold_seed_batched_chunk(
                args,
                reads,
                quals,
                barcode_umis,
                &chunk_read_ids,
                cold_index,
                cold_scorer,
                regroup_by_locus,
                read_score_stats.as_deref_mut(),
            );
            merge_cold_fallback_chunk(
                chunk_result,
                stats,
                score_stats,
                stage_times,
                function_timings,
                bucket_sizes,
                read_candidate_stats,
                assignments,
                cold_to_primary_gene,
                &mut summary,
            );
        }
        if args.cold_evict_interval > 0 && summary.shard_groups % args.cold_evict_interval == 0 {
            let advise_started = Instant::now();
            if cold_index.advise_dontneed().is_ok() {
                function_timings.add_elapsed("cold_index_advise_dontneed", advise_started);
            }
        }
        group_start = group_end;
    }
    if args.cold_evict_interval > 0 {
        let advise_started = Instant::now();
        if cold_index.advise_dontneed().is_ok() {
            function_timings.add_elapsed("cold_index_advise_dontneed", advise_started);
        }
    }
    summary
}

#[allow(clippy::too_many_arguments)]
fn map_cold_seed_batched_chunk(
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
        args.cold_max_seeds_per_read.max(1),
        args.max_postings_per_seed,
        args.search_reverse_complement,
        seed_planner(args.seed_planner),
        CandidateSearchMode::SparseProbeNoFallback(SparseProbeConfig {
            stride: args.sparse_probe_stride,
            max_seeds: args.cold_max_seeds_per_read.max(1),
            min_seed_hits: 1,
        }),
        candidate_pruning(args),
    );
    stage_times.candidate_generation_seconds += candidate_started.elapsed().as_secs_f64();
    function_timings.add_elapsed(
        "cold_generate_candidate_hits_seed_batched",
        candidate_started,
    );

    let bucket_started = Instant::now();
    let buckets = if regroup_by_locus {
        make_candidate_locus_buckets(hits, args.candidate_bin_size)
    } else {
        make_single_hit_buckets(hits, args.candidate_bin_size)
    };
    stage_times.bucket_build_seconds += bucket_started.elapsed().as_secs_f64();
    if regroup_by_locus {
        function_timings.add_elapsed("cold_make_candidate_locus_buckets", bucket_started);
    } else {
        function_timings.add_elapsed("cold_make_single_hit_buckets", bucket_started);
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
    function_timings.add_elapsed("cold_assign_streamed_chunk", assignment_started);

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

fn should_try_cold_fallback(assignment: &Assignment) -> bool {
    matches!(
        assignment.assignment_type,
        AssignmentType::Unmapped | AssignmentType::AmbiguousGene | AssignmentType::AntisenseGene
    )
}

fn cold_shard_key(kmer_code: u64, k: u8, prefix_bits: u8) -> u64 {
    if kmer_code == u64::MAX || prefix_bits == 0 {
        return kmer_code;
    }
    let total_bits = k.saturating_mul(2);
    kmer_code >> total_bits.saturating_sub(prefix_bits)
}

fn first_sparse_seed_key(seq: &[u8], k: u8, stride: u32, include_reverse: bool) -> u64 {
    let stride = stride.max(1);
    let encoded = encode_acgt(seq);
    let forward = iter_kmers_2bit(&encoded, k)
        .find(|kmer| kmer.pos % stride == 0)
        .map_or(u64::MAX, |kmer| kmer.code);
    if !include_reverse {
        return forward;
    }
    let reverse_encoded = reverse_complement(&encoded);
    let reverse = iter_kmers_2bit(&reverse_encoded, k)
        .find(|kmer| kmer.pos % stride == 0)
        .map_or(u64::MAX, |kmer| kmer.code);
    forward.min(reverse)
}

#[allow(clippy::too_many_arguments)]
fn merge_cold_fallback_chunk(
    chunk_result: StreamedChunkResult,
    stats: &mut CandidateGenerationStats,
    score_stats: &mut ScoreFailureStats,
    stage_times: &mut MappingStageTimes,
    function_timings: &mut FunctionTimings,
    bucket_sizes: &mut Vec<usize>,
    read_candidate_stats: &mut [CandidateGenerationStats],
    assignments: &mut [Assignment],
    cold_to_primary_gene: Option<&[u32]>,
    summary: &mut ColdFallbackSummary,
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
            add_stats(slot, read_stats);
        }
    }
    for (read_id, mut assignment) in chunk_result.assignments {
        if !is_rescue_replacement(&assignment) {
            continue;
        }
        if let Some(map) = cold_to_primary_gene {
            let Some(cold_gene_id) = assignment.gene_id else {
                continue;
            };
            let Some(&primary_gene_id) = map.get(cold_gene_id as usize) else {
                continue;
            };
            assignment.gene_id = Some(primary_gene_id);
            assignment.transcript_id = None;
        }
        assignment.flags |= SCORE_FLAG_RESCUED;
        if let Some(slot) = assignments.get_mut(read_id as usize) {
            *slot = assignment;
            summary.replaced_reads += 1;
        }
    }
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
    let (hits, mut stats, mut per_read_stats) = generate_candidate_hits_seed_batched(
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
    let mut bucket_sizes: Vec<_> = buckets.iter().map(|bucket| bucket.hits.len()).collect();

    let scoring_started = Instant::now();
    let (scored, mut score_stats) = score_buckets_flat(
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
    let mut assignments = make_streamed_chunk_assignments(chunk_read_ids, &scored, barcode_umis);
    stage_times.assignment_seconds += assignment_started.elapsed().as_secs_f64();
    function_timings.add_elapsed("assign_streamed_chunk", assignment_started);

    if args.rescue_failed_reads || args.rescue_score_failed_reads {
        let rescue_read_ids: Vec<_> = assignments
            .iter()
            .filter_map(|(read_id, assignment)| {
                let read_stats = per_read_stats
                    .iter()
                    .find(|(stats_read_id, _)| stats_read_id == read_id)
                    .map(|(_, stats)| *stats)
                    .unwrap_or_default();
                should_rescue_assignment(args, assignment, read_stats).then_some(*read_id)
            })
            .collect();
        if !rescue_read_ids.is_empty() {
            let rescue_result = rescue_seed_batched_chunk(
                args,
                reads,
                quals,
                barcode_umis,
                &rescue_read_ids,
                index,
                scorer,
                regroup_by_locus,
                None,
            );
            add_stats(&mut stats, rescue_result.stats);
            add_score_stats(&mut score_stats, rescue_result.score_stats);
            stage_times.candidate_generation_seconds +=
                rescue_result.stage_times.candidate_generation_seconds;
            stage_times.bucket_build_seconds += rescue_result.stage_times.bucket_build_seconds;
            stage_times.scoring_seconds += rescue_result.stage_times.scoring_seconds;
            stage_times.assignment_seconds += rescue_result.stage_times.assignment_seconds;
            function_timings.add_all(rescue_result.function_timings);
            bucket_sizes.extend(rescue_result.bucket_sizes);
            for (read_id, read_stats) in rescue_result.per_read_stats {
                if let Some((_, slot)) = per_read_stats
                    .iter_mut()
                    .find(|(slot_read_id, _)| *slot_read_id == read_id)
                {
                    add_stats(slot, read_stats);
                }
            }
            for (read_id, mut rescued) in rescue_result.assignments {
                if !is_rescue_replacement(&rescued) {
                    continue;
                }
                rescued.flags |= SCORE_FLAG_RESCUED;
                if let Some((_, assignment)) = assignments
                    .iter_mut()
                    .find(|(assignment_read_id, _)| *assignment_read_id == read_id)
                {
                    *assignment = rescued;
                }
            }
        }
    }

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

#[allow(clippy::too_many_arguments)]
fn rescue_seed_batched_chunk(
    args: &MapArgs,
    reads: &[(u64, Vec<u8>)],
    quals: &[(u64, Vec<u8>)],
    barcode_umis: &[(String, String)],
    rescue_read_ids: &[u64],
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
        rescue_read_ids,
        args.min_seed_quality,
        args.max_seeds_per_read,
        args.max_postings_per_seed,
        args.search_reverse_complement || args.rescue_reverse_complement,
        seed_planner(args.seed_planner),
        CandidateSearchMode::Full,
        CandidatePruningMode::None,
    );
    stage_times.candidate_generation_seconds += candidate_started.elapsed().as_secs_f64();
    function_timings.add_elapsed(
        "rescue_generate_candidate_hits_seed_batched",
        candidate_started,
    );

    let bucket_started = Instant::now();
    let buckets = if regroup_by_locus {
        make_candidate_locus_buckets(hits, args.candidate_bin_size)
    } else {
        make_single_hit_buckets(hits, args.candidate_bin_size)
    };
    stage_times.bucket_build_seconds += bucket_started.elapsed().as_secs_f64();
    if regroup_by_locus {
        function_timings.add_elapsed("rescue_make_candidate_locus_buckets", bucket_started);
    } else {
        function_timings.add_elapsed("rescue_make_single_hit_buckets", bucket_started);
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
    let assignments = make_streamed_chunk_assignments(rescue_read_ids, &scored, barcode_umis);
    stage_times.assignment_seconds += assignment_started.elapsed().as_secs_f64();
    function_timings.add_elapsed("rescue_assign_streamed_chunk", assignment_started);

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

fn should_rescue_assignment(
    args: &MapArgs,
    assignment: &Assignment,
    read_stats: CandidateGenerationStats,
) -> bool {
    (args.rescue_failed_reads && matches!(assignment.assignment_type, AssignmentType::Unmapped))
        || (args.rescue_score_failed_reads
            && matches!(assignment.assignment_type, AssignmentType::Unmapped)
            && assignment.candidate_count == 0
            && read_stats.candidate_hits > 0)
        || (args.rescue_antisense_reads
            && matches!(assignment.assignment_type, AssignmentType::AntisenseGene))
}

fn is_rescue_replacement(assignment: &Assignment) -> bool {
    matches!(
        assignment.assignment_type,
        AssignmentType::UniqueGene | AssignmentType::AmbiguousTranscriptSameGene
    )
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

fn per_read_u64(count: u64, reads: u64) -> f64 {
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

fn target_assignment_rate(assignments: &[Assignment], flag: u16) -> f64 {
    assignments
        .iter()
        .filter(|assignment| assignment.flags & flag != 0)
        .count() as f64
        / assignments.len().max(1) as f64
}

fn mixed_target_assignment_rate(assignments: &[Assignment]) -> f64 {
    assignments
        .iter()
        .filter(|assignment| {
            let count = [
                SCORE_FLAG_TARGET_EXON,
                SCORE_FLAG_TARGET_INTRON,
                SCORE_FLAG_TARGET_GENE_BODY,
            ]
            .iter()
            .filter(|&&flag| assignment.flags & flag != 0)
            .count();
            count > 1
        })
        .count() as f64
        / assignments.len().max(1) as f64
}

fn unknown_target_assignment_rate(assignments: &[Assignment]) -> f64 {
    assignments
        .iter()
        .filter(|assignment| {
            matches!(
                assignment.assignment_type,
                AssignmentType::UniqueGene
                    | AssignmentType::AmbiguousTranscriptSameGene
                    | AssignmentType::AmbiguousGene
                    | AssignmentType::AntisenseGene
            ) && assignment.flags
                & (SCORE_FLAG_TARGET_EXON | SCORE_FLAG_TARGET_INTRON | SCORE_FLAG_TARGET_GENE_BODY)
                == 0
        })
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
        "read_id\treason\tseq_len\tnum_valid_kmers\tselected_seeds\tseed_hits\tcandidate_hits\tscored_candidates\tscore_seen\tscore_oob\tscore_mismatch\tscore_min_len\tscore_trimmed_oob\tscore_full_length_mismatch_only\tsparse_probe_attempted\tsparse_probe_accepted\tsparse_probe_fallback\tflags\n",
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
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n",
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
            stats.sparse_probe_attempted_reads,
            stats.sparse_probe_accepted_reads,
            stats.sparse_probe_fallback_reads,
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

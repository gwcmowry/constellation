use crate::dna::{encode_acgt, iter_kmers_2bit, reverse_complement};
use crate::index::IndexAccess;
use crate::{GeneId, ReadId, TranscriptId};
use rayon::prelude::*;
use rustc_hash::FxHashMap;
use smallvec::SmallVec;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct CandidateHit {
    pub read_id: ReadId,
    pub transcript_id: TranscriptId,
    pub gene_id: GeneId,
    pub pos: u32,
    pub strand: u8,
    pub seed_count: u16,
    pub seed_score: u16,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
pub struct CandidateBucketKey {
    pub transcript_id: TranscriptId,
    pub pos_bin: u32,
    pub strand: u8,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct CandidateLocusBucket {
    pub key: CandidateBucketKey,
    pub hits: Vec<CandidateHit>,
}

#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub enum SeedPlanner {
    #[default]
    RawFrequency,
    GeneIdf,
}

#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub struct CandidateGenerationStats {
    pub seed_lookups: u64,
    pub candidate_hits: u64,
    pub postings_skipped_due_to_frequency: u64,
    pub seeds_skipped_due_to_quality: u64,
    pub seed_candidates_considered: u64,
    pub early_stopped_reads: u64,
    pub seed_lookups_saved_by_early_stop: u64,
    pub selected_seeds: u64,
    pub reads_with_no_valid_kmers: u64,
    pub reads_with_no_selected_seeds: u64,
    pub reads_with_no_seed_postings: u64,
    pub reads_all_selected_seeds_over_frequency_cap: u64,
    pub reads_with_candidate_hits: u64,
    pub selected_seeds_absent_from_index: u64,
    pub selected_seeds_over_frequency_cap: u64,
    pub selected_seeds_with_postings: u64,
    pub query_seed_occurrences: u64,
    pub distinct_seed_lists_loaded: u64,
    pub candidate_votes: u64,
    pub seed_groups_skipped_due_to_frequency: u64,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct QuerySeedChoice {
    pub seed_pos: u32,
    pub kmer_code: u64,
    pub read_strand: u8,
    pub raw_postings: u32,
    pub transcript_df: u32,
    pub gene_df: u32,
    pub seed_score: u16,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
pub struct QuerySeed {
    pub kmer_code: u64,
    pub read_id: ReadId,
    pub seed_pos: u32,
    pub read_strand: u8,
    pub raw_postings: u32,
    pub seed_score: u16,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
pub struct CandidateVote {
    pub read_id: ReadId,
    pub transcript_id: TranscriptId,
    pub candidate_start: u32,
    pub strand: u8,
    pub seed_score: u16,
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct CandidateEarlyStopConfig {
    pub posterior_threshold: f64,
    pub prior_alpha: f64,
    pub min_seed_lookups: usize,
    pub min_top_seed_count: u16,
}

impl CandidateEarlyStopConfig {
    pub fn is_enabled(self) -> bool {
        self.posterior_threshold.is_finite()
            && self.posterior_threshold > 0.0
            && self.posterior_threshold < 1.0
            && self.prior_alpha.is_finite()
            && self.prior_alpha > 0.0
    }
}

pub fn generate_candidate_hits(
    index: &dyn IndexAccess,
    read_id: ReadId,
    seq: &[u8],
    max_seeds_per_read: usize,
    max_postings_per_seed: usize,
) -> Vec<CandidateHit> {
    generate_candidate_hits_with_stats(
        index,
        read_id,
        seq,
        max_seeds_per_read,
        max_postings_per_seed,
    )
    .0
}

pub fn generate_candidate_hits_with_stats(
    index: &dyn IndexAccess,
    read_id: ReadId,
    seq: &[u8],
    max_seeds_per_read: usize,
    max_postings_per_seed: usize,
) -> (Vec<CandidateHit>, CandidateGenerationStats) {
    generate_candidate_hits_with_quality_stats(
        index,
        read_id,
        seq,
        None,
        0,
        max_seeds_per_read,
        max_postings_per_seed,
        false,
        None,
        SeedPlanner::RawFrequency,
    )
}

#[allow(clippy::too_many_arguments)]
pub fn generate_candidate_hits_with_quality_stats(
    index: &dyn IndexAccess,
    read_id: ReadId,
    seq: &[u8],
    qual: Option<&[u8]>,
    min_seed_quality: u8,
    max_seeds_per_read: usize,
    max_postings_per_seed: usize,
    search_reverse_complement: bool,
    early_stop: Option<CandidateEarlyStopConfig>,
    seed_planner: SeedPlanner,
) -> (Vec<CandidateHit>, CandidateGenerationStats) {
    let encoded = encode_acgt(seq);
    let seed_limit = max_seeds_per_read.max(1);
    let mut stats = CandidateGenerationStats::default();
    let mut seed_choices: SmallVec<[QuerySeedChoice; 256]> = SmallVec::new();

    collect_seed_choices(
        index,
        &encoded,
        qual,
        min_seed_quality,
        0,
        seed_planner,
        &mut stats,
        &mut seed_choices,
    );
    if search_reverse_complement {
        let reverse_encoded = reverse_complement(&encoded);
        collect_seed_choices(
            index,
            &reverse_encoded,
            None,
            min_seed_quality,
            1,
            seed_planner,
            &mut stats,
            &mut seed_choices,
        );
    }

    sort_seed_choices(&mut seed_choices, seed_planner);

    let mut grouped: FxHashMap<u64, (u16, u16)> =
        FxHashMap::with_capacity_and_hasher(seed_limit.saturating_mul(4), Default::default());

    let selected_seed_count = seed_choices.len().min(seed_limit);
    stats.selected_seeds = selected_seed_count as u64;
    if stats.seed_candidates_considered == 0 {
        stats.reads_with_no_valid_kmers = 1;
    }
    if selected_seed_count == 0 {
        stats.reads_with_no_selected_seeds = 1;
    }
    let mut selected_over_frequency_cap = 0_u64;
    for (selected_idx, seed) in seed_choices.into_iter().take(seed_limit).enumerate() {
        stats.seed_lookups += 1;
        if seed.raw_postings == 0 {
            stats.selected_seeds_absent_from_index += 1;
            continue;
        }
        if seed.raw_postings as usize > max_postings_per_seed {
            stats.postings_skipped_due_to_frequency += seed.raw_postings as u64;
            stats.selected_seeds_over_frequency_cap += 1;
            selected_over_frequency_cap += 1;
            continue;
        }
        stats.selected_seeds_with_postings += 1;
        for posting in index.seed_postings(seed.kmer_code) {
            if posting.pos < seed.seed_pos {
                continue;
            }
            let candidate_start = posting.pos - seed.seed_pos;
            let key = pack_candidate_key(posting.transcript_id, candidate_start, seed.read_strand);
            let entry = grouped.entry(key).or_default();
            entry.0 = entry.0.saturating_add(1);
            entry.1 = entry.1.saturating_add(seed.seed_score);
        }
        if should_stop_early(&grouped, selected_idx + 1, selected_seed_count, early_stop) {
            stats.early_stopped_reads = 1;
            stats.seed_lookups_saved_by_early_stop =
                (selected_seed_count - selected_idx - 1) as u64;
            break;
        }
    }

    let hits: Vec<_> = grouped
        .into_iter()
        .map(|(key, (seed_count, seed_score))| {
            let (transcript_id, pos, strand) = unpack_candidate_key(key);
            CandidateHit {
                read_id,
                transcript_id,
                gene_id: index.transcript_gene_id(transcript_id).unwrap_or(0),
                pos,
                strand,
                seed_count,
                seed_score,
            }
        })
        .collect();
    stats.candidate_hits = hits.len() as u64;
    if stats.candidate_hits > 0 {
        stats.reads_with_candidate_hits = 1;
    } else if selected_seed_count > 0 && stats.selected_seeds_with_postings == 0 {
        if selected_over_frequency_cap == selected_seed_count as u64 {
            stats.reads_all_selected_seeds_over_frequency_cap = 1;
        } else {
            stats.reads_with_no_seed_postings = 1;
        }
    }
    (hits, stats)
}

#[allow(clippy::too_many_arguments)]
pub fn generate_candidate_hits_seed_batched(
    index: &dyn IndexAccess,
    reads: &[(ReadId, Vec<u8>)],
    quals: &[(ReadId, Vec<u8>)],
    read_ids: &[ReadId],
    min_seed_quality: u8,
    max_seeds_per_read: usize,
    max_postings_per_seed: usize,
    search_reverse_complement: bool,
    seed_planner: SeedPlanner,
) -> (
    Vec<CandidateHit>,
    CandidateGenerationStats,
    Vec<(ReadId, CandidateGenerationStats)>,
) {
    let seed_limit = max_seeds_per_read.max(1);
    let mut total_stats = CandidateGenerationStats::default();
    let mut read_stats = Vec::with_capacity(read_ids.len());
    let mut query_seeds = Vec::new();

    for &read_id in read_ids {
        let Some((_, seq)) = reads.get(read_id as usize) else {
            continue;
        };
        let qual = quals.get(read_id as usize).map(|(_, qual)| qual.as_slice());
        let (selected, mut stats) = select_seed_choices_for_read(
            index,
            seq,
            qual,
            min_seed_quality,
            seed_limit,
            search_reverse_complement,
            seed_planner,
        );
        stats.query_seed_occurrences = selected.len() as u64;
        for seed in selected {
            query_seeds.push(QuerySeed {
                kmer_code: seed.kmer_code,
                read_id,
                seed_pos: seed.seed_pos,
                read_strand: seed.read_strand,
                raw_postings: seed.raw_postings,
                seed_score: seed.seed_score,
            });
        }
        add_generation_stats(&mut total_stats, stats);
        read_stats.push((read_id, stats));
    }
    let read_stat_offsets: FxHashMap<ReadId, usize> = read_stats
        .iter()
        .enumerate()
        .map(|(idx, (read_id, _))| (*read_id, idx))
        .collect();

    query_seeds.sort_unstable();
    let mut votes = Vec::new();
    let mut group_start = 0;
    while group_start < query_seeds.len() {
        let kmer_code = query_seeds[group_start].kmer_code;
        let mut group_end = group_start + 1;
        while group_end < query_seeds.len() && query_seeds[group_end].kmer_code == kmer_code {
            group_end += 1;
        }
        let postings_len = query_seeds[group_start].raw_postings as usize;
        if postings_len == 0 {
            update_seed_stats_for_group(
                &query_seeds[group_start..group_end],
                &mut read_stats,
                &read_stat_offsets,
                |stats| stats.selected_seeds_absent_from_index += 1,
            );
            group_start = group_end;
            continue;
        }
        if postings_len > max_postings_per_seed {
            total_stats.postings_skipped_due_to_frequency +=
                postings_len as u64 * (group_end - group_start) as u64;
            total_stats.seed_groups_skipped_due_to_frequency += 1;
            update_seed_stats_for_group(
                &query_seeds[group_start..group_end],
                &mut read_stats,
                &read_stat_offsets,
                |stats| {
                    stats.selected_seeds_over_frequency_cap += 1;
                    stats.postings_skipped_due_to_frequency += postings_len as u64;
                },
            );
            group_start = group_end;
            continue;
        }
        total_stats.distinct_seed_lists_loaded += 1;
        update_seed_stats_for_group(
            &query_seeds[group_start..group_end],
            &mut read_stats,
            &read_stat_offsets,
            |stats| stats.selected_seeds_with_postings += 1,
        );
        for posting in index.seed_postings(kmer_code) {
            for seed in &query_seeds[group_start..group_end] {
                if posting.pos < seed.seed_pos {
                    continue;
                }
                votes.push(CandidateVote {
                    read_id: seed.read_id,
                    transcript_id: posting.transcript_id,
                    candidate_start: posting.pos - seed.seed_pos,
                    strand: seed.read_strand,
                    seed_score: seed.seed_score,
                });
            }
        }
        group_start = group_end;
    }
    total_stats.candidate_votes = votes.len() as u64;

    votes.sort_unstable();
    let mut hits = Vec::new();
    let mut idx = 0;
    while idx < votes.len() {
        let first = votes[idx];
        let mut seed_count = 1_u16;
        let mut seed_score = first.seed_score;
        idx += 1;
        while idx < votes.len()
            && votes[idx].read_id == first.read_id
            && votes[idx].transcript_id == first.transcript_id
            && votes[idx].candidate_start == first.candidate_start
            && votes[idx].strand == first.strand
        {
            seed_count = seed_count.saturating_add(1);
            seed_score = seed_score.saturating_add(votes[idx].seed_score);
            idx += 1;
        }
        hits.push(CandidateHit {
            read_id: first.read_id,
            transcript_id: first.transcript_id,
            gene_id: index.transcript_gene_id(first.transcript_id).unwrap_or(0),
            pos: first.candidate_start,
            strand: first.strand,
            seed_count,
            seed_score,
        });
    }
    total_stats.candidate_hits = hits.len() as u64;

    for (read_id, stats) in &mut read_stats {
        let start = hits.partition_point(|hit| hit.read_id < *read_id);
        let end = start + hits[start..].partition_point(|hit| hit.read_id == *read_id);
        stats.candidate_hits = (end - start) as u64;
        if stats.candidate_hits > 0 {
            stats.reads_with_candidate_hits = 1;
        } else if stats.selected_seeds > 0 && stats.selected_seeds_with_postings == 0 {
            if stats.selected_seeds_over_frequency_cap == stats.selected_seeds {
                stats.reads_all_selected_seeds_over_frequency_cap = 1;
            } else {
                stats.reads_with_no_seed_postings = 1;
            }
        }
    }
    for (_, stats) in &read_stats {
        total_stats.reads_with_candidate_hits += stats.reads_with_candidate_hits;
        total_stats.reads_with_no_seed_postings += stats.reads_with_no_seed_postings;
        total_stats.reads_all_selected_seeds_over_frequency_cap +=
            stats.reads_all_selected_seeds_over_frequency_cap;
        total_stats.selected_seeds_absent_from_index += stats.selected_seeds_absent_from_index;
        total_stats.selected_seeds_over_frequency_cap += stats.selected_seeds_over_frequency_cap;
        total_stats.selected_seeds_with_postings += stats.selected_seeds_with_postings;
    }

    (hits, total_stats, read_stats)
}

fn collect_seed_choices(
    index: &dyn IndexAccess,
    encoded: &crate::dna::EncodedSeq,
    qual: Option<&[u8]>,
    min_seed_quality: u8,
    read_strand: u8,
    seed_planner: SeedPlanner,
    stats: &mut CandidateGenerationStats,
    seed_choices: &mut SmallVec<[QuerySeedChoice; 256]>,
) {
    for seed in iter_kmers_2bit(encoded, index.k()) {
        if let Some(qual) = qual {
            let start = seed.pos as usize;
            let end = start + index.k() as usize;
            if end <= qual.len()
                && qual[start..end]
                    .iter()
                    .any(|&q| q.saturating_sub(33) < min_seed_quality)
            {
                stats.seeds_skipped_due_to_quality += 1;
                continue;
            }
        }
        stats.seed_candidates_considered += 1;
        let stats_for_seed = index.kmer_stats(seed.code);
        let (raw_postings, transcript_df, gene_df) = stats_for_seed
            .map(|stats| (stats.raw_postings, stats.transcript_df, stats.gene_df))
            .unwrap_or((0, 0, 0));
        seed_choices.push(QuerySeedChoice {
            seed_pos: seed.pos,
            kmer_code: seed.code,
            read_strand,
            raw_postings,
            transcript_df,
            gene_df,
            seed_score: seed_score(index.num_genes(), gene_df, seed_planner),
        });
    }
}

fn select_seed_choices_for_read(
    index: &dyn IndexAccess,
    seq: &[u8],
    qual: Option<&[u8]>,
    min_seed_quality: u8,
    seed_limit: usize,
    search_reverse_complement: bool,
    seed_planner: SeedPlanner,
) -> (SmallVec<[QuerySeedChoice; 256]>, CandidateGenerationStats) {
    let encoded = encode_acgt(seq);
    let mut stats = CandidateGenerationStats::default();
    let mut seed_choices = SmallVec::new();
    collect_seed_choices(
        index,
        &encoded,
        qual,
        min_seed_quality,
        0,
        seed_planner,
        &mut stats,
        &mut seed_choices,
    );
    if search_reverse_complement {
        let reverse_encoded = reverse_complement(&encoded);
        collect_seed_choices(
            index,
            &reverse_encoded,
            None,
            min_seed_quality,
            1,
            seed_planner,
            &mut stats,
            &mut seed_choices,
        );
    }
    sort_seed_choices(&mut seed_choices, seed_planner);
    let selected_seed_count = seed_choices.len().min(seed_limit);
    stats.selected_seeds = selected_seed_count as u64;
    stats.seed_lookups = selected_seed_count as u64;
    if stats.seed_candidates_considered == 0 {
        stats.reads_with_no_valid_kmers = 1;
    }
    if selected_seed_count == 0 {
        stats.reads_with_no_selected_seeds = 1;
    }
    seed_choices.truncate(selected_seed_count);
    (seed_choices, stats)
}

fn update_seed_stats_for_group(
    group: &[QuerySeed],
    read_stats: &mut [(ReadId, CandidateGenerationStats)],
    read_stat_offsets: &FxHashMap<ReadId, usize>,
    mut update: impl FnMut(&mut CandidateGenerationStats),
) {
    for seed in group {
        if let Some(&idx) = read_stat_offsets.get(&seed.read_id) {
            let (_, stats) = &mut read_stats[idx];
            update(stats);
        }
    }
}

fn add_generation_stats(total: &mut CandidateGenerationStats, next: CandidateGenerationStats) {
    total.seed_lookups += next.seed_lookups;
    total.seed_candidates_considered += next.seed_candidates_considered;
    total.seeds_skipped_due_to_quality += next.seeds_skipped_due_to_quality;
    total.selected_seeds += next.selected_seeds;
    total.reads_with_no_valid_kmers += next.reads_with_no_valid_kmers;
    total.reads_with_no_selected_seeds += next.reads_with_no_selected_seeds;
    total.query_seed_occurrences += next.query_seed_occurrences;
}

fn sort_seed_choices(seed_choices: &mut [QuerySeedChoice], seed_planner: SeedPlanner) {
    match seed_planner {
        SeedPlanner::RawFrequency => seed_choices.sort_unstable_by_key(|seed| {
            (
                seed.raw_postings == 0,
                nonzero_or_max(seed.raw_postings),
                seed.seed_pos,
                seed.kmer_code,
                seed.read_strand,
            )
        }),
        SeedPlanner::GeneIdf => seed_choices.sort_unstable_by_key(|seed| {
            (
                seed.raw_postings == 0,
                nonzero_or_max(seed.gene_df),
                std::cmp::Reverse(seed.seed_score),
                nonzero_or_max(seed.raw_postings),
                seed.seed_pos,
                seed.kmer_code,
                seed.read_strand,
            )
        }),
    }
}

fn nonzero_or_max(value: u32) -> u32 {
    if value == 0 {
        u32::MAX
    } else {
        value
    }
}

fn seed_score(num_genes: usize, gene_df: u32, seed_planner: SeedPlanner) -> u16 {
    match seed_planner {
        SeedPlanner::RawFrequency => 1,
        SeedPlanner::GeneIdf => {
            if gene_df == 0 {
                return 0;
            }
            let idf = ((num_genes as f64 + 1.0) / (gene_df as f64 + 1.0)).ln();
            (1.0 + idf * 64.0).round().clamp(1.0, u16::MAX as f64) as u16
        }
    }
}

fn should_stop_early(
    grouped: &FxHashMap<u64, (u16, u16)>,
    seed_lookups: usize,
    selected_seed_count: usize,
    early_stop: Option<CandidateEarlyStopConfig>,
) -> bool {
    let Some(config) = early_stop.filter(|config| config.is_enabled()) else {
        return false;
    };
    if seed_lookups < config.min_seed_lookups || grouped.len() < 2 {
        return false;
    }
    let mut top = 0_u16;
    let mut second = 0_u16;
    for &(count, _) in grouped.values() {
        if count > top {
            second = top;
            top = count;
        } else if count > second {
            second = count;
        }
    }
    if top < config.min_top_seed_count {
        return false;
    }
    let posterior_top_vs_runner_up =
        (top as f64 + config.prior_alpha) / (top as f64 + second as f64 + 2.0 * config.prior_alpha);
    let remaining = selected_seed_count.saturating_sub(seed_lookups) as u16;
    posterior_top_vs_runner_up >= config.posterior_threshold
        && top > second.saturating_add(remaining)
}

#[inline]
fn pack_candidate_key(transcript_id: TranscriptId, pos: u32, strand: u8) -> u64 {
    ((transcript_id as u64) << 32) | ((pos as u64) << 1) | u64::from(strand & 1)
}

#[inline]
fn unpack_candidate_key(key: u64) -> (TranscriptId, u32, u8) {
    (
        (key >> 32) as TranscriptId,
        ((key >> 1) & 0x7fff_ffff) as u32,
        (key & 1) as u8,
    )
}

pub fn sort_candidate_hits(hits: &mut [CandidateHit], candidate_bin_size: u32) {
    let bin_size = candidate_bin_size.max(1);
    hits.par_sort_by_key(|hit| {
        (
            hit.transcript_id,
            hit.pos / bin_size,
            hit.strand,
            hit.read_id,
            hit.pos,
        )
    });
}

pub fn make_candidate_locus_buckets(
    mut hits: Vec<CandidateHit>,
    candidate_bin_size: u32,
) -> Vec<CandidateLocusBucket> {
    let bin_size = candidate_bin_size.max(1);
    sort_candidate_hits(&mut hits, bin_size);
    let mut buckets: Vec<CandidateLocusBucket> = Vec::new();

    for hit in hits {
        let key = CandidateBucketKey {
            transcript_id: hit.transcript_id,
            pos_bin: hit.pos / bin_size,
            strand: hit.strand,
        };
        match buckets.last_mut() {
            Some(bucket) if bucket.key == key => bucket.hits.push(hit),
            _ => buckets.push(CandidateLocusBucket {
                key,
                hits: vec![hit],
            }),
        }
    }

    buckets
}

pub fn make_single_hit_buckets(
    hits: Vec<CandidateHit>,
    candidate_bin_size: u32,
) -> Vec<CandidateLocusBucket> {
    let bin_size = candidate_bin_size.max(1);
    hits.into_iter()
        .map(|hit| CandidateLocusBucket {
            key: CandidateBucketKey {
                transcript_id: hit.transcript_id,
                pos_bin: hit.pos / bin_size,
                strand: hit.strand,
            },
            hits: vec![hit],
        })
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::index_build::build_transcript_index;
    use std::io::Write;

    #[test]
    fn candidate_locus_sorting_is_stable() {
        let hits = vec![
            CandidateHit {
                read_id: 2,
                transcript_id: 1,
                gene_id: 1,
                pos: 130,
                strand: 0,
                seed_count: 1,
                seed_score: 1,
            },
            CandidateHit {
                read_id: 1,
                transcript_id: 1,
                gene_id: 1,
                pos: 64,
                strand: 0,
                seed_count: 1,
                seed_score: 1,
            },
        ];
        let buckets = make_candidate_locus_buckets(hits, 64);
        assert_eq!(buckets.len(), 2);
        assert_eq!(buckets[0].key.pos_bin, 1);
        assert_eq!(buckets[1].key.pos_bin, 2);
    }

    #[test]
    fn gene_idf_planner_prefers_low_gene_df_seed() {
        let mut fasta = tempfile::NamedTempFile::new().unwrap();
        writeln!(
            fasta,
            ">tx0|gene=GENE_A\nAAAGGG\n>tx1|gene=GENE_A\nAAATTT\n>tx2|gene=GENE_A\nAAAGAA\n>tx3|gene=GENE_B\nGGGCCC\n>tx4|gene=GENE_C\nTTTCCC"
        )
        .unwrap();
        let index = build_transcript_index(fasta.path(), 3, 256).unwrap();

        let (raw_hits, _) = generate_candidate_hits_with_quality_stats(
            &index,
            0,
            b"AAACCC",
            None,
            0,
            1,
            256,
            false,
            None,
            SeedPlanner::RawFrequency,
        );
        let (idf_hits, _) = generate_candidate_hits_with_quality_stats(
            &index,
            0,
            b"AAACCC",
            None,
            0,
            1,
            256,
            false,
            None,
            SeedPlanner::GeneIdf,
        );

        assert!(raw_hits.iter().all(|hit| hit.gene_id != 0));
        assert!(idf_hits.iter().all(|hit| hit.gene_id == 0));
        assert!(idf_hits.iter().all(|hit| hit.seed_score > hit.seed_count));
    }

    #[test]
    fn seed_batched_matches_per_read_exact_reads() {
        let mut fasta = tempfile::NamedTempFile::new().unwrap();
        writeln!(
            fasta,
            ">tx0|gene=GENE_A\nAAACCCGGG\n>tx1|gene=GENE_B\nTTTAAACCC"
        )
        .unwrap();
        let index = build_transcript_index(fasta.path(), 3, 256).unwrap();
        let reads = vec![(0, b"AAACCC".to_vec()), (1, b"CCCGGG".to_vec())];
        let quals = vec![(0, b"IIIIII".to_vec()), (1, b"IIIIII".to_vec())];
        let mut per_read = Vec::new();
        for (read_id, seq) in &reads {
            per_read.extend(
                generate_candidate_hits_with_quality_stats(
                    &index,
                    *read_id,
                    seq,
                    Some(&quals[*read_id as usize].1),
                    0,
                    4,
                    256,
                    false,
                    None,
                    SeedPlanner::RawFrequency,
                )
                .0,
            );
        }
        let (mut batched, stats, _) = generate_candidate_hits_seed_batched(
            &index,
            &reads,
            &quals,
            &[0, 1],
            0,
            4,
            256,
            false,
            SeedPlanner::RawFrequency,
        );
        sort_hits_for_test(&mut per_read);
        sort_hits_for_test(&mut batched);
        assert_eq!(batched, per_read);
        assert!(stats.distinct_seed_lists_loaded < stats.query_seed_occurrences);
        assert!(stats.candidate_votes >= stats.candidate_hits);
    }

    fn sort_hits_for_test(hits: &mut [CandidateHit]) {
        hits.sort_unstable_by_key(|hit| {
            (
                hit.read_id,
                hit.transcript_id,
                hit.pos,
                hit.strand,
                hit.seed_count,
                hit.seed_score,
            )
        });
    }
}

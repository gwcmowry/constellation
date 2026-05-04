use serde::{Deserialize, Serialize};
use std::collections::BTreeMap;
use std::time::Duration;

#[derive(Debug, Clone, Default, PartialEq, Serialize, Deserialize)]
pub struct MapMetrics {
    pub mode: String,
    pub score_mode: String,
    pub num_reads: u64,
    pub wall_seconds: f64,
    pub index_load_seconds: f64,
    pub fastq_load_seconds: f64,
    pub preprocess_seconds: f64,
    pub mapping_seconds: f64,
    pub candidate_generation_seconds: f64,
    pub bucket_build_seconds: f64,
    pub scoring_seconds: f64,
    pub assignment_seconds: f64,
    pub write_seconds: f64,
    pub reads_per_second: f64,
    pub bases_per_second: f64,
    pub seed_lookups_per_read: f64,
    pub seed_candidates_considered_per_read: f64,
    pub candidate_hits_per_read: f64,
    pub candidate_buckets_per_read: f64,
    pub postings_skipped_due_to_frequency: u64,
    pub seeds_skipped_due_to_quality: u64,
    pub early_stopped_reads: u64,
    pub seed_lookups_saved_by_early_stop: u64,
    pub scored_candidates_per_read: f64,
    pub mean_candidate_bucket_size: f64,
    pub median_candidate_bucket_size: f64,
    pub max_candidate_bucket_size: u64,
    pub reads_with_no_valid_kmers: u64,
    pub reads_with_no_selected_seeds: u64,
    pub reads_with_no_seed_postings: u64,
    pub reads_all_selected_seeds_over_frequency_cap: u64,
    pub reads_with_candidate_hits: u64,
    pub selected_seeds_absent_from_index: u64,
    pub selected_seeds_over_frequency_cap: u64,
    pub selected_seeds_with_postings: u64,
    pub query_seed_occurrences_per_read: f64,
    pub distinct_seed_lists_loaded: u64,
    pub posting_list_reuse_factor: f64,
    pub candidate_votes_per_read: f64,
    pub seed_groups_skipped_due_to_frequency: u64,
    pub sparse_probe_attempted_reads: u64,
    pub sparse_probe_accepted_reads: u64,
    pub sparse_probe_fallback_reads: u64,
    pub candidate_hits_pruned: u64,
    pub score_candidates_seen: u64,
    pub score_candidates_out_of_bounds: u64,
    pub score_candidates_failed_mismatch: u64,
    pub score_candidates_failed_min_scored_len: u64,
    pub score_candidates_failed_trimmed_out_of_bounds: u64,
    pub score_candidates_full_length_mismatch_only: u64,
    pub score_candidates_passed_full_length: u64,
    pub score_candidates_passed_trimmed_or_softclipped: u64,
    pub unique_gene_rate: f64,
    pub ambiguous_gene_rate: f64,
    pub same_gene_multitranscript_rate: f64,
    pub multi_gene_ambiguous_rate: f64,
    pub antisense_gene_rate: f64,
    pub gene_countable_rate: f64,
    pub low_complexity_rate: f64,
    pub low_quality_rate: f64,
    pub unmapped_rate: f64,
    pub function_seconds: BTreeMap<String, f64>,
}

#[derive(Debug, Clone, Default, PartialEq)]
pub struct PerfSummary {
    pub counters: BTreeMap<String, f64>,
}

impl PerfSummary {
    pub fn l1_miss_rate(&self) -> Option<f64> {
        ratio(
            self.counters.get("L1-dcache-load-misses").copied(),
            self.counters.get("L1-dcache-loads").copied(),
        )
    }

    pub fn llc_miss_rate(&self) -> Option<f64> {
        ratio(
            self.counters.get("LLC-load-misses").copied(),
            self.counters.get("LLC-loads").copied(),
        )
    }
}

pub fn parse_perf_stat(text: &str) -> PerfSummary {
    let mut counters = BTreeMap::new();
    for line in text.lines() {
        let trimmed = line.trim();
        if trimmed.is_empty() || trimmed.starts_with('#') {
            continue;
        }
        if let Some((value, event)) =
            parse_csv_perf_line(trimmed).or_else(|| parse_text_perf_line(trimmed))
        {
            counters.insert(event, value);
        }
    }
    PerfSummary { counters }
}

fn parse_csv_perf_line(line: &str) -> Option<(f64, String)> {
    let fields: Vec<_> = line.split(',').collect();
    if fields.len() < 3 {
        return None;
    }
    let value = parse_counter(fields[0])?;
    let event = fields[2].trim().to_owned();
    if event.is_empty() {
        None
    } else {
        Some((value, event))
    }
}

fn parse_text_perf_line(line: &str) -> Option<(f64, String)> {
    let mut fields = line.split_whitespace();
    let value = parse_counter(fields.next()?)?;
    let event = fields.next()?.to_owned();
    Some((value, event))
}

fn parse_counter(text: &str) -> Option<f64> {
    let normalized = text.replace(',', "");
    normalized.parse().ok()
}

fn ratio(num: Option<f64>, den: Option<f64>) -> Option<f64> {
    let (Some(num), Some(den)) = (num, den) else {
        return None;
    };
    if den == 0.0 {
        None
    } else {
        Some(num / den)
    }
}

impl MapMetrics {
    pub fn finish_rates(&mut self, elapsed: Duration, total_bases: u64) {
        self.wall_seconds = elapsed.as_secs_f64();
        if self.wall_seconds > 0.0 {
            self.reads_per_second = self.num_reads as f64 / self.wall_seconds;
            self.bases_per_second = total_bases as f64 / self.wall_seconds;
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn parses_perf_stat_text_rates() {
        let perf = parse_perf_stat(
            "1,000 L1-dcache-loads\n100 L1-dcache-load-misses\n200 LLC-loads\n10 LLC-load-misses\n",
        );
        assert_eq!(perf.l1_miss_rate(), Some(0.1));
        assert_eq!(perf.llc_miss_rate(), Some(0.05));
    }
}

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
pub struct CandidateGenerationStats {
    pub seed_lookups: u64,
    pub candidate_hits: u64,
    pub postings_skipped_due_to_frequency: u64,
    pub seeds_skipped_due_to_quality: u64,
    pub seed_candidates_considered: u64,
    pub early_stopped_reads: u64,
    pub seed_lookups_saved_by_early_stop: u64,
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
) -> (Vec<CandidateHit>, CandidateGenerationStats) {
    let encoded = encode_acgt(seq);
    let seed_limit = max_seeds_per_read.max(1);
    let mut stats = CandidateGenerationStats::default();
    let mut seed_choices: SmallVec<[(usize, u32, u64, u8); 256]> = SmallVec::new();

    collect_seed_choices(
        index,
        &encoded,
        qual,
        min_seed_quality,
        0,
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
            &mut stats,
            &mut seed_choices,
        );
    }

    seed_choices.sort_unstable();

    let mut grouped: FxHashMap<u64, u16> =
        FxHashMap::with_capacity_and_hasher(seed_limit.saturating_mul(4), Default::default());

    let selected_seed_count = seed_choices.len().min(seed_limit);
    for (selected_idx, (postings_len, seed_pos, seed_code, read_strand)) in
        seed_choices.into_iter().take(seed_limit).enumerate()
    {
        stats.seed_lookups += 1;
        if postings_len == usize::MAX {
            continue;
        }
        if postings_len > max_postings_per_seed {
            stats.postings_skipped_due_to_frequency += postings_len as u64;
            continue;
        }
        for posting in index.seed_postings(seed_code) {
            if posting.pos < seed_pos {
                continue;
            }
            let candidate_start = posting.pos - seed_pos;
            let key = pack_candidate_key(posting.transcript_id, candidate_start, read_strand);
            *grouped.entry(key).or_default() += 1;
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
        .map(|(key, seed_count)| {
            let (transcript_id, pos, strand) = unpack_candidate_key(key);
            CandidateHit {
                read_id,
                transcript_id,
                gene_id: index.transcript_gene_id(transcript_id).unwrap_or(0),
                pos,
                strand,
                seed_count,
                seed_score: seed_count,
            }
        })
        .collect();
    stats.candidate_hits = hits.len() as u64;
    (hits, stats)
}

fn collect_seed_choices(
    index: &dyn IndexAccess,
    encoded: &crate::dna::EncodedSeq,
    qual: Option<&[u8]>,
    min_seed_quality: u8,
    read_strand: u8,
    stats: &mut CandidateGenerationStats,
    seed_choices: &mut SmallVec<[(usize, u32, u64, u8); 256]>,
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
        let postings_len = index.posting_count(seed.code);
        let frequency_rank = if postings_len == 0 {
            usize::MAX
        } else {
            postings_len
        };
        seed_choices.push((frequency_rank, seed.pos, seed.code, read_strand));
    }
}

fn should_stop_early(
    grouped: &FxHashMap<u64, u16>,
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
    for &count in grouped.values() {
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
}

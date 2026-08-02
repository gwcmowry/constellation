use std::collections::{BTreeMap, HashMap};

use crate::dna::{kmer_keys, normalize_seq};
use crate::index::{intersect_pathsets, pathset_size, CertificateStore, PathCertIndex};
use crate::io::fastq::FastqRecord;

pub const FLAG_AMBIGUOUS_PATHSET: u32 = 1;
pub const FLAG_MULTI_COMPONENT: u32 = 1 << 1;

#[derive(Clone, Debug)]
pub struct MappingConfig {
    pub max_postings_per_seed: usize,
    pub min_seed_hits: usize,
}

impl Default for MappingConfig {
    fn default() -> Self {
        Self {
            max_postings_per_seed: 256,
            min_seed_hits: 1,
        }
    }
}

#[derive(Clone, Debug, Eq, PartialEq)]
pub enum MappingStatus {
    Mapped,
    MultiComponent,
    Unmapped,
}

#[derive(Clone, Debug, Eq, PartialEq)]
pub struct MappingRecord {
    pub read_id: String,
    pub status: MappingStatus,
    pub certificate_ids: Vec<u32>,
}

#[derive(Clone, Debug, Default, Eq, PartialEq)]
pub struct MappingSummary {
    pub reads: u64,
    pub mapped: u64,
    pub multi_component: u64,
    pub unmapped: u64,
    pub seed_lookups: u64,
    pub candidate_components: u64,
    pub pathset_intersections: u64,
}

pub fn map_records(
    index: &PathCertIndex,
    records: &[FastqRecord],
    config: &MappingConfig,
) -> (Vec<MappingRecord>, CertificateStore, MappingSummary) {
    let mut certificates = CertificateStore::default();
    let mut summary = MappingSummary::default();
    let mut output = Vec::with_capacity(records.len());

    for record in records {
        let mapping = map_one(index, record, config, &mut certificates, &mut summary);
        output.push(mapping);
    }

    (output, certificates, summary)
}

pub fn map_records_component_bucketed(
    index: &PathCertIndex,
    records: &[FastqRecord],
    config: &MappingConfig,
) -> (Vec<MappingRecord>, CertificateStore, MappingSummary) {
    let mut summary = MappingSummary {
        reads: records.len() as u64,
        ..MappingSummary::default()
    };
    let mut component_buckets: BTreeMap<u32, Vec<(usize, u64)>> = BTreeMap::new();

    for (read_idx, record) in records.iter().enumerate() {
        let seq = normalize_seq(&record.seq);
        for key in kmer_keys(&seq, index.k) {
            summary.seed_lookups += 1;
            let postings = index.seed_index.lookup(key);
            if postings.is_empty() || postings.len() > config.max_postings_per_seed {
                continue;
            }
            for posting in postings {
                component_buckets
                    .entry(posting.component_id)
                    .or_default()
                    .push((read_idx, posting.pathset));
            }
        }
    }

    let mut resolved_by_read: Vec<Vec<(u32, u64)>> = vec![Vec::new(); records.len()];
    for (component_id, work_items) in component_buckets {
        let Some(component) = index.component(component_id) else {
            continue;
        };
        let all_paths = component.all_paths_mask();
        let mut local: HashMap<usize, (u64, usize)> = HashMap::new();
        for (read_idx, pathset) in work_items {
            let entry = local.entry(read_idx).or_insert((all_paths, 0));
            entry.0 = intersect_pathsets(entry.0, pathset);
            entry.1 += 1;
            summary.pathset_intersections += 1;
        }
        for (read_idx, (pathset, seed_hits)) in local {
            if pathset != 0 && seed_hits >= config.min_seed_hits {
                resolved_by_read[read_idx].push((component_id, pathset));
            }
        }
    }

    let mut certificates = CertificateStore::default();
    let mut mappings = Vec::with_capacity(records.len());
    for (read_idx, resolved) in resolved_by_read.into_iter().enumerate() {
        summary.candidate_components += resolved.len() as u64;
        if resolved.is_empty() {
            summary.unmapped += 1;
            mappings.push(MappingRecord {
                read_id: records[read_idx].id.clone(),
                status: MappingStatus::Unmapped,
                certificate_ids: Vec::new(),
            });
            continue;
        }

        let multi_component = resolved.len() > 1;
        if multi_component {
            summary.multi_component += 1;
        } else {
            summary.mapped += 1;
        }

        let mut certificate_ids = Vec::with_capacity(resolved.len());
        for (component_id, pathset) in resolved {
            let mut flags = 0;
            if pathset_size(pathset) > 1 {
                flags |= FLAG_AMBIGUOUS_PATHSET;
            }
            if multi_component {
                flags |= FLAG_MULTI_COMPONENT;
            }
            certificate_ids.push(certificates.intern(component_id, pathset, flags));
        }
        mappings.push(MappingRecord {
            read_id: records[read_idx].id.clone(),
            status: if multi_component {
                MappingStatus::MultiComponent
            } else {
                MappingStatus::Mapped
            },
            certificate_ids,
        });
    }

    (mappings, certificates, summary)
}

pub fn map_one(
    index: &PathCertIndex,
    record: &FastqRecord,
    config: &MappingConfig,
    certificates: &mut CertificateStore,
    summary: &mut MappingSummary,
) -> MappingRecord {
    summary.reads += 1;
    let seq = normalize_seq(&record.seq);
    let seed_keys = kmer_keys(&seq, index.k);
    if seed_keys.is_empty() {
        summary.unmapped += 1;
        return MappingRecord {
            read_id: record.id.clone(),
            status: MappingStatus::Unmapped,
            certificate_ids: Vec::new(),
        };
    }

    let mut by_component: HashMap<u32, Vec<u64>> = HashMap::new();
    for key in seed_keys {
        summary.seed_lookups += 1;
        let postings = index.seed_index.lookup(key);
        if postings.is_empty() || postings.len() > config.max_postings_per_seed {
            continue;
        }
        for posting in postings {
            by_component
                .entry(posting.component_id)
                .or_default()
                .push(posting.pathset);
        }
    }

    let mut resolved = BTreeMap::new();
    for (component_id, pathsets) in by_component {
        if pathsets.len() < config.min_seed_hits {
            continue;
        }
        let Some(component) = index.component(component_id) else {
            continue;
        };
        let mut candidate = component.all_paths_mask();
        for pathset in pathsets {
            summary.pathset_intersections += 1;
            candidate = intersect_pathsets(candidate, pathset);
            if candidate == 0 {
                break;
            }
        }
        if candidate != 0 {
            resolved.insert(component_id, candidate);
        }
    }

    summary.candidate_components += resolved.len() as u64;
    if resolved.is_empty() {
        summary.unmapped += 1;
        return MappingRecord {
            read_id: record.id.clone(),
            status: MappingStatus::Unmapped,
            certificate_ids: Vec::new(),
        };
    }

    let multi_component = resolved.len() > 1;
    if multi_component {
        summary.multi_component += 1;
    } else {
        summary.mapped += 1;
    }

    let mut certificate_ids = Vec::new();
    for (component_id, pathset) in resolved {
        let mut flags = 0;
        if pathset_size(pathset) > 1 {
            flags |= FLAG_AMBIGUOUS_PATHSET;
        }
        if multi_component {
            flags |= FLAG_MULTI_COMPONENT;
        }
        certificate_ids.push(certificates.intern(component_id, pathset, flags));
    }

    MappingRecord {
        read_id: record.id.clone(),
        status: if multi_component {
            MappingStatus::MultiComponent
        } else {
            MappingStatus::Mapped
        },
        certificate_ids,
    }
}

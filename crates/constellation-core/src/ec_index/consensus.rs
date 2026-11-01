//! Conservative forward-strand assignment from agreement among distinct seeds.
use super::*;

#[derive(Debug, Default)]
pub struct EcConsensusScratch {
    counts: Vec<u32>,
    classes: Vec<u8>,
    touched: Vec<GeneId>,
    keys: Vec<u64>,
    states: Vec<SplicingState>,
}

impl EcConsensusScratch {
    fn reset(&mut self, genes: usize) {
        for gene in self.touched.drain(..) {
            self.counts[gene as usize] = 0;
            self.classes[gene as usize] = 0;
        }
        self.counts.resize(genes, 0);
        self.classes.resize(genes, 0);
        self.keys.clear();
    }
}

impl LoadedEcIndex {
    /// Intersect the gene sets of all usable, distinct matched seeds. At least
    /// two different matched k-mers are required. Missing
    /// seeds are tolerated (sequencing errors); contradictory matched seeds are
    /// not. Class annotations never resolve a tie between different genes.
    ///
    /// This is a forward-strand policy; callers must not pool reverse-complement
    /// evidence into the same intersection. It needs no new index format.
    pub fn map_encoded_consensus_into(
        &self,
        encoded: &EncodedSeq,
        stride: u32,
        max_gene_df: u32,
        scratch: &mut EcConsensusScratch,
        mapping: &mut EcResolvedMapping,
    ) -> EcMappingStats {
        self.map_consensus_impl::<false>(encoded, stride, max_gene_df, scratch, mapping)
    }

    /// Preserve class compatibility separately for every surviving gene.
    pub fn map_encoded_splicing_into(
        &self,
        encoded: &EncodedSeq,
        stride: u32,
        max_gene_df: u32,
        scratch: &mut EcConsensusScratch,
        mapping: &mut EcResolvedMapping,
    ) -> EcMappingStats {
        self.map_consensus_impl::<true>(encoded, stride, max_gene_df, scratch, mapping)
    }

    fn map_consensus_impl<const KEEP_STATES: bool>(
        &self,
        encoded: &EncodedSeq,
        stride: u32,
        max_gene_df: u32,
        scratch: &mut EcConsensusScratch,
        mapping: &mut EcResolvedMapping,
    ) -> EcMappingStats {
        scratch.reset(self.num_genes());
        if KEEP_STATES {
            scratch
                .states
                .resize(self.num_genes(), SplicingState::Unknown);
        }
        let mut top = std::mem::take(&mut mapping.top_gene_ids);
        top.clear();
        *mapping = EcResolvedMapping {
            top_gene_ids: top,
            ..Default::default()
        };
        let mut stats = EcMappingStats::default();
        for_each_selected_kmer(encoded, self.k(), self.sampling(), stride.max(1), |kmer| {
            stats.seed_lookups += 1;
            let Some(hit) = self.lookup_hit(kmer.code) else {
                return;
            };
            if hit.gene_df > max_gene_df.max(1) {
                stats.seeds_skipped_high_gene_df += 1;
                return;
            }
            if scratch.keys.contains(&kmer.code) {
                return;
            }
            scratch.keys.push(kmer.code);
            self.visit_lookup_hit(hit, max_gene_df.max(1), &mut stats, |gene, _, class| {
                let i = gene as usize;
                if scratch.counts[i] == 0 {
                    scratch.touched.push(gene);
                    if KEEP_STATES {
                        scratch.states[i] = SplicingState::from_class_bits(class);
                    }
                } else if KEEP_STATES {
                    scratch.states[i] =
                        scratch.states[i].intersect(SplicingState::from_class_bits(class));
                }
                scratch.counts[i] += 1;
                scratch.classes[i] |= class;
            });
        });
        mapping.gene_candidate_count = scratch.touched.len() as u32;
        if scratch.keys.len() < 2 {
            return stats;
        }
        let required = scratch.keys.len() as u32;
        for &gene in &scratch.touched {
            if scratch.counts[gene as usize] == required {
                mapping.top_gene_ids.push(gene);
                let class = scratch.classes[gene as usize];
                mapping.has_exon_support |= class & EC_CLASS_EXON != 0;
                mapping.has_intron_support |= class & EC_CLASS_INTRON != 0;
                mapping.has_antisense_support |= class & EC_CLASS_ANTISENSE != 0;
            }
        }
        mapping.top_gene_ids.sort_unstable();
        if !mapping.top_gene_ids.is_empty() {
            mapping.score = required.min(u16::MAX as u32) as u16;
            mapping.best_score = mapping.score;
            stats.reads_with_gene_support = 1;
            if mapping.top_gene_ids.len() == 1 {
                mapping.gene_id = Some(mapping.top_gene_ids[0]);
            }
        }
        stats
    }
}

impl EcConsensusScratch {
    /// Only query genes in the result of the most recent splicing-aware call.
    pub fn splicing_state(&self, gene: GeneId) -> SplicingState {
        self.states[gene as usize]
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Write;

    #[test]
    fn shared_sequence_stays_ambiguous_and_scratch_does_not_leak() {
        let mut fasta = tempfile::NamedTempFile::new().unwrap();
        writeln!(fasta, ">a|gene=A exon\nACGTGCAATGACTGGTACAGCTAG\n>b|gene=B intron\nACGTGCAATGACTGGTACAGCTAG\n>c|gene=C exon\nTGCCTAGGATCCGATGTCAGGTCA").unwrap();
        let index = LoadedEcIndex::Owned(EcIndex::build_from_fasta(fasta.path(), 7).unwrap());
        let mut scratch = EcConsensusScratch::default();
        let mut output = EcResolvedMapping::default();
        for _ in 0..2 {
            index.map_encoded_consensus_into(
                &encode_acgt(b"ACGTGCAATGACTGGTACAGCTAG"),
                1,
                256,
                &mut scratch,
                &mut output,
            );
            assert!(output.gene_id.is_none());
            assert_eq!(output.top_gene_ids.len(), 2);
            assert!(output.has_exon_support && output.has_intron_support);
            index.map_encoded_consensus_into(
                &encode_acgt(b"TGCCTAGGATCCGATGTCAGGTCA"),
                1,
                256,
                &mut scratch,
                &mut output,
            );
            assert_eq!(output.top_gene_ids.len(), 1);
            assert_eq!(index.gene_name(output.gene_id.unwrap()), Some("C"));
            index.map_encoded_consensus_into(
                &encode_acgt(b"ACGTGCANNNNNNNNNNNNNNN"),
                1,
                256,
                &mut scratch,
                &mut output,
            );
            assert_eq!(output.best_score, 0);
            assert!(output.top_gene_ids.is_empty());
            assert!(!output.has_exon_support && !output.has_intron_support);
        }
    }

    #[test]
    fn independent_intersection_oracle_across_layouts_and_repeated_seeds() {
        let mut fasta = tempfile::NamedTempFile::new().unwrap();
        writeln!(fasta, ">a|gene=A exon\nACGTGCAATGACTGGTACAGCTAG\n>b|gene=B intron\nACGTGCAATGACTGGTACAGCTAG\n>r|gene=R exon\nAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA").unwrap();
        let index = EcIndex::build_from_fasta(fasta.path(), 7).unwrap();
        let path = tempfile::NamedTempFile::new().unwrap();
        index.save_mmap(path.path()).unwrap();
        for loaded in [
            LoadedEcIndex::Owned(index),
            LoadedEcIndex::load(path.path()).unwrap(),
        ] {
            let mut scratch = EcConsensusScratch::default();
            let mut output = EcResolvedMapping::default();
            for seq in [
                b"ACGTGCAATGACTGGTACAGCTAG".as_slice(),
                b"AAAAAAAANNNACGTGCAATGACTGGTACAGCTAG",
                b"AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA",
                b"NN",
                b"",
            ] {
                let encoded = encode_acgt(seq);
                let mut keys = std::collections::HashSet::new();
                let mut sets = Vec::new();
                let mut positions = Vec::new();
                for kmer in iter_kmers_2bit(&encoded, 7) {
                    if let Some(record) = loaded.lookup(kmer.code) {
                        if keys.insert(kmer.code) {
                            let mut genes = std::collections::HashSet::new();
                            loaded.for_each_ec_target(record.ec_id, |gene, _| {
                                genes.insert(gene);
                            });
                            sets.push(genes);
                            positions.push(kmer.pos);
                        }
                    }
                }
                let mut expected = Vec::new();
                if sets.len() >= 2 {
                    expected = sets[0]
                        .iter()
                        .copied()
                        .filter(|g| sets.iter().all(|s| s.contains(g)))
                        .collect();
                    expected.sort_unstable();
                }
                loaded.map_encoded_consensus_into(&encoded, 1, 256, &mut scratch, &mut output);
                assert_eq!(output.top_gene_ids, expected, "{seq:?}");
            }
        }
    }
}

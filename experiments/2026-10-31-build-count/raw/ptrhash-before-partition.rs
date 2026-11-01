//! Direct v9 construction avoids materializing and decoding a v4 prefix index.
use super::*;
use rayon::prelude::*;

impl EcIndex {
    pub fn save_ptrhash(&self, path: impl AsRef<Path>) -> Result<(), EcIndexError> {
        if !(1..=32).contains(&self.k) {
            return Err(EcIndexError::InvalidK(self.k));
        }
        self.sampling.validate(self.k)?;
        // PtrHash cannot construct an empty key set; v4 is an exact supported
        // empty-index representation and needs no lookup acceleration.
        if self.kmers.is_empty() {
            return self.save_mmap(path);
        }
        let started = std::time::Instant::now();
        let keys: Vec<_> = self.kmers.par_iter().map(|r| r.kmer_code).collect();
        let ptrhash = EcPtrHash::try_new(&keys, PtrHashParams::default_fast())
            .ok_or_else(|| EcIndexError::PtrHash("construction failed after ten seeds".into()))?;
        build_phase("ptrhash_construct", started);
        let started = std::time::Instant::now();
        let mut slots = vec![EcPtrHashSlot::default(); ptrhash.max_index()];
        let payload = |record: &EcKmerRecord| -> Result<EcPtrHashSlot, EcIndexError> {
            let start = self.ec_offsets[record.ec_id as usize];
            let end = self.ec_offsets[record.ec_id as usize + 1];
            let len = end
                .checked_sub(start)
                .filter(|&n| n > 0)
                .ok_or(EcIndexError::Truncated)?;
            Ok(EcPtrHashSlot {
                kmer_code: record.kmer_code,
                target_start: start,
                gene_df: u16::try_from(record.gene_df)
                    .map_err(|_| EcIndexError::TooLarge("EC v9 gene df"))?,
                target_len: u16::try_from(len)
                    .map_err(|_| EcIndexError::TooLarge("EC v9 target length"))?,
            })
        };
        ptrhash.index_stream::<32, _>(&keys).enumerate().fold(
            Ok(()),
            |result, (i, slot)| -> Result<(), EcIndexError> {
                result?;
                let target = slots.get_mut(slot).ok_or(EcIndexError::Truncated)?;
                if target.target_len != 0 {
                    return Err(EcIndexError::PtrHash(
                        "two retained keys mapped to one slot".into(),
                    ));
                }
                *target = payload(&self.kmers[i])?;
                Ok(())
            },
        )?;
        build_phase("ptrhash_fill", started);
        let started = std::time::Instant::now();
        self.kmers
            .par_chunks(65536)
            .try_for_each(|chunk| -> Result<(), EcIndexError> {
                ptrhash
                    .index_stream::<32, _>(chunk.iter().map(|r| r.kmer_code))
                    .enumerate()
                    .fold(Ok(()), |result, (i, slot)| -> Result<(), EcIndexError> {
                        result?;
                        if slots[slot] != payload(&chunk[i])? {
                            return Err(EcIndexError::PtrHash(
                                "positive-key verification failed".into(),
                            ));
                        }
                        Ok(())
                    })
            })?;
        // Bounded even for tiny, fully populated k-mer universes.
        let mask = if self.k == 32 {
            u64::MAX
        } else {
            (1u64 << (2 * self.k)) - 1
        };
        let mut state = 0x6a09_e667_f3bc_c909u64;
        for _ in 0..1_000_000 {
            state = splitmix64(state);
            let key = state & mask;
            if keys.binary_search(&key).is_err() {
                let slot = slots[ptrhash.index(&key)];
                if slot.target_len != 0 && slot.kmer_code == key {
                    return Err(EcIndexError::PtrHash(
                        "negative-key verification failed".into(),
                    ));
                }
            }
        }
        build_phase("ptrhash_verify", started);
        let started = std::time::Instant::now();
        let mut names = Vec::new();
        let mut metas = Vec::new();
        for name in &self.genes {
            metas.push(EcGeneMeta {
                name_start: u32::try_from(names.len())
                    .map_err(|_| EcIndexError::TooLarge("gene names"))?,
                name_len: u32::try_from(name.len())
                    .map_err(|_| EcIndexError::TooLarge("gene name"))?,
            });
            names.extend_from_slice(name.as_bytes());
        }
        let mut targets = Vec::new();
        write_v4_targets(&mut targets, &self.ec_gene_ids, &self.ec_class_bits)?;
        write_ptrhash_index(
            path,
            self.k,
            self.sampling,
            keys.len(),
            &ptrhash,
            &slots,
            &metas,
            &names,
            &self.ec_offsets,
            &targets,
        )?;
        build_phase("ptrhash_write", started);
        Ok(())
    }
}

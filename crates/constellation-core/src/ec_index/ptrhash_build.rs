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
        // Permute record references into disjoint slot ranges. This replaces
        // one latency-bound scatter over the entire index with independent
        // smaller scatters, without unsafe aliased mutable access.
        if self.kmers.len() > u32::MAX as usize || slots.len() > u32::MAX as usize {
            return Err(EcIndexError::TooLarge("PtrHash partition record/slot IDs"));
        }
        let width = slots.len().div_ceil(256);
        let parts = slots.len().div_ceil(width);
        let chunks: Vec<Vec<Vec<u64>>> = self
            .kmers
            .par_chunks(1 << 20)
            .enumerate()
            .map(|(chunk_id, chunk)| {
                let mut buckets: Vec<_> = (0..parts)
                    .map(|_| Vec::with_capacity(chunk.len().div_ceil(parts)))
                    .collect();
                ptrhash
                    .index_stream::<32, _>(chunk.iter().map(|r| r.kmer_code))
                    .enumerate()
                    .for_each(|(i, slot)| {
                        let record_id = chunk_id * (1 << 20) + i;
                        buckets[slot / width].push(((slot as u64) << 32) | record_id as u64);
                    });
                buckets
            })
            .collect();
        let mut buckets: Vec<Vec<u64>> = (0..parts)
            .map(|i| Vec::with_capacity(chunks.iter().map(|c| c[i].len()).sum()))
            .collect();
        for chunk in chunks {
            for (dest, mut source) in buckets.iter_mut().zip(chunk) {
                dest.append(&mut source);
            }
        }
        build_phase("ptrhash_partition", started);
        let started = std::time::Instant::now();
        slots
            .par_chunks_mut(width)
            .zip(buckets.into_par_iter())
            .enumerate()
            .try_for_each(|(part, (slots, records))| -> Result<(), EcIndexError> {
                for packed in records {
                    let slot = (packed >> 32) as usize - part * width;
                    let target = slots.get_mut(slot).ok_or(EcIndexError::Truncated)?;
                    if target.target_len != 0 {
                        return Err(EcIndexError::PtrHash(
                            "two retained keys mapped to one slot".into(),
                        ));
                    }
                    *target = payload(&self.kmers[packed as u32 as usize])?;
                }
                Ok(())
            })?;
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

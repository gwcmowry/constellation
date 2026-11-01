//! Read-only diagnostic of the v9 index and a bounded real-read trace.
//! Run with: INDEX R1.fastq R2.fastq READ_COUNT OUTPUT.json
//! Microbenchmarks exclude FASTQ I/O and are not end-to-end throughput claims.
use constellation_core::dna::{encode_acgt, iter_closed_syncmers_2bit, EncodedSeq};
use constellation_core::ec_index::{EcScoreScratch, LoadedEcIndex};
use constellation_core::sketch::is_low_complexity_read;
use epserde::prelude::*;
use memmap2::Mmap;
use ptr_hash::{bucket_fn::Linear, hash::Xxh3Int, PtrHash};
use rayon::prelude::*;
use rustc_hash::{FxHashMap, FxHashSet};
use serde_json::json;
use std::fs::File;
use std::hint::black_box;
use std::io::{BufRead, BufReader, Cursor};
use std::time::Instant;

type Hash = PtrHash<u64, Linear, Vec<u32>, Xxh3Int, Vec<u8>, false, false>;

#[repr(C)]
#[derive(Clone, Copy, bytemuck::Pod, bytemuck::Zeroable)]
struct Slot {
    key: u64,
    start: u32,
    df: u16,
    len: u16,
}

// This diagnostic explicitly supports the existing little-endian v9 layout.
// Keep all experimental layout decoding outside the production mapper.
struct View {
    mmap: Mmap,
    hash: Hash,
    slots: usize,
    nslots: usize,
    target_offset: usize,
    target_count: usize,
    k: u8,
    s: u8,
}
impl View {
    fn load(path: &str) -> Self {
        let mmap = unsafe { Mmap::map(&File::open(path).unwrap()).unwrap() };
        assert_eq!(&mmap[..8], b"CSTLEC2\0");
        let u32_at = |p| u32::from_le_bytes(mmap[p..p + 4].try_into().unwrap());
        let u64_at = |p| u64::from_le_bytes(mmap[p..p + 8].try_into().unwrap()) as usize;
        assert_eq!(u32_at(8), 9);
        assert_eq!(&mmap[16..20], b"SYNH");
        let hash_start = u64_at(64);
        let hash_len = u64_at(72);
        let hash = unsafe {
            Hash::deserialize_full(&mut Cursor::new(&mmap[hash_start..hash_start + hash_len]))
                .unwrap()
        };
        let view = Self {
            hash,
            slots: u64_at(88),
            nslots: u64_at(40),
            target_offset: u64_at(104),
            target_count: u64_at(56),
            k: u32_at(12) as u8,
            s: u32_at(20) as u8,
            mmap,
        };
        assert_eq!(view.hash.max_index(), view.nslots);
        let _ = view.slot_slice();
        assert!(view.target_offset + view.target_count * 5 <= view.mmap.len());
        view
    }
    fn slot_slice(&self) -> &[Slot] {
        bytemuck::cast_slice(&self.mmap[self.slots..self.slots + self.nslots * 16])
    }
    fn lookup(&self, key: u64) -> Option<Slot> {
        let slot = self.slot_slice()[self.hash.index(&key)];
        (slot.key == key && slot.len > 0).then_some(slot)
    }
    fn signature(slot: Slot) -> u64 {
        ((slot.start as u64) << 32) | ((slot.df as u64) << 16) | slot.len as u64
    }
}

struct Read {
    seq: Vec<u8>,
    encoded: EncodedSeq,
    eligible: bool,
}

fn record(reader: &mut impl BufRead) -> Option<[Vec<u8>; 4]> {
    let mut lines = std::array::from_fn(|_| Vec::new());
    for (i, line) in lines.iter_mut().enumerate() {
        if reader.read_until(b'\n', line).unwrap() == 0 {
            assert_eq!(i, 0, "truncated FASTQ");
            return None;
        }
        while matches!(line.last(), Some(b'\r' | b'\n')) {
            line.pop();
        }
    }
    assert!(lines[0].starts_with(b"@") && lines[2].starts_with(b"+"));
    assert_eq!(lines[1].len(), lines[3].len());
    Some(lines)
}

fn timed(mut run: impl FnMut() -> u64) -> serde_json::Value {
    let warmup = Instant::now();
    let expected = black_box(run()); // warmup
    let iterations = (0.3 / warmup.elapsed().as_secs_f64())
        .ceil()
        .clamp(1.0, 64.0) as usize;
    let mut seconds = Vec::new();
    let mut sample_seconds = Vec::new();
    for _ in 0..3 {
        let start = Instant::now();
        for _ in 0..iterations {
            assert_eq!(black_box(run()), expected);
        }
        let elapsed = start.elapsed().as_secs_f64();
        sample_seconds.push(elapsed);
        seconds.push(elapsed / iterations as f64);
    }
    json!({"seconds": seconds, "sample_seconds": sample_seconds,
           "iterations_per_sample": iterations, "checksum": expected})
}

fn slot_checksum(key: u64, slot: Slot) -> u64 {
    if slot.key == key && slot.len > 0 {
        View::signature(slot).wrapping_add(key)
    } else {
        0
    }
}

#[inline]
fn prefetch_slot(slots: &[Slot], idx: usize) {
    // idx comes from the validated PtrHash with max_index == slots.len().
    #[cfg(target_arch = "x86_64")]
    unsafe {
        std::arch::x86_64::_mm_prefetch(
            slots.as_ptr().add(idx).cast(),
            std::arch::x86_64::_MM_HINT_T0,
        );
    }
}

fn main() {
    let args: Vec<_> = std::env::args().collect();
    assert_eq!(
        args.len(),
        6,
        "INDEX R1.fastq R2.fastq READ_COUNT OUTPUT.json"
    );
    let limit: usize = args[4].parse().unwrap();
    let index = LoadedEcIndex::load(&args[1]).unwrap();
    let view = View::load(&args[1]);
    let mut r1 = BufReader::new(File::open(&args[2]).unwrap());
    let mut r2 = BufReader::new(File::open(&args[3]).unwrap());
    let mut reads = Vec::with_capacity(limit);
    for _ in 0..limit {
        let a = record(&mut r1).expect("R1 shorter than requested trace");
        let b = record(&mut r2).expect("R2 shorter than requested trace");
        let seq = b[1].to_ascii_uppercase();
        let valid_r1 = a[1].len() >= 28
            && a[1][..28]
                .iter()
                .all(|b| matches!(b.to_ascii_uppercase(), b'A' | b'C' | b'G' | b'T' | b'U'));
        let quality = b[3]
            .iter()
            .map(|q| q.saturating_sub(33) as u64)
            .sum::<u64>() as f64
            / b[3].len().max(1) as f64;
        let eligible = valid_r1 && quality >= 10.0 && !is_low_complexity_read(&seq, view.k);
        reads.push(Read {
            encoded: encode_acgt(&seq),
            seq,
            eligible,
        });
    }
    let mut keys = Vec::new();
    let mut frequencies = FxHashMap::<u64, u32>::default();
    let mut training_frequencies = FxHashMap::<u64, u32>::default();
    let mut sequences = FxHashSet::default();
    let mut batch_sequences = FxHashSet::default();
    let mut batch_unique_sequences = 0_u64;
    let mut training_occurrences = 0_u64;
    let mut hits = 0_u64;
    let mut capped = 0_u64;
    let mut singleton = 0_u64;
    let mut votes = 0_u64;
    let mut run_votes = 0_u64;
    let mut unique_ec_votes = 0_u64;
    let mut runs = 0_u64;
    let mut unique_ecs = 0_u64;
    let mut eligible_count = 0_u64;
    let mut one_ec_reads = 0_u64;
    let mut supported_reads = 0_u64;
    for (read_idx, read) in reads.iter().enumerate() {
        if read_idx % 65_536 == 0 {
            batch_unique_sequences += batch_sequences.len() as u64;
            batch_sequences.clear();
        }
        if !read.eligible {
            continue;
        }
        eligible_count += 1;
        sequences.insert(read.seq.as_slice());
        batch_sequences.insert(read.seq.as_slice());
        let mut last = None;
        let mut seen = FxHashSet::default();
        for kmer in iter_closed_syncmers_2bit(&read.encoded, view.k, view.s) {
            keys.push(kmer.code);
            *frequencies.entry(kmer.code).or_default() += 1;
            if read_idx < limit / 2 {
                *training_frequencies.entry(kmer.code).or_default() += 1;
                training_occurrences += 1;
            }
            let slot = view.lookup(kmer.code);
            assert_eq!(slot.is_some(), index.lookup_hit(kmer.code).is_some());
            let Some(slot) = slot else {
                last = None;
                continue;
            };
            if slot.df > 256 {
                capped += 1;
                last = None;
                continue;
            }
            hits += 1;
            votes += slot.len as u64;
            singleton += u64::from(slot.len == 1);
            let signature = View::signature(slot);
            if last != Some(signature) {
                runs += 1;
                run_votes += slot.len as u64;
            }
            if seen.insert(signature) {
                unique_ecs += 1;
                unique_ec_votes += slot.len as u64;
            }
            last = Some(signature);
        }
        supported_reads += u64::from(!seen.is_empty());
        one_ec_reads += u64::from(seen.len() == 1);
    }
    batch_unique_sequences += batch_sequences.len() as u64;
    let mut training: Vec<_> = training_frequencies
        .iter()
        .map(|(&key, &count)| (key, count))
        .collect();
    training.sort_unstable_by_key(|&(key, count)| (std::cmp::Reverse(count), key));
    let holdout_coverage: Vec<_> = [1_000, 10_000, 100_000, 1_000_000].into_iter().map(|n| {
        json!({"keys": n, "holdout_occurrences": training.iter().take(n).map(|&(key, count)| (frequencies[&key] - count) as u64).sum::<u64>()})
    }).collect();
    for chunk in keys.chunks(4096) {
        let checked = view
            .hash
            .index_stream::<32, _>(chunk)
            .fold(0_usize, |i, pos| {
                assert_eq!(pos, view.hash.index(&chunk[i]));
                i + 1
            });
        assert_eq!(checked, chunk.len());
    }
    let mut counts: Vec<_> = frequencies.values().copied().collect();
    counts.sort_unstable_by(|a, b| b.cmp(a));
    let coverage: Vec<_> = [1_000, 10_000, 100_000, 1_000_000].into_iter().map(|n| {
        json!({"keys": n, "occurrences": counts.iter().take(n).map(|&x| x as u64).sum::<u64>()})
    }).collect();
    let trace = json!({
        "reads": limit, "eligible_reads": eligible_count, "supported_reads": supported_reads,
        "unique_eligible_sequences": sequences.len(), "seed_occurrences": keys.len(),
        "unique_eligible_sequences_within_65536_read_batches": batch_unique_sequences,
        "unique_keys": counts.len(), "singleton_keys": counts.iter().filter(|&&x| x == 1).count(),
        "accepted_hits": hits, "capped_hits": capped,
        "absent_keys": keys.len() as u64 - hits - capped, "singleton_target_hits": singleton,
        "gene_votes": votes, "same_ec_runs": runs, "run_gene_visits": run_votes,
        "distinct_ecs_per_read_sum": unique_ecs, "distinct_ec_gene_visits": unique_ec_votes,
        "reads_with_one_distinct_accepted_ec": one_ec_reads, "hot_key_coverage": coverage,
        "training_reads": limit / 2, "holdout_seed_occurrences": keys.len() as u64 - training_occurrences,
        "hot_key_holdout_coverage": holdout_coverage, "streaming_indices_verified": keys.len(),
        "index_bytes": view.mmap.len(), "index_keys": view.hash.n(), "index_slots": view.nslots,
        "slot_bytes": view.nslots * 16, "target_bytes": view.target_count * 5,
    });
    eprintln!("trace complete: {} reads, {} seeds", limit, keys.len());
    std::fs::write(
        &args[5],
        serde_json::to_string_pretty(&json!({"trace": trace})).unwrap() + "\n",
    )
    .unwrap();
    // Free frequency analysis state before timing; the input trace stays immutable.
    drop(sequences);
    drop(batch_sequences);
    drop(frequencies);
    drop(counts);
    drop(training_frequencies);
    drop(training);
    let mut benches = serde_json::Map::new();
    for threads in [1, 12] {
        let pool = rayon::ThreadPoolBuilder::new()
            .num_threads(threads)
            .build()
            .unwrap();
        let mut row = serde_json::Map::new();
        row.insert(
            "encode".into(),
            timed(|| {
                pool.install(|| {
                    reads
                        .par_iter()
                        .map(|r| black_box(encode_acgt(&r.seq)).len as u64)
                        .sum()
                })
            }),
        );
        row.insert(
            "select_syncmers".into(),
            timed(|| {
                pool.install(|| {
                    reads
                        .par_iter()
                        .filter(|r| r.eligible)
                        .map(|r| {
                            iter_closed_syncmers_2bit(&r.encoded, view.k, view.s)
                                .map(|k| black_box(k).pos as u64 + 1)
                                .sum::<u64>()
                        })
                        .sum()
                })
            }),
        );
        row.insert(
            "map_encoded".into(),
            timed(|| {
                pool.install(|| {
                    reads
                        .par_chunks(4096)
                        .map(|chunk| {
                            let mut scratch = EcScoreScratch::new(index.num_genes());
                            chunk
                                .iter()
                                .filter(|r| r.eligible)
                                .map(|r| {
                                    let (m, stats) = index
                                        .map_encoded_sparse_resolved_with_scratch(
                                            &r.encoded,
                                            false,
                                            12,
                                            256,
                                            false,
                                            &mut scratch,
                                        );
                                    black_box(
                                        m.gene_id.unwrap_or(u32::MAX) as u64
                                            + m.score as u64
                                            + stats.gene_votes,
                                    )
                                })
                                .sum::<u64>()
                        })
                        .sum()
                })
            }),
        );
        row.insert(
            "lookup_api_preselected".into(),
            timed(|| {
                pool.install(|| {
                    keys.par_chunks(4096)
                        .map(|chunk| {
                            chunk
                                .iter()
                                .map(|&key| u64::from(black_box(index.lookup_hit(key)).is_some()))
                                .sum::<u64>()
                        })
                        .sum()
                })
            }),
        );
        row.insert(
            "hash_only".into(),
            timed(|| {
                pool.install(|| {
                    keys.par_chunks(4096)
                        .map(|chunk| {
                            chunk
                                .iter()
                                .map(|key| black_box(view.hash.index(key)) as u64)
                                .sum::<u64>()
                        })
                        .sum()
                })
            }),
        );
        row.insert(
            "lookup_scalar_preselected".into(),
            timed(|| {
                pool.install(|| {
                    keys.par_chunks(4096)
                        .map(|chunk| {
                            chunk.iter().fold(0_u64, |sum, &key| {
                                sum.wrapping_add(slot_checksum(
                                    key,
                                    view.slot_slice()[view.hash.index(&key)],
                                ))
                            })
                        })
                        .reduce(|| 0, u64::wrapping_add)
                })
            }),
        );
        row.insert(
            "lookup_stream_prefetch_preselected".into(),
            timed(|| {
                pool.install(|| {
                    keys.par_chunks(4096)
                        .map(|chunk| {
                            let slots = view.slot_slice();
                            let mut positions = [0_usize; 32];
                            let mut sum = 0_u64;
                            let mut i = 0_usize;
                            // PtrHash 2.0.2's streaming iterator supports fold/for_each, not next.
                            view.hash.index_stream::<32, _>(chunk).for_each(|next| {
                                prefetch_slot(slots, next);
                                if i >= 32 {
                                    sum = sum.wrapping_add(slot_checksum(
                                        chunk[i - 32],
                                        slots[positions[i % 32]],
                                    ));
                                }
                                positions[i % 32] = next;
                                i += 1;
                            });
                            assert_eq!(i, chunk.len());
                            for j in chunk.len().saturating_sub(32)..chunk.len() {
                                sum = sum.wrapping_add(slot_checksum(
                                    chunk[j],
                                    slots[positions[j % 32]],
                                ));
                            }
                            sum
                        })
                        .reduce(|| 0, u64::wrapping_add)
                })
            }),
        );
        assert_eq!(
            row["lookup_scalar_preselected"]["checksum"],
            row["lookup_stream_prefetch_preselected"]["checksum"]
        );
        benches.insert(threads.to_string(), json!(row));
        eprintln!("microbenchmarks complete: {threads} threads");
    }
    std::fs::write(
        &args[5],
        serde_json::to_string_pretty(&json!({"trace": trace, "microbenchmarks": benches})).unwrap()
            + "\n",
    )
    .unwrap();
}

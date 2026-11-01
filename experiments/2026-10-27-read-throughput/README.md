# Read throughput diagnostic and algorithm directions

The next substantial improvement should address both the read representation and
the number of independent index operations. The current PtrHash mapper spends
about 28% of elapsed time assembling and destroying per-read objects, outside its
previously reported mapping stage. A controlled allocator substitution recovered
17.5% throughput with identical output. For a larger algorithmic step, the best
candidate is an exact sequence-aware index that follows verified stretches with
constant gene evidence, instead of hashing every selected seed independently.

These are diagnostic results and proposed experiments. Production mapping code
and dependencies are unchanged. The diagnostic example, timing patches, scripts,
and raw measurements are retained here. The old overview documents describe p27;
the actual retained implementation at commit `3c3fd49` is PtrHash v9.

Measurements were collected on 2026-10-27 on the Ryzen 9 7900X, with CPU affinity
0–11 and 12 Rayon workers, release defaults, warm file cache, batch size 65,536,
and the first 10 million PBMC 10k v3 L001 pairs. Each controlled configuration had
one warmup and three serial measured runs. The inputs were recreated from the
retained public corpus; their hashes are in [inputs.json](raw/inputs.json).
Gzip level 1 was used for all new compressed-input comparisons, so exact timing
comparisons should use these fresh controls. All non-timing aggregate metrics
also match the August PtrHash result exactly.

Cargo-flamegraph was attempted first. The host has `perf_event_paranoid=4`, and
sudo requires a password; both sampling and hardware counters were unavailable.
The permission logs are saved in `raw/`. The new evidence comes from explicit
wall timers, an allocation audit, controlled pipeline ablations, and replay of
real queries. No new dTLB miss rate or sampled CPU percentage is claimed.

| Configuration | Median wall time | Range | Throughput | Median max RSS |
|---|---:|---:|---:|---:|
| Original binary, gzip, no output | 8.121 s | 7.942–8.131 s | 1.231M reads/s | 4.358 GiB |
| Timing instrumented, gzip, no output | 8.066 s | 7.917–8.167 s | 1.240M reads/s | 4.363 GiB |
| Same timing binary, TCMalloc | 6.863 s | 6.786–6.866 s | 1.457M reads/s | 4.346 GiB |
| Timing instrumented, uncompressed input | 7.788 s | 7.729–7.971 s | 1.284M reads/s | 4.355 GiB |

All measured runs had zero major page faults. The allocator comparison is an
ablation of allocation/ownership cost, not adoption of TCMalloc as the proposed
architecture. It improved throughput 17.5% against the identical instrumented
binary and 18.3% against the original binary. Uncompressed input improved 3.6%
against the instrumented gzip control; decompression is not the main opportunity.
The raw table and equality checks are in
[pipeline-summary.json](raw/pipeline-summary.json).

| Instrumented wall-time component | System allocator | TCMalloc |
|---|---:|---:|
| Prepare and map reads | 5.568 s | 5.686 s |
| Assemble paired batches | 0.980 s | 0.344 s |
| Destroy FASTQ batch objects | 1.019 s | 0.454 s |
| Destroy mapped read objects | 0.290 s | 0.228 s |
| Drain assignments and statistics | 0.075 s | 0.076 s |
| Wait for input batches | 0.046 s | 0.035 s |
| Load index | 0.035 s | 0.035 s |

Entries are separate medians and need not sum to median wall time. The first
three newly exposed ownership costs total 2.289 seconds, 28.4% of wall time.
Batch assembly includes dropping R1 quality vectors; batch destruction drops
the three surviving per-read sequence/quality vectors. Barcode and UMI strings
are allocated and cloned during mapping, then freed with mapped objects. The
existing input timing only measured channel waiting, so it concealed assembly
and destruction work.

A separate audit counted **324,616 score scratch initializations** for 10M reads
and 58,395 genes, with 5.768 seconds of summed CPU wall time inside constructors.
This is not 5.768 seconds of process elapsed time. Each constructor requests
about 1.168 MB of score/epoch arrays; their cumulative logical size is about
379 GB, which must not be confused with measured DRAM traffic. `map_init`
initializes per Rayon job, not once per worker, as documented by
[Rayon](https://docs.rs/rayon/latest/rayon/iter/trait.ParallelIterator.html#method.map_init).
The scratch is reused within a job but repeatedly discarded across small jobs
and batches. See [scratch-audit.log](raw/scratch-audit.log).

Full RAD output was compared byte-for-byte for the original binary and the
diagnostic binary with TCMalloc: all 10M records matched, totaling 369,461,060
bytes, SHA-256 `b7c5b5ad39e5b2bca42cf253d1d9e10ef9ac5e120f752384681a14fe4cb3b99a`.
The retained assignment rates are 78.50747% countable, 9.49510% ambiguous,
11.26884% unmapped. See [output-verification.json](raw/output-verification.json).

The first 1M reads were also traced using the real filters and sampled index.
Of these, 992,365 passed the prefilters. The diagnostic distinguishes actual
absent keys from high-frequency hits: the existing aggregate field
`selected_seeds_absent_from_index` includes both because it subtracts accepted
hits from attempted lookups.

| Trace property | Observation | Consequence |
|---|---:|---|
| Selected seed occurrences | 10,243,437 | About ten independent probes per input read |
| Accepted hits | 7,546,061 | Most probes find useful indexed evidence |
| Truly absent keys | 2,190,638, 21.4% | A negative filter can avoid only this fraction of probes |
| Hits above the gene-frequency cap | 506,738 | Preserve this distinct outcome and its diagnostics |
| Consecutive runs of identical accepted EC payload | 2,395,167 | 3.15 accepted hits per run on average |
| Gene score updates | 45,610,670 | Large gene sets make accumulation substantial |
| Gene visits after ideal per-run aggregation | 38,207,244 | Only 16.2% fewer visits; aggregation alone is limited |
| Accepted hits with one target | 66.9% | Inline singleton payloads could eliminate a dependent fetch |
| Exact eligible R2 duplicates across 1M reads | 26.8% | Whole-read memoization has potential, with a sizable working set |
| Exact eligible R2 duplicates within 65,536-read batches | 12.9% | Bounded deduplication has a smaller, measurable ceiling |

The same-EC run count is an opportunity measure, not proof that those sequences
form one unbranched reference path. Equal endpoint payloads do not authorize
skipping intervening sequence. Even an idealized scheme that replaces every
accepted run with one lookup while retaining all absent/capped queries reduces
total random probes by about 50.3%, not 68.3%; the latter percentage applies only
to accepted hits. Real extension and fallback work must be included.

The index has 259,516,542 keys, 262,137,600 slots, and a total size of
4,427,647,629 bytes. Full-key/payload slots occupy 4,194,201,600 bytes, **94.7%**
of the index. The target array is only 133,827,300 bytes. Compressing gene sets
alone therefore misses the dominant storage cost. Rebuilding p27 took 124.66 s
and 4.90 GiB peak RSS; conversion with the existing exhaustive PtrHash audit
took 60.73 s and 9.73 GiB. Input preparation overlapped parts of rebuilding, so
these build times are provenance, not a controlled build-speed comparison.

Hot-key admission looks more promising on this sampled workload than the old
dense trace suggested. A table selected using only the first 500k reads covered
the following fractions of seed occurrences in the next 500k reads:

| Admitted keys | Holdout coverage |
|---:|---:|
| 1,000 | 13.35% |
| 10,000 | 31.08% |
| 100,000 | 43.63% |
| 1,000,000 | 52.55% |

This supports a shared immutable frequency-admitted table as an experiment. It
does not validate the rejected per-worker cache: cold-key admission, writes,
and twelve copies would change the cost. The holdout is adjacent data from one
library, not evidence of portability across tissues or sequencing runs.

The isolated replay measurements below are seconds per complete 1M-read trace.
Each row has a warmup and three samples; fast cases repeat the trace within each
sample to reduce timer noise. Mapping uses pre-encoded eligible reads and
4,096-read scratch chunks. Lookup rows use already-materialized seed keys.

| Isolated operation | 1 worker | 12 workers |
|---|---:|---:|
| Encode R2 | 0.511 | 0.0447 |
| Select closed syncmers | 0.749 | 0.0681 |
| Existing encoded-read mapping and resolution | 2.760 | 0.2496 |
| Existing lookup API over preselected keys | 0.545 | 0.0490 |
| PtrHash addressing only | 0.139 | 0.0140 |
| Direct exact slot lookup, scalar | 0.324 | 0.0294 |
| Direct exact slot lookup, streamed with prefetch | 0.104 | 0.0149 |

The streamed prototype prefetches both PtrHash pilot metadata and the resulting
payload slots. It is 3.12x faster than the direct scalar kernel at one worker and
1.97x at twelve. It verifies full keys, and all 10,243,437 streaming slot indices
were checked against scalar addressing. The kernels produce matching payload
checksums. These gains exclude seed collection, scoring, and output, and the
direct view hoists validation/dispatch work compared with the generic API.
They must not be reported as mapper speedups. Timings of separate kernels are
also not an additive CPU profile: cache state, scheduling, and instruction
overlap differ. The robust conclusion is that selection/encoding are material
and the dependency chain in lookup can be overlapped. Full raw timings and
trace counts are in [trace.json](raw/trace.json).

The next experiments should be pursued in this order:

1. **Replace per-read ownership with reusable packed batches.** Parse sequences
   and qualities into contiguous buffers with offsets. Retain barcode and UMI
   as packed integers through mapping and binary output. Use inline short-read
   encoding with an explicit fallback for longer reads. Give each worker or
   bounded large chunk reusable score state, and return compact mapping records.
   Recycle the buffers as units. This removes allocations from the data model
   instead of optimizing individual allocations. Preserve read order, quality
   filtering, diagnostics, scores, and output exactly. The first gate is at
   least 20% higher end-to-end throughput against the original binary, in both
   no-output and normal RAD+zstd modes, without increased memory. Compare against
   the TCMalloc control too. The measured allocator win is evidence for this
   direction, not a measured result for the packed design.

2. **Build a sequence-aware index with constant-evidence runs.** Store a compact
   representation of reference strings/unitigs and boundaries where the complete
   gene/class payload changes. A seed lookup returns a sequence location; verify
   the subsequent read bases against that sequence and reuse the payload while
   the match and color interval continue. Fall back at mismatches, branches,
   sequence boundaries, or evidence changes. Preserve the existing selected
   seeds and their multiplicities; aggregate their nonnegative weighted votes
   with the same saturation and exon/intron/antisense tie rules. This attacks
   both the 4.19 GB slot representation and the number of independent probes.
   [SSHash](https://github.com/jermp/sshash) supplies exact associative streaming
   queries; its [Rust implementation](https://github.com/COMBINE-lab/sshash-rs)
   is a practical candidate to evaluate. The color-boundary idea is supported
   by [Themisto's key-k-mer construction](https://academic.oup.com/bioinformatics/article/39/Supplement_1/i260/7210444).
   These structures require an adapter and a new build path: they are not
   drop-in replacements. In particular, the current selector canonicalizes
   s-mer order but stores directional k-mer codes; a canonical dictionary must
   return orientation and retain the appropriate directional payloads.

   For this larger change, require at least 1.5x total throughput over the
   current baseline, an index no larger than 5 GB, mapping RSS below 6 GiB,
   construction below 32 GiB peak RSS, and byte-identical output. Count the
   complete index/build cost, including all unitig and color structures. The
   trace does not yet establish that these gates are achievable. A target near
   2M reads/s is an experiment goal, not a prediction. A rank/select-based
   [RSHash](https://github.com/jonsmcode/rshash) dictionary is another structural
   alternative worth comparing after establishing the sequence/color adapter.

3. **Try a bounded two-stage lookup pipeline as a lower-cost bridge.** Keep
   32–64 outstanding seeds in a small worker-local ring, prefetch pilots and
   payload slots, then accumulate in read order. Avoid global seed sorting or
   scattered per-read score tables. The new premise is the v9 direct-slot
   dependency chain, which did not exist in the rejected v4 seed-batching
   experiments. Include materialization, partial batches, and scoring in the
   benchmark; require at least 10% total improvement and identical output.
   This can improve latency hiding but cannot replace the algorithmic benefit
   of avoiding probes.

4. **Evaluate a small shared hot-key table, then exact read deduplication.**
   A 10k–100k-key table has a favorable holdout hit rate and a much smaller
   footprint than a million-key cache. Build it from an initial mapping sample,
   store full keys and exact hit/miss/capped outcomes, freeze it, and count its
   training/build cost in runtime. Preserve all lookup statistics and scores.
   Require at least 10% end-to-end improvement including warmup. Exact R2
   deduplication can share mapping results after each read independently passes
   quality/barcode filtering, but its 12.9% within-batch duplicate fraction is
   only a work-count ceiling. A shared barcode/UMI alone is not proof that reads
   have identical mapping evidence. Retain no cache implementation that fails
   the total-runtime gate.

5. **Vectorize selection after the representation changes.** Fuse compatible
   per-base work and evaluate a regular sliding-window minimum algorithm across
   several reads, avoiding the branch-heavy monotonic queue where profitable.
   [SimdMinimizers](https://drops.dagstuhl.de/entities/document/10.4230/LIPIcs.SEA.2025.20)
   provides relevant SIMD algorithm work, and the current library exposes
   [closed syncmers](https://docs.rs/simd-minimizers/latest/simd_minimizers/fn.closed_syncmers.html).
   Its default hasher differs from our canonical mixed-hash rule, so substituting
   the API without adapting the hasher would require a new index and a separate
   accuracy experiment. Keep seed identity, ties, strand handling, and N handling
   exact for the first prototype. A smaller exact v9 experiment is to inline the
   singleton gene/class in the existing 16-byte slot; this can remove a fetch
   for 66.9% of accepted hits, but it is a secondary optimization.

Fewer absolute-position probes, broader prefix layouts, blind early stopping,
global sorting, and negative-only filters are lower priorities given the
previous failures and this trace. Aggregate EC counting alone removes only
16.2–16.8% of gene visits. GPU offload would also leave the measured host-side
ownership bottleneck intact; it is not the first experiment here. No GPU
inspection or computation was performed.

An Amdahl calculation explains why the two architectural changes belong
together. Using the instrumented medians, even halving the entire 5.568-second
mapping stage leaves about 5.282 seconds total: only 1.53x throughput. Conversely,
removing all 2.289 seconds of exposed batch ownership work, unrealistically,
would give about 1.40x. These are optimistic arithmetic bounds, not forecasts;
lookup alone is only part of the mapping stage, and faster mapping may expose
input waiting. A substantial advance requires reducing both per-read overhead
and work within mapping.

Reproduction uses the existing workspace and local corpus. The example is built
with `cargo build --release -p constellation-core --example read_path_diagnostic`.
Its positional arguments are `INDEX R1.fastq R2.fastq READ_COUNT OUTPUT.json`.
Run it under `RAYON_NUM_THREADS=12 taskset -c 0-11`. It accepts only the current
v9 layout; the raw input fixture intentionally avoids timing decompression.

`prepare_inputs.py` recreates raw/gzip subsets and hashes. `bench_pipeline.py`
takes `LABEL BINARY system|tcmalloc gzip|raw` and records a warmup plus three
serial runs, guarded by a shared benchmark lock. `summarize.py` checks all
non-timing metrics against each other and the historical control.
`verify_output.py` compares complete RAD output. Run Python scripts with UV,
for example `UV_CACHE_DIR=/tmp/uv-cache uv run --python 3.12 SCRIPT.py`.

Large inputs, indexes, output files, and preserved baseline/diagnostic binaries
are in the git-ignored `benchdata/throughput-2026-10-27/`. The original binary is
`constellation.baseline`; `constellation.diagnostic` includes the extra scratch
counters. For the exact three-run timing instrumentation, apply
[pipeline-timing.patch](pipeline-timing.patch) to the baseline source and build;
[scratch-timing.patch](scratch-timing.patch) is an alternative, cumulative patch
for the allocation audit. Neither patch remains applied to production source.
The initial allocator series whose last repetition overlapped a raw-input
warmup was excluded and rerun in full; those excluded artifacts and an initial
PtrHash iterator API error log are preserved under the git-ignored
`excluded-runs/` subdirectory. Only the final serial runs enter the tables.

Final validation passed: all 76 workspace tests, workspace formatting checks,
Python script parsing, the full-query streaming/scalar comparison, and the
10M-read output comparison. A one-read fixture with no eligible seeds also
completed, exercising the empty-stream case. The production source and dependency
manifests have no diff. Future candidate implementations should be removed if
their stated gates fail, while preserving their measurements in the experiment
ledger.

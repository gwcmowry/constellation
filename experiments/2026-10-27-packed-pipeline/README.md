# Packed pipeline and exact read throughput experiments

Retained: the normal gzip-to-RAD+zstd pipeline is **3.00x faster**, with slightly
lower memory use, the same index, and byte-identical decoded output. The changes
are enabled automatically for normal read-at-a-time gene-EC mapping; no index
rebuild or new command-line option is required for existing v9 indexes.

| Mode | Original median | Optimized median | Speedup | Optimized throughput |
|---|---:|---:|---:|---:|
| No output | 7.955 s | 2.746 s | 2.90x | 3.641M reads/s |
| RAD with zstd level 3 | 8.949 s | 2.979 s | 3.00x | 3.357M reads/s |
| TSV, single measured run | 9.069 s | 2.947 s | 3.08x | 3.393M reads/s |

No-output and RAD+zstd results are three-run medians after a warmup. Their
optimized ranges were 2.736–2.764 s and 2.965–3.111 s, respectively. These are
10M PBMC pairs from the diagnostic corpus, warm cache, gzip level 1 inputs,
CPUs 0–11 and 12 Rayon workers; all measured runs had zero major faults.
The original binary was retained from commit `3c3fd49`; both binaries used
the same inputs and the same 4,427,647,629-byte PtrHash index. Zstd-output
peak RSS fell from a median 4.372 GiB to 4.265 GiB. Process CPU time, including
the compressor, fell from 65.73 to 27.00 seconds. Headline timings use the
application's end-to-end wall timer, including index load and output completion;
the full `/usr/bin/time -v` reports are retained too.

The normal pipeline now has three recyclable contiguous input slabs per mate,
one persistent score/encoding/resolution state per Rayon worker, and ordered
output chunks independent of worker scheduling. On this run the score scratch
is initialized 12 times, instead of the diagnostic's 324,616 times. A safe,
normally uncontended mutex associates each state with its Rayon worker; it is
acquired once per chunk, with no nested parallel work while held. Input errors
and early cancellation disconnect both channels before joining reader threads.

Mapping produces statistics and output bytes directly, eliminating owned
per-read assignments, barcode/UMI string clones, and intermediate mapped-read
vectors. ASCII barcode/UMI fields are borrowed and encoded on output; invalid
UTF-8 boundaries retain the prior lossy-string behavior. Gene-tie lists and
encoded sequence buffers retain capacity across reads, including changing read
lengths and N masks. Filtering and encoding share one lookup-table-based pass;
U and invalid-base handling exactly match the old encoding and complexity rules.

V9 lookup gathers at most 32 seeds from one read at a time, pipelines PtrHash
pilot reads eight queries ahead, prefetches the resulting slots, then checks
full keys and scores in the original order. This is bounded memory and does not
sort or mix reads. Empty streams, partial groups, long reads spanning multiple
groups, reverse complements, and saturated scores are checked against the
independent source-v4 scalar mapper. Other index versions retain their existing
lookup paths.

Once those changes made mapping faster, gzip decoding became the bottleneck.
Enabling flate2's `zlib-rs` backend removed most of that wait. The dependency
adds one locked Rust package and keeps streaming concatenated-gzip support and
CRC/truncation validation. The backend choice follows the supported
[flate2 configuration](https://github.com/rust-lang/flate2-rs#backends); the
measurement here, rather than an upstream benchmark, establishes its benefit
for this workload. No allocator replacement or CPU-specific compiler flags are
needed by the retained implementation.

The intermediate measurements were used to select implementations, not as the
headline benchmark. Each row below is one measured no-output run after warmup:

| Incremental candidate | Total | Preparation and mapping | Input wait |
|---|---:|---:|---:|
| Reusable slabs, scratch and fused output | 5.250 s | 5.123 s | 0.067 s |
| Add fused encoding and complexity check | 4.595 s | 4.193 s | 0.341 s |
| Add bounded v9 prefetch | 4.196 s | 2.704 s | 1.435 s |
| Add zlib-rs input decoding | 2.852 s | 2.768 s | 0.028 s |
| Give each worker one persistent state | 2.644 s | 2.422 s | 0.177 s |

An additional exact aggregation of repeated EC payloads took 2.888 s versus
2.852 s for its immediate control. It showed no useful total-time gain and was
removed. Its code is preserved in `rejected-run-aggregation.patch`; its raw runs
have the `runs-` prefix. The smaller prototype results remain single-run
diagnostics; the final configuration was then retested in the controlled series
above. A new sequence/unitig index and a learned hot-key table were not needed
to obtain the retained speedup and are not part of this change.

All non-timing metrics match the baseline in every benchmark run. Complete
10M-read RAD output and decompressed zstd output match the original 369,461,060
bytes (SHA-256 `b7c5b5ad39e5b2bca42cf253d1d9e10ef9ac5e120f752384681a14fe4cb3b99a`).
Complete TSV output also matches: 674,434,328 bytes, SHA-256
`2e2a622111f13194c9019103af8eef005fb8820ee4706a18a8b5c79c254ab702`.
Both compressed files are 168,459,890 bytes. Countable, ambiguous and unmapped
rates remain 78.50747%, 9.49510% and 11.26884%. See
[summary.json](raw/summary.json) and [output-equality.json](raw/output-equality.json).

The remaining dominant stage is the fused preparation, seed selection, lookup,
scoring and record-formatting work: 2.687 s of the 2.979 s RAD+zstd run, with
0.219 s in output and 0.014 s in input waiting (separate medians). Those worker
timings include output construction, so they should not be interpreted as
lookup-only CPU time. The former roughly 2.3 seconds of batch assembly and
per-read object destruction is gone. Diagnostics, locality sorting, and older
experimental traversals deliberately use the original richer object path;
their throughput is not represented by the fast-path benchmark.

Validation adds differential CLI coverage for TSV/RAD, dense/v4/v9 indexes,
reverse complements and trimming, short/invalid/U/Unicode barcodes, N/lowercase
DNA, empty and long reads, and quality filtering. Unit tests cover input slab
recycling/cancellation, concatenated and corrupt/truncated gzip input, fused
encoding with changing lengths and arbitrary byte values, and independent
v4-versus-v9 score/statistics/tie-list comparisons with reusable output and long
repetitive reads. The workspace suite has 79 tests. Formatting and Clippy are
also checked; the remaining Clippy warnings are in pre-existing code.

To reproduce, build with `CARGO_BUILD_JOBS=12 RAYON_NUM_THREADS=12 taskset -c 0-11
cargo build --release --workspace`, preserving the original and optimized
binaries in the diagnostic data directory as `constellation.baseline` and
`constellation.optimized`. Run `final_suite.py` with UV and then `summarize.py`.
The suite serializes benchmark processes, runs fresh warmups, checks all
non-timing metrics, and compares complete outputs. Raw logs are in `raw/`;
large binaries and outputs stay in git-ignored
`benchdata/throughput-2026-10-27/`. The fixture creation, source hashes and
earlier diagnostic are in `../2026-10-27-read-throughput/`.

The following paragraphs record the predeclared experiment scope and gates.

The experiment replaces the default gene-EC per-read ownership pipeline with
bounded reusable FASTQ slabs, persistent worker scratch/encoding/resolution
storage, and fused mapping, statistics, and output formatting. Scores, filters,
gene/class tie ordering, diagnostics, and byte-level output are the oracle.
The existing path remains available for diagnostics and historical traversal
options. No new index or change in evidence selection is required.

The initial acceptance gate is at least 20% greater end-to-end throughput than
the retained PtrHash binary, at no greater mapping RSS, in both no-output and
RAD+zstd modes. Run one warmup and three serial warm-cache repetitions on the
same 10M PBMC pairs, batch size 65,536, CPUs 0–11, and 12 Rayon workers. Compare
all non-timing metrics and complete RAD bytes; test TSV, malformed and truncated
input, invalid and short barcodes, U/N/lowercase bases, read length changes,
reverse complements, trimming, dense and sampled indexes, and recycled batches.

After this gate, separately measure bounded prefetched lookup, fused per-base
preparation, and exact run aggregation. The new premise for prefetching is the
PtrHash direct-slot dependency chain, with read-local scoring retained; this
does not repeat global seed sorting/prefix tiling. Require a material total-time
win, including materialization and output, before retention. Remove failed
production experiments while preserving their measurements here.

Hardware sampling and counters were attempted in the preceding diagnostic and
are still unavailable (`perf_event_paranoid=4`, sudo requires a password).
Explicit stage timings and controlled ablations are retained as raw JSON/time
files under `raw/`; no new hardware-counter percentages are claimed.

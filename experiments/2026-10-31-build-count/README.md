# Faster index construction and state counting; explicit USA assignment policy

The build and counting changes preserve the conservative output while reducing
wall time. A separate optional USA resolver increases assigned molecules, with a
measured source-origin precision tradeoff. It remains opt-in and does not recover
sequence information absent from the sampled index.

## Measured result

Three serial repetitions per stage, same ten-million-pair workload:

| Stage | Previous median | New median | Speedup |
|---|---:|---:|---:|
| Build usable v9 index, including old conversion | 185.047s | 24.407s | 7.58× |
| FASTQ → gene/state assignments | 4.782s | 4.385s | 1.09× |
| Assignments → all count matrices | 7.881s | 1.538s | 5.13× |

A separate complete chained run took **196.02s → 30.76s (6.37× faster)**,
including index creation on that run. The sums of stage medians are 197.71s and
30.33s; those are a separate statistic. With a prebuilt index, the corresponding
map-plus-count median sums are 12.66s and 5.92s (2.14×). Mapping's algorithm is
unchanged; the small observed mapping difference is not an intended optimization
and may reflect build/layout/compiler/cache variation.

The three optimized builds range from 24.37–24.57s; the old build-plus-conversion
runs range from 184.70–185.81s. Build peak RSS rises from a median **9.84 GiB to
14.35 GiB**. Counter RSS falls from **1.17 GiB to 0.48 GiB**. Mapping RSS remains
about 4.27 GiB. Final index size is unchanged at **4,427,647,629 bytes**.

All 259,516,542 retained-key payloads match the baseline; one million further
probes agree. The ten-million-read CSTRAD2 output and every conservative matrix
and axis file are byte-identical, including the separate chained runs. All
**91 workspace tests** pass. Index binary hashes differ; semantic index identity,
rather than artifact byte identity, is the validation criterion.

Raw measurements: [summary.json](raw/summary.json), [validation.json](raw/validation.json),
[index audit](raw/index-validation.json), and [tests](raw/tests.log).

## Workload and boundaries

Same prepared Ensembl 93 splici FASTA/t2g, k31/s20 closed syncmers, 58,395 genes,
259,516,542 retained keys; ten million original PBMC read pairs and the published
11,769-cell barcode whitelist. Cell calling is not performed; the whitelist is
supplied to both versions. Same machine, CPUs 0–11, 12 Rayon workers. CPU
benchmarks run serially. No GPU or additional hardware is used.

Build timing begins with the prepared reference FASTA and ends with a usable v9
index. The old path includes BOTH the v4 build and conversion to v9. The new path
builds v9 directly. Reference FASTA generation from genome/GTF, downloads and
compilation are outside these timings. Mapping includes index loading, compressed
FASTQ input, state-preserving CSTRAD2 output and teardown. Counting includes
barcode correction, resolution, directional UMI correction and all seven sparse
matrices. Timings include process startup/teardown, with normal buffered output;
they are not cold-cache or fsync-durability measurements.

Stage timings use three repetitions. The report separately records their median
sums and a complete serial build → map → count execution for each version. The
latter includes index construction on that invocation; it is not a fused command
or a prebuilt-index run. The speed table uses the unchanged conservative policy;
native-compatible assignment changes are evaluated separately.

## Why the build is faster

The baseline instrumented run spends 46.13 seconds generating/sorting/writing
runs, 75.84 seconds merging and interning ECs, and about 3.5 seconds writing v4.
Conversion then costs about 60 seconds. This made the serial external merge the
first target.

The partitioned algorithm:

1. Reads bounded batches of reference targets and extracts seeds in parallel.
   Assigns every full key to a high-bit partition; writes compact temporary
   records without a per-transcript hash-set insertion for every seed.
2. Sorts partitions independently, combines duplicate key/gene records with
   exact class-bit union, and interns ECs with a hash table. Key ranges do not
   overlap, so there is no serial heap merge of all occurrences.
3. Combines the partition-local EC dictionaries in key order, preserving the
   old gene IDs, EC ordering and annotation semantics.
4. Constructs PtrHash directly from the retained keys, without writing and
   decoding a v4 prefix index.
5. Partitions slot references into disjoint ranges before filling the v9 table
   in parallel. Every collision is checked. A serial scattered fill took about
   18 seconds in the first direct-build prototype; partition plus parallel fill
   takes about 4 seconds. This changes construction, not the query layout.
6. Verifies every retained key/payload before writing the index, plus bounded
   random absent-key checks. The standalone regression audit additionally checks
   the final index against the previous builder.

New construction uses roughly 5.45 GB of temporary partition data, removed as
partitions are consumed, and more peak memory than the old path. The default CLI
builder is partitioned; `--build-mode streaming` retains the previous algorithm.
The old lower-peak path remains available by combining streaming construction
with v4 output and subsequent conversion. The direct writer itself retains its
input tables while constructing the lookup table. The final index format and
size are unchanged. `--format ptrhash` is explicit;
`--format mmap` and legacy formats remain available.

## Why counting is faster

The old counter spent about 2.39 seconds ingesting/global hashing, 1.91 seconds
resolving molecules, 2.44 seconds correcting UMIs/building sparse trees, and 0.92
seconds writing. It allocated one candidate vector for every distinct CB/UMI.

The new counter keeps barcodes/UMIs packed, caches barcode correction by integer
code, and appends evidence to contiguous buffers per corrected cell. Cells are
independent: each sorts its own UMIs, resolves genes/states and corrects UMIs in
Rayon. Counts are gathered once into a flat sorted array, instead of repeatedly
inserting millions of entries into B-trees. Original barcode order and the old
UMI ordering/tie rules are preserved, so conservative outputs remain identical.

The first new count run spends about 0.88 seconds ingesting, 0.16 seconds
resolving/correcting, 0.09 seconds gathering/sorting, and 0.43 seconds writing.
Input parsing/collation and output serialization now dominate. A future fused
map/count pipeline could avoid the intermediate assignment file; that is not
implemented here.

## Assignment-rate result

The new explicit mode is:

```sh
constellation count-splicing --assignments reads.cstrad2 \
  --barcode-whitelist barcodes.txt --resolution usa --no-correct-umis \
  --out-prefix counts --emit-metrics counts.metrics.json
```

It applies the exact-UMI `cr-like` / `WinnerTakeAll` USA rule tested against
alevin-fry 0.19.0: maximum-supported gene/class labels, followed by the native
spliced preference in cross-gene ties of at most ten labels. A votes for both S
and U. Ties between S and U of one gene remain A. This currently requires exact
UMIs and does not implement native UMI edit correction. Groups containing our
C/X annotations fall back to conservative intersection and are counted in metrics.

| Same Constellation read evidence, exact UMIs | Intersection | USA |
|---|---:|---:|
| Spliced | 2,661,927 | 2,823,561 |
| Unspliced | 2,011,089 | 2,011,271 |
| Ambiguous, including C/X | 852,738 | 894,749 |
| Total | 5,525,754 | 5,729,581 |

Net increase: **203,827 molecules / 3.69%**, under matched exact-UMI handling.
There are 203,836 newly assigned groups and nine previously assigned groups
lost under USA rules. The policy uses spliced tie-breaking in 203,516 groups;
1,548 C/X groups fall back to intersection. No previously uniquely resolved gene
is changed in this input; 228 retained groups change their state. These counts
should not be described as 203,827 newly proven correct molecules.

With piscem evidence, the new mode produces S=2,959,526, U=1,994,063, A=841,163,
total=5,794,752, matching native alevin-fry. The validation compares **every
nonzero matrix coordinate and value**, not only totals. The reference,
orientation, unique Hamming-1 barcode correction, whitelist and exact-UMI
handling are matched; see the preceding state-comparison experiment for the
lossless target-set-to-CSTRAD2 adapter.

For comparison, our ordinary directional UMI correction plus intersection gives
5,514,020 molecules. Comparing USA directly with that total also includes the
effect of disabling UMI edit correction; it is not the clean policy-only gain.

## Accuracy tradeoff

The frozen USA tie rule was evaluated on the existing 91bp and 150bp source-origin
holdouts. Tagged UMIs make each read a distinct synthetic molecule. These are
360,000 positive reads per length plus 30,000 intergenic and 30,000 random-DNA
controls, using the previously saved mappings. This isolates the resolver's
single-read behavior; real multi-read molecules have different evidence.

| Length | Policy | Unique source-gene precision | Correct unique recall |
|---|---|---:|---:|
| 91bp | Intersection | 99.945% | 85.181% |
| 91bp | USA | 99.657% | 87.304% |
| 150bp | Intersection | 99.984% | 88.977% |
| 150bp | USA | 99.768% | 90.725% |

The additional 91bp calls contain 7,642 correct and 913 wrong source-gene calls;
909 of the latter are intronic origins. At 150bp there are 6,293 correct and 710
wrong additional calls, all 710 intronic. Shared sequence can be fundamentally
ambiguous: an annotation preference raises the assignment rate without supplying
more evidence. Random-DNA controls remain unassigned. Raw results also report
unique calls on intergenic controls; these are different from any-candidate
assignment rates in the prior mapper report.

Therefore the default remains conservative. Improving recall while preserving
its precision requires targeted sequence/path verification, especially for weak
short reads and shared S/U evidence. The prior full-reference audit showed that
most excess state ambiguity comes from unsampled discriminating k-mers. The new
build and counting algorithms do not change that sampling rule or silently
force those states.

## Validation, provenance and reproduction

- Full workspace tests cover k=1,3,15,31,32; dense/sampled builds; class union,
  repeats, lowercase/ambiguous bases, empty indexes; state resolution, packed
  barcode/UMI ordering, compression and malformed CSTRAD2.
- `verify_index.rs` checks every retained key/payload in the full human index,
  gene dictionaries, occupied-slot counts and one million deterministic probes.
  Binary artifact bytes are allowed to differ; query evidence must agree.
- `verify.py` compares all seven conservative matrices and their axes byte for
  byte, the full mapping output, and all native USA matrix entries.
- Hardware sampling was attempted first with cargo-flamegraph, then with sudo
  through an approved elevated call. Host `perf_event_paranoid=4` and unavailable
  noninteractive sudo prevent sampling. The saved errors are not profiles;
  phase timers and process CPU/RSS statistics supply the diagnostic measurements.
- Raw commands, timings, phase dumps, validation and quality results are in
  `raw/`. The original executable and an instrumented baseline are preserved in
  ignored `benchdata/build-count-2026-10-31/`. Prior unrelated changes in the
  workspace are retained. No dependency was added.

Use the existing UV-managed NumPy/PyArrow environment at
`benchdata/comparison-2026-10-27/venv/bin/python`. `benchmark.py` records serial
pinned commands; `baseline`, `optimized` and `repeat-build` reproduce the saved
suites (the original baseline binary must be retained). `check_quality.py`
reproduces the policy audit. The `pipeline` suite measures complete chained runs. The exact standalone Rust
audit build command is saved in `raw/verify-index-build.json`. Redundant index copies from repeated
builds are removed by `finalize.py` after validation; the baseline final index,
final optimized index, binaries, count outputs and raw logs are retained. The
cleanup manifest records exactly which experiment-created files were removed.

Recommended construction, once `cargo build --release` is complete:

```sh
RAYON_NUM_THREADS=12 taskset -c 0-11 target/release/constellation index-ec \
  --transcripts reference.fa --t2g-map t2g.tsv --k 31 --syncmer-s 20 \
  --format ptrhash --out reference.ptrhash.ecidx
```

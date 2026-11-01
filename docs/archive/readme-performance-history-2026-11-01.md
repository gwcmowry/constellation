# Historical README performance notes

Archived on 2026-11-01 when the README was reorganized around the current
CLI. This preserves older measurements and descriptions; statements about the
“current” index, missing counting support, and host restrictions below belong
to their original experiments and are superseded by the current README.

### Earlier benchmarks and implementation history

The [build-to-matrix experiment](../../experiments/2026-10-31-build-count/README.md)
measures index construction, mapping and counting together. Direct partitioned
index construction and cell-wise counting preserve the conservative outputs.
The optional USA resolver increases assignment rate with an explicit spliced
preference; its source-origin accuracy tradeoff is measured separately.

For the measured comparison against piscem 0.24.0 and kallisto 0.52.0, including
accuracy, negative controls, two read lengths, and two PBMC lanes, see the
[2026-10-30 comparison](../../experiments/2026-10-30-consensus/README.md).
The new opt-in `map --gene-ec-consensus` mode intersects the gene sets of at least
two distinct matched seeds and preserves cross-gene ambiguity. It substantially
improves specificity, at a cost in speed and short-read recall. It uses existing
indexes; its score is the distinct supporting-seed count. It supports forward
read-at-a-time mapping, with the ordinary packed pipeline. The default weighted
policy and its historical outputs remain available.

Controlled 10M PBMC comparison on 2026-10-27, three serial warm-cache runs,
12 physical cores and 12 Rayon workers, with the same existing 4.428GB PtrHash
index:

| Mode | Original wall time | Current wall time | Speedup | Current throughput |
|---|---:|---:|---:|---:|
| No output | 7.955 s | 2.746 s | 2.90x | 3.641M reads/s |
| RAD + zstd level 3 | 8.949 s | 2.979 s | 3.00x | 3.357M reads/s |

Complete RAD and TSV outputs are byte-identical to the original. Countable,
ambiguous and unmapped rates are unchanged at 78.50747%, 9.49510% and 11.26884%.
Median max RSS with zstd output fell from 4.372 to 4.265 GiB. Reusable input
slabs, persistent worker state, fused preparation/output construction, bounded
v9 prefetching and faster gzip decoding account for the gain. No index rebuild
is needed for an existing v9 index. Diagnostics and historical traversal modes
retain the richer object pipeline. Raw measurements, validation and limitations
are in [the packed pipeline experiment](../../experiments/2026-10-27-packed-pipeline/README.md).

Historical sampled-versus-dense comparison (2026-08-02, three serial runs):

```text
path             index       median throughput   countable   ambiguous   unmapped   max RSS
dense v4 p24     17.892 GB   1.048994M reads/s   78.53896%   9.26827%    11.46418%  16.823 GiB
syncmer s20 p27   3.538 GB   1.068055M reads/s   78.50747%   9.49510%    11.26884%   3.964 GiB
```

The retained sampled index is 5.06x smaller, uses 4.24x less mapping RSS, and was
1.82% faster by median throughput. Its countable rate differs by -0.03149
percentage points; unmapped improves by 0.19534 points while multi-gene ambiguity
increases by 0.22683 points. The workload is 10M PBMC 10k v3 read pairs on a
Ryzen 9 7900X, release build, warm cache, and a hard cap of 12 threads/cores.

The direct resolver scans reusable epoch-tagged scores without allocating and
sorting an owned score vector for every read. Differential tests cover the old
and new resolution paths; aggregate dense assignments are unchanged. The
sampled index is self-describing, uses a canonical mixed-hash closed-syncmer
selector, and ignores the old absolute query stride so sampling is not compounded.

The older controlled zstd-output result from the dense index was:

```text
zstd-output internal wall         10.13 s
zstd-output throughput            0.987M reads/s
zstd-output mapping               7.09 s
assignment/record construction    0.21 s
zstd/file write                   0.15 s
compressed output size            160M
```

Raw sampled-index runs are archived in
`experiments/archive/2026-08-02-syncmer-index/`; the preceding locality audit is
in `experiments/archive/2026-08-02-memory-locality-audit/`.

Historical hardware-counter profile from the same v4 lookup family:

```text
task-clock                        105.18 s, 11.1 CPUs utilized
instructions                      513.3B
cycles                            551.6B
IPC                               0.93
L1D load miss rate                2.22%
branch miss rate                  4.70%
dTLB load miss rate               32.78%
```

Function-specific perf samples:

```text
LoadedEcIndex::lookup             34.2% cycles, 55.6% dTLB-load-misses
encode_acgt                       9.1% cycles, 28.0% branch-misses
sketch_read                       7.0% cycles, 17.7% branch-misses
accumulate_ec_seq                 6.3% cycles, 9.7% dTLB-load-misses
KmerIter::next                    5.7% cycles
parse_tenx_3p_v3_r1               2.9% cycles
gzip inflate                      2.6% cycles
sort/dedup                        2.5% cycles
zstd compression                  0.9% cycles
```

The remaining dominant stage is fused preparation, selection, lookup, scoring
and output formatting. New hardware-counter profiling remains unavailable, so
the historical counters above are not measurements of the current pipeline.
Global seed batching, prefix tiling, larger v5-v7 layouts, higher absolute probe
strides, per-worker lookup caching, and mmap huge-page advice failed the combined
speed/quality/memory bar. The retained v9 prefetching is bounded within one read
and does not sort or regroup reads; see the failure ledger before revisiting
older approaches.

## Constellation vs simpleaf/alevin-fry

100M read-pair run on the same PBMC 10k v3 L001 subset:

```text
tool/path                         wall time   runtime split                                  max RSS
Constellation gene-EC TSV         1:56.86     116.19s internal; 72.48s map; 22.71s assign    25.3G
Constellation gene-EC RAD-like    1:58.01     117.25s internal; 73.54s map; 22.51s assign    25.3G
simpleaf full quant pipeline      2:24.20     134.77s map; 1.23s GPL; 3.89s collate; 4.30s quant 3.93G
```

On this benchmark, Constellation is faster wall-clock than the full simpleaf pipeline, but it uses much more memory. The current prefix24 mmap EC index keeps an explicit posting structure resident through the OS page cache and reached about `25.3G` RSS, while simpleaf/alevin-fry was about `3.93G`.

The comparison is not yet fully output-equivalent. Constellation currently emits per-read gene/EC assignments, while simpleaf produces RAD plus downstream UMI-resolution/count outputs. The 100M aggregate gene-total comparison against simpleaf was about `0.94` Pearson/Spearman on log common genes; a tighter comparison still needs a molecule-level output/counting path that matches alevin-fry semantics more closely.

Output size on the 100M run:

```text
TSV                               6.44 GiB
TSV + zstd -3                     1.81 GiB
RAD-like binary                   3.13 GiB
RAD-like + zstd -3                1.39 GiB
RAD-like + zstd -10               1.20 GiB
```

Zstd is cheap and useful, but it does not make per-read output 10x smaller. Getting there likely requires molecule/EC aggregation instead of one variable-length record per read.

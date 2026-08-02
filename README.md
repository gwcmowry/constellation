# Constellation

`constellation` is a Rust scRNA-seq mapper/counting prototype focused on high-throughput 10x-style gene assignment. The current best path is the gene-EC mapper backed by a prefix24 mmap EC index. Superseded plans are archived under `docs/archive/`; rejected approaches and their revisit gates are recorded in `docs/failed_performance_experiments.md`.

## Current Path

Build a splici-style gene-EC index:

```bash
cargo run --release -p constellation-cli -- index-ec \
  --transcripts /path/to/splici_or_gene_ec_reference.fa \
  --t2g-map /path/to/t2g.tsv \
  --k 31 \
  --format mmap \
  --out /tmp/reference.prefix24.mmap.ecidx
```

Map paired 10x FASTQs:

```bash
cargo run --release -p constellation-cli -- map \
  --index /tmp/reference.prefix24.mmap.ecidx \
  --r1 /path/to/R1.fastq.gz \
  --r2 /path/to/R2.fastq.gz \
  --batch-size 65536 \
  --output-format gene-ec-rad \
  --output-compression zstd \
  --zstd-level 3 \
  --emit-metrics /tmp/constellation.metrics.json \
  --out /tmp/constellation.cstrad.zst
```

The RAD-like output is Constellation-specific, not currently alevin-fry-compatible RAD. It stores per-read barcode/UMI, primary gene, assignment type, score, flags, and top tied gene IDs.

## Current Performance

Latest controlled configuration (2026-08-02):

```text
index                 17GB human splici v4 prefix24 mmap EC index
input                 10M PBMC 10k v3 read pairs
mapper                gene-EC prefix24 mmap path
output                Constellation RAD-like binary, optionally zstd level 3
host                  Ryzen 9 7900X; release build; hard cap of 12 threads/cores
```

10M read-pair warm-cache profile:

```text
no-output internal wall           9.50-9.57 s
no-output throughput              1.044-1.052M reads/s
zstd-output internal wall         10.13 s
zstd-output throughput            0.987M reads/s
zstd-output mapping               7.09 s
assignment/record construction    0.21 s
zstd/file write                   0.15 s
unique gene rate                  78.54%
ambiguous gene rate               9.27%
unmapped rate                     11.46%
compressed output size            160M
max RSS                           16.8 GiB
```

The accepted dense epoch-score accumulator and cheaper direct preprocessing are
byte-identical to the previous compressed output and improve the controlled
no-output path by roughly 1–2%. Raw runs are archived in
`experiments/archive/2026-08-02-memory-locality-audit/`.

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

The main remaining bottleneck is EC index lookup locality and address translation pressure. Output compression is not currently a material runtime cost. Straight seed batching, prefix tiling, larger v5-v7 layouts, higher absolute probe strides, per-worker lookup caching, and mmap huge-page advice all failed the combined speed/quality/memory bar; do not retry them without reading the failure ledger.

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

## Development

```bash
taskset -c 0-11 cargo fmt --all
CARGO_BUILD_JOBS=12 RAYON_NUM_THREADS=12 taskset -c 0-11 cargo test --workspace
CARGO_BUILD_JOBS=12 RAYON_NUM_THREADS=12 taskset -c 0-11 cargo build --release --workspace
```

Useful profiling commands:

```bash
perf stat -d -d -d -- target/release/constellation map ...
perf record -F 999 -g -- target/release/constellation map ...
perf report --stdio --no-children --sort comm,dso,symbol
```

On the current host, hardware-counter collection is blocked by
`perf_event_paranoid=4`; the commands above require the host policy to change.
Keep all mapper benchmarks at or below 12 threads.

Generated benchmark data, large references, and local profiler outputs should stay outside git.

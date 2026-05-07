# Constellation

`constellation` is an experimental Rust RNA-seq/scRNA-seq mapper focused on cache-aware batching. Instead of processing each read independently, it sketches reads into similarity buckets, performs seed lookup, regroups candidate hits by transcript/reference locus, and scores many read-candidate pairs against cache-resident sequence windows. The project uses scalar reference implementations for correctness and `pulp`-based SIMD kernels for regular scoring filters. The first target is transcriptome-level 10x-style gene-expression assignment with synthetic and public PBMC benchmarks.

This repository was initialized from the original `cachebatch` project brief, with the project and CLI renamed to `constellation`.

## Current MVP

- Cargo workspace with `constellation-core` and `constellation-cli`.
- 2-bit DNA encoding, decoding, reverse complement, k-mer iteration, and tests.
- Transcript FASTA parser and in-memory k-mer index builder.
- Compact mmap-backed `.cbidx` serialization, with JSON loading retained for prototype compatibility.
- Streaming compact-index builder for large FASTA targets that avoid holding all postings in memory.
- `constellation index`, `inspect-index`, `simulate`, `map`, and `bench-report` commands.
- Deterministic synthetic paired FASTQ generation.
- Read-at-a-time, sketch-bucket, and candidate-locus mapping modes for cache-locality comparisons.
- Gzip and plain FASTQ input support.
- `needletail` FASTA parsing, `rayon` sorting, and mmap-backed index loading.
- Optional `--gtf` transcript-to-gene mapping during indexing.
- GTF/genome-derived target generation for exon transcripts, gene bodies, introns-only, and exon-plus-gene-body targets.
- Low-quality and low-complexity read classification.
- Rare-anchor seed selection that prefers lower-frequency k-mers before candidate lookup.
- Gene-level k-mer document frequency in compact v3+ indexes, stored as a compact side array for `gene-idf` seed planning.
- Compact v4 indexes class-partition each k-mer posting list by exon, intron, gene-body, then unknown target class while preserving v2/v3 read compatibility.
- Seed-batched retrieval mode that groups query seed occurrences by k-mer and loads each selected posting list once per batch.
- Sparse-probe candidate search with fallback for faster average-case lookup.
- Transitional exon-first candidate search that probes exon postings first from v4 class-partitioned indexes, then falls back to the normal sparse/full path when exon evidence is insufficient.
- Streamed seed-batched mapping chunks so candidate generation, scoring, and assignment run on bounded working sets.
- Gene-level WAND-lite candidate pruning with configurable score ratio and max loci per retained gene.
- Reverse-complement seed lookup generated on the fly from encoded query k-mers, without building a larger reverse-complement index.
- Target-class-aware scoring and assignment for exon transcript, intron, and gene-body hits.
- Conservative exon-priority assignment that resolves same-gene exon/gene-body ties without resolving multi-gene ties.
- Optional failed-read rescue pass that reruns full unpruned seed-batched candidate search for unmapped or antisense-only reads and only replaces the call when it becomes gene-countable.
- Hot candidate-search tables use 16-byte k-mer entries and 8-byte postings so reference slices can fit more cleanly into L1/L2 cache.
- Scalar scoring plus a `pulp` feature-gated SIMD Hamming path with scalar-vs-SIMD differential tests.
- Configurable mismatch, poly-A/poly-T trimming, TSO trimming, low-quality tail trimming, and softclip fallback scoring.
- Assignment TSV output, JSON metrics, and truth-aware benchmark reporting for simulated reads.
- Streaming `count` output as Matrix Market plus barcode/features TSVs, with optional 10x whitelist barcode correction, feature filtering, conservative molecule-gene resolution, and directional UMI correction.
- Split assignment metrics for unique gene, same-gene multi-transcript, multi-gene ambiguous, antisense-gene diagnostics, and gene-countable rates.
- Per-read unmapped diagnostics and aggregate score-failure decomposition.

## Current Performance Snapshot

Latest fast path: gene-EC mapping against the Ensembl 93 splici-style prefix24 mmap EC index:

```text
/tmp/human_ensembl93_splici_r91.k31.t2g.prefix24.mmap.ecidx
```

The current recommended mapper configuration is `--mode gene-ec`, `--output-format gene-ec-rad`, optionally `--output-compression zstd --zstd-level 3`. The gene-EC path uses the prefix24 mmap lookup table, sparse read probing, streamed paired-FASTQ decode, RAD-like binary records, and optional streaming zstd output.

Recent Constellation-only profile on 10M PBMC 10k v3 read pairs:

```text
input                           /tmp/pbmc_10k_v3_10M_L001_R{1,2}.fastq.gz
output                          gene-EC RAD-like + zstd level 3
perf stat elapsed               9.45 s
internal wall time              8.69 s
throughput                      1.15M reads/s
mapping/candidate generation    4.16 s  (47.9%)
assignment/record construction  2.28 s  (26.2%)
zstd/file write                 0.20 s  (2.3%)
FASTQ batch wait                0.04 s  (hidden by pipeline)
unique gene rate                78.54%
ambiguous gene rate             9.27%
unmapped rate                   11.46%
compressed output size          160M
```

Perf hardware-counter profile from the same run:

```text
task-clock                      105.18 s, 11.1 CPUs utilized
instructions                    513.3B
cycles                          551.6B
IPC                             0.93
L1D load miss rate              2.22%
branch miss rate                4.70%
dTLB load miss rate             32.78%
```

Function-specific perf samples show the remaining mapper bottleneck is EC index locality, not output compression:

```text
LoadedEcIndex::lookup           34.2% cycles, 55.6% dTLB-load-misses
encode_acgt                     9.1% cycles, 28.0% branch-misses
sketch_read                     7.0% cycles, 17.7% branch-misses
accumulate_ec_seq               6.3% cycles, 9.7% dTLB-load-misses
KmerIter::next                  5.7% cycles
parse_tenx_3p_v3_r1             2.9% cycles
gzip inflate                    2.6% cycles
sort/dedup                      2.5% cycles
zstd compression                0.9% cycles
```

100M read-pair comparison against `simpleaf quant`/alevin-fry on the same L001 subset:

```text
tool/path                       wall time   mapper/runtime split                       max RSS
Constellation gene-EC TSV       1:56.86     116.19s internal; 72.48s map; 22.71s assign 25.3G
Constellation gene-EC RAD-like  1:58.01     117.25s internal; 73.54s map; 22.51s assign 25.3G
simpleaf full quant pipeline    2:24.20     134.77s map; 1.23s GPL; 3.89s collate; 4.30s quant 3.93G
```

On this benchmark, Constellation is faster wall-clock than the full simpleaf pipeline but uses much more resident memory because the current EC index is a large explicit mmap posting structure. The comparison is not yet output-equivalent: Constellation emits per-read gene/EC assignments, while simpleaf produces RAD plus downstream UMI-resolution/count outputs. Aggregate gene-total correlation against simpleaf on the 100M run was about `0.94` Pearson/Spearman on log common genes, but the exact UMI/count pipeline still needs a closer apples-to-apples comparison.

Output size status:

```text
100M TSV                        6.44 GiB
100M TSV + zstd -3              1.81 GiB
100M RAD-like binary            3.13 GiB
100M RAD-like + zstd -3         1.39 GiB
100M RAD-like + zstd -10        1.20 GiB
```

Zstd is useful and cheap, but it is not enough for a 10x reduction on per-read records. Getting to that level likely requires molecule/EC aggregation rather than one variable-length record per read.

Historical compact-index benchmark: 100k read pairs from 10x PBMC 10k v3, lane L001 subset, mapped with release build against the Ensembl 93 exon-plus-gene-body compact v3 index:

```text
benchdata/reference/ensembl93/human_ensembl93_exon_plus_gene_body.k21.genedf.cbidx
```

The current fast path is seed-batched, sparse-probe candidate search, streamed chunks, gene-level WAND-lite pruning, reverse-complement search, TSO/polyA/polyT/low-quality trimming, and softclip scoring. Use the internal metrics JSON wall-clock fields for timing; process-level `/usr/bin/time` includes mmap/load/teardown effects that are not comparable with the mapper hot path.

```text
config                         max_mismatches=6, right softclip=16, RC + TSO/polyA/polyT/lowQ trim
internal wall time             0.669 s
mapping hot-path time          0.512 s
throughput                     149.5k reads/s
unique gene rate               42.77%
same-gene multi-transcript     33.53%
multi-gene ambiguous           10.07%
gene-countable rate            76.30%
unmapped rate                  12.68%
low-complexity rate            0.88%
```

Experimental low-memory two-tier mapping is available with `--hot-index`.
The hot index is an exon-only compact index and `--index` remains the broad exon-plus-gene-body cold index. Cold fallback is scheduled by raw sparse seed prefix, uses no-fallback sparse probing by default, and calls `MADV_DONTNEED` every `--cold-evict-interval` shard groups to cap resident cold pages. On the same 100k subset:

```text
hot index                      human_ensembl93_exon_transcripts.k21.compact.cbidx (4.5G)
cold index                     human_ensembl93_exon_plus_gene_body.targetmeta.k21.genedf.cbidx (42G)
cold shard prefix bits         10
cold evict interval            32
max RSS                        14.4G
process wall time              13.0 s
internal mapping time          12.7 s
gene-countable rate            81.63%
cold fallback reads            48,207
cold fallback replacements     30,790
```

With RSS uncapped (`--cold-evict-interval 0`) and the broad 42G cold index, the same two-tier strategy reached 32.1G max RSS, 4.0s process wall, 2.92s internal mapping, and the same 81.63% gene-countable rate. Setting `--cold-shard-prefix-bits 0` kept RSS similar but degraded wall time to 50s, so shard-prefix scheduling remains important even when the cold mmap is allowed to stay resident.

A compact splici-inspired flanked-intron cold path was also built:

```text
target FASTA                    Homo_sapiens.GRCh38.93.intron_flanks86.targetmeta.fa (48M)
cold index                      human_ensembl93_intron_flanks86.targetmeta.k21.genedf.cbidx (1.1G)
target bases                    45.1M
index build RSS                 8.6G
two-tier max RSS                5.8G
process wall time               2.22s
internal mapping time           1.81s
gene-countable rate             51.72%
cold fallback replacements      877
```

The flanked-intron index is memory-friendly but currently too sparse to recover the broad gene-body recall. The failed all-introns target is 1.6G of FASTA sequence and the indexing attempt was killed at about 56G RSS, so full introns are not a practical cold path in the current explicit-postings index format.

Reverse-complement search is a large recall gain but costs throughput. On the same 100k subset:

```text
mode             gene-countable   unmapped   internal throughput
forward only     68.03%           20.41%     172.4k reads/s
RC enabled       73.98%           15.52%     165.0k reads/s
RC + TSO trim    76.30%           12.68%     149.5k reads/s
```

Class-partitioned compact v4 indexes are now available for testing a unified-index alternative to the two-tier hot/cold path:

```text
index                           human_ensembl93_exon_plus_gene_body.targetmeta.k21.genedf.v4.classpart.cbidx (42G)
plain sparse-probe v4           0.513s mapping, 144k reads/s, 76.40% gene-countable
sparse exon-first v4            0.569s mapping, 134k reads/s, 77.91% gene-countable
sparse exon-first + score rescue 0.737s mapping, 109k reads/s, 78.35% gene-countable
```

The v4 layout does not materially change index size because it reorders existing postings rather than adding a posting-boundary side table. The first exon-first attempt used full seed selection and was rejected as too slow; the current implementation uses sparse exon probes and only falls back to normal sparse/full candidate search for reads without enough exon evidence.

Cell Ranger comparison is currently only against `filtered_feature_bc_matrix`, which is a cell/gene UMI matrix rather than read-level mapping truth. On the 100k PBMC subset, about 91.5% of Constellation UMIs land in Cell Ranger filtered barcodes, all Cell Ranger genes are present in the Ensembl 93 reference, and about 5.6% of Constellation UMIs are assigned to genes outside the Cell Ranger feature set. Log gene-total correlation against the filtered matrix is about 0.51 Pearson and 0.56 Spearman on this tiny read subset, so use it as a rough sanity check rather than a mapping-rate target.

## Handoff Status

Implemented recall and speed work from the Codex mapping brief:

```text
done   gene-countable same-gene transcript ambiguity
done   split assignment metrics and score-failure diagnostics
done   exon-plus-gene-body target generation for intron/pre-mRNA support
done   compact v3 gene-DF index support and gene-idf seed planning
done   sparse-probe candidate search with full fallback
done   streamed seed-batched candidate generation/scoring/assignment
done   WAND-lite and gene-level WAND-lite pruning
done   on-the-fly reverse-complement lookup
done   compact v4 target-class posting partitioning
done   sparse exon-first candidate search over unified v4 indexes
done   TSO trimming and antisense diagnostic assignment class
done   target-class metrics and conservative exon-priority assignment
done   optional post-score rescue for unmapped/antisense reads
done   exon-plus-introns-only target generation mode
done   count-stage 10x barcode whitelist correction, feature whitelist filtering, molecule-gene resolution, and directional UMI correction
done   experimental two-tier hot/cold mapping with shard-sorted cold fallback and configurable cold mmap eviction
```

Open architecture questions for the next handoff:

```text
1. Cell Ranger comparability:
   use the implemented `count --feature-whitelist` mode when comparing against 10x matrices, and decide whether filtered matrices or raw feature references are the right comparator for each benchmark.

2. Target class priority:
   benchmark the implemented exon-priority rule and compare exon-plus-gene-body against exon-plus-introns-only indexes before choosing the default reference architecture.

3. Barcode and UMI correction:
   benchmark the implemented count-stage corrections on full-size datasets and tune whether molecule-gene resolution should skip or retain weak multi-gene molecules.

4. Rescue alignment:
   add a bounded fallback for reads with weak exact-k seed support, such as shorter/spaced seeds plus local alignment, instead of only loosening the full fast path.

5. Reference architecture:
   evaluate a Cell Ranger-like feature reference or splici-style reference so high genome mapping rates do not come from reads that are not countable transcriptome evidence.
```

## Example

```bash
cargo run -p constellation-cli -- index \
  --transcripts tests/tiny_transcriptome.fa \
  --k 15 \
  --out /tmp/tiny.cbidx

cargo run -p constellation-cli -- inspect-index --index /tmp/tiny.cbidx

cargo run -p constellation-cli -- simulate \
  --transcripts tests/tiny_transcriptome.fa \
  --num-reads 1000 \
  --read-len 40 \
  --out-prefix /tmp/sim

cargo run -p constellation-cli -- map \
  --index /tmp/tiny.cbidx \
  --r1 /tmp/sim_R1.fastq \
  --r2 /tmp/sim_R2.fastq \
  --chemistry tenx-3p-v3 \
  --mode candidate-locus \
  --seed-planner raw-frequency \
  --retrieval-mode per-read \
  --out /tmp/assignments.tsv \
  --emit-metrics /tmp/metrics.json

cargo run -p constellation-cli -- bench-report \
  --assignments /tmp/assignments.tsv \
  --truth /tmp/sim_truth.tsv \
  --index /tmp/tiny.cbidx \
  --metrics /tmp/metrics.json

cargo run -p constellation-cli -- count \
  --assignments /tmp/assignments.tsv \
  --index /tmp/tiny.cbidx \
  --barcode-whitelist /path/to/3M-february-2018.txt \
  --feature-whitelist /path/to/features.tsv.gz \
  --emit-metrics /tmp/counts.metrics.json \
  --out-prefix /tmp/counts
```

To compile and run the `pulp` SIMD scoring path through the CLI:

```bash
cargo run -p constellation-cli --features simd -- map \
  --index /tmp/tiny.cbidx \
  --r1 /tmp/sim_R1.fastq \
  --r2 /tmp/sim_R2.fastq \
  --mode candidate-locus \
  --score-mode pulp \
  --out /tmp/assignments.pulp.tsv
```

Mapping modes:

```text
read-at-a-time   baseline mode; scores each read's candidate hits independently
sketch-bucket    sketch-sorted read order without candidate-locus regrouping
candidate-locus  sketch-sorted reads plus candidate hit regrouping by transcript/locus
```

Seed planners:

```text
raw-frequency  prefer selected seeds with fewer raw postings
gene-idf       prefer low gene-document-frequency seeds and weight candidate votes by gene-level IDF
```

Retrieval modes:

```text
per-read      load selected posting lists independently for each read
seed-batched  sort selected query seeds by k-mer, load each posting list once, then reduce candidate votes
```

Target kinds:

```bash
cargo run -p constellation-cli -- build-transcriptome-target \
  --genome /path/to/genome.fa \
  --gtf /path/to/annotation.gtf \
  --target-kind exon-plus-gene-body \
  --out /tmp/target.fa
```

```text
exon-transcripts
gene-bodies
introns-only
intron-flanks
exon-plus-gene-body
exon-plus-introns-only
```

Simulation scenarios:

```text
exact
low-complexity
high-expression-skew
```

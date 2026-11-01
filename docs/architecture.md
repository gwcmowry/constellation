# Architecture

`quant` provides the sample-level workflow: it launches `map --splicing-state`
and `count-splicing` sequentially, checks their record totals, and writes a run
manifest. Separate processes release mapper memory before counting. This is
orchestration of the existing stages, not fused mapping/counting; CSTRAD2 is
still materialized and removed on success unless retention is requested.
`--threads` controls each stage's Rayon pool. The conservative resolver is the
default; `quant --resolution usa` selects exact-UMI USA automatically.

The opt-in `--splicing-state` path specializes consensus mapping to intersect
reference classes per candidate gene, preserving shared, conflicting and unknown
evidence separately. State storage and updates compile out of the ordinary
consensus path; the packed worker dispatches the optional specialization once
per chunk. CSTRAD2 adds per-candidate state bytes and an embedded gene dictionary.
The standalone `count-splicing` command intersects gene/state evidence across
CB/UMI groups and merges state constraints during directional UMI correction.
See [splicing-state output](splicing_state.md) for the exact policy and limits.

The optional `--gene-ec-consensus` policy uses reusable per-worker count and class
arrays to intersect all usable distinct matched-seed gene sets. It requires two
different matched k-mers, ignores absent seeds, and rejects contradictory matched
evidence. Exon/intron annotations describe surviving candidates but never break
a cross-gene tie. This policy keeps the existing index and packed I/O pipeline;
its emitted score counts distinct supporting seeds. It is currently restricted
to forward read-at-a-time mapping, with incompatible options rejected explicitly.
See the [comparison and limitations](../experiments/2026-10-30-consensus/README.md).

The current architecture is the gene-EC mapper. Older positional compact-index modes remain only as historical code paths and should not be used for current performance comparisons.

## Mapper Flow

```text
paired FASTQ streams
  -> zlib-rs gzip decode into three reusable slabs per mate
  -> borrowed 10x barcode/UMI fields and R2 quality checks
  -> fused encoding and low-complexity filtering
  -> sequence-intrinsic closed-syncmer probes
  -> exact PtrHash v9 lookup, at most 32 seeds buffered within a read
  -> per-worker dense epoch-tagged gene score accumulation
  -> direct resolution with reusable tied-gene storage
  -> ordered reusable RAD-like binary or TSV output buffers
  -> optional streaming zstd compression
```

FASTQ decoding runs ahead of mapping through bounded channels. Input buffers
are recycled rather than allocating sequence/quality vectors for every read.
Each Rayon worker owns persistent encoding, scoring and resolution state;
separate ordered output chunks preserve record order even when scheduling
changes. A chunk acquires its worker's state once, with no nested parallel work.
Mapping aggregates metrics and formats records directly, avoiding an owned
assignment object for each read. Diagnostics and historical traversal modes
retain the old richer object path.

## Index

Build directly with `constellation index-ec --format ptrhash --syncmer-s 20`,
supplying the reference FASTA and transcript-to-gene/class table. The default
partitioned builder writes v9 without an intermediate conversion. The previous
streaming builder and mmap layout remain available explicitly; the bare
`index-ec` defaults still select dense mmap for compatibility.

The hashed closed-syncmer rule is sequence- and strand-intrinsic: reference and
query select the same exact k-mer independently of absolute transcript/read
position. A rolling canonical s-mer code and fixed monotonic queue make selection
linear and allocation-free. Sampling metadata is stored in the v4 reserved header
fields; incompatible lexicographic prototype indexes are rejected rather than
silently queried with different semantics.

V9 stores full keys and direct EC payloads in non-minimal perfect-hash slots,
so negative queries are rejected exactly. Pilot metadata is prefetched eight
queries ahead, then payload slots are prefetched before scoring. Buffering is
bounded to 32 seeds within a single read and supports arbitrarily long reads
through successive groups. It does not change sampling, evidence order, scores,
or index format. Other formats retain their lookup implementations.

The measured human v9 index is 4.428GB. In the October 27 weighted-mapping test
on 10M PBMC pairs, the 12-core warm-cache
median is 3.641M reads/s with output disabled and 3.357M reads/s with RAD+zstd,
at 4.265GiB max RSS for compressed output. Both modes preserve all assignment
metrics and complete output bytes. Measurements are in
`experiments/2026-10-27-packed-pipeline/`. The previous p27 sampled layout and
dense indexes remain compatible; rejected densities and rescue policies are
recorded in [the failure ledger](failed_performance_experiments.md). These
weighted-mapper timings are distinct from the newer consensus/state-aware
[100M/500M complete workflow comparison](../experiments/2026-11-01-large-scale/README.md).

The direct mapper reuses one dense score/epoch array per Rayon worker. Reset cost
is proportional to genes touched by the read, and score updates avoid a per-read
hash table. Direct resolution scans touched scores in place; tied top-gene IDs
are collected only for output modes that need them, with capacity reused across
reads.

## Output

The preferred mapper output is `--output-format gene-ec-rad`. It is a Constellation-specific compact per-read format, not alevin-fry-compatible RAD. Streaming zstd is available with `--output-compression zstd`; it is cheap enough for normal benchmark runs but does not remove the need for molecule-level aggregation.

## Current Constraints

- The state counter ingests and partitions evidence serially before resolving cells with Rayon. At 500M pairs, ingestion takes 49.19 of 60.46 default counting seconds.
- Counting retains evidence across the sample; default peak RSS grows from 3.63 GiB at 100M to 14.29 GiB at 500M. Mapping remains near 4.27 GiB.
- The intermediate CSTRAD2 still requires a write/read round trip and about 16.2 GiB at 500M pairs. A fused or bounded-memory path remains future work.
- Native comparisons against piscem/alevin-fry and kallisto/bustools are complete at 100M and 500M. Mapping evidence and molecule policies differ; matched-evidence USA controls establish policy parity within their tested scope, not broad biological equivalence.
- Gene/state compatibility is not transcript-path alignment or isoform quantification. Cell calling, independent-dataset scaling and broader accuracy validation remain separate work.

# Architecture

Constellation follows the dataflow described in the project brief:

```text
FASTQ reads
  -> compact/sketched reads
  -> sketch-sorted read buckets
  -> seed lookup
  -> candidate-locus buckets
  -> scalar/SIMD scoring
  -> gene/transcript assignments
```

The current implementation is correctness-first. The scalar scorer is the oracle, and `PulpScorer` currently shares the same trait boundary while SIMD kernels are added behind it.

Implemented mapping modes:

```text
read-at-a-time
sketch-bucket
candidate-locus
```

`candidate-locus` is the default because it directly exercises the central cache-locality hypothesis: seed lookups are generated from sketch-sorted reads, candidate hits are regrouped by transcript/locus bins, and scoring then works over locus-local buckets.

Current filtering before seed lookup:

```text
low_quality      mean Phred quality below --min-mean-quality
low_complexity   insufficient valid kmers or homopolymer-dominated sequence
```

FASTQ input supports plain text and `.gz` files.

Indexing currently supports:

```text
FASTA transcript sequences through needletail
optional GTF transcript_id -> gene_id mapping
mmap-backed index loading for the prototype .cbidx format
```

Candidate seed selection considers valid k-mers, filters low-quality seeds when qualities are available, and probes lower-frequency k-mers first. This is the first implementation step toward rare-anchor/IDF seed scheduling.

Counting currently emits:

```text
<prefix>_matrix.mtx
<prefix>_barcodes.tsv
<prefix>_features.tsv
```

Counts are exact UMI-deduplicated over `unique_gene` assignments.

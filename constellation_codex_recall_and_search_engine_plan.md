# Constellation: Codex Implementation Plan for Mapping Recall and Search-Engine Candidate Retrieval

**Audience:** Codex / implementation agent  
**Project:** `constellation`, a Rust cache-batched scRNA-seq mapper/quantifier  
**Current problem:** current early prototype processes roughly **170–175k read pairs/s**, but the tested PBMC subset has roughly **50% unmapped reads** and a large combined ambiguity bucket.  
**Primary goal:** improve biologically useful gene assignment and explain every lost read before doing deeper SIMD or WAND-style pruning.

This file is intended to be dropped into the repository as something like:

```text
docs/codex_recall_and_search_engine_plan.md
```

Codex should treat this as a prioritized implementation spec. Implement tasks in order unless the user explicitly asks otherwise.

---

## 0. Current repo facts to preserve

The current repo already has the essential skeleton:

```text
crates/constellation-core/src/
  assign.rs
  batch.rs
  candidate.rs
  chemistry.rs
  dna.rs
  fastq.rs
  gtf.rs
  index.rs
  index_build.rs
  lookup.rs
  metrics.rs
  minimizer.rs
  read_store.rs
  score.rs
  score_scalar.rs
  score_simd.rs
  simulate.rs
  sketch.rs
  transcriptome_target.rs

crates/constellation-cli/src/
  cmd_bench_report.rs
  cmd_build_target.rs
  cmd_count.rs
  cmd_index.rs
  cmd_map.rs
  cmd_simulate.rs
  main.rs
```

Important current behavior:

1. `candidate.rs` already uses the correct seed-offset diagonal for candidate start:

   ```rust
   candidate_start = posting.pos - seed_pos
   ```

   Do not spend time “fixing” candidate start unless tests show a real issue. Add a regression test for this behavior instead.

2. `cmd_map.rs` currently supports three mapping modes:

   ```text
   read-at-a-time
   sketch-bucket
   candidate-locus
   ```

3. `candidate-locus` currently sorts reads by sketch and then regroups candidate hits by transcript/position bin. That improves scoring locality, but candidate generation is still effectively read-centric: each read chooses seeds and loads posting lists independently.

4. `score_scalar.rs` currently uses full-length Hamming scoring with a default of `max_mismatches = 4`; candidates that extend beyond the transcript window are skipped.

5. `assign.rs` already distinguishes:

   ```text
   unique_gene
   ambiguous_gene
   ambiguous_transcript_same_gene
   unmapped
   low_complexity
   low_quality
   ```

6. `cmd_count.rs` currently only counts `unique_gene` rows. This likely discards reads marked `ambiguous_transcript_same_gene`, even though those are gene-countable when `gene_id` is present.

7. `metrics.rs` currently combines `AmbiguousGene` and `AmbiguousTranscriptSameGene` into a single `ambiguous_gene_rate`. That hides a biologically important distinction.

8. `transcriptome_target.rs` currently only builds `TranscriptomeTargetKind::ExonTranscripts`. For modern 10x-style whole-transcriptome GEX, exon-only targets likely leave many usable intronic/pre-mRNA reads unmapped. 10x states that intronic mapped reads can account for about 20–40% of reads and that Cell Ranger v7.0+ counts intronic reads by default for whole-transcriptome gene expression data: https://www.10xgenomics.com/support/software/cell-ranger/latest/miscellaneous/cr-intron-mode-rec

---

## 1. Global rules for Codex

### 1.1 Priorities

Implement in this order:

```text
1. Explain and recover mapping recall.
2. Preserve correctness and conservative ambiguity.
3. Improve biologically useful gene-countable assignment.
4. Add search-engine-inspired candidate retrieval.
5. Only then add WAND/block-pruning optimizations.
6. SIMD is not the bottleneck until recall and candidate generation are understood.
```

### 1.2 Do not do yet

Do **not** implement any of these before diagnostics and recall work:

```text
- aggressive early stopping that can reduce recall;
- expression priors as final assignment evidence;
- dense vector embedding or HNSW-style read-neighbor search;
- full spliced genomic alignment;
- BAM/CRAM output;
- hand-written AVX-512 kernels;
- large unsafe rewrites of the index format without tests.
```

### 1.3 Testing rule

Every change must include at least one of:

```text
- unit test;
- integration test;
- synthetic benchmark validation;
- metrics comparison command.
```

Silent correctness bugs are more dangerous than crashes in this project.

---

## 2. Highest-level diagnosis

The current 49–50% unmapped rate is unlikely to be solved by SIMD. Likely contributors:

```text
1. Exon-only target misses intronic/pre-mRNA reads.
2. Full-length Hamming with max 4 mismatches rejects reads with poly-A tails, adapter/low-quality suffixes, or short inserts.
3. Same-gene transcript ambiguity may be counted as “ambiguous” and then discarded by `count`, reducing useful gene-level signal.
4. Raw posting-frequency seed filtering can throw away gene-informative k-mers duplicated across isoforms of the same gene.
5. Candidate generation is still read-centric, so posting lists are repeatedly loaded instead of being reused search-engine style.
```

The immediate goal is to turn:

```text
unmapped_rate = 49.96%
```

into a decomposition like:

```text
likely_intronic_or_target_missing              18.2%
score_rejected_full_length_mismatch_or_tail     9.4%
all_selected_seeds_over_frequency_cap           6.1%
no_seed_with_postings                           4.0%
candidate_out_of_bounds                         2.2%
reverse_complement_only                         0.3%
other                                           ...
```

Do not optimize blind.

---

## 3. Task 1 — Split assignment semantics and count same-gene transcript ambiguity

### Motivation

`ambiguous_transcript_same_gene` is usually gene-countable for scRNA-seq gene-expression quantification when `gene_id` is present. It should not be merged into true multi-gene ambiguity in metrics, and `count` should not discard it.

### Files likely affected

```text
crates/constellation-core/src/assign.rs
crates/constellation-core/src/metrics.rs
crates/constellation-cli/src/cmd_map.rs
crates/constellation-cli/src/cmd_count.rs
crates/constellation-cli/src/cmd_bench_report.rs
crates/constellation-cli/tests/cli_smoke.rs
tests/expected_assignments.tsv
```

### Required changes

Add or expose these rates in `MapMetrics`:

```rust
pub unique_gene_rate: f64,
pub same_gene_multitranscript_rate: f64,
pub multi_gene_ambiguous_rate: f64,
pub gene_countable_rate: f64,
pub low_complexity_rate: f64,
pub low_quality_rate: f64,
pub unmapped_rate: f64,
```

Keep backward compatibility if helpful by retaining `ambiguous_gene_rate`, but redefine or document it clearly. Preferred:

```text
ambiguous_gene_rate = true multi-gene ambiguity only
same_gene_multitranscript_rate = same-gene transcript ambiguity
```

Update `cmd_map.rs` metrics construction:

```rust
let unique_gene = rate(&assignments, AssignmentType::UniqueGene);
let same_gene_multitranscript = rate(&assignments, AssignmentType::AmbiguousTranscriptSameGene);
let multi_gene_ambiguous = rate(&assignments, AssignmentType::AmbiguousGene);
let gene_countable = unique_gene + same_gene_multitranscript;
```

Update `cmd_count.rs` to count both:

```text
unique_gene
ambiguous_transcript_same_gene with present gene_id
```

Still exclude:

```text
ambiguous_gene
unmapped
low_complexity
low_quality
```

Update `cmd_bench_report.rs` to report correct/evaluable gene assignments for both `unique_gene` and `ambiguous_transcript_same_gene` when `gene_id` is present.

### Acceptance criteria

```text
- `cargo test` passes.
- Metrics JSON contains `same_gene_multitranscript_rate` and `gene_countable_rate`.
- `count` includes same-gene multi-transcript rows with valid gene_id.
- True multi-gene ambiguous rows remain excluded from count matrix.
- Bench report prints separate lines for unique, same-gene transcript ambiguous, multi-gene ambiguous, and gene-countable assignments.
```

### Tests to add

In `assign.rs`:

```rust
#[test]
fn same_gene_multi_transcript_assignment_is_gene_countable() {
    // Two best candidates with same gene, different transcripts.
    // Expected assignment_type: AmbiguousTranscriptSameGene
    // Expected gene_id: Some(gene)
    // Expected transcript_id: None
}
```

In `cmd_count.rs` tests or CLI smoke tests:

```text
assignment_type=ambiguous_transcript_same_gene with gene_id=1 should increment matrix count.
assignment_type=ambiguous_gene should not increment matrix count.
```

---

## 4. Task 2 — Add mapping failure diagnostics before changing algorithms

### Motivation

The current `unmapped_rate` is too coarse. The mapper needs to explain whether reads are lost because seeds are absent, seeds are over frequency cap, candidates fail score, candidates are out of bounds, or reads are filtered before candidate generation.

### Files likely affected

```text
crates/constellation-core/src/candidate.rs
crates/constellation-core/src/score.rs
crates/constellation-core/src/score_scalar.rs
crates/constellation-core/src/score_simd.rs
crates/constellation-core/src/metrics.rs
crates/constellation-cli/src/cmd_map.rs
crates/constellation-cli/src/main.rs
```

### Add data structures

In `candidate.rs`, extend `CandidateGenerationStats` or add a separate struct:

```rust
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub struct CandidateFailureStats {
    pub reads_with_no_valid_kmers: u64,
    pub reads_with_no_selected_seeds: u64,
    pub reads_with_no_seed_postings: u64,
    pub reads_all_selected_seeds_over_frequency_cap: u64,
    pub reads_with_candidate_hits: u64,
    pub selected_seeds_absent_from_index: u64,
    pub selected_seeds_over_frequency_cap: u64,
    pub selected_seeds_with_postings: u64,
}
```

In scoring, add:

```rust
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub struct ScoreFailureStats {
    pub candidates_seen: u64,
    pub candidates_out_of_bounds: u64,
    pub candidates_failed_mismatch: u64,
    pub candidates_passed_full_length: u64,
    pub candidates_passed_trimmed_or_softclipped: u64,
}
```

Keep stats additive and easy to reduce across threads.

### Required changes

1. Track the number of reads that have zero valid k-mers after `iter_kmers_2bit`.
2. Track selected seeds that have zero postings.
3. Track selected seeds skipped due to frequency cap.
4. Track reads where all selected seeds were absent or over cap.
5. Track reads that produced candidate hits but later yielded no scored candidate.
6. Track scoring failures separately from candidate-generation failures.

Add a CLI option:

```bash
--emit-unmapped-diagnostics path.tsv
```

The TSV should be sampled if full output is too large. Suggested columns:

```tsv
read_id	reason	seq_len	num_valid_kmers	selected_seeds	seed_hits	candidate_hits	scored_candidates	flags
```

Suggested `reason` values:

```text
low_quality
low_complexity
invalid_barcode_umi
no_valid_kmers
no_selected_seeds
no_seed_with_postings
all_selected_seeds_over_frequency_cap
candidate_hits_but_all_out_of_bounds
candidate_hits_but_score_failed
unmapped_unknown
```

### Acceptance criteria

```text
- Metrics JSON reports counts for all failure reasons.
- `--emit-unmapped-diagnostics` writes a TSV with header and at least one row for unmapped reads.
- Existing output assignment TSV remains unchanged unless diagnostics flag is used.
- `cargo test` passes.
```

### Benchmark command

```bash
cargo run -p constellation-cli --release -- map \
  --index /path/to/index.cbidx \
  --r1 /path/to/R1.fastq.gz \
  --r2 /path/to/R2.fastq.gz \
  --chemistry tenx-3p-v3 \
  --mode candidate-locus \
  --out /tmp/assignments.tsv \
  --emit-metrics /tmp/metrics.json \
  --emit-unmapped-diagnostics /tmp/unmapped.tsv
```

---

## 5. Task 3 — Add soft-clipped / trimmed scoring

### Motivation

The current scalar scorer performs full-length Hamming and rejects reads when the full read extends beyond the transcript or has more than 4 mismatches. This is too brittle for real 10x data with terminal poly-A/T, adapter-like suffixes, low-quality tails, or short inserts.

### Files likely affected

```text
crates/constellation-core/src/score.rs
crates/constellation-core/src/score_scalar.rs
crates/constellation-core/src/score_simd.rs
crates/constellation-cli/src/main.rs
crates/constellation-cli/src/cmd_map.rs
crates/constellation-core/src/metrics.rs
```

### Add scoring config

Create a config passed into both scalar and SIMD scorers:

```rust
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct ScoreConfig {
    pub max_mismatches: u32,
    pub max_right_softclip: u16,
    pub max_left_softclip: u16,
    pub trim_poly_a: bool,
    pub trim_poly_t: bool,
    pub trim_low_quality_tail: bool,
    pub min_scored_len: u16,
    pub min_tail_phred: u8,
}

impl Default for ScoreConfig {
    fn default() -> Self {
        Self {
            max_mismatches: 4,
            max_right_softclip: 0,
            max_left_softclip: 0,
            trim_poly_a: false,
            trim_poly_t: false,
            trim_low_quality_tail: false,
            min_scored_len: 35,
            min_tail_phred: 10,
        }
    }
}
```

Add CLI args:

```bash
--max-mismatches 4
--max-right-softclip 16
--max-left-softclip 0
--trim-poly-a
--trim-poly-t
--trim-low-quality-tail
--min-scored-length 35
--min-tail-phred 10
```

### Scoring algorithm

For each candidate:

1. Try full-length Hamming exactly as today.
2. If it passes, emit a normal scored candidate.
3. If it fails and trimming/soft-clipping is enabled:
   - compute a scoring slice for the read by trimming terminal low-quality bases if enabled;
   - trim terminal homopolymer A/T if enabled;
   - allow right soft clip up to `max_right_softclip`;
   - optionally allow left soft clip up to `max_left_softclip`, but default should be 0;
   - require `scored_len >= min_scored_len`;
   - penalize clipping so a full-length clean match wins over a shorter clipped match.

Suggested score:

```rust
score = scored_len - mismatches - softclip_penalty
```

Use saturating arithmetic. Keep score as `u16` for now.

### Flags

Add candidate flags in `score.rs`:

```rust
pub const SCORE_FLAG_FULL_LENGTH: u16 = 1 << 0;
pub const SCORE_FLAG_TRIMMED_POLY_A: u16 = 1 << 1;
pub const SCORE_FLAG_TRIMMED_POLY_T: u16 = 1 << 2;
pub const SCORE_FLAG_TRIMMED_LOW_QUALITY: u16 = 1 << 3;
pub const SCORE_FLAG_RIGHT_SOFTCLIP: u16 = 1 << 4;
pub const SCORE_FLAG_LEFT_SOFTCLIP: u16 = 1 << 5;
```

Consider extending `ScoredCandidate` with:

```rust
pub scored_len: u16,
pub right_softclip: u16,
pub left_softclip: u16,
```

If this causes too much churn in TSV output, keep these as internal metrics/flags first.

### Acceptance criteria

```text
- Existing exact synthetic tests still pass with default config.
- With softclip options disabled, output should match current scalar output.
- New tests recover reads with terminal poly-A tails.
- New tests recover reads with low-quality right tails when low-quality trimming is enabled.
- Metrics report full-length pass vs trimmed/soft-clipped pass counts.
- SIMD scorer either supports the same config or explicitly falls back to scalar for non-default clipping mode.
```

### Tests to add

```text
score_scalar::poly_a_tail_can_softclip
score_scalar::low_quality_tail_can_trim
score_scalar::default_scoring_matches_old_behavior
score_scalar::full_length_match_beats_softclipped_match
```

---

## 6. Task 4 — Add gene-body and intron target support

### Motivation

The current target builder only supports spliced exon transcripts. For 10x whole-transcriptome GEX, intronic reads can be usable and are counted by default in modern Cell Ranger. Exon-only targets likely explain a major part of the unmapped rate.

### Files likely affected

```text
crates/constellation-core/src/transcriptome_target.rs
crates/constellation-core/src/gtf.rs
crates/constellation-core/src/index_build.rs
crates/constellation-core/src/index.rs
crates/constellation-cli/src/cmd_build_target.rs
crates/constellation-cli/src/main.rs
```

### Add target kinds

Update:

```rust
pub enum TranscriptomeTargetKind {
    ExonTranscripts,
    GeneBodies,
    IntronsOnly,
    ExonPlusGeneBody,
}
```

Add CLI:

```bash
constellation build-transcriptome-target \
  --genome genome.fa.gz \
  --gtf genes.gtf.gz \
  --target-kind exon-transcripts \
  --out target.fa
```

Supported values:

```text
exon-transcripts
gene-bodies
introns-only
exon-plus-gene-body
```

### Header format

Use consistent headers and make the index parser accept all older variants:

```text
>{target_id}|gene:{gene_id}|target:{target_kind}|contig:{seqname}|strand:{+|-}|start:{start}|end:{end}
```

Backward-compatible gene parsing must accept:

```text
|gene:GENE
|gene=GENE
gene_id "GENE"
gene_id=GENE
gene_id:GENE
```

### Gene body behavior

For `GeneBodies`:

1. For each gene, derive a genomic interval.
2. Prefer explicit `gene` features in GTF when present.
3. If no gene feature exists, derive gene body as min start / max end over all exons for the gene.
4. Emit one gene-level target per gene.
5. Reverse-complement minus-strand gene bodies so sequence is in transcript orientation, consistent with exon transcript behavior.

### Introns-only behavior

For `IntronsOnly`:

1. Build merged exon intervals per gene.
2. Build gene body interval.
3. Subtract merged exons from gene body to get intronic intervals.
4. Emit either:
   - one concatenated intron target per gene; or
   - multiple intron segment targets per gene.
5. Prefer one concatenated intron target per gene for first implementation to keep assignment simple.

### ExonPlusGeneBody behavior

For `ExonPlusGeneBody`:

Emit the existing exon transcript targets plus gene-body targets. Mark `target:exon_transcript` vs `target:gene_body` in headers.

### Assignment implications

Initially, gene-body/intron targets can map directly to `gene_id`. Do not require transcript-level precision for intronic reads.

If possible, add target class metadata to the index:

```rust
pub enum TargetClass {
    ExonTranscript,
    GeneBody,
    Intron,
}
```

If too much churn, encode target class in transcript name and add metrics later.

### Acceptance criteria

```text
- `build-transcriptome-target --target-kind exon-transcripts` produces output matching current behavior, aside from optional header metadata.
- `--target-kind gene-bodies` emits one gene-level target per gene.
- `--target-kind introns-only` emits intronic sequence and excludes exon bases in a simple synthetic GTF.
- `--target-kind exon-plus-gene-body` emits both classes.
- Synthetic intronic reads map to the correct gene-level target.
- Metrics can compare exon-only vs gene-body/exon-plus-gene-body mapping rate.
```

### Tests to add

```text
transcriptome_target::builds_gene_body_targets_from_gene_feature
transcriptome_target::builds_gene_body_targets_from_exon_bounds_without_gene_feature
transcriptome_target::builds_introns_only_by_subtracting_merged_exons
transcriptome_target::exon_plus_gene_body_contains_both_target_classes
index_build::parses_gene_colon_and_gene_equals_headers
```

---

## 7. Task 5 — Add gene-level document frequency and IDF-weighted seed planning

### Motivation

Current seed planning uses raw posting count. That can misclassify useful k-mers as too frequent when they appear in many isoforms of the same gene. For gene-level scRNA-seq assignment, the relevant search-engine analogue is **gene document frequency**, not raw transcript-position frequency.

### Files likely affected

```text
crates/constellation-core/src/index.rs
crates/constellation-core/src/index_build.rs
crates/constellation-core/src/candidate.rs
crates/constellation-cli/src/cmd_index.rs
crates/constellation-cli/src/cmd_map.rs
```

### Add k-mer stats

Add:

```rust
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub struct KmerStats {
    pub raw_postings: u32,
    pub transcript_df: u32,
    pub gene_df: u32,
}
```

Extend `IndexAccess`:

```rust
fn kmer_stats(&self, kmer_code: u64) -> Option<KmerStats>;
```

Maintain backward compatibility:

```text
- Compact v2 indexes can return KmerStats { raw_postings: posting_count, transcript_df: posting_count, gene_df: posting_count }.
- New compact v3 indexes should store real transcript_df and gene_df.
```

### Serialization approach

Preferred simple approach:

```text
- Add KmerStats fields to serialized KmerEntry / CompactKmerEntry.
- Bump compact version to 3.
- Keep loader support for v2 by detecting version and filling stats from postings_len.
```

Alternative if preserving 16-byte hot k-mer entries is important:

```text
- Keep CompactKmerEntry at 16 bytes.
- Add a parallel compact stats array aligned by kmer index.
- Access stats only in seed planning.
```

Codex may choose the simpler route first, but document the tradeoff in code comments.

### Compute stats during index build

For each k-mer:

```text
raw_postings = postings.len()
transcript_df = distinct transcript_id count
gene_df = distinct gene_id count
```

Use sorted postings to compute distinct counts without huge temporary sets where possible.

### Seed planning changes

Replace the tuple-based seed choice:

```rust
SmallVec<[(usize, u32, u64, u8); 256]>
```

with a struct:

```rust
#[derive(Debug, Clone, Copy)]
pub struct QuerySeedChoice {
    pub kmer_code: u64,
    pub read_offset: u32,
    pub read_strand: u8,
    pub raw_postings: u32,
    pub transcript_df: u32,
    pub gene_df: u32,
    pub idf_q8: u16,
    pub min_phred: u8,
}
```

Compute approximate IDF:

```rust
idf = ln((num_genes + 1) / (gene_df + 1))
idf_q8 = clamp(round(idf * 256.0), 0, u16::MAX)
```

Sort seed choices by:

```text
1. absent seeds last;
2. low gene_df first;
3. high idf first;
4. lower raw_postings first;
5. positional diversity if implemented;
6. lower read_offset as deterministic tiebreaker.
```

### Fallback passes

Avoid a single hard frequency cap that kills recall. Implement staged seed selection:

```text
Pass 1: gene_df <= 32, max M seeds
Pass 2: if no candidate, gene_df <= 256 or 1024, max additional M seeds
Pass 3: diagnostic-only, allow common seeds with cap
```

Add CLI options:

```bash
--seed-planner raw-frequency|gene-idf
--max-gene-df-pass1 32
--max-gene-df-pass2 1024
--seed-fallback-pass2
```

Default can remain current raw-frequency until tests pass, then switch to gene-idf.

### Candidate scoring changes

Use seed score as weighted support:

```rust
seed_score += idf_q8
seed_count += 1
```

Keep assignment conservative. Do not make high IDF alone override a candidate that fails sequence scoring.

### Acceptance criteria

```text
- `inspect-index` reports raw posting, transcript_df, and gene_df histograms.
- Synthetic isoform-duplication test shows same-gene multi-isoform k-mers have high raw postings but low gene_df.
- Gene-IDF planner does not reduce synthetic exact-read recall.
- PBMC subset metrics compare raw-frequency vs gene-idf seed planner.
```

### Tests to add

```text
index_build::computes_gene_df_lower_than_raw_postings_for_isoforms
candidate::gene_idf_planner_prefers_low_gene_df_seed
candidate::gene_idf_fallback_recovers_when_rare_seeds_absent
```

---

## 8. Task 6 — Implement search-engine style seed-occurrence batching

### Motivation

This is the highest-value cache/search-engine architecture change. Current candidate generation still performs posting-list lookup per read. A search engine loads a posting list once and applies it to all queries containing that term. The analogous change here is:

```text
per-read seed lookup
```

becomes:

```text
extract seed occurrences for many reads -> sort by kmer_code -> load each posting list once -> emit candidate votes -> reduce votes to candidate hits
```

### Files likely affected

```text
crates/constellation-core/src/candidate.rs
crates/constellation-core/src/batch.rs
crates/constellation-core/src/metrics.rs
crates/constellation-cli/src/cmd_map.rs
crates/constellation-cli/src/main.rs
benches/lookup.rs
benches/end_to_end_synthetic.rs
```

### Add retrieval mode

Add CLI:

```bash
--retrieval-mode per-read|seed-batched
```

Default should be `per-read` until tests prove equivalence, then consider switching.

### New structs

```rust
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
pub struct QuerySeed {
    pub kmer_code: u64,
    pub read_id: ReadId,
    pub read_offset: u32,
    pub read_strand: u8,
    pub weight_q8: u16,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
pub struct CandidateVote {
    pub read_id: ReadId,
    pub transcript_id: TranscriptId,
    pub gene_id: GeneId,
    pub candidate_start: u32,
    pub strand: u8,
    pub weight_q8: u16,
}
```

### Seed-batched algorithm

Within each sketch bucket or map batch:

```text
1. For every non-filtered read, extract selected QuerySeed records.
2. Sort QuerySeed by kmer_code.
3. For each kmer_code group:
   a. Load posting list once.
   b. If posting list is absent, record diagnostics and continue.
   c. If posting list exceeds pass-specific cap, record diagnostics and continue.
   d. For each query seed occurrence and each posting:
        if posting.pos >= read_offset:
            candidate_start = posting.pos - read_offset
            emit CandidateVote
4. Sort CandidateVote by (read_id, transcript_id, candidate_start, strand).
5. Reduce adjacent votes into CandidateHit:
   - seed_count = number of distinct supporting seeds/votes
   - seed_score = sum weight_q8, saturating to u16
6. Continue existing candidate-locus bucketing and scoring.
```

### Important guardrails

The cross-product of seed occurrences and postings can explode. Prevent this with:

```text
- IDF/gene_df seed selection before batching;
- max postings per seed;
- max seed occurrences per kmer group if necessary;
- diagnostic counters for skipped groups;
- fallback to per-read mode for problematic groups if helpful.
```

Do not allow seed-batched mode to silently skip more reads than per-read mode without metrics.

### Metrics to add

```rust
pub seed_occurrences_per_read: f64,
pub distinct_seed_lists_loaded: u64,
pub posting_list_loads_per_read: f64,
pub posting_list_reuse_factor: f64,
pub candidate_votes_per_read: f64,
pub seed_batched_skipped_large_seed_groups: u64,
```

Definitions:

```text
seed_occurrences_per_read = query_seed_count / read_count
posting_list_reuse_factor = query_seed_count / distinct_seed_lists_loaded
candidate_votes_per_read = candidate_vote_count / read_count
```

### Acceptance criteria

```text
- On exact synthetic reads, seed-batched mode produces the same assignments as per-read mode.
- On current tiny tests, output is deterministic.
- Metrics show fewer distinct posting-list loads than query seed occurrences.
- Candidate generation speed is benchmarked for per-read vs seed-batched mode.
- If seed-batched mode changes recall on public PBMC subset, unmapped diagnostics explain why.
```

### Tests to add

```text
candidate::seed_batched_matches_per_read_exact_reads
candidate::seed_batched_handles_mid_read_seed_offsets
candidate::seed_batched_skips_too_frequent_seed_with_metric
candidate::candidate_vote_reduction_sums_seed_score_and_count
```

### Benchmark command

```bash
cargo run -p constellation-cli --release -- map \
  --index /path/to/index.cbidx \
  --r1 /path/to/R1.fastq.gz \
  --r2 /path/to/R2.fastq.gz \
  --chemistry tenx-3p-v3 \
  --mode candidate-locus \
  --retrieval-mode seed-batched \
  --seed-planner gene-idf \
  --emit-metrics /tmp/seed_batched.metrics.json \
  --out /tmp/seed_batched.assignments.tsv
```

---

## 9. Task 7 — Add WAND-like pruning only after recall is healthy

### Motivation

Search engines use upper bounds to avoid scoring candidates that cannot win. This is valuable, but premature pruning can worsen the current mapping-rate problem. Implement only after Tasks 1–6 have good diagnostics and recall.

### Files likely affected

```text
crates/constellation-core/src/candidate.rs
crates/constellation-core/src/index.rs
crates/constellation-core/src/index_build.rs
crates/constellation-core/src/metrics.rs
```

### Upgrade current early stop

Current early stop is count-based top-vs-runner-up. Replace or supplement it with weighted gene-level upper bounds:

```text
current_best_gene_score
current_second_gene_score
remaining_possible_idf_score
```

Stop only when:

```text
best_gene_score > second_gene_score + remaining_possible_idf_score
AND best has enough distinct supporting seeds
AND best support is diagonal-consistent
AND at least min_seed_lookups have been performed
```

### Block-level metadata later

For long posting lists, add block metadata:

```rust
pub struct PostingBlockMeta {
    pub kmer_code: u64,
    pub postings_start: u32,
    pub postings_len: u32,
    pub max_weight_q8: u16,
    pub min_gene_id: GeneId,
    pub max_gene_id: GeneId,
}
```

Do not implement full block-max WAND before seed-batched candidate generation is stable.

### Acceptance criteria

```text
- WAND/early-stop mode is off by default.
- Enabling it cannot change synthetic exact-read assignment unless explicitly allowed by test scenario.
- Metrics report early-stopped reads and seed lookups saved.
- PBMC run reports speed and assignment-rate delta with early stop on/off.
```

---

## 10. Required experiment matrix

After each major task, run a small matrix. Keep the same 100k PBMC subset so deltas are comparable.

### 10.1 Baseline

```bash
cargo run -p constellation-cli --release -- map \
  --index exon_transcripts.cbidx \
  --r1 PBMC_R1.fastq.gz \
  --r2 PBMC_R2.fastq.gz \
  --chemistry tenx-3p-v3 \
  --mode candidate-locus \
  --out /tmp/base.assignments.tsv \
  --emit-metrics /tmp/base.metrics.json
```

### 10.2 Sensitivity: relaxed candidate generation

```bash
cargo run -p constellation-cli --release -- map \
  --index exon_transcripts.cbidx \
  --r1 PBMC_R1.fastq.gz \
  --r2 PBMC_R2.fastq.gz \
  --chemistry tenx-3p-v3 \
  --mode candidate-locus \
  --max-seeds-per-read 24 \
  --max-postings-per-seed 4096 \
  --min-seed-quality 0 \
  --early-stop-posterior 1.0 \
  --out /tmp/relaxed.assignments.tsv \
  --emit-metrics /tmp/relaxed.metrics.json
```

Interpretation:

```text
mapping jumps      -> seed cap/frequency cap/quality filter is too aggressive
mapping unchanged  -> target/scoring/biology issue is more likely
```

### 10.3 Orientation test

```bash
cargo run -p constellation-cli --release -- map \
  --index exon_transcripts.cbidx \
  --r1 PBMC_R1.fastq.gz \
  --r2 PBMC_R2.fastq.gz \
  --chemistry tenx-3p-v3 \
  --mode candidate-locus \
  --search-reverse-complement \
  --out /tmp/rc.assignments.tsv \
  --emit-metrics /tmp/rc.metrics.json
```

Interpretation:

```text
large jump -> strand/orientation/index convention issue
small/no jump -> orientation probably not primary issue
```

### 10.4 Scoring sensitivity

After Task 3:

```bash
cargo run -p constellation-cli --release -- map \
  --index exon_transcripts.cbidx \
  --r1 PBMC_R1.fastq.gz \
  --r2 PBMC_R2.fastq.gz \
  --chemistry tenx-3p-v3 \
  --mode candidate-locus \
  --max-mismatches 6 \
  --max-right-softclip 16 \
  --trim-poly-a \
  --trim-low-quality-tail \
  --min-scored-length 35 \
  --out /tmp/softclip.assignments.tsv \
  --emit-metrics /tmp/softclip.metrics.json
```

### 10.5 Target sensitivity

After Task 4:

```bash
constellation build-transcriptome-target \
  --genome genome.fa.gz \
  --gtf genes.gtf.gz \
  --target-kind exon-plus-gene-body \
  --out exon_plus_gene_body.fa

constellation index \
  --transcripts exon_plus_gene_body.fa \
  --k 21 \
  --out exon_plus_gene_body.cbidx

constellation map \
  --index exon_plus_gene_body.cbidx \
  --r1 PBMC_R1.fastq.gz \
  --r2 PBMC_R2.fastq.gz \
  --chemistry tenx-3p-v3 \
  --mode candidate-locus \
  --out /tmp/exon_plus.assignments.tsv \
  --emit-metrics /tmp/exon_plus.metrics.json
```

### 10.6 Search-engine candidate retrieval

After Task 6:

```bash
constellation map \
  --index exon_plus_gene_body.cbidx \
  --r1 PBMC_R1.fastq.gz \
  --r2 PBMC_R2.fastq.gz \
  --chemistry tenx-3p-v3 \
  --mode candidate-locus \
  --seed-planner gene-idf \
  --retrieval-mode seed-batched \
  --out /tmp/seed_batched.assignments.tsv \
  --emit-metrics /tmp/seed_batched.metrics.json
```

---

## 11. Success metrics

### 11.1 Minimum success for this phase

```text
- Unmapped reads are decomposed into actionable failure reasons.
- same_gene_multitranscript is counted as gene-countable.
- `count` output uses same-gene transcript ambiguous reads correctly.
- Soft clipping/trimming and/or gene-body targets improve gene_countable_rate without obvious false unique inflation.
```

### 11.2 Strong success

```text
- gene_countable_rate increases materially on PBMC subset.
- unmapped_rate decreases materially after exon+gene-body target and soft-clipping.
- multi_gene_ambiguous_rate does not explode.
- synthetic exact and error benchmarks preserve correctness.
- seed-batched retrieval reduces repeated posting-list loads and improves candidate-generation time or cache counters.
```

### 11.3 Performance metrics to watch

```text
reads_per_second
candidate_generation_seconds
scoring_seconds
seed_lookups_per_read
seed_occurrences_per_read
posting_list_loads_per_read
posting_list_reuse_factor
candidate_hits_per_read
candidate_votes_per_read
scored_candidates_per_read
mean_candidate_bucket_size
LLC_miss_rate if perf is available
```

### 11.4 Biological metrics to watch

```text
unique_gene_rate
same_gene_multitranscript_rate
multi_gene_ambiguous_rate
gene_countable_rate
unmapped_rate
low_complexity_rate
low_quality_rate
```

---

## 12. Suggested implementation order as Codex prompts

Use these as sequential Codex prompts.

### Prompt A — Assignment metrics and count semantics

```text
Implement Task 1 from docs/codex_recall_and_search_engine_plan.md.
Focus only on splitting assignment metrics and making `count` include `ambiguous_transcript_same_gene` rows with valid gene_id. Add tests. Do not change candidate generation or scoring.
```

### Prompt B — Failure diagnostics

```text
Implement Task 2 from docs/codex_recall_and_search_engine_plan.md.
Add mapping failure reason counters and `--emit-unmapped-diagnostics`. Add tests with tiny synthetic cases. Do not change assignment decisions except for reporting diagnostics.
```

### Prompt C — Soft-clipped scoring

```text
Implement Task 3 from docs/codex_recall_and_search_engine_plan.md.
Add ScoreConfig and CLI options for max mismatches, right soft clip, poly-A/T trimming, low-quality tail trimming, and min scored length. Preserve old behavior when options are disabled. Add scalar tests. SIMD can fall back to scalar for non-default clipping mode if needed.
```

### Prompt D — Gene-body and intron targets

```text
Implement Task 4 from docs/codex_recall_and_search_engine_plan.md.
Extend build-transcriptome-target with target-kind values: exon-transcripts, gene-bodies, introns-only, exon-plus-gene-body. Add tests on a small GTF/genome. Preserve current exon-transcripts behavior.
```

### Prompt E — Gene_df / IDF seed planner

```text
Implement Task 5 from docs/codex_recall_and_search_engine_plan.md.
Add KmerStats with raw_postings, transcript_df, and gene_df. Add IndexAccess::kmer_stats. Implement a gene-idf seed planner behind `--seed-planner`. Preserve raw-frequency mode. Add inspect-index histograms and tests.
```

### Prompt F — Seed-occurrence batching

```text
Implement Task 6 from docs/codex_recall_and_search_engine_plan.md.
Add `--retrieval-mode per-read|seed-batched`. Implement QuerySeed extraction, kmer_code sorting, posting-list reuse, CandidateVote emission, and vote reduction into CandidateHit. Synthetic outputs must match per-read mode before public data benchmarking.
```

### Prompt G — WAND-like pruning

```text
Only after Tasks 1–6 are green: implement Task 7 from docs/codex_recall_and_search_engine_plan.md. Keep WAND/weighted early stop off by default and prove synthetic recall is preserved.
```

---

## 13. Notes on expected failure modes

### 13.1 Gene-body target can increase ambiguity

Adding gene bodies/introns may recover reads but can also increase multi-gene ambiguity in paralogous regions. That is acceptable if it is reported honestly. Prefer conservative `ambiguous_gene` over false `unique_gene`.

### 13.2 Soft clipping can create false positives

Keep `min_scored_len` conservative. A 20 bp clipped fragment is not enough. Start with `min_scored_len = 35` or higher.

### 13.3 Gene-IDF can bias against repeats but must not erase real biology

Use fallback passes. Do not hard-drop all high-frequency k-mers unless diagnostics prove it is safe.

### 13.4 Seed-batched mode can explode candidate votes

The cross-product of query seed occurrences and posting lists can become huge. Always track `candidate_votes_per_read` and skipped large groups.

---

## 14. Definition of done for this document

This implementation phase is done when the repo can produce a comparison table like:

```text
configuration                  reads/s   gene_countable   unique_gene   same_gene_tx   multi_gene_ambig   unmapped   notes
baseline exon-only              ...        ...             ...           ...            ...                ...        current behavior
+ split metrics/count           ...        ...             ...           ...            ...                ...        count semantics fixed
+ softclip/trim                 ...        ...             ...           ...            ...                ...        recover tail failures
+ exon+gene-body target         ...        ...             ...           ...            ...                ...        recover intronic/gene-body
+ gene-idf seed planner         ...        ...             ...           ...            ...                ...        better seed selection
+ seed-batched retrieval        ...        ...             ...           ...            ...                ...        posting-list reuse
```

and an unmapped diagnostics table like:

```text
reason                                      count       fraction
low_quality                                 ...         ...
low_complexity                              ...         ...
invalid_barcode_umi                         ...         ...
no_valid_kmers                              ...         ...
no_seed_with_postings                       ...         ...
all_selected_seeds_over_frequency_cap       ...         ...
candidate_hits_but_all_out_of_bounds        ...         ...
candidate_hits_but_score_failed             ...         ...
unmapped_unknown                            ...         ...
```

Only after those tables exist should the project focus on deeper pruning, compressed posting-list block metadata, or more SIMD.

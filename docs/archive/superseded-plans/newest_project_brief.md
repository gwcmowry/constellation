According to a document from 2026-05-04, yes: mapping/assignment quality should remain the primary objective for the next phase, but with a narrow definition of “quality.” You are no longer at the “why is half the dataset unmapped?” stage. The current README shows 76.30% gene-countable, 12.68% unmapped, 10.07% multi-gene ambiguous, and 149.5k reads/s on the PBMC 100k subset using seed-batched retrieval, sparse-probe search, WAND-lite, RC search, TSO/polyA/polyT/lowQ trimming, and softclip scoring.  ￼  ￼

So I would frame the next phase as:

First: make the mapping semantics biologically comparable to Cell Ranger/STARsolo-style outputs. Then: optimize the hot path under those semantics.

Not just “increase mapping rate at all costs.”

1. Should mapping rate come before throughput?

Yes — but I would be precise: optimize gene-countable, biologically defensible assignment, not raw mapping rate.

A higher raw mapping rate can be bad if it comes from permissive gene-body matching, over-softclipping, antisense leakage, paralog/pseudogene confusion, or assigning reads that should be multi-gene ambiguous. The repo already notes the open problems that matter here: feature whitelist/evaluation mode, target class priority, barcode/UMI correction, bounded rescue alignment, and better reference architecture.  ￼

The next recall goal should probably be:

unmapped:                 12.68% → <8–10%
multi-gene ambiguous:     10.07% → ideally lower, not higher
gene-countable:           76.30% → 80–85%+
same-gene multi-tx:       33.53% → better explained by target class

The key is to avoid turning the remaining 12.68% unmapped into 8% more ambiguous garbage. The most important thing is to classify the remaining unmapped reads by cause:

target/reference missing
seed support too weak
sparse-probe accepted too early
scoring rejected
RC-only / orientation issue
polyA/TSO/adapter issue
feature-set mismatch
true low complexity
true multi-gene/repetitive

The project plan explicitly says to “explain and recover mapping recall” before aggressive WAND/block pruning or SIMD work, and to avoid implementing aggressive early stopping before diagnostics and recall work.  ￼ That still seems right.

2. What I would do next for mapping quality

The highest-leverage next steps are not generic “loosen thresholds.” They are semantic improvements.

A. Target class priority

Your current best reference is exon-plus-gene-body. That likely helps intronic/pre-mRNA recall, but it can also duplicate exonic signal: the same read may hit both transcript targets and gene-body targets for the same gene. That can inflate the 33.53% same-gene multi-transcript bucket.

I would implement:

score exonic transcript candidates first
if confident same-gene exonic assignment exists, skip gene-body/intron candidates
only score intronic/gene-body candidates for unresolved reads

This could improve both quality and speed.

B. Exon-plus-introns-only / splici-style target

Instead of:

exon transcripts + full gene body

try:

exon transcripts + introns-only gene-level targets

or a splici-style reference. The repo already lists exon, gene-body, introns-only, and exon-plus-gene-body target kinds, but the README flags reference architecture as an open question.  ￼

This may reduce redundant same-gene hits while preserving intronic read recovery.

C. Post-score rescue

The sparse-probe path may be fast but can still lose reads if the sparse candidate set had seed support but then failed scoring. A good rescue path is:

run sparse-probe
score candidates
assign reads
for unmapped / weak / antisense-only reads:
    rerun full candidate search
    rescore
    replace only if better/gene-countable

This is better than loosening the fast path globally.

D. Feature whitelist mode

The README says about 5.6% of Constellation UMIs are assigned to genes outside the Cell Ranger feature set, and the current comparison is only against filtered_feature_bc_matrix, not read-level truth.  ￼ That means some “disagreement” is probably reference/feature-set mismatch, not mapping failure.

Add a mode that labels assignments as:

inside_10x_feature_set
outside_10x_feature_set
feature_set_missing

before trying to optimize concordance.

E. Barcode/UMI correction

For publication/competitiveness, this matters more than squeezing another few percent of read-level assignment. Without 10x-style barcode correction and molecule-level UMI correction, count-matrix comparisons will remain noisy.

3. Where the biggest efficiency gains probably remain

Based on the README, the current fast path already implemented several major performance ideas: seed-batched retrieval, streamed chunks, WAND-lite pruning, compact v3 gene-DF index, on-the-fly RC lookup, and 16-byte k-mer entries / 8-byte postings for cache-friendlier candidate search tables.  ￼

So the next efficiency gains are likely in four places:

1. avoid doing expensive work for reads already confidently resolved;
2. reduce candidate explosion from gene-body/exon duplication;
3. make scoring regular enough for SIMD again;
4. further reduce random index/posting access.

4. Best SIMD opportunities

The current SIMD path is probably not helping much in the highest-recall configuration, because trimming and softclip scoring push you toward scalar fallback. The README says the project has scalar scoring plus a pulp Hamming path, but the current best config uses softclip scoring and TSO/polyA/polyT/lowQ trimming.  ￼

The best strategy is not to SIMD the full softclip scorer first. Instead, use SIMD as a fast rejection/acceptance layer.

A. Two-stage scoring

For each candidate:

Stage 1: SIMD fixed-window Hamming on pre-trimmed read view
Stage 2: scalar softclip fallback only if Stage 1 fails but seed evidence is strong

Metrics to add:

fast_score_attempts
fast_score_passes
softclip_fallback_attempts
softclip_fallback_passes
softclip_avoided_candidates
scoring_ns_per_candidate

This lets SIMD handle the common easy cases while scalar softclip rescues the hard cases.

B. SIMD packed 2-bit mismatch counting

If reads and reference windows are 2-bit encoded, a very good kernel is:

xor read_word/ref_word
collapse 2-bit fields to 1 mismatch bit per base
popcount mismatch bits
compare to threshold

For ~90–100 bp reads, this is only a few u64 words per candidate. SIMD can process multiple candidates at once if the score block is structure-of-arrays.

The target layout should be:

read_word0[N], read_word1[N], read_word2[N], read_word3[N]
ref_word0[N],  ref_word1[N],  ref_word2[N],  ref_word3[N]
candidate_meta[N]

not:

Vec<Candidate { read_id, transcript_id, pos, ... }>

C. SIMD quality/trim masks

Precompute masks once per read:

valid_base_mask
low_quality_mask
trimmed_scoring_len
polyA/polyT trim point
TSO trim point

Then scoring becomes:

mismatches &= valid_scoring_mask

rather than rescanning qualities and suffixes per candidate.

D. SIMD seed-vote reduction only if data are sorted

Candidate-vote reduction is potentially vectorizable only after sorting by candidate key. If you have arrays like:

candidate_key[N]
seed_weight[N]

then SIMD can help compare adjacent keys and reduce weights, but this is secondary. Sorting and memory traffic probably dominate.

5. Best prefetching opportunities

Manual prefetching only helps when you know the next memory address early enough and the access pattern is predictable. It will not fix random lookup by itself. Use it in sorted, streaming phases.

A. Posting-list prefetch in seed-batched retrieval

Seed-batched retrieval sorts query seed occurrences by k-mer and loads each selected posting list once per batch. That is exactly where prefetch can work.  ￼

When iterating sorted QuerySeed groups:

current_seed → current posting list
next_seed    → next posting list metadata
next_next    → maybe posting-list start address

Prefetch:

next KmerEntry
next posting-list first cache line
maybe next posting-list block metadata

This is much more plausible than prefetching in read-at-a-time mode.

B. Reference-window prefetch in candidate-locus buckets

After candidate hits are sorted/regrouped by transcript/locus, the next few buckets are known. Prefetch:

next transcript metadata
next reference sequence window
next read-store sequence words

This is a clean fit because the project’s original hypothesis is candidate-locus regrouping for cache-resident scoring blocks.  ￼

C. Avoid prefetching tiny or high-frequency posting lists blindly

For very small lists, prefetch overhead can exceed benefit. For giant repetitive lists, prefetching may just pollute cache. Gate it:

if postings_len between small_min and large_max:
    prefetch

Something like:

prefetch only 64 <= postings_len <= 4096

then tune empirically.

D. Use block metadata first

A better long-term prefetch strategy is to split posting lists into blocks:

PostingBlockMeta {
    start,
    len,
    max_gene_score,
    min/max gene id or target class,
    maybe block-level target-class counts
}

Then prefetch/decode one block at a time. This also enables block-skipping.

6. Best low-level caching / memory-layout opportunities

A. Target-class priority reduces cache pressure

The fastest candidate is the one you never score. If an exonic hit confidently explains the read, do not touch gene-body/intron reference windows. This reduces irregular reference access and candidate sorting volume.

This is probably a bigger win than micro-optimizing Hamming.

B. RC fallback instead of always-on RC

The README shows the tradeoff clearly:

forward only:   68.03% gene-countable, 20.41% unmapped, 172.4k reads/s
RC enabled:     73.98% gene-countable, 15.52% unmapped, 165.0k reads/s
RC + TSO trim:  76.30% gene-countable, 12.68% unmapped, 149.5k reads/s

RC is recall-positive but costs throughput.  ￼

Make RC conditional:

forward pass first
RC only for unmapped / weak / antisense-only reads

This should recover most RC benefit while reducing lookup and scoring work.

C. Precompute scoring views per read

Do not recompute these per candidate:

reverse complement
TSO trim
polyA/polyT trim
low-quality tail trim
scoring length
quality mask

Create:

PreparedReadView {
    fwd_seq_words,
    fwd_qual_mask,
    fwd_scoring_start,
    fwd_scoring_len,
    rc_seq_words optional,
    rc_qual_mask optional,
    trim_flags,
}

Then candidate scoring only loads the prepared view.

D. Deduplicate identical R2 sequences

The project brief already flags exact or near-exact read deduplication as a later extension: many scRNA-seq R2 sequences repeat, and one could map unique R2 sequences once and propagate assignments back to CB/UMI records while preserving metadata.  ￼

This may be very powerful for abundant genes. Start with exact R2 dedup per chunk:

hash trimmed R2 sequence
map unique sequence once
copy assignment to all read IDs

Caveat: preserve quality-dependent behavior. For a first pass, only dedup when sequence and relevant quality-mask class match.

E. Reduce allocation churn in candidate votes

Seed-batched retrieval can explode candidate votes. The implementation plan explicitly warns that candidate votes can explode and says to track candidate_votes_per_read and skipped large groups.  ￼

Use reusable per-thread arenas:

query_seeds.clear()
candidate_votes.clear()
candidate_hits.clear()
score_blocks.clear()

Avoid allocating new Vecs per read, per seed group, or per bucket.

F. Replace comparison sort with radix sort for fixed-width keys

Many of your hot keys are fixed-width integers:

kmer_code
candidate_key = transcript_id | pos_bin | strand
read_id
gene_id

Rayon comparison sort is convenient, but radix sort can be much faster and more cache-predictable for u64/u128 packed keys. This is worth testing for:

QuerySeed sort by kmer_code
CandidateVote sort by read_id/candidate
CandidateHit sort by candidate locus

Do not rewrite everything at once; benchmark one sort stage.

7. Irregular memory access: where to attack it

The main irregular-access path is probably still:

selected seed → k-mer entry lookup → posting-list slice → candidate vote emission → reference window lookup

The current architecture has already improved this by seed-batching and compact index layout. The next steps are:

A. Convert lookup from “random per seed” to “ordered seed group stream”

You already have seed-batched retrieval, but I would verify with metrics:

posting_list_loads_per_read
posting_list_reuse_factor
distinct_seed_lists_loaded
candidate_votes_per_read

The implementation plan says these are the right metrics to watch.  ￼

If posting-list reuse is low, the read sketch/batch grouping is not co-locating similar seeds well enough.

B. Shard by target class / gene block

Instead of one global stream of candidates, group work into:

exon transcript shard
intron/gene-body shard
high-frequency/repetitive shard
mitochondrial/ribosomal shard

Then each shard has more predictable access patterns and different pruning thresholds.

C. Hot-gene cache

In PBMC/scRNA-seq, expression is skewed. A small number of genes will be responsible for many reads. Track top genes/loci in a run and cache their reference windows / candidate metadata more aggressively.

A simple version:

within chunk, count candidate genes
promote top N genes to hot set
process hot genes in larger locus batches

This is a search-engine-style “hot documents” tier.

D. Blocked posting lists

Long-term, split posting lists into cache-line-aligned blocks:

KmerEntry -> [PostingBlockMeta]
PostingBlockMeta -> compressed postings

Block metadata enables:

skip blocks by target class
skip blocks by gene/WAND bound
prefetch next block
decode block into L1-sized scratch

This is the best direction for reducing random memory traffic after recall semantics stabilize.

8. Where not to spend time yet

I would not prioritize:

hand AVX-512 kernels
HNSW / dense vector read embeddings
expression priors as final evidence
full splice-aware genomic alignment
global WAND aggressiveness tuning before rescue/diagnostics

The repo’s own plan says SIMD is not the bottleneck until recall and candidate generation are understood, and warns against aggressive early stopping before diagnostics.  ￼ That remains correct.

9. Practical next milestone

I would set the next milestone as:

“Constellation 0.2: biologically comparable assignment semantics”

Acceptance targets:

unmapped_rate < 10% on PBMC 100k
gene_countable_rate > 80%
multi_gene_ambiguous does not increase materially
target-class breakdown reported
feature-whitelist breakdown reported
RC fallback implemented
post-score rescue implemented
exon-first/intron-second scoring implemented

Then set the performance target:

recover throughput from 149.5k reads/s back toward 170–200k reads/s
without losing the improved recall

The best performance work for that milestone is:

1. RC fallback instead of always-on RC.
2. target-class priority to avoid scoring unnecessary gene-body/intron hits.
3. precomputed read scoring views.
4. two-stage SIMD/scalar scoring.
5. radix sorting or arena reuse in seed/candidate vote paths.
6. measured prefetching in seed-batched posting-list traversal.

Bottom line

Your instinct is right: do not chase throughput first. At this stage, raw throughput optimization can easily lock in the wrong reference semantics or pruning behavior.

But you do not need to stop all performance work. Focus on performance changes that also clarify or preserve mapping quality:

Recall/semantics first:
  target class priority
  feature whitelist mode
  post-score rescue
  exon+introns-only/splici reference
  barcode/UMI correction
Performance next, but compatible:
  RC fallback
  precomputed scoring views
  two-stage SIMD scoring
  seed/posting prefetch
  candidate-vote arena reuse
  radix sort for packed keys
  blocked posting lists

The deepest efficiency win is still the original thesis: make the computation more ordered before trying to make each operation faster. SIMD and prefetching will work best once the seed/posting/candidate/reference accesses are already sorted into predictable streams.
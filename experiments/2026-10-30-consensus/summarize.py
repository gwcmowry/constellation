"""Render the report from saved measurements; no benchmark processes are run."""
import json
from pathlib import Path
import statistics

HERE=Path(__file__).parent
RAW=HERE/"raw"
TOOLS=["constellation","constellation-consensus","piscem-forward","kallisto-forward"]
NAMES=["Constellation weighted (original)","Constellation consensus (new)","piscem 0.24.0 / rs 0.11.1","kallisto 0.52.0"]
timing=json.loads((RAW/"timing-summary.json").read_text())
truth={c:json.loads((RAW/f"{c}-accuracy.json").read_text()) for c in ["heldout","heldout150"]}
real=json.loads((RAW/"lane2-concordance.json").read_text())
perf=[]
for tool,name in zip(TOOLS,NAMES):
    t=timing["pbmc10m/"+tool]
    perf.append(f"| {name} | {t['median_seconds']:.3f} | {t['million_pairs_per_second']:.3f} | {t['median_rss_gib']:.2f} | {t['median_cpu_seconds']:.2f} |")
accuracy={}
for corpus in truth:
    rows=[]
    for tool,name in zip(TOOLS,NAMES):
        d=truth[corpus][tool];a=d["positive"]
        rows.append(f"| {name} | {100*a['unique_precision']:.4f}% | {100*a['unique_recall']:.3f}% | {100*a['candidate_set_recall']:.3f}% | {100*d['intergenic/exact']['assignment_rate']:.3f}% |")
    accuracy[corpus]="\n".join(rows)
ka=[json.loads((RAW/f"pbmc10m-kallisto-forward-{i}.run.json").read_text()) for i in [1,2,3]]
ka_loop=statistics.median(d["pseudoalignment_loop_seconds"] for d in ka)
ka_before=statistics.median(d["before_pseudoalignment_seconds"] for d in ka)
ka_after=statistics.median(d["after_pseudoalignment_seconds"] for d in ka)
speedup=timing['pbmc10m/piscem-forward']['median_seconds']/timing['pbmc10m/constellation-consensus']['median_seconds']
old=real['read_concordance']['constellation_vs_piscem-forward']
new=real['read_concordance']['constellation-consensus_vs_piscem-forward']
counts=real['common_count_policy']['comparisons']
text=f"""# Competitor comparison and consensus assignment

The new optional consensus policy brings Constellation close to piscem/kallisto
unique-assignment precision on this suite while retaining a **{speedup:.2f}x
whole-process throughput advantage over piscem**. It is not an unconditional
improvement: it is slower than the original weighted policy and loses some
91bp recall. At 150bp its correct unique recall slightly exceeds both peers.
These results support a useful fast, conservative mapping option, not a universal
SOTA claim.

## Throughput and memory

10 million original PBMC L001 FASTQ pairs, gzip level 1, the same Ensembl 93 splici
FASTA/t2g and k=31; CPUs 0–11 (12 physical cores), 12 workers/execution slots,
serial runs, warmup followed by three measured repetitions with rotating tool
order. Median whole-process time includes index loading, FASTQ decoding, mapping,
native output and teardown. Index construction, barcode correction, UMI counting,
and matrix generation are excluded. Constellation writes uncompressed CSTRAD,
piscem writes native RAD, and kallisto writes BUS; their metadata differ.

| Tool/policy | Total seconds | Million pairs/s | Peak RSS GiB | CPU seconds |
|---|---:|---:|---:|---:|
{chr(10).join(perf)}

Kallisto spends a median **{ka_before:.2f}s before its read-processing loop,
{ka_loop:.2f}s in that loop, and {ka_after:.2f}s afterwards** with this unusually
large splici index. Its loop processes {10/ka_loop:.3f}M pairs/s. Therefore its
whole-process ratio is not a pure mapping-algorithm speed ratio. Loop markers
were timestamped externally and checked against the release's ProcessBUSReads
source; raw events are saved. Other stage boundaries differ slightly.

Full ranges, CPU times, memory and page faults: [timing-summary.json](raw/timing-summary.json).
Every run has its command, executable SHA-256, environment, process timer,
native metrics, stderr and `/usr/bin/time -v` record under `raw/`.

## Independent controlled truth

The generator reconstructs sequences from the genome and GTF without using
mapper code. Each holdout has 90,000 newly sampled source windows, balanced
between mature transcripts, splice junctions, and introns. Each window produces
four reads: exact, 0.5% substitutions, 2% substitutions, and one internal 1bp
indel. Thus each has 360,000 positive reads plus 30,000 intergenic genomic and
30,000 random-DNA controls. Both holdouts were generated after freezing the
two-seed consensus rule (seeds 20261030 and 20261031); they were not used to tune
it. Source transcripts and intronic windows were checked against the shared
reference. These are separately sampled datasets, not a new biological sample
or an external standardized benchmark. The error variants share source windows.

Precision = correct unique source-gene calls / unique calls on positive reads.
Recall = correct unique source-gene calls / all positive reads. Candidate recall
allows a multi-gene result containing the source gene. Intergenic assignment is
the fraction of negative-origin reads assigned any gene. Shared sequences make
some source origins unidentifiable; preserving ambiguity is appropriate. These
are controlled source-origin scores, not estimated PBMC biological error rates.

### 91bp holdout

| Tool/policy | Unique precision | Correct unique recall | Candidate recall | Intergenic assigned |
|---|---:|---:|---:|---:|
{accuracy['heldout']}

### 150bp holdout

| Tool/policy | Unique precision | Correct unique recall | Candidate recall | Intergenic assigned |
|---|---:|---:|---:|---:|
{accuracy['heldout150']}

All tools assigned zero random-DNA controls. Per-class and error-condition
results are in [91bp accuracy](raw/heldout-accuracy.json) and
[150bp accuracy](raw/heldout150-accuracy.json). Their manifests record exact
inputs, seeds, sizes, hashes and source-validation counts.

## Correcting orientation and validating joins

The initial October 27 comparison allowed both competitor orientations while
Constellation searched forward reference sequences. That overstated real-read
disagreement and distorted unique-recall comparisons. This report supersedes it:
use forward piscem hits and `kallisto bus --fr-stranded`. Piscem still writes its
normal orientation-tagged RAD; filtering is part of evaluation, not its timer.
The earlier unstranded files remain archived explicitly for audit.

Only accuracy inputs replace R1 with a fixed valid barcode and a unique encoded
UMI. R2 sequence and quality remain unchanged for real reads. This permits exact
read joins: piscem's printed record ordinal is not an input read ID. Every BUS
input-read flag is cross-checked against its UMI, and all decoded counts are
checked against native tool summaries. The RAD parser was independently checked
with alevin-fry on unchanged record bytes after bridging the v2 header. The
Constellation parser was checked against 100,000 TSV records. Timed PBMC inputs
retain original barcodes, UMIs and qualities.

## Real reads and downstream sanity check

One million reads from each of PBMC lanes L001 and L002 were compared. These are
two sequencing lanes of the same biological library, not independent biological
replicates. Only reads whose original 28-base R1 barcode/UMI are valid enter the
pairwise results. There is no independent truth for these real reads.

On L002, agreement with piscem among reads both call uniquely rises from
**{100*old['unique_gene_agreement']:.4f}% to {100*new['unique_gene_agreement']:.4f}%**
({new['both_unique']:,} joint unique calls for the new policy).
Under the same simple counting rule, molecule-set Jaccard rises from
**{counts['constellation_vs_piscem-forward']['molecule_jaccard']:.4f} to
{counts['constellation-consensus_vs_piscem-forward']['molecule_jaccard']:.4f}**,
and log1p total gene-count Pearson correlation rises from
**{counts['constellation_vs_piscem-forward']['gene_total_log1p']['pearson']:.4f} to
{counts['constellation-consensus_vs_piscem-forward']['gene_total_log1p']['pearson']:.4f}**.
L001 shows the same direction. These comparisons do not establish correctness
for every additional or missing molecule.

The common counting rule uses original valid CB/UMIs, barcodes observed in at
least 20 input reads, unique genes, and exact `(barcode, UMI, gene)` deduplication.
It does not perform cell calling, barcode/UMI correction or ambiguous-molecule
resolution. Correlations conditioned on the union of nonzero matrix entries are
also saved, but are sensitive to the many count-one entries and are not treated
as biological validation. See [L001](raw/real-concordance.json) and
[L002](raw/lane2-concordance.json).

## Implemented change

`map --gene-ec-consensus` intersects the gene sets of all usable distinct matched
seeds and requires at least two different k-mers. Missing seeds are tolerated;
contradictory matched seeds reject the call. Shared evidence stays ambiguous;
exon/intron class support no longer forces a single gene. Repeated copies of one
k-mer cannot masquerade as independent supporting sequences. Reusable worker
arrays avoid per-read gene-map allocation. Existing indexes and packed input/
output remain usable. The score in this mode is distinct supporting-seed count.

The original lookup already checked full keys exactly. Its weakness was using
weighted partial evidence and annotation-based tie resolution as a confident
read-level gene decision. Exact seed membership alone does not establish that
the entire read came from that gene.

A first development candidate required three distinct seeds spread across at
least k/2 start positions. It reached 99.977% precision but only 80.068% unique
recall on development positives and was rejected as too restrictive. The retained
two-seed rule reached 99.951% precision and 84.904% recall before holdouts.
The stricter candidate's development summaries are archived under `strict-development-*`.
The original policy remains the default because changing it trades assignment
rate and speed for specificity. Use the flag explicitly when that tradeoff fits
the analysis. Reverse-complement searching, alternative traversal/locality modes,
adaptive rescue and unmapped-diagnostic mode are explicitly rejected with this
option; they do not silently fall back to another decision rule.

## Remaining weaknesses and next direction

Short erroneous reads can have too few intact sampled seeds. Consensus improves
precision by abstaining, so 91bp recall remains below the competitors. Gene-set
agreement also does not enforce one consistent transcript/locus path. The next
useful algorithmic step is targeted sequence/path verification for weak and
conflicting cases, to recover true reads without reverting to forced assignment.
A genome-decoy comparison is also needed: this shared reference deliberately has
no poison/D-list and therefore does not represent every tool's best recommended
reference setup. Intergenic sequences can genuinely share exact sequence with
annotated genes.

The fast path still has room for improvement, including batched lookup in the
new consensus implementation. Native stage metrics and process CPU usage are
saved. cargo-flamegraph was attempted again, but this machine requires a sudo
password for sampling; [the failure is recorded](raw/flamegraph.log). No hardware
counter or flamegraph-based localization is claimed.

This suite does not cover other organisms/assays, GPUs, aligner-level outputs,
native full quantification workflows, or a full precision/recall parameter sweep.
It cannot establish universal SOTA. Kallisto and piscem were verified as the same
latest public releases on October 30; downloaded release metadata is retained.

## Reproduce

Install the two pinned analysis dependencies with UV into the isolated benchmark
venv, and build/test the Rust workspace with `CARGO_BUILD_JOBS=12`,
`RAYON_NUM_THREADS=12`, and `taskset -c 0-11`. The benchmark uses the downloaded
[piscem 0.24.0](https://github.com/COMBINE-lab/piscem/releases/tag/v0.24.0) and
[kallisto 0.52.0](https://github.com/pachterlab/kallisto/releases/tag/v0.52.0)
executables and the shared indexes under `benchdata/comparison-2026-10-27`.
The weighted baseline is the retained `constellation.optimized` executable;
consensus uses `target/release/constellation`. Complete commands and binary hashes
are in each `*.run.json`. The original FASTA, GTF and genome hashes are in the
October 27 input manifest; new holdout hashes are saved here.

```sh
UV_CACHE_DIR=/tmp/constellation-comparison-uv uv pip install \\
  --python benchdata/comparison-2026-10-27/venv/bin/python \\
  -r experiments/2026-10-30-consensus/requirements.txt
taskset -c 0-11 benchdata/comparison-2026-10-27/venv/bin/python \\
  experiments/2026-10-30-consensus/validate_and_benchmark.py
benchdata/comparison-2026-10-27/venv/bin/python \\
  experiments/2026-10-30-consensus/summarize.py
```

Large FASTQs, indexes, per-read Parquet truth/results and binaries are git-ignored.
The generator, native-output decoder, evaluations, commands and raw measurements
are retained. Tests cover evidence intersection, contradictory genes, repeated
seeds, ambiguous annotations, empty/N-containing reads, scratch reuse, sampled
v4/v9 agreement, and incompatible CLI settings, alongside the existing regression
suite. See `raw/tests.log` and the final validation records.
"""
(HERE/"README.md").write_text(text)
print("Report written",flush=True)

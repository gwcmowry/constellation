# Index construction and native counting stage comparison

This experiment compares Constellation, piscem + alevin-fry, and kallisto +
bustools on the same prepared reference and ten-million-pair input. It changes
no production code. Index construction is timed from scratch; counting starts
from the saved native mapping output, with its full correction/grouping pipeline.

## Measured result

Same 10M input pairs, prepared splici reference, twelve-core cap. Build numbers
are one fresh run per tool; counting numbers are medians of three runs after
warmup. All index builds start from the reference FASTA.

| Builder | Build wall time | Peak RSS | Final index size |
|---|---:|---:|---:|
| Constellation direct PtrHash | 24.39s | 14.18 GiB | 4.428 GB |
| piscem 0.24.0 | 190.71s | 7.40 GiB | 2.976 GB |
| kallisto 0.52.0 | 1099.38s | 38.95 GiB | 6.013 GB |

| Post-mapping workflow | Complete count wall time | CPU-seconds | Peak RSS |
|---|---:|---:|---:|
| Our mappings, conservative + UMI edit correction | 1.585s | 2.96 | 0.474 GiB |
| Our mappings, exact-UMI USA | 1.463s | 2.19 | 0.565 GiB |
| Piscem + alevin-fry USA | 1.263s | 9.49 | 0.535 GiB |
| Kallisto + bustools gene-UMI S/U/A | 19.880s | 20.35 | 2.035 GiB |

Native stage breakdowns (medians):

| Workflow | Barcode correction / permit list | Sort / collate | Final molecule/matrix step |
|---|---:|---:|---:|
| alevin-fry | 0.288s | 0.210s | 0.766s |
| bustools | 0.454s | 1.056s | 18.395s |

With **identical piscem evidence and exact-UMI USA results**, our counter takes
1.422s versus alevin-fry's complete 1.263s pipeline.
Our aggregate CPU cost is 2.13 CPU-seconds versus
9.49; alevin-fry has a modest wall-time advantage while
using more CPU in this test. The projection adapter is excluded from this
counter-only control. Default Constellation additionally corrects UMI errors,
so its first row is not an identical-resolution comparison.

The earlier unpinned, twelve-thread builds were 190.88s for piscem and 1,099.60s
for kallisto; the fresh pinned results above are the headline measurements.
Our earlier three-build median was 24.41s. Index construction is a reusable
one-time reference cost; counting is repeated for each sample.

Raw data: [summary.json](raw/summary.json), [validated count outputs](raw/validation.json),
and [binary provenance](raw/tools.json).

Native alevin-fry's default `--small-thresh 100` was also tested in three full
pipelines: median 1.306s, with exactly the same matrices.
This checks that forcing `--small-thresh 0` in the controlled comparison did not
create the main speed ranking. [Default-path control](raw/af-default-control.json).

Bustools' S/U/A count takes **16.914s even on a zero-record BUS file**
with the same EC metadata and split list. Much of its roughly 18.4s final count
time is therefore setup for this EC dictionary, rather than processing the 10M
input pairs. This fixed component makes linear extrapolation to larger inputs
inappropriate. [Empty-input control](raw/bus-empty-split.run.json).



## What the count phase does

Mapping turns each read pair into a cell-barcode/UMI tag and compatible reference
targets or genes, with splicing-class evidence. Counting turns those records into
raw gene-by-cell **molecule counts**:

1. Match cell barcodes to the supplied list, correcting uniquely resolvable
   one-substitution errors and excluding unsupported barcodes.
2. Group read evidence by corrected cell barcode and UMI. Repeated PCR copies
   of a captured molecule can generate several reads with the same tags.
3. Resolve compatible genes and splicing state across those reads. Policies
   differ in how they handle conflicting evidence, cross-gene ties and UMI
   collisions; this is an inference step, not just file-format conversion.
4. Optionally correct UMI sequence errors to avoid counting one original
   molecule more than once. Constellation's default directional correction is
   enabled in its default row; the exact-UMI rows disable it explicitly.
5. Sum resolved molecules into spliced, unspliced and ambiguous gene/cell matrices
   and write axes/metadata. Constellation additionally writes the total and
   ambiguity-component matrices.

For example, five compatible reads with the same corrected cell barcode, UMI
and gene normally contribute one molecule to that cell/gene, not five. Different
UMIs normally contribute separate molecules. All timing boundaries here end at
raw matrices. Cell calling is excluded: every method receives the same published
11,769-barcode list.

## Controls and timing boundaries

- Prepared Ensembl 93 splici FASTA and transcript-to-gene/state annotation, k=31.
  Exactly the same FASTA file is supplied to all builders.
- Constellation uses k31/s20 closed-syncmer sampling and the direct v9 builder.
  Piscem uses its native k31/m19 defaults, with a 24 GiB sort-memory ceiling.
  Kallisto uses native k31 defaults, choosing m23. Its log reports native poly-A
  trimming and replacement of ambiguous reference bases; these were not changed.
- Builds are one fresh, successful run per tool, sequential on CPUs 0–11 with
  twelve threads where supported. This is not a three-repeat index benchmark.
  Previous single-run competitor build measurements are retained for context.
- Count inputs are native output from the original ten million PBMC pairs:
  CSTRAD2, piscem RAD, and kallisto BUS respectively. Mapping is excluded from
  the count timers. Read evidence and numbers of mapped records differ across
  mappers, as documented in the preceding comparison.
- Counting has a warmup and three measured runs per configuration, with rotating
  order. All processes are serial across compared workflows. Native tool thread
  settings are capped at twelve and processes are pinned to the same CPUs.
  Bustools `count` is single-threaded in this release; twelve available cores do
  not make that phase parallel.
- Timings include native dictionary loading, intermediate output and final
  matrix output, process startup and teardown. Reference preparation from
  genome/GTF, downloading tools, compilation and mapping are excluded from these
  stage timings. Normal warm filesystem caching is used; no cache flush/fsync
  durability benchmark is implied.

The indexes do different work. Ours retains 259.5 million sampled forward keys
and gene/class sets; the other two build much denser compacted graphs and retain
reference-target/path information. Piscem reports 1.451 billion canonical k-mers.
Thus the faster build is a useful whole-tool result, not an equal-information
algorithm race. Nor does a faster build imply a smaller index or lower memory.

## Native pipelines and policies

**Constellation default:** `count-splicing` performs barcode correction,
cell partitioning, conservative gene intersection, directional UMI correction,
state resolution and all seven sparse matrices in one process.

**Constellation USA:** `--resolution usa --no-correct-umis` applies the exact-UMI
winner-take-all/spliced-preference policy. Its own-mapping row remains a native
workflow; the separate piscem-evidence row is a control described below.

**Piscem + alevin-fry 0.19.0:**

- `generate-permit-list -d fw --valid-bc ... --cell-bc-correction unique
  --cell-bc-neighborhood hamming-1`
- `collate`
- `quant -r cr-like --small-thresh 0`, native winner-take-all USA model and exact
  UMIs. No molecule-table dump, bootstraps or EM are requested.

The recorded metadata verifies `WinnerTakeAll`, `crlike_umi_edit=0`, and no
molecule dump. All three runs reproduce the previous native matrix values
exactly. The earlier comparison enabled a molecule Parquet dump for validation;
that optional output is excluded from these timing runs.

**Kallisto 0.52.0 + bustools 0.45.1:**

- `correct` against the same barcode list;
- `sort -t 12 -m 2G --no-flags`;
- `count --genecounts --umi-gene --split intronic-targets.txt`.

The split list contains the 73,193 U-class targets from the common reference.
`.mtx` holds spliced/mature counts, `.2.mtx` unspliced/nascent counts, and
`.ambiguous.mtx` the ambiguous layer. UMIs are matched exactly. Its collision
handling and spliced-priority rules differ from both other resolvers; totals
are not expected to be identical. The native source and documented split-output
behavior support this configuration. See the [release notes](https://github.com/BUStools/bustools/releases/tag/v0.43.0),
[spliced-priority update](https://github.com/BUStools/bustools/releases/tag/v0.43.2),
and [v0.45.1 counting implementation](https://github.com/BUStools/bustools/blob/v0.45.1/src/bustools_count.cpp).

## Checks against misleading count comparisons

Timing only alevin-fry `quant` would omit permit-list/correction and collation;
timing only bustools `count` would omit correction and sorting. The primary
comparison includes all of those stages.

A separate bustools control removes `--split` while keeping the same corrected,
sorted BUS input. It takes about 1.47 seconds for the final gene-only count step,
versus about 18.4 seconds for native S/U/A output. Those produce different output
and assignment semantics, so the gene-only number is not used in the main
state-output comparison. Source inspection shows a linear search through the
intronic-target list during split EC preprocessing; this is a plausible contributor
to the difference, not a sampled profile attributing every extra second.

We also streamed bustools `correct` into `sort` (one plus eleven worker slots).
The three upstream runs take about 1.40 seconds and produce a byte-identical
sorted BUS file. Materialized correction plus sorting takes about 1.51 seconds.
Avoiding that intermediate file cannot explain the much longer S/U/A count
stage. Native binaries are unchanged.

The matched-evidence control runs our USA counter on the previously generated
lossless gene/state projection of piscem's mappings. Every normalized S/U/A
coordinate and value agrees with native alevin-fry. Projection-adapter time is
excluded, so this is an isolated counter comparison, not a timed hybrid
piscem-to-Constellation production workflow. The input representations differ:
our input is already gene-projected, while alevin-fry reads native target IDs.

## Output validation and interpretation

`validate_counts.py` normalizes gene/barcode axes and verifies every nonzero
coordinate across all three repetitions. All values are finite positive integers.
The matched piscem-evidence control and all native alevin-fry runs exactly match
the earlier native USA matrix. Normalized matrices are saved in Parquet, with
canonical axes under `raw/`.

| Evidence + resolver | S | U | A, including our diagnostic categories | Total |
|---|---:|---:|---:|---:|
| Ours + default conservative | 2,656,486 | 2,007,141 | 850,393 | 5,514,020 |
| Ours + exact-UMI USA | 2,823,561 | 2,011,271 | 894,749 | 5,729,581 |
| Piscem + alevin-fry / our matched USA control | 2,959,526 | 1,994,063 | 841,163 | 5,794,752 |
| Kallisto + bustools gene-UMI split | 2,807,312 | 2,011,001 | 1,002,480 | 5,820,793 |

The total is not an accuracy ranking. Mapping evidence, cross-gene resolution,
state rules, collision treatment and UMI edit correction all affect it. The
preceding [assignment-policy experiment](../2026-10-31-build-count/README.md)
measures the source-origin tradeoff in choosing additional assignments.

## Recommended larger-data follow-up

Start with the **full PBMC 10k v3 dataset already on disk**: 638,901,019 read
pairs, 11,769 called cells, 47.30 GB of compressed R1/R2 files across two lanes.
The read/cell totals come from the [official summary](https://cf.10xgenomics.com/samples/cell-exp/3.0.0/pbmc_10k_v3/pbmc_10k_v3_web_summary.html);
local file sizes and paths are recorded in [larger-data-plan.json](raw/larger-data-plan.json).
Run 100 million pairs first to measure memory, then all 638.9 million. Preserve
paired reads and original CB/UMI identities, sample across both lanes, and count
both lanes jointly because they contain the same cell library.

For independent validation, use [20k Human PBMCs, 3′ HT v3.1, Chromium X](https://www.10xgenomics.com/datasets/20-k-human-pbm-cs-3-ht-v-3-1-chromium-x-3-1-high-6-1-0):
23,837 cells and 818,107,363 read pairs, verified in its public web summary.
Start with fixed 100M and 500M subsets. Its R1 contains a 16-base barcode and
12-base UMI; R2 is 90 bases. The FASTQs have not been downloaded in this turn.

A subsequent biological challenge could use [10k Adult Mouse Brain Nuclei,
Next GEM v3.1](https://www.10xgenomics.com/datasets/10k-Mouse-Brain-CNIK-3p-nextgem):
9,684 nuclei at approximately 60,000 pairs/nucleus, or roughly 580M pairs from
the published depth. This requires a matched mouse reference and tests a
nuclei/antisense-rich workload beyond PBMCs. Begin with 100M pairs.

Read-volume scaling and cell-number scaling test different things. The full
existing PBMC data deepens the same cell library; the independent 20k library
increases cell diversity/cardinality as well. For a future multi-library atlas,
process libraries independently and prefix cell identities when merging matrices;
the same 16-base barcode can occur in different libraries. Raw FASTQs are needed
for mapper benchmarks; expression-matrix downloads alone cannot test mapping.

More reads do not change the cost of building an index for the same reference.
Keep reference and policies fixed while measuring read-volume scaling. Preserve
native gzip inputs, record decoding/mapping/counting separately, and reuse
scratch outputs across repetitions. The current disk has roughly 125 GB free;
retain only necessary intermediates. Full-data count RAM needs measurement at
100M before committing to 500M/full; the 10M counter is not evidence of constant
memory as read volume grows.

## Reproduction

Use the existing UV-managed environment at
`benchdata/comparison-2026-10-27/venv/bin/python` (NumPy/PyArrow).
`run.py counts` performs warmup plus three measured counting runs;
`validate_counts.py` checks matrices; `controls.py` performs the gene-only and
streamed-BUS controls; `run.py builds` builds every index from scratch.
Do not overlap these CPU workloads. Rebuild to a fresh output directory when
repeating native index construction.

`raw/` retains exact commands, process times, RSS/CPU counters, build-stage event
timestamps, binary hashes, release metadata and complete native logs. Downloaded
tools, sources and larger outputs live in git-ignored `benchdata/`.
Only the newly needed bustools binary/source was downloaded; existing pinned
piscem, kallisto and alevin-fry binaries were reused. The input reads and reference
were unchanged. No production mapper, builder or counter code was edited.

# Detailed gene/splicing-state comparison

Constellation's optional state output largely agrees with piscem and kallisto,
but retains extra ambiguity. The downstream matrices are not interchangeable:
read mapping and molecule-resolution rules both contribute to differences.
This experiment changes no production mapping/counting behavior.

## Methods and controls

- Constellation `--splicing-state`, piscem 0.24.0 (piscem-rs 0.11.1), kallisto
  0.52.0, and the downloaded alevin-fry 0.19.0 release.
- Same Ensembl 93 splici FASTA/t2g, k=31, and forward strand interpretation.
  Piscem target hits are filtered to forward orientation; kallisto's saved runs
  used `--fr-stranded`. Target IDs are projected to `(gene, S/U)` after mapping,
  preserving all compatible states of each gene.
- Read comparison: one million identity-tagged pairs from each PBMC lane, plus
  the existing independent 420k-read 91bp and 150bp simulation holdouts. Original
  R2 sequences/qualities are retained. Only originally valid R1s enter real-read
  agreement metrics. Tagged UMIs identify the original read; BUS `-n` flags provide
  an additional independent identity check.
- Counting comparison: the original ten-million-pair PBMC inputs, the same
  published 11,769-barcode whitelist, and unique Hamming-1 barcode correction.
  Native piscem mappings and our existing full state output are reused.
- The common-counter test converts piscem's final forward target sets to CSTRAD2
  gene/state sets and runs the same `count-splicing` executable on both mappers'
  evidence. This projection retains the target mapper's per-read gene/state
  compatibility; it does not replace its mapping algorithm with ours. Adapter
  record order is not original FASTQ order; it is used only for counting. Native
  scores are unspecified (zero) and ignored by the counter.
- Separately, native alevin-fry GPL/collation/USA quantification runs on the same
  piscem RAD. The primary native configuration is `cr-like`, `WinnerTakeAll`,
  `--small-thresh 0`, with no UMI edit correction. Exact-UMI common-counter controls
  separate this correction difference from the resolution policy itself.
- Processes are run serially on CPUs 0–11. These are functional assignment/count
  comparisons, not a new repeated runtime benchmark. Native molecule-table output
  is enabled for checking the counts.

## Per-read agreement

S/U labels describe compatibility with reference classes, not directly observed
biological maturity. These real-read results measure agreement, not truth.

| Real-data subset | L001 | L002 |
|---|---:|---:|
| Originally valid R1 pairs | 999,029 | 999,158 |
| Same unique gene in both Constellation and piscem | 655,454 | 637,422 |
| Same state in that shared-unique-gene subset | 98.636% | 98.611% |
| Ours A, piscem S | 8,372 | 8,332 |
| Ours A, piscem U | 567 | 519 |
| S versus U reversals in that subset | 0 | 0 |
| Exact full gene/state sets, among reads mapped by both | 95.832% | 95.716% |

Every state disagreement in the shared-unique-gene subset is extra ambiguity in
our output. The same pattern holds when examining every shared candidate gene.
This does not say that every read has an agreed gene: unmapped and differing
gene-set results are retained separately in the raw comparisons.

Kallisto gives essentially the same result: its state agreement with piscem on
shared unique genes is 99.995% on both lanes. Our corresponding agreement with
kallisto is 98.631% and 98.608%. See [read-summary.json](raw/read-summary.json)
for complete denominators, gene-set differences, confusion tables and synthetic
holdout results. Per-read gene/state lists are saved as Parquet in `benchdata/`.

## Independent reference audit and diagnosis

For each read length, 1,500 error-free synthetic reads were selected randomly,
balanced across mature-transcript, junction and intronic origins. Their entire
sequences were searched directly in the reference targets belonging to the known
source gene, separately for S and U. This checks full-sequence reference
compatibility without using any mapper to define truth.

| Balanced audit | 91bp | 150bp |
|---|---:|---:|
| Our exact source-gene state set | 1,484 / 1,500 | 1,494 / 1,500 |
| Our extra state ambiguity | 14 | 5 |
| Source gene absent from our result | 2 | 1 |
| Piscem exact source-gene state set | 1,499 / 1,500 | 1,499 / 1,500 |
| Kallisto exact source-gene state set | 1,500 / 1,500 | 1,500 / 1,500 |

An additional diagnostic sample deliberately selected 200 exact-read
disagreements per length. The full-reference check supported the competitors'
more specific state in all 400. This enriched sample is for diagnosis, not a
population error-rate estimate.

The cause was then checked rather than inferred from agreement:

- **382 cases:** every class-discriminating 31-mer was absent from the sampled
  index because it was not selected by the closed-syncmer rule.
- **6 cases:** at least one discriminating k-mer was selected but rejected by
  the gene-frequency cap; other discriminators can also be unsampled.
- **12 cases:** even the union of all 31-mers for the gene/class is insufficient
  to disambiguate. Those k-mers occur in the extra class across reference targets,
  but no single target contains the full read. Target-path identity matters here.

In total, the read-only Rust index audit checked 1,699 discriminating k-mers:
1,678 unselected and 21 filtered by gene frequency. There were no unexpected
lookup/class inconsistencies or usable discriminators that our scorer ignored.
Thus most sampled disagreements are a loss of discriminating sequence evidence,
not an inverted S/U encoding or an output-parser mistake.

Sources: [reference-audit.json](raw/reference-audit.json),
[sampling-audit-summary.json](raw/sampling-audit-summary.json),
[sampling-audit.log](raw/sampling-audit.log). Audited sequences and cases are
stored in `reference-audit-cases.parquet`; diagnostic selection is explicit.

## Counting under the same rules

Both sets of evidence use our conservative gene intersection and state-aware
directional UMI correction here. Same input pairs, gene/barcode axes and whitelist:

| Final layer | Constellation mappings | Piscem mappings | Our difference |
|---|---:|---:|---:|
| S | 2,656,486 | 2,797,507 | -5.04% |
| U | 2,007,141 | 1,989,631 | +0.88% |
| A, including conflicting/unknown evidence | 850,393 | 802,337 | +5.99% |
| Total | 5,514,020 | 5,589,475 | -1.35% |

These aggregate differences include reads/genes one mapper accepts and the other
does not. On the full sparse gene-by-barcode matrices (including implicit zeros),
Pearson correlations are 0.9884 for S, 0.9588 for U, and 0.9818 for A. Correlation
alone is insufficient: the report also gives L1 discrepancies and weighted
Jaccard overlap. Log1p gene-total correlations across all 58,395 genes are 0.9969, 0.9896 and 0.9937.

## Native alevin-fry USA counts

Native `cr-like` with its default winner-take-all state model produces:

| Layer | Piscem + native alevin-fry |
|---|---:|
| S | 2,959,526 |
| U | 1,994,063 |
| A | 841,163 |
| Total | 5,794,752 |

The native molecule Parquet contains 6,475,696 rows, but **680,944 are unresolved
multi-gene groups**. Only its 5,794,752 assigned rows contribute to the matrix;
this was checked explicitly. Treating every Parquet row as an assigned molecule
would overstate the native total.

Controlled decomposition of the 280,732-total-count difference from our default:

| Change, in sequence | Total | Increase |
|---|---:|---:|
| Our mappings and our counter, UMI correction on | 5,514,020 | — |
| Disable our UMI edit correction | 5,525,754 | 11,734 |
| Use piscem mappings with that same exact-UMI counter | 5,602,921 | 77,167 |
| Use native alevin-fry resolution on the same piscem mappings | 5,794,752 | 191,831 |

The two exact-UMI counters receive the same **7,204,991** accepted piscem read
records and form the same **6,475,696** barcode/UMI groups. Our rule requires a
single gene surviving intersection; native cr-like retains maximally supported
gene/state alternatives and applies its USA resolution rules. The different
resolution policies, rather than input identity or barcode filtering, explain
the final 191,831-count increment. Extra assignments are not automatically more
correct; there is no real-data molecule ground truth here.

We also tested native `--sa-model prefer-ambig` on the same collated input. It
produces S=2,946,628, U=1,993,275, A=961,264, total=5,901,167. This reinforces that
native USA model settings materially affect the result. Both native settings and
the exact-UMI controls are retained; no setting was selected to maximize agreement.

See [count-comparison.json](raw/count-comparison.json) and
[native model comparison](raw/native-ambiguity-model-comparison.json). Every
normalized sparse matrix is also saved as Parquet with canonical gene/cell IDs.

## Implications

The option preserves useful state information and has high agreement, but its
layers should not be treated as interchangeable with native alevin-fry output.
Two distinct next steps are indicated:

1. Verify ambiguous states using additional class-discriminating sequence or a
   compact reference-context/path index. More gene-score tie-breaking cannot
   recover evidence that the sampled index omitted.
2. Offer an explicit native-compatible UMI resolution policy if matching
   alevin-fry matrices is the objective. Keep the current conservative policy
   available; changing resolution silently would change the measurement.

The counter is also a remaining throughput opportunity. The native GPL,
collation and default quantification functional runs took about 0.31s, 0.20s and
1.00s at 12 threads, including the molecule dump. Our earlier state-count run was
7.86s and is largely sequential. These single runs with different policies/output
sets are not a fair algorithm speed ratio, but motivate a controlled study of
cell-wise parallel counting. No production code was changed in this comparison.

## Reproduce and provenance

Use the existing UV-managed Python environment at
`benchdata/comparison-2026-10-27/venv/bin/python` (NumPy and PyArrow). Run serially
under `taskset -c 0-11`:

```sh
python experiments/2026-10-30-state-comparison/compare_reads.py
python experiments/2026-10-30-state-comparison/run_native.py
python experiments/2026-10-30-state-comparison/compare_counts.py
python experiments/2026-10-30-state-comparison/audit_reference.py
```

The native commands, flags, timing logs and release metadata are saved under
`raw/af-*.run.json` and `raw/alevin-fry-release.json`. Run the native stages before
`compare_counts.py`; normalize the published barcode list by removing its `-1`
suffix. The alternative ambiguity-model comparison uses a separate output
directory. `sampling_audit.rs` is a read-only diagnostic linked against the
existing release library; its exact build command is saved in
`raw/sampling-audit-build.json`. Input mappings/reference files come from the
previous frozen comparisons. Large outputs, released binaries/source and
Parquet data are git-ignored under `benchdata/state-comparison-2026-10-30`.

Primary implementation references: [alevin-fry v0.19.0](https://github.com/COMBINE-lab/alevin-fry/releases/tag/v0.19.0),
[released resolution implementation](https://github.com/COMBINE-lab/alevin-fry/blob/v0.19.0/src/pugutils.rs),
and [USA quantification documentation](https://alevin-fry.readthedocs.io/en/latest/quant.html).
For exact options/defaults, the saved release's help, source and `quant.json`
take precedence over moving online documentation.

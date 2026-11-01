# Optional gene-by-splicing-state output

For a complete sample workflow, `constellation quant --index ... --r1 ...
--r2 ... --barcode-whitelist ... --out-dir ... --threads 12` runs both stages
below with state output enabled. Its default is conservative intersection and
directional UMI correction. `quant --resolution usa` selects exact-UMI USA
automatically; use `--keep-assignments` to preserve its CSTRAD2 for re-counting.
See the [README quickstart](../README.md#quickstart) for a runnable example.

`map --splicing-state` is opt-in. It enables forward-strand consensus mapping and
preserves a separate reference-class compatibility state for every candidate
gene. The default weighted mapper and ordinary `--gene-ec-consensus` mode retain
their existing formats and behavior. No index rebuild is needed when the index
was built with complete spliced/intronic annotations, such as the existing splici
reference with its three-column transcript-to-gene table.

```sh
constellation map --threads 12 \
  --splicing-state \
  --index /path/to/reference.ptrhash.ecidx \
  --r1 /path/to/R1.fastq.gz --r2 /path/to/R2.fastq.gz \
  --out reads.cstrad2

constellation count-splicing --threads 12 \
  --assignments reads.cstrad2 \
  --barcode-whitelist /path/to/accepted_barcodes.tsv.gz \
  --out-prefix counts --emit-metrics counts.metrics.json
```

The counting command reads the gene dictionary from CSTRAD2; it does not need an
index. Plain, gzip and zstd CSTRAD2 are accepted. Mapping can write zstd directly
with `--output-compression zstd`; zstd input/output requires the `zstd` executable.
`--output-format tsv` instead adds a `gene_states` column with entries such as
`12:S;47:A`. This is an inspection format; `count-splicing` requires CSTRAD2.

## Meaning of the states

| Label | Stored code | Meaning for a candidate gene |
|---|---:|---|
| S | 1 | Spliced-reference class remains compatible |
| U | 2 | Intronic/unspliced-reference class remains compatible |
| A | 3 | Both classes remain compatible |
| C | 4 | Matched evidence has no common S/U class |
| X | 8 | Stored reference annotations provide no S/U constraint |

Only genes passing the ordinary two-distinct-seed consensus rule are returned.
For each gene, the S/U class sets are intersected across its matched seeds. An
unknown class is uninformative. A conflict remains a conflict; OR-ing S and U
into an apparently compatible A would discard this distinction. C/X candidates
remain gene candidates, so adding state output does not change consensus gene
assignments. Neither class annotation nor splicing state breaks cross-gene ties.

These are **reference-class compatibility labels**, not proof that a physical
molecule is mature or unspliced. Splici intron targets include exon flanks, and
some exonic sequence can occur in both classes. Index construction/read length
and incomplete annotation affect the interpretation. Existing unannotated
indexes cannot recover missing annotations; X marks missing stored constraints.
This remains gene-level seed compatibility, not full transcript/path alignment
or isoform quantification. The default conservative counting policy differs
from alevin-fry's native USA resolution; the optional USA policy below was
validated on matched evidence. Neither is a complete RNA-velocity workflow.

## Resolving molecules

1. Apply the existing unique one-mismatch barcode correction when a whitelist is
   supplied. `--no-correct-barcodes` keeps only exact whitelist members. Without a
   whitelist, observed valid barcodes are used; this command does not call cells.
2. Group mapped reads by corrected barcode and exact UMI. Intersect candidate
   genes across reads. A group with exactly one surviving gene is retained;
   unresolved or conflicting gene groups are excluded and counted in metrics.
   A true UMI collision involving different genes can therefore be excluded.
3. Intersect the selected gene's state constraints across the grouped reads.
   A plus S resolves to S; A plus U resolves to U; S plus U becomes C.
4. Within each cell/gene, optionally merge one-mismatch UMI neighbors using the
   deterministic one-hop directional rule: visit roots in descending read count,
   breaking ties by UMI code, and absorb unassigned neighbors when
   `root_count >= 2 * neighbor_count - 1`. Intersect states across each merge.
   `--no-correct-umis` retains exact UMIs. Correction is performed jointly across
   states, avoiding double-counting one corrected molecule in both S and U.

This is a conservative, explicit resolution policy. It preserves uncertainty
instead of choosing a gene or state from an unresolved tie. The default mapper's
previous counting command is unchanged.

## Optional USA resolution

`count-splicing --resolution usa --no-correct-umis` implements the exact-UMI
`cr-like` / `WinnerTakeAll` USA resolution tested against alevin-fry 0.19.0:

- Each read votes once for every compatible gene/class label. A votes for both
  S and U. Keep the labels with maximal read support within a corrected CB/UMI.
- One winning label produces S or U; tied S/U labels of one gene produce A.
- Cross-gene ties with 2–10 labels prefer the sole spliced gene, if present.
  Include A if that gene's U label also ties. Other cross-gene ties are excluded.
- C/X evidence has no equivalent target label in native USA. Any CB/UMI group
  containing it falls back to conservative intersection, reported in metrics.

This is an explicit biological/reference-class preference, not additional
sequence evidence. In the controlled single-read-molecule source-origin audit,
it improves recall but lowers precision, especially for intronic reads shared
with another gene's spliced reference. `intersection` remains the default.
USA currently requires `--no-correct-umis`; it does not emulate alevin-fry's
optional UMI edit-correction model. The ordinary state-aware directional
correction remains available with the default intersection policy.

The [full experiment](../experiments/2026-10-31-build-count/README.md) checks native
matrix identity on piscem evidence and quantifies the tradeoff on our evidence.
Counting now partitions packed evidence by corrected cell, resolves/corrects
cells with Rayon, and writes deterministic sorted matrices. Use `--threads N`
or `RAYON_NUM_THREADS` to limit worker parallelism. The default conservative matrices
remain byte-identical on the ten-million-pair regression input.

## Count outputs

All matrices use the same `_features.tsv` and `_barcodes.tsv` axes:

- `_spliced.mtx`: S molecules.
- `_unspliced.mtx`: U molecules.
- `_ambiguous.mtx`: A, C and X molecules together.
- `_matrix.mtx`: the total, exactly spliced + unspliced + ambiguous.
- `_compatible_both.mtx`, `_conflicting.mtx`, `_unknown.mtx`: separate components
  of the ambiguous matrix, for inspection or alternative downstream filtering.

Each accepted, corrected molecule contributes once to the total. Metrics report
excluded gene groups, barcode correction, pre/post-correction molecule counts,
and the five state totals. Conflicting/unknown counts remain separately available
per gene and barcode; their inclusion in the aggregate ambiguous matrix does not
mean they were compatible with both reference classes.

## Binary layout

CSTRAD2 is a new self-contained format; existing CSTRAD1 consumers must reject it
unless updated. All integers are little-endian. The header contains:

- Eight bytes `CSTRAD2\0`.
- Five u32 values: version 2, fixed record width 28, implicit ordered-read-ID flag
  1, candidate entry width 5, and number of genes.
- For each gene, a u16 UTF-8 byte length followed by the name. Its ordinal is the
  gene ID. A dictionary is required even for an empty input stream.

The fixed 28-byte record is the CSTRAD1 layout: packed barcode/UMI u64, valid-base
mask u32, primary gene u32 (`u32::MAX` when absent), candidate count u32, score
u16, class/QC flags u16, assignment type u8, reserved zero u8, and returned-gene
count u16. Each returned gene then occupies five bytes: gene ID u32 and one
state code u8. Candidate genes are sorted and unique. The score counts distinct
matched supporting seeds. Original sequence/quality strings are not stored.
Barcode and UMI lengths are fixed at 16 and 12 bases. Bases use A=0, C=1, G=2,
T=3, with the first barcode base in the least-significant two bits and the UMI
following the barcode; valid-mask bit i describes base i.

See the [measurement and validation record](../experiments/2026-10-30-splicing/README.md)
for the cost of this option and the regression checks.

The [detailed comparison with piscem/kallisto and native alevin-fry USA](../experiments/2026-10-30-state-comparison/README.md)
finds about 98.6% state agreement when both tools uniquely call the same gene on
two PBMC lanes. Our remaining disagreements in that subset retain extra ambiguity.
Sparse sampling and molecule-resolution policy both affect the final count layers;
they are not interchangeable with native alevin-fry output.

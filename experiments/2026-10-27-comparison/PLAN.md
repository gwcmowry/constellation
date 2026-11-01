# Matched mapper comparison

Evaluation protocol written before collecting comparison results. Production
mapping/scoring code and parameters are frozen at the retained packed-pipeline
implementation. This experiment does not tune accuracy on the evaluation reads.

## Scope and controls

- Current downloaded public releases: piscem 0.24.0 and kallisto 0.52.0. Save
  release metadata, executable checksums, CLI help, commands and raw logs.
- Identical Ensembl 93 splici FASTA and transcript-to-gene table. All use k=31;
  Constellation retains its existing closed-syncmer s=20 index. Native index
  structures differ by design. Reuse the existing piscem index if compatible.
- Serial processes restricted to CPUs 0–11 and 12 workers/execution slots.
  Warmup then three measured runs. Report whole-process wall time and RSS,
  including index loading and native mapping output. Index construction is
  separate. Output is Constellation RAD, piscem RAD, or kallisto BUS; these
  formats carry different metadata. No barcode correction, cell calling or
  UMI-resolution time is included in mapping throughput.
- Throughput corpus: the same first 10M PBMC L001 pairs used for optimization,
  gzip level 1, with original barcodes/UMIs. Report compressed input sizes.
- Use common unstranded candidate-gene interpretation for read-level comparison.
  Also retain Constellation's stricter countable assignment labels.

## Accuracy and concordance

Create a separate deterministic simulator from the genome and GTF, without
calling Constellation's simulator or sampling its index. Select mature transcript,
junction-enriched, and intronic source sequences and emit matched 91bp reads with
no errors, 0.5% substitutions, 2% substitutions, or one internal 1bp indel. Sample
30,000 independent source windows per class and apply each error condition to
each window. Sample 30,000 intergenic genomic and 30,000 random DNA negative
controls. Save source gene, coordinates, conditions and read IDs in Parquet.
Sampling is broad across annotated transcripts, not an expression-matched PBMC
model. Intergenic-origin sequences can have homologs in annotated genes.

Unique-gene precision is correct unique assignments / all unique assignments
on positive reads; unique-gene recall is correct unique assignments / positive
reads. Also report wrong unique assignments / positive reads, candidate-set
recall (source gene appears in the returned set), ambiguity, and negative-control
assignment rates. Identical sequences shared between genes make source-gene truth
partly unidentifiable; report ambiguous behavior instead of treating assignment
rate as accuracy. Error variants of a source window are correlated, not independent
replicates. Synthetic negatives are summarized separately, not mixed into an
arbitrary overall precision score.

Real-data comparison: first 1M PBMC pairs, with unchanged R2 sequence/qualities
and artificial valid R1 identity tags. Encode read ordinal in the 12bp UMI so
piscem RAD can be joined reliably: its printed ID is a mapped-record ordinal,
not an input-read ID. This is an accuracy-only corpus, never a timing input.
Record original R1 validity and report agreement/discordance, not real-data
accuracy. Collapse transcript candidates to Ensembl gene sets using the common
t2g table. Do not equate correlation of total gene counts with per-read truth.

As a downstream sanity check, restore original real-read CB/UMIs and apply one
common counting policy to all outputs: valid barcodes observed in at least 20
input reads, unique-gene assignments, and exact `(barcode, UMI, gene)` duplicate
removal. Report molecule overlap and sparse gene-by-barcode count concordance.
This threshold is mapper-independent, not a biological cell caller; this check
does not validate barcode/UMI correction, ambiguous molecule resolution or the
tools' full downstream quantification pipelines.

## Interpretation

This is a controlled single-machine comparison, not a field-wide SOTA proof.
The shared reference has no genome decoy/D-list; that choice isolates mapper
behavior but is not every tool's best recommended reference configuration.
Downstream biological equivalence remains separate from raw mapping throughput.
If results expose quality losses, report them even if speed improves.

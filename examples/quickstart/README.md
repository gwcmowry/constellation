# Tiny runnable example

These synthetic reads exercise the installation and the full gene/state path;
they are not a biological accuracy or throughput benchmark. The
[main quickstart](../../README.md#try-the-bundled-example) has copyable commands.

The reference has two genes and five targets. One sequence occurs in both the
spliced and intronic classes of GENE_A. Seven read pairs cover two cell barcodes;
two reads are PCR duplicates sharing a barcode, UMI and gene.

Expected default/USA results: **6 molecules = 4 spliced + 1 unspliced + 1
ambiguous**, across two genes and two cells. The CLI integration test checks
these totals and compares all matrix/axis files to the separate native stages.

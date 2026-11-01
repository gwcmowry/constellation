# Initial competitor comparison

Use the [corrected comparison and tested consensus implementation](../2026-10-30-consensus/README.md)
for conclusions. These scripts, frozen inputs, and raw runs are retained for
reproduction and audit.

The initial piscem and kallisto results allowed both mapping orientations,
whereas the Constellation index/query used forward reference orientation. That
made the initial real-read agreement and unique-gene recall comparisons
inappropriate for ranking accuracy. The corrected experiment uses forward
piscem hits and `kallisto bus --fr-stranded`, and separately records native
mapping time. It keeps the original observations rather than silently replacing
the original unstranded results.

The simulator reconstructs transcripts and introns from the genome and GTF,
checks them against the shared reference, and does not import mapper code.
`read_outputs.py` handles Constellation's primary decision separately from its
broader score-tie list. Its piscem reader supports the release's RAD v2 format;
`validate_parsers.py` checks unchanged record bytes with alevin-fry's independent
legacy reader after bridging the header. Every validation BUS read number is
checked against its identity-tagged UMI. Large inputs and per-read outputs are
stored in git-ignored `benchdata/comparison-2026-10-27`.

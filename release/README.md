# Constellation @VERSION@

Portable Linux x86_64 executable, statically linked with musl and built for the
baseline x86-64 CPU target. Rust, Python, glibc and a GPU are not needed to run
this binary. Plain and gzip input are built in; optional zstd streams require
the `zstd` executable. Reference preparation and cell calling remain separate.

## Install

Verify the downloaded archive before extracting it:

```sh
sha256sum -c @ARCHIVE@.sha256
tar -xzf @ARCHIVE@
cd @DIRECTORY@
sh install.sh
```

The installer verifies the extracted files and installs into `~/.local/bin`.
Add that directory to your `PATH`, or run `~/.local/bin/constellation` directly.
Use `sh install.sh --prefix /your/prefix` for a different location. Installing
replaces the existing `constellation` executable at that destination.
You can also run `./constellation` directly without installing.

## Try the bundled example

From this extracted directory:

```sh
mkdir -p demo-output
./constellation index-ec --threads 2 \
  --transcripts examples/quickstart/reference.fa \
  --t2g-map examples/quickstart/t2g.tsv \
  --format ptrhash --syncmer-s 20 --out demo-output/reference.ecidx
./constellation quant --threads 2 \
  --index demo-output/reference.ecidx \
  --r1 examples/quickstart/R1.fastq --r2 examples/quickstart/R2.fastq \
  --barcode-whitelist examples/quickstart/barcodes.txt \
  --out-dir demo-output/sample
```

Expected: 7 read pairs, 2 genes, 2 cells, **6 molecules** (4 spliced, 1 unspliced,
1 ambiguous). Select a new output directory on a repeat run. `quant --resolution
usa` enables the optional exact-UMI USA resolver; conservative intersection with
directional UMI correction is the default.

R1 and R2 are the two mates of paired-end sequencing. For the supported 10x 3′
v3 layout, R1 contains the 16-base cell barcode and 12-base UMI, and R2 is the
RNA-derived sequence to map. These roles are specific to the library protocol.

## Provenance and notices

`BUILD-INFO.json` records the source commit, source fingerprint, dirty-tree flag,
toolchain, target, binary hash and smoke-test result. A filename containing
`-dev-` is a local working-tree snapshot, not a tagged release.
`THIRD-PARTY.json` inventories build/runtime dependencies; license/copyright
files and applicable source copies are under `third-party/`.
Constellation is licensed under MIT OR Apache-2.0, as declared in its Cargo
manifests. The notices for bundled dependencies retain their respective terms.

Full documentation: https://github.com/gwcmowry/constellation

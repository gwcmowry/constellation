# Optional splicing-state output and counting

Implemented opt-in `map --splicing-state` and a corresponding `count-splicing`
command. The option enables consensus mapping and preserves reference-class
compatibility for each candidate gene. It does not change consensus gene calls.
Existing default mapping and ordinary consensus output remain CSTRAD1; state
output is self-contained CSTRAD2 (or extended TSV for inspection).

See [usage, state meanings, molecule policy and binary layout](../../docs/splicing_state.md).

## Mapping cost

Same original 10M PBMC L001 FASTQ pairs, prebuilt 4.428GB PtrHash k31/s20 index,
gzip level 1 inputs, CPUs 0–11 and 12 Rayon workers. One warmup and three serial
measurements per mode, rotating order. Whole-process times include loading,
decoding, mapping, native uncompressed output, and exit. Counting is separate.

| Mode | Median seconds | Million pairs/s | Output bytes |
|---|---:|---:|---:|
| Weighted, before | 3.067 | 3.260 | 369,461,060 |
| Weighted, after | 3.083 | 3.244 | 369,461,060 |
| Consensus, before | 4.384 | 2.281 | 334,172,392 |
| Consensus, after | 4.391 | 2.277 | 334,172,392 |
| Consensus with splicing states | 4.562 | 2.192 | 348,708,203 |

State preservation adds **3.90% elapsed mapping time** relative to consensus in
the same build. Default-mode timing changes by 0.5%, within the observed run
variation (before 3.063–3.193s; after 3.022–3.192s). The difference from weighted
mapping includes the previously implemented consensus algorithm; it should not
all be attributed to state annotations.

The optional state code is specialized at compile time and selected once per
worker chunk. Default workers do not allocate or update the state array. On the
middle measured runs, the fused preparation/mapping/formatting stage increased
from 4.184s to 4.368s; input waiting remained below 0.008s and output writing was
about 0.11s. These are coarse stage timings, not isolated kernel measurements.
Splicing-mode RSS was 4,474,044 KiB on that run. Hardware profiling was attempted
with cargo-flamegraph, but sudo required a password; the failure is saved in
[flamegraph.log](raw/flamegraph.log).

Commands, executable hashes, all repetitions and sizes:
[summary.json](raw/summary.json), `raw/*.run.json`, native metrics and time logs.

## Exact regression and output checks

Across all 10M pairs:

- Weighted output before/after is byte-identical, SHA-256
  `b7c5b5ad39e5b2bca42cf253d1d9e10ef9ac5e120f752384681a14fe4cb3b99a`.
- Consensus output before/after is byte-identical, SHA-256
  `3d59110cda44a37dca8d6186d85f3342b124ecb0db06b7b97df31093f0b5a53c`.
- Removing the embedded dictionary and state bytes from CSTRAD2 reproduces every
  CSTRAD1 consensus record exactly, including gene sets, scores and flags.
- The state dictionary contains 58,395 unique gene names. Every gene/state entry
  is checked for a valid ID and state code.

These checks verify added information and compatibility, not biological truth
for splicing state. Full results: [output-checks.json](raw/output-checks.json).

## Real-data counting check

One complete count run on the 10M-pair state output used the existing published
PBMC filtered barcode list, unique one-mismatch barcode correction, conservative
gene intersection within CB/UMI groups, and one-hop directional UMI correction.
It took **7.861s whole-process**, including output matrices, with peak RSS
1,230,888 KiB. This is one functional timing, not a repeated performance study.

| Final molecule class | Count |
|---|---:|
| Spliced | 2,656,486 |
| Unspliced | 2,007,141 |
| Compatible with both | 849,400 |
| Conflicting reference-class evidence | 993 |
| Unknown annotation | 0 |
| Total | 5,514,020 |

The aggregate ambiguous matrix contains 850,393 molecules. Its three components
are also emitted separately. All 11,769 barcode columns share the same gene
axis. Per-entry streaming checks confirm both `total = S + U + ambiguous` and
`ambiguous = compatible_both + conflicting + unknown`; matrix sums also agree
with native metrics.

The policy excluded 1,036,441 unresolved multi-gene CB/UMI groups and 2,495 groups
with conflicting gene evidence. These exclusions, including possible true UMI
collisions, are explicit limitations of the conservative policy. No equivalence
to alevin-fry's native molecule-resolution algorithm is claimed. This is not a
complete biological/RNA-velocity validation. See
[count metrics](raw/count.metrics.json) and [matrix checks](raw/matrix-checks.json).

## Verification and reproduction

The workspace's 86 tests pass. New tests cover state intersection algebra,
cross-gene ambiguity resolved by another read, S/U conflicts within and across
reads, unknown annotation, barcode correction, joint state/UMI correction,
sampled v4/v9 equivalence, extended TSV, plain/gzip/zstd CSTRAD2, corruption,
truncation, empty inputs and matrix conservation. Formatting passes; Clippy
completes with the existing warnings and no new splicing-module warnings.

The before binary was preserved at
`benchdata/comparison-2026-10-27/constellation.before-splicing` before editing.
Build with `CARGO_BUILD_JOBS=12 RAYON_NUM_THREADS=12 taskset -c 0-11 cargo build
--release --workspace`, then run:

```sh
taskset -c 0-11 benchdata/comparison-2026-10-27/venv/bin/python \
  experiments/2026-10-30-splicing/benchmark.py
```

The Python environment was installed with UV during the preceding comparison;
this script itself uses the standard library. Large outputs and matrices are
git-ignored under `benchdata/splicing-2026-10-30`. Tests and profiler-attempt logs
are saved in `raw/`.

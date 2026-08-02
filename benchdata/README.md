# Local benchmark corpus

This directory is intentionally retained outside git. As of 2026-08-02 it
contains approximately 152 GiB:

- `pbmc_10k_v3/` (approximately 93 GiB): FASTQs, assignment outputs, metrics,
  and historical experiment logs.
- `reference/` (approximately 59 GiB): reference inputs and generated indexes.

Do not delete this directory during repository cleanup. It contains the raw
evidence for older experiments. New small metrics and profiler summaries should
go in a dated directory under `experiments/archive/`; large inputs, indexes, and
assignment files should remain here or in `/tmp` and be referenced by a manifest.

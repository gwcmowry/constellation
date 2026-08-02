# Benchmark Plan

Initial measurements should compare input-order and component-bucketed mapping
on the same index and read batch.

Collect:

- reads per second
- mapped, ambiguous, multi-component, and unmapped fragments
- mean seed lookups per read
- mean candidate components per mapped read
- mean pathset intersections per mapped read
- index size

Linux hardware counters:

```bash
perf stat -e cycles,instructions,cache-misses,LLC-load-misses,branches,branch-misses \
  target/release/pathcert map --index index.pcidx --reads reads.fq --out out.ec.tsv
```

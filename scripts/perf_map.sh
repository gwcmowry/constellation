#!/usr/bin/env bash
set -euo pipefail

perf stat \
  -x, \
  -e cycles,instructions,cache-references,cache-misses,branches,branch-misses,L1-dcache-loads,L1-dcache-load-misses,LLC-loads,LLC-load-misses \
  "$@"

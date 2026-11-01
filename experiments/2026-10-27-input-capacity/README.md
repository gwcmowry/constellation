# Input capacity after the packed pipeline changes

The current input pipeline can supply **4.337 million read pairs/s** in isolation
on the existing PBMC fixture. Zlib-rs decoding alone supplies **6.230 million
pairs/s worth of decompressed data**. The current reader is therefore a plausible
next bottleneck, but 4.3M pairs/s is not an inherent limit of the decoder or
hardware. The next experiment should separate large-block decoding from record
parsing and validation, or reduce their cost with a byte-oriented parser.

These are direct release-mode benchmarks of the existing `FastqStreamReader`,
`Slab::fill`, and paired `Pipeline`, without index lookup, mapping or output.
Every case has a full warmup and three serial measurements. CPU affinity was
0–11; individual mate tests use one input thread, and paired tests decode/read
the two mates concurrently. There are no concurrent mapping workers in these
measurements, so their capacity is an isolated upper bound, not guaranteed
end-to-end throughput. The process reported zero filesystem inputs, confirming
that the data was cached. No GPU or hardware-counter profiling was involved.

The fixture has 10M pairs, 28-base R1 and 91-base R2, with gzip level 1 inputs.
R1 is 234,052,182 compressed bytes and 1,154,400,627 raw FASTQ bytes. R2 is
554,509,759 compressed bytes and 2,414,400,627 raw FASTQ bytes. Input hashes and
provenance are in `../2026-10-27-read-throughput/raw/inputs.json`. The decoder
is zlib-rs 0.6.7 through flate2 1.1.9. Rates in GB/s use decimal bytes and include
FASTQ headers, separators and qualities, not just sequence bases.

| Operation on all 10M records or pairs | Median seconds | Rate |
|---|---:|---:|
| R1 gzip decode only, one thread | 0.621 | 16.091M R1 record equivalents/s; 1.858 GB/s |
| R2 gzip decode only, one thread | 1.596 | 6.267M R2 record equivalents/s; 1.513 GB/s |
| Both gzip decoders concurrently | 1.605 | 6.230M pair equivalents/s; 2.223 GB/s combined |
| R1 decode and parse/validate | 1.116 | 8.964M R1 records/s |
| R2 decode and parse/validate | 2.149 | 4.653M R2 records/s |
| R1 decode, parse and pack slab | 1.185 | 8.437M R1 records/s |
| R2 decode, parse and pack slab | 2.261 | 4.423M R2 records/s |
| Full paired gzip input pipeline, including recycling | 2.306 | 4.337M pairs/s |
| R1 raw FASTQ parse and pack | 0.621 | 16.116M R1 records/s |
| R2 raw FASTQ parse and pack | 0.759 | 13.174M R2 records/s |
| Full paired raw FASTQ input pipeline | 0.797 | 12.548M pairs/s |

Decode-only tests consume all bytes through a reusable 1 MiB output buffer and
do not parse/count records; the record-equivalent rates use the known fixture
count, and decoded byte counts are checked against the input manifest. Parsed
tests all verify 10,000,000 records, with paired order/count checks and 65,536-read
batches. Raw controls assume already-decompressed cached files; they exclude
the cost of producing or reading those files from uncached storage.

The paired gzip pipeline varied from 2.299 to 2.313 seconds. Its throughput is
29.2% above the retained gzip-to-RAD+zstd mapper's 3.357M pairs/s, and 19.1% above
the no-output mapper's 3.641M pairs/s. CPU/cache competition with mapping can
reduce that headroom. The previously reported 0.014-second input wait in the
RAD+zstd run measures only consumer starvation; decoding executes concurrently
and takes much longer. That wait metric cannot be used to calculate decoder
throughput.

R2 determines the input rate. Both mates are already processed concurrently,
and R1's full reader can run almost twice as fast as R2's. Within each mate's
thread, however, decode, parse/validate, and copy/uppercase are sequential. The
R2 differences suggest roughly 0.55 seconds beyond inflation for parsing and
validation, another 0.11 seconds for packing, and about 0.05 seconds for paired
pipeline coordination. These differences are informative, not an additive
hardware profile. Increasing buffer/channel depth alone cannot remove the
dominant per-record work.

A useful next design is to let a reader decode large raw FASTQ blocks, identify
complete-record boundaries efficiently, and let mapping workers validate and
consume sequence/quality slices. This can overlap parsing with the next block's
inflation and avoid copying every identifier and separator into temporary
strings. A simpler initial prototype could replace repeated `read_line` calls
with bulk newline scanning while preserving the existing ownership model.
Preserve UTF-8/error behavior, sequence/quality length checks, mate pairing,
CRLF, concatenated gzip streams, CRC errors and truncated-record handling.
An input target of 5–6M pairs/s is supported by the measured decoder capacity,
but is not an implemented or benchmarked result. Mapper CPU time and memory
traffic must be included in any acceptance test.

Further scaling can use independently compressed input shards or blocks, or a
parallel decoder for ordinary gzip. DEFLATE back-references may cross block
boundaries and require up to 32 KiB of prior output, so arbitrary file offsets
are not independent starting points. [RFC 1951](https://www.rfc-editor.org/rfc/rfc1951)
defines that dependency. [BGZF](https://www.htslib.org/doc/bgzip.html) packages
independent gzip blocks with explicit lengths; [rapidgzip](https://github.com/mxmlnkn/rapidgzip)
implements parallel decompression for ordinary gzip. These are candidates for
evaluation, not measured improvements here. Count preprocessing/index creation,
CRC validation, buffering and the CPU cores taken from mapping in total runtime.

At the current fixture's compression ratio, 4.337M pairs/s requires about
342 MB/s of compressed input, and 6.230M pairs/s about 491 MB/s. Cold disk or
network throughput must sustain that rate too. Warm-cache results do not test
that constraint, and other read lengths or compression settings can change
both codec speed and the conversion from bytes/s to pairs/s.

The diagnostic was temporarily added as an ignored release test using the
actual private production reader, then removed. Production code and the normal
test suite are unchanged by this investigation. To reproduce, apply
`input-diagnostic.patch`, then build the release test binary:

```sh
CARGO_BUILD_JOBS=12 RAYON_NUM_THREADS=12 taskset -c 0-11 \
  cargo test --release -p constellation-cli --bin constellation --no-run \
  --message-format=json --offline > experiments/2026-10-27-input-capacity/raw/build.jsonl
UV_CACHE_DIR=/tmp/uv-cache uv run --python 3.12 experiments/2026-10-27-input-capacity/run.py
UV_CACHE_DIR=/tmp/uv-cache uv run --python 3.12 experiments/2026-10-27-input-capacity/summarize.py
```

The runner shares the earlier benchmark lock to prevent overlap and saves the
command, test log, `/usr/bin/time` report and raw samples. Results are in
[measurements.json](raw/measurements.json) and [summary.json](raw/summary.json).

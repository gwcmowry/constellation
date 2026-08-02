# PathCert Architecture

The hot path is:

```text
FASTQ batch -> seed extraction -> sorted global seed lookup
            -> component buckets -> local pathset intersection
            -> certificate ID -> EC count
```

The prototype avoids per-read transcript-list expansion. A seed lookup returns
component-local pathsets. The resolver intersects those pathsets inside each
candidate component and interns the final compatible transcript set as a
certificate.

The first implementation stores transcript sequences directly so correctness
tests can focus on certificate behavior. The public module boundaries are
intended to survive a later structure-of-arrays graph store:

- `graph`: component and transcript-path construction
- `index`: sorted seed index, pathset operations, certificate store
- `mapper`: batch mapping, component bucketing, and EC counts
- `io`: FASTA, FASTQ, and small GTF parsing

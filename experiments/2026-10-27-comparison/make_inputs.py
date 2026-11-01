"""Independent genome/GTF truth and identity-tagged real reads (no mapper imports)."""
from __future__ import annotations

import bisect
import collections
import gzip
import hashlib
import json
from pathlib import Path
import random
import re

import pyarrow as pa
import pyarrow.parquet as pq

ROOT = Path(__file__).resolve().parents[2]
DATA = ROOT / "benchdata/comparison-2026-10-27"
REF = ROOT / "benchdata/reference/ensembl93"
RAW = Path(__file__).parent / "raw"
SEED = 20261027
READ_LENGTH = 91
MASK = 0x5A3C69
CB = "ACGTGTCAGTACGATC"
COMPLEMENT = bytes.maketrans(b"ACGTN", b"TGCAN")


def tag(read_id):
    value = read_id ^ MASK
    return "".join("ACGT"[(value >> shift) & 3] for shift in range(22, -1, -2))


def rc(seq):
    return seq.translate(COMPLEMENT)[::-1]


def write_pair(r1, r2, read_id, sequence, quality=None):
    header = f"@comparison_{read_id}\n".encode()
    r1.write(header + (CB + tag(read_id) + "\n+\n" + "I" * 28 + "\n").encode())
    r2.write(header + sequence + b"\n+\n" + (quality or b"I" * len(sequence)) + b"\n")


def load_annotation():
    transcripts = {}
    genes = collections.defaultdict(list)
    gene_pattern = re.compile(r'gene_id "([^"]+)"')
    tx_pattern = re.compile(r'transcript_id "([^"]+)"')
    with (REF / "Homo_sapiens.GRCh38.93.gtf").open() as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            chrom, _, feature, start, end, _, strand, _, attrs = line.rstrip().split("\t")
            if feature not in ("gene", "exon"):
                continue
            gene = gene_pattern.search(attrs).group(1)
            start, end = int(start) - 1, int(end)
            if feature == "gene":
                genes[chrom].append((start, end))
            else:
                tx = tx_pattern.search(attrs).group(1)
                row = transcripts.setdefault(tx, [gene, chrom, strand, []])
                assert row[:3] == [gene, chrom, strand]
                row[3].append((start, end))
    for row in transcripts.values():
        row[3].sort()
    return transcripts, genes


def load_genome():
    genome = {}
    name, chunks = None, []
    with gzip.open(REF / "Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz", "rb") as handle:
        for line in handle:
            if line.startswith(b">"):
                if name is not None:
                    genome[name] = b"".join(chunks).upper()
                name, chunks = line[1:].split()[0].decode(), []
            else:
                chunks.append(line.strip())
        genome[name] = b"".join(chunks).upper()
    return genome


def synthetic():
    L = READ_LENGTH
    rng = random.Random(SEED)
    transcripts, gene_intervals = load_annotation()
    print(f"annotation: {len(transcripts)} transcripts", flush=True)
    genome = load_genome()
    print(f"genome: {len(genome)} sequences", flush=True)
    t2g = {}
    with (REF / "simpleaf_splici_r91_k31/ref/t2g_3col.tsv").open() as handle:
        for line in handle:
            name, gene, state = line.split()
            if state == "S":
                t2g[name] = gene
    # Restrict positives to transcripts the common reference says it represents.
    txs = {tx: row for tx, row in transcripts.items()
           if tx in t2g and row[1] in genome and sum(b-a for a, b in row[3]) >= (L + 1)}
    assert all(t2g[tx] == row[0] for tx, row in txs.items())
    mature_ids = sorted(txs)
    junction_ids = [tx for tx in mature_ids if len(txs[tx][3]) > 1]
    introns = sorted({(row[0], row[1], a[1], b[0], row[2])
                      for row in txs.values() for a, b in zip(row[3], row[3][1:])
                      if b[0] - a[1] >= (L + 1)})
    cache = {}

    def transcript(tx):
        if tx not in cache:
            _, chrom, strand, exons = txs[tx]
            seq = b"".join(genome[chrom][a:b] for a, b in exons)
            cache[tx] = seq if strand == "+" else rc(seq)
        return cache[tx]

    def source(kind):
        while True:
            if kind == "intronic":
                gene, chrom, start, end, strand = rng.choice(introns)
                position = rng.randrange(start, end - L)
                seq = genome[chrom][position:position + (L + 1)]
                if strand == "-":
                    seq = rc(seq)
                tx, coordinate, junction = "", position, -1
            else:
                tx = rng.choice(mature_ids if kind == "mature" else junction_ids)
                gene, chrom, strand, exons = txs[tx]
                whole = transcript(tx)
                junction = -1
                if kind == "junction":
                    lengths = [b-a for a,b in (exons if strand == "+" else exons[::-1])]
                    junctions, total = [], 0
                    for length in lengths[:-1]:
                        total += length
                        if 20 <= total <= len(whole) - 21:
                            junctions.append(total)
                    if not junctions:
                        continue
                    junction = rng.choice(junctions)
                    lo, hi = max(0, junction - (L - 20)), min(junction - 20, len(whole) - (L + 1))
                    if lo > hi:
                        continue
                    coordinate = rng.randint(lo, hi)
                else:
                    coordinate = rng.randrange(len(whole) - L)
                seq = whole[coordinate:coordinate + (L + 1)]
            if len(seq) == (L + 1) and all(b in b"ACGT" for b in seq):
                return seq, dict(gene=gene, chromosome=chrom, strand=strand,
                                 transcript=tx, source_position=coordinate,
                                 junction_position=junction)

    def perturb(seq, error):
        if error == "indel1":
            position = rng.randint(15, (L - 16))
            if rng.randrange(2):
                return seq[:position] + bytes([rng.choice(b"ACGT")]) + seq[position:(L - 1)], 1
            return seq[:position] + seq[position+1:(L + 1)], 1
        out, changes = bytearray(seq[:L]), 0
        rate = {"exact": 0, "sub005": .005, "sub02": .02}[error]
        for i, base in enumerate(out):
            if rng.random() < rate:
                out[i] = rng.choice([b for b in b"ACGT" if b != base])
                changes += 1
        return bytes(out), changes

    rows = []
    intron_windows = collections.defaultdict(dict)
    with gzip.open(DATA / "synthetic_R1.fastq.gz", "wb", compresslevel=1) as r1, \
         gzip.open(DATA / "synthetic_R2.fastq.gz", "wb", compresslevel=1) as r2:
        source_id = 0
        for kind in ("mature", "junction", "intronic"):
            for _ in range(30_000):
                clean, meta = source(kind)
                if kind == "intronic":
                    intron_windows[meta["gene"]][source_id] = clean[:L]
                for error in ("exact", "sub005", "sub02", "indel1"):
                    seq, changes = perturb(clean, error)
                    read_id = len(rows)
                    assert len(seq) == L
                    write_pair(r1, r2, read_id, seq)
                    rows.append(dict(read_id=read_id, source_id=source_id, source_class=kind,
                                     error=error, changes=changes, **meta))
                source_id += 1
            print(f"generated {kind}", flush=True)

        # Complement of the union of ALL GTF gene intervals on primary chromosomes.
        gaps = []
        for chrom in [str(i) for i in range(1,23)] + ["X", "Y"]:
            end = 0
            for a,b in sorted(gene_intervals[chrom]):
                if a - end >= (L + 1):
                    gaps.append((chrom, end, a))
                end = max(end,b)
            if len(genome[chrom]) - end >= (L + 1):
                gaps.append((chrom, end, len(genome[chrom])))
        cumulative, total = [], 0
        for _, a,b in gaps:
            total += b-a-(L - 1)
            cumulative.append(total)
        for kind in ("intergenic", "random"):
            for _ in range(30_000):
                chrom, coordinate, strand = "", -1, "+"
                if kind == "random":
                    seq = bytes(rng.choices(b"ACGT", k=L))
                else:
                    while True:
                        draw = rng.randrange(total)
                        idx = bisect.bisect_right(cumulative, draw)
                        chrom, a, _ = gaps[idx]
                        coordinate = a + draw - (cumulative[idx-1] if idx else 0)
                        seq = genome[chrom][coordinate:coordinate+L]
                        if all(b in b"ACGT" for b in seq):
                            break
                    if rng.randrange(2):
                        seq, strand = rc(seq), "-"
                read_id = len(rows)
                write_pair(r1, r2, read_id, seq)
                rows.append(dict(read_id=read_id, source_id=source_id, source_class=kind,
                                 error="exact", changes=0, gene="", chromosome=chrom,
                                 strand=strand, transcript="", source_position=coordinate,
                                 junction_position=-1))
                source_id += 1
    # Independently check represented positives against the common FASTA after
    # sampling. This checks the generator; it does not call or filter on a mapper.
    all_t2g = {}
    with (REF / "simpleaf_splici_r91_k31/ref/t2g_3col.tsv").open() as f:
        for line in f:
            name, gene, _ = line.split()
            all_t2g[name] = gene
    checked_transcripts = set()
    represented_introns = set()

    def check_reference(name, chunks):
        if name is None:
            return
        seq = b"".join(chunks).upper()
        if name in cache:
            assert seq == cache[name], f"genome/GTF transcript differs from reference: {name}"
            checked_transcripts.add(name)
        gene = all_t2g[name]
        for sid, window in intron_windows.get(gene, {}).items():
            if sid not in represented_introns and (window in seq or rc(window) in seq):
                represented_introns.add(sid)

    with (REF / "simpleaf_splici_r91_k31/ref/roers_ref.fa").open("rb") as f:
        name, chunks = None, []
        for line in f:
            if line.startswith(b">"):
                check_reference(name, chunks)
                name, chunks = line[1:].split()[0].decode(), []
            else:
                chunks.append(line.strip())
        check_reference(name, chunks)
    assert checked_transcripts == set(cache)
    assert len(represented_introns) == 30_000, f"missing intronic sources: {30_000-len(represented_introns)}"
    pq.write_table(pa.Table.from_pylist(rows), DATA / "synthetic-truth.parquet", compression="zstd")
    return dict(seed=SEED, num_reads=len(rows), positive_sources=90_000,
                negative_sources=60_000, eligible_transcripts=len(mature_ids),
                eligible_introns=len(introns), reference_checked_transcripts=len(checked_transcripts),
                reference_checked_intronic_sources=len(represented_introns), source_counts=dict(collections.Counter(
                    row["source_class"] for row in rows)))


def real():
    source = ROOT / "benchdata/throughput-2026-10-27"
    rows = []
    with gzip.open(source / "pbmc10m_R1.fastq.gz", "rb") as old1, \
         gzip.open(source / "pbmc10m_R2.fastq.gz", "rb") as old2, \
         gzip.open(DATA / "real1m_R1.fastq.gz", "wb", compresslevel=1) as r1, \
         gzip.open(DATA / "real1m_R2.fastq.gz", "wb", compresslevel=1) as r2:
        for read_id in range(1_000_000):
            one, two = [old1.readline() for _ in range(4)], [old2.readline() for _ in range(4)]
            assert one[0].split()[0] == two[0].split()[0] and two[-1]
            seq1, seq2, qual2 = one[1].strip(), two[1].strip(), two[3].strip()
            assert len(seq2) == len(qual2)
            write_pair(r1, r2, read_id, seq2, qual2)
            rows.append(dict(read_id=read_id, original_name=two[0].strip().decode(),
                             original_cb=seq1[:16].decode(), original_umi=seq1[16:28].decode(),
                             original_r1_valid=len(seq1)>=28 and all(b in b"ACGT" for b in seq1[:28])))
    pq.write_table(pa.Table.from_pylist(rows), DATA / "real1m-truth.parquet", compression="zstd")
    return dict(num_reads=len(rows), original_r1_valid=sum(r["original_r1_valid"] for r in rows))


if __name__ == "__main__":
    DATA.mkdir(parents=True, exist_ok=True)
    manifest = dict(synthetic=synthetic(), real1m=real(), umi_xor_mask=MASK,
                    fixed_barcode=CB, read_length=91, quality_model="Phred40 for synthetic; original R2 for real")
    manifest["inputs"] = {}
    for path in [REF / "Homo_sapiens.GRCh38.93.gtf",
                 REF / "Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz",
                 REF / "simpleaf_splici_r91_k31/ref/roers_ref.fa",
                 REF / "simpleaf_splici_r91_k31/ref/t2g_3col.tsv",
                 *sorted(DATA.glob("*.fastq.gz"))]:
        with path.open("rb") as handle:
            digest = hashlib.file_digest(handle, "sha256").hexdigest()
        manifest["inputs"][str(path.relative_to(ROOT))] = dict(bytes=path.stat().st_size, sha256=digest)
    (RAW / "input-manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    print(json.dumps({k:v for k,v in manifest.items() if k!="inputs"}), flush=True)

"""Decode native outputs, joining by encoded UMI, never mapped-record ordinal.

RAD v2 schema checked against piscem-rs 0.11.1 src/io/rad.rs. Legacy RAD parsing
is retained for cross-checking with the installed alevin-fry view implementation.
Only the exact uncompressed single-barcode schema used here is accepted.
"""
from __future__ import annotations

import argparse
import collections
import json
from pathlib import Path
import struct

import pyarrow as pa
import pyarrow.parquet as pq

from make_inputs import CB, DATA, MASK, RAW, REF, ROOT


def pack_dna(s):
    out = 0
    for base in s:
        out = (out << 2) | "ACGT".index(base)
    return out


def unpack_dna(value, length):
    return "".join("ACGT"[(value >> shift) & 3] for shift in range(length*2-2, -1, -2))


def read_exact(handle, size):
    data = handle.read(size)
    assert len(data) == size, (len(data), size)
    return data


def read_int(handle, fmt):
    return struct.unpack("<"+fmt, read_exact(handle, struct.calcsize("<"+fmt)))[0]


def read_string(handle):
    return read_exact(handle, read_int(handle, "H")).decode()


def constellation_genes():
    with (ROOT / "benchdata/throughput-2026-10-27/reference.ptrhash.ecidx").open("rb") as f:
        header = read_exact(f,128)
        assert header[:8] == b"CSTLEC2\0" and struct.unpack_from("<I",header,8)[0] == 9
        n = struct.unpack_from("<Q",header,24)[0]
        genes_offset = struct.unpack_from("<Q",header,80)[0]
        names_offset, names_length = struct.unpack_from("<QQ",header,112)
        f.seek(genes_offset)
        records = read_exact(f,n*8)
        f.seek(names_offset)
        names = read_exact(f,names_length)
    return [names[a:a+b].decode() for a,b in struct.iter_unpack("<II",records)]


def gene_maps():
    original = constellation_genes()
    genes = sorted(original)
    assert len(genes) == len(set(genes))
    ids = {name:i for i,name in enumerate(genes)}
    t2g = {}
    with (REF / "simpleaf_splici_r91_k31/ref/t2g_3col.tsv").open() as f:
        for line in f:
            target, gene, _ = line.split()
            t2g[target] = ids[gene]
    assert set(t2g.values()) == set(range(len(genes)))
    (RAW / "gene-names.json").write_text(json.dumps(genes)+"\n")
    return genes, [ids[name] for name in original], t2g


def constellation_records(path):
    with path.open("rb") as f:
        header = read_exact(f,24)
        assert header == b"CSTRAD1\0" + struct.pack("<IIII",1,28,1,0)
        while raw := f.read(28):
            code, valid, primary, candidate_count, score, flags, state, reserved, ngenes = struct.unpack("<QIIIHHBBH",raw)
            genes = struct.unpack("<"+"I"*ngenes,read_exact(f,ngenes*4))
            bc = "".join("ACGT"[(code >> (i*2))&3] for i in range(16))
            umi = "".join("ACGT"[(code >> (i*2))&3] for i in range(16,28))
            yield dict(bc=pack_dna(bc), umi=pack_dna(umi), genes=genes, state=state,
                       primary=primary, valid=valid, score=score, flags=flags)


def rad_records(path, metadata=None):
    with path.open("rb") as f:
        version = 1
        if read_exact(f,8) == b"RAD_FILE":
            major, minor, extlen = struct.unpack("<BBI",read_exact(f,6))
            assert (major,minor) == (2,0)
            read_exact(f,extlen)
            version = 2
        else:
            f.seek(0)
        header_start=f.tell()
        paired, nrefs = struct.unpack("<BQ",read_exact(f,9))
        names = [read_string(f) for _ in range(nrefs)]
        chunk_count_offset=f.tell()
        nchunks = read_int(f,"Q")
        sections = []
        for _ in range(3):
            tags = []
            for _ in range(read_int(f,"H")):
                name, typeid = read_string(f), read_int(f,"B")
                assert typeid in (1,2,3,4,8), (name,typeid)
                role = None
                if version == 2:
                    role_code, plen = struct.unpack("<BB",read_exact(f,2))
                    role = (role_code, list(read_exact(f,plen)))
                tags.append((name,typeid,role))
            sections.append(tags)
        widths = {1:"B",2:"H",3:"I",4:"Q"}
        file_values = {name:read_string(f) if kind==8 else read_int(f,widths[kind])
                       for name,kind,_ in sections[0]}
        assert [(name,kind) for name,kind,_ in sections[1]] == [("b",3),("u",3)]
        assert len(sections[2])==1 and sections[2][0][1]==3, sections[2]
        if version==2:
            assert sections[1][0][2] == (1,[0,16])
            assert sections[1][1][2] == (2,[12])
            assert file_values["known_rad_type"] == "sc_rna_basic"
        if metadata is not None:
            metadata.update(version=version, nrefs=nrefs, nchunks=nchunks,
                            tags=sections, file_values=file_values, header_start=header_start,
                            chunk_count_offset=chunk_count_offset, chunks_offset=f.tell())
        for _ in range(nchunks):
            chunk_bytes, nrecords = struct.unpack("<II",read_exact(f,8))
            chunk = read_exact(f,chunk_bytes-8)
            position = 0
            for _ in range(nrecords):
                na, bc, umi = struct.unpack_from("<III",chunk,position)
                position += 12
                hits = struct.unpack_from("<"+"I"*na,chunk,position)
                position += na*4
                assert all((hit&0x7fffffff)<nrefs for hit in hits)
                yield bc, umi, [(names[h&0x7fffffff],bool(h&0x80000000)) for h in hits]
            assert position==len(chunk), (position,len(chunk))
        assert not f.read(1), "unparsed bytes after declared chunks"


def bus_records(directory, t2g):
    names = (directory / "transcripts.txt").read_text().splitlines()
    ecs = {}
    with (directory / "matrix.ec").open() as f:
        for line in f:
            fields = line.split()
            ecs[int(fields[0])] = tuple(sorted({t2g[names[int(i)]] for i in fields[1].split(",")}))
    with (directory / "output.bus").open("rb") as f:
        magic, version, bclen, umilen, textlen = struct.unpack("<4sIIII",read_exact(f,20))
        assert (magic,version,bclen,umilen)==(b"BUS\0",1,16,12)
        read_exact(f,textlen)
        while raw := f.read(32):
            bc, umi, ec, count, flags, padding = struct.unpack("<QQiIII",raw)
            assert count==1 and padding==0, (count,padding)
            # Validation runs also request kallisto's native input-read number.
            assert flags == umi^MASK, (flags,umi^MASK)
            yield bc,umi,ecs[ec]


def normalize(corpus, tool, nreads, original_ids, t2g):
    directory = DATA / f"{corpus}-{tool}"
    gene_sets = [()] * nreads
    evidence_sets = [()] * nreads
    state = [0] * nreads
    seen = bytearray(nreads)
    intern = {():()}
    count = 0
    meta = {}
    if tool.startswith("constellation"):
        def rows():
            for ordinal,row in enumerate(constellation_records(directory / "map.cstrad")):
                assert row["valid"]==(1<<28)-1
                assert row["umi"]^MASK == ordinal
                genes = tuple(sorted(original_ids[i] for i in row["genes"]))
                evidence = genes
                if row["state"] in (1,3):
                    assert original_ids[row["primary"]] in genes
                    # Constellation resolves some total-score ties using exon/
                    # intron support. Its primary assignment is the final call;
                    # retain the broader tied-score evidence separately.
                    genes = (original_ids[row["primary"]],)
                yield row["bc"],row["umi"],genes,row["state"],evidence
    elif tool.startswith("piscem"):
        def rows():
            for bc,umi,hits in rad_records(directory / "map.rad",meta):
                genes = tuple(sorted({t2g[target] for target,fw in hits if fw or not tool.endswith("forward")}))
                yield bc,umi,genes,1 if len(genes)==1 else (2 if genes else 0),genes
    else:
        def rows():
            for bc,umi,genes in bus_records(directory,t2g):
                yield bc,umi,genes,1 if len(genes)==1 else (2 if genes else 0),genes
    expected_bc = pack_dna(CB)
    for bc,umi,genes,label,evidence in rows():
        assert bc==expected_bc, (bc,expected_bc)
        read_id = umi^MASK
        assert 0 <= read_id < nreads and not seen[read_id], (read_id,nreads)
        assert len(genes)==len(set(genes))
        seen[read_id] = 1
        gene_sets[read_id] = intern.setdefault(genes,genes)
        evidence_sets[read_id] = intern.setdefault(evidence,evidence)
        state[read_id] = label
        count += 1
    if tool.startswith("constellation"):
        assert count==nreads
    elif tool.startswith("piscem"):
        info = json.loads((directory/"map_info.json").read_text())
        assert count==info["num_mapped"] and info["num_reads"]==nreads
    else:
        info = json.loads((directory/"run_info.json").read_text())
        assert count==info["n_pseudoaligned"] and info["n_processed"]==nreads
    table = pa.table({"read_id":pa.array(range(nreads),type=pa.uint32()),
                      "genes":pa.array(gene_sets,type=pa.list_(pa.uint32())),
                      "evidence_genes":pa.array(evidence_sets,type=pa.list_(pa.uint32())),
                      "assignment_type":pa.array(state,type=pa.uint8())})
    pq.write_table(table,DATA/f"{corpus}-{tool}.parquet",compression="zstd")
    meta.update(records=count, returned_gene_sets=sum(bool(g) for g in gene_sets),
                unique_gene=sum(len(g)==1 for g in gene_sets),
                unique_calls_from_multiple_score_ties=sum(len(g)==1 and len(e)>1 for g,e in zip(gene_sets,evidence_sets)),
                assignment_type_counts=dict(collections.Counter(state)))
    (RAW/f"{corpus}-{tool}-parse.json").write_text(json.dumps(meta,indent=2)+"\n")
    print(json.dumps(dict(corpus=corpus,tool=tool,**{k:v for k,v in meta.items() if k!="tags"})),flush=True)


if __name__=="__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--corpus",required=True,choices=["synthetic","real1m"])
    parser.add_argument("--tools",nargs="+",default=["constellation","piscem","kallisto"])
    args = parser.parse_args()
    names,original,t2g=gene_maps()
    nreads=pq.read_metadata(DATA/f"{args.corpus}-truth.parquet").num_rows
    for tool in args.tools:
        normalize(args.corpus,tool,nreads,original,t2g)

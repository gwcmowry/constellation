"""Serial 12-core cost comparison, exact output projection, and real-data counting."""
import collections
import hashlib
import heapq
import itertools
import json
import os
from pathlib import Path
import statistics
import struct
import subprocess
import time

ROOT=Path(__file__).resolve().parents[2]
RAW=Path(__file__).parent/"raw"
DATA=ROOT/"benchdata/splicing-2026-10-30"
INPUT=ROOT/"benchdata/throughput-2026-10-27"
DATA.mkdir(exist_ok=True)
BEFORE=ROOT/"benchdata/comparison-2026-10-27/constellation.before-splicing"
AFTER=ROOT/"target/release/constellation"
ENV=dict(os.environ,RAYON_NUM_THREADS="12",OMP_NUM_THREADS="12")
MODES={"before-weighted":(BEFORE,[]),"after-weighted":(AFTER,[]),
       "before-consensus":(BEFORE,["--gene-ec-consensus"]),
       "after-consensus":(AFTER,["--gene-ec-consensus"]),
       "splicing":(AFTER,["--splicing-state"])}


def digest(path):
    with path.open("rb") as f:return hashlib.file_digest(f,"sha256").hexdigest()


def run(label,repeat,cmd):
    prefix=RAW/f"{label}-{repeat}"
    command=["/usr/bin/time","-v","-o",str(prefix)+".time","taskset","-c","0-11",*map(str,cmd)]
    start=time.perf_counter()
    with Path(str(prefix)+".log").open("w") as log:
        result=subprocess.run(command,env=ENV,stdout=log,stderr=log,cwd=ROOT)
    elapsed=time.perf_counter()-start
    row=dict(label=label,repeat=repeat,seconds=elapsed,command=command,exit_code=result.returncode)
    Path(str(prefix)+".run.json").write_text(json.dumps(row,indent=2)+"\n")
    assert result.returncode==0,prefix
    print(label,repeat,round(elapsed,3),flush=True)
    return row


for i,repeat in enumerate(["warmup",1,2,3]):
    order=list(MODES)[i:]+list(MODES)[:i]
    for label in order:
        binary,extra=MODES[label]
        run(label,repeat,[binary,"map",*extra,"--index",INPUT/"reference.ptrhash.ecidx",
                         "--r1",INPUT/"pbmc10m_R1.fastq.gz","--r2",INPUT/"pbmc10m_R2.fastq.gz",
                         "--out",DATA/f"{label}.rad","--emit-metrics",RAW/f"{label}-{repeat}.metrics.json"])
summary={}
for label in MODES:
    times=[json.loads((RAW/f"{label}-{i}.run.json").read_text())["seconds"] for i in [1,2,3]]
    summary[label]=dict(seconds=times,median_seconds=statistics.median(times),
                        million_pairs_per_second=10/statistics.median(times),output_bytes=(DATA/f"{label}.rad").stat().st_size)
summary["binaries"]={str(p):digest(p) for p in [BEFORE,AFTER]}
(RAW/"summary.json").write_text(json.dumps(summary,indent=2)+"\n")

checks={}
for mode in ["weighted","consensus"]:
    a,b=[digest(DATA/f"{phase}-{mode}.rad") for phase in ["before","after"]]
    assert a==b
    checks[mode]=dict(byte_identical=True,sha256=a)
states=collections.Counter()
with (DATA/"splicing.rad").open("rb") as s,(DATA/"after-consensus.rad").open("rb") as c:
    h=s.read(28);assert h[:8]==b"CSTRAD2\0"
    ngenes=struct.unpack_from("<I",h,24)[0]
    names=[]
    for _ in range(ngenes):names.append(s.read(struct.unpack("<H",s.read(2))[0]).decode())
    assert len(names)==len(set(names))
    assert c.read(24)==b"CSTRAD1\0"+struct.pack("<IIII",1,28,1,0)
    reads=0
    while fixed:=s.read(28):
        assert fixed==c.read(28)
        n=struct.unpack_from("<H",fixed,26)[0]
        entries=s.read(n*5);assert len(entries)==n*5
        genes=[]
        for gene,state in struct.iter_unpack("<IB",entries):
            assert gene<ngenes and state in [1,2,3,4,8]
            genes.append(gene);states[state]+=1
        assert c.read(n*4)==struct.pack("<"+"I"*n,*genes)
        reads+=1
    assert not c.read(1) and reads==10_000_000
checks["splicing_projection"]=dict(identical_gene_assignments=True,reads=reads,genes=ngenes,candidate_state_counts=dict(states))
(RAW/"output-checks.json").write_text(json.dumps(checks,indent=2)+"\n")

run("count-splicing","validation",[AFTER,"count-splicing","--assignments",DATA/"splicing.rad",
       "--barcode-whitelist",ROOT/"benchdata/pbmc_10k_v3/filtered_feature_bc_matrix/barcodes.tsv.gz",
       "--out-prefix",DATA/"counts","--emit-metrics",RAW/"count.metrics.json"])


def matrix(label):
    with (DATA/f"counts_{label}.mtx").open() as f:
        assert next(f).startswith("%%MatrixMarket")
        genes,cells,n=map(int,next(f).split());assert genes==ngenes
        previous=None;seen=0
        for line in f:
            g,c,v=map(int,line.split())
            assert 1<=g<=genes and 1<=c<=cells and v>0
            assert previous is None or previous<(g,c)
            previous=(g,c);seen+=1
            yield g,c,v
        assert seen==n


def combined(labels):
    for (g,c),entries in itertools.groupby(heapq.merge(*(matrix(l) for l in labels)),key=lambda r:r[:2]):
        yield g,c,sum(v for _,_,v in entries)


for labels,total in [(["spliced","unspliced","ambiguous"],"matrix"),(["compatible_both","conflicting","unknown"],"ambiguous")]:
    for a,b in itertools.zip_longest(combined(labels),matrix(total)):assert a==b
m=json.loads((RAW/"count.metrics.json").read_text())
assert sum(v for _,_,v in matrix("matrix"))==m["molecules_after_umi_correction"]
(RAW/"matrix-checks.json").write_text(json.dumps(dict(total_equals_s_plus_u_plus_a=True,ambiguous_equals_both_plus_conflict_plus_unknown=True,total_molecules=m["molecules_after_umi_correction"]),indent=2)+"\n")
print("All output and molecule conservation checks passed",flush=True)

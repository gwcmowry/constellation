"""Frozen-policy holdouts and serial comparison. Run with the UV-managed venv."""
import contextlib
import gzip
import hashlib
import json
from pathlib import Path
import sys

ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/"experiments/2026-10-27-comparison"))
import make_inputs as m
import run_mapping as r
import read_outputs as o
import evaluate as e
import pyarrow as pa
import pyarrow.parquet as pq

RAW=Path(__file__).parent/"raw"
r.RAW=e.RAW=RAW
e.TOOLS=["constellation","constellation-consensus","piscem-forward","kallisto-forward"]


def generate():
    data=m.DATA
    for corpus,length,seed in [("heldout",91,20261030),("heldout150",150,20261031)]:
        m.DATA=data/corpus
        m.DATA.mkdir(exist_ok=True)
        m.SEED=seed
        m.READ_LENGTH=length
        with (RAW/f"{corpus}-generation.log").open("w") as log,contextlib.redirect_stdout(log):
            manifest=m.synthetic()
        manifest.update(read_length=length,seed=seed,files={})
        for suffix in ["_R1.fastq.gz","_R2.fastq.gz","-truth.parquet"]:
            source=m.DATA/("synthetic"+suffix)
            dest=data/(corpus+suffix)
            source.replace(dest)
            with dest.open("rb") as f: digest=hashlib.file_digest(f,"sha256").hexdigest()
            manifest["files"][dest.name]=dict(bytes=dest.stat().st_size,sha256=digest)
        (RAW/f"{corpus}-manifest.json").write_text(json.dumps(manifest,indent=2)+"\n")
        print(corpus,"generated and reference-validated",flush=True)
    m.DATA=data
    rows=[]
    original=ROOT/"benchdata/pbmc_10k_v3/pbmc_10k_v3_fastqs"
    with gzip.open(original/"pbmc_10k_v3_S1_L002_R1_001.fastq.gz","rb") as a, \
         gzip.open(original/"pbmc_10k_v3_S1_L002_R2_001.fastq.gz","rb") as b, \
         gzip.open(data/"lane2_R1.fastq.gz","wb",compresslevel=1) as x, \
         gzip.open(data/"lane2_R2.fastq.gz","wb",compresslevel=1) as y:
        for i in range(1_000_000):
            one,two=[a.readline() for _ in range(4)],[b.readline() for _ in range(4)]
            assert one[0].split()[0]==two[0].split()[0] and two[-1]
            seq=one[1].strip()
            m.write_pair(x,y,i,two[1].strip(),two[3].strip())
            rows.append(dict(read_id=i,original_cb=seq[:16].decode(),original_umi=seq[16:28].decode(),
                             original_r1_valid=len(seq)>=28 and all(c in b"ACGT" for c in seq[:28])))
    pq.write_table(pa.Table.from_pylist(rows),data/"lane2-truth.parquet",compression="zstd")
    print("Second lane generated",flush=True)


def validate():
    _,original,t2g=o.gene_maps()
    for corpus,n in [("heldout",420000),("heldout150",420000),("lane2",1000000)]:
        for tool in e.TOOLS:
            r.run(tool,corpus,"validation")
            with (RAW/f"{corpus}-{tool}-parse.log").open("w") as log,contextlib.redirect_stdout(log):
                o.normalize(corpus,tool,n,original,t2g)
        with (RAW/f"{corpus}-evaluation.log").open("w") as log,contextlib.redirect_stdout(log):
            if corpus=="lane2": e.real(corpus)
            else: e.synthetic(corpus)
        print(corpus,"evaluated",flush=True)


def benchmark():
    for i,repeat in enumerate(["warmup",1,2,3]):
        order=e.TOOLS[i:]+e.TOOLS[:i]
        for tool in order:r.run(tool,"pbmc10m",repeat)
    with (RAW/"timing-summary.log").open("w") as log,contextlib.redirect_stdout(log):e.timing()


if __name__=="__main__":
    generate()
    validate()
    benchmark()

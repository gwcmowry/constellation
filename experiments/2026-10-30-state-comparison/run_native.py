"""Reproduce native USA runs using the downloaded, pinned alevin-fry release."""
import gzip
import json
import os
import subprocess
import time
from pathlib import Path

ROOT=Path(__file__).resolve().parents[2]
DATA=ROOT/'benchdata/state-comparison-2026-10-30'
RAW=Path(__file__).parent/'raw'
AF=DATA/'alevin-fry-x86_64-unknown-linux-gnu/alevin-fry'
MAP=ROOT/'benchdata/comparison-2026-10-27/pbmc10m-piscem-forward'
G=DATA/'af-collation'
W=DATA/'barcodes.txt'
with gzip.open(ROOT/'benchdata/pbmc_10k_v3/filtered_feature_bc_matrix/barcodes.tsv.gz','rt') as f:
    W.write_text(''.join(line.strip().split('-')[0]+'\n' for line in f))
commands=[('generate-permit-list',['generate-permit-list','-i',MAP,'-o',G,'-d','fw','-t','12','--valid-bc',W,'--cell-bc-correction','unique','--cell-bc-neighborhood','hamming-1']),
          ('collate',['collate','-i',G,'-r',MAP,'-t','12'])]
for label,destination,extra in [('quant','af-usa',[]),('prefer-ambiguity','af-usa-prefer-ambiguity',['--sa-model','prefer-ambig'])]:
    commands.append((label,['quant','-i',G,'-o',DATA/destination,'-m',ROOT/'benchdata/reference/ensembl93/simpleaf_splici_r91_k31/ref/t2g_3col.tsv',
                            '-r','cr-like','-t','12','--small-thresh','0','--dump-molecules',*extra]))
for label,args in commands:
    cmd=['/usr/bin/time','-v','-o',str(RAW/f'af-{label}.time'),'taskset','-c','0-11',str(AF),*map(str,args)]
    start=time.perf_counter()
    with (RAW/f'af-{label}.log').open('w') as log:
        result=subprocess.run(cmd,env=dict(os.environ,RAYON_NUM_THREADS='12'),stdout=log,stderr=log,cwd=ROOT)
    info=dict(command=cmd,seconds=time.perf_counter()-start,exit_code=result.returncode)
    (RAW/f'af-{label}.run.json').write_text(json.dumps(info,indent=2)+'\n')
    print(label,info['seconds'],result.returncode,flush=True)
    assert result.returncode==0

"""Sparse matrix agreement on canonical axes, after all timed native jobs."""
import json
from pathlib import Path
import numpy as np
import pyarrow.parquet as pq
ROOT=Path(__file__).resolve().parents[2];RAW=Path(__file__).parent/'raw';DATA=ROOT/'benchdata/large-scale-2026-11-01/matrices'
GENES=json.loads((ROOT/'experiments/2026-11-01-competitor-stages/raw/canonical-genes.json').read_text());CELLS=json.loads((ROOT/'experiments/2026-11-01-competitor-stages/raw/canonical-cells.json').read_text());NC=len(CELLS);N=len(GENES)*NC

def load(name,state):
    t=pq.read_table(DATA/f'{name}-{state}.parquet');g=t['gene_id'].combine_chunks().to_numpy().astype(np.uint64);c=t['cell_id'].combine_chunks().to_numpy().astype(np.uint64)
    values=t['count'].combine_chunks().to_numpy();keys=g*NC+c
    assert len(keys)==0 or (keys[1:]>keys[:-1]).all()
    return keys,values

def compare(a,b):
    ka,x=a;kb,y=b;sx,sy=int(x.sum()),int(y.sum());xx,yy=int(np.dot(x,x)),int(np.dot(y,y))
    common,ia,ib=np.intersect1d(ka,kb,assume_unique=True,return_indices=True)
    xy=int(np.dot(x[ia],y[ib]));minimum=int(np.minimum(x[ia],y[ib]).sum());l1=sx+sy-2*minimum
    gx=np.bincount((ka//NC).astype(np.int64),weights=x.astype(np.float64),minlength=len(GENES));gy=np.bincount((kb//NC).astype(np.int64),weights=y.astype(np.float64),minlength=len(GENES))
    return {'a_total':sx,'b_total':sy,'nonzero_union':len(ka)+len(kb)-len(common),'l1':l1,'l1_over_sum':l1/(sx+sy),'weighted_jaccard':minimum/(sx+sy-minimum),'full_matrix_pearson':(xy-sx*sy/N)/np.sqrt((xx-sx*sx/N)*(yy-sy*sy/N)),'gene_total_log1p_pearson':float(np.corrcoef(np.log1p(gx),np.log1p(gy))[0,1])}
report={}
for size in [100,500]:
    ref=f'{size}m-piscem-count'
    for a in [f'{size}m-ours-count',f'{size}m-ours-usa',f'{size}m-kallisto-count']:
        label=a+'_vs_'+ref;report[label]={}
        for state in 'SUA':
            report[label][state]=compare(load(a,state),load(ref,state))
            print(label,state,report[label][state]['full_matrix_pearson'],flush=True)
        (RAW/'matrix-comparison.json').write_text(json.dumps(report,indent=2)+'\n')

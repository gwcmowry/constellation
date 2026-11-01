"""Check the native default small-cell fast path against the controlled USA run."""
import json,sys,time
from run import run,AF,OLD,REF,STATE,DATA,RAW,ROOT
import numpy as np
sys.path.insert(0,str(ROOT/'experiments/2026-10-30-state-comparison'))
from compare_counts import native_af
expected=native_af()
results=[]
for i in range(1,4):
    name=f'af-default-{i}';base=DATA/name;base.mkdir(exist_ok=True);col=base/'collation';q=base/'quant'
    start=time.perf_counter()
    run(name+'-permit',[AF,'generate-permit-list','-i',OLD/'pbmc10m-piscem-forward','-o',col,'-d','fw','-t','12','--valid-bc',STATE/'barcodes.txt','--cell-bc-correction','unique','--cell-bc-neighborhood','hamming-1'])
    run(name+'-collate',[AF,'collate','-i',col,'-r',OLD/'pbmc10m-piscem-forward','-t','12'])
    run(name+'-quant',[AF,'quant','-i',col,'-o',q,'-m',REF/'t2g_3col.tsv','-r','cr-like','-t','12'])
    elapsed=time.perf_counter()-start
    actual=native_af(q/'alevin')
    for state in 'SUA':assert all(np.array_equal(a,b) for a,b in zip(expected[state],actual[state]))
    metadata=json.loads((q/'quant.json').read_text())
    results.append({'seconds':elapsed,'small_thresh':metadata['quant_options']['small_thresh'],'tiny_cell_resolved_count':len(metadata['tiny_cell_resolved_cell_numbers']),'matrix_matches_control':True})
(RAW/'af-default-control.json').write_text(json.dumps(results,indent=2)+'\n')

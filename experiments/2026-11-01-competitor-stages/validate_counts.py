"""Normalize native axes; validate repeated outputs and exact matched-evidence control."""
import hashlib,json,sys
from pathlib import Path
import numpy as np
import pyarrow as pa
import pyarrow.parquet as pq
ROOT=Path(__file__).resolve().parents[2]
RAW=Path(__file__).parent/'raw';DATA=ROOT/'benchdata/competitor-stages-2026-11-01'
sys.path.insert(0,str(ROOT/'experiments/2026-10-30-state-comparison'))
from compare_counts import native_ours,native_af,read_mtx,sorted_arrays,IDS,CB_IDS,NCELLS,GENES,BARCODES

def native_bus(folder):
    genes=(folder/'counts.genes.txt').read_text().splitlines()
    cells=(folder/'counts.barcodes.txt').read_text().splitlines()
    gids=np.array([IDS[g] for g in genes],dtype=np.uint64);cids=np.array([CB_IDS[c] for c in cells],dtype=np.uint64)
    output={}
    for state,suffix in [('S','.mtx'),('U','.2.mtx'),('A','.ambiguous.mtx')]:
        keys=[];values=[]
        for c,g,value in read_mtx(folder/('counts'+suffix)):
            if value>0:keys.append(gids[g]*NCELLS+cids[c]);values.append(value)
            else:assert value==0
        output[state]=sorted_arrays(keys,values)
    return output

def digest(states):
    h=hashlib.sha256()
    for s,(keys,counts) in states.items():
        assert np.isfinite(counts).all() and (counts>0).all() and (counts==np.round(counts)).all()
        h.update(s.encode());h.update(keys.tobytes());h.update(counts.tobytes())
    return h.hexdigest()

report={'totals':{},'matrix_digests':{},'repeat_identity':{},'native_metadata':{}}
gold=native_af()
for name,loader in [('ours-default',lambda i:native_ours(DATA/f'ours-ours-default-{i}')),
    ('ours-usa',lambda i:native_ours(DATA/f'ours-ours-usa-{i}')),
    ('ours-piscem-usa',lambda i:native_ours(DATA/f'ours-piscem-usa-{i}')),
    ('alevin-fry',lambda i:native_af(DATA/f'af-{i}/quant/alevin')),
    ('bustools',lambda i:native_bus(DATA/f'bus-{i}'))]:
    hashes=[]
    for i in range(1,4):
        states=loader(i);hashes.append(digest(states))
        if name in ['ours-piscem-usa','alevin-fry']:
            for s in 'SUA':assert all(np.array_equal(a,b) for a,b in zip(states[s],gold[s]))
        if name=='alevin-fry':
            info=json.loads((DATA/f'af-{i}/quant/quant.json').read_text())['quant_options']
            assert info['sa_model']=='WinnerTakeAll' and info['crlike_umi_edit']==0 and not info['dump_molecules'] and info['small_thresh']==0
            report['native_metadata'][str(i)]=info
        if i==1:
            report['totals'][name]={s:int(count.sum()) for s,(_,count) in states.items()}
            report['totals'][name]['total']=sum(report['totals'][name].values())
            for s,(keys,count) in states.items():
                pq.write_table(pa.table({'gene_id':pa.array(keys//NCELLS,type=pa.uint32()),'cell_id':pa.array(keys%NCELLS,type=pa.uint32()),'count':pa.array(count,type=pa.uint64())}),DATA/f'{name}-{s}.parquet',compression='zstd')
    assert len(set(hashes))==1,name
    report['matrix_digests'][name]=hashes[0];report['repeat_identity'][name]=True
report['matched_piscem_evidence_usa_matrices_exactly_match_native']=True
report['all_counts_positive_finite_integers']=True
(RAW/'canonical-genes.json').write_text(json.dumps(GENES)+'\n');(RAW/'canonical-cells.json').write_text(json.dumps(BARCODES)+'\n')
(RAW/'validation.json').write_text(json.dumps(report,indent=2)+'\n')
print(json.dumps(report['totals'],indent=2))

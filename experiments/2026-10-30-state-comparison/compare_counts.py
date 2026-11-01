"""Project piscem target/state compatibility losslessly into the common counter.

The adapter is counting-only: record order is piscem's output order, not original
FASTQ order. Score/flags are unknown and set to zero; the counter does not use them.
"""
import collections
import hashlib
import json
import os
from pathlib import Path
import struct
import subprocess
import time

import numpy as np
import pyarrow as pa
import pyarrow.parquet as pq

from compare_reads import ROOT,RAW,DATA,OLD,GENES,IDS,T2G,constellation_genes,rad_records,unpack_dna

BARCODES=sorted((DATA/'barcodes.txt').read_text().splitlines())
CB_IDS={cb:i for i,cb in enumerate(BARCODES)}
NGENES,NCELLS=len(GENES),len(BARCODES)


def adapter():
    names=constellation_genes();original={g:i for i,g in enumerate(names)}
    targets={t:(original[GENES[g]],s) for t,(g,s) in T2G.items()}
    reverse=[sum(((b>>(2*i))&3)<<(2*(3-i)) for i in range(4)) for b in range(256)]
    def little(value,nbytes):
        out=0
        for j in range(nbytes):out=(out<<8)|reverse[(value>>(8*j))&255]
        return out
    destination=DATA/'piscem-common-counter.cstrad2'
    total=mapped=0
    with destination.open('wb',buffering=1<<20) as f:
        f.write(b'CSTRAD2\0'+struct.pack('<5I',2,28,1,5,len(names)))
        for name in names:b=name.encode();f.write(struct.pack('<H',len(b))+b)
        for bc,umi,hits in rad_records(OLD/'pbmc10m-piscem-forward/map.rad'):
            by_gene={}
            for target,fw in hits:
                if fw:
                    g,s=targets[target];by_gene[g]=by_gene.get(g,0)|s
            pairs=sorted(by_gene.items());n=len(pairs)
            code=little(bc,4)|(little(umi,3)<<32)
            if total<10_000:
                decoded=''.join('ACGT'[(code>>(2*i))&3] for i in range(28))
                assert decoded==unpack_dna(bc,16)+unpack_dna(umi,12)
            f.write(struct.pack('<QIIIHHBBH',code,(1<<28)-1,pairs[0][0] if n==1 else 0xffffffff,n,0,0,1 if n==1 else 2 if n else 5,0,n))
            for g,s in pairs:f.write(struct.pack('<IB',g,s))
            total+=1;mapped+=bool(n)
    info=json.loads((OLD/'pbmc10m-piscem-forward/map_info.json').read_text());assert total==info['num_mapped']
    assert mapped==7_754_302 # Independently reported by native GPL with -d fw.
    (RAW/'adapter.json').write_text(json.dumps(dict(input_rad_records=total,forward_mapped=mapped,barcode_umi_roundtrips_checked=10000,output=str(destination),score_and_flags='not supplied; zero; ignored by counter'),indent=2)+'\n')
    print('piscem converted',total,mapped,flush=True)


def count(label,path,correct):
    cmd=['/usr/bin/time','-v','-o',str(RAW/f'{label}.time'),'taskset','-c','0-11',str(ROOT/'target/release/constellation'),
         'count-splicing','--assignments',str(path),'--barcode-whitelist',str(DATA/'barcodes.txt'),
         '--out-prefix',str(DATA/label),'--emit-metrics',str(RAW/f'{label}.metrics.json')]
    if not correct:cmd+=['--no-correct-umis']
    start=time.perf_counter()
    with (RAW/f'{label}.log').open('w') as log:subprocess.run(cmd,env=dict(os.environ,RAYON_NUM_THREADS='12'),stdout=log,stderr=log,check=True)
    (RAW/f'{label}.run.json').write_text(json.dumps(dict(command=cmd,seconds=time.perf_counter()-start),indent=2)+'\n')
    print(label,'counted',flush=True)


def read_mtx(path):
    with path.open() as f:
        line=next(f);assert line.startswith('%%MatrixMarket')
        while True:
            line=next(f)
            if not line.startswith('%'):break
        nr,nc,n=map(int,line.split())
        for line in f:
            a,b,v=line.split();yield int(a)-1,int(b)-1,float(v)


def native_ours(prefix):
    names=[line.split('\t')[1] for line in Path(str(prefix)+'_features.tsv').read_text().splitlines()]
    cells=Path(str(prefix)+'_barcodes.tsv').read_text().splitlines()
    gene_ids=np.array([IDS[g] for g in names],dtype=np.uint64)
    cell_ids=np.array([CB_IDS[c] for c in cells],dtype=np.uint64)
    output={}
    for state,label in [('S','spliced'),('U','unspliced'),('A','ambiguous')]:
        keys=[];values=[]
        for g,c,v in read_mtx(Path(str(prefix)+f'_{label}.mtx')):
            keys.append(gene_ids[g]*NCELLS+cell_ids[c]);values.append(v)
        output[state]=sorted_arrays(keys,values)
    return output


def sorted_arrays(keys,values):
    keys=np.asarray(keys,dtype=np.uint64);values=np.asarray(values,dtype=np.float64)
    order=np.argsort(keys);keys,values=keys[order],values[order]
    assert len(keys)==0 or (keys[1:]>keys[:-1]).all()
    return keys,values


def native_af(folder=None):
    folder=folder or DATA/'af-usa/alevin'
    names=(folder/'quants_mat_cols.txt').read_text().splitlines()
    cells=(folder/'quants_mat_rows.txt').read_text().splitlines()
    assert len(names)==3*NGENES
    base=names[:NGENES]
    assert names[NGENES:2*NGENES]==[g+'-U' for g in base]
    assert names[2*NGENES:]==[g+'-A' for g in base]
    gids=np.array([IDS[g] for g in base],dtype=np.uint64)
    cids=np.array([CB_IDS[c] for c in cells],dtype=np.uint64)
    rows={s:([],[]) for s in 'SUA'}
    for c,feature,v in read_mtx(folder/'quants_mat.mtx'):
        state='SUA'[feature//NGENES];g=feature%NGENES
        rows[state][0].append(gids[g]*NCELLS+cids[c]);rows[state][1].append(v)
    return {s:sorted_arrays(*v) for s,v in rows.items()}


def compare(a,b):
    ka,va=a;kb,vb=b
    keys=np.union1d(ka,kb);x=np.zeros(len(keys));y=x.copy()
    x[np.searchsorted(keys,ka)]=va;y[np.searchsorted(keys,kb)]=vb
    sx,sy=x.sum(),y.sum();xx,yy=np.dot(x,x),np.dot(y,y);xy=np.dot(x,y);n=NGENES*NCELLS
    l1=np.abs(x-y).sum()
    gx=np.bincount((keys//NCELLS).astype(np.int64),weights=x,minlength=NGENES)
    gy=np.bincount((keys//NCELLS).astype(np.int64),weights=y,minlength=NGENES)
    return dict(a_total=float(sx),b_total=float(sy),l1=float(l1),l1_over_sum=float(l1/(sx+sy)),
                weighted_jaccard=float((sx+sy-l1)/(sx+sy+l1)),
                full_matrix_pearson=float((xy-sx*sy/n)/np.sqrt((xx-sx*sx/n)*(yy-sy*sy/n))),
                gene_total_log1p_pearson=float(np.corrcoef(np.log1p(gx),np.log1p(gy))[0,1]),
                nonzero_union=len(keys))


def summarize():
    (RAW/'cell-barcodes.json').write_text(json.dumps(BARCODES)+'\n')
    datasets={'ours_corrected':native_ours(ROOT/'benchdata/splicing-2026-10-30/counts'),
              'ours_exact_umi':native_ours(DATA/'ours-exact-umi'),
              'piscem_common_corrected':native_ours(DATA/'piscem-common-corrected'),
              'piscem_common_exact_umi':native_ours(DATA/'piscem-common-exact-umi'),
              'piscem_native_af':native_af()}
    report={'universe':dict(genes=NGENES,barcodes=NCELLS),'totals':{},'comparisons':{}}
    for name,states in datasets.items():
        report['totals'][name]={s:float(v.sum()) for s,(k,v) in states.items()}
        report['totals'][name]['total']=sum(report['totals'][name].values())
        for s,(keys,values) in states.items():
            pq.write_table(pa.table({'gene_id':pa.array(keys//NCELLS,type=pa.uint32()),'cell_id':pa.array(keys%NCELLS,type=pa.uint32()),'count':values}),DATA/f'{name}-{s}.parquet',compression='zstd')
    for a,b in [('ours_corrected','piscem_common_corrected'),('ours_exact_umi','piscem_common_exact_umi'),
                ('ours_corrected','piscem_native_af'),('ours_exact_umi','piscem_native_af'),
                ('piscem_common_exact_umi','piscem_native_af')]:
        report['comparisons'][a+'_vs_'+b]={s:compare(datasets[a][s],datasets[b][s]) for s in 'SUA'}
    # Native molecule-table rows include unresolved groups; only assigned rows
    # contribute to the integer matrix. Validate that interpretation independently.
    native=pq.read_table(DATA/'af-usa/alevin/molecules.parquet',columns=['status'])
    native=pa.table({'status':native['status'].cast(pa.string())})
    status=native.group_by('status').aggregate([('status','count')]).to_pylist()
    report['native_molecule_status_counts']=status
    assigned=sum(r['status_count'] for r in status if r['status']=='assigned')
    assert assigned==report['totals']['piscem_native_af']['total']
    report['native_assigned_rows_match_matrix']=True
    (RAW/'count-comparison.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(report['totals'],indent=2),flush=True)


def analyze_alternative():
    prefer=native_af(DATA/'af-usa-prefer-ambiguity/alevin')
    report={'prefer_ambiguity_totals':{s:float(v.sum()) for s,(k,v) in prefer.items()},'comparisons':{}}
    for name in ['ours_exact_umi','piscem_common_exact_umi','piscem_native_af']:
        result={}
        for state in 'SUA':
            t=pq.read_table(DATA/f'{name}-{state}.parquet').to_pydict()
            keys=np.asarray(t['gene_id'],dtype=np.uint64)*NCELLS+np.asarray(t['cell_id'],dtype=np.uint64)
            values=np.asarray(t['count'])
            result[state]=compare((keys,values),prefer[state])
        report['comparisons'][name+'_vs_af_prefer_ambiguity']=result
    (RAW/'native-ambiguity-model-comparison.json').write_text(json.dumps(report,indent=2)+'\n')


if __name__=='__main__':
    adapter()
    count('piscem-common-corrected',DATA/'piscem-common-counter.cstrad2',True)
    count('piscem-common-exact-umi',DATA/'piscem-common-counter.cstrad2',False)
    count('ours-exact-umi',ROOT/'benchdata/splicing-2026-10-30/splicing.rad',False)
    for label in ['piscem-common-corrected','piscem-common-exact-umi']:
        assert json.loads((RAW/f'{label}.metrics.json').read_text())['mapped_valid_rows']==7_204_991
    summarize()
    analyze_alternative()

"""Compare native read-level gene/state compatibility; identities are UMI-tagged."""
import collections
import contextlib
import itertools
import json
import os
from pathlib import Path
import struct
import subprocess
import sys
import time

import pyarrow as pa
import pyarrow.parquet as pq

ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/"experiments/2026-10-27-comparison"))
from read_outputs import rad_records,read_exact,read_string,pack_dna,unpack_dna,constellation_genes
from make_inputs import MASK,CB

RAW=Path(__file__).parent/"raw"
OLD=ROOT/"benchdata/comparison-2026-10-27"
DATA=ROOT/"benchdata/state-comparison-2026-10-30"
REFERENCE=ROOT/"benchdata/reference/ensembl93/simpleaf_splici_r91_k31/ref"
GENES=sorted(constellation_genes())
IDS={g:i for i,g in enumerate(GENES)}
T2G={}
for line in (REFERENCE/"t2g_3col.tsv").read_text().splitlines():
    target,gene,state=line.split();assert state in ["S","U"]
    T2G[target]=(IDS[gene],1 if state=="S" else 2)
LABELS={0:"missing",1:"S",2:"U",3:"A",4:"C",8:"X"}


def cstrad2(path):
    with path.open("rb") as f:
        h=read_exact(f,28);assert h[:8]==b"CSTRAD2\0"
        assert struct.unpack_from("<4I",h,8)==(2,28,1,5)
        names=[read_string(f) for _ in range(struct.unpack_from("<I",h,24)[0])]
        assert len(names)==len(set(names))
        while fixed:=f.read(28):
            assert len(fixed)==28
            code,valid,primary,candidates,score,flags,kind,reserved,n=struct.unpack("<QIIIHHBBH",fixed)
            entries=struct.unpack("<"+"IB"*n,read_exact(f,n*5))
            pairs=tuple(sorted((IDS[names[entries[j]]],entries[j+1]) for j in range(0,2*n,2)))
            assert all(s in [1,2,3,4,8] for _,s in pairs)
            umi=0
            for i in range(16,28):umi=(umi<<2)|((code>>(2*i))&3)
            yield umi^MASK,pairs


def collapse(targets):
    by_gene={}
    for target in targets:
        g,s=T2G[target];by_gene[g]=by_gene.get(g,0)|s
    return tuple(sorted(by_gene.items()))


def load_piscem(corpus,n):
    rows=[()]*n;seen=set();metadata={}
    folder=OLD/f"{corpus}-piscem-forward"
    for bc,umi,hits in rad_records(folder/"map.rad",metadata):
        i=umi^MASK;assert 0<=i<n and i not in seen and bc==pack_dna(CB)
        seen.add(i);rows[i]=collapse(target for target,fw in hits if fw)
    info=json.loads((folder/"map_info.json").read_text())
    assert len(seen)==info["num_mapped"] and info["num_reads"]==n
    return rows


def load_kallisto(corpus,n):
    folder=OLD/f"{corpus}-kallisto-forward"
    targets=(folder/"transcripts.txt").read_text().splitlines()
    ecs={}
    for line in (folder/"matrix.ec").read_text().splitlines():
        ec,values=line.split();ecs[int(ec)]=collapse(targets[int(i)] for i in values.split(','))
    rows=[()]*n;seen=set()
    with (folder/"output.bus").open('rb') as f:
        magic,v,b,u,length=struct.unpack('<4s4I',read_exact(f,20));assert (magic,v,b,u)==(b'BUS\0',1,16,12)
        read_exact(f,length)
        while raw:=f.read(32):
            bc,umi,ec,count,flags,pad=struct.unpack('<QQiIII',raw)
            i=umi^MASK;assert i==flags and 0<=i<n and i not in seen and count==1 and bc==pack_dna(CB)
            seen.add(i);rows[i]=ecs[ec]
    info=json.loads((folder/'run_info.json').read_text());assert len(seen)==info['n_pseudoaligned'] and info['n_processed']==n
    return rows


def read_comparison(a,b,eligible):
    c=collections.Counter();conf=collections.Counter();shared=collections.Counter()
    for x,y,valid in zip(a,b,eligible):
        if not valid:continue
        c['reads']+=1;c['exact_gene_state_sets']+=x==y
        c['a_mapped']+=bool(x);c['b_mapped']+=bool(y)
        if x and y:
            c['both_mapped']+=1;c['exact_sets_both_mapped']+=x==y
        dx,dy=dict(x),dict(y)
        if dx.keys()==dy.keys() and dx:c['same_nonempty_gene_set']+=1
        for g in dx.keys()&dy.keys():shared[(LABELS[dx[g]],LABELS[dy[g]])]+=1
        if len(x)==len(y)==1 and x[0][0]==y[0][0]:
            c['same_unique_gene']+=1
            conf[(LABELS[x[0][1]],LABELS[y[0][1]])]+=1
            c['same_state_same_unique_gene']+=x[0][1]==y[0][1]
    c=dict(c)
    c['state_agreement_given_same_unique_gene']=c.get('same_state_same_unique_gene',0)/max(1,c.get('same_unique_gene',0))
    c['confusion_same_unique_gene']={f'{a}->{b}':v for (a,b),v in sorted(conf.items())}
    c['confusion_all_shared_genes']={f'{a}->{b}':v for (a,b),v in sorted(shared.items())}
    return c


def save(corpus,tool,rows):
    pq.write_table(pa.table({'read_id':pa.array(range(len(rows)),type=pa.uint32()),
                            'genes':pa.array([[g for g,s in r] for r in rows],type=pa.list_(pa.uint32())),
                            'states':pa.array([[s for g,s in r] for r in rows],type=pa.list_(pa.uint8()))}),
                   DATA/f'{corpus}-{tool}.parquet',compression='zstd')


def main():
    (RAW/'gene-names.json').write_text(json.dumps(GENES)+'\n')
    results={}
    for corpus,n in [('real1m',1000000),('lane2',1000000),('heldout',420000),('heldout150',420000)]:
        out=DATA/f'{corpus}-ours.cstrad2'
        cmd=['taskset','-c','0-11',str(ROOT/'target/release/constellation'),'map','--splicing-state',
             '--index',str(ROOT/'benchdata/throughput-2026-10-27/reference.ptrhash.ecidx'),
             '--r1',str(OLD/f'{corpus}_R1.fastq.gz'),'--r2',str(OLD/f'{corpus}_R2.fastq.gz'),'--out',str(out)]
        with (RAW/f'{corpus}-map.log').open('w') as log:subprocess.run(cmd,env=dict(os.environ,RAYON_NUM_THREADS='12'),stdout=log,stderr=log,check=True)
        (RAW/f'{corpus}-map-command.json').write_text(json.dumps(cmd,indent=2)+'\n')
        ours=[]
        for i,row in cstrad2(out):assert i==len(ours);ours.append(row)
        assert len(ours)==n
        rows={'ours':ours,'piscem':load_piscem(corpus,n),'kallisto':load_kallisto(corpus,n)}
        for tool,value in rows.items():save(corpus,tool,value)
        truth=pq.read_table(OLD/f'{corpus}-truth.parquet').to_pydict()
        eligible=truth.get('original_r1_valid',[True]*n)
        report={'input_reads':n,'eligible_reads':sum(eligible),'tools':{},'comparisons':{}}
        for tool,data in rows.items():
            counter=collections.Counter()
            for r,valid in zip(data,eligible):
                if not valid:continue
                if not r:counter['unmapped']+=1
                elif len(r)>1:counter['multi_gene']+=1
                else:counter[LABELS[r[0][1]]]+=1
            report['tools'][tool]=dict(counter)
        for a,b in itertools.combinations(rows,2):report['comparisons'][f'{a}_vs_{b}']=read_comparison(rows[a],rows[b],eligible)
        results[corpus]=report
        (RAW/f'{corpus}-read-comparison.json').write_text(json.dumps(report,indent=2)+'\n')
        print(corpus,{k:v['state_agreement_given_same_unique_gene'] for k,v in report['comparisons'].items()},flush=True)
    (RAW/'read-summary.json').write_text(json.dumps(results,indent=2)+'\n')


if __name__=='__main__':main()

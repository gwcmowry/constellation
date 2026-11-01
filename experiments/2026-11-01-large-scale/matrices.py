"""Normalize native sparse matrices with Arrow's compiled parser; keep Parquet."""
import hashlib,json,sys
from pathlib import Path
import numpy as np
import pyarrow as pa
import pyarrow.csv as csv
import pyarrow.parquet as pq
ROOT=Path(__file__).resolve().parents[2];RAW=Path(__file__).parent/'raw';DATA=ROOT/'benchdata/large-scale-2026-11-01'
PREV=ROOT/'experiments/2026-11-01-competitor-stages/raw'
GENES=json.loads((PREV/'canonical-genes.json').read_text());CELLS=json.loads((PREV/'canonical-cells.json').read_text())
GIDS={g:i for i,g in enumerate(GENES)};CIDS={c:i for i,c in enumerate(CELLS)};NC=len(CELLS)

def mtx(path):
    skip=0
    with path.open() as f:
        while True:
            line=f.readline();skip+=1
            if not line.startswith('%'):shape=tuple(map(int,line.split()));break
    if shape[2]==0:return tuple(np.array([],dtype=np.uint64) for _ in range(3)),shape
    table=csv.read_csv(path,read_options=csv.ReadOptions(column_names=['a','b','n'],skip_rows=skip,block_size=8<<20,use_threads=True),parse_options=csv.ParseOptions(delimiter=' '),convert_options=csv.ConvertOptions(column_types={'a':pa.uint32(),'b':pa.uint32(),'n':pa.float64()}))
    assert table.num_rows==shape[2],(path,shape,table.num_rows)
    counts=table['n'].combine_chunks().to_numpy()
    assert np.isfinite(counts).all() and (counts>=0).all() and (counts==np.round(counts)).all()
    return (table['a'].combine_chunks().to_numpy().astype(np.uint64)-1,table['b'].combine_chunks().to_numpy().astype(np.uint64)-1,counts.astype(np.uint64)),shape

def write_state(label,state,g,c,n):
    keep=n>0;g,c,n=g[keep],c[keep],n[keep]
    keys=g*NC+c;order=np.argsort(keys);keys,n=keys[order],n[order]
    assert len(keys)==0 or (keys[1:]>keys[:-1]).all(),(label,state,'duplicate matrix coordinates')
    table=pa.table({'gene_id':pa.array(keys//NC,type=pa.uint32()),'cell_id':pa.array(keys%NC,type=pa.uint32()),'count':pa.array(n,type=pa.uint64())})
    out=DATA/'matrices';out.mkdir(exist_ok=True)
    pq.write_table(table,out/f'{label}-{state}.parquet',compression='zstd',row_group_size=1_000_000)
    h=hashlib.sha256();h.update(memoryview(keys));h.update(memoryview(n));digest=h.hexdigest()
    return {'total':int(n.sum()),'nonzero':len(n),'sha256_canonical_arrays':digest,'genes_with_counts':int(len(np.unique(keys//NC))),'cells_with_counts':int(len(np.unique(keys%NC)))}

def normalize(label,kind,path):
    report={}
    if kind=='ours':
        genes=[line.split('\t')[1] for line in Path(str(path)+'_features.tsv').read_text().splitlines()];cells=Path(str(path)+'_barcodes.tsv').read_text().splitlines()
        gids=np.array([GIDS[g] for g in genes],dtype=np.uint64);cids=np.array([CIDS[c] for c in cells],dtype=np.uint64)
        for state,suffix in [('S','spliced'),('U','unspliced'),('A','ambiguous')]:
            (a,b,n),shape=mtx(Path(str(path)+f'_{suffix}.mtx'));assert shape[:2]==(len(genes),len(cells))
            report[state]=write_state(label,state,gids[a],cids[b],n)
        metrics=json.loads((RAW/f'{label}.metrics.json').read_text());assert sum(x['total'] for x in report.values())==metrics['molecules_after_umi_correction']
    elif kind=='af':
        path=path/'alevin';names=(path/'quants_mat_cols.txt').read_text().splitlines();cells=(path/'quants_mat_rows.txt').read_text().splitlines();ng=len(names)//3
        genes=names[:ng];assert names[ng:2*ng]==[g+'-U' for g in genes] and names[2*ng:]==[g+'-A' for g in genes]
        gids=np.array([GIDS[g] for g in genes],dtype=np.uint64);cids=np.array([CIDS[c] for c in cells],dtype=np.uint64)
        (a,b,n),shape=mtx(path/'quants_mat.mtx');assert shape[:2]==(len(cells),len(names))
        for k,state in enumerate('SUA'):
            mask=b//ng==k;report[state]=write_state(label,state,gids[b[mask]%ng],cids[a[mask]],n[mask])
    else:
        genes=Path(str(path)+'.genes.txt').read_text().splitlines();cells=Path(str(path)+'.barcodes.txt').read_text().splitlines()
        gids=np.array([GIDS[g] for g in genes],dtype=np.uint64);cids=np.array([CIDS[c] for c in cells],dtype=np.uint64)
        for state,suffix in [('S','.mtx'),('U','.2.mtx'),('A','.ambiguous.mtx')]:
            (a,b,n),shape=mtx(Path(str(path)+suffix));assert shape[:2]==(len(cells),len(genes));report[state]=write_state(label,state,gids[b],cids[a],n)
    report['total']=sum(x['total'] for x in report.values());report['barcodes']=len(cells);report['genes']=len(genes)
    (RAW/f'{label}.matrix-summary.json').write_text(json.dumps(report,indent=2)+'\n')
    print(label,report['total'],flush=True)
    return report
if __name__=='__main__':normalize(sys.argv[1],sys.argv[2],Path(sys.argv[3]))

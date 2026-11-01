"""Isolate the USA spliced tie preference on single-read synthetic molecules.

This intentionally measures the cost of resolving sequence ambiguity; source
origin is not always identifiable from shared sequence. No parameters are tuned.
"""
import collections,json
from pathlib import Path
import pyarrow.parquet as pq
ROOT=Path(__file__).resolve().parents[2]
GENES=json.loads((ROOT/'experiments/2026-10-30-state-comparison/raw/gene-names.json').read_text())
report={}
for corpus in ['heldout','heldout150']:
    truth=pq.read_table(ROOT/f'benchdata/comparison-2026-10-27/{corpus}-truth.parquet').to_pylist()
    rows=pq.read_table(ROOT/f'benchdata/state-comparison-2026-10-30/{corpus}-ours.parquet').to_pylist()
    counts=collections.Counter(); by_class=collections.defaultdict(collections.Counter)
    for row in rows:
        tr=truth[row['read_id']]; positive=tr['source_class'] in ['mature','junction','intronic']
        pairs=list(zip(row['genes'],row['states']))
        old=pairs[0][0] if len(pairs)==1 else None
        call=old
        if len(pairs)>1 and all(s in [1,2,3] for _,s in pairs):
            labels=[(g,v) for g,s in pairs for v in (1,2) if s&v]
            spliced=[g for g,s in labels if s==1]
            if 2<=len(labels)<=10 and len(spliced)==1:call=spliced[0]
        if positive:
            counts['positives']+=1
            for policy,gene in [('intersection',old),('usa',call)]:
                if gene is not None:
                    counts[policy+'_assigned']+=1
                    counts[policy+'_correct']+=GENES[gene]==tr['gene']
            if call is not None and old is None:
                tag='correct' if GENES[call]==tr['gene'] else 'wrong_source_gene'
                counts['additional_'+tag]+=1;by_class[tr['source_class']][tag]+=1
        else:
            counts[tr['source_class']+'_n']+=1
            counts[tr['source_class']+'_intersection_assigned']+=old is not None
            counts[tr['source_class']+'_usa_assigned']+=call is not None
    for policy in ['intersection','usa']:
        counts[policy+'_precision']=counts[policy+'_correct']/counts[policy+'_assigned']
        counts[policy+'_recall']=counts[policy+'_correct']/counts['positives']
    report[corpus]={'counts':dict(counts),'additional_by_origin':{k:dict(v) for k,v in by_class.items()}}
Path(__file__).with_name('raw').joinpath('quality.json').write_text(json.dumps(report,indent=2)+'\n')
print(json.dumps(report,indent=2))

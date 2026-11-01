"""Independent exact-read class audit by full-sequence searches in the FASTA.

The balanced random audit estimates neither biological maturity nor full-genome
origin. It checks compatibility with the supplied reference for the source gene.
The enriched discordance sample is diagnostic and is summarized separately.
"""
import collections
import gzip
import json
import random
from pathlib import Path

import pyarrow as pa
import pyarrow.parquet as pq

from compare_reads import RAW,DATA,OLD,REFERENCE,GENES,IDS,T2G,LABELS


def main():
    rng=random.Random(20261030)
    queries=[]; predictions={};truths={}
    for corpus in ['heldout','heldout150']:
        truth=pq.read_table(OLD/f'{corpus}-truth.parquet').to_pydict();truths[corpus]=truth
        rows={}
        for tool in ['ours','piscem','kallisto']:
            t=pq.read_table(DATA/f'{corpus}-{tool}.parquet').to_pydict()
            rows[tool]=[dict(zip(g,s)) for g,s in zip(t['genes'],t['states'])]
        predictions[corpus]=rows
        selected={}
        for source in ['mature','junction','intronic']:
            eligible=[i for i,(s,e) in enumerate(zip(truth['source_class'],truth['error'])) if s==source and e=='exact']
            for i in rng.sample(eligible,500):selected[i]='balanced'
        discordant=[]
        for i,(gene,error) in enumerate(zip(truth['gene'],truth['error'])):
            if not gene or error!='exact':continue
            g=IDS[gene];a,b=rows['ours'][i],rows['piscem'][i]
            if a.get(g)==3 and b.get(g) in (1,2) and len(a)==len(b)==1:discordant.append(i)
        for i in rng.sample(discordant,min(200,len(discordant))):
            selected[i]='balanced+discordance' if i in selected else 'discordance'
        with gzip.open(OLD/f'{corpus}_R2.fastq.gz','rb') as f:
            i=0
            while header:=f.readline():
                seq=f.readline().strip();assert f.readline().startswith(b'+');assert f.readline()
                if i in selected:
                    queries.append(dict(corpus=corpus,read_id=i,gene=IDS[truth['gene'][i]],sequence=seq,
                                        sample=selected[i],source_class=truth['source_class'][i],expected=0))
                i+=1
        assert i==420000
    by_gene=collections.defaultdict(list)
    for q in queries:by_gene[q['gene']].append(q)
    def consume(target,pieces):
        if target is None:return
        gene,state=T2G[target]
        if gene not in by_gene:return
        seq=b''.join(pieces).upper()
        for q in by_gene[gene]:
            if not q['expected']&state and q['sequence'] in seq:q['expected']|=state
    with (REFERENCE/'roers_ref.fa').open('rb') as f:
        target=None;pieces=[]
        for line in f:
            if line.startswith(b'>'):
                consume(target,pieces);target=line[1:].split()[0].decode();pieces=[]
            else:pieces.append(line.strip())
        consume(target,pieces)
    report={};cases=[]
    for corpus in ['heldout','heldout150']:
        report[corpus]={}
        for sample in ['balanced','discordance']:
            subset=[q for q in queries if q['corpus']==corpus and sample in q['sample']]
            out={'n':len(subset),'reference_masks':dict(collections.Counter(LABELS[q['expected']] for q in subset)),'tools':{}}
            for tool in predictions[corpus]:
                stats=collections.Counter()
                for q in subset:
                    expected=q['expected'];actual=predictions[corpus][tool][q['read_id']].get(q['gene'],0)
                    if not expected:label='source_not_present_forward_in_reference'
                    elif not actual:label='source_gene_not_reported'
                    elif actual==expected:label='exact_state_set'
                    elif actual==3:label='extra_state_ambiguity'
                    elif actual==4:label='conflicting_evidence'
                    elif actual==8:label='unknown'
                    elif expected==3:label='lost_reference_ambiguity'
                    else:label='opposite_definite_state'
                    stats[label]+=1
                out['tools'][tool]=dict(stats)
            report[corpus][sample]=out
    # Diagnose up to 200 exact discordances per read length: did sparse sampling
    # miss a class-discriminating k-mer, or is the gene/state k-mer union itself
    # too broad even if EVERY k-mer is checked?
    targets_needed=collections.defaultdict(list)
    for q in queries:
        if 'discordance' in q['sample'] and q['expected'] in (1,2):
            q['extra_state']=3^q['expected'];q['uncovered_kmers']=set(q['sequence'][i:i+31] for i in range(len(q['sequence'])-30))
            targets_needed[q['gene']].append(q)
    def examine(target,pieces):
        if target is None:return
        g,s=T2G[target]
        if g not in targets_needed:return
        seq=b''.join(pieces).upper()
        for q in targets_needed[g]:
            if s==q['extra_state']:
                q['uncovered_kmers'].difference_update([k for k in q['uncovered_kmers'] if k in seq])
    with (REFERENCE/'roers_ref.fa').open('rb') as f:
        target=None;pieces=[]
        for line in f:
            if line.startswith(b'>'):
                examine(target,pieces);target=line[1:].split()[0].decode();pieces=[]
            else:pieces.append(line.strip())
        examine(target,pieces)
    diagnosis=collections.defaultdict(collections.Counter)
    for q in queries:
        d={k:v for k,v in q.items() if k not in ['uncovered_kmers','sequence']}
        d.setdefault('extra_state',None)
        d['diagnosis']=None;d['extra_class_absent_kmers']=None
        d['gene_name']=GENES[q['gene']]
        d['sequence']=q['sequence'].decode()
        for tool in predictions[q['corpus']]:d[tool]=predictions[q['corpus']][tool][q['read_id']].get(q['gene'],0)
        if 'uncovered_kmers' in q:
            reason='discriminating_kmers_not_enforced' if q['uncovered_kmers'] else 'target_path_identity_lost_even_with_dense_kmers'
            diagnosis[q['corpus']][reason]+=1;d['diagnosis']=reason;d['extra_class_absent_kmers']=len(q['uncovered_kmers'])
        cases.append(d)
    with (DATA/'discriminating-kmers.tsv').open('w') as f:
        for q in queries:
            if 'uncovered_kmers' in q:
                for kmer in sorted(q['uncovered_kmers']):
                    f.write(f"{q['corpus']}:{q['read_id']}\t{GENES[q['gene']]}\t{q['extra_state']}\t{kmer.decode()}\n")
    report['discordance_diagnosis']={k:dict(v) for k,v in diagnosis.items()}
    (RAW/'reference-audit.json').write_text(json.dumps(report,indent=2)+'\n')
    pq.write_table(pa.Table.from_pylist(cases),DATA/'reference-audit-cases.parquet',compression='zstd')
    print(json.dumps(report,indent=2),flush=True)


if __name__=='__main__':main()

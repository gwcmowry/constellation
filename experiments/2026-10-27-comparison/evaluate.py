"""Truth accuracy, real-read agreement, common exact-UMI count sanity check."""
import collections
import itertools
import json
from pathlib import Path
import statistics

import numpy as np
import pyarrow.parquet as pq

from make_inputs import DATA, RAW

TOOLS = ["constellation", "piscem", "kallisto"]


def read_sets(corpus, tool):
    table = pq.read_table(DATA/f"{corpus}-{tool}.parquet")
    assert table["read_id"].to_pylist() == list(range(len(table)))
    return [tuple(g) for g in table["genes"].to_pylist()], table["assignment_type"].to_numpy()


def ratio(a,b):
    return a/b if b else None


def synthetic(corpus="synthetic"):
    truth = pq.read_table(DATA/f"{corpus}-truth.parquet").to_pydict()
    names = json.loads((RAW/"gene-names.json").read_text())
    ids = {g:i for i,g in enumerate(names)}
    truth_genes = np.array([ids.get(g,-1) for g in truth["gene"]])
    groups = collections.defaultdict(list)
    for i,(source,error) in enumerate(zip(truth["source_class"],truth["error"])):
        groups[f"{source}/{error}"].append(i)
        groups["positive" if truth_genes[i]>=0 else "negative"].append(i)
        if truth_genes[i]>=0:
            groups[f"positive/{error}"].append(i)
    report = {}
    for tool in TOOLS:
        sets, state = read_sets(corpus,tool)
        assert len(sets)==len(truth_genes)
        sizes=np.fromiter((len(s) for s in sets),dtype=np.int32)
        correct=np.fromiter((len(s)==1 and s[0]==g for s,g in zip(sets,truth_genes)),dtype=bool)
        retained=np.fromiter((g in s for s,g in zip(sets,truth_genes)),dtype=bool)
        evidence=pq.read_table(DATA/f"{corpus}-{tool}.parquet",columns=["evidence_genes"])["evidence_genes"].to_pylist()
        evidence_retained=np.fromiter((g in s for s,g in zip(evidence,truth_genes)),dtype=bool)
        tool_report={}
        for name, rows in groups.items():
            ix=np.array(rows)
            n=len(ix)
            unique=int((sizes[ix]==1).sum())
            right=int(correct[ix].sum())
            mapped=int((sizes[ix]>0).sum())
            pos=bool((truth_genes[ix]>=0).all())
            tool_report[name]=dict(n=n,unique=unique,ambiguous=int((sizes[ix]>1).sum()),
                                   unmapped=int((sizes[ix]==0).sum()),mapped=mapped,
                                   assignment_rate=mapped/n,unique_assignment_rate=unique/n)
            if pos:
                tool_report[name].update(correct_unique=right,wrong_unique=unique-right,
                                         unique_precision=ratio(right,unique),unique_recall=right/n,
                                         wrong_unique_fraction=(unique-right)/n,
                                         candidate_set_recall=float(retained[ix].mean()),
                                         evidence_set_recall=float(evidence_retained[ix].mean()))
            if tool=="constellation":
                countable=np.isin(state[ix],[1,3])
                tool_report[name]["native_countable"]=int(countable.sum())
                tool_report[name]["native_countable_correct"]=int((countable&correct[ix]).sum())
        report[tool]=tool_report
    (RAW/("accuracy.json" if corpus=="synthetic" else f"{corpus}-accuracy.json")).write_text(json.dumps(report,indent=2)+"\n")
    for tool in TOOLS:
        print(tool,json.dumps({k:report[tool][k] for k in ["positive","intergenic/exact","random/exact"]}),flush=True)


def sparse_correlation(a,b,log=False):
    keys=set(a)|set(b)
    x=np.fromiter((a.get(k,0) for k in keys),dtype=float)
    y=np.fromiter((b.get(k,0) for k in keys),dtype=float)
    if log:
        x,y=np.log1p(x),np.log1p(y)
    return dict(pearson=float(np.corrcoef(x,y)[0,1]),nonzero_union=len(keys))


def real(corpus="real1m"):
    truth=pq.read_table(DATA/f"{corpus}-truth.parquet").to_pydict()
    sets={}; states={}
    for tool in TOOLS:
        sets[tool],states[tool]=read_sets(corpus,tool)
    n=len(truth["read_id"])
    valid=truth["original_r1_valid"]
    summary={tool:dict(n=n,unique=sum(len(g)==1 for g in sets[tool]),
                       ambiguous=sum(len(g)>1 for g in sets[tool]),
                       unmapped=sum(not g for g in sets[tool]),
                       original_r1_valid_n=sum(valid),
                       valid_unique=sum(len(g)==1 and v for g,v in zip(sets[tool],valid)),
                       valid_ambiguous=sum(len(g)>1 and v for g,v in zip(sets[tool],valid)),
                       valid_unmapped=sum(not g and v for g,v in zip(sets[tool],valid))) for tool in TOOLS}
    comparisons={}
    for a,b in itertools.combinations(TOOLS,2):
        c=collections.Counter()
        for ga,gb,v in zip(sets[a],sets[b],valid):
            c["all_reads"]+=1
            c["equal_gene_sets"]+=ga==gb
            if not v:
                continue
            c["original_r1_valid"]+=1
            c["valid_equal_gene_sets"]+=ga==gb
            c["either_mapped"]+=bool(ga or gb)
            c["both_mapped"]+=bool(ga and gb)
            c["equal_nonempty_sets"]+=bool(ga and ga==gb)
            if ga and gb:
                c["mapped_disjoint"]+=not set(ga).intersection(gb)
            if len(ga)==len(gb)==1:
                c["both_unique"]+=1
                c["same_unique_gene"]+=ga==gb
            c["a_only_mapped"]+=bool(ga and not gb)
            c["b_only_mapped"]+=bool(gb and not ga)
            c["a_unique_b_ambiguous_contains_a"]+=bool(len(ga)==1 and len(gb)>1 and ga[0] in gb)
            c["b_unique_a_ambiguous_contains_b"]+=bool(len(gb)==1 and len(ga)>1 and gb[0] in ga)
        comparisons[f"{a}_vs_{b}"]=dict(c,unique_gene_agreement=ratio(c["same_unique_gene"],c["both_unique"]),
                                      disjoint_among_both_mapped=ratio(c["mapped_disjoint"],c["both_mapped"]))

    # Fixed, mapper-independent barcode inclusion; this is NOT a cell caller.
    cb_frequency=collections.Counter(cb for cb,v in zip(truth["original_cb"],valid) if v)
    common_barcodes={cb for cb,count in cb_frequency.items() if count>=20}
    molecules={}; matrices={}; gene_counts={}
    for tool in TOOLS:
        molecules[tool]={(cb,umi,g[0]) for g,cb,umi,v,state in zip(
            sets[tool],truth["original_cb"],truth["original_umi"],valid,states[tool])
            if len(g)==1 and v and cb in common_barcodes and (tool!="constellation" or state in (1,3))}
        matrices[tool]=collections.Counter((cb,gene) for cb,umi,gene in molecules[tool])
        gene_counts[tool]=collections.Counter(gene for cb,umi,gene in molecules[tool])
    counts=dict(policy="Original valid CB/UMI; CB >=20 input reads; unique genes; exact UMI deduplication; no correction or multi-gene resolution",
                retained_barcodes=len(common_barcodes),
                molecule_counts={t:len(m) for t,m in molecules.items()},comparisons={})
    for a,b in itertools.combinations(TOOLS,2):
        keys=set(matrices[a])|set(matrices[b])
        counts["comparisons"][f"{a}_vs_{b}"]=dict(
            shared_molecules=len(molecules[a]&molecules[b]),
            molecule_jaccard=len(molecules[a]&molecules[b])/len(molecules[a]|molecules[b]),
            matrix_log1p=sparse_correlation(matrices[a],matrices[b],log=True),
            gene_total_log1p=sparse_correlation(gene_counts[a],gene_counts[b],log=True),
            matrix_l1_over_total= sum(abs(matrices[a].get(k,0)-matrices[b].get(k,0)) for k in keys)
                                  /(len(molecules[a])+len(molecules[b])))
    report=dict(assignment_summary=summary,read_concordance=comparisons,common_count_policy=counts)
    (RAW/("real-concordance.json" if corpus=="real1m" else f"{corpus}-concordance.json")).write_text(json.dumps(report,indent=2)+"\n")
    print(json.dumps(report),flush=True)


def timing():
    report={}
    for corpus in ("pbmc10m","pbmc10m-raw"):
        for tool in TOOLS:
            results=[]
            for run in (1,2,3):
                prefix=RAW/f"{corpus}-{tool}-{run}"
                if not Path(str(prefix)+".run.json").exists():
                    continue
                row=json.loads(Path(str(prefix)+".run.json").read_text())
                assert row["exit_code"]==0
                resource={}
                for line in Path(str(prefix)+".time").read_text().splitlines():
                    for key in ["Maximum resident set size (kbytes)","User time (seconds)","System time (seconds)","Major (requiring I/O) page faults"]:
                        if line.strip().startswith(key+":"):
                            resource[key]=float(line.split(":",1)[1].strip())
                results.append(dict(wall=row["process_wall_seconds"],rss_kib=resource["Maximum resident set size (kbytes)"],
                                    cpu_seconds=resource["User time (seconds)"]+resource["System time (seconds)"],
                                    major_faults=resource["Major (requiring I/O) page faults"]))
            if results:
                wall=statistics.median(r["wall"] for r in results)
                report[f"{corpus}/{tool}"]=dict(repeats=results,median_seconds=wall,
                                                million_pairs_per_second=10/wall,
                                                median_rss_gib=statistics.median(r["rss_kib"] for r in results)/(1024**2),
                                                median_cpu_seconds=statistics.median(r["cpu_seconds"] for r in results))
    (RAW/"timing-summary.json").write_text(json.dumps(report,indent=2)+"\n")
    print(json.dumps(report),flush=True)


if __name__=="__main__":
    synthetic()
    real()
    timing()

"""Attribute assignment differences without changing or tuning any mapper."""
import collections
import json

import numpy as np
import pyarrow as pa
import pyarrow.parquet as pq

from make_inputs import DATA, RAW


def load(corpus,tool):
    return pq.read_table(DATA/f"{corpus}-{tool}.parquet").to_pydict()


def main():
    names=json.loads((RAW/"gene-names.json").read_text())
    ids={name:i for i,name in enumerate(names)}
    truth=pq.read_table(DATA/"synthetic-truth.parquet").to_pydict()
    c,p,k=[load("synthetic",tool) for tool in ["constellation","piscem","kallisto"]]
    attribution=collections.Counter()
    conservative=collections.Counter()
    records=[]
    for i,(gene,source,error) in enumerate(zip(truth["gene"],truth["source_class"],truth["error"])):
        expected=ids.get(gene,-1)
        final,evidence=c["genes"][i],c["evidence_genes"][i]
        if expected>=0:
            conservative["positive_reads"]+=1
            if len(evidence)==1:
                conservative["unique"]+=1
                conservative["correct_unique"]+=evidence[0]==expected
            conservative["source_retained_in_score_ties"]+=expected in evidence
            if len(final)==1 and final[0]!=expected:
                attribution["wrong_unique"]+=1
                attribution["wrong_unique_source_in_score_ties"]+=expected in evidence
                attribution["wrong_unique_source_absent_from_score_ties"]+=expected not in evidence
                if error=="exact":
                    attribution["wrong_unique_exact"]+=1
                    attribution["wrong_unique_exact_alternative_in_both_competitor_sets"]+=final[0] in p["genes"][i] and final[0] in k["genes"][i]
                    attribution["wrong_unique_exact_source_unique_in_both_competitors"]+=p["genes"][i]==k["genes"][i]==[expected]
        elif source=="intergenic":
            conservative["intergenic_reads"]+=1
            conservative["intergenic_unique"]+=len(evidence)==1
        if len(final)==1 and final[0]!=expected:
            records.append(dict(read_id=i,source_class=source,error=error,source_gene=gene,
                                chosen_gene=names[final[0]],
                                score_tie_genes=[names[j] for j in evidence],
                                piscem_genes=[names[j] for j in p["genes"][i]],
                                kallisto_genes=[names[j] for j in k["genes"][i]]))
    conservative=dict(conservative)
    conservative["unique_precision"]=conservative["correct_unique"]/conservative["unique"]
    conservative["unique_recall"]=conservative["correct_unique"]/conservative["positive_reads"]
    conservative["description"]="Offline diagnostic only: treat every multiple-total-score tie as ambiguous, ignoring exon/intron tie resolution. Not a new timed mapper or tuned production change."
    pq.write_table(pa.Table.from_pylist(records),DATA/"source-misassignment-cases.parquet",compression="zstd")

    # Error variants from a source whose EXACT read is uniquely and correctly
    # called by both comparison tools: a diagnostic consensus stratum, not a new
    # independent ground-truth definition or a replacement for the full result.
    consensus_sources={truth["source_id"][i] for i,gene in enumerate(truth["gene"])
                       if gene and truth["error"][i]=="exact" and p["genes"][i]==k["genes"][i]==[ids[gene]]}
    strata={}
    for label,rows in [("constellation",c),("piscem",p),("kallisto",k)]:
        counts=collections.Counter()
        for i,gene in enumerate(truth["gene"]):
            if gene and truth["source_id"][i] in consensus_sources:
                counts["positive_reads"]+=1
                genes=rows["genes"][i]
                if len(genes)==1:
                    counts["unique"]+=1
                    counts["correct_unique"]+=genes[0]==ids[gene]
        strata[label]=dict(counts,unique_precision=counts["correct_unique"]/counts["unique"],
                           unique_recall=counts["correct_unique"]/counts["positive_reads"])

    c,p,k=[load("real1m",tool) for tool in ["constellation","piscem","kallisto"]]
    valid=pq.read_table(DATA/"real1m-truth.parquet",columns=["original_r1_valid"])["original_r1_valid"].to_pylist()
    real=collections.Counter()
    for i,v in enumerate(valid):
        if not v:
            continue
        cg,pg,kg=c["genes"][i],p["genes"][i],k["genes"][i]
        if len(cg)==len(pg)==len(kg)==1 and pg==kg and cg!=pg:
            real["constellation_unique_disagrees_with_both_same_unique"]+=1
            real["consensus_gene_present_in_constellation_score_ties"]+=pg[0] in c["evidence_genes"][i]
    report=dict(synthetic_error_attribution=dict(attribution),conservative_tie_diagnostic=conservative,
                consensus_unique_exact_source_stratum=strata,real_disagreement_attribution=dict(real))
    (RAW/"quality-diagnostics.json").write_text(json.dumps(report,indent=2)+"\n")
    print(json.dumps(report,indent=2),flush=True)


if __name__=="__main__":
    main()

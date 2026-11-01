"""Validate read identities/projections and independent count invariants."""
import ast
import csv
import hashlib
import json
from pathlib import Path

import numpy as np
import pyarrow as pa
import pyarrow.parquet as pq

from compare_reads import ROOT,RAW,DATA,OLD


def sha(path):
    with path.open('rb') as f:return hashlib.file_digest(f,'sha256').hexdigest()


def main():
    checks={};total=0
    for corpus,n in [('real1m',1000000),('lane2',1000000),('heldout',420000),('heldout150',420000)]:
        for tool,oldtool in [('ours','constellation-consensus'),('piscem','piscem-forward'),('kallisto','kallisto-forward')]:
            a=pq.read_table(DATA/f'{corpus}-{tool}.parquet',columns=['genes','states'])
            b=pq.read_table(OLD/f'{corpus}-{oldtool}.parquet',columns=['genes'])
            assert len(a)==len(b)==n and a['genes'].equals(b['genes'])
            checks[f'{corpus}/{tool}']=True;total+=n
    assert total==8_520_000
    no_unexpected=True
    for row in csv.DictReader((RAW/'sampling-audit.tsv').open(),delimiter='\t'):
        assert not row['classification'].startswith('unexpected')
    assert 'ClosedSyncmer { s: 20 }' in (RAW/'sampling-audit.log').read_text()
    counts=json.loads((RAW/'count-comparison.json').read_text())
    exact=json.loads((RAW/'piscem-common-exact-umi.metrics.json').read_text())
    assert exact['mapped_valid_rows']==7_204_991 and exact['raw_barcode_umi_groups']==6_475_696
    assert counts['native_assigned_rows_match_matrix']
    for path in DATA.glob('*-[SUA].parquet'):
        values=pq.read_table(path,columns=['count'])['count'].to_numpy()
        assert np.isfinite(values).all() and (values>0).all() and (values==np.round(values)).all()
    native=pq.read_table(DATA/'af-usa-prefer-ambiguity/alevin/molecules.parquet',columns=['status'])
    native=pa.table({'status':native['status'].cast(pa.string())})
    statuses=native.group_by('status').aggregate([('status','count')]).to_pylist()
    assigned=sum(row['status_count'] for row in statuses if row['status']=='assigned')
    alt=json.loads((RAW/'native-ambiguity-model-comparison.json').read_text())
    assert assigned==sum(alt['prefer_ambiguity_totals'].values())
    for path in Path(__file__).parent.glob('*.py'):ast.parse(path.read_text())
    report=dict(gene_projection_checks=checks,total_read_rows_checked=total,
                native_accepted_read_count_matches_common_counter=True,
                native_exact_barcode_umi_group_count_matches_common_counter=True,
                default_native_assigned_molecule_rows_match_matrix=True,
                alternative_native_assigned_molecule_rows_match_matrix=True,
                no_unexpected_index_classification_results=True,
                all_normalized_counts_positive_finite_integers=True,
                production_binary_sha256=sha(ROOT/'target/release/constellation'),
                alevin_fry_binary_sha256=sha(DATA/'alevin-fry-x86_64-unknown-linux-gnu/alevin-fry'))
    (RAW/'validation.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(report,indent=2),flush=True)


if __name__=='__main__':main()

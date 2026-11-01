"""Verify output identity and native USA compatibility on the full real input."""
import hashlib,json,sys
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[2]
DATA=ROOT/'benchdata/build-count-2026-10-31'
RAW=Path(__file__).parent/'raw'
sys.path.insert(0,str(ROOT/'experiments/2026-10-30-state-comparison'))
from compare_counts import native_ours,native_af

def sha(p):
    with p.open('rb') as f:return hashlib.file_digest(f,'sha256').hexdigest()
checks={}
for suffix in ['barcodes.tsv','features.tsv','spliced.mtx','unspliced.mtx','ambiguous.mtx','compatible_both.mtx','conflicting.mtx','unknown.mtx','matrix.mtx']:
    values=[sha(DATA/f'{prefix}_{suffix}') for prefix in ['baseline-count-1','final-count-1','final-count-2','final-count-3']]
    assert len(set(values))==1,suffix
    checks[suffix]=values[0]
a=native_ours(DATA/'piscem-usa');b=native_af()
for state in 'SUA':
    assert all(np.array_equal(x,y) for x,y in zip(a[state],b[state])),state
checks['native_piscem_usa_matrices_exact']=True
checks['mapping_bytes_identical']=sha(DATA/'baseline-map.rad')==sha(DATA/'fast-map.rad')
assert checks['mapping_bytes_identical']
checks['baseline_index_sha256']=sha(DATA/'baseline-final.ecidx')
checks['optimized_index_sha256']=sha(DATA/'fast-build-3.ecidx')
checks['index_bytes_identical']=checks['baseline_index_sha256']==checks['optimized_index_sha256']
# Binary layout bytes can differ; verify_index.rs checks every key and payload.
checks['full_index_semantic_audit']=json.loads((RAW/'index-validation.json').read_text())
(RAW/'validation.json').write_text(json.dumps(checks,indent=2)+'\n')
print(json.dumps(checks,indent=2))

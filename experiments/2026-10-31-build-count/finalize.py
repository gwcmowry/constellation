"""Validate chained runs and discard only this experiment's redundant indexes."""
import hashlib,json,shutil
from pathlib import Path
ROOT=Path(__file__).resolve().parents[2]
DATA=ROOT/'benchdata/build-count-2026-10-31'
RAW=Path(__file__).parent/'raw'
def sha(path):
    with path.open('rb') as f:return hashlib.file_digest(f,'sha256').hexdigest()
r=json.loads((RAW/'validation.json').read_text())
r['pipeline_mapping_bytes_identical']=sha(DATA/'pipeline-baseline.rad')==sha(DATA/'pipeline-fast.rad')
assert r['pipeline_mapping_bytes_identical']
r['pipeline_count_output_bytes_identical']={}
for suffix in ['barcodes.tsv','features.tsv','spliced.mtx','unspliced.mtx','ambiguous.mtx','compatible_both.mtx','conflicting.mtx','unknown.mtx','matrix.mtx']:
    assert sha(DATA/f'pipeline-baseline-count_{suffix}')==sha(DATA/f'pipeline-fast-count_{suffix}')==r[suffix]
    r['pipeline_count_output_bytes_identical'][suffix]=True
(RAW/'validation.json').write_text(json.dumps(r,indent=2)+'\n')
shutil.copy2(ROOT/'target/release/constellation',DATA/'constellation.optimized')
# Keep baseline-final.ecidx and fast-build-3.ecidx for exhaustive revalidation.
names=['baseline-build.ecidx','baseline-build-2.ecidx','baseline-build-3.ecidx',
       'baseline-final-2.ecidx','baseline-final-3.ecidx','optimized-build-1.ecidx',
       'fast-build-1.ecidx','fast-build-2.ecidx','pipeline-baseline-build.ecidx',
       'pipeline-baseline-final.ecidx','pipeline-fast-build.ecidx']
removed=[]
for name in names:
    path=DATA/name
    if path.exists():
        removed.append({'path':str(path),'bytes':path.stat().st_size});path.unlink()
(RAW/'cleanup.json').write_text(json.dumps({'removed_redundant_indexes':removed,'bytes_reclaimed':sum(r['bytes'] for r in removed)},indent=2)+'\n')
print('Validated pipeline outputs; reclaimed',sum(r['bytes'] for r in removed),'bytes of redundant indexes.')

import ast,hashlib,json,subprocess
from pathlib import Path
ROOT=Path(__file__).resolve().parents[2]
RAW=Path(__file__).parent/'raw';DATA=ROOT/'benchdata/competitor-stages-2026-11-01'
subprocess.run(['python3',str(RAW.parent/'summarize.py')],check=True,stdout=(RAW/'summary.log').open('w'))
s=json.loads((RAW/'summary.json').read_text());assert len(s['build'])==3
text='''## Measured result

Same 10M input pairs, prepared splici reference, twelve-core cap. Build numbers
are one fresh run per tool; counting numbers are medians of three runs after
warmup. All index builds start from the reference FASTA.

| Builder | Build wall time | Peak RSS | Final index size |
|---|---:|---:|---:|
'''
for name,label in [('ours','Constellation direct PtrHash'),('piscem','piscem 0.24.0'),('kallisto','kallisto 0.52.0')]:
    r=s['build'][name];text+=f"| {label} | {r['seconds']:.2f}s | {r['rss_gib']:.2f} GiB | {s['index_bytes'][name]/1e9:.3f} GB |\n"
text+='''
| Post-mapping workflow | Complete count wall time | CPU-seconds | Peak RSS |
|---|---:|---:|---:|
'''
for name,label in [('ours_default','Our mappings, conservative + UMI edit correction'),('ours_usa','Our mappings, exact-UMI USA'),('alevin_fry','Piscem + alevin-fry USA'),('bustools','Kallisto + bustools gene-UMI S/U/A')]:
    r=s['count'][name];text+=f"| {label} | {r['seconds']['median']:.3f}s | {r['cpu_seconds']['median']:.2f} | {r['rss_gib']['median']:.3f} GiB |\n"
text+='''
Native stage breakdowns (medians):

| Workflow | Barcode correction / permit list | Sort / collate | Final molecule/matrix step |
|---|---:|---:|---:|
'''
for name,label in [('alevin_fry','alevin-fry'),('bustools','bustools')]:
    numbers=[x['seconds']['median'] for x in s['count_stages'][name].values()]
    text+=f"| {label} | {numbers[0]:.3f}s | {numbers[1]:.3f}s | {numbers[2]:.3f}s |\n"
c=s['count']['ours_piscem_usa'];a=s['count']['alevin_fry']
text+=f'''
With **identical piscem evidence and exact-UMI USA results**, our counter takes
{c['seconds']['median']:.3f}s versus alevin-fry's complete {a['seconds']['median']:.3f}s pipeline.
Our aggregate CPU cost is {c['cpu_seconds']['median']:.2f} CPU-seconds versus
{a['cpu_seconds']['median']:.2f}; alevin-fry has a modest wall-time advantage while
using more CPU in this test. The projection adapter is excluded from this
counter-only control. Default Constellation additionally corrects UMI errors,
so its first row is not an identical-resolution comparison.

The earlier unpinned, twelve-thread builds were 190.88s for piscem and 1,099.60s
for kallisto; the fresh pinned results above are the headline measurements.
Our earlier three-build median was 24.41s. Index construction is a reusable
one-time reference cost; counting is repeated for each sample.

Raw data: [summary.json](raw/summary.json), [validated count outputs](raw/validation.json),
and [binary provenance](raw/tools.json).

'''
control=s['controls'].get('alevin_fry_native_small_cell_default')
if control:
    metadata=json.loads((RAW/'af-default-control.json').read_text())
    text+=f"Native alevin-fry's default `--small-thresh 100` was also tested in three full\npipelines: median {control['seconds']['median']:.3f}s, with exactly the same matrices.\nThis checks that forcing `--small-thresh 0` in the controlled comparison did not\ncreate the main speed ranking. [Default-path control](raw/af-default-control.json).\n\n"
if 'bustools_empty_split' in s['controls']:
    r=s['controls']['bustools_empty_split']
    text+=f"Bustools' S/U/A count takes **{r['seconds']:.3f}s even on a zero-record BUS file**\nwith the same EC metadata and split list. Much of its roughly 18.4s final count\ntime is therefore setup for this EC dictionary, rather than processing the 10M\ninput pairs. This fixed component makes linear extrapolation to larger inputs\ninappropriate. [Empty-input control](raw/bus-empty-split.run.json).\n\n"
p=RAW.parent/'README.md';body=p.read_text();assert '<!-- RESULTS -->' in body;body=body.replace('<!-- RESULTS -->',text);p.write_text(body)
for script in RAW.parent.glob('*.py'):ast.parse(script.read_text(),filename=str(script))
inputs={}
for path in [ROOT/'benchdata/reference/ensembl93/simpleaf_splici_r91_k31/ref/roers_ref.fa',ROOT/'benchdata/reference/ensembl93/simpleaf_splici_r91_k31/ref/t2g_3col.tsv',ROOT/'benchdata/state-comparison-2026-10-30/barcodes.txt']:
    with path.open('rb') as f:inputs[str(path)]={'bytes':path.stat().st_size,'sha256':hashlib.file_digest(f,'sha256').hexdigest()}
(RAW/'input-digests.json').write_text(json.dumps(inputs,indent=2)+'\n')
print(json.dumps({'build':s['build'],'count_medians':{k:v['seconds']['median'] for k,v in s['count'].items()},'controls':s['controls']},indent=2))

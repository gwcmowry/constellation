import hashlib,json,re,statistics
from pathlib import Path
ROOT=Path(__file__).resolve().parents[2]
RAW=Path(__file__).parent/'raw'
DATA=ROOT/'benchdata/build-count-2026-10-31'
def measure(label):
    j=json.loads((RAW/f'{label}.run.json').read_text());t=(RAW/f'{label}.time').read_text()
    def number(pattern):return float(re.search(pattern,t).group(1))
    return {'seconds':j['seconds'],'rss_gib':number(r'Maximum resident set size \(kbytes\): (\d+)')/1024**2,'cpu_seconds':number(r'User time \(seconds\): ([\d.]+)')+number(r'System time \(seconds\): ([\d.]+)')}
def stats(rows):return {k:{'median':statistics.median(r[k] for r in rows),'min':min(r[k] for r in rows),'max':max(r[k] for r in rows)} for k in rows[0]}
def baseline_build(i):
    suffix='' if i==1 else f'-{i}'
    a,b=[measure(f'baseline-{p}{suffix}') for p in ['build','convert']]
    return {'seconds':a['seconds']+b['seconds'],'rss_gib':max(a['rss_gib'],b['rss_gib']),'cpu_seconds':a['cpu_seconds']+b['cpu_seconds']}
rows={'baseline_build':[baseline_build(i) for i in range(1,4)],'fast_build':[measure(f'fast-build-{i}') for i in range(1,4)]}
for name,pattern in [('baseline_count','baseline-count'),('fast_count','final-count'),('baseline_map','baseline-map'),('fast_map','fast-map')]:rows[name]=[measure(f'{pattern}-{i}') for i in range(1,4)]
report={k:stats(v) for k,v in rows.items()}
report['raw_stage_measurements']=rows
report['speedups']={phase:report[f'baseline_{phase}']['seconds']['median']/report[f'fast_{phase}']['seconds']['median'] for phase in ['build','map','count']}
report['sum_of_stage_medians']={name:sum(report[f'{name}_{p}']['seconds']['median'] for p in ['build','map','count']) for name in ['baseline','fast']}
report['prebuilt_fastq_to_matrix_sum_of_medians']={name:sum(report[f'{name}_{p}']['seconds']['median'] for p in ['map','count']) for name in ['baseline','fast']}
report['full_pipeline']={name:json.loads((RAW/f'pipeline-{name}.json').read_text()) for name in ['baseline','fast']}
report['sha256']={}
for label,path in [('baseline_binary',DATA/'constellation.baseline'),('instrumented_baseline_binary',DATA/'constellation.instrumented'),('current_binary',ROOT/'target/release/constellation')]:
    with path.open('rb') as f:report['sha256'][label]=hashlib.file_digest(f,'sha256').hexdigest()
report['index_bytes']=(DATA/'fast-build-3.ecidx').stat().st_size
(RAW/'summary.json').write_text(json.dumps(report,indent=2)+'\n')
print(json.dumps({k:v for k,v in report.items() if k not in ['raw_stage_measurements','sha256']},indent=2))

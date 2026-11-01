import json,re,statistics
from pathlib import Path
ROOT=Path(__file__).resolve().parents[2]
RAW=Path(__file__).parent/'raw';DATA=ROOT/'benchdata/competitor-stages-2026-11-01'
def run_data(label):
    j=json.loads((RAW/f'{label}.run.json').read_text());t=(RAW/f'{label}.time').read_text()
    def n(pattern):return float(re.search(pattern,t).group(1))
    assert j['exit_code']==0
    return dict(seconds=j['seconds'],cpu_seconds=n(r'User time \(seconds\): ([\d.]+)')+n(r'System time \(seconds\): ([\d.]+)'),rss_gib=n(r'Maximum resident set size \(kbytes\): (\d+)')/1024**2)
def summary(rows):return {k:dict(median=statistics.median(r[k] for r in rows),min=min(r[k] for r in rows),max=max(r[k] for r in rows)) for k in rows[0]}
report={'count':{},'count_stages':{},'build':{}}
for name,prefix in [('ours_default','ours-ours-default'),('ours_usa','ours-ours-usa'),('ours_piscem_usa','ours-piscem-usa')]:
    report['count'][name]=summary([run_data(f'{prefix}-{i}') for i in range(1,4)])
for name,prefix,stages in [('alevin_fry','af',['permit','collate','quant']),('bustools','bus',['correct','sort','count'])]:
    rows=[]
    for i in range(1,4):
        parts=[run_data(f'{prefix}-{i}-{s}') for s in stages]
        row={'seconds':json.loads((RAW/f'{prefix}-{i}-pipeline.json').read_text())['seconds'],'cpu_seconds':sum(p['cpu_seconds'] for p in parts),'rss_gib':max(p['rss_gib'] for p in parts)}
        rows.append(row)
    report['count'][name]=summary(rows)
    report['count_stages'][name]={s:summary([run_data(f'{prefix}-{i}-{s}') for i in range(1,4)]) for s in stages}
report['controls']={'bustools_gene_only_final_count':run_data('bus-gene-only'),'bustools_streamed_correct_sort':summary([run_data(f'bus-stream-{i}') for i in range(1,4)])}
for tool in ['ours','piscem','kallisto']:
    if (RAW/f'build-{tool}.run.json').exists():report['build'][tool]=run_data(f'build-{tool}')
if (RAW/'af-default-control.json').exists():
    values=json.loads((RAW/'af-default-control.json').read_text());rows=[]
    for i,value in enumerate(values,1):
        parts=[run_data(f'af-default-{i}-{stage}') for stage in ['permit','collate','quant']]
        rows.append({'seconds':value['seconds'],'cpu_seconds':sum(p['cpu_seconds'] for p in parts),'rss_gib':max(p['rss_gib'] for p in parts)})
    report['controls']['alevin_fry_native_small_cell_default']=summary(rows)
if (RAW/'bus-empty-split.run.json').exists():report['controls']['bustools_empty_split']=run_data('bus-empty-split')
if len(report['build'])==3:
    files={'ours':[DATA/'constellation.ecidx'],'piscem':sorted(DATA.glob('piscem-splici*')),'kallisto':[DATA/'kallisto-splici.idx']}
    report['index_files']={name:[{'path':str(p),'bytes':p.stat().st_size} for p in paths if p.is_file()] for name,paths in files.items()}
    report['index_bytes']={name:sum(p['bytes'] for p in paths) for name,paths in report['index_files'].items()}
report['count_output_validation']=json.loads((RAW/'validation.json').read_text())['totals']
(RAW/'summary.json').write_text(json.dumps(report,indent=2)+'\n')
print(json.dumps(report,indent=2))

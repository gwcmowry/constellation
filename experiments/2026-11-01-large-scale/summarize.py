import json,re
from pathlib import Path
ROOT=Path(__file__).resolve().parents[2];RAW=Path(__file__).parent/'raw'

def stats(label,result):
    text=(RAW/f'{label}.time').read_text()
    def n(pattern):return float(re.search(pattern,text).group(1))
    return {'seconds':result['seconds'],'cpu_seconds':n(r'User time \(seconds\): ([\d.]+)')+n(r'System time \(seconds\): ([\d.]+)'),
      'rss_gib':n(r'Maximum resident set size \(kbytes\): (\d+)')/1024**2,
      'major_page_faults':int(n(r'Major \(requiring I/O\) page faults: (\d+)')),
      'filesystem_input_blocks':int(n(r'File system inputs: (\d+)')),
      'filesystem_output_blocks':int(n(r'File system outputs: (\d+)'))}

report={'cases':{},'speedups':{},'scaling':{},'protocol':{'replicates':1,'cores':12,'prebuilt_indexes':True,'combined_times':'sums of measured native stages; extra USA count and validation excluded'}}
for size in (100,500):
    for tool in ('ours','piscem','kallisto'):
        tag=f'{size}m-{tool}';path=RAW/f'{tag}.pipeline.json'
        if not path.exists():continue
        case=json.loads(path.read_text());out={key:case[key] for key in ['pairs','map_seconds','count_seconds','pipeline_seconds','matrix_totals','count_policy']}
        out['stages']={name:stats(tag+'-'+name,result) for name,result in case['steps'].items()}
        out['million_pairs_per_second']=case['pairs']/case['map_seconds']/1e6
        out['pipeline_peak_rss_gib']=max(s['rss_gib'] for name,s in out['stages'].items() if name!='usa')
        if tool=='ours':
            for mode in ('count','usa'):
                out[mode+'_metrics']=json.loads((RAW/f'{tag}-{mode}.metrics.json').read_text())
            before=out['count_metrics']['molecules_before_umi_correction'];usa=out['usa_metrics']['molecules_after_umi_correction']
            out['usa_policy_net_gain_vs_exact_intersection']={'molecules':usa-before,'fraction':(usa-before)/before}
            out['usa_pipeline_seconds']=case['usa_pipeline_seconds']
            m=json.loads((RAW/f'{tag}-map.metrics.json').read_text());out['mapping_metrics']=m
            out['forward_assigned_reads']=round(case['pairs']*(m['unique_gene_rate']+m['ambiguous_gene_rate']))
        elif tool=='piscem':
            m=json.loads((RAW/f'{tag}-collation-generate_permit_list.json').read_text());out['barcode_correction']=m['correction_stats']
            out['forward_assigned_reads']=sum(m['correction_stats'][s] for s in ['exact_reads','corrected_reads','ambiguous_reads','not_found_reads'])
            out['mapping_metrics']=json.loads((RAW/f'{tag}-map-map_info.json').read_text())
        else:
            m=json.loads((RAW/f'{tag}-map-run_info.json').read_text());out['mapping_metrics']=m;out['forward_assigned_reads']=m['n_pseudoaligned']
        out['forward_assignment_fraction']=out['forward_assigned_reads']/case['pairs']
        report['cases'][tag]=out
    if all(f'{size}m-{t}' in report['cases'] for t in ('ours','piscem','kallisto')):
        our=report['cases'][f'{size}m-ours'];report['speedups'][str(size)]={t:{'mapping':report['cases'][f'{size}m-{t}']['map_seconds']/our['map_seconds'],'map_and_count':report['cases'][f'{size}m-{t}']['pipeline_seconds']/our['pipeline_seconds']} for t in ['piscem','kallisto']}
for tool in ('ours','piscem','kallisto'):
    if all(f'{n}m-{tool}' in report['cases'] for n in [100,500]):
        a,b=[report['cases'][f'{n}m-{tool}'] for n in [100,500]]
        report['scaling'][tool]={key:b[key]/a[key] for key in ['map_seconds','count_seconds','pipeline_seconds','pipeline_peak_rss_gib']}
(RAW/'summary.json').write_text(json.dumps(report,indent=2)+'\n')
print(json.dumps({k:v for k,v in report.items() if k!='cases'},indent=2))
for name,c in report['cases'].items():print(name,{k:c[k] for k in ['map_seconds','count_seconds','pipeline_seconds','pipeline_peak_rss_gib']})

import json,os,subprocess
from pathlib import Path
ROOT=Path(__file__).resolve().parents[2];RAW=Path(__file__).parent/'raw';D=ROOT/'benchdata/large-scale-2026-11-01/smoke'
C=ROOT/'target/release/constellation';P=ROOT/'benchdata/comparison-2026-10-27/tools/piscem-x86_64-unknown-linux-gnu/piscem';K=ROOT/'benchdata/comparison-2026-10-27/tools/kallisto/kallisto'
def run(label,args):
 cmd=['taskset','-c','0-11',*map(str,args)]
 with (RAW/f'smoke-{label}.log').open('w') as f:subprocess.run(cmd,env=dict(os.environ,RAYON_NUM_THREADS='12'),stdout=f,stderr=f,check=True)
 return cmd
cmds=[]
cmds.append(run('index-ours',[C,'index-ec','--transcripts',D/'ref.fa','--t2g-map',D/'t2g.tsv','--k','31','--syncmer-s','20','--format','ptrhash','--out',D/'ours.idx']))
cmds.append(run('index-piscem',[P,'build','-s',D/'ref.fa','-k','31','-m','19','-t','12','-o',D/'piscem','-w',D/'p-work']))
cmds.append(run('index-kallisto',[K,'index','-t','12','-k','31','-i',D/'kallisto.idx','-T',D/'k-work',D/'ref.fa']))
for label,n in [('small',2000),('large',10000)]:
 r1,r2=[D/f'{label}_R{i}.fastq.gz' for i in (1,2)]
 cmds.append(run(label+'-ours',[C,'map','--splicing-state','--index',D/'ours.idx','--r1',r1,'--r2',r2,'--out',D/f'{label}.cstrad2','--emit-metrics',D/f'{label}.json']))
 assert json.loads((D/f'{label}.json').read_text())['num_reads']==n
 cmds.append(run(label+'-piscem',[P,'map-sc','-i',D/'piscem','-g','chromium_v3','-t','12','-1',r1,'-2',r2,'-o',D/f'{label}-piscem']))
 assert json.loads((D/f'{label}-piscem/map_info.json').read_text())['num_reads']==n
 cmds.append(run(label+'-kallisto',[K,'bus','-i',D/'kallisto.idx','-x','10xv3','--fr-stranded','-t','12','-o',D/f'{label}-kallisto',r1,r2]))
 assert json.loads((D/f'{label}-kallisto/run_info.json').read_text())['n_processed']==n
(RAW/'smoke-validation.json').write_text(json.dumps({'all_three_tools_count_both_gzip_member_subsets_correctly':True,'commands':cmds},indent=2)+'\n')
print('All three native readers processed exactly 2,000 and 10,000 paired records across gzip members.')

"""Serial, pinned, subprocess wall/RSS measurements; raw phase logs are retained."""
import json, os, subprocess, sys, time
from pathlib import Path
ROOT=Path(__file__).resolve().parents[2]
DATA=ROOT/'benchdata/build-count-2026-10-31'
RAW=Path(__file__).parent/'raw'
REF=ROOT/'benchdata/reference/ensembl93/simpleaf_splici_r91_k31/ref'
ENV=dict(os.environ,RAYON_NUM_THREADS='12',CONSTELLATION_PROFILE_BUILD='1')
def run(label,args):
    cmd=['/usr/bin/time','-v','-o',str(RAW/f'{label}.time'),'taskset','-c','0-11',*map(str,args)]
    started=time.perf_counter()
    with (RAW/f'{label}.log').open('w') as log:
        p=subprocess.run(cmd,cwd=ROOT,env=ENV,stdout=log,stderr=log)
    report=dict(command=cmd,seconds=time.perf_counter()-started,exit_code=p.returncode)
    (RAW/f'{label}.run.json').write_text(json.dumps(report,indent=2)+'\n')
    print(label,round(report['seconds'],3),p.returncode,flush=True)
    if p.returncode:raise RuntimeError(label)
def count(binary,label,extra=(),assignments=None):
    run(label,[binary,'count-splicing','--assignments',assignments or ROOT/'benchdata/splicing-2026-10-30/splicing.rad','--barcode-whitelist',ROOT/'benchdata/state-comparison-2026-10-30/barcodes.txt','--out-prefix',DATA/label,'--emit-metrics',RAW/f'{label}.metrics.json',*extra])
def build(binary,label,extra=()):
    run(label,[binary,'index-ec','--transcripts',REF/'roers_ref.fa','--t2g-map',REF/'t2g_3col.tsv','--k','31','--syncmer-s','20','--out',DATA/f'{label}.ecidx',*extra])
if __name__=='__main__':
    if sys.argv[1]=='baseline':
        b=DATA/'constellation.instrumented'
        for i in range(3):count(b,f'baseline-count-{i+1}')
        build(b,'baseline-build')
        ENV['CONSTELLATION_EC_PTRHASH']='1'
        run('baseline-convert',[b,'convert-ec-index','--index',DATA/'baseline-build.ecidx','--out',DATA/'baseline-final.ecidx'])
    elif sys.argv[1]=='optimized':
        b=ROOT/'target/release/constellation'
        for i in range(3):count(b,f'optimized-count-{i+1}')
        count(b,'ours-usa',['--resolution','usa','--no-correct-umis'])
        count(b,'piscem-usa',['--resolution','usa','--no-correct-umis'],ROOT/'benchdata/state-comparison-2026-10-30/piscem-common-counter.cstrad2')
        build(b,'optimized-build-1',['--format','ptrhash'])
    elif sys.argv[1]=='repeat-build':
        b=ROOT/'target/release/constellation'
        old=DATA/'constellation.instrumented'
        for i in range(1,4):
            build(b,f'fast-build-{i}',['--format','ptrhash'])
            if i<3:
                build(old,f'baseline-build-{i+1}')
                ENV['CONSTELLATION_EC_PTRHASH']='1'
                run(f'baseline-convert-{i+1}',[old,'convert-ec-index','--index',DATA/f'baseline-build-{i+1}.ecidx','--out',DATA/f'baseline-final-{i+1}.ecidx'])
        for i in range(3):
            count(b,f'final-count-{i+1}')
            for name,binary,index in [('baseline',old,DATA/'baseline-final.ecidx'),('fast',b,DATA/'fast-build-3.ecidx')]:
                run(f'{name}-map-{i+1}',[binary,'map','--splicing-state','--index',index,
                    '--r1',ROOT/'benchdata/throughput-2026-10-27/pbmc10m_R1.fastq.gz',
                    '--r2',ROOT/'benchdata/throughput-2026-10-27/pbmc10m_R2.fastq.gz',
                    '--out',DATA/f'{name}-map.rad','--emit-metrics',RAW/f'{name}-map-{i+1}.metrics.json'])
    elif sys.argv[1]=='pipeline':
        for name,binary,extra in [('baseline',DATA/'constellation.instrumented',[]),('fast',ROOT/'target/release/constellation',['--format','ptrhash'])]:
            began=time.perf_counter()
            label=f'pipeline-{name}'
            build(binary,label+'-build',extra)
            index=DATA/f'{label}-build.ecidx'
            if name=='baseline':
                ENV['CONSTELLATION_EC_PTRHASH']='1'
                converted=DATA/f'{label}-final.ecidx'
                run(label+'-convert',[binary,'convert-ec-index','--index',index,'--out',converted])
                index=converted
            assignments=DATA/f'{label}.rad'
            run(label+'-map',[binary,'map','--splicing-state','--index',index,
                '--r1',ROOT/'benchdata/throughput-2026-10-27/pbmc10m_R1.fastq.gz',
                '--r2',ROOT/'benchdata/throughput-2026-10-27/pbmc10m_R2.fastq.gz',
                '--out',assignments,'--emit-metrics',RAW/f'{label}-map.metrics.json'])
            count(binary,label+'-count',assignments=assignments)
            report={'seconds':time.perf_counter()-began,'includes':['reference_FASTA_to_index','index_conversion_if_needed','paired_FASTQ_to_CSTRAD2','CSTRAD2_to_all_count_matrices'],'excludes':['reference_preparation','compilation','downloads'],'policy':'intersection_with_directional_UMI_correction'}
            (RAW/f'{label}.json').write_text(json.dumps(report,indent=2)+'\n')
            print(label,round(report['seconds'],3),flush=True)

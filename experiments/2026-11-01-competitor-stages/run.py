"""Matched stage measurements, serial on CPUs 0–11; no mapping in count timers."""
import hashlib,json,os,subprocess,sys,time
from pathlib import Path
ROOT=Path(__file__).resolve().parents[2]
RAW=Path(__file__).parent/'raw'
DATA=ROOT/'benchdata/competitor-stages-2026-11-01'
OLD=ROOT/'benchdata/comparison-2026-10-27'
REF=ROOT/'benchdata/reference/ensembl93/simpleaf_splici_r91_k31/ref'
STATE=ROOT/'benchdata/state-comparison-2026-10-30'
BIN=ROOT/'target/release/constellation'
PSC=OLD/'tools/piscem-x86_64-unknown-linux-gnu/piscem'
KAL=OLD/'tools/kallisto/kallisto'
BUS=DATA/'tools/bustools/bustools'
AF=STATE/'alevin-fry-x86_64-unknown-linux-gnu/alevin-fry'
ENV=dict(os.environ,RAYON_NUM_THREADS='12',OMP_NUM_THREADS='12',OPENBLAS_NUM_THREADS='1')

def run(label,args):
    cmd=['/usr/bin/time','-v','-o',str(RAW/f'{label}.time'),'taskset','-c','0-11',*map(str,args)]
    t=time.perf_counter()
    with (RAW/f'{label}.log').open('w') as log:
        if label.startswith('build-'):
            with (RAW/f'{label}.events.jsonl').open('w') as events:
                p=subprocess.Popen(cmd,cwd=ROOT,env=ENV,stdout=subprocess.PIPE,stderr=subprocess.STDOUT,text=True)
                for line in p.stdout:
                    log.write(line);log.flush()
                    events.write(json.dumps({'seconds':time.perf_counter()-t,'line':line.rstrip()})+'\n');events.flush()
                p.wait()
        else:
            p=subprocess.run(cmd,cwd=ROOT,env=ENV,stdout=log,stderr=log)
    result=dict(command=cmd,seconds=time.perf_counter()-t,exit_code=p.returncode)
    (RAW/f'{label}.run.json').write_text(json.dumps(result,indent=2)+'\n')
    print(label,round(result['seconds'],3),p.returncode,flush=True)
    if p.returncode:raise RuntimeError(label)
    return result['seconds']

def prepare():
    lines=[line.split() for line in (REF/'t2g_3col.tsv').read_text().splitlines()]
    assert all(s in ['S','U'] for _,_,s in lines)
    (DATA/'t2g.tsv').write_text(''.join(f'{t}\t{g}\n' for t,g,_ in lines))
    (DATA/'intronic-targets.txt').write_text(''.join(t+'\n' for t,_,s in lines if s=='U'))
    tools={}
    for name,path in [('constellation',BIN),('piscem',PSC),('kallisto',KAL),('bustools',BUS),('alevin-fry',AF)]:
        with path.open('rb') as f:tools[name]={'path':str(path),'sha256':hashlib.file_digest(f,'sha256').hexdigest()}
    (RAW/'tools.json').write_text(json.dumps(tools,indent=2)+'\n')
    (RAW/'reference-inputs.json').write_text(json.dumps({'fasta':str(REF/'roers_ref.fa'),'targets':len(lines),'intronic_targets':sum(s=='U' for _,_,s in lines),'whitelist':str(STATE/'barcodes.txt')},indent=2)+'\n')

def ours(i,policy='default',source='ours'):
    label=f'ours-{source}-{policy}-{i}'
    assignments=ROOT/'benchdata/splicing-2026-10-30/splicing.rad' if source=='ours' else STATE/'piscem-common-counter.cstrad2'
    extra=[] if policy=='default' else ['--resolution','usa','--no-correct-umis']
    run(label,[BIN,'count-splicing','--assignments',assignments,'--barcode-whitelist',STATE/'barcodes.txt','--out-prefix',DATA/label,'--emit-metrics',RAW/f'{label}.metrics.json',*extra])

def af(i):
    stem=f'af-{i}';out=DATA/stem;col=out/'collation';dest=out/'quant';out.mkdir(exist_ok=True)
    started=time.perf_counter();stage=[]
    stage.append(run(stem+'-permit',[AF,'generate-permit-list','-i',OLD/'pbmc10m-piscem-forward','-o',col,'-d','fw','-t','12','--valid-bc',STATE/'barcodes.txt','--cell-bc-correction','unique','--cell-bc-neighborhood','hamming-1']))
    stage.append(run(stem+'-collate',[AF,'collate','-i',col,'-r',OLD/'pbmc10m-piscem-forward','-t','12']))
    stage.append(run(stem+'-quant',[AF,'quant','-i',col,'-o',dest,'-m',REF/'t2g_3col.tsv','-r','cr-like','-t','12','--small-thresh','0']))
    (RAW/f'{stem}-pipeline.json').write_text(json.dumps({'seconds':time.perf_counter()-started,'stages_seconds':stage,'policy':'cr-like WTA exact UMI, no molecule dump'},indent=2)+'\n')

def bus(i):
    stem=f'bus-{i}';out=DATA/stem;out.mkdir(exist_ok=True);source=OLD/'pbmc10m-kallisto-forward'
    started=time.perf_counter();stage=[]
    stage.append(run(stem+'-correct',[BUS,'correct','-w',STATE/'barcodes.txt','-o',out/'corrected.bus',source/'output.bus']))
    stage.append(run(stem+'-sort',[BUS,'sort','-t','12','-m','2G','-T',out/'sort-tmp','--no-flags','-o',out/'sorted.bus',out/'corrected.bus']))
    stage.append(run(stem+'-count',[BUS,'count','-o',out/'counts','-g',DATA/'t2g.tsv','-e',source/'matrix.ec','-t',source/'transcripts.txt','--genecounts','--umi-gene','--split',DATA/'intronic-targets.txt',out/'sorted.bus']))
    (RAW/f'{stem}-pipeline.json').write_text(json.dumps({'seconds':time.perf_counter()-started,'stages_seconds':stage,'policy':'gene UMI collapsing, split S/U/A, no UMI edit correction'},indent=2)+'\n')

def builds():
    (DATA/'piscem-work').mkdir(exist_ok=True);(DATA/'kallisto-tmp').mkdir(exist_ok=True)
    run('build-ours',[BIN,'index-ec','--transcripts',REF/'roers_ref.fa','--t2g-map',REF/'t2g_3col.tsv','--k','31','--syncmer-s','20','--format','ptrhash','--out',DATA/'constellation.ecidx'])
    run('build-piscem',[PSC,'build','-s',REF/'roers_ref.fa','-k','31','-m','19','-t','12','-o',DATA/'piscem-splici','-w',DATA/'piscem-work','--ram-limit-gib','24'])
    run('build-kallisto',[KAL,'index','-t','12','-k','31','-i',DATA/'kallisto-splici.idx','-T',DATA/'kallisto-tmp',REF/'roers_ref.fa'])

if __name__=='__main__':
    prepare()
    if sys.argv[1]=='counts':
        funcs=[lambda i:ours(i),lambda i:ours(i,'usa'),af,bus,lambda i:ours(i,'usa','piscem')]
        for f in funcs:f('warmup')
        for i in range(1,4):
            for f in funcs[i:]+funcs[:i]:f(i)
    elif sys.argv[1]=='builds':builds()

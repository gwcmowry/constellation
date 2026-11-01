"""Serial 100M/500M paired-read, native-output scaling experiment."""
import hashlib,json,os,shutil,signal,subprocess,sys,threading,time
from pathlib import Path
ROOT=Path(__file__).resolve().parents[2];RAW=Path(__file__).parent/'raw';DATA=ROOT/'benchdata/large-scale-2026-11-01'
PRIOR=ROOT/'benchdata/competitor-stages-2026-11-01';OLD=ROOT/'benchdata/comparison-2026-10-27';STATE=ROOT/'benchdata/state-comparison-2026-10-30';REF=ROOT/'benchdata/reference/ensembl93/simpleaf_splici_r91_k31/ref'
C=DATA/'constellation.frozen';P=OLD/'tools/piscem-x86_64-unknown-linux-gnu/piscem';K=OLD/'tools/kallisto/kallisto';B=PRIOR/'tools/bustools/bustools';A=STATE/'alevin-fry-x86_64-unknown-linux-gnu/alevin-fry'
ENV=dict(os.environ,RAYON_NUM_THREADS='12',OMP_NUM_THREADS='12',OPENBLAS_NUM_THREADS='1')
LIMIT=46*1024**3

def write_progress(**fields):
    (RAW/'progress.json').write_text(json.dumps(dict(utc=time.strftime('%Y-%m-%dT%H:%M:%SZ',time.gmtime()),**fields),indent=2)+'\n')

def rss_tree(pid):
    pending=[pid];total=0
    while pending:
        p=pending.pop()
        try:
            status=Path(f'/proc/{p}/status').read_text()
            for line in status.splitlines():
                if line.startswith('VmRSS:'):total+=int(line.split()[1])*1024;break
            pending.extend(map(int,Path(f'/proc/{p}/task/{p}/children').read_text().split()))
        except (FileNotFoundError,ProcessLookupError):pass
    return total

def run(label,args):
    if shutil.disk_usage(DATA).free<8*1024**3:raise RuntimeError('Less than 8 GiB free before '+label)
    cmd=['/usr/bin/time','-v','-o',str(RAW/f'{label}.time'),'taskset','-c','0-11',*map(str,args)]
    saved=RAW/f'{label}.run.json'
    if saved.exists():
        prior=json.loads(saved.read_text())
        if prior['exit_code']==0:
            assert prior['command']==cmd,'Completed stage command changed: '+label
            print(label,'reuse completed stage',flush=True);return prior
    started=time.perf_counter();events=[];sampled_peak=0;guard=False;stopped=threading.Event()
    write_progress(stage=label,status='running',command=cmd)
    with (RAW/f'{label}.log').open('w') as log,(RAW/f'{label}.events.jsonl').open('w') as event_file:
        p=subprocess.Popen(cmd,cwd=ROOT,env=ENV,stdout=subprocess.PIPE,stderr=subprocess.STDOUT,text=True,start_new_session=True)
        def monitor():
            nonlocal sampled_peak,guard
            with (RAW/f'{label}.memory.jsonl').open('w') as f:
                while p.poll() is None:
                    resident=rss_tree(p.pid);sampled_peak=max(sampled_peak,resident)
                    f.write(json.dumps({'seconds':time.perf_counter()-started,'rss_bytes':resident})+'\n');f.flush()
                    if resident>LIMIT:
                        guard=True;os.killpg(p.pid,signal.SIGTERM);return
                    stopped.wait(2)
        watcher=threading.Thread(target=monitor,daemon=True);watcher.start()
        for line in p.stdout:
            log.write(line);log.flush()
            event={'seconds':time.perf_counter()-started,'line':line.rstrip()};events.append(event)
            event_file.write(json.dumps(event)+'\n');event_file.flush()
        rc=p.wait();elapsed=time.perf_counter()-started;stopped.set();watcher.join()
    result={'command':cmd,'seconds':elapsed,'exit_code':rc,'sampled_peak_rss_bytes':sampled_peak,'memory_guard_gib':46,'stopped_by_memory_guard':guard}
    if '-map' in label and 'kallisto' in label:
        begin=end=None
        for e in events:
            if '[quant] finding pseudoalignments for the reads' in e['line']:begin=e['seconds']
            if begin is not None and end is None and e['line'].strip()=='done':end=e['seconds']
        if begin is not None and end is not None:result.update(before_read_loop_seconds=begin,read_loop_seconds=end-begin,after_read_loop_seconds=elapsed-end)
    (RAW/f'{label}.run.json').write_text(json.dumps(result,indent=2)+'\n')
    write_progress(stage=label,status='complete' if rc==0 else 'failed',seconds=elapsed,exit_code=rc)
    print(label,round(elapsed,3),rc,flush=True)
    if rc:raise RuntimeError(f'{label} failed; memory guard={guard}')
    return result

def normalize(label,kind,path):
    write_progress(stage=label+'-validation',status='running')
    with (RAW/f'{label}.validation.log').open('w') as f:
        subprocess.run(['taskset','-c','0-11',sys.executable,str(RAW.parent/'matrices.py'),label,kind,str(path)],env=ENV,stdout=f,stderr=f,check=True)
    return json.loads((RAW/f'{label}.matrix-summary.json').read_text())

def archive_json(prefix,directory):
    for p in directory.rglob('*.json'):
        name=str(p.relative_to(directory)).replace('/','-')
        shutil.copy2(p,RAW/f'{prefix}-{name}')

def cleanup(label,work):
    assert work.resolve().is_relative_to((DATA/'work').resolve())
    files=[{'path':str(p.relative_to(work)),'bytes':p.stat().st_size} for p in work.rglob('*') if p.is_file()]
    (RAW/f'{label}.artifacts.json').write_text(json.dumps({'removed_scratch_root':str(work),'files':files,'bytes':sum(x['bytes'] for x in files)},indent=2)+'\n')
    shutil.rmtree(work)

def warm_index(tool):
    files={'ours':[PRIOR/'constellation.ecidx'],'piscem':[p for p in sorted(PRIOR.glob('piscem-splici*')) if p.is_file()],'kallisto':[PRIOR/'kallisto-splici.idx']}[tool]
    started=time.perf_counter();count=0;buffer=bytearray(8<<20)
    for path in files:
        with path.open('rb',buffering=0) as f:
            while n:=f.readinto(buffer):count+=n
    return {'seconds':time.perf_counter()-started,'bytes':count,'files':list(map(str,files)),'included_in_mapping_time':False}

def pipeline(size,tool):
    n=size*1_000_000;tag=f'{size}m-{tool}';work=DATA/'work'/tag
    if (RAW/f'{tag}.pipeline.json').exists():
        print(tag,'already complete',flush=True);return
    work.mkdir(parents=True,exist_ok=True)
    data='small' if size==100 else 'large';r1,r2=[DATA/f'inputs/{data}_R{i}.fastq.gz' for i in (1,2)]
    result={'pairs':n,'tool':tool,'steps':{},'count_policy':{},'matrix_totals':{},'index_warmup':warm_index(tool)}
    if tool=='ours':
        result['steps']['map']=run(tag+'-map',[C,'map','--splicing-state','--index',PRIOR/'constellation.ecidx','--r1',r1,'--r2',r2,'--out',work/'map.cstrad2','--emit-metrics',RAW/f'{tag}-map.metrics.json'])
        assert json.loads((RAW/f'{tag}-map.metrics.json').read_text())['num_reads']==n
        for mode,extra in [('count',[]),('usa',['--resolution','usa','--no-correct-umis'])]:
            label=tag+'-'+mode
            result['steps'][mode]=run(label,[C,'count-splicing','--assignments',work/'map.cstrad2','--barcode-whitelist',STATE/'barcodes.txt','--out-prefix',work/mode,'--emit-metrics',RAW/f'{label}.metrics.json',*extra])
            metrics=json.loads((RAW/f'{label}.metrics.json').read_text());assert metrics['rows']==n
            assert metrics['molecules_after_umi_correction']<=metrics['molecules_before_umi_correction']<=metrics['raw_barcode_umi_groups']<=metrics['mapped_valid_rows']<=n
            result['matrix_totals'][mode]=normalize(label,'ours',work/mode)
            result['count_policy'][mode]=metrics['resolution']
        count_steps=['count']
    elif tool=='piscem':
        result['steps']['map']=run(tag+'-map',[P,'map-sc','-i',PRIOR/'piscem-splici','-g','chromium_v3','-t','12','-1',r1,'-2',r2,'-o',work/'map'])
        assert json.loads((work/'map/map_info.json').read_text())['num_reads']==n
        col=work/'collation';quant=work/'quant'
        result['steps']['permit']=run(tag+'-permit',[A,'generate-permit-list','-i',work/'map','-o',col,'-d','fw','-t','12','--valid-bc',STATE/'barcodes.txt','--cell-bc-correction','unique','--cell-bc-neighborhood','hamming-1'])
        result['steps']['collate']=run(tag+'-collate',[A,'collate','-i',col,'-r',work/'map','-t','12'])
        result['steps']['quant']=run(tag+'-quant',[A,'quant','-i',col,'-o',quant,'-m',REF/'t2g_3col.tsv','-r','cr-like','-t','12','--small-thresh','0'])
        q=json.loads((quant/'quant.json').read_text())['quant_options'];assert q['crlike_umi_edit']==0 and q['sa_model']=='WinnerTakeAll' and not q['dump_molecules']
        result['matrix_totals']['count']=normalize(tag+'-count','af',quant);result['count_policy']['count']='cr-like WTA exact UMI'
        count_steps=['permit','collate','quant']
    else:
        native=work/'map'
        result['steps']['map']=run(tag+'-map',[K,'bus','-i',PRIOR/'kallisto-splici.idx','-x','10xv3','--fr-stranded','-t','12','-o',native,r1,r2])
        assert json.loads((native/'run_info.json').read_text())['n_processed']==n
        result['steps']['correct']=run(tag+'-correct',[B,'correct','-w',STATE/'barcodes.txt','-o',work/'corrected.bus',native/'output.bus'])
        result['steps']['sort']=run(tag+'-sort',[B,'sort','-t','12','-m','2G','-T',work/'sort-tmp','--no-flags','-o',work/'sorted.bus',work/'corrected.bus'])
        result['steps']['count']=run(tag+'-count',[B,'count','-o',work/'counts','-g',PRIOR/'t2g.tsv','-e',native/'matrix.ec','-t',native/'transcripts.txt','--genecounts','--umi-gene','--split',PRIOR/'intronic-targets.txt',work/'sorted.bus'])
        result['matrix_totals']['count']=normalize(tag+'-count','bus',work/'counts');result['count_policy']['count']='gene UMI split, exact UMIs'
        count_steps=['correct','sort','count']
    for values in result['matrix_totals'].values():assert values['total']<=n
    result['count_seconds']=sum(result['steps'][s]['seconds'] for s in count_steps)
    result['map_seconds']=result['steps']['map']['seconds'];result['pipeline_seconds']=result['map_seconds']+result['count_seconds']
    if tool=='ours':result['usa_pipeline_seconds']=result['map_seconds']+result['steps']['usa']['seconds']
    (RAW/f'{tag}.pipeline.json').write_text(json.dumps(result,indent=2)+'\n');archive_json(tag,work);cleanup(tag,work)
    write_progress(stage=tag,status='pipeline-complete',pipeline_seconds=result['pipeline_seconds'])

def prepare_environment():
    p=DATA/'inputs/paired-validation.json';j=json.loads(p.read_text());assert j['large_pairs']==500_000_000 and j['small_pairs']==100_000_000 and j['r1_r2_ids_match']
    assert p.stat().st_mtime>=max(f.stat().st_mtime for f in (DATA/'inputs').glob('*.gz'))
    shutil.copy2(p,RAW/'input-validation.json');shutil.copy2(ROOT/'target/release/constellation',C)
    binaries={}
    for name,path in [('ours',C),('piscem',P),('kallisto',K),('bustools',B),('alevin-fry',A)]:
        with path.open('rb') as f:binaries[name]={'path':str(path),'sha256':hashlib.file_digest(f,'sha256').hexdigest()}
    inputs={str(f):{'bytes':f.stat().st_size} for f in (DATA/'inputs').glob('*.gz')}
    (RAW/'provenance.json').write_text(json.dumps({'binaries':binaries,'inputs':inputs,'cpus':'0-11','threads':12,'replicates':1,'memory_guard_gib':46,'index_build_included':False,'gzip_level':1},indent=2)+'\n')
    (DATA/'work').mkdir(exist_ok=True)

if __name__=='__main__':
    if sys.argv[1]=='init':prepare_environment()
    elif sys.argv[1]=='one':pipeline(int(sys.argv[2]),sys.argv[3])
    elif sys.argv[1]=='all':
        prepare_environment()
        for size in (100,500):
            for tool in ('ours','piscem','kallisto'):
                try:pipeline(size,tool)
                except Exception as e:
                    (RAW/f'{size}m-{tool}.failure.json').write_text(json.dumps({'error':repr(e)},indent=2)+'\n')
                    write_progress(stage=f'{size}m-{tool}',status='failed',error=repr(e));raise
        write_progress(stage='all',status='complete')

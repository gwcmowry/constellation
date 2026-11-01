import hashlib,json,sys
from run import run,BUS,OLD,DATA,RAW
base=OLD/'pbmc10m-kallisto-forward'
run('bus-gene-only',[BUS,'count','-o',DATA/'bus-gene-only','-g',DATA/'t2g.tsv','-e',base/'matrix.ec','-t',base/'transcripts.txt','--genecounts','--umi-gene',DATA/'bus-1/sorted.bus'])
def sha(p):
    with p.open('rb') as f:return hashlib.file_digest(f,'sha256').hexdigest()
expected=sha(DATA/'bus-1/sorted.bus');hashes={}
for i in range(1,4):
    label=f'bus-stream-{i}'
    run(label,[sys.executable,str(RAW.parent/'stream_bus.py'),label])
    hashes[label]=sha(DATA/label/'sorted.bus');assert hashes[label]==expected
(RAW/'stream-validation.json').write_text(json.dumps({'sorted_bus_sha256':expected,'streamed_output_sha256':hashes},indent=2)+'\n')

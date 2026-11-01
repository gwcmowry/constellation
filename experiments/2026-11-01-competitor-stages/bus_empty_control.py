"""Separate fixed split-metadata setup from processing actual molecule records."""
import json,struct
from run import run,BUS,OLD,DATA,RAW
source=OLD/'pbmc10m-kallisto-forward'
with (DATA/'bus-1/sorted.bus').open('rb') as f:
    header=f.read(20);magic,v,b,u,n=struct.unpack('<4s4I',header)
    assert (magic,v,b,u)==(b'BUS\0',1,16,12)
    header+=f.read(n)
empty=DATA/'empty.bus';empty.write_bytes(header)
run('bus-empty-split',[BUS,'count','-o',DATA/'bus-empty','-g',DATA/'t2g.tsv','-e',source/'matrix.ec','-t',source/'transcripts.txt','--genecounts','--umi-gene','--split',DATA/'intronic-targets.txt',empty])
(RAW/'empty-control.json').write_text(json.dumps({'bus_records':0,'same_ecmap':str(source/'matrix.ec'),'same_split_list':str(DATA/'intronic-targets.txt')},indent=2)+'\n')

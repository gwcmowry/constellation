"""Stream correct -> sort within a 12-worker budget (1 + 11)."""
import subprocess,sys
from run import BUS,OLD,STATE,DATA,RAW
label=sys.argv[1];out=DATA/label;out.mkdir(exist_ok=True)
with (RAW/f'{label}-correct.log').open('w') as c, (RAW/f'{label}-sort.log').open('w') as s:
    first=subprocess.Popen([str(BUS),'correct','-w',str(STATE/'barcodes.txt'),'-p',str(OLD/'pbmc10m-kallisto-forward/output.bus')],stdout=subprocess.PIPE,stderr=c)
    second=subprocess.Popen([str(BUS),'sort','-t','11','-m','2G','-T',str(out/'tmp'),'--no-flags','-o',str(out/'sorted.bus'),'-'],stdin=first.stdout,stderr=s)
    first.stdout.close()
    a=second.wait();b=first.wait()
    assert a==b==0,(a,b)

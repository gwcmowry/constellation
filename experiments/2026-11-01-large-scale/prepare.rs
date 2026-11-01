//! Bounded parallel gzip preparation, preserving every selected FASTQ byte.
use flate2::{read::MultiGzDecoder,write::GzEncoder,Compression};
use sha2::{Digest,Sha256};
use serde_json::json;
use std::{collections::BTreeMap,fs::File,io::{BufRead,BufReader,BufWriter,Write},path::{Path,PathBuf},sync::{mpsc,Arc,Mutex},thread,time::Instant};
struct Job { mate:usize, order:u64, small:bool, bytes:Vec<u8> }
struct Done { mate:usize, order:u64, small:bool, bytes:Vec<u8> }
fn hex(hash:Sha256)->String {hash.finalize().iter().map(|x|format!("{x:02x}")).collect()}
fn read_mate(source:&Path,out:&Path,mate:usize,large:u64,small:u64,block_size:usize,send:mpsc::SyncSender<Job>) {
    let mut order=0;let mut all_hash=Sha256::new();let mut small_hash=Sha256::new();let mut lanes=Vec::new();
    let started=Instant::now();
    for lane in 1..=2 {
        let path=source.join(format!("pbmc_10k_v3_S1_L00{lane}_R{}_001.fastq.gz",mate+1));
        let mut input=BufReader::with_capacity(1<<20,MultiGzDecoder::new(BufReader::with_capacity(1<<20,File::open(&path).unwrap())));
        let mut bytes=Vec::with_capacity(block_size+1024);let mut block_small=true;
        let mut first=String::new();let mut last=String::new();
        for i in 0..large {
            if i==small && !bytes.is_empty() {
                send.send(Job{mate,order,small:block_small,bytes}).unwrap();order+=1;
                bytes=Vec::with_capacity(block_size+1024);block_small=false;
            }
            let start=bytes.len();
            assert!(input.read_until(b'\n',&mut bytes).unwrap()>0,"source exhausted: {} after {i}",path.display());
            assert_eq!(bytes[start],b'@');
            let header_end=bytes.len();
            let id_end=start+bytes[start..header_end].iter().position(|b|b.is_ascii_whitespace()).unwrap_or(header_end-start);
            all_hash.update(&bytes[start..id_end]);all_hash.update([0]);
            if i<small {small_hash.update(&bytes[start..id_end]);small_hash.update([0]);}
            if i==0 {first=String::from_utf8(bytes[start..id_end].to_vec()).unwrap();}
            if i+1==large {last=String::from_utf8(bytes[start..id_end].to_vec()).unwrap();}
            assert!(input.read_until(b'\n',&mut bytes).unwrap()>0);let seq_end=bytes.len();
            assert!(input.read_until(b'\n',&mut bytes).unwrap()>0);let plus_end=bytes.len();
            assert_eq!(bytes[seq_end],b'+');
            assert!(input.read_until(b'\n',&mut bytes).unwrap()>0);
            assert_eq!(seq_end-header_end,bytes.len()-plus_end,"FASTQ sequence/quality length mismatch");
            if bytes.len()>=block_size {
                send.send(Job{mate,order,small:block_small,bytes}).unwrap();order+=1;
                bytes=Vec::with_capacity(block_size+1024);block_small=i+1<small;
            }
            if (i+1)%10_000_000==0 {eprintln!("prepared R{} L00{lane} {} records {:.1}s",mate+1,i+1,started.elapsed().as_secs_f64());}
        }
        if !bytes.is_empty(){send.send(Job{mate,order,small:block_small,bytes}).unwrap();order+=1;}
        lanes.push(json!({"path":path,"records":large,"first_id":first,"last_id":last}));
    }
    let report=json!({"mate":mate+1,"large_records":2*large,"small_records":2*small,"large_id_sha256":hex(all_hash),"small_id_sha256":hex(small_hash),"gzip_members":order,"lanes":lanes,"seconds":started.elapsed().as_secs_f64()});
    std::fs::write(out.join(format!("prepare-R{}.json",mate+1)),serde_json::to_vec_pretty(&report).unwrap()).unwrap();
}
fn main(){
    let a:Vec<_>=std::env::args().collect();assert_eq!(a.len(),6,"source output large_per_lane small_per_lane block_bytes");
    let source=PathBuf::from(&a[1]);let out=PathBuf::from(&a[2]);std::fs::create_dir_all(&out).unwrap();
    let large:u64=a[3].parse().unwrap();let small:u64=a[4].parse().unwrap();let block_size:usize=a[5].parse().unwrap();assert!(small>0 && small<=large);
    let (send,recv)=mpsc::sync_channel::<Job>(20);let recv=Arc::new(Mutex::new(recv));
    let (result,result_recv)=mpsc::sync_channel::<Done>(20);
    let pool=rayon::ThreadPoolBuilder::new().num_threads(10).build().unwrap();
    for _ in 0..10 {let recv=recv.clone();let result=result.clone();pool.spawn(move|| {
        loop {let job=match recv.lock().unwrap().recv(){Ok(j)=>j,Err(_)=>break};
            let mut encoder=GzEncoder::new(Vec::with_capacity(job.bytes.len()/3),Compression::fast());
            encoder.write_all(&job.bytes).unwrap();let bytes=encoder.finish().unwrap();
            result.send(Done{mate:job.mate,order:job.order,small:job.small,bytes}).unwrap();
        }
    });}
    drop(result);
    let readers:Vec<_>=(0..2).map(|mate|{let source=source.clone();let out=out.clone();let send=send.clone();thread::spawn(move||read_mate(&source,&out,mate,large,small,block_size,send))}).collect();drop(send);
    let mut outputs:Vec<_>=(0..2).map(|mate|[BufWriter::with_capacity(1<<20,File::create(out.join(format!("small_R{}.fastq.gz",mate+1))).unwrap()),BufWriter::with_capacity(1<<20,File::create(out.join(format!("large_R{}.fastq.gz",mate+1))).unwrap())]).collect();
    let mut expected=[0u64;2];let mut pending:[BTreeMap<u64,Done>;2]=std::array::from_fn(|_|BTreeMap::new());
    for done in result_recv {
        let mate=done.mate;assert!(pending[mate].insert(done.order,done).is_none());
        while let Some(done)=pending[mate].remove(&expected[mate]) {
            outputs[mate][1].write_all(&done.bytes).unwrap();if done.small {outputs[mate][0].write_all(&done.bytes).unwrap();} expected[mate]+=1;
        }
    }
    for reader in readers {reader.join().unwrap();}drop(pool);
    assert!(pending.iter().all(BTreeMap::is_empty));for pair in &mut outputs {for file in pair {file.flush().unwrap();}}
    let reports:Vec<serde_json::Value>=(1..=2).map(|m|serde_json::from_reader(File::open(out.join(format!("prepare-R{m}.json"))).unwrap()).unwrap()).collect();
    for key in ["large_records","small_records","large_id_sha256","small_id_sha256"]{assert_eq!(reports[0][key],reports[1][key],"mate identity/count mismatch {key}");}
    std::fs::write(out.join("paired-validation.json"),serde_json::to_vec_pretty(&json!({"r1_r2_ids_match":true,"large_pairs":2*large,"small_pairs":2*small,"selection":"first N pairs per lane, both lanes jointly; small subset nested in large","gzip_level":1,"reports":reports})).unwrap()).unwrap();
}

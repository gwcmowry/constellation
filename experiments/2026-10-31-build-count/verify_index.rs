// Exhaustive retained-key comparison plus one million deterministic probes.
use constellation_core::ec_index::LoadedEcIndex;
use memmap2::Mmap;
use rayon::prelude::*;
use std::fs::File;
fn word(bytes:&[u8], offset:usize)->u64 {u64::from_le_bytes(bytes[offset..offset+8].try_into().unwrap())}
fn main() {
    let args:Vec<_>=std::env::args().collect();
    let old=LoadedEcIndex::load(&args[1]).unwrap();
    let new=LoadedEcIndex::load(&args[2]).unwrap();
    assert_eq!(old.num_genes(),new.num_genes());
    for id in 0..old.num_genes(){assert_eq!(old.gene_name(id as u32),new.gene_name(id as u32));}
    let file=File::open(&args[1]).unwrap();
    let bytes=unsafe {Mmap::map(&file).unwrap()};
    assert_eq!(u32::from_le_bytes(bytes[8..12].try_into().unwrap()),9);
    let expected=word(&bytes,32) as usize;
    let count=word(&bytes,40) as usize;
    let offset=word(&bytes,88) as usize;
    let checked:usize=bytes[offset..offset+16*count].par_chunks(16*65536).map(|chunk|{
        let mut n=0;
        for record in chunk.chunks_exact(16) {
            if u16::from_le_bytes(record[14..16].try_into().unwrap()) == 0 {continue;}
            let code=word(record,0);
            assert_eq!(old.lookup_hit(code),new.lookup_hit(code),"key {code}"); n+=1;
        }
        n
    }).sum();
    assert_eq!(checked,expected);
    let new_file=File::open(&args[2]).unwrap();
    let new_bytes=unsafe {Mmap::map(&new_file).unwrap()};
    assert_eq!(word(&new_bytes,32) as usize,expected);
    let new_count=word(&new_bytes,40) as usize;
    let new_offset=word(&new_bytes,88) as usize;
    let occupied:usize=new_bytes[new_offset..new_offset+16*new_count].par_chunks(16*65536)
        .map(|chunk| chunk.chunks_exact(16).filter(|r| u16::from_le_bytes(r[14..16].try_into().unwrap())>0).count()).sum();
    assert_eq!(occupied,expected);
    let mut state=0xdeadbeef12345678u64;
    let k=u32::from_le_bytes(bytes[12..16].try_into().unwrap());
    let mask=if k==32{u64::MAX}else{(1u64<<(2*k))-1};
    for _ in 0..1_000_000 {
        state^=state<<13; state^=state>>7; state^=state<<17;
        let key=state&mask; assert_eq!(old.lookup_hit(key),new.lookup_hit(key));
    }
    println!("{{\"retained_keys_verified\":{checked},\"extra_probes\":1000000,\"gene_dictionaries_equal\":true}}");
}

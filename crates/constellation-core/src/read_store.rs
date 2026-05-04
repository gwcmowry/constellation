use crate::dna::{decode_to_vec, encode_acgt, EncodedSeq};
use crate::ReadId;

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct CompactRead {
    pub read_id: ReadId,
    pub seq_offset: u64,
    pub len: u16,
    pub qual_offset: u64,
    pub cb: u64,
    pub umi: u64,
    pub flags: u16,
}

#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub struct ReadStore {
    pub reads: Vec<CompactRead>,
    pub seq_arena: Vec<u8>,
    pub qual_arena: Vec<u8>,
}

impl ReadStore {
    pub fn push(
        &mut self,
        read_id: ReadId,
        seq: &[u8],
        qual: &[u8],
        cb: u64,
        umi: u64,
        flags: u16,
    ) {
        let seq_offset = self.seq_arena.len() as u64;
        let qual_offset = self.qual_arena.len() as u64;
        self.seq_arena
            .extend(seq.iter().map(|b| b.to_ascii_uppercase()));
        self.qual_arena.extend(qual);
        self.reads.push(CompactRead {
            read_id,
            seq_offset,
            len: seq.len().min(u16::MAX as usize) as u16,
            qual_offset,
            cb,
            umi,
            flags,
        });
    }

    pub fn seq(&self, read: &CompactRead) -> &[u8] {
        let start = read.seq_offset as usize;
        let end = start + read.len as usize;
        &self.seq_arena[start..end]
    }

    pub fn qual(&self, read: &CompactRead) -> &[u8] {
        let start = read.qual_offset as usize;
        let end = start + read.len as usize;
        &self.qual_arena[start..end]
    }

    pub fn encoded_seq(&self, read: &CompactRead) -> EncodedSeq {
        encode_acgt(self.seq(read))
    }
}

pub fn packed_roundtrip(seq: &[u8]) -> Vec<u8> {
    decode_to_vec(&encode_acgt(seq))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn read_store_roundtrip() {
        let mut store = ReadStore::default();
        store.push(7, b"acgtn", b"IIIII", 1, 2, 0);
        let read = &store.reads[0];
        assert_eq!(store.seq(read), b"ACGTN");
        assert_eq!(store.qual(read), b"IIIII");
        assert_eq!(packed_roundtrip(store.seq(read)), b"ACGTN");
    }
}

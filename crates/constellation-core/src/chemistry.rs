use crate::dna::Base;
use thiserror::Error;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Chemistry {
    Tenx3pV3,
}

impl Chemistry {
    pub fn parse(name: &str) -> Result<Self, ChemistryError> {
        match name {
            "tenx-3p-v3" | "10x-3p-v3" => Ok(Self::Tenx3pV3),
            other => Err(ChemistryError::Unknown(other.to_owned())),
        }
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct BarcodeUmi {
    pub cell_barcode: u64,
    pub umi: u64,
    pub cell_barcode_seq: String,
    pub umi_seq: String,
}

#[derive(Debug, Error)]
pub enum ChemistryError {
    #[error("unknown chemistry: {0}")]
    Unknown(String),
    #[error("R1 is too short for 10x 3' v3: expected at least 28 bases, got {0}")]
    TooShort(usize),
    #[error("R1 barcode/UMI contains a non-ACGT base")]
    InvalidBase,
}

pub fn parse_tenx_3p_v3_r1(seq: &[u8]) -> Result<BarcodeUmi, ChemistryError> {
    if seq.len() < 28 {
        return Err(ChemistryError::TooShort(seq.len()));
    }
    let cb = &seq[..16];
    let umi = &seq[16..28];
    Ok(BarcodeUmi {
        cell_barcode: pack_acgt(cb)?,
        umi: pack_acgt(umi)?,
        cell_barcode_seq: String::from_utf8(cb.to_ascii_uppercase()).expect("ACGT is UTF-8"),
        umi_seq: String::from_utf8(umi.to_ascii_uppercase()).expect("ACGT is UTF-8"),
    })
}

fn pack_acgt(seq: &[u8]) -> Result<u64, ChemistryError> {
    let mut out = 0_u64;
    for &b in seq {
        out = (out << 2) | Base::from_ascii(b).ok_or(ChemistryError::InvalidBase)? as u64;
    }
    Ok(out)
}

#[cfg(test)]
mod tests {
    #[test]
    fn parse_tenx_3p_v3_r1() {
        let parsed = super::parse_tenx_3p_v3_r1(b"ACGTACGTACGTACGTTTTTCCCCAAAAGGGG").unwrap();
        assert_eq!(parsed.cell_barcode_seq, "ACGTACGTACGTACGT");
        assert_eq!(parsed.umi_seq, "TTTTCCCCAAAA");
    }
}

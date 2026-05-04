use crate::fastq::{write_fastq, FastqRecord};
use crate::index_build::read_fasta;
use std::fs;
use std::path::{Path, PathBuf};
use thiserror::Error;

#[derive(Debug, Clone)]
pub struct SimulateConfig {
    pub transcripts: PathBuf,
    pub num_reads: usize,
    pub read_len: usize,
    pub error_rate: f64,
    pub scenario: SimulationScenario,
    pub out_prefix: PathBuf,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SimulationScenario {
    Exact,
    LowComplexity,
    HighExpressionSkew,
}

#[derive(Debug, Error)]
pub enum SimulateError {
    #[error("simulation input error: {0}")]
    Build(#[from] crate::index_build::BuildIndexError),
    #[error("simulation I/O error: {0}")]
    Io(#[from] std::io::Error),
    #[error("no transcript is long enough for read length {0}")]
    NoEligibleTranscript(usize),
}

pub fn simulate(config: &SimulateConfig) -> Result<(), SimulateError> {
    let transcripts = read_fasta(&config.transcripts)?;
    if !transcripts
        .iter()
        .any(|record| record.seq.len() >= config.read_len)
    {
        return Err(SimulateError::NoEligibleTranscript(config.read_len));
    }

    let mut rng = Lcg::new(0xC057_E11A_710A_u64);
    let mut r1 = Vec::with_capacity(config.num_reads);
    let mut r2 = Vec::with_capacity(config.num_reads);
    let mut truth = String::from("read_id\ttranscript\tgene\tpos\n");
    let eligible: Vec<_> = transcripts
        .iter()
        .filter(|record| record.seq.len() >= config.read_len)
        .collect();

    for read_id in 0..config.num_reads {
        let tx_idx = match config.scenario {
            SimulationScenario::HighExpressionSkew if eligible.len() > 1 && read_id % 10 < 8 => 0,
            _ => read_id % eligible.len(),
        };
        let tx = eligible[tx_idx];
        let max_start = tx.seq.len() - config.read_len;
        let start = if max_start == 0 {
            0
        } else {
            (rng.next_u64() as usize) % (max_start + 1)
        };
        let mut read_seq =
            if config.scenario == SimulationScenario::LowComplexity && read_id % 5 == 0 {
                vec![b'A'; config.read_len]
            } else {
                tx.seq[start..start + config.read_len].to_vec()
            };
        for base in &mut read_seq {
            if rng.next_f64() < config.error_rate {
                *base = substitute(*base, &mut rng);
            }
        }
        let id = format!("sim_read_{read_id}");
        r1.push(FastqRecord {
            id: id.clone(),
            seq: random_acgt(28, &mut rng),
            qual: vec![b'I'; 28],
        });
        r2.push(FastqRecord {
            id,
            seq: read_seq,
            qual: vec![b'I'; config.read_len],
        });
        truth.push_str(&format!("{read_id}\t{}\t{}\t{start}\n", tx.name, tx.gene));
    }

    write_fastq(with_suffix(&config.out_prefix, "_R1.fastq"), &r1)?;
    write_fastq(with_suffix(&config.out_prefix, "_R2.fastq"), &r2)?;
    fs::write(with_suffix(&config.out_prefix, "_truth.tsv"), truth)?;
    Ok(())
}

impl SimulationScenario {
    pub fn parse(value: &str) -> Option<Self> {
        match value {
            "exact" => Some(Self::Exact),
            "low-complexity" => Some(Self::LowComplexity),
            "high-expression-skew" => Some(Self::HighExpressionSkew),
            _ => None,
        }
    }
}

fn with_suffix(prefix: &Path, suffix: &str) -> PathBuf {
    PathBuf::from(format!("{}{}", prefix.display(), suffix))
}

fn random_acgt(len: usize, rng: &mut Lcg) -> Vec<u8> {
    (0..len)
        .map(|_| match rng.next_u64() & 0b11 {
            0 => b'A',
            1 => b'C',
            2 => b'G',
            _ => b'T',
        })
        .collect()
}

fn substitute(base: u8, rng: &mut Lcg) -> u8 {
    let choices = match base.to_ascii_uppercase() {
        b'A' => [b'C', b'G', b'T'],
        b'C' => [b'A', b'G', b'T'],
        b'G' => [b'A', b'C', b'T'],
        b'T' => [b'A', b'C', b'G'],
        _ => [b'A', b'C', b'G'],
    };
    choices[(rng.next_u64() as usize) % choices.len()]
}

#[derive(Debug, Clone, Copy)]
struct Lcg(u64);

impl Lcg {
    fn new(seed: u64) -> Self {
        Self(seed)
    }

    fn next_u64(&mut self) -> u64 {
        self.0 = self
            .0
            .wrapping_mul(6364136223846793005)
            .wrapping_add(1442695040888963407);
        self.0
    }

    fn next_f64(&mut self) -> f64 {
        (self.next_u64() >> 11) as f64 / ((1_u64 << 53) as f64)
    }
}

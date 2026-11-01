//! Versioned, self-contained gene/state compatibility records (CSTRAD2).
use anyhow::{anyhow, bail, ensure, Context, Result};
use constellation_core::ec_index::{LoadedEcIndex, SplicingState};
use flate2::read::MultiGzDecoder;
use std::fs::File;
use std::io::{BufReader, Read, Seek, Write};
use std::path::Path;
use std::process::{Child, Command, Stdio};

const MAGIC: &[u8; 8] = b"CSTRAD2\0";

pub fn write_header(writer: &mut impl Write, index: &LoadedEcIndex) -> Result<()> {
    writer.write_all(MAGIC)?;
    for n in [2u32, 28, 1, 5, u32::try_from(index.num_genes())?] {
        writer.write_all(&n.to_le_bytes())?;
    }
    for gene in 0..index.num_genes() {
        let name = index.gene_name(gene as u32).context("missing gene name")?;
        writer.write_all(&u16::try_from(name.len())?.to_le_bytes())?;
        writer.write_all(name.as_bytes())?;
    }
    Ok(())
}

#[derive(Default)]
pub struct StateRecord {
    pub barcode_umi: u64,
    pub valid_mask: u32,
    pub assignment_type: u8,
    pub candidates: Vec<(u32, SplicingState)>,
}

pub struct StateReader {
    reader: BufReader<Box<dyn Read>>,
    child: Option<Child>,
    pub genes: Vec<String>,
}

impl StateReader {
    pub fn open(path: &Path) -> Result<Self> {
        let mut file = File::open(path)?;
        let mut magic = [0; 4];
        file.read_exact(&mut magic)?;
        file.rewind()?;
        let mut child = None;
        let input: Box<dyn Read> = if magic == [0x28, 0xb5, 0x2f, 0xfd] {
            let mut process = Command::new("zstd")
                .args(["-q", "-d", "-c", "--"])
                .arg(path)
                .stdin(Stdio::null())
                .stdout(Stdio::piped())
                .spawn()
                .context("failed to spawn zstd for CSTRAD2 input")?;
            let stdout = process.stdout.take().context("missing zstd stdout")?;
            child = Some(process);
            Box::new(stdout)
        } else if magic[..2] == [0x1f, 0x8b] {
            Box::new(MultiGzDecoder::new(file))
        } else {
            Box::new(file)
        };
        let mut result = Self {
            reader: BufReader::with_capacity(1 << 20, input),
            child,
            genes: Vec::new(),
        };
        let mut header = [0; 28];
        result.reader.read_exact(&mut header)?;
        ensure!(
            &header[..8] == MAGIC,
            "count-splicing requires CSTRAD2 from map --splicing-state (not CSTRAD1 or TSV)"
        );
        for (offset, expected) in [(8, 2), (12, 28), (16, 1), (20, 5)] {
            ensure!(
                u32_at(&header, offset) == expected,
                "unsupported CSTRAD2 header"
            );
        }
        let n = u32_at(&header, 24);
        for _ in 0..n {
            let mut size = [0; 2];
            result.reader.read_exact(&mut size)?;
            let mut bytes = vec![0; u16::from_le_bytes(size) as usize];
            result.reader.read_exact(&mut bytes)?;
            let name = String::from_utf8(bytes)?;
            ensure!(
                !name.contains(['\t', '\n', '\r']),
                "invalid gene name in CSTRAD2"
            );
            result.genes.push(name);
        }
        Ok(result)
    }

    pub fn next_into(&mut self, record: &mut StateRecord) -> Result<bool> {
        let mut fixed = [0; 28];
        if self.reader.read(&mut fixed[..1])? == 0 {
            if let Some(mut child) = self.child.take() {
                ensure!(child.wait()?.success(), "zstd decompression failed");
            }
            return Ok(false);
        }
        self.reader
            .read_exact(&mut fixed[1..])
            .context("truncated CSTRAD2 record")?;
        record.barcode_umi = u64::from_le_bytes(fixed[..8].try_into().unwrap());
        record.valid_mask = u32_at(&fixed, 8);
        record.assignment_type = fixed[24];
        ensure!(
            (1..=7).contains(&record.assignment_type) && fixed[25] == 0,
            "invalid CSTRAD2 assignment type/reserved byte"
        );
        let n = u16::from_le_bytes(fixed[26..28].try_into().unwrap()) as usize;
        ensure!(
            n <= self.genes.len(),
            "CSTRAD2 candidate count exceeds gene dictionary"
        );
        record.candidates.clear();
        for _ in 0..n {
            let mut bytes = [0; 5];
            self.reader
                .read_exact(&mut bytes)
                .context("truncated CSTRAD2 gene/state")?;
            let gene = u32_at(&bytes, 0);
            ensure!(
                (gene as usize) < self.genes.len(),
                "CSTRAD2 gene ID out of bounds"
            );
            ensure!(
                record
                    .candidates
                    .last()
                    .is_none_or(|&(previous, _)| previous < gene),
                "CSTRAD2 genes must be sorted and unique"
            );
            let state = SplicingState::from_code(bytes[4])
                .ok_or_else(|| anyhow!("invalid CSTRAD2 splicing state"))?;
            record.candidates.push((gene, state));
        }
        let primary = u32_at(&fixed, 12);
        match record.assignment_type {
            1 | 3 => ensure!(
                n == 1 && primary == record.candidates[0].0,
                "inconsistent unique-gene record"
            ),
            2 => ensure!(
                n > 1 && primary == u32::MAX,
                "inconsistent ambiguous-gene record"
            ),
            5..=7 => ensure!(
                n == 0 && primary == u32::MAX,
                "filtered record contains gene candidates"
            ),
            _ => bail!("unsupported assignment type for forward CSTRAD2"),
        }
        Ok(true)
    }
}

impl Drop for StateReader {
    fn drop(&mut self) {
        if let Some(mut child) = self.child.take() {
            let _ = child.kill();
            let _ = child.wait();
        }
    }
}

fn u32_at(bytes: &[u8], offset: usize) -> u32 {
    u32::from_le_bytes(bytes[offset..offset + 4].try_into().unwrap())
}

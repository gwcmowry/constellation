use std::collections::HashMap;
use std::fs::File;
use std::io::{self, BufRead, BufReader, BufWriter, Write};
use std::path::Path;

use crate::dna::encode_kmer;
use crate::graph::{Component, TranscriptPath};

#[derive(Clone, Debug, Eq, PartialEq)]
pub struct SeedPosting {
    pub key: u64,
    pub component_id: u32,
    pub pathset: u64,
}

#[derive(Clone, Debug)]
pub struct SeedIndex {
    pub k: usize,
    postings: Vec<SeedPosting>,
}

impl SeedIndex {
    pub fn build(components: &[Component], k: usize) -> Self {
        let mut by_seed_component: HashMap<(u64, u32), u64> = HashMap::new();
        for component in components {
            for (tx_idx, transcript) in component.transcripts.iter().enumerate() {
                if transcript.sequence.len() < k {
                    continue;
                }
                let bit = 1u64 << tx_idx;
                for window in transcript.sequence.windows(k) {
                    if let Some(key) = encode_kmer(window) {
                        *by_seed_component.entry((key, component.id)).or_insert(0) |= bit;
                    }
                }
            }
        }

        let mut postings: Vec<_> = by_seed_component
            .into_iter()
            .map(|((key, component_id), pathset)| SeedPosting {
                key,
                component_id,
                pathset,
            })
            .collect();
        postings.sort_by_key(|posting| (posting.key, posting.component_id));
        Self { k, postings }
    }

    pub fn lookup(&self, key: u64) -> &[SeedPosting] {
        let start = self.postings.partition_point(|posting| posting.key < key);
        let end = self.postings[start..].partition_point(|posting| posting.key == key) + start;
        &self.postings[start..end]
    }

    pub fn len(&self) -> usize {
        self.postings.len()
    }
}

#[derive(Clone, Debug, Eq, PartialEq)]
pub struct Certificate {
    pub id: u32,
    pub component_id: u32,
    pub pathset: u64,
    pub flags: u32,
    pub count: u64,
}

#[derive(Clone, Debug, Default)]
pub struct CertificateStore {
    by_key: HashMap<(u32, u64, u32), u32>,
    certificates: Vec<Certificate>,
}

impl CertificateStore {
    pub fn intern(&mut self, component_id: u32, pathset: u64, flags: u32) -> u32 {
        if let Some(&id) = self.by_key.get(&(component_id, pathset, flags)) {
            if let Some(cert) = self.certificates.get_mut(id as usize) {
                cert.count += 1;
            }
            return id;
        }
        let id = self.certificates.len() as u32;
        self.by_key.insert((component_id, pathset, flags), id);
        self.certificates.push(Certificate {
            id,
            component_id,
            pathset,
            flags,
            count: 1,
        });
        id
    }

    pub fn certificates(&self) -> &[Certificate] {
        &self.certificates
    }
}

#[derive(Clone, Debug)]
pub struct PathCertIndex {
    pub k: usize,
    pub components: Vec<Component>,
    pub seed_index: SeedIndex,
}

impl PathCertIndex {
    pub fn new(k: usize, mut components: Vec<Component>) -> io::Result<Self> {
        if k == 0 || k > 31 {
            return Err(io::Error::new(
                io::ErrorKind::InvalidInput,
                "k must be between 1 and 31",
            ));
        }
        for (idx, component) in components.iter_mut().enumerate() {
            component.id = idx as u32;
            if component.transcripts.len() > 64 {
                return Err(io::Error::new(
                    io::ErrorKind::InvalidData,
                    format!(
                        "component {} has more than 64 transcripts",
                        component.gene_id
                    ),
                ));
            }
        }
        let seed_index = SeedIndex::build(&components, k);
        Ok(Self {
            k,
            components,
            seed_index,
        })
    }

    pub fn component(&self, component_id: u32) -> Option<&Component> {
        self.components.get(component_id as usize)
    }

    pub fn write_text(&self, path: impl AsRef<Path>) -> io::Result<()> {
        let mut writer = BufWriter::new(File::create(path)?);
        writeln!(writer, "PATHCERT_TEXT_INDEX\t1")?;
        writeln!(writer, "K\t{}", self.k)?;
        for component in &self.components {
            writeln!(
                writer,
                "COMP\t{}\t{}\t{}",
                component.id,
                escape(&component.gene_id),
                component.transcripts.len()
            )?;
            for transcript in &component.transcripts {
                writeln!(
                    writer,
                    "TX\t{}\t{}",
                    escape(&transcript.name),
                    String::from_utf8_lossy(&transcript.sequence)
                )?;
            }
            writeln!(writer, "ENDCOMP")?;
        }
        Ok(())
    }

    pub fn read_text(path: impl AsRef<Path>) -> io::Result<Self> {
        let file = File::open(path)?;
        let reader = BufReader::new(file);
        let mut lines = reader.lines();

        let magic = next_line(&mut lines, "magic")?;
        if magic != "PATHCERT_TEXT_INDEX\t1" {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                "not a PathCert text index",
            ));
        }
        let k_line = next_line(&mut lines, "k")?;
        let k = k_line
            .strip_prefix("K\t")
            .ok_or_else(|| invalid_data("missing K line"))?
            .parse::<usize>()
            .map_err(invalid_data)?;

        let mut components = Vec::new();
        let mut current: Option<Component> = None;
        for line in lines {
            let line = line?;
            if line == "ENDCOMP" {
                if let Some(component) = current.take() {
                    components.push(component);
                }
                continue;
            }
            let fields: Vec<&str> = line.split('\t').collect();
            match fields.as_slice() {
                ["COMP", id, gene_id, _count] => {
                    if current.is_some() {
                        return Err(invalid_data("nested component in index"));
                    }
                    current = Some(Component {
                        id: id.parse::<u32>().map_err(invalid_data)?,
                        gene_id: unescape(gene_id),
                        transcripts: Vec::new(),
                    });
                }
                ["TX", name, seq] => {
                    let Some(component) = current.as_mut() else {
                        return Err(invalid_data("transcript outside component"));
                    };
                    component.transcripts.push(TranscriptPath {
                        name: unescape(name),
                        sequence: seq.as_bytes().to_vec(),
                    });
                }
                _ => return Err(invalid_data(format!("bad index line: {line}"))),
            }
        }

        if current.is_some() {
            return Err(invalid_data("unterminated component"));
        }

        Self::new(k, components)
    }
}

pub fn intersect_pathsets(a: u64, b: u64) -> u64 {
    a & b
}

pub fn pathset_size(pathset: u64) -> u32 {
    pathset.count_ones()
}

fn next_line(
    lines: &mut impl Iterator<Item = io::Result<String>>,
    label: &str,
) -> io::Result<String> {
    lines.next().ok_or_else(|| {
        io::Error::new(
            io::ErrorKind::UnexpectedEof,
            format!("missing index {label} line"),
        )
    })?
}

fn invalid_data(error: impl ToString) -> io::Error {
    io::Error::new(io::ErrorKind::InvalidData, error.to_string())
}

fn escape(text: &str) -> String {
    text.replace('\\', "\\\\").replace('\t', "\\t")
}

fn unescape(text: &str) -> String {
    let mut output = String::new();
    let mut chars = text.chars();
    while let Some(ch) = chars.next() {
        if ch == '\\' {
            match chars.next() {
                Some('t') => output.push('\t'),
                Some('\\') => output.push('\\'),
                Some(other) => {
                    output.push('\\');
                    output.push(other);
                }
                None => output.push('\\'),
            }
        } else {
            output.push(ch);
        }
    }
    output
}

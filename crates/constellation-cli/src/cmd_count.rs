use crate::CountArgs;
use anyhow::{anyhow, Context, Result};
use constellation_core::index::{IndexAccess, LoadedIndex};
use flate2::read::MultiGzDecoder;
use rustc_hash::{FxHashMap, FxHashSet};
use serde_json::json;
use std::cmp::Reverse;
use std::collections::{BTreeMap, BTreeSet};
use std::fs::{self, File};
use std::io::{BufRead, BufReader, Read};
use std::path::{Path, PathBuf};
use std::time::Instant;

const COUNTABLE_ASSIGNMENTS: [&str; 2] = ["unique_gene", "ambiguous_transcript_same_gene"];
const CELL_BARCODE_LEN: u8 = 16;
const UMI_INVALID_FLAG: u64 = 1_u64 << 57;
const UMI_CODE_MASK: u64 = UMI_INVALID_FLAG - 1;

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
struct MoleculeKey {
    cell_id: u32,
    umi: u64,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
struct MoleculeGeneKey {
    cell_id: u32,
    umi: u64,
    gene_id: u32,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, PartialOrd, Ord)]
struct CellGeneKey {
    cell_id: u32,
    gene_id: u32,
}

#[derive(Debug, Default)]
struct CountMetrics {
    rows: u64,
    countable_rows: u64,
    rows_rejected_by_barcode: u64,
    rows_rejected_by_feature: u64,
    barcode_exact: u64,
    barcode_corrected: u64,
    raw_molecule_gene_pairs: u64,
    accepted_molecules: u64,
    ambiguous_molecules_skipped: u64,
    raw_gene_umi_pairs: u64,
    corrected_gene_umi_pairs: u64,
}

#[derive(Debug)]
struct BarcodeCorrector {
    whitelist: FxHashSet<u64>,
}

#[derive(Debug, Default)]
struct CellInterner {
    ids_by_seq: FxHashMap<String, u32>,
    ids_by_code: FxHashMap<u64, u32>,
    seqs: Vec<String>,
}

#[derive(Debug, Clone, Copy)]
enum BarcodeResolution {
    Exact(u32),
    Corrected(u32),
    Unchecked(u32),
}

impl CellInterner {
    fn get_or_insert_seq(&mut self, seq: &str) -> Result<u32> {
        if let Some(&id) = self.ids_by_seq.get(seq) {
            return Ok(id);
        }
        let id = u32::try_from(self.seqs.len()).context("too many cell barcodes")?;
        self.ids_by_seq.insert(seq.to_owned(), id);
        self.seqs.push(seq.to_owned());
        Ok(id)
    }

    fn get_or_insert_code(&mut self, code: u64) -> Result<u32> {
        if let Some(&id) = self.ids_by_code.get(&code) {
            return Ok(id);
        }
        let seq = decode_acgt(code, CELL_BARCODE_LEN);
        let id = self.get_or_insert_seq(&seq)?;
        self.ids_by_code.insert(code, id);
        Ok(id)
    }
}

pub fn run_count(args: CountArgs) -> Result<()> {
    let started = Instant::now();
    let index = LoadedIndex::load(&args.index)?;
    let barcode_corrector = args
        .barcode_whitelist
        .as_ref()
        .map(|path| load_barcode_whitelist(path))
        .transpose()?;
    let feature_filter = args
        .feature_whitelist
        .as_ref()
        .map(|path| load_feature_filter(path, &index))
        .transpose()?;

    let file = File::open(&args.assignments)
        .with_context(|| format!("failed to open {}", args.assignments.display()))?;
    let reader = BufReader::new(file);
    let mut lines = reader.lines();
    let Some(header) = lines.next() else {
        return Err(anyhow!("assignment file is empty"));
    };
    validate_assignment_header(&header?)?;

    let mut metrics = CountMetrics::default();
    let mut cells = CellInterner::default();
    let mut molecule_gene_counts: FxHashMap<MoleculeGeneKey, u32> = FxHashMap::default();
    let mut barcode_cache: FxHashMap<String, Option<BarcodeResolution>> = FxHashMap::default();

    for line in lines {
        let line = line?;
        if line.trim().is_empty() {
            continue;
        }
        metrics.rows += 1;
        let fields: Vec<_> = line.split('\t').collect();
        if fields.len() < 9 {
            return Err(anyhow!(
                "assignment line {} has wrong field count",
                metrics.rows + 1
            ));
        }
        if !COUNTABLE_ASSIGNMENTS.contains(&fields[3]) || fields[4] == "." {
            continue;
        }
        metrics.countable_rows += 1;
        let gene_id: u32 = fields[4].parse()?;
        if feature_filter
            .as_ref()
            .is_some_and(|filter| !filter.contains(&gene_id))
        {
            metrics.rows_rejected_by_feature += 1;
            continue;
        }
        let Some(cell_id) = corrected_cell_id(
            fields[1],
            &barcode_corrector,
            args.correct_barcodes,
            &mut cells,
            &mut barcode_cache,
            &mut metrics,
        )?
        else {
            metrics.rows_rejected_by_barcode += 1;
            continue;
        };
        let Some(umi) = pack_umi(fields[2]) else {
            continue;
        };
        *molecule_gene_counts
            .entry(MoleculeGeneKey {
                cell_id,
                umi,
                gene_id,
            })
            .or_default() += 1;
    }

    metrics.raw_molecule_gene_pairs = molecule_gene_counts.len() as u64;
    let gene_umi_counts = resolve_molecules(&molecule_gene_counts, &args, &mut metrics);
    let counts = collapse_umis(&gene_umi_counts, args.correct_umis, args.umi_edit_distance);
    metrics.corrected_gene_umi_pairs = counts.values().map(|&count| count as u64).sum();

    write_outputs(&args, &index, &cells.seqs, &counts)?;
    if let Some(path) = args.emit_metrics.as_ref() {
        fs::write(
            path,
            serde_json::to_string_pretty(&json!({
                "rows": metrics.rows,
                "countable_rows": metrics.countable_rows,
                "rows_rejected_by_barcode": metrics.rows_rejected_by_barcode,
                "rows_rejected_by_feature": metrics.rows_rejected_by_feature,
                "barcode_exact": metrics.barcode_exact,
                "barcode_corrected": metrics.barcode_corrected,
                "raw_molecule_gene_pairs": metrics.raw_molecule_gene_pairs,
                "accepted_molecules": metrics.accepted_molecules,
                "ambiguous_molecules_skipped": metrics.ambiguous_molecules_skipped,
                "raw_gene_umi_pairs": metrics.raw_gene_umi_pairs,
                "corrected_gene_umi_pairs": metrics.corrected_gene_umi_pairs,
                "cells": cells.seqs.len(),
                "nonzero_gene_cell_entries": counts.len(),
                "count_seconds": started.elapsed().as_secs_f64(),
                "barcode_whitelist": args.barcode_whitelist.as_ref().map(|path| path.display().to_string()),
                "feature_whitelist": args.feature_whitelist.as_ref().map(|path| path.display().to_string()),
                "correct_barcodes": args.correct_barcodes,
                "correct_umis": args.correct_umis,
                "umi_edit_distance": args.umi_edit_distance,
                "resolve_molecule_genes": args.resolve_molecule_genes,
                "molecule_gene_ratio": args.molecule_gene_ratio,
            }))? + "\n",
        )?;
    }
    Ok(())
}

fn validate_assignment_header(header: &str) -> Result<()> {
    let fields: Vec<_> = header.split('\t').collect();
    let expected = [
        "read_id",
        "cell_barcode",
        "umi",
        "assignment_type",
        "gene_id",
        "transcript_id",
        "candidate_count",
        "score",
        "flags",
    ];
    if fields.len() < expected.len() || fields[..expected.len()] != expected {
        return Err(anyhow!("unexpected assignment TSV header"));
    }
    Ok(())
}

fn corrected_cell_id(
    barcode: &str,
    corrector: &Option<BarcodeCorrector>,
    correct_barcodes: bool,
    cells: &mut CellInterner,
    cache: &mut FxHashMap<String, Option<BarcodeResolution>>,
    metrics: &mut CountMetrics,
) -> Result<Option<u32>> {
    if let Some(cached) = cache.get(barcode) {
        return Ok(record_barcode_resolution(*cached, metrics));
    }
    let resolved = if let Some(corrector) = corrector {
        let Some(code) = pack_acgt(barcode) else {
            cache.insert(barcode.to_owned(), None);
            return Ok(None);
        };
        if corrector.whitelist.contains(&code) {
            Some(BarcodeResolution::Exact(cells.get_or_insert_code(code)?))
        } else if correct_barcodes {
            let corrected = unique_one_mismatch_whitelist_neighbor(code, &corrector.whitelist);
            if let Some(corrected) = corrected {
                Some(BarcodeResolution::Corrected(
                    cells.get_or_insert_code(corrected)?,
                ))
            } else {
                None
            }
        } else {
            None
        }
    } else {
        Some(BarcodeResolution::Unchecked(
            cells.get_or_insert_seq(barcode)?,
        ))
    };
    cache.insert(barcode.to_owned(), resolved);
    Ok(record_barcode_resolution(resolved, metrics))
}

fn record_barcode_resolution(
    resolution: Option<BarcodeResolution>,
    metrics: &mut CountMetrics,
) -> Option<u32> {
    match resolution {
        Some(BarcodeResolution::Exact(id)) => {
            metrics.barcode_exact += 1;
            Some(id)
        }
        Some(BarcodeResolution::Corrected(id)) => {
            metrics.barcode_corrected += 1;
            Some(id)
        }
        Some(BarcodeResolution::Unchecked(id)) => Some(id),
        None => None,
    }
}

fn resolve_molecules(
    molecule_gene_counts: &FxHashMap<MoleculeGeneKey, u32>,
    args: &CountArgs,
    metrics: &mut CountMetrics,
) -> FxHashMap<CellGeneKey, Vec<(u64, u32)>> {
    let mut entries: Vec<_> = molecule_gene_counts
        .iter()
        .map(|(key, &count)| (*key, count))
        .collect();
    entries.sort_unstable_by_key(|(key, count)| (key.cell_id, key.umi, Reverse(*count)));

    let mut out: FxHashMap<CellGeneKey, Vec<(u64, u32)>> = FxHashMap::default();
    let mut idx = 0;
    while idx < entries.len() {
        let molecule = MoleculeKey {
            cell_id: entries[idx].0.cell_id,
            umi: entries[idx].0.umi,
        };
        let start = idx;
        while idx < entries.len()
            && entries[idx].0.cell_id == molecule.cell_id
            && entries[idx].0.umi == molecule.umi
        {
            idx += 1;
        }
        let group = &entries[start..idx];
        let selected = if !args.resolve_molecule_genes || group.len() == 1 {
            None
        } else {
            dominant_gene(group, args.molecule_gene_ratio.max(1))
        };
        if let Some((gene_id, count)) = selected {
            out.entry(CellGeneKey {
                cell_id: molecule.cell_id,
                gene_id,
            })
            .or_default()
            .push((molecule.umi, count));
            metrics.accepted_molecules += 1;
        } else if group.len() == 1 || !args.resolve_molecule_genes {
            for (key, count) in group {
                out.entry(CellGeneKey {
                    cell_id: key.cell_id,
                    gene_id: key.gene_id,
                })
                .or_default()
                .push((key.umi, *count));
                metrics.accepted_molecules += 1;
            }
        } else {
            metrics.ambiguous_molecules_skipped += 1;
        }
    }
    metrics.raw_gene_umi_pairs = out.values().map(|umis| umis.len() as u64).sum();
    out
}

fn dominant_gene(group: &[(MoleculeGeneKey, u32)], ratio: u32) -> Option<(u32, u32)> {
    let mut ranked: Vec<_> = group
        .iter()
        .map(|(key, count)| (key.gene_id, *count))
        .collect();
    ranked.sort_unstable_by_key(|&(_, count)| Reverse(count));
    let (top_gene, top_count) = ranked[0];
    let second_count = ranked.get(1).map_or(0, |&(_, count)| count);
    if second_count == 0 || top_count >= ratio.saturating_mul(second_count) {
        Some((top_gene, top_count))
    } else {
        None
    }
}

fn collapse_umis(
    gene_umi_counts: &FxHashMap<CellGeneKey, Vec<(u64, u32)>>,
    correct_umis: bool,
    umi_edit_distance: u8,
) -> BTreeMap<(u32, usize), u32> {
    let mut counts = BTreeMap::new();
    for (cell_gene, umis) in gene_umi_counts {
        let molecule_count = if correct_umis && umi_edit_distance > 0 {
            directional_umi_count(umis, umi_edit_distance)
        } else {
            umis.len() as u32
        };
        if molecule_count > 0 {
            counts.insert(
                (cell_gene.gene_id, cell_gene.cell_id as usize),
                molecule_count,
            );
        }
    }
    counts
}

fn directional_umi_count(umis: &[(u64, u32)], max_edit_distance: u8) -> u32 {
    if max_edit_distance == 0 || umis.len() <= 1 {
        return umis.len() as u32;
    }
    let counts: FxHashMap<_, _> = umis.iter().copied().collect();
    let mut ordered = umis.to_vec();
    ordered.sort_unstable_by_key(|&(umi, count)| (Reverse(count), umi));
    let mut assigned = FxHashSet::default();
    let mut roots = 0_u32;
    for (umi, count) in ordered {
        if assigned.contains(&umi) {
            continue;
        }
        roots += 1;
        assigned.insert(umi);
        if max_edit_distance == 1 {
            for neighbor in one_mismatch_umis(umi) {
                if assigned.contains(&neighbor) {
                    continue;
                }
                if let Some(&neighbor_count) = counts.get(&neighbor) {
                    if count >= neighbor_count.saturating_mul(2).saturating_sub(1) {
                        assigned.insert(neighbor);
                    }
                }
            }
        }
    }
    roots
}

fn write_outputs(
    args: &CountArgs,
    index: &LoadedIndex,
    cell_list: &[String],
    counts: &BTreeMap<(u32, usize), u32>,
) -> Result<()> {
    fs::write(
        with_suffix(&args.out_prefix, "_barcodes.tsv"),
        cell_list.join("\n") + "\n",
    )?;
    let mut features = String::new();
    for idx in 0..index.num_genes() {
        let gene = index.gene_name(idx as u32).unwrap_or(".");
        features.push_str(&format!("{idx}\t{gene}\tGene Expression\n"));
    }
    fs::write(with_suffix(&args.out_prefix, "_features.tsv"), features)?;

    let mut matrix = String::from("%%MatrixMarket matrix coordinate integer general\n");
    matrix.push_str(&format!(
        "{} {} {}\n",
        index.num_genes(),
        cell_list.len(),
        counts.len()
    ));
    for ((gene_id, cell_idx), count) in counts {
        matrix.push_str(&format!("{} {} {}\n", gene_id + 1, cell_idx + 1, count));
    }
    fs::write(with_suffix(&args.out_prefix, "_matrix.mtx"), matrix)?;
    Ok(())
}

fn load_barcode_whitelist(path: &Path) -> Result<BarcodeCorrector> {
    let mut whitelist = FxHashSet::default();
    for line in read_lines_maybe_gz(path)? {
        let barcode = line?;
        let barcode = barcode
            .split_once('-')
            .map_or(barcode.as_str(), |(base, _)| base)
            .trim();
        if barcode.is_empty() {
            continue;
        }
        let Some(code) = pack_acgt(barcode) else {
            continue;
        };
        whitelist.insert(code);
    }
    if whitelist.is_empty() {
        return Err(anyhow!("barcode whitelist is empty: {}", path.display()));
    }
    Ok(BarcodeCorrector { whitelist })
}

fn load_feature_filter(path: &Path, index: &LoadedIndex) -> Result<FxHashSet<u32>> {
    let mut feature_names = BTreeSet::new();
    for line in read_lines_maybe_gz(path)? {
        let line = line?;
        let fields: Vec<_> = line.split('\t').collect();
        if let Some(gene_id) = fields.first() {
            feature_names.insert((*gene_id).to_owned());
        }
        if let Some(gene_name) = fields.get(1) {
            feature_names.insert((*gene_name).to_owned());
        }
    }
    let mut filter = FxHashSet::default();
    for idx in 0..index.num_genes() {
        if let Some(gene) = index.gene_name(idx as u32) {
            if feature_names.contains(gene) {
                filter.insert(idx as u32);
            }
        }
    }
    if filter.is_empty() {
        return Err(anyhow!(
            "feature whitelist did not match any index genes: {}",
            path.display()
        ));
    }
    Ok(filter)
}

fn read_lines_maybe_gz(path: &Path) -> Result<Box<dyn Iterator<Item = std::io::Result<String>>>> {
    let file = File::open(path).with_context(|| format!("failed to open {}", path.display()))?;
    let reader: Box<dyn Read> = if path.extension().is_some_and(|ext| ext == "gz") {
        Box::new(MultiGzDecoder::new(file))
    } else {
        Box::new(file)
    };
    Ok(Box::new(BufReader::new(reader).lines()))
}

fn unique_one_mismatch_whitelist_neighbor(code: u64, whitelist: &FxHashSet<u64>) -> Option<u64> {
    let mut found = None;
    for neighbor in one_mismatch_codes(code, CELL_BARCODE_LEN) {
        if !whitelist.contains(&neighbor) {
            continue;
        }
        if found.replace(neighbor).is_some() {
            return None;
        }
    }
    found
}

fn one_mismatch_umis(umi: u64) -> Vec<u64> {
    if umi & UMI_INVALID_FLAG != 0 {
        return Vec::new();
    }
    let len = (umi >> 58) as u8;
    let code = umi & UMI_CODE_MASK;
    one_mismatch_codes(code, len)
        .into_iter()
        .map(|code| pack_with_len(code, len))
        .collect()
}

fn one_mismatch_codes(code: u64, len: u8) -> Vec<u64> {
    let mut out = Vec::with_capacity(len as usize * 3);
    for pos in 0..len {
        let shift = 2 * (len - pos - 1);
        let base = (code >> shift) & 0b11;
        for alt in 0..4 {
            if alt != base {
                out.push((code & !(0b11 << shift)) | (alt << shift));
            }
        }
    }
    out
}

fn pack_umi(seq: &str) -> Option<u64> {
    let len = u8::try_from(seq.len()).ok()?;
    if len == 0 || len > 28 {
        return None;
    }
    Some(match pack_acgt(seq) {
        Some(code) => pack_with_len(code, len),
        None => ((len as u64) << 58) | UMI_INVALID_FLAG | stable_umi_hash(seq),
    })
}

fn pack_with_len(code: u64, len: u8) -> u64 {
    ((len as u64) << 58) | code
}

fn stable_umi_hash(seq: &str) -> u64 {
    let mut hash = 0xcbf29ce484222325_u64;
    for byte in seq.bytes() {
        hash ^= byte as u64;
        hash = hash.wrapping_mul(0x100000001b3);
    }
    hash & UMI_CODE_MASK
}

fn pack_acgt(seq: &str) -> Option<u64> {
    let mut out = 0_u64;
    for byte in seq.bytes() {
        out = (out << 2)
            | match byte {
                b'A' | b'a' => 0,
                b'C' | b'c' => 1,
                b'G' | b'g' => 2,
                b'T' | b't' => 3,
                _ => return None,
            };
    }
    Some(out)
}

fn decode_acgt(code: u64, len: u8) -> String {
    let mut out = String::with_capacity(len as usize);
    for pos in 0..len {
        let shift = 2 * (len - pos - 1);
        out.push(match (code >> shift) & 0b11 {
            0 => 'A',
            1 => 'C',
            2 => 'G',
            _ => 'T',
        });
    }
    out
}

fn with_suffix(prefix: &Path, suffix: &str) -> PathBuf {
    PathBuf::from(format!("{}{}", prefix.display(), suffix))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn directional_umi_correction_collapses_one_mismatch_low_count_neighbor() {
        let parent = pack_umi("AAAAAAAAAAAA").unwrap();
        let child = pack_umi("AAAAAAAAAAAT").unwrap();
        let far = pack_umi("CCCCCCCCCCCC").unwrap();
        assert_eq!(
            directional_umi_count(&[(parent, 4), (child, 1), (far, 1)], 1),
            2
        );
    }

    #[test]
    fn barcode_correction_requires_unique_neighbor() {
        let a = pack_acgt("AAAAAAAAAAAAAAAA").unwrap();
        let c = pack_acgt("CAAAAAAAAAAAAAAA").unwrap();
        let query = pack_acgt("TAAAAAAAAAAAAAAA").unwrap();
        let mut whitelist = FxHashSet::default();
        whitelist.insert(a);
        assert_eq!(
            unique_one_mismatch_whitelist_neighbor(query, &whitelist),
            Some(a)
        );
        whitelist.insert(c);
        assert_eq!(
            unique_one_mismatch_whitelist_neighbor(query, &whitelist),
            None
        );
    }
}

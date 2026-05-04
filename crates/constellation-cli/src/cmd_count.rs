use crate::CountArgs;
use anyhow::{anyhow, Result};
use constellation_core::index::{IndexAccess, LoadedIndex};
use std::collections::{BTreeMap, BTreeSet};
use std::fs;
use std::path::{Path, PathBuf};

pub fn run_count(args: CountArgs) -> Result<()> {
    let index = LoadedIndex::load(&args.index)?;
    let text = fs::read_to_string(args.assignments)?;
    let mut cells = BTreeSet::new();
    let mut dedup = BTreeSet::new();

    for (line_idx, line) in text.lines().enumerate().skip(1) {
        if line.trim().is_empty() {
            continue;
        }
        let fields: Vec<_> = line.split('\t').collect();
        if fields.len() != 9 {
            return Err(anyhow!(
                "assignment line {} has wrong field count",
                line_idx + 1
            ));
        }
        if !matches!(fields[3], "unique_gene" | "ambiguous_transcript_same_gene")
            || fields[4] == "."
        {
            continue;
        }
        let gene_id: u32 = fields[4].parse()?;
        let cb = fields[1].to_owned();
        let umi = fields[2].to_owned();
        cells.insert(cb.clone());
        dedup.insert((gene_id, cb, umi));
    }

    let cell_list: Vec<_> = cells.into_iter().collect();
    let cell_ids: BTreeMap<_, _> = cell_list
        .iter()
        .enumerate()
        .map(|(idx, cell)| (cell.as_str(), idx))
        .collect();
    let mut counts: BTreeMap<(u32, usize), u32> = BTreeMap::new();
    for (gene_id, cell, _umi) in dedup {
        let Some(&cell_idx) = cell_ids.get(cell.as_str()) else {
            continue;
        };
        *counts.entry((gene_id, cell_idx)).or_default() += 1;
    }

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

fn with_suffix(prefix: &Path, suffix: &str) -> PathBuf {
    PathBuf::from(format!("{}{}", prefix.display(), suffix))
}

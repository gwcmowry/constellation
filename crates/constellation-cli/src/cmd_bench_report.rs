use crate::BenchReportArgs;
use anyhow::{anyhow, Result};
use constellation_core::index::{IndexAccess, LoadedIndex};
use constellation_core::metrics::parse_perf_stat;
use std::collections::BTreeMap;
use std::fs;

pub fn run_bench_report(args: BenchReportArgs) -> Result<()> {
    let assignments = fs::read_to_string(args.assignments)?;
    let assignment_rows = parse_assignments(&assignments)?;
    println!("assignments\t{}", assignment_rows.len());

    if let Some(metrics) = args.metrics {
        let metrics = fs::read_to_string(metrics)?;
        println!("metrics_json\t{}", metrics.replace('\n', ""));
    }

    if let Some(perf_stat) = args.perf_stat {
        let perf = parse_perf_stat(&fs::read_to_string(perf_stat)?);
        if let Some(rate) = perf.l1_miss_rate() {
            println!("L1_miss_rate\t{rate:.6}");
        }
        if let Some(rate) = perf.llc_miss_rate() {
            println!("LLC_miss_rate\t{rate:.6}");
        }
    }

    if let Some(truth) = args.truth {
        let truth_rows = parse_truth(&fs::read_to_string(truth)?)?;
        println!("truth\t{}", truth_rows.len());
        let index = args.index.map(LoadedIndex::load).transpose()?;
        report_truth_comparison(&assignment_rows, &truth_rows, index.as_ref());
    }

    Ok(())
}

#[derive(Debug, Clone)]
struct AssignmentRow {
    read_id: u64,
    assignment_type: String,
    gene_id: Option<u32>,
}

#[derive(Debug, Clone)]
struct TruthRow {
    read_id: u64,
    gene: String,
}

fn parse_assignments(text: &str) -> Result<Vec<AssignmentRow>> {
    let mut rows = Vec::new();
    for (line_idx, line) in text.lines().enumerate().skip(1) {
        if line.trim().is_empty() {
            continue;
        }
        let fields: Vec<_> = line.split('\t').collect();
        if fields.len() != 9 {
            return Err(anyhow!(
                "assignment line {} has {} fields, expected 9",
                line_idx + 1,
                fields.len()
            ));
        }
        rows.push(AssignmentRow {
            read_id: fields[0].parse()?,
            assignment_type: fields[3].to_owned(),
            gene_id: parse_optional_u32(fields[4])?,
        });
    }
    Ok(rows)
}

fn parse_truth(text: &str) -> Result<Vec<TruthRow>> {
    let mut rows = Vec::new();
    for (line_idx, line) in text.lines().enumerate().skip(1) {
        if line.trim().is_empty() {
            continue;
        }
        let fields: Vec<_> = line.split('\t').collect();
        if fields.len() < 4 {
            return Err(anyhow!(
                "truth line {} has {} fields, expected at least 4",
                line_idx + 1,
                fields.len()
            ));
        }
        rows.push(TruthRow {
            read_id: fields[0].parse()?,
            gene: fields[2].to_owned(),
        });
    }
    Ok(rows)
}

fn parse_optional_u32(text: &str) -> Result<Option<u32>> {
    if text == "." {
        Ok(None)
    } else {
        Ok(Some(text.parse()?))
    }
}

fn report_truth_comparison(
    assignments: &[AssignmentRow],
    truth: &[TruthRow],
    index: Option<&LoadedIndex>,
) {
    let truth_by_read: BTreeMap<_, _> = truth
        .iter()
        .map(|row| (row.read_id, row.gene.as_str()))
        .collect();
    let mut evaluable = 0_u64;
    let mut correct = 0_u64;
    let mut unique = 0_u64;
    let mut unique_correct = 0_u64;
    let mut missing_truth = 0_u64;

    for assignment in assignments {
        let Some(truth_gene) = truth_by_read.get(&assignment.read_id) else {
            missing_truth += 1;
            continue;
        };
        let Some(assigned_gene) = assignment
            .gene_id
            .and_then(|gene_id| resolve_gene_name(gene_id, index))
        else {
            continue;
        };
        evaluable += 1;
        if assigned_gene == *truth_gene {
            correct += 1;
        }
        if assignment.assignment_type == "unique_gene" {
            unique += 1;
            if assigned_gene == *truth_gene {
                unique_correct += 1;
            }
        }
    }

    println!("truth_missing_assignments\t{missing_truth}");
    println!("evaluable_gene_assignments\t{evaluable}");
    println!(
        "correct_gene_assignment_rate\t{:.6}",
        ratio(correct, evaluable)
    );
    println!("unique_gene_assignments\t{unique}");
    println!(
        "unique_gene_correct_rate\t{:.6}",
        ratio(unique_correct, unique)
    );
    if index.is_none() {
        println!(
            "truth_note\tprovide --index to compare numeric assignment gene IDs against gene names"
        );
    }
}

fn resolve_gene_name(gene_id: u32, index: Option<&LoadedIndex>) -> Option<String> {
    match index {
        Some(index) => index.gene_name(gene_id).map(str::to_owned),
        None => Some(gene_id.to_string()),
    }
}

fn ratio(num: u64, den: u64) -> f64 {
    if den == 0 {
        0.0
    } else {
        num as f64 / den as f64
    }
}

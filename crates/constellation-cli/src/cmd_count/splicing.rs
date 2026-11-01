//! Conservative CB/UMI gene intersection, followed by state-aware UMI correction.
use super::*;
use crate::splicing_io::{StateReader, StateRecord};
use crate::SplicingCountArgs;
use constellation_core::ec_index::SplicingState;

struct Molecule {
    reads: u32,
    genes: Vec<(u32, SplicingState)>,
}

impl Molecule {
    fn add(&mut self, genes: &[(u32, SplicingState)]) {
        self.reads = self.reads.saturating_add(1);
        let mut j = 0;
        self.genes.retain_mut(|(gene, state)| {
            while j < genes.len() && genes[j].0 < *gene {
                j += 1;
            }
            if j < genes.len() && genes[j].0 == *gene {
                *state = state.intersect(genes[j].1);
                true
            } else {
                false
            }
        });
    }
}

#[derive(Default)]
struct CellEvidence {
    reads: Vec<ReadEvidence>,
    candidates: Vec<(u32, SplicingState)>,
}
struct ReadEvidence {
    umi: u64,
    start: u32,
    end: u32,
}
#[derive(Default)]
struct CellCounts {
    counts: Vec<(u32, [u32; 5])>,
    raw: usize,
    unresolved: usize,
    conflicting: usize,
    before: usize,
    rescued: usize,
    changed_genes: usize,
    changed_states: usize,
    spliced_ties: usize,
    fallback: usize,
}

// CSTRAD2 uses the first nucleotide in the low bits; the existing barcode and
// UMI correction routines use it in the high bits. Preserve their tie ordering.
fn reverse_pairs(code: u64, len: u8) -> u64 {
    let reversed = code.reverse_bits();
    (((reversed & 0x5555_5555_5555_5555) << 1) | ((reversed & 0xaaaa_aaaa_aaaa_aaaa) >> 1))
        >> (64 - 2 * len)
}

fn correct_packed_barcode(
    code: u64,
    corrector: &Option<BarcodeCorrector>,
    correct: bool,
    cells: &mut CellInterner,
    cache: &mut FxHashMap<u64, Option<BarcodeResolution>>,
    metrics: &mut CountMetrics,
) -> Result<Option<u32>> {
    let resolution = if let Some(&cached) = cache.get(&code) {
        cached
    } else {
        let resolution = if let Some(corrector) = corrector {
            if corrector.whitelist.contains(&code) {
                Some(BarcodeResolution::Exact(cells.get_or_insert_code(code)?))
            } else if correct {
                unique_one_mismatch_whitelist_neighbor(code, &corrector.whitelist)
                    .map(|c| {
                        cells
                            .get_or_insert_code(c)
                            .map(BarcodeResolution::Corrected)
                    })
                    .transpose()?
            } else {
                None
            }
        } else {
            Some(BarcodeResolution::Unchecked(
                cells.get_or_insert_code(code)?,
            ))
        };
        cache.insert(code, resolution);
        resolution
    };
    Ok(record_barcode_resolution(resolution, metrics))
}

pub fn run(args: SplicingCountArgs) -> Result<()> {
    use crate::SplicingResolutionArg;
    use rayon::prelude::*;
    anyhow::ensure!(
        args.resolution != SplicingResolutionArg::Usa || !args.correct_umis,
        "--resolution usa currently requires --no-correct-umis (exact UMI matching)"
    );
    let started = Instant::now();
    let mut input = StateReader::open(&args.assignments)?;
    let corrector = args
        .barcode_whitelist
        .as_ref()
        .map(|path| load_barcode_whitelist(path))
        .transpose()?;
    let mut cells = CellInterner::default();
    let mut cache = FxHashMap::default();
    let mut metrics = CountMetrics::default();
    let mut evidence: Vec<CellEvidence> = Vec::new();
    let mut record = StateRecord::default();
    let mut invalid_barcodes = 0u64;
    while input.next_into(&mut record)? {
        metrics.rows += 1;
        if record.candidates.is_empty() {
            continue;
        }
        if record.valid_mask != (1 << 28) - 1 {
            invalid_barcodes += 1;
            continue;
        }
        let barcode = reverse_pairs(record.barcode_umi & 0xffff_ffff, 16);
        let Some(cell_id) = correct_packed_barcode(
            barcode,
            &corrector,
            args.correct_barcodes,
            &mut cells,
            &mut cache,
            &mut metrics,
        )?
        else {
            metrics.rows_rejected_by_barcode += 1;
            continue;
        };
        let umi = pack_with_len(
            reverse_pairs((record.barcode_umi >> 32) & 0xff_ffff, 12),
            12,
        );
        metrics.countable_rows += 1;
        if evidence.len() <= cell_id as usize {
            evidence.resize_with(cell_id as usize + 1, CellEvidence::default);
        }
        let cell = &mut evidence[cell_id as usize];
        let start = u32::try_from(cell.candidates.len())
            .context("too many candidate entries in one cell")?;
        cell.candidates.extend_from_slice(&record.candidates);
        let end = u32::try_from(cell.candidates.len())
            .context("too many candidate entries in one cell")?;
        cell.reads.push(ReadEvidence { umi, start, end });
    }
    let ingest_seconds = started.elapsed().as_secs_f64();
    let phase = Instant::now();
    let resolved: Vec<_> = evidence
        .into_par_iter()
        .map(|cell| resolve_cell(cell, args.resolution, args.correct_umis))
        .collect();
    let resolve_seconds = phase.elapsed().as_secs_f64();
    let phase = Instant::now();
    let mut totals = CellCounts::default();
    let mut state_totals = [0u64; 5];
    let mut counts = Vec::new();
    for (cell, result) in resolved.into_iter().enumerate() {
        totals.raw += result.raw;
        totals.unresolved += result.unresolved;
        totals.conflicting += result.conflicting;
        totals.before += result.before;
        totals.rescued += result.rescued;
        totals.changed_genes += result.changed_genes;
        totals.changed_states += result.changed_states;
        totals.spliced_ties += result.spliced_ties;
        totals.fallback += result.fallback;
        for (gene, values) in result.counts {
            for (total, count) in state_totals.iter_mut().zip(values) {
                *total += count as u64;
            }
            counts.push((gene, cell as u32, values));
        }
    }
    counts.par_sort_unstable_by_key(|&(gene, cell, _)| (gene, cell));
    let gather_seconds = phase.elapsed().as_secs_f64();
    let phase = Instant::now();
    fs::write(
        with_suffix(&args.out_prefix, "_barcodes.tsv"),
        if cells.seqs.is_empty() {
            String::new()
        } else {
            cells.seqs.join("\n") + "\n"
        },
    )?;
    let mut features = String::new();
    for (id, name) in input.genes.iter().enumerate() {
        features.push_str(&format!("{id}\t{name}\tGene Expression\n"));
    }
    fs::write(with_suffix(&args.out_prefix, "_features.tsv"), features)?;
    for (label, from, to) in [
        ("spliced", 0, 1),
        ("unspliced", 1, 2),
        ("compatible_both", 2, 3),
        ("conflicting", 3, 4),
        ("unknown", 4, 5),
        ("ambiguous", 2, 5),
        ("matrix", 0, 5),
    ] {
        write_cell_matrix(
            &args.out_prefix,
            label,
            input.genes.len(),
            cells.seqs.len(),
            &counts,
            from..to,
        )?;
    }
    if let Some(path) = args.emit_metrics {
        fs::write(
            path,
            serde_json::to_vec_pretty(&json!({
                "rows":metrics.rows,"mapped_valid_rows":metrics.countable_rows,
                "invalid_barcode_rows":invalid_barcodes,"barcode_rejected_rows":metrics.rows_rejected_by_barcode,
                "barcode_exact":metrics.barcode_exact,"barcode_corrected":metrics.barcode_corrected,
                "raw_barcode_umi_groups":totals.raw,"unresolved_gene_groups":totals.unresolved,
                "conflicting_gene_groups":totals.conflicting,"molecules_before_umi_correction":totals.before,
                "molecules_after_umi_correction":state_totals.iter().sum::<u64>(),
                "spliced":state_totals[0],"unspliced":state_totals[1],"compatible_both":state_totals[2],
                "conflicting_state":state_totals[3],"unknown_state":state_totals[4],
                "ambiguous_including_conflict_unknown":state_totals[2..].iter().sum::<u64>(),
                "cells":cells.seqs.len(),"correct_umis":args.correct_umis,
                "resolution":args.resolution.label(),
                "policy_rescued_groups":totals.rescued,"policy_changed_genes":totals.changed_genes,
                "policy_changed_states":totals.changed_states,"usa_spliced_tiebreak_groups":totals.spliced_ties,
                "usa_unknown_conflict_fallback_groups":totals.fallback,
                "policy":if args.resolution == SplicingResolutionArg::Intersection {
                    "intersection of candidate genes within corrected barcode/exact UMI; skip unresolved/conflicting genes; deterministic one-hop directional UMI correction within cell/gene; intersect states across merged reads/UMIs"
                } else { "exact UMI, maximum-support gene/state labels; USA spliced preference for cross-gene ties of at most ten labels; C/X groups use conservative intersection" },
                "phase_seconds":{"ingest":ingest_seconds,"resolve_and_collapse":resolve_seconds,"gather":gather_seconds,"write":phase.elapsed().as_secs_f64()},
                "count_seconds":started.elapsed().as_secs_f64()
            }))?,
        )?;
    }
    Ok(())
}

fn resolve_cell(
    mut cell: CellEvidence,
    policy: crate::SplicingResolutionArg,
    correct: bool,
) -> CellCounts {
    use crate::SplicingResolutionArg;
    cell.reads.sort_unstable_by_key(|r| r.umi);
    let mut result = CellCounts::default();
    let mut molecules = Vec::new();
    let mut votes = Vec::new();
    let mut intersection = Molecule {
        reads: 0,
        genes: Vec::new(),
    };
    let mut i = 0;
    while i < cell.reads.len() {
        let begin = i;
        let umi = cell.reads[i].umi;
        while i < cell.reads.len() && cell.reads[i].umi == umi {
            i += 1;
        }
        result.raw += 1;
        intersection.reads = 0;
        intersection.genes.clear();
        votes.clear();
        let mut known = true;
        for read in &cell.reads[begin..i] {
            let candidates = &cell.candidates[read.start as usize..read.end as usize];
            if intersection.reads == 0 {
                intersection.genes.extend_from_slice(candidates);
                intersection.reads = 1;
            } else {
                intersection.add(candidates);
            }
            if policy == SplicingResolutionArg::Usa {
                for &(gene, state) in candidates {
                    match state {
                        SplicingState::Spliced => votes.push((gene, 1u8)),
                        SplicingState::Unspliced => votes.push((gene, 2)),
                        SplicingState::Both => {
                            votes.push((gene, 1));
                            votes.push((gene, 2));
                        }
                        _ => known = false,
                    }
                }
            }
        }
        let conservative = (intersection.genes.len() == 1).then(|| intersection.genes[0]);
        let selected = if policy == SplicingResolutionArg::Usa && known {
            let (selected, tie) = resolve_usa(&mut votes);
            result.spliced_ties += usize::from(tie);
            if let Some((gene, state)) = selected {
                if let Some((old_gene, old_state)) = conservative {
                    result.changed_genes += usize::from(gene != old_gene);
                    result.changed_states += usize::from(gene == old_gene && state != old_state);
                } else {
                    result.rescued += 1;
                }
            }
            selected
        } else {
            if policy == SplicingResolutionArg::Usa {
                result.fallback += 1;
            }
            conservative
        };
        if let Some((gene, state)) = selected {
            molecules.push((gene, umi, intersection.reads, state));
        } else if intersection.genes.is_empty() {
            result.conflicting += 1;
        } else {
            result.unresolved += 1;
        }
    }
    result.before = molecules.len();
    molecules.sort_unstable_by_key(|&(gene, umi, _, _)| (gene, umi));
    let mut i = 0;
    while i < molecules.len() {
        let begin = i;
        let gene = molecules[i].0;
        while i < molecules.len() && molecules[i].0 == gene {
            i += 1;
        }
        let mut counts = [0u32; 5];
        let umis = molecules[begin..i]
            .iter()
            .map(|&(_, u, n, s)| (u, n, s))
            .collect();
        for state in collapse_states(umis, correct) {
            let category = match state {
                SplicingState::Spliced => 0,
                SplicingState::Unspliced => 1,
                SplicingState::Both => 2,
                SplicingState::Conflict => 3,
                SplicingState::Unknown => 4,
            };
            counts[category] += 1;
        }
        result.counts.push((gene, counts));
    }
    result
}

/// The cr-like WinnerTakeAll USA rule, with no UMI edit correction. Both votes
/// once for each class; tied S/U of one gene remains Both. The cross-gene S
/// preference is explicit and opt-in, and does not imply stronger evidence.
fn resolve_usa(votes: &mut [(u32, u8)]) -> (Option<(u32, SplicingState)>, bool) {
    votes.sort_unstable();
    let mut best = Vec::new();
    let mut maximum = 0;
    let mut i = 0;
    while i < votes.len() {
        let label = votes[i];
        let begin = i;
        while i < votes.len() && votes[i] == label {
            i += 1;
        }
        match (i - begin).cmp(&maximum) {
            std::cmp::Ordering::Greater => {
                maximum = i - begin;
                best.clear();
                best.push(label);
            }
            std::cmp::Ordering::Equal => best.push(label),
            _ => (),
        }
    }
    if best.len() == 1 {
        return (
            Some((best[0].0, SplicingState::from_code(best[0].1).unwrap())),
            false,
        );
    }
    if best.len() == 2 && best[0].0 == best[1].0 {
        return (Some((best[0].0, SplicingState::Both)), false);
    }
    if (2..=10).contains(&best.len()) {
        let mut spliced = best.iter().filter(|x| x.1 == 1);
        if let Some(&(gene, _)) = spliced.next() {
            if spliced.next().is_none() {
                let state = if best.contains(&(gene, 2)) {
                    SplicingState::Both
                } else {
                    SplicingState::Spliced
                };
                return (Some((gene, state)), true);
            }
        }
    }
    (None, false)
}

fn write_cell_matrix(
    prefix: &Path,
    label: &str,
    genes: usize,
    cells: usize,
    counts: &[(u32, u32, [u32; 5])],
    states: std::ops::Range<usize>,
) -> Result<()> {
    use std::io::Write;
    let count = |values: &[u32; 5]| values[states.clone()].iter().sum::<u32>();
    let nonzero = counts
        .iter()
        .filter(|(_, _, values)| count(values) > 0)
        .count();
    let mut out = std::io::BufWriter::with_capacity(
        1 << 20,
        File::create(with_suffix(prefix, &format!("_{label}.mtx")))?,
    );
    writeln!(
        out,
        "%%MatrixMarket matrix coordinate integer general\n{genes} {cells} {nonzero}"
    )?;
    for (gene, cell, values) in counts {
        let n = count(values);
        if n > 0 {
            writeln!(out, "{} {} {n}", gene + 1, cell + 1)?;
        }
    }
    out.flush()?;
    Ok(())
}

fn collapse_states(mut umis: Vec<(u64, u32, SplicingState)>, correct: bool) -> Vec<SplicingState> {
    if !correct || umis.len() <= 1 {
        return umis.into_iter().map(|(_, _, s)| s).collect();
    }
    let evidence: FxHashMap<_, _> = umis.iter().map(|&(u, n, s)| (u, (n, s))).collect();
    umis.sort_unstable_by_key(|&(u, n, _)| (Reverse(n), u));
    let mut assigned = FxHashSet::default();
    let mut out = Vec::new();
    for (umi, count, mut state) in umis {
        if !assigned.insert(umi) {
            continue;
        }
        for neighbor in one_mismatch_umis(umi) {
            if assigned.contains(&neighbor) {
                continue;
            }
            if let Some(&(n, s)) = evidence.get(&neighbor) {
                if count >= n.saturating_mul(2).saturating_sub(1) {
                    assigned.insert(neighbor);
                    state = state.intersect(s);
                }
            }
        }
        out.push(state);
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;
    use SplicingState::*;
    #[test]
    fn packed_endianness_preserves_barcode_and_umi_ordering() {
        for sequence in ["ACGTGTCAGTACGATC", "AAAAAAAAAAAAAAAA", "TGCATGCATGCATGCA"] {
            let little = sequence.bytes().enumerate().fold(0u64, |v, (i, b)| {
                v | ((b"ACGT".iter().position(|&x| x == b).unwrap() as u64) << (2 * i))
            });
            assert_eq!(reverse_pairs(little, 16), pack_acgt(sequence).unwrap());
            assert_eq!(
                pack_with_len(reverse_pairs(little & 0xff_ffff, 12), 12),
                pack_umi(&sequence[..12]).unwrap()
            );
        }
    }
    #[test]
    fn usa_keeps_ties_explicit_and_matches_native_spliced_preference() {
        assert_eq!(resolve_usa(&mut [(0, 1), (0, 2)]), (Some((0, Both)), false));
        assert_eq!(
            resolve_usa(&mut [(0, 1), (1, 2)]),
            (Some((0, Spliced)), true)
        );
        assert_eq!(resolve_usa(&mut [(0, 1), (1, 1)]), (None, false));
        assert_eq!(resolve_usa(&mut [(0, 2), (1, 2)]), (None, false));
        assert_eq!(
            resolve_usa(&mut [(0, 1), (0, 2), (1, 2)]),
            (Some((0, Both)), true)
        );
        assert_eq!(
            resolve_usa(&mut [(0, 1), (0, 2), (1, 2), (1, 2)]),
            (Some((1, Unspliced)), false)
        );
        let mut many: Vec<_> = (0..11).map(|g| (g, if g == 0 { 1 } else { 2 })).collect();
        assert_eq!(resolve_usa(&mut many), (None, false));
    }
    #[test]
    fn usa_falls_back_for_conflicting_and_unknown_constraints() {
        let cell = CellEvidence {
            reads: vec![ReadEvidence {
                umi: 1,
                start: 0,
                end: 2,
            }],
            candidates: vec![(0, Unknown), (1, Spliced)],
        };
        let result = resolve_cell(cell, crate::SplicingResolutionArg::Usa, false);
        assert_eq!(result.fallback, 1);
        assert_eq!(result.unresolved, 1);
        assert_eq!(result.before, 0);
    }
    #[test]
    fn ambiguous_reads_resolve_genes_and_state_without_forcing_conflicts() {
        let mut molecule = Molecule {
            reads: 1,
            genes: vec![(0, Both), (1, Unspliced)],
        };
        molecule.add(&[(0, Spliced)]);
        assert_eq!(molecule.genes, vec![(0, Spliced)]);
        molecule.add(&[(0, Unspliced)]);
        assert_eq!(molecule.genes, vec![(0, Conflict)]);
        molecule.add(&[(1, Unspliced)]);
        assert!(molecule.genes.is_empty());
    }
    #[test]
    fn umi_correction_merges_state_constraints_not_separate_state_counts() {
        let a = pack_umi("AAAAAAAAAAAA").unwrap();
        let b = pack_umi("AAAAAAAAAAAC").unwrap();
        assert_eq!(
            collapse_states(vec![(a, 5, Both), (b, 2, Spliced)], true),
            vec![Spliced]
        );
        assert_eq!(
            collapse_states(vec![(a, 5, Spliced), (b, 2, Unspliced)], true),
            vec![Conflict]
        );
        assert_eq!(
            collapse_states(vec![(a, 5, Both), (b, 2, Spliced)], false),
            vec![Both, Spliced]
        );
    }
}

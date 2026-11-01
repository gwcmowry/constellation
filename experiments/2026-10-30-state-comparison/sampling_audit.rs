//! Read-only check of exact class-discriminating k-mers against the actual index.
use constellation_core::dna::{encode_acgt, iter_closed_syncmers_2bit, iter_kmers_2bit};
use constellation_core::ec_index::{LoadedEcIndex, EcMappingStats};
use std::collections::{BTreeMap, HashMap};
use std::io::{BufRead, Write};

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let args: Vec<_> = std::env::args().collect();
    let index = LoadedEcIndex::load(&args[1])?;
    assert_eq!(index.k(), 31);
    println!("sampling: {:?}", index.sampling());
    let genes: HashMap<_, _> = (0..index.num_genes())
        .map(|g| (index.gene_name(g as u32).unwrap().to_owned(), g as u32)).collect();
    let mut counts = BTreeMap::new();
    let mut total = BTreeMap::new();
    for line in std::io::BufReader::new(std::fs::File::open(&args[2])?).lines() {
        let line = line?;
        let fields: Vec<_> = line.split('\t').collect();
        assert_eq!(fields.len(), 4);
        let gene = genes[fields[1]];
        let extra: u8 = fields[2].parse()?;
        let encoded = encode_acgt(fields[3].as_bytes());
        assert_eq!(fields[3].len(), 31);
        let code = iter_kmers_2bit(&encoded, 31).next().unwrap().code;
        let selected = iter_closed_syncmers_2bit(&encoded, 31, 20).next().is_some();
        let hit = index.lookup_hit(code);
        let classification = match (selected, hit) {
            (false, None) => "not_selected",
            (false, Some(_)) => "unexpected_unselected_present",
            (true, None) => "unexpected_selected_missing",
            (true, Some(hit)) => {
                let df = index.lookup(code).unwrap().gene_df;
                let mut scores = Default::default();
                index.accumulate_lookup_hit(hit, u32::MAX, &mut scores, &mut EcMappingStats::default());
                let mask = scores.get(&gene).map(|s| u8::from(s.exon_score>0) | (u8::from(s.intron_score>0)<<1)).unwrap_or(0);
                if mask & extra != 0 { "unexpected_extra_class_present" }
                else if mask == 0 { "unexpected_source_gene_absent" }
                else if df > 256 { "filtered_by_gene_df" }
                else { "unexpected_usable_discriminator" }
            }
        };
        *counts.entry((fields[0].to_owned(), classification)).or_insert(0u32) += 1;
        *total.entry(classification).or_insert(0u32) += 1;
    }
    let mut out = std::io::BufWriter::new(std::fs::File::create(&args[3])?);
    writeln!(out, "case\tclassification\tkmers")?;
    for ((case, classification), n) in counts { writeln!(out, "{case}\t{classification}\t{n}")?; }
    for (classification, n) in total { println!("{classification}\t{n}"); }
    Ok(())
}

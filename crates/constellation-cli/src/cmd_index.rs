use crate::{IndexArgs, InspectIndexArgs};
use anyhow::Result;
use constellation_core::gtf::transcript_gene_map;
use constellation_core::index::{IndexAccess, LoadedIndex};
use constellation_core::index_build::{
    build_transcript_index, build_transcript_index_with_gene_map,
};

pub fn run_index(args: IndexArgs) -> Result<()> {
    let index = if let Some(gtf) = args.gtf {
        let gene_map = transcript_gene_map(gtf)?;
        build_transcript_index_with_gene_map(
            args.transcripts,
            args.k,
            args.max_kmer_frequency,
            Some(&gene_map),
        )?
    } else {
        build_transcript_index(args.transcripts, args.k, args.max_kmer_frequency)?
    };
    index.save_auto(args.out)?;
    Ok(())
}

pub fn run_inspect_index(args: InspectIndexArgs) -> Result<()> {
    let index = LoadedIndex::load(args.index)?;
    let stats = index.stats();
    println!("num_transcripts\t{}", stats.num_transcripts);
    println!("num_genes\t{}", stats.num_genes);
    println!("num_distinct_kmers\t{}", stats.num_distinct_kmers);
    println!("num_postings\t{}", stats.num_postings);
    println!("k\t{}", stats.k);
    println!("max_postings_per_kmer\t{}", stats.max_postings_per_kmer);
    println!("high_frequency_kmers\t{}", stats.high_frequency_kmers);
    Ok(())
}

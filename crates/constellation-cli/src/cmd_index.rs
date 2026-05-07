use crate::{ConvertEcIndexArgs, EcIndexFormatArg, IndexArgs, IndexEcArgs, InspectIndexArgs};
use anyhow::Result;
use constellation_core::ec_index::{EcIndex, LoadedEcIndex};
use constellation_core::gtf::transcript_gene_map;
use constellation_core::index::{IndexAccess, LoadedIndex};
use constellation_core::index_build::{
    build_compact_transcript_index_streaming, build_transcript_index,
    build_transcript_index_with_gene_map,
};

pub fn run_index(args: IndexArgs) -> Result<()> {
    let use_streaming_compact = args.gtf.is_none()
        && args
            .out
            .extension()
            .is_none_or(|extension| extension != "json")
        && std::fs::metadata(&args.transcripts)
            .map(|metadata| metadata.len() > 1_000_000_000)
            .unwrap_or(false);
    if use_streaming_compact {
        return build_compact_transcript_index_streaming(
            args.transcripts,
            args.k,
            args.max_kmer_frequency,
            args.out,
        )
        .map_err(Into::into);
    }

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

pub fn run_index_ec(args: IndexEcArgs) -> Result<()> {
    let index =
        EcIndex::build_from_fasta_with_t2g(args.transcripts, args.t2g_map.as_ref(), args.k)?;
    match args.format {
        EcIndexFormatArg::Mmap => index.save_mmap(args.out)?,
        EcIndexFormatArg::Legacy => index.save(args.out)?,
    }
    Ok(())
}

pub fn run_convert_ec_index(args: ConvertEcIndexArgs) -> Result<()> {
    let index = LoadedEcIndex::load(args.index)?;
    index.save_mmap(args.out)?;
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
    print_histogram("raw_postings", &stats.raw_postings_histogram);
    print_histogram("transcript_df", &stats.transcript_df_histogram);
    print_histogram("gene_df", &stats.gene_df_histogram);
    Ok(())
}

fn print_histogram(name: &str, hist: &constellation_core::index::KmerHistogram) {
    println!("{name}_hist_0\t{}", hist.zero);
    println!("{name}_hist_1\t{}", hist.one);
    println!("{name}_hist_2_4\t{}", hist.two_to_four);
    println!("{name}_hist_5_16\t{}", hist.five_to_sixteen);
    println!("{name}_hist_17_64\t{}", hist.seventeen_to_sixty_four);
    println!("{name}_hist_65_256\t{}", hist.sixty_five_to_256);
    println!("{name}_hist_gt_256\t{}", hist.above_256);
}

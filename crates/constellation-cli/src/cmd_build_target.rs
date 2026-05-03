use crate::BuildTranscriptomeTargetArgs;
use anyhow::Result;
use constellation_core::transcriptome_target::{
    build_transcriptome_target, TranscriptomeTargetKind,
};

pub fn run_build_transcriptome_target(args: BuildTranscriptomeTargetArgs) -> Result<()> {
    let stats = build_transcriptome_target(
        args.genome,
        args.gtf,
        args.out,
        TranscriptomeTargetKind::ExonTranscripts,
    )?;
    println!("num_genome_contigs\t{}", stats.num_genome_contigs);
    println!("num_transcripts\t{}", stats.num_transcripts);
    println!("num_exons\t{}", stats.num_exons);
    println!("total_bases\t{}", stats.total_bases);
    Ok(())
}

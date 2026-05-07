use crate::BuildTranscriptomeTargetArgs;
use anyhow::Result;
use constellation_core::transcriptome_target::{
    build_transcriptome_target, TranscriptomeTargetKind,
};

pub fn run_build_transcriptome_target(args: BuildTranscriptomeTargetArgs) -> Result<()> {
    let kind = match args.target_kind {
        crate::TargetKindArg::ExonTranscripts => TranscriptomeTargetKind::ExonTranscripts,
        crate::TargetKindArg::GeneBodies => TranscriptomeTargetKind::GeneBodies,
        crate::TargetKindArg::IntronsOnly => TranscriptomeTargetKind::IntronsOnly,
        crate::TargetKindArg::IntronFlanks => TranscriptomeTargetKind::IntronFlanks,
        crate::TargetKindArg::ExonPlusGeneBody => TranscriptomeTargetKind::ExonPlusGeneBody,
        crate::TargetKindArg::ExonPlusIntronsOnly => TranscriptomeTargetKind::ExonPlusIntronsOnly,
    };
    let stats = build_transcriptome_target(args.genome, args.gtf, args.out, kind)?;
    println!("num_genome_contigs\t{}", stats.num_genome_contigs);
    println!("num_transcripts\t{}", stats.num_transcripts);
    println!("num_exons\t{}", stats.num_exons);
    println!("total_bases\t{}", stats.total_bases);
    Ok(())
}

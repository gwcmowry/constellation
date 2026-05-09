use anyhow::Result;
use clap::{ArgAction, Parser, Subcommand, ValueEnum};
use std::path::PathBuf;

mod cmd_bench_report;
mod cmd_build_target;
mod cmd_count;
mod cmd_index;
mod cmd_map;
mod cmd_simulate;

#[derive(Debug, Parser)]
#[command(name = "constellation")]
#[command(about = "Cache-aware RNA-seq/scRNA-seq mapper prototype")]
struct Cli {
    #[command(subcommand)]
    command: Commands,
}

#[derive(Debug, Subcommand)]
enum Commands {
    Index(IndexArgs),
    IndexEc(IndexEcArgs),
    ConvertEcIndex(ConvertEcIndexArgs),
    InspectIndex(InspectIndexArgs),
    Simulate(SimulateArgs),
    BuildTranscriptomeTarget(BuildTranscriptomeTargetArgs),
    Map(MapArgs),
    Count(CountArgs),
    BenchReport(BenchReportArgs),
}

#[derive(Debug, Parser)]
struct IndexArgs {
    #[arg(long)]
    transcripts: PathBuf,
    #[arg(long)]
    gtf: Option<PathBuf>,
    #[arg(long, default_value_t = 21)]
    k: u8,
    #[arg(long, default_value_t = 256)]
    max_kmer_frequency: u32,
    #[arg(long)]
    out: PathBuf,
}

#[derive(Debug, Parser)]
struct IndexEcArgs {
    #[arg(long)]
    transcripts: PathBuf,
    #[arg(long)]
    t2g_map: Option<PathBuf>,
    #[arg(long, default_value_t = 31)]
    k: u8,
    #[arg(long, value_enum, default_value_t = EcIndexFormatArg::Mmap)]
    format: EcIndexFormatArg,
    #[arg(long)]
    out: PathBuf,
}

#[derive(Debug, Parser)]
struct ConvertEcIndexArgs {
    #[arg(long)]
    index: PathBuf,
    #[arg(long)]
    out: PathBuf,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, ValueEnum)]
enum EcIndexFormatArg {
    Mmap,
    Legacy,
}

#[derive(Debug, Parser)]
struct InspectIndexArgs {
    #[arg(long)]
    index: PathBuf,
}

#[derive(Debug, Parser)]
struct SimulateArgs {
    #[arg(long)]
    transcripts: PathBuf,
    #[arg(long, default_value_t = 1000)]
    num_reads: usize,
    #[arg(long, default_value_t = 91)]
    read_len: usize,
    #[arg(long, default_value_t = 0.0)]
    error_rate: f64,
    #[arg(long, default_value = "exact")]
    scenario: String,
    #[arg(long, default_value = "tenx-3p-v3")]
    chemistry: String,
    #[arg(long)]
    out_prefix: PathBuf,
}

#[derive(Debug, Parser)]
struct BuildTranscriptomeTargetArgs {
    #[arg(long)]
    genome: PathBuf,
    #[arg(long)]
    gtf: PathBuf,
    #[arg(long, value_enum, default_value_t = TargetKindArg::ExonTranscripts)]
    target_kind: TargetKindArg,
    #[arg(long)]
    out: PathBuf,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, ValueEnum)]
enum TargetKindArg {
    ExonTranscripts,
    GeneBodies,
    IntronsOnly,
    IntronFlanks,
    ExonPlusGeneBody,
    ExonPlusIntronsOnly,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, ValueEnum)]
enum ScoreMode {
    Scalar,
    Pulp,
}

impl ScoreMode {
    fn as_str(self) -> &'static str {
        match self {
            Self::Scalar => "scalar",
            Self::Pulp => "pulp",
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, ValueEnum)]
enum MapMode {
    ReadAtATime,
    SketchBucket,
    CandidateLocus,
    GeneEc,
}

impl MapMode {
    fn as_str(self) -> &'static str {
        match self {
            Self::ReadAtATime => "read-at-a-time",
            Self::SketchBucket => "sketch-bucket",
            Self::CandidateLocus => "candidate-locus",
            Self::GeneEc => "gene-ec",
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, ValueEnum)]
enum SeedPlannerArg {
    RawFrequency,
    GeneIdf,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, ValueEnum)]
enum RetrievalModeArg {
    PerRead,
    SeedBatched,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, ValueEnum)]
enum CandidateSearchArg {
    Full,
    ExonFirst,
    SparseProbe,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, ValueEnum)]
enum ReadChunkingArg {
    Sketch,
    BestSeed,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, ValueEnum)]
enum CandidatePruningArg {
    None,
    WandLite,
    GeneWandLite,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, ValueEnum)]
enum MapOutputFormatArg {
    Tsv,
    GeneEcRad,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, ValueEnum)]
enum OutputCompressionArg {
    None,
    Zstd,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, ValueEnum)]
enum LibraryStrandArg {
    Unstranded,
    Forward,
    Reverse,
}

#[derive(Debug, Parser)]
struct MapArgs {
    #[arg(long)]
    index: PathBuf,
    #[arg(long, hide = true)]
    hot_index: Option<PathBuf>,
    #[arg(long)]
    r1: PathBuf,
    #[arg(long)]
    r2: PathBuf,
    #[arg(long, default_value = "tenx-3p-v3")]
    chemistry: String,
    #[arg(long, default_value_t = 65_536)]
    batch_size: usize,
    #[arg(long, hide = true)]
    cold_batch_size: Option<usize>,
    #[arg(long, default_value_t = 1, hide = true)]
    cold_max_seeds_per_read: usize,
    #[arg(long, default_value_t = 10, hide = true)]
    cold_shard_prefix_bits: u8,
    #[arg(long, default_value_t = 32, hide = true)]
    cold_evict_interval: u64,
    #[arg(long = "no-cold-shard-scheduling", action = ArgAction::SetFalse, default_value_t = true, hide = true)]
    cold_shard_scheduling: bool,
    #[arg(long, default_value_t = false, hide = true)]
    preload_cold_index: bool,
    #[arg(long, value_enum, default_value_t = ScoreMode::Scalar, hide = true)]
    score_mode: ScoreMode,
    #[arg(long, value_enum, default_value_t = MapMode::GeneEc, hide = true)]
    mode: MapMode,
    #[arg(long, default_value_t = 8, hide = true)]
    max_seeds_per_read: usize,
    #[arg(long, value_enum, default_value_t = SeedPlannerArg::RawFrequency, hide = true)]
    seed_planner: SeedPlannerArg,
    #[arg(long, value_enum, default_value_t = RetrievalModeArg::PerRead, hide = true)]
    retrieval_mode: RetrievalModeArg,
    #[arg(long, value_enum, default_value_t = CandidateSearchArg::Full, hide = true)]
    candidate_search: CandidateSearchArg,
    #[arg(long, value_enum, default_value_t = CandidatePruningArg::GeneWandLite, hide = true)]
    candidate_pruning: CandidatePruningArg,
    #[arg(long, default_value_t = 2, hide = true)]
    wand_min_top_seed_count: u16,
    #[arg(long, default_value_t = 20, hide = true)]
    wand_score_ratio_percent: u16,
    #[arg(long, default_value_t = 4, hide = true)]
    wand_max_loci_per_gene: usize,
    #[arg(long, value_enum, default_value_t = ReadChunkingArg::Sketch, hide = true)]
    read_chunking: ReadChunkingArg,
    #[arg(long, default_value_t = 12)]
    sparse_probe_stride: u32,
    #[arg(long, default_value_t = 3, hide = true)]
    sparse_probe_max_seeds: usize,
    #[arg(long, default_value_t = 2, hide = true)]
    sparse_probe_min_seed_hits: u16,
    #[arg(long, default_value_t = 256)]
    max_postings_per_seed: usize,
    #[arg(long, default_value_t = 64, hide = true)]
    candidate_bin_size: u32,
    #[arg(long, default_value_t = 10.0)]
    min_mean_quality: f64,
    #[arg(long, default_value_t = 10, hide = true)]
    min_seed_quality: u8,
    #[arg(long, default_value_t = 6, hide = true)]
    max_mismatches: u32,
    #[arg(long, default_value_t = 0, hide = true)]
    max_right_softclip: u16,
    #[arg(long, default_value_t = 0, hide = true)]
    max_left_softclip: u16,
    #[arg(long, default_value_t = false, hide = true)]
    trim_poly_a: bool,
    #[arg(long, default_value_t = false, hide = true)]
    trim_poly_t: bool,
    #[arg(long, default_value_t = false, hide = true)]
    trim_low_quality_tail: bool,
    #[arg(long, default_value_t = false)]
    trim_tso: bool,
    #[arg(long, default_value_t = 3)]
    max_tso_mismatches: u8,
    #[arg(long, default_value_t = 10)]
    min_tso_match_len: u8,
    #[arg(long, value_enum, default_value_t = LibraryStrandArg::Unstranded, hide = true)]
    library_strand: LibraryStrandArg,
    #[arg(long, default_value_t = 35, hide = true)]
    min_scored_length: u16,
    #[arg(long, default_value_t = 10, hide = true)]
    min_tail_phred: u8,
    #[arg(long, default_value_t = false)]
    search_reverse_complement: bool,
    #[arg(long, default_value_t = false, hide = true)]
    rescue_failed_reads: bool,
    #[arg(long, default_value_t = false, hide = true)]
    rescue_score_failed_reads: bool,
    #[arg(long, default_value_t = true, hide = true)]
    rescue_antisense_reads: bool,
    #[arg(long, default_value_t = true, hide = true)]
    rescue_reverse_complement: bool,
    #[arg(long, default_value_t = 1.0, hide = true)]
    early_stop_posterior: f64,
    #[arg(long, default_value_t = 0.25, hide = true)]
    early_stop_prior_alpha: f64,
    #[arg(long, default_value_t = 4, hide = true)]
    early_stop_min_seed_lookups: usize,
    #[arg(long, default_value_t = 4, hide = true)]
    early_stop_min_top_seed_count: u16,
    #[arg(long)]
    emit_metrics: Option<PathBuf>,
    #[arg(long)]
    emit_unmapped_diagnostics: Option<PathBuf>,
    #[arg(long, default_value_t = false)]
    skip_assignments: bool,
    #[arg(long, value_enum, default_value_t = MapOutputFormatArg::GeneEcRad)]
    output_format: MapOutputFormatArg,
    #[arg(long, value_enum, default_value_t = OutputCompressionArg::None)]
    output_compression: OutputCompressionArg,
    #[arg(long, default_value_t = 3)]
    zstd_level: i32,
    #[arg(long)]
    out: PathBuf,
}

#[derive(Debug, Parser)]
struct BenchReportArgs {
    #[arg(long)]
    assignments: PathBuf,
    #[arg(long)]
    truth: Option<PathBuf>,
    #[arg(long)]
    index: Option<PathBuf>,
    #[arg(long)]
    metrics: Option<PathBuf>,
    #[arg(long)]
    perf_stat: Option<PathBuf>,
}

#[derive(Debug, Parser)]
struct CountArgs {
    #[arg(long)]
    assignments: PathBuf,
    #[arg(long)]
    index: PathBuf,
    #[arg(long)]
    barcode_whitelist: Option<PathBuf>,
    #[arg(long = "no-correct-barcodes", action = ArgAction::SetFalse, default_value_t = true)]
    correct_barcodes: bool,
    #[arg(long = "no-correct-umis", action = ArgAction::SetFalse, default_value_t = true)]
    correct_umis: bool,
    #[arg(long, default_value_t = 1)]
    umi_edit_distance: u8,
    #[arg(long = "no-resolve-molecule-genes", action = ArgAction::SetFalse, default_value_t = true)]
    resolve_molecule_genes: bool,
    #[arg(long, default_value_t = 2)]
    molecule_gene_ratio: u32,
    #[arg(long)]
    feature_whitelist: Option<PathBuf>,
    #[arg(long)]
    emit_metrics: Option<PathBuf>,
    #[arg(long)]
    out_prefix: PathBuf,
}

fn main() -> Result<()> {
    let cli = Cli::parse();
    match cli.command {
        Commands::Index(args) => cmd_index::run_index(args),
        Commands::IndexEc(args) => cmd_index::run_index_ec(args),
        Commands::ConvertEcIndex(args) => cmd_index::run_convert_ec_index(args),
        Commands::InspectIndex(args) => cmd_index::run_inspect_index(args),
        Commands::Simulate(args) => cmd_simulate::run_simulate(args),
        Commands::BuildTranscriptomeTarget(args) => {
            cmd_build_target::run_build_transcriptome_target(args)
        }
        Commands::Map(args) => cmd_map::run_map(args),
        Commands::Count(args) => cmd_count::run_count(args),
        Commands::BenchReport(args) => cmd_bench_report::run_bench_report(args),
    }
}

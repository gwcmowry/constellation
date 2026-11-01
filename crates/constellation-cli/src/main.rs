use anyhow::Result;
use clap::{ArgAction, Parser, Subcommand, ValueEnum};
use std::num::NonZeroUsize;
use std::path::PathBuf;

mod cmd_bench_report;
mod cmd_build_target;
mod cmd_count;
mod cmd_index;
mod cmd_map;
mod cmd_quant;
mod cmd_simulate;
mod splicing_io;

#[derive(Debug, Parser)]
#[command(name = "constellation", version)]
#[command(about = "Fast gene and splicing-state mapping/counting for 10x 3' v3 reads")]
struct Cli {
    /// Rayon workers; otherwise use RAYON_NUM_THREADS or the available CPU count.
    /// Input/output helper threads are additional; this does not set CPU affinity.
    #[arg(short = 't', long, global = true)]
    threads: Option<NonZeroUsize>,
    #[command(subcommand)]
    command: Commands,
}

#[derive(Debug, Subcommand)]
enum Commands {
    /// Build a legacy positional index (not used by the current gene-EC mapper).
    Index(IndexArgs),
    /// Build a gene/state index; recommended: --format ptrhash --syncmer-s 20.
    IndexEc(IndexEcArgs),
    ConvertEcIndex(ConvertEcIndexArgs),
    InspectIndex(InspectIndexArgs),
    Simulate(SimulateArgs),
    BuildTranscriptomeTarget(BuildTranscriptomeTargetArgs),
    /// Map paired FASTQs to gene assignments; --splicing-state retains S/U evidence.
    Map(MapArgs),
    /// Map and count a sample into gene-by-cell splicing-state matrices.
    Quant(QuantArgs),
    Count(CountArgs),
    /// Resolve gene/splicing compatibility through UMI counting from CSTRAD2.
    CountSplicing(SplicingCountArgs),
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
    /// Build a closed-syncmer sampled index using this s-mer length.
    #[arg(long)]
    syncmer_s: Option<u8>,
    #[arg(long, value_enum, default_value_t = EcIndexFormatArg::Mmap)]
    format: EcIndexFormatArg,
    /// Partitioned parallel construction, or the lower-memory streaming builder.
    #[arg(long, value_enum, default_value_t = EcBuildModeArg::Partitioned)]
    build_mode: EcBuildModeArg,
    #[arg(long)]
    out: PathBuf,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, ValueEnum)]
enum EcBuildModeArg {
    Partitioned,
    Streaming,
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
    /// Build the fast exact PtrHash lookup directly (no intermediate index).
    Ptrhash,
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
enum GeneEcTraversalArg {
    ReadAtATime,
    SeedBatched,
    PrefixTiled,
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
    /// Preserve per-gene S/U compatibility. Enables forward-strand consensus mapping.
    #[arg(long)]
    splicing_state: bool,
    /// Require agreement of at least two distinct matched seeds; preserve cross-gene ties.
    /// Uses forward-strand read-at-a-time mapping; score is supporting seed count.
    #[arg(long)]
    gene_ec_consensus: bool,
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
    #[arg(long, value_enum, default_value_t = GeneEcTraversalArg::ReadAtATime, hide = true)]
    gene_ec_traversal: GeneEcTraversalArg,
    #[arg(long, default_value_t = 12)]
    sparse_probe_stride: u32,
    #[arg(long, default_value_t = false, hide = true)]
    gene_ec_adaptive_rescue: bool,
    #[arg(long, hide = true)]
    gene_ec_rescue_stride: Option<u32>,
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

#[derive(Debug, Clone, Copy, PartialEq, Eq, ValueEnum)]
enum SplicingResolutionArg {
    Intersection,
    Usa,
}
impl SplicingResolutionArg {
    fn label(self) -> &'static str {
        match self {
            Self::Intersection => "intersection",
            Self::Usa => "usa",
        }
    }
}

#[derive(Debug, Parser)]
struct SplicingCountArgs {
    /// CSTRAD2 produced by map --splicing-state (plain, gzip, or zstd).
    #[arg(long)]
    assignments: PathBuf,
    #[arg(long)]
    barcode_whitelist: Option<PathBuf>,
    #[arg(long = "no-correct-barcodes", action = ArgAction::SetFalse, default_value_t = true)]
    correct_barcodes: bool,
    #[arg(long = "no-correct-umis", action = ArgAction::SetFalse, default_value_t = true)]
    correct_umis: bool,
    /// Molecule resolver; usa applies alevin-fry cr-like/WTA spliced tie preference.
    /// Currently requires --no-correct-umis.
    #[arg(long, value_enum, default_value_t = SplicingResolutionArg::Intersection)]
    resolution: SplicingResolutionArg,
    #[arg(long)]
    emit_metrics: Option<PathBuf>,
    #[arg(long)]
    out_prefix: PathBuf,
}

#[derive(Debug, Parser)]
struct QuantArgs {
    /// Reusable gene-EC index with spliced/intronic reference annotations.
    #[arg(long)]
    index: PathBuf,
    #[arg(long)]
    r1: PathBuf,
    #[arg(long)]
    r2: PathBuf,
    /// Accepted cell barcodes, one per line; plain or gzip, optional -1 suffix.
    /// This workflow uses the supplied cells; it does not perform cell calling.
    #[arg(long, alias = "barcodes")]
    barcode_whitelist: PathBuf,
    /// Create a new output directory; an existing path is never overwritten.
    #[arg(short = 'o', long)]
    out_dir: PathBuf,
    #[arg(long, default_value = "tenx-3p-v3")]
    chemistry: String,
    /// Conservative intersection by default; usa automatically uses exact UMIs.
    #[arg(long, value_enum, default_value_t = SplicingResolutionArg::Intersection)]
    resolution: SplicingResolutionArg,
    /// Disable directional UMI correction for intersection (implicit with usa).
    #[arg(long)]
    no_correct_umis: bool,
    /// Accept only exact whitelist barcodes.
    #[arg(long)]
    no_correct_barcodes: bool,
    /// Retain assignments.cstrad2 after successful counting for later re-counting.
    #[arg(long)]
    keep_assignments: bool,
}

fn main() -> Result<()> {
    let cli = Cli::parse();
    // Quant launches the same native stages in separate processes, so mapping
    // memory is released before counting. Only its children need worker pools.
    if !matches!(&cli.command, Commands::Quant(_)) {
        if let Some(threads) = cli.threads {
            rayon::ThreadPoolBuilder::new()
                .num_threads(threads.get())
                .build_global()?;
        }
    }
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
        Commands::Quant(args) => cmd_quant::run(args, cli.threads),
        Commands::Count(args) => cmd_count::run_count(args),
        Commands::CountSplicing(args) => cmd_count::splicing::run(args),
        Commands::BenchReport(args) => cmd_bench_report::run_bench_report(args),
    }
}

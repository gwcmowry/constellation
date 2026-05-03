use anyhow::Result;
use clap::{Parser, Subcommand, ValueEnum};
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
    #[arg(long)]
    out: PathBuf,
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
}

impl MapMode {
    fn as_str(self) -> &'static str {
        match self {
            Self::ReadAtATime => "read-at-a-time",
            Self::SketchBucket => "sketch-bucket",
            Self::CandidateLocus => "candidate-locus",
        }
    }
}

#[derive(Debug, Parser)]
struct MapArgs {
    #[arg(long)]
    index: PathBuf,
    #[arg(long)]
    r1: PathBuf,
    #[arg(long)]
    r2: PathBuf,
    #[arg(long, default_value = "tenx-3p-v3")]
    chemistry: String,
    #[arg(long, default_value_t = 2048)]
    batch_size: usize,
    #[arg(long, value_enum, default_value_t = ScoreMode::Scalar)]
    score_mode: ScoreMode,
    #[arg(long, value_enum, default_value_t = MapMode::CandidateLocus)]
    mode: MapMode,
    #[arg(long, default_value_t = 8)]
    max_seeds_per_read: usize,
    #[arg(long, default_value_t = 256)]
    max_postings_per_seed: usize,
    #[arg(long, default_value_t = 64)]
    candidate_bin_size: u32,
    #[arg(long, default_value_t = 10.0)]
    min_mean_quality: f64,
    #[arg(long, default_value_t = 10)]
    min_seed_quality: u8,
    #[arg(long, default_value_t = false)]
    search_reverse_complement: bool,
    #[arg(long, default_value_t = 1.0)]
    early_stop_posterior: f64,
    #[arg(long, default_value_t = 0.25)]
    early_stop_prior_alpha: f64,
    #[arg(long, default_value_t = 4)]
    early_stop_min_seed_lookups: usize,
    #[arg(long, default_value_t = 4)]
    early_stop_min_top_seed_count: u16,
    #[arg(long)]
    emit_metrics: Option<PathBuf>,
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
    out_prefix: PathBuf,
}

fn main() -> Result<()> {
    let cli = Cli::parse();
    match cli.command {
        Commands::Index(args) => cmd_index::run_index(args),
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

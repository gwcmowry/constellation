use crate::SimulateArgs;
use anyhow::{anyhow, Result};
use constellation_core::chemistry::Chemistry;
use constellation_core::simulate::{simulate, SimulateConfig, SimulationScenario};

pub fn run_simulate(args: SimulateArgs) -> Result<()> {
    let _chemistry = Chemistry::parse(&args.chemistry)?;
    simulate(&SimulateConfig {
        transcripts: args.transcripts,
        num_reads: args.num_reads,
        read_len: args.read_len,
        error_rate: args.error_rate,
        scenario: SimulationScenario::parse(&args.scenario)
            .ok_or_else(|| anyhow!("unknown simulation scenario: {}", args.scenario))?,
        out_prefix: args.out_prefix,
    })?;
    Ok(())
}

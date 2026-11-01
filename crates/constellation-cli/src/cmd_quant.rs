//! Sample-level orchestration of the existing native mapping/counting stages.
use crate::{cmd_count, QuantArgs, SplicingResolutionArg};
use anyhow::{ensure, Context, Result};
use constellation_core::chemistry::Chemistry;
use serde_json::{json, Value};
use std::{fs, num::NonZeroUsize, path::Path, process::Command, time::Instant};

fn read_json(path: &Path) -> Result<Value> {
    Ok(serde_json::from_slice(&fs::read(path)?)?)
}

fn stage(command: &mut Command, name: &str) -> Result<f64> {
    eprintln!("constellation quant: {name}");
    let start = Instant::now();
    let status = command
        .status()
        .with_context(|| format!("could not start {name}"))?;
    ensure!(status.success(), "{name} failed ({status})");
    Ok(start.elapsed().as_secs_f64())
}

pub fn run(args: QuantArgs, threads: Option<NonZeroUsize>) -> Result<()> {
    // Validate cheap inputs before spending time mapping or creating outputs.
    Chemistry::parse(&args.chemistry)?;
    for (name, path) in [
        ("index", &args.index),
        ("R1", &args.r1),
        ("R2", &args.r2),
        ("barcode whitelist", &args.barcode_whitelist),
    ] {
        ensure!(
            fs::metadata(path)
                .with_context(|| format!("cannot access {name}: {}", path.display()))?
                .is_file(),
            "{name} must be a regular file: {}",
            path.display()
        );
        fs::File::open(path).with_context(|| format!("cannot open {name}: {}", path.display()))?;
    }
    cmd_count::load_barcode_whitelist(&args.barcode_whitelist)?;
    if let Some(parent) = args.out_dir.parent().filter(|p| !p.as_os_str().is_empty()) {
        fs::create_dir_all(parent)?;
    }
    // create_dir, rather than create_dir_all, also rejects a concurrent run
    // claiming the same directory. Partial outputs are preserved on failure.
    fs::create_dir(&args.out_dir).with_context(|| {
        format!("cannot create output directory {}; choose a new --out-dir (existing paths are not overwritten)", args.out_dir.display())
    })?;

    let executable = std::env::current_exe()?;
    let command = || {
        let mut c = Command::new(&executable);
        if let Some(n) = threads {
            c.arg("--threads").arg(n.to_string());
        }
        c
    };
    let assignments = args.out_dir.join("assignments.cstrad2");
    let map_metrics = args.out_dir.join("map.metrics.json");
    let count_metrics = args.out_dir.join("count.metrics.json");
    let prefix = args.out_dir.join("counts");
    let correct_umis = args.resolution != SplicingResolutionArg::Usa && !args.no_correct_umis;
    let started = Instant::now();
    let mut map_seconds = None;
    let mut count_seconds = None;
    let result = (|| -> Result<()> {
        map_seconds = Some(stage(
            command()
                .arg("map")
                .arg("--splicing-state")
                .arg("--index")
                .arg(&args.index)
                .arg("--r1")
                .arg(&args.r1)
                .arg("--r2")
                .arg(&args.r2)
                .arg("--chemistry")
                .arg(&args.chemistry)
                .arg("--out")
                .arg(&assignments)
                .arg("--emit-metrics")
                .arg(&map_metrics),
            "mapping",
        )?);
        let mut count = command();
        count
            .arg("count-splicing")
            .arg("--assignments")
            .arg(&assignments)
            .arg("--barcode-whitelist")
            .arg(&args.barcode_whitelist)
            .arg("--resolution")
            .arg(args.resolution.label())
            .arg("--out-prefix")
            .arg(&prefix)
            .arg("--emit-metrics")
            .arg(&count_metrics);
        if !correct_umis {
            count.arg("--no-correct-umis");
        }
        if args.no_correct_barcodes {
            count.arg("--no-correct-barcodes");
        }
        count_seconds = Some(stage(&mut count, "counting")?);
        let mapping = read_json(&map_metrics)?;
        let counting = read_json(&count_metrics)?;
        let reads = mapping["num_reads"]
            .as_u64()
            .context("missing mapping read count")?;
        ensure!(
            counting["rows"].as_u64() == Some(reads),
            "mapping/counting record totals differ"
        );
        if !args.keep_assignments {
            fs::remove_file(&assignments)
                .context("could not remove completed mapping intermediate")?;
        }
        Ok(())
    })();
    let manifest = json!({
        "version": env!("CARGO_PKG_VERSION"),
        "status": if result.is_ok() { "complete" } else { "failed" },
        "error": result.as_ref().err().map(|e| format!("{e:#}")),
        "command": std::env::args_os().map(|s| s.to_string_lossy().into_owned()).collect::<Vec<_>>(),
        "resolution": args.resolution.label(), "correct_umis": correct_umis,
        "correct_barcodes": !args.no_correct_barcodes,
        "threads_requested": threads.map(NonZeroUsize::get),
        "rayon_num_threads_env": std::env::var("RAYON_NUM_THREADS").ok(),
        "keep_assignments": args.keep_assignments,
        "map_seconds": map_seconds, "count_seconds": count_seconds,
        "total_seconds": started.elapsed().as_secs_f64()
    });
    fs::write(
        args.out_dir.join("run.json"),
        serde_json::to_vec_pretty(&manifest)?,
    )?;
    result.with_context(|| {
        format!(
            "incomplete quant outputs retained in {}",
            args.out_dir.display()
        )
    })?;
    eprintln!(
        "constellation quant: complete; matrices and metrics in {}",
        args.out_dir.display()
    );
    Ok(())
}

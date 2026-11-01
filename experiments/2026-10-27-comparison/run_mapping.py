"""Serial, affinity-limited native-output mapping benchmarks and accuracy runs."""
import argparse
import fcntl
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import time

ROOT = Path(__file__).resolve().parents[2]
DATA = ROOT / "benchdata/comparison-2026-10-27"
RAW = Path(__file__).parent / "raw"
OLD = ROOT / "benchdata/throughput-2026-10-27"
REF = ROOT / "benchdata/reference/ensembl93/simpleaf_splici_r91_k31"
BINARIES = {
    "constellation": OLD / "constellation.optimized",
    "piscem": DATA / "tools/piscem-x86_64-unknown-linux-gnu/piscem",
    "kallisto": DATA / "tools/kallisto/kallisto",
}
BINARIES["constellation-consensus"] = ROOT / "target/release/constellation"
BINARIES["piscem-forward"] = BINARIES["piscem"]
BINARIES["kallisto-forward"] = BINARIES["kallisto"]


def run(tool, corpus, repeat):
    prefix = RAW / f"{corpus}-{tool}-{repeat}"
    out = DATA / f"{corpus}-{tool}"
    out.mkdir(exist_ok=True)
    if corpus in ("pbmc10m", "pbmc10m-raw"):
        suffix = ".fastq" if corpus.endswith("-raw") else ".fastq.gz"
        r1, r2 = [OLD / f"pbmc10m_R{i}{suffix}" for i in (1,2)]
        count = 10_000_000
    elif corpus == "smoke":
        r1, r2 = [ROOT / f"benchdata/pbmc_10k_v3/pbmc_10k_v3_100k_L001_R{i}.fastq.gz" for i in (1,2)]
        count = 100_000
    else:
        r1, r2 = [DATA / f"{corpus}_R{i}.fastq.gz" for i in (1,2)]
        count = 420_000 if corpus in ("synthetic", "heldout", "heldout150") else 1_000_000
    cmd = [str(BINARIES[tool])]
    if tool.startswith("constellation"):
        cmd += ["map", "--index", str(OLD / "reference.ptrhash.ecidx"),
                "--r1", str(r1), "--r2", str(r2), "--batch-size", "65536",
                "--output-format", "gene-ec-rad", "--out", str(out / "map.cstrad"),
                "--emit-metrics", str(prefix) + ".metrics.json"]
        if tool.endswith("consensus"):
            cmd += ["--gene-ec-consensus"]
    elif tool.startswith("piscem"):
        cmd += ["map-sc", "-i", str(DATA / "piscem-splici"),
                "-g", "chromium_v3", "-t", "12", "-1", str(r1), "-2", str(r2),
                "-o", str(out)]
    else:
        cmd += ["bus", "-i", str(DATA / "kallisto-splici.idx"),
                "-x", "10xv3", "--fr-stranded" if tool.endswith("forward") else "--unstranded", "-t", "12", "-o", str(out),
                *(["-n"] if corpus in ("synthetic", "heldout", "heldout150", "real1m", "lane2") else []), str(r1), str(r2)]
    timed = ["/usr/bin/time", "-v", "-o", str(prefix) + ".time", "taskset", "-c", "0-11", *cmd]
    env = dict(os.environ, RAYON_NUM_THREADS="12", OMP_NUM_THREADS="12")
    start = time.perf_counter()
    events = []
    with Path(str(prefix) + ".log").open("wb") as log:
        with subprocess.Popen(timed, cwd=ROOT, env=env, stdout=subprocess.PIPE,
                              stderr=subprocess.STDOUT) as process:
            while chunk := os.read(process.stdout.fileno(), 65536):
                events.append(dict(seconds=time.perf_counter()-start,text=chunk.decode(errors="replace")))
                log.write(chunk)
            returncode = process.wait()
    elapsed = time.perf_counter() - start
    metadata = dict(tool=tool, corpus=corpus, repeat=repeat, command=timed,
                    environment={k: env[k] for k in ("RAYON_NUM_THREADS", "OMP_NUM_THREADS")},
                    process_wall_seconds=elapsed, pairs=count,
                    pairs_per_second=count/elapsed, exit_code=returncode)
    Path(str(prefix)+".events.json").write_text(json.dumps(events,indent=2)+"\n")
    if tool.startswith("kallisto"):
        # std::cerr is unbuffered. The release's ProcessBUSReads prints these
        # markers immediately around MP.processReads(); retain raw timestamps.
        start_phase = end_phase = None
        accumulated = ""
        for event in events:
            accumulated += event["text"]
            marker = "[quant] finding pseudoalignments for the reads ..."
            if start_phase is None and marker in accumulated:
                start_phase = event["seconds"]
            if start_phase is not None and end_phase is None and " done\n" in accumulated.split(marker,1)[1]:
                end_phase = event["seconds"]
        if start_phase is not None and end_phase is not None:
            metadata["pseudoalignment_loop_seconds"] = end_phase-start_phase
            metadata["before_pseudoalignment_seconds"] = start_phase
            metadata["after_pseudoalignment_seconds"] = elapsed-end_phase
    with BINARIES[tool].open("rb") as handle:
        metadata["executable_sha256"] = hashlib.file_digest(handle,"sha256").hexdigest()
    for candidate in ("map_info.json", "run_info.json"):
        if (out / candidate).exists():
            shutil.copyfile(out / candidate, str(prefix) + "." + candidate)
    Path(str(prefix) + ".run.json").write_text(json.dumps(metadata, indent=2) + "\n")
    print(json.dumps({k:metadata[k] for k in ("tool","corpus","repeat","process_wall_seconds","pairs_per_second","exit_code")}), flush=True)
    if returncode:
        raise RuntimeError(f"{tool} failed: see {prefix}.log")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--corpus", required=True,
                        choices=["smoke","pbmc10m","pbmc10m-raw","synthetic","real1m"])
    parser.add_argument("--tools", nargs="+", default=list(BINARIES), choices=list(BINARIES))
    parser.add_argument("--repeats", type=int, default=0,
                        help="0 is a single accuracy/smoke run; positive includes warmup")
    args = parser.parse_args()
    with (OLD / "benchmark.lock").open("w") as lock:
        fcntl.flock(lock, fcntl.LOCK_EX)
        # Interleave tool order across repeats to avoid assigning all drift to one tool.
        rounds = ["warmup", *range(1,args.repeats+1)] if args.repeats else ["validation"]
        for index, repeat in enumerate(rounds):
            order = args.tools[index % len(args.tools):] + args.tools[:index % len(args.tools)]
            for tool in order:
                run(tool, args.corpus, repeat)

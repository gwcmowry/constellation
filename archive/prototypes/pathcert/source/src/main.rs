use std::env;
use std::fs::File;
use std::hint::black_box;
use std::io::{self, BufWriter, Write};
use std::path::PathBuf;
use std::time::Instant;

use pathcert::graph::{
    components_from_genome_gtf, components_from_transcript_fasta, Component, TranscriptPath,
};
use pathcert::index::PathCertIndex;
use pathcert::io::fastq::read_fastq;
use pathcert::io::fastq::FastqRecord;
use pathcert::mapper::{
    map_records, map_records_component_bucketed, MappingConfig, MappingStatus, MappingSummary,
};

fn main() {
    if let Err(error) = run() {
        eprintln!("error: {error}");
        std::process::exit(1);
    }
}

fn run() -> io::Result<()> {
    let mut args = env::args().skip(1);
    let Some(command) = args.next() else {
        print_help();
        return Ok(());
    };
    let args: Vec<String> = args.collect();

    match command.as_str() {
        "index" => index_cmd(&args),
        "map" => map_cmd(&args),
        "inspect" => inspect_cmd(&args),
        "bench" => bench_cmd(&args),
        "bench-fastq" => bench_fastq_cmd(&args),
        "help" | "--help" | "-h" => {
            print_help();
            Ok(())
        }
        other => Err(invalid_input(format!("unknown command {other}"))),
    }
}

fn index_cmd(args: &[String]) -> io::Result<()> {
    let transcripts = opt_path(args, "--transcripts");
    let genome = opt_path(args, "--genome");
    let gtf = opt_path(args, "--gtf");
    let out = required_path(args, "--out")?;
    let k = opt_value(args, "--k")
        .unwrap_or("21")
        .parse::<usize>()
        .map_err(invalid_input)?;

    let components = if let Some(transcripts) = transcripts {
        components_from_transcript_fasta(transcripts)?
    } else {
        let genome = genome.ok_or_else(|| invalid_input("missing --genome or --transcripts"))?;
        let gtf = gtf.ok_or_else(|| invalid_input("missing --gtf"))?;
        components_from_genome_gtf(genome, gtf)?
    };

    let index = PathCertIndex::new(k, components)?;
    index.write_text(out)?;
    eprintln!(
        "wrote PathCert index: components={} seed_postings={} k={}",
        index.components.len(),
        index.seed_index.len(),
        index.k
    );
    Ok(())
}

fn map_cmd(args: &[String]) -> io::Result<()> {
    let index_path = required_path(args, "--index")?;
    let reads_path = required_path(args, "--reads")?;
    let out_path = required_path(args, "--out")?;
    let summary_path = opt_path(args, "--summary");
    let max_postings_per_seed = opt_value(args, "--max-postings-per-seed")
        .unwrap_or("256")
        .parse::<usize>()
        .map_err(invalid_input)?;

    let index = PathCertIndex::read_text(index_path)?;
    let reads = read_fastq(reads_path)?;
    let config = MappingConfig {
        max_postings_per_seed,
        ..MappingConfig::default()
    };
    let (mappings, certificates, summary) = map_records(&index, &reads, &config);

    let mut writer = BufWriter::new(File::create(out_path)?);
    writeln!(
        writer,
        "read_id\tstatus\tecert_id\tcomponent_id\tgene_id\tpathset\ttranscripts\tflags"
    )?;
    for mapping in mappings {
        if mapping.status == MappingStatus::Unmapped {
            writeln!(writer, "{}\tunmapped\t.\t.\t.\t.\t.\t0", mapping.read_id)?;
            continue;
        }
        for cert_id in mapping.certificate_ids {
            let cert = &certificates.certificates()[cert_id as usize];
            let component = index.component(cert.component_id).ok_or_else(|| {
                invalid_input(format!(
                    "certificate references component {}",
                    cert.component_id
                ))
            })?;
            let names = component.transcript_names_for_mask(cert.pathset).join(",");
            writeln!(
                writer,
                "{}\t{:?}\t{}\t{}\t{}\t{}\t{}\t{}",
                mapping.read_id,
                mapping.status,
                cert.id,
                cert.component_id,
                component.gene_id,
                cert.pathset,
                names,
                cert.flags
            )?;
        }
    }

    if let Some(path) = summary_path {
        let mut writer = BufWriter::new(File::create(path)?);
        writeln!(writer, "metric\tvalue")?;
        writeln!(writer, "reads\t{}", summary.reads)?;
        writeln!(writer, "mapped\t{}", summary.mapped)?;
        writeln!(writer, "multi_component\t{}", summary.multi_component)?;
        writeln!(writer, "unmapped\t{}", summary.unmapped)?;
        writeln!(writer, "seed_lookups\t{}", summary.seed_lookups)?;
        writeln!(
            writer,
            "candidate_components\t{}",
            summary.candidate_components
        )?;
        writeln!(
            writer,
            "pathset_intersections\t{}",
            summary.pathset_intersections
        )?;
        writeln!(
            writer,
            "certificates\t{}",
            certificates.certificates().len()
        )?;
    }

    Ok(())
}

fn bench_fastq_cmd(args: &[String]) -> io::Result<()> {
    let index_path = required_path(args, "--index")?;
    let reads_path = required_path(args, "--reads")?;
    let max_postings_per_seed = parse_opt(args, "--max-postings-per-seed", 256usize)?;

    let index = PathCertIndex::read_text(index_path)?;
    let reads = read_fastq(reads_path)?;
    let config = MappingConfig {
        max_postings_per_seed,
        ..MappingConfig::default()
    };

    println!(
        "dataset\treads\t{}\tk\t{}\tcomponents\t{}\tseed_postings\t{}",
        reads.len(),
        index.k,
        index.components.len(),
        index.seed_index.len()
    );
    println!("mode\tseconds\treads_per_sec\tmapped\tmulti_component\tunmapped\tseed_lookups\tcandidate_components\tpathset_intersections\tcertificates");

    let start = Instant::now();
    let (_mappings, certs, summary) = map_records(&index, black_box(&reads), &config);
    let seconds = start.elapsed().as_secs_f64();
    print_bench_row(
        "input_order",
        seconds,
        reads.len(),
        &summary,
        certs.certificates().len(),
    );

    let start = Instant::now();
    let (_mappings, certs, summary) =
        map_records_component_bucketed(&index, black_box(&reads), &config);
    let seconds = start.elapsed().as_secs_f64();
    print_bench_row(
        "component_bucketed",
        seconds,
        reads.len(),
        &summary,
        certs.certificates().len(),
    );

    Ok(())
}

fn bench_cmd(args: &[String]) -> io::Result<()> {
    let component_count = parse_opt(args, "--components", 256usize)?;
    let transcripts_per_component = parse_opt(args, "--transcripts-per-component", 2usize)?;
    let read_count = parse_opt(args, "--reads", 200_000usize)?;
    let read_len = parse_opt(args, "--read-len", 48usize)?;
    let k = parse_opt(args, "--k", 15usize)?;
    let seq_len = parse_opt(args, "--seq-len", 192usize)?;

    if transcripts_per_component == 0 || transcripts_per_component > 64 {
        return Err(invalid_input("--transcripts-per-component must be 1..64"));
    }
    if read_len == 0 || read_len > seq_len {
        return Err(invalid_input("--read-len must be 1..seq-len"));
    }

    let components = synthetic_components(component_count, transcripts_per_component, seq_len);
    let reads = synthetic_reads(&components, read_count, read_len);
    let index = PathCertIndex::new(k, components)?;
    let config = MappingConfig::default();

    println!(
        "dataset\tcomponents\t{}\ttranscripts_per_component\t{}\treads\t{}\tread_len\t{}\tk\t{}\tseed_postings\t{}",
        component_count,
        transcripts_per_component,
        read_count,
        read_len,
        k,
        index.seed_index.len()
    );
    println!("mode\tseconds\treads_per_sec\tmapped\tmulti_component\tunmapped\tseed_lookups\tcandidate_components\tpathset_intersections\tcertificates");

    let start = Instant::now();
    let (_mappings, certs, summary) = map_records(&index, black_box(&reads), &config);
    let seconds = start.elapsed().as_secs_f64();
    print_bench_row(
        "input_order",
        seconds,
        read_count,
        &summary,
        certs.certificates().len(),
    );

    let start = Instant::now();
    let (_mappings, certs, summary) =
        map_records_component_bucketed(&index, black_box(&reads), &config);
    let seconds = start.elapsed().as_secs_f64();
    print_bench_row(
        "component_bucketed",
        seconds,
        read_count,
        &summary,
        certs.certificates().len(),
    );

    Ok(())
}

fn print_bench_row(
    mode: &str,
    seconds: f64,
    read_count: usize,
    summary: &MappingSummary,
    certificate_count: usize,
) {
    println!(
        "{}\t{:.6}\t{:.0}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
        mode,
        seconds,
        read_count as f64 / seconds.max(f64::MIN_POSITIVE),
        summary.mapped,
        summary.multi_component,
        summary.unmapped,
        summary.seed_lookups,
        summary.candidate_components,
        summary.pathset_intersections,
        certificate_count
    );
}

fn synthetic_components(
    component_count: usize,
    transcripts_per_component: usize,
    seq_len: usize,
) -> Vec<Component> {
    let mut components = Vec::with_capacity(component_count);
    for component_id in 0..component_count {
        let shared_prefix = synthetic_seq(component_id as u64, 0, seq_len / 4);
        let shared_suffix = synthetic_seq(component_id as u64, 9_999, seq_len / 4);
        let unique_len = seq_len.saturating_sub(shared_prefix.len() + shared_suffix.len());
        let mut transcripts = Vec::with_capacity(transcripts_per_component);
        for tx_idx in 0..transcripts_per_component {
            let mut sequence = Vec::with_capacity(seq_len);
            sequence.extend_from_slice(&shared_prefix);
            sequence.extend_from_slice(&synthetic_seq(
                component_id as u64,
                tx_idx as u64 + 1,
                unique_len,
            ));
            sequence.extend_from_slice(&shared_suffix);
            transcripts.push(TranscriptPath {
                name: format!("TX{component_id}_{tx_idx}"),
                sequence,
            });
        }
        components.push(Component {
            id: component_id as u32,
            gene_id: format!("GENE{component_id}"),
            transcripts,
        });
    }
    components
}

fn synthetic_reads(
    components: &[Component],
    read_count: usize,
    read_len: usize,
) -> Vec<FastqRecord> {
    let mut reads = Vec::with_capacity(read_count);
    for read_idx in 0..read_count {
        let component = &components[read_idx % components.len()];
        let transcript =
            &component.transcripts[(read_idx / components.len()) % component.transcripts.len()];
        let max_start = transcript.sequence.len() - read_len;
        let start = if max_start == 0 {
            0
        } else {
            lcg(read_idx as u64 + 17) as usize % (max_start + 1)
        };
        let seq = transcript.sequence[start..start + read_len].to_vec();
        reads.push(FastqRecord {
            id: format!("read{read_idx}"),
            qual: vec![b'I'; seq.len()],
            seq,
        });
    }
    reads
}

fn synthetic_seq(component_seed: u64, transcript_seed: u64, len: usize) -> Vec<u8> {
    let mut state = component_seed
        .wrapping_mul(0x9E37_79B9_7F4A_7C15)
        .wrapping_add(transcript_seed.wrapping_mul(0xD1B5_4A32_D192_ED03))
        .wrapping_add(0xA076_1D64_78BD_642F);
    let alphabet = *b"ACGT";
    let mut seq = Vec::with_capacity(len);
    for _ in 0..len {
        state = lcg(state);
        seq.push(alphabet[(state as usize >> 32) & 3]);
    }
    seq
}

fn lcg(value: u64) -> u64 {
    value
        .wrapping_mul(6_364_136_223_846_793_005)
        .wrapping_add(1_442_695_040_888_963_407)
}

fn inspect_cmd(args: &[String]) -> io::Result<()> {
    let index_path = required_path(args, "--index")?;
    let component_filter = opt_value(args, "--component");
    let index = PathCertIndex::read_text(index_path)?;
    println!("k\t{}", index.k);
    println!("components\t{}", index.components.len());
    println!("seed_postings\t{}", index.seed_index.len());
    for component in &index.components {
        if let Some(filter) = component_filter {
            if component.gene_id != filter && component.id.to_string() != filter {
                continue;
            }
        }
        println!(
            "component\t{}\t{}\ttranscripts={}",
            component.id,
            component.gene_id,
            component.transcripts.len()
        );
        for (idx, transcript) in component.transcripts.iter().enumerate() {
            println!(
                "transcript\t{}\t{}\tlen={}",
                idx,
                transcript.name,
                transcript.sequence.len()
            );
        }
    }
    Ok(())
}

fn opt_value<'a>(args: &'a [String], name: &str) -> Option<&'a str> {
    args.windows(2)
        .find(|window| window[0] == name)
        .map(|window| window[1].as_str())
}

fn opt_path(args: &[String], name: &str) -> Option<PathBuf> {
    opt_value(args, name).map(PathBuf::from)
}

fn parse_opt<T>(args: &[String], name: &str, default: T) -> io::Result<T>
where
    T: std::str::FromStr,
    T::Err: ToString,
{
    opt_value(args, name)
        .map(|value| value.parse::<T>().map_err(invalid_input))
        .unwrap_or(Ok(default))
}

fn required_path(args: &[String], name: &str) -> io::Result<PathBuf> {
    opt_path(args, name).ok_or_else(|| invalid_input(format!("missing {name}")))
}

fn invalid_input(error: impl ToString) -> io::Error {
    io::Error::new(io::ErrorKind::InvalidInput, error.to_string())
}

fn print_help() {
    eprintln!(
        "PathCert commands:\n  index --transcripts transcripts.fa --k 21 --out index.pcidx\n  index --genome genome.fa --gtf annotation.gtf --k 21 --out index.pcidx\n  map --index index.pcidx --reads reads.fq --out ec.tsv [--summary summary.tsv]\n  inspect --index index.pcidx [--component GENE]\n  bench [--components 256] [--reads 200000] [--k 15]\n  bench-fastq --index index.pcidx --reads reads.fq"
    );
}

use std::fs;
use std::path::{Path, PathBuf};
use std::process::Command;

use pathcert::graph::{
    components_from_genome_gtf, components_from_transcript_fasta, toy_exon_skipping_component,
};
use pathcert::index::PathCertIndex;
use pathcert::io::fastq::FastqRecord;
use pathcert::mapper::{map_records, MappingConfig, MappingStatus};

#[test]
fn toy_reads_emit_ambiguous_and_distinguishing_certificates() {
    let component = toy_exon_skipping_component();
    let index = PathCertIndex::new(5, vec![component]).unwrap();
    let records = vec![
        fastq("shared_exon1", "AAAACCCCGGGG"),
        fastq("include_junction", "GGGGTTTTGGGG"),
        fastq("skip_junction", "GGGGACACACAC"),
    ];

    let (mappings, certificates, _summary) =
        map_records(&index, &records, &MappingConfig::default());

    assert_eq!(mappings[0].status, MappingStatus::Mapped);
    assert_eq!(mappings[1].status, MappingStatus::Mapped);
    assert_eq!(mappings[2].status, MappingStatus::Mapped);

    let shared = &certificates.certificates()[mappings[0].certificate_ids[0] as usize];
    let include = &certificates.certificates()[mappings[1].certificate_ids[0] as usize];
    let skip = &certificates.certificates()[mappings[2].certificate_ids[0] as usize];
    assert_eq!(shared.pathset, 0b11);
    assert_eq!(include.pathset, 0b01);
    assert_eq!(skip.pathset, 0b10);
}

#[test]
fn exact_simulated_reads_include_source_transcript() {
    let component = toy_exon_skipping_component();
    let index = PathCertIndex::new(5, vec![component.clone()]).unwrap();
    let mut records = Vec::new();
    let mut expected = Vec::new();

    for (tx_idx, transcript) in component.transcripts.iter().enumerate() {
        for (read_idx, window) in transcript.sequence.windows(12).step_by(3).enumerate() {
            records.push(fastq(
                &format!("{}_{}", transcript.name, read_idx),
                std::str::from_utf8(window).unwrap(),
            ));
            expected.push(1u64 << tx_idx);
        }
    }

    let (mappings, certificates, _summary) =
        map_records(&index, &records, &MappingConfig::default());
    for (mapping, expected_bit) in mappings.iter().zip(expected) {
        assert_ne!(
            mapping.status,
            MappingStatus::Unmapped,
            "{}",
            mapping.read_id
        );
        let mut includes_source = false;
        for cert_id in &mapping.certificate_ids {
            let cert = &certificates.certificates()[*cert_id as usize];
            includes_source |= cert.pathset & expected_bit != 0;
        }
        assert!(
            includes_source,
            "source transcript excluded for {}",
            mapping.read_id
        );
    }
}

#[test]
fn reuses_existing_tiny_transcriptome_fixture() {
    let components = components_from_transcript_fasta(concat!(
        env!("CARGO_MANIFEST_DIR"),
        "/../../../../tests/tiny_transcriptome.fa"
    ))
    .unwrap();
    let index = PathCertIndex::new(7, components).unwrap();
    let reads = vec![fastq("fixture_read0_cdna", "ACGTACGTACGTACGTACGT")];

    let (mappings, certificates, summary) = map_records(&index, &reads, &MappingConfig::default());
    assert_eq!(summary.reads, 1);
    assert_ne!(mappings[0].status, MappingStatus::Unmapped);

    let mut genes = Vec::new();
    for cert_id in &mappings[0].certificate_ids {
        let cert = &certificates.certificates()[*cert_id as usize];
        genes.push(index.component(cert.component_id).unwrap().gene_id.clone());
    }
    assert!(genes.iter().any(|gene| gene == "GENE_A"));
}

#[test]
fn gtf_genome_builder_constructs_transcript_sequences() {
    let dir = temp_dir("pathcert_gtf");
    let genome = dir.join("genome.fa");
    let gtf = dir.join("genes.gtf");
    fs::write(&genome, ">chr1\nAAAACCCCGGGGTTTT\n").unwrap();
    fs::write(
        &gtf,
        concat!(
            "chr1\ttest\texon\t1\t4\t.\t+\t.\tgene_id \"G1\"; transcript_id \"T1\";\n",
            "chr1\ttest\texon\t9\t12\t.\t+\t.\tgene_id \"G1\"; transcript_id \"T1\";\n",
            "chr1\ttest\texon\t1\t4\t.\t+\t.\tgene_id \"G1\"; transcript_id \"T2\";\n",
            "chr1\ttest\texon\t13\t16\t.\t+\t.\tgene_id \"G1\"; transcript_id \"T2\";\n",
        ),
    )
    .unwrap();

    let components = components_from_genome_gtf(&genome, &gtf).unwrap();
    assert_eq!(components.len(), 1);
    assert_eq!(components[0].gene_id, "G1");
    assert_eq!(components[0].transcripts.len(), 2);
    let sequences: Vec<String> = components[0]
        .transcripts
        .iter()
        .map(|tx| String::from_utf8(tx.sequence.clone()).unwrap())
        .collect();
    assert!(sequences.contains(&"AAAAGGGG".to_string()));
    assert!(sequences.contains(&"AAAATTTT".to_string()));
}

#[test]
fn cli_indexes_maps_and_writes_summary() {
    let dir = temp_dir("pathcert_cli");
    let index = dir.join("tiny.pcidx");
    let out = dir.join("tiny.ec.tsv");
    let summary = dir.join("tiny.summary.tsv");

    assert_success(
        pathcert_bin()
            .args([
                "index",
                "--transcripts",
                concat!(
                    env!("CARGO_MANIFEST_DIR"),
                    "/../../../../tests/tiny_transcriptome.fa"
                ),
                "--k",
                "7",
                "--out",
                path(&index),
            ])
            .output()
            .unwrap(),
    );
    assert_success(
        pathcert_bin()
            .args([
                "map",
                "--index",
                path(&index),
                "--reads",
                concat!(
                    env!("CARGO_MANIFEST_DIR"),
                    "/../../../../tests/tiny_reads_R2.fastq"
                ),
                "--out",
                path(&out),
                "--summary",
                path(&summary),
            ])
            .output()
            .unwrap(),
    );

    let rows = fs::read_to_string(out).unwrap();
    assert!(rows.starts_with("read_id\tstatus\tecert_id"));
    assert!(rows.contains("GENE_A"));
    let metrics = fs::read_to_string(summary).unwrap();
    assert!(metrics.contains("reads\t1"));
    assert!(metrics.contains("certificates\t"));
}

fn fastq(id: &str, seq: &str) -> FastqRecord {
    FastqRecord {
        id: id.to_string(),
        seq: seq.as_bytes().to_vec(),
        qual: vec![b'I'; seq.len()],
    }
}

fn path(path: &Path) -> &str {
    path.to_str().unwrap()
}

fn pathcert_bin() -> Command {
    Command::new(env!("CARGO_BIN_EXE_pathcert"))
}

fn assert_success(output: std::process::Output) {
    assert!(
        output.status.success(),
        "status: {:?}\nstdout:\n{}\nstderr:\n{}",
        output.status.code(),
        String::from_utf8_lossy(&output.stdout),
        String::from_utf8_lossy(&output.stderr)
    );
}

fn temp_dir(name: &str) -> PathBuf {
    let mut dir = std::env::temp_dir();
    dir.push(format!(
        "{}_{}_{}",
        name,
        std::process::id(),
        std::thread::current().name().unwrap_or("test")
    ));
    let _ = fs::remove_dir_all(&dir);
    fs::create_dir_all(&dir).unwrap();
    dir
}

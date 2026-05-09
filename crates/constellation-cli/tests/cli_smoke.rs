use std::fs;
use std::path::Path;
use std::process::Command;

fn constellation() -> Command {
    Command::new(env!("CARGO_BIN_EXE_constellation"))
}

#[test]
fn cli_index_simulate_gene_ec_map_and_report() {
    let tmp = tempfile::tempdir().unwrap();
    let fasta = tmp.path().join("tiny.fa");
    fs::write(
        &fasta,
        concat!(
            ">A1|gene=GENE_A\n",
            "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT\n",
            ">B1|gene=GENE_B\n",
            "TGCATGCATGCATGCATGCATGCATGCATGCATGCATGCA\n",
        ),
    )
    .unwrap();
    let index = tmp.path().join("tiny.cbidx");
    let ec_index = tmp.path().join("tiny.ecidx");

    assert_success(
        constellation()
            .args([
                "index",
                "--transcripts",
                path(&fasta),
                "--k",
                "11",
                "--out",
                path(&index),
            ])
            .output()
            .unwrap(),
    );

    let inspect = constellation()
        .args(["inspect-index", "--index", path(&index)])
        .output()
        .unwrap();
    assert_success(inspect);

    assert_success(
        constellation()
            .args([
                "index-ec",
                "--transcripts",
                path(&fasta),
                "--k",
                "11",
                "--out",
                path(&ec_index),
            ])
            .output()
            .unwrap(),
    );

    let prefix = tmp.path().join("sim");
    assert_success(
        constellation()
            .args([
                "simulate",
                "--transcripts",
                path(&fasta),
                "--num-reads",
                "12",
                "--read-len",
                "20",
                "--scenario",
                "exact",
                "--out-prefix",
                path(&prefix),
            ])
            .output()
            .unwrap(),
    );

    let assignments = tmp.path().join("gene-ec.assignments.tsv");
    let metrics = tmp.path().join("gene-ec.metrics.json");
    let diagnostics = tmp.path().join("gene-ec.unmapped.tsv");
    assert_success(
        constellation()
            .args([
                "map",
                "--index",
                path(&ec_index),
                "--r1",
                path(&tmp.path().join("sim_R1.fastq")),
                "--r2",
                path(&tmp.path().join("sim_R2.fastq")),
                "--mode",
                "gene-ec",
                "--output-format",
                "tsv",
                "--emit-metrics",
                path(&metrics),
                "--emit-unmapped-diagnostics",
                path(&diagnostics),
                "--out",
                path(&assignments),
            ])
            .output()
            .unwrap(),
    );
    let metrics_text = fs::read_to_string(&metrics).unwrap();
    assert!(metrics_text.contains("\"mode\": \"gene-ec\""));
    assert!(metrics_text.contains("\"same_gene_multitranscript_rate\""));
    assert!(metrics_text.contains("\"gene_countable_rate\""));
    let diagnostics_text = fs::read_to_string(&diagnostics).unwrap();
    assert!(diagnostics_text.starts_with("read_id\treason\tseq_len\t"));

    let report = constellation()
        .args([
            "bench-report",
            "--assignments",
            path(&assignments),
            "--truth",
            path(&tmp.path().join("sim_truth.tsv")),
            "--metrics",
            path(&metrics),
        ])
        .output()
        .unwrap();
    assert_success(report);

    let count_prefix = tmp.path().join("counts");
    assert_success(
        constellation()
            .args([
                "count",
                "--assignments",
                path(&assignments),
                "--index",
                path(&index),
                "--out-prefix",
                path(&count_prefix),
            ])
            .output()
            .unwrap(),
    );
    assert!(tmp.path().join("counts_matrix.mtx").exists());
    assert!(tmp.path().join("counts_barcodes.tsv").exists());
    assert!(tmp.path().join("counts_features.tsv").exists());
}

fn path(path: &Path) -> &str {
    path.to_str().unwrap()
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

#[test]
fn cli_golden_assignment_types() {
    let tmp = tempfile::tempdir().unwrap();
    let fasta = tmp.path().join("golden.fa");
    fs::write(
        &fasta,
        concat!(
            ">A1|gene=GENE_A\n",
            "ACGTACGTACGTACGT\n",
            ">B1|gene=GENE_B\n",
            "TTGCGATCGATCGGTA\n",
            ">C1|gene=GENE_C\n",
            "ACGTACGTACGTACGT\n",
        ),
    )
    .unwrap();
    let index = tmp.path().join("golden.cbidx");
    let ec_index = tmp.path().join("golden.ecidx");
    assert_success(
        constellation()
            .args([
                "index",
                "--transcripts",
                path(&fasta),
                "--k",
                "5",
                "--out",
                path(&index),
            ])
            .output()
            .unwrap(),
    );
    assert_success(
        constellation()
            .args([
                "index-ec",
                "--transcripts",
                path(&fasta),
                "--k",
                "5",
                "--out",
                path(&ec_index),
            ])
            .output()
            .unwrap(),
    );

    let r1 = tmp.path().join("R1.fastq");
    let r2 = tmp.path().join("R2.fastq");
    let r1_seq = "ACGTACGTACGTACGTTTTTCCCCAAAA";
    fs::write(
        &r1,
        [
            ("ambiguous", r1_seq, "IIIIIIIIIIIIIIIIIIIIIIIIIIII"),
            ("unique", r1_seq, "IIIIIIIIIIIIIIIIIIIIIIIIIIII"),
            ("unmapped", r1_seq, "IIIIIIIIIIIIIIIIIIIIIIIIIIII"),
            ("low_complexity", r1_seq, "IIIIIIIIIIIIIIIIIIIIIIIIIIII"),
            ("low_quality", r1_seq, "IIIIIIIIIIIIIIIIIIIIIIIIIIII"),
        ]
        .into_iter()
        .map(|(id, seq, qual)| format!("@{id}\n{seq}\n+\n{qual}\n"))
        .collect::<String>(),
    )
    .unwrap();
    fs::write(
        &r2,
        [
            ("ambiguous", "ACGTACGTACGT", "IIIIIIIIIIII"),
            ("unique", "TTGCGATCGATC", "IIIIIIIIIIII"),
            ("unmapped", "GATACGCTAGTC", "IIIIIIIIIIII"),
            ("low_complexity", "AAAAAAAAAAAA", "IIIIIIIIIIII"),
            ("low_quality", "TTGCGATCGATC", "!!!!!!!!!!!!"),
        ]
        .into_iter()
        .map(|(id, seq, qual)| format!("@{id}\n{seq}\n+\n{qual}\n"))
        .collect::<String>(),
    )
    .unwrap();

    let assignments = tmp.path().join("assignments.tsv");
    assert_success(
        constellation()
            .args([
                "map",
                "--index",
                path(&ec_index),
                "--r1",
                path(&r1),
                "--r2",
                path(&r2),
                "--max-postings-per-seed",
                "256",
                "--output-format",
                "tsv",
                "--out",
                path(&assignments),
            ])
            .output()
            .unwrap(),
    );
    let rows = fs::read_to_string(assignments).unwrap();
    assert!(rows.contains("\tambiguous_gene\t"));
    assert!(rows.contains("\tunique_gene\t"));
    assert!(rows.contains("\tunmapped\t"));
    assert!(rows.contains("\tlow_complexity\t"));
    assert!(rows.contains("\tlow_quality\t"));

    let manual_assignments = tmp.path().join("manual_assignments.tsv");
    fs::write(
        &manual_assignments,
        concat!(
            "read_id\tcell_barcode\tumi\tassignment_type\tgene_id\ttranscript_id\tcandidate_count\tscore\tflags\n",
            "0\tCB\tUMI1\tunique_gene\t0\t0\t1\t12\t0\n",
            "1\tCB\tUMI2\tambiguous_transcript_same_gene\t0\t.\t2\t12\t0\n",
            "2\tCB\tUMI3\tambiguous_gene\t.\t.\t2\t12\t0\n",
        ),
    )
    .unwrap();
    let count_prefix = tmp.path().join("manual_counts");
    assert_success(
        constellation()
            .args([
                "count",
                "--assignments",
                path(&manual_assignments),
                "--index",
                path(&index),
                "--out-prefix",
                path(&count_prefix),
            ])
            .output()
            .unwrap(),
    );
    let matrix = fs::read_to_string(tmp.path().join("manual_counts_matrix.mtx")).unwrap();
    assert!(matrix.contains("\n1 1 2\n"));

    let corrected_assignments = tmp.path().join("corrected_assignments.tsv");
    fs::write(
        &corrected_assignments,
        concat!(
            "read_id\tcell_barcode\tumi\tassignment_type\tgene_id\ttranscript_id\tcandidate_count\tscore\tflags\n",
            "0\tAAAAAAAAAAAAAAAA\tAAAAAAAAAAAA\tunique_gene\t0\t0\t1\t12\t0\n",
            "1\tAAAAAAAAAAAAAAAA\tAAAAAAAAAAAA\tunique_gene\t0\t0\t1\t12\t0\n",
            "2\tCAAAAAAAAAAAAAAA\tAAAAAAAAAAAT\tunique_gene\t0\t0\t1\t12\t0\n",
            "3\tAAAAAAAAAAAAAAAA\tCCCCCCCCCCCC\tunique_gene\t0\t0\t1\t12\t0\n",
        ),
    )
    .unwrap();
    let whitelist = tmp.path().join("barcodes.tsv");
    fs::write(&whitelist, "AAAAAAAAAAAAAAAA\n").unwrap();
    let corrected_prefix = tmp.path().join("corrected_counts");
    let corrected_metrics = tmp.path().join("corrected_counts.metrics.json");
    assert_success(
        constellation()
            .args([
                "count",
                "--assignments",
                path(&corrected_assignments),
                "--index",
                path(&index),
                "--barcode-whitelist",
                path(&whitelist),
                "--emit-metrics",
                path(&corrected_metrics),
                "--out-prefix",
                path(&corrected_prefix),
            ])
            .output()
            .unwrap(),
    );
    let corrected_matrix =
        fs::read_to_string(tmp.path().join("corrected_counts_matrix.mtx")).unwrap();
    assert!(corrected_matrix.contains("\n1 1 2\n"));
    let corrected_barcodes =
        fs::read_to_string(tmp.path().join("corrected_counts_barcodes.tsv")).unwrap();
    assert_eq!(corrected_barcodes, "AAAAAAAAAAAAAAAA\n");
    let corrected_metrics = fs::read_to_string(corrected_metrics).unwrap();
    assert!(corrected_metrics.contains("\"barcode_corrected\": 1"));
}

#[test]
fn cli_build_transcriptome_target_kinds() {
    let tmp = tempfile::tempdir().unwrap();
    let genome = tmp.path().join("genome.fa");
    fs::write(&genome, ">1\nAACCGGTTAACCGGTT\n").unwrap();
    let gtf = tmp.path().join("genes.gtf");
    fs::write(
        &gtf,
        concat!(
            "1\ttest\tgene\t2\t15\t.\t+\t.\tgene_id \"GENE1\";\n",
            "1\ttest\texon\t2\t4\t.\t+\t.\tgene_id \"GENE1\"; transcript_id \"TX1\";\n",
            "1\ttest\texon\t9\t12\t.\t+\t.\tgene_id \"GENE1\"; transcript_id \"TX1\";\n",
        ),
    )
    .unwrap();

    for (kind, target_marker) in [
        ("exon-transcripts", "target:exon_transcript"),
        ("gene-bodies", "target:gene_body"),
        ("introns-only", "target:intron"),
        ("intron-flanks", "target:intron"),
        ("exon-plus-gene-body", "target:gene_body"),
    ] {
        let out = tmp.path().join(format!("{kind}.fa"));
        assert_success(
            constellation()
                .args([
                    "build-transcriptome-target",
                    "--genome",
                    path(&genome),
                    "--gtf",
                    path(&gtf),
                    "--target-kind",
                    kind,
                    "--out",
                    path(&out),
                ])
                .output()
                .unwrap(),
        );
        let text = fs::read_to_string(out).unwrap();
        assert!(text.contains(target_marker));
    }
}

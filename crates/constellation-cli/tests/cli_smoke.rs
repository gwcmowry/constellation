use std::fs;
use std::path::Path;
use std::process::Command;

fn constellation() -> Command {
    Command::new(env!("CARGO_BIN_EXE_constellation"))
}

#[test]
fn cli_index_simulate_map_modes_and_report() {
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

    for mode in ["read-at-a-time", "sketch-bucket", "candidate-locus"] {
        let assignments = tmp.path().join(format!("{mode}.assignments.tsv"));
        let metrics = tmp.path().join(format!("{mode}.metrics.json"));
        assert_success(
            constellation()
                .args([
                    "map",
                    "--index",
                    path(&index),
                    "--r1",
                    path(&tmp.path().join("sim_R1.fastq")),
                    "--r2",
                    path(&tmp.path().join("sim_R2.fastq")),
                    "--mode",
                    mode,
                    "--out",
                    path(&assignments),
                    "--emit-metrics",
                    path(&metrics),
                ])
                .output()
                .unwrap(),
        );
        let metrics_text = fs::read_to_string(&metrics).unwrap();
        assert!(metrics_text.contains(&format!("\"mode\": \"{mode}\"")));

        let report = constellation()
            .args([
                "bench-report",
                "--assignments",
                path(&assignments),
                "--truth",
                path(&tmp.path().join("sim_truth.tsv")),
                "--index",
                path(&index),
                "--metrics",
                path(&metrics),
            ])
            .output()
            .unwrap();
        assert_success(report);
    }

    let pulp_assignments = tmp.path().join("candidate-locus-pulp.assignments.tsv");
    let pulp_metrics = tmp.path().join("candidate-locus-pulp.metrics.json");
    assert_success(
        constellation()
            .args([
                "map",
                "--index",
                path(&index),
                "--r1",
                path(&tmp.path().join("sim_R1.fastq")),
                "--r2",
                path(&tmp.path().join("sim_R2.fastq")),
                "--mode",
                "candidate-locus",
                "--score-mode",
                "pulp",
                "--out",
                path(&pulp_assignments),
                "--emit-metrics",
                path(&pulp_metrics),
            ])
            .output()
            .unwrap(),
    );
    let pulp_metrics_text = fs::read_to_string(&pulp_metrics).unwrap();
    assert!(pulp_metrics_text.contains("\"score_mode\": \"pulp\""));

    let count_prefix = tmp.path().join("counts");
    assert_success(
        constellation()
            .args([
                "count",
                "--assignments",
                path(&pulp_assignments),
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
                path(&index),
                "--r1",
                path(&r1),
                "--r2",
                path(&r2),
                "--max-postings-per-seed",
                "256",
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
}

use std::{
    fs,
    path::{Path, PathBuf},
    process::{Command, Output},
};

fn fixture(name: &str) -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("../../examples/quickstart")
        .join(name)
}
fn p(path: &Path) -> &str {
    path.to_str().unwrap()
}
fn run(args: &[&str]) -> Output {
    Command::new(env!("CARGO_BIN_EXE_constellation"))
        .env("RAYON_NUM_THREADS", "1")
        .args(args)
        .output()
        .unwrap()
}
fn ok(args: &[&str]) {
    let out = run(args);
    assert!(
        out.status.success(),
        "{}",
        String::from_utf8_lossy(&out.stderr)
    );
}
fn json(path: &Path) -> serde_json::Value {
    serde_json::from_slice(&fs::read(path).unwrap()).unwrap()
}
fn index(path: &Path) {
    ok(&[
        "index-ec",
        "--threads",
        "2",
        "--transcripts",
        p(&fixture("reference.fa")),
        "--t2g-map",
        p(&fixture("t2g.tsv")),
        "--format",
        "ptrhash",
        "--syncmer-s",
        "20",
        "--out",
        p(path),
    ]);
}

#[test]
fn quant_matches_native_stages_and_demo_counts_for_both_policies() {
    let temp = tempfile::tempdir().unwrap();
    let idx = temp.path().join("reference index.ecidx");
    let r1 = fixture("R1.fastq");
    let r2 = fixture("R2.fastq");
    let barcodes = fixture("barcodes.txt");
    index(&idx);
    let mapped = temp.path().join("native.cstrad2");
    ok(&[
        "--threads",
        "2",
        "map",
        "--splicing-state",
        "--index",
        p(&idx),
        "--r1",
        p(&fixture("R1.fastq")),
        "--r2",
        p(&fixture("R2.fastq")),
        "--out",
        p(&mapped),
    ]);
    for policy in ["intersection", "usa"] {
        let output = temp.path().join(format!("sample {policy}"));
        let mut quant = vec![
            "quant",
            "--index",
            p(&idx),
            "--r1",
            p(&r1),
            "--r2",
            p(&r2),
            "--barcodes",
            p(&barcodes),
            "--out-dir",
            p(&output),
            "--threads",
            "2",
            "--resolution",
            policy,
        ];
        if policy == "usa" {
            quant.push("--keep-assignments");
        }
        ok(&quant);
        let prefix = temp.path().join(format!("native-{policy}"));
        let mut count = vec![
            "count-splicing",
            "--assignments",
            p(&mapped),
            "--barcode-whitelist",
            p(&barcodes),
            "--out-prefix",
            p(&prefix),
            "--resolution",
            policy,
        ];
        if policy == "usa" {
            count.push("--no-correct-umis");
        }
        ok(&count);
        for suffix in [
            "_matrix.mtx",
            "_spliced.mtx",
            "_unspliced.mtx",
            "_ambiguous.mtx",
            "_compatible_both.mtx",
            "_conflicting.mtx",
            "_unknown.mtx",
            "_features.tsv",
            "_barcodes.tsv",
        ] {
            assert_eq!(
                fs::read(output.join(format!("counts{suffix}"))).unwrap(),
                fs::read(format!("{}{suffix}", prefix.display())).unwrap(),
                "{policy}: {suffix}"
            );
        }
        let counts = json(&output.join("count.metrics.json"));
        assert_eq!(counts["rows"], 7);
        assert_eq!(counts["molecules_after_umi_correction"], 6);
        assert_eq!(counts["spliced"], 4);
        assert_eq!(counts["unspliced"], 1);
        assert_eq!(counts["compatible_both"], 1);
        assert_eq!(counts["correct_umis"], policy == "intersection");
        let manifest = json(&output.join("run.json"));
        assert_eq!(manifest["status"], "complete");
        assert_eq!(manifest["threads_requested"], 2);
        if policy == "usa" {
            assert_eq!(
                fs::read(output.join("assignments.cstrad2")).unwrap(),
                fs::read(&mapped).unwrap()
            );
        } else {
            assert!(!output.join("assignments.cstrad2").exists());
        }
        let before = fs::read(output.join("run.json")).unwrap();
        assert!(
            !run(&quant).status.success(),
            "existing output must be rejected"
        );
        assert_eq!(fs::read(output.join("run.json")).unwrap(), before);
    }
}

#[test]
fn quant_rejects_bad_inputs_and_preserves_failed_mapping() {
    let temp = tempfile::tempdir().unwrap();
    let idx = temp.path().join("ref.idx");
    index(&idx);
    let empty_barcodes = temp.path().join("empty.txt");
    fs::write(&empty_barcodes, "").unwrap();
    let output = temp.path().join("bad-whitelist");
    let result = run(&[
        "quant",
        "--index",
        p(&idx),
        "--r1",
        p(&fixture("R1.fastq")),
        "--r2",
        p(&fixture("R2.fastq")),
        "--barcodes",
        p(&empty_barcodes),
        "-o",
        p(&output),
    ]);
    assert!(!result.status.success());
    assert!(
        !output.exists(),
        "bad whitelist must fail before mapping/output creation"
    );
    let truncated = temp.path().join("truncated.fastq");
    fs::write(&truncated, "@broken\nACGT\n+\n").unwrap();
    let failed = temp.path().join("failed-map");
    assert!(!run(&[
        "quant",
        "--index",
        p(&idx),
        "--r1",
        p(&fixture("R1.fastq")),
        "--r2",
        p(&truncated),
        "--barcodes",
        p(&fixture("barcodes.txt")),
        "-o",
        p(&failed)
    ])
    .status
    .success());
    assert_eq!(json(&failed.join("run.json"))["status"], "failed");
    assert!(!failed.join("counts_matrix.mtx").exists());
    let zero_threads = run(&[
        "--threads",
        "0",
        "index-ec",
        "--transcripts",
        p(&fixture("reference.fa")),
        "--out",
        p(&temp.path().join("zero.idx")),
    ]);
    assert!(!zero_threads.status.success());
    assert!(String::from_utf8_lossy(&zero_threads.stderr).contains("invalid value"));
    let version = run(&["--version"]);
    assert!(version.status.success());
    assert!(String::from_utf8_lossy(&version.stdout).contains(env!("CARGO_PKG_VERSION")));
}

use std::{fs, process::Command};

fn run(args: &[&str]) {
    let result = Command::new(env!("CARGO_BIN_EXE_constellation"))
        .env("RAYON_NUM_THREADS", "2")
        .env("CONSTELLATION_EC_PTRHASH", "1")
        .args(args)
        .output()
        .unwrap();
    assert!(
        result.status.success(),
        "{}",
        String::from_utf8_lossy(&result.stderr)
    );
}

#[test]
fn consensus_preserves_ties_across_sampled_layouts_and_rejects_incompatible_options() {
    let tmp = tempfile::tempdir().unwrap();
    let fa = tmp.path().join("ref.fa");
    let v4 = tmp.path().join("v4.idx");
    let v9 = tmp.path().join("v9.idx");
    let r1 = tmp.path().join("r1.fq");
    let r2 = tmp.path().join("r2.fq");
    let seq = "ACGTGCAATGACTGGTACAGCTAGCTGACTGACATGCATCG";
    fs::write(
        &fa,
        format!(">a|gene=A exon\n{seq}\n>b|gene=B intron\n{seq}\n"),
    )
    .unwrap();
    fs::write(
        &r1,
        "@r\nACGTACGTACGTACGTACGTACGTACGT\n+\nIIIIIIIIIIIIIIIIIIIIIIIIIIII\n",
    )
    .unwrap();
    fs::write(&r2, format!("@r\n{seq}\n+\n{}\n", "I".repeat(seq.len()))).unwrap();
    run(&[
        "index-ec",
        "--transcripts",
        fa.to_str().unwrap(),
        "--k",
        "7",
        "--syncmer-s",
        "4",
        "--out",
        v4.to_str().unwrap(),
    ]);
    run(&[
        "convert-ec-index",
        "--index",
        v4.to_str().unwrap(),
        "--out",
        v9.to_str().unwrap(),
    ]);
    let mut outputs = Vec::new();
    for index in [&v4, &v9] {
        let out = tmp.path().join("out.tsv");
        let args = [
            "map",
            "--gene-ec-consensus",
            "--index",
            index.to_str().unwrap(),
            "--r1",
            r1.to_str().unwrap(),
            "--r2",
            r2.to_str().unwrap(),
            "--output-format",
            "tsv",
            "--out",
            out.to_str().unwrap(),
        ];
        run(&args);
        let text = fs::read_to_string(&out).unwrap();
        assert!(text.contains("ambiguous_gene\t.\t."), "{text}");
        outputs.push(text);
        let result = Command::new(env!("CARGO_BIN_EXE_constellation"))
            .args(args)
            .arg("--search-reverse-complement")
            .output()
            .unwrap();
        assert!(!result.status.success());
        assert!(String::from_utf8_lossy(&result.stderr).contains("forward-strand"));
    }
    assert_eq!(outputs[0], outputs[1]);
}

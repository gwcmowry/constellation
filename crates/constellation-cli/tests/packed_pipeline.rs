use std::fs;
use std::path::Path;
use std::process::Command;

fn run(args: &[&str]) {
    let output = Command::new(env!("CARGO_BIN_EXE_constellation"))
        .env("RAYON_NUM_THREADS", "2")
        .args(args)
        .output()
        .unwrap();
    assert!(
        output.status.success(),
        "{}",
        String::from_utf8_lossy(&output.stderr)
    );
}

fn path(p: &Path) -> &str {
    p.to_str().unwrap()
}

#[test]
fn packed_and_diagnostic_paths_have_identical_output() {
    let tmp = tempfile::tempdir().unwrap();
    let fasta = tmp.path().join("ref.fa");
    fs::write(&fasta, ">a|gene=G1 target:exon\nACGTGCAATGACTGGTACAGCTAGCTGACTGACATGCATCG\n>b|gene=G2 target:intron\nACGTGCAATGACTGGTACAGCTAGCTGACTGACATGCATCG\n>c|gene=G3\nTGCCAGGATCTGAACTCGATGACGTTAGCATGCAACCTGGA\n").unwrap();
    let r1 = tmp.path().join("r1.fq");
    let r2 = tmp.path().join("r2.fq");
    let barcodes = [
        "ACGTACGTACGTACGTTTTTCCCCAAAA",
        "",
        "acgt",
        "NNNNACGTACGTACGTCCCCAAAATTTT",
        "UUUUACGTACGTACGTCCCCAAAATTTT",
        "éACGTACGTACGTACGTCCCCAAAATTTT",
    ];
    let sequences = [
        "ACGTGCAATGACTGGTACAGCTAGCTGACTGACATGCATCG".to_owned(),
        "TGCCAGGATCTGAACTCGATGACGTTAGCATGCAACCTGGA".to_owned(),
        "ACGTNNNNGCAATGACTGGTACAGCTAGCTGACTGACATGCATCG".to_owned(),
        "ACGTGCAATGACTGGTACAGCTAGCTGACTGACATGCATCG".repeat(8),
        "acgtucaatgactggtacagctagctgactgacatgcatcg".to_owned(),
        "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA".to_owned(),
        "".to_owned(),
        "ACG".to_owned(),
        "AAGCAGTGGTATCAACGCAGAGTACATGGGACGTGCAATGACTGGTACAGCTAGCTGACTGACATGCATCG".to_owned(),
    ];
    let mut left = String::new();
    let mut right = String::new();
    for i in 0..300 {
        let cb = barcodes[i % barcodes.len()];
        let seq = &sequences[(i / barcodes.len()) % sequences.len()];
        left += &format!("@{i}\r\n{cb}\r\n+\r\n{}\r\n", "I".repeat(cb.len()));
        right += &format!(
            "@{i}\r\n{seq}\r\n+\r\n{}\r\n",
            if i % 7 == 0 { "!" } else { "I" }.repeat(seq.len())
        );
    }
    fs::write(&r1, left).unwrap();
    fs::write(&r2, right).unwrap();
    for sampled in [false, true] {
        let index = tmp.path().join(if sampled {
            "sampled.ecidx"
        } else {
            "dense.ecidx"
        });
        let mut args = vec![
            "index-ec",
            "--transcripts",
            path(&fasta),
            "--k",
            "7",
            "--out",
            path(&index),
        ];
        if sampled {
            args.extend(["--syncmer-s", "4"]);
        }
        run(&args);
        let mut indexes = vec![index.clone()];
        if sampled {
            let v9 = tmp.path().join("ptrhash.ecidx");
            let output = Command::new(env!("CARGO_BIN_EXE_constellation"))
                .env("RAYON_NUM_THREADS", "2")
                .env("CONSTELLATION_EC_PTRHASH", "1")
                .args([
                    "convert-ec-index",
                    "--index",
                    path(&index),
                    "--out",
                    path(&v9),
                ])
                .output()
                .unwrap();
            assert!(
                output.status.success(),
                "{}",
                String::from_utf8_lossy(&output.stderr)
            );
            indexes.push(v9);
        }
        for index in indexes {
            for format in ["tsv", "gene-ec-rad"] {
                for reverse_trim in [false, true] {
                    let fast = tmp.path().join("fast.out");
                    let slow = tmp.path().join("slow.out");
                    let diag = tmp.path().join("diagnostics.tsv");
                    let mut args = vec![
                        "map",
                        "--index",
                        path(&index),
                        "--r1",
                        path(&r1),
                        "--r2",
                        path(&r2),
                        "--output-format",
                        format,
                        "--out",
                        path(&fast),
                    ];
                    if reverse_trim {
                        args.extend(["--search-reverse-complement", "--trim-tso"]);
                    }
                    run(&args);
                    let out_pos = args.iter().position(|&s| s == "--out").unwrap() + 1;
                    args[out_pos] = path(&slow);
                    args.extend(["--emit-unmapped-diagnostics", path(&diag)]);
                    run(&args);
                    assert_eq!(
                        fs::read(&fast).unwrap(),
                        fs::read(&slow).unwrap(),
                        "{index:?}, {format}, reverse_trim={reverse_trim}"
                    );
                }
            }
        }
    }
}

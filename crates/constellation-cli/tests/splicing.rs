use std::{fs, io::Write, path::Path, process::Command};

fn run(args: &[&str]) -> std::process::Output {
    Command::new(env!("CARGO_BIN_EXE_constellation"))
        .env("RAYON_NUM_THREADS", "2")
        .env("CONSTELLATION_EC_PTRHASH", "1")
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
fn p(path: &Path) -> &str {
    path.to_str().unwrap()
}
fn dna(mut seed: u64) -> String {
    (0..96)
        .map(|_| {
            seed = seed.wrapping_mul(6364136223846793005).wrapping_add(1);
            b"ACGT"[(seed >> 62) as usize] as char
        })
        .collect()
}
fn umi(n: usize) -> String {
    let word: String = [(n >> 4) & 3, (n >> 2) & 3, n & 3]
        .into_iter()
        .map(|i| b"ACGT"[i] as char)
        .collect();
    word.repeat(4)
}

#[test]
fn state_output_and_molecule_counts_across_layouts_compression_and_corruption() {
    let tmp = tempfile::tempdir().unwrap();
    let fa = tmp.path().join("ref.fa");
    let t2g = tmp.path().join("t2g.tsv");
    let v4 = tmp.path().join("v4.idx");
    let v9 = tmp.path().join("v9.idx");
    let r1 = tmp.path().join("r1.fq");
    let r2 = tmp.path().join("r2.fq");
    let whitelist = tmp.path().join("whitelist.txt");
    let seqs: Vec<_> = (1..=6).map(dna).collect();
    let targets = [
        ("s", "G", "S", 0),
        ("u", "G", "U", 1),
        ("as", "G", "S", 2),
        ("au", "G", "U", 2),
        ("x", "G", "?", 3),
        ("gs", "G", "S", 4),
        ("hu", "H", "U", 4),
        ("h", "H", "U", 5),
    ];
    let mut fasta = String::new();
    let mut annotations = String::new();
    for (tx, gene, state, seq) in targets {
        fasta += &format!(">{tx}\n{}\n", seqs[seq]);
        annotations += &format!("{tx}\t{gene}\t{state}\n");
    }
    fs::write(&fa, fasta).unwrap();
    fs::write(&t2g, annotations).unwrap();
    let cb = "ACGTGTCAGTACGATC";
    fs::write(&whitelist, format!("{cb}\n")).unwrap();
    let cases = [
        (2, 0),
        (0, 0),
        (1, 1),
        (2, 2),
        (0, 3),
        (1, 3),
        (3, 4),
        (4, 5),
        (0, 5),
        (0, 6),
        (5, 6),
        (5, 7),
    ];
    let mut left = String::new();
    let mut right = String::new();
    for (i, (seq, id)) in cases.into_iter().enumerate() {
        let barcode = if i == 1 { "CCGTGTCAGTACGATC" } else { cb };
        left += &format!("@{i}\n{barcode}{}\n+\n{}\n", umi(id), "I".repeat(28));
        right += &format!("@{i}\n{}\n+\n{}\n", seqs[seq], "I".repeat(96));
    }
    let chimeric = format!("{}{}", &seqs[0][..45], &seqs[1][..45]);
    left += &format!("@12\n{cb}{}\n+\n{}\n", umi(8), "I".repeat(28));
    right += &format!("@12\n{chimeric}\n+\n{}\n", "I".repeat(90));
    fs::write(&r1, left).unwrap();
    fs::write(&r2, right).unwrap();
    ok(&[
        "index-ec",
        "--transcripts",
        p(&fa),
        "--t2g-map",
        p(&t2g),
        "--k",
        "15",
        "--syncmer-s",
        "9",
        "--out",
        p(&v4),
    ]);
    ok(&["convert-ec-index", "--index", p(&v4), "--out", p(&v9)]);
    let mut canonical = None;
    for index in [&v4, &v9] {
        let out = tmp.path().join("map.rad");
        ok(&[
            "map",
            "--splicing-state",
            "--index",
            p(index),
            "--r1",
            p(&r1),
            "--r2",
            p(&r2),
            "--out",
            p(&out),
        ]);
        let bytes = fs::read(&out).unwrap();
        assert_eq!(&bytes[..8], b"CSTRAD2\0");
        if let Some(ref expected) = canonical {
            assert_eq!(&bytes, expected);
        } else {
            canonical = Some(bytes.clone());
        }
        let text = tmp.path().join("map.tsv");
        ok(&[
            "map",
            "--splicing-state",
            "--index",
            p(index),
            "--r1",
            p(&r1),
            "--r2",
            p(&r2),
            "--out",
            p(&text),
            "--output-format",
            "tsv",
        ]);
        let tsv = fs::read_to_string(text).unwrap();
        let rows: Vec<_> = tsv.lines().collect();
        assert!(rows[0].ends_with("gene_states"));
        for (i, expected) in [
            "0:A", "0:S", "0:U", "0:A", "0:S", "0:U", "0:X", "0:S;1:U", "0:S", "0:S", "1:U", "1:U",
            "0:C",
        ]
        .into_iter()
        .enumerate()
        {
            assert_eq!(rows[i + 1].split('\t').next_back().unwrap(), expected);
        }
        let gz = tmp.path().join("map.rad.gz");
        let mut encoder = flate2::write::GzEncoder::new(Vec::new(), flate2::Compression::fast());
        encoder.write_all(&bytes).unwrap();
        let compressed = encoder.finish().unwrap();
        fs::write(&gz, &compressed).unwrap();
        let corrupt = tmp.path().join("corrupt.gz");
        let mut bad_gzip = compressed.clone();
        let crc = bad_gzip.len() - 8;
        bad_gzip[crc] ^= 1;
        fs::write(&corrupt, bad_gzip).unwrap();
        assert!(!run(&[
            "count-splicing",
            "--assignments",
            p(&corrupt),
            "--out-prefix",
            p(&tmp.path().join("corrupt"))
        ])
        .status
        .success());
        assert!(!tmp.path().join("corrupt_matrix.mtx").exists());
        let mut paths = vec![out.clone(), gz];
        if Command::new("zstd").arg("--version").output().is_ok() {
            let zst = tmp.path().join("map.rad.zst");
            ok(&[
                "map",
                "--splicing-state",
                "--index",
                p(index),
                "--r1",
                p(&r1),
                "--r2",
                p(&r2),
                "--out",
                p(&zst),
                "--output-compression",
                "zstd",
            ]);
            paths.push(zst);
        }
        for input in paths {
            let prefix = tmp.path().join("counts");
            let metrics = tmp.path().join("counts.json");
            ok(&[
                "count-splicing",
                "--assignments",
                p(&input),
                "--barcode-whitelist",
                p(&whitelist),
                "--emit-metrics",
                p(&metrics),
                "--out-prefix",
                p(&prefix),
            ]);
            let m: serde_json::Value = serde_json::from_slice(&fs::read(metrics).unwrap()).unwrap();
            for (key, n) in [
                ("rows", 13),
                ("spliced", 2),
                ("unspliced", 2),
                ("compatible_both", 1),
                ("conflicting_state", 2),
                ("unknown_state", 1),
                ("conflicting_gene_groups", 1),
                ("molecules_after_umi_correction", 8),
                ("barcode_corrected", 1),
            ] {
                assert_eq!(m[key], n, "{key}: {m}");
            }
            let sum = |name: &str| -> u32 {
                fs::read_to_string(tmp.path().join(format!("counts_{name}.mtx")))
                    .unwrap()
                    .lines()
                    .skip(2)
                    .map(|l| l.split_whitespace().last().unwrap().parse::<u32>().unwrap())
                    .sum()
            };
            assert_eq!(
                (
                    sum("spliced"),
                    sum("unspliced"),
                    sum("ambiguous"),
                    sum("matrix")
                ),
                (2, 2, 4, 8)
            );
            assert_eq!(
                (sum("compatible_both"), sum("conflicting"), sum("unknown")),
                (1, 2, 1)
            );
        }
        for bad in [
            bytes[..bytes.len() - 1].to_vec(),
            {
                let mut b = bytes.clone();
                *b.last_mut().unwrap() = 255;
                b
            },
            {
                let mut b = bytes.clone();
                b[6] = b'1';
                b
            },
        ] {
            let path = tmp.path().join("bad.rad");
            fs::write(&path, bad).unwrap();
            let result = run(&[
                "count-splicing",
                "--assignments",
                p(&path),
                "--out-prefix",
                p(&tmp.path().join("bad")),
            ]);
            assert!(!result.status.success());
            assert!(!tmp.path().join("bad_matrix.mtx").exists());
        }
        // A header-only stream has zero molecules and zero barcode rows.
        let mut header_end = 28;
        for _ in 0..2 {
            let len =
                u16::from_le_bytes(bytes[header_end..header_end + 2].try_into().unwrap()) as usize;
            header_end += 2 + len;
        }
        let empty = tmp.path().join("empty.rad");
        fs::write(&empty, &bytes[..header_end]).unwrap();
        ok(&[
            "count-splicing",
            "--assignments",
            p(&empty),
            "--out-prefix",
            p(&tmp.path().join("empty")),
        ]);
        assert!(fs::read(tmp.path().join("empty_barcodes.tsv"))
            .unwrap()
            .is_empty());
        assert!(fs::read_to_string(tmp.path().join("empty_matrix.mtx"))
            .unwrap()
            .contains("2 0 0"));
    }
}

use criterion::{Criterion, criterion_group, criterion_main};
use std::hint::black_box;
use std::path::PathBuf;
use std::process::Command;

fn bench_bam_fingerprint(c: &mut Criterion) {
    let bin = env!("CARGO_BIN_EXE_rsomics-bam-fingerprint");
    let manifest = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
    let bam = manifest.join("tests/golden/small.bam");
    let out = tempfile::NamedTempFile::new().unwrap();

    c.bench_function("rsomics-bam-fingerprint golden", |b| {
        b.iter(|| {
            let status = Command::new(black_box(bin))
                .args([
                    "--bam",
                    bam.to_str().unwrap(),
                    "--out-fingerprint",
                    out.path().to_str().unwrap(),
                ])
                .status()
                .unwrap();
            assert!(status.success());
        });
    });
}

criterion_group!(benches, bench_bam_fingerprint);
criterion_main!(benches);

//! ChIP-enrichment fingerprint, matching deeptools `plotFingerprint` defaults.
//!
//! Sampling, chunkSize derivation, per-bin coverage arithmetic, and output
//! format all match deeptools at one BAM / one process — see `sampling.rs`
//! and `coverage.rs` for the upstream-cited formulas. Output is the data tables
//! (raw counts, fingerprint curve, quality metrics), not the PNG.

#![allow(clippy::cast_precision_loss)]

use std::num::NonZero;
use std::path::Path;

use rsomics_common::{Result, RsomicsError};

mod coverage;
mod metrics;
mod output;
mod sampling;

pub use metrics::QualityMetrics;
pub use output::{write_fingerprint, write_quality_metrics, write_raw_counts};

#[derive(Debug, Clone)]
pub struct FingerprintOpts {
    /// Window size in bp (deeptools default: 500).
    pub bin_size: u32,
    /// Target number of sampled bins across the genome (deeptools default: 500000).
    pub number_of_samples: u64,
    /// Minimum mapping quality (deeptools default: 0 = no filter).
    pub min_mapq: u8,
    /// Skip reads whose FLAG has any of these bits set (deeptools `samFlagExclude`,
    /// default None → 0). Unmapped (0x4) reads are always skipped regardless.
    pub sam_flag_exclude: u16,
    /// Keep only reads whose FLAG has all of these bits set (deeptools
    /// `samFlagInclude`, default None → 0 = no requirement).
    pub sam_flag_include: u16,
}

impl Default for FingerprintOpts {
    fn default() -> Self {
        Self {
            bin_size: 500,
            number_of_samples: 500_000,
            min_mapq: 0,
            sam_flag_exclude: 0,
            sam_flag_include: 0,
        }
    }
}

pub struct Fingerprint {
    pub counts: Vec<u64>,
}

pub fn fingerprint(
    input: &Path,
    opts: &FingerprintOpts,
    workers: NonZero<usize>,
) -> Result<Fingerprint> {
    let chroms = sampling::read_chrom_geom(input, workers)?;
    let genome_size: u64 = chroms.iter().map(|c| c.length).sum();
    if genome_size == 0 {
        return Err(RsomicsError::InvalidInput(
            "BAM header declares zero genome length".into(),
        ));
    }

    let mapped = sampling::count_mapped(input, workers)?;
    if mapped == 0 {
        return Err(RsomicsError::InvalidInput(
            "no mapped reads found — check MAPPING quality / flag filters".into(),
        ));
    }

    let step_size = (genome_size / opts.number_of_samples).max(1);

    let mean_chrom_len = genome_size as f64 / chroms.len() as f64;
    if mean_chrom_len < step_size as f64 {
        let min_samples = (genome_size as f64 / mean_chrom_len) as u64;
        return Err(RsomicsError::InvalidInput(format!(
            "--number-of-samples has to be bigger than {min_samples}"
        )));
    }

    let chunk_size = sampling::chunk_length(step_size, opts.bin_size, mapped, genome_size, 1);
    let bin_size = u64::from(opts.bin_size);

    let (regions, n_bins) = sampling::layout_regions(&chroms, step_size, bin_size, chunk_size);
    if n_bins == 0 {
        return Err(RsomicsError::InvalidInput(
            "no bins were sampled — decrease --bin-size or --number-of-samples".into(),
        ));
    }

    let counts = coverage::accumulate_coverage(
        input, opts, workers, &chroms, &regions, n_bins, step_size, bin_size,
    )?;
    Ok(Fingerprint { counts })
}

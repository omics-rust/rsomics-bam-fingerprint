//! ChIP-enrichment fingerprint, matching deeptools `plotFingerprint` default semantics.
//!
//! The genome is sampled at evenly spaced windows (bins). Each bin's value is
//! the per-base read coverage summed over the bin. Sorting those per-bin values
//! and plotting their normalised cumulative sum against the normalised bin rank
//! gives a Lorenz curve: a diagonal means coverage is uniform (no enrichment),
//! a sharp elbow near the right means a few bins hold most of the signal (strong
//! ChIP enrichment).
//!
//! ## Deterministic sampling (source parity with deeptools)
//!
//! deeptools does **not** pick bins at random. With the whole genome it sets
//! `stepSize = max(genomeSize / numberOfSamples, 1)` and walks each chromosome
//! at `stepSize` intervals, taking one bin `[i, i + binSize)` per step. The walk
//! happens inside genome "chunks" (`mapReduce` partitions) of length `chunkSize`,
//! and a bin that would cross a chunk boundary (`i + binSize > chunkEnd`) is
//! dropped — so chunk geometry, not just stepSize, decides which bins exist.
//! `chunkSize` itself is derived from the BAM's mapped-read count:
//!
//! ```text
//! reads_per_bp = max_mapped / genomeSize
//! chunkSize    = floor(stepSize * 1000 / (reads_per_bp * nBams))
//! chunkSize    = max(chunkSize, stepSize, binSize)
//! ```
//!
//! Reproducing all of this makes our sampled-bin set identical to deeptools'
//! when run with one BAM at `--numberOfProcessors 1` (no task shuffle).
//!
//! ## Per-bin coverage rule (`SumCoveragePerBin`)
//!
//! For each bin `[reg0, reg1)` (a single-tile region, `tileSize = binSize`):
//! every read with `FLAG & 0x4 == 0` is decomposed into aligned blocks
//! (`get_blocks()`: reference-consuming M/=/X runs, split by D and N). Each block
//! contributes `min(blockEnd, reg1) - max(blockStart, reg0)` bases (clamped to
//! `[0, binSize]`) to the bin. Default filters: skip unmapped only — no MAPQ
//! floor, no FLAG include/exclude, duplicates kept, reads not extended
//! (`extendReads=False` → `defaultFragmentLength == "read length"`).
//!
//! ## Output (`--out-raw-counts`)
//!
//! Mirrors `plotFingerprint --outRawCounts`: a `#plotFingerprint --outRawCounts`
//! header, a quoted label line, then one integer per sampled bin (the per-base
//! coverage, formatted `%d` exactly as deeptools casts the float column). Bins
//! appear in genome order (chromosome header order, ascending position).
//!
//! ## Fingerprint + summary (`--out-fingerprint`, `--out-quality-metrics`)
//!
//! The fingerprint is the sorted, cumulative, max-normalised coverage vs the
//! normalised rank — the data behind deeptools' PNG. The quality metrics
//! (AUC, X-intercept, elbow) use deeptools' exact formulas over the sorted
//! cumulative curve.
//!
//! ## Scope
//!
//! We emit the data tables (raw counts, fingerprint curve, summary metrics), not
//! the PNG plot — the same split as `rsomics-bam-signal`, which emits bedGraph
//! rather than bigWig. Synthetic/JSD/CHANCE columns (only produced with
//! `--JSDsample` against a reference BAM) are out of scope for this
//! single-input crate.

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

/// The fingerprint result: per-bin per-base coverage in genome order.
pub struct Fingerprint {
    pub counts: Vec<u64>,
}

/// Compute the per-bin per-base coverage fingerprint for `input`.
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

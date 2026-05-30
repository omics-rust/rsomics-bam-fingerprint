// Genome sampling: chromosome geometry, chunkSize derivation (deeptools
// `get_chunk_length`), and region layout (`mapReduce` + `count_reads_in_region`).

use std::num::NonZero;
use std::path::Path;

use rsomics_bamio::raw::{self, RawRecord};
use rsomics_common::{Result, RsomicsError};

pub(crate) struct ChromGeom {
    pub(crate) tid: usize,
    pub(crate) length: u64,
}

/// A coverage region as deeptools constructs it in `count_reads_in_region`.
///
/// When `stepSize != binSize` each region is one bin (`n_bins == 1`); when
/// `stepSize == binSize` a whole genome chunk becomes one tiled region of
/// `n_bins` contiguous `bin_size` tiles (the partial tail bin is dropped, as
/// deeptools floor-divides). `out_lo` is this region's offset into the global
/// genome-ordered counts vector.
pub(crate) struct Region {
    pub(crate) tid: usize,
    pub(crate) start: u64,
    /// Exclusive end of the region (`reg[1]` in deeptools): the bin end for a
    /// single-bin region, the chunk end for a tiled region (may exceed
    /// `start + n_bins * bin_size` when the partial tail bin was floor-dropped).
    pub(crate) end: u64,
    pub(crate) n_bins: u64,
    pub(crate) out_lo: usize,
}

pub(crate) fn read_chrom_geom(input: &Path, workers: NonZero<usize>) -> Result<Vec<ChromGeom>> {
    let mut reader = rsomics_bamio::open_with_workers(input, workers)?;
    let header = reader.read_header().map_err(RsomicsError::Io)?;
    Ok(header
        .reference_sequences()
        .iter()
        .enumerate()
        .map(|(tid, (_, seq))| ChromGeom {
            tid,
            length: usize::from(seq.length()) as u64,
        })
        .collect())
}

/// Counts reads with `FLAG & 0x4 == 0`, matching pysam `.mapped` — the
/// idxstats total deeptools uses to derive `chunkSize`. No MAPQ or
/// include/exclude filtering enters this count.
pub(crate) fn count_mapped(input: &Path, workers: NonZero<usize>) -> Result<u64> {
    let mut reader = rsomics_bamio::open_with_workers(input, workers)?;
    reader.read_header().map_err(RsomicsError::Io)?;
    let mut record = RawRecord::default();
    let mut mapped: u64 = 0;
    while raw::read_record(reader.get_mut(), &mut record)? != 0 {
        if record.flags() & 0x4 == 0 && record.reference_sequence_id() >= 0 {
            mapped += 1;
        }
    }
    Ok(mapped)
}

/// deeptools `get_chunk_length` with `len(bamFilesHandles) == n_bams`.
pub(crate) fn chunk_length(
    step_size: u64,
    bin_size: u32,
    max_mapped: u64,
    genome_size: u64,
    n_bams: u64,
) -> u64 {
    let reads_per_bp = max_mapped as f64 / genome_size as f64;
    let mut chunk = (step_size as f64 * 1e3 / (reads_per_bp * n_bams as f64)) as u64;
    if chunk < step_size {
        chunk = step_size;
    }
    if u64::from(bin_size) > 0 && chunk < u64::from(bin_size) {
        chunk = u64::from(bin_size);
    }
    chunk
}

/// Builds deeptools' coverage regions in genome order. Mirrors the `mapReduce`
/// chunk walk (`range(0, size, chunkSize)`) crossed with `count_reads_in_region`:
///
/// - `stepSize == binSize` → each chunk is one tiled region spanning
///   `nBins = (chunkEnd - chunkStart) // binSize` contiguous tiles (partial tail
///   dropped, deeptools floor-divides).
/// - `stepSize != binSize` → one single-bin region per `range(chunkStart,
///   chunkEnd, stepSize)` step, skipping a bin that crosses the chunk end.
///
/// Returns the regions plus the total bin count (the length of the output).
pub(crate) fn layout_regions(
    chroms: &[ChromGeom],
    step_size: u64,
    bin_size: u64,
    chunk_size: u64,
) -> (Vec<Region>, usize) {
    let tiled = step_size == bin_size;
    let mut regions = Vec::new();
    let mut out_lo = 0usize;
    for chrom in chroms {
        let mut chunk_start = 0u64;
        while chunk_start < chrom.length {
            let chunk_end = (chunk_start + chunk_size).min(chrom.length);
            if tiled {
                let n_bins = (chunk_end - chunk_start) / bin_size;
                if n_bins > 0 {
                    regions.push(Region {
                        tid: chrom.tid,
                        start: chunk_start,
                        end: chunk_end,
                        n_bins,
                        out_lo,
                    });
                    out_lo += n_bins as usize;
                }
            } else {
                let mut i = chunk_start;
                while i < chunk_end {
                    if i + bin_size > chunk_end {
                        break;
                    }
                    regions.push(Region {
                        tid: chrom.tid,
                        start: i,
                        end: i + bin_size,
                        n_bins: 1,
                        out_lo,
                    });
                    out_lo += 1;
                    i += step_size;
                }
            }
            chunk_start += chunk_size;
        }
    }
    (regions, out_lo)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn chunk_length_floors_to_step_and_bin() {
        // reads_per_bp = 1.0; chunk = step*1000 / (1*1) = 500000, dominates.
        assert_eq!(chunk_length(500, 500, 1_000_000, 1_000_000, 1), 500_000);
        // dense: reads_per_bp large → chunk floored to step.
        let c = chunk_length(50, 50, 1_000_000, 1_000, 1);
        assert!(c >= 50);
    }
}

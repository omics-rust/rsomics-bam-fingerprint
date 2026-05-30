// Per-read coverage accumulation: CIGAR → aligned blocks → bin counts.
// Mirrors deeptools `SumCoveragePerBin.get_coverage_of_region`.

use std::num::NonZero;
use std::path::Path;

use rsomics_bamio::raw::{self, RawRecord};
use rsomics_common::{Result, RsomicsError};

use crate::FingerprintOpts;
use crate::sampling::Region;

const CIGAR_MATCH: u8 = 0;
const CIGAR_INSERTION: u8 = 1;
const CIGAR_DELETION: u8 = 2;
const CIGAR_SKIP: u8 = 3;
const CIGAR_SOFT_CLIP: u8 = 4;
const CIGAR_SEQ_MATCH: u8 = 7;
const CIGAR_SEQ_MISMATCH: u8 = 8;

/// pysam `get_blocks()`: aligned blocks of reference-consuming CIGAR. M/=/X
/// extend the current block; D, N **and I** all break it (pysam emits a fresh
/// block at every insertion even though I consumes no reference, so the two
/// blocks abut). Soft/hard clips and padding are ignored. Fills `out` (cleared
/// first) so the caller can reuse one allocation across records.
pub(crate) fn aligned_blocks(start0: u64, record: &RawRecord, out: &mut Vec<(u64, u64)>) {
    blocks_from_cigar(start0, record.cigar_ops(), out);
}

pub(crate) fn blocks_from_cigar(
    start0: u64,
    cigar: impl Iterator<Item = (u8, u32)>,
    out: &mut Vec<(u64, u64)>,
) {
    out.clear();
    let mut pos = start0;
    let mut block_start = start0;
    let mut in_block = false;
    for (kind, len) in cigar {
        let len = u64::from(len);
        match kind {
            CIGAR_MATCH | CIGAR_SEQ_MATCH | CIGAR_SEQ_MISMATCH => {
                if !in_block {
                    block_start = pos;
                    in_block = true;
                }
                pos += len;
            }
            // I, D and N all break the current block (pysam emits a fresh block
            // at every insertion too, even though I consumes no reference).
            CIGAR_INSERTION | CIGAR_DELETION | CIGAR_SKIP => {
                if in_block {
                    out.push((block_start, pos));
                    in_block = false;
                }
                if kind != CIGAR_INSERTION {
                    pos += len;
                }
            }
            CIGAR_SOFT_CLIP => {}
            _ => {}
        }
    }
    if in_block && pos > block_start {
        out.push((block_start, pos));
    }
}

/// Routes one read's aligned blocks into its chromosome's regions, applying
/// deeptools' `SumCoveragePerBin.get_coverage_of_region` arithmetic per region.
/// `chrom_regions` is ascending by start; the read can touch the (single) tiled
/// region of its chunk, or — in the sparse single-bin layout — at most a handful
/// of consecutive single-bin regions.
pub(crate) fn add_read(
    counts: &mut [i64],
    chrom_regions: &[Region],
    blocks: &[(u64, u64)],
    step_size: u64,
    bin_size: u64,
) {
    let Some(&(read_start, _)) = blocks.first() else {
        return;
    };
    let read_end = blocks.last().unwrap().1;

    let _ = step_size;
    // First region whose end is past read_start.
    let mut r = chrom_regions.partition_point(|reg| reg.end <= read_start);
    while r < chrom_regions.len() {
        let reg = &chrom_regions[r];
        if reg.start >= read_end {
            break;
        }
        cover_region(
            &mut counts[reg.out_lo..reg.out_lo + reg.n_bins as usize],
            reg,
            blocks,
            bin_size,
        );
        r += 1;
    }
}

/// deeptools `SumCoveragePerBin.get_coverage_of_region` for one region. `cov` is
/// the bin slice for this region (`tileSize == bin_size`, `nRegBins == n_bins`).
/// Per-base coverage is spread across tiles with deeptools' exact integer
/// arithmetic — including the `ceil`/`eIdx` clamp and `last_eIdx` block guard
/// that the upstream tiled path uses (and that make wide reads over-count a tile
/// rather than clipping to true per-base overlap; matched for byte parity).
pub(crate) fn cover_region(cov: &mut [i64], reg: &Region, blocks: &[(u64, u64)], bin_size: u64) {
    let reg0 = reg.start as i64;
    let reg1 = reg.end as i64;
    let tile = bin_size as i64;
    let n_reg_bins = reg.n_bins as i64;
    let len_cov = cov.len() as i64;
    let cov_end = reg0 + len_cov * tile;

    let mut last_eidx: Option<i64> = None;
    for &(bs, be) in blocks {
        let mut frag_start = bs as i64;
        let mut frag_end = be as i64;
        if frag_end - frag_start == 0 {
            continue;
        }
        if frag_end <= reg0 || frag_start >= reg1 {
            continue;
        }
        if frag_start < reg0 {
            frag_start = reg0;
        }
        if frag_end > cov_end {
            frag_end = cov_end;
        }

        let mut s_idx = ((frag_start - reg0) / tile).max(0);
        let mut e_idx = (div_ceil_i64(frag_end - reg0, tile)).min(n_reg_bins);
        if e_idx >= len_cov {
            e_idx = len_cov - 1;
        }
        if let Some(last) = last_eidx {
            s_idx = last.max(s_idx);
            if s_idx >= e_idx {
                continue;
            }
        }

        // First bin: partial overlap from frag_start to the bin's end.
        let first = if frag_end < reg0 + (s_idx + 1) * tile {
            frag_end - frag_start
        } else {
            reg0 + (s_idx + 1) * tile - frag_start
        };
        cov[s_idx as usize] += first.min(tile);

        let mut k = s_idx + 1;
        while k < e_idx {
            cov[k as usize] += tile;
            k += 1;
        }
        while e_idx - s_idx >= n_reg_bins {
            e_idx -= 1;
        }
        if e_idx > s_idx {
            let mut last = frag_end - (reg0 + e_idx * tile);
            if last > tile {
                last = tile;
            } else if last < 0 {
                last = 0;
            }
            cov[e_idx as usize] += last;
        }
        last_eidx = Some(e_idx);
    }
}

fn div_ceil_i64(a: i64, b: i64) -> i64 {
    (a + b - 1) / b
}

/// Single streaming pass: every kept read is routed to its chromosome's regions
/// and its aligned blocks are spread across the bins via deeptools'
/// `SumCoveragePerBin.get_coverage_of_region` arithmetic.
#[allow(clippy::too_many_arguments)]
pub(crate) fn accumulate_coverage(
    input: &Path,
    opts: &FingerprintOpts,
    workers: NonZero<usize>,
    chroms: &[crate::sampling::ChromGeom],
    regions: &[Region],
    n_bins: usize,
    step_size: u64,
    bin_size: u64,
) -> Result<Vec<u64>> {
    // Per-chromosome slice [lo, hi) into `regions`, indexed by tid. Regions are
    // grouped by chromosome in header order, ascending start within each.
    let mut chrom_span: Vec<(usize, usize)> = vec![(0, 0); chroms.len()];
    {
        let mut idx = 0usize;
        while idx < regions.len() {
            let tid = regions[idx].tid;
            let lo = idx;
            while idx < regions.len() && regions[idx].tid == tid {
                idx += 1;
            }
            chrom_span[tid] = (lo, idx);
        }
    }

    let mut counts = vec![0i64; n_bins];

    let mut reader = rsomics_bamio::open_with_workers(input, workers)?;
    reader.read_header().map_err(RsomicsError::Io)?;
    let mut record = RawRecord::default();
    let mut blocks: Vec<(u64, u64)> = Vec::new();

    while raw::read_record(reader.get_mut(), &mut record)? != 0 {
        let flags = record.flags();
        if flags & 0x4 != 0 {
            continue;
        }
        let tid = record.reference_sequence_id();
        if tid < 0 {
            continue;
        }
        if opts.sam_flag_include != 0 && (flags & opts.sam_flag_include) != opts.sam_flag_include {
            continue;
        }
        if opts.sam_flag_exclude != 0 && (flags & opts.sam_flag_exclude) != 0 {
            continue;
        }
        if opts.min_mapq > 0 && record.mapping_quality() < opts.min_mapq {
            continue;
        }

        let tid = tid as usize;
        let (lo, hi) = match chrom_span.get(tid) {
            Some(&span) if span.0 < span.1 => span,
            _ => continue,
        };
        let chrom_regions = &regions[lo..hi];

        let start0 = record.alignment_start() as u64;
        aligned_blocks(start0, &record, &mut blocks);
        add_read(&mut counts, chrom_regions, &blocks, step_size, bin_size);
    }

    Ok(counts.into_iter().map(|c| c.max(0) as u64).collect())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::sampling::Region;

    fn blocks(record_blocks: &[(u64, u64)]) -> Vec<(u64, u64)> {
        record_blocks.to_vec()
    }

    #[test]
    fn div_ceil_matches_python_ceil() {
        assert_eq!(div_ceil_i64(5, 200), 1);
        assert_eq!(div_ceil_i64(200, 200), 1);
        assert_eq!(div_ceil_i64(201, 200), 2);
        assert_eq!(div_ceil_i64(1046, 200), 6);
    }

    #[test]
    fn single_bin_region_is_clamped_per_base_overlap() {
        // One bin [1000,1200), one read block [1096,1146): contributes 50 bp.
        let reg = Region {
            tid: 0,
            start: 1000,
            end: 1200,
            n_bins: 1,
            out_lo: 0,
        };
        let mut cov = vec![0i64];
        cover_region(&mut cov, &reg, &blocks(&[(1096, 1146)]), 200);
        assert_eq!(cov[0], 50);
    }

    #[test]
    fn tiled_region_spreads_full_middle_tile() {
        // deeptools tiled quirk: a read [959,1009) with sIdx=4 eIdx=6 dumps a
        // full 200 into the middle tile (bin 5) even though it only reaches 1009.
        let reg = Region {
            tid: 0,
            start: 0,
            end: 10000,
            n_bins: 50,
            out_lo: 0,
        };
        let mut cov = vec![0i64; 50];
        cover_region(&mut cov, &reg, &blocks(&[(959, 1009)]), 200);
        assert_eq!(cov[5], 200);
        assert_eq!(cov[4], 41); // 1000 - 959
    }

    #[test]
    fn blocks_split_on_insertion_deletion_skip() {
        let mut out = Vec::new();
        // 30M5I25M from 1999 → pysam [(1999,2029),(2029,2054)] (split on I).
        blocks_from_cigar(
            1999,
            [(CIGAR_MATCH, 30), (CIGAR_INSERTION, 5), (CIGAR_MATCH, 25)].into_iter(),
            &mut out,
        );
        assert_eq!(out, vec![(1999, 2029), (2029, 2054)]);

        // 40M5D20M from 1199 → [(1199,1239),(1244,1264)] (split on D, gap 5).
        blocks_from_cigar(
            1199,
            [(CIGAR_MATCH, 40), (CIGAR_DELETION, 5), (CIGAR_MATCH, 20)].into_iter(),
            &mut out,
        );
        assert_eq!(out, vec![(1199, 1239), (1244, 1264)]);

        // 30M100N30M from 499 → [(499,529),(629,659)] (split on N, gap 100).
        blocks_from_cigar(
            499,
            [(CIGAR_MATCH, 30), (CIGAR_SKIP, 100), (CIGAR_MATCH, 30)].into_iter(),
            &mut out,
        );
        assert_eq!(out, vec![(499, 529), (629, 659)]);

        // 5S55M from 1599 → soft-clip ignored, [(1599,1654)].
        blocks_from_cigar(
            1599,
            [(CIGAR_SOFT_CLIP, 5), (CIGAR_MATCH, 55)].into_iter(),
            &mut out,
        );
        assert_eq!(out, vec![(1599, 1654)]);
    }
}

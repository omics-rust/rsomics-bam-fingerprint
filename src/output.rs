// Table writers for --out-raw-counts, --out-fingerprint, --out-quality-metrics.
// Format mirrors deeptools `plotFingerprint --outRawCounts`.

use std::io::{BufWriter, Write};

use rsomics_common::{Result, RsomicsError};

use crate::Fingerprint;

/// Write the `--outRawCounts` table: deeptools header, quoted label, one `%d`
/// per sampled bin in genome order.
pub fn write_raw_counts(out: &mut dyn Write, label: &str, fp: &Fingerprint) -> Result<()> {
    let mut w = BufWriter::new(out);
    writeln!(w, "#plotFingerprint --outRawCounts").map_err(RsomicsError::Io)?;
    writeln!(w, "'{label}'").map_err(RsomicsError::Io)?;
    for &c in &fp.counts {
        writeln!(w, "{c}").map_err(RsomicsError::Io)?;
    }
    w.flush().map_err(RsomicsError::Io)
}

/// Write the fingerprint curve: `rank<TAB>fraction`, one row per sampled bin.
pub fn write_fingerprint(out: &mut dyn Write, fp: &Fingerprint) -> Result<()> {
    let mut w = BufWriter::new(out);
    writeln!(w, "#rank\tfraction").map_err(RsomicsError::Io)?;
    for (x, y) in fp.cumulative_curve() {
        writeln!(w, "{x}\t{y}").map_err(RsomicsError::Io)?;
    }
    w.flush().map_err(RsomicsError::Io)
}

/// Write the quality-metrics summary: deeptools' AUC / X-intercept / elbow.
pub fn write_quality_metrics(out: &mut dyn Write, label: &str, fp: &Fingerprint) -> Result<()> {
    let m = fp.quality_metrics();
    let mut w = BufWriter::new(out);
    writeln!(w, "Sample\tAUC\tX-intercept\tElbow Point").map_err(RsomicsError::Io)?;
    writeln!(w, "{}\t{}\t{}\t{}", label, m.auc, m.x_int, m.elbow).map_err(RsomicsError::Io)?;
    w.flush().map_err(RsomicsError::Io)
}

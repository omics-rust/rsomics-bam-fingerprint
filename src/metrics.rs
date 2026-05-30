// Quality metrics over the sorted cumulative coverage curve.
// Formulas mirror deeptools `plotFingerprint.main`.

use crate::Fingerprint;

pub struct QualityMetrics {
    pub auc: f64,
    pub x_int: f64,
    pub elbow: f64,
}

impl Fingerprint {
    /// Sorted, cumulative, max-normalised curve (the y-axis of the plot).
    /// Returned paired with the normalised rank (x-axis), both length N.
    pub fn cumulative_curve(&self) -> Vec<(f64, f64)> {
        let n = self.counts.len();
        let mut sorted = self.counts.clone();
        sorted.sort_unstable();
        let total: u64 = sorted.iter().sum();
        let mut acc: u64 = 0;
        sorted
            .iter()
            .enumerate()
            .map(|(i, &v)| {
                acc += v;
                let x = i as f64 / n as f64;
                let y = if total == 0 {
                    0.0
                } else {
                    acc as f64 / total as f64
                };
                (x, y)
            })
            .collect()
    }

    pub fn quality_metrics(&self) -> QualityMetrics {
        let n = self.counts.len();
        let mut sorted = self.counts.clone();
        sorted.sort_unstable();
        let total: u64 = sorted.iter().sum();

        let mut counts = Vec::with_capacity(n);
        let mut acc: u64 = 0;
        for &v in &sorted {
            acc += v;
            counts.push(if total == 0 {
                0.0
            } else {
                acc as f64 / total as f64
            });
        }

        let auc = counts.iter().sum::<f64>() / n as f64;

        let x_int = counts
            .iter()
            .position(|&c| c > 0.0)
            .map_or(0.0, |k| (k + 1) as f64 / n as f64);

        // Maximum vertical distance from diagonal — the Lorenz-curve elbow.
        let elbow = if n <= 1 {
            0.0
        } else {
            let mut best_k = 0usize;
            let mut best = f64::NEG_INFINITY;
            for (i, &c) in counts.iter().enumerate() {
                let line = i as f64 / (n as f64 - 1.0);
                let diff = line - c;
                if diff > best {
                    best = diff;
                    best_k = i;
                }
            }
            (best_k + 1) as f64 / n as f64
        };

        QualityMetrics { auc, x_int, elbow }
    }
}

#[cfg(test)]
mod tests {
    use crate::Fingerprint;

    #[test]
    fn quality_metrics_uniform_curve() {
        // Uniform counts → AUC ≈ 0.5+, elbow near 1, x-int = 1/n.
        let fp = Fingerprint {
            counts: vec![10; 100],
        };
        let m = fp.quality_metrics();
        assert!((m.x_int - 0.01).abs() < 1e-12);
        assert!(m.auc > 0.0 && m.auc < 1.0);
    }

    #[test]
    fn quality_metrics_all_zero_is_safe() {
        let fp = Fingerprint {
            counts: vec![0; 10],
        };
        let m = fp.quality_metrics();
        assert_eq!(m.auc, 0.0);
        assert_eq!(m.x_int, 0.0);
    }
}

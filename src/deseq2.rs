use crate::{LocalFit, LocalRegressionConfig, LocfitError, PredictionMethod};

/// DESeq2-oriented local dispersion trend.
///
/// Inputs and outputs are normal-scale means and dispersions. Internally, the
/// local fit is performed on `ln(mean)` and `ln(dispersion)`, with normal-scale
/// means used as prior weights. Points with `dispersion < min_disp * 10` are
/// excluded from the local regression.
#[derive(Clone, Debug)]
pub struct Deseq2LocalDispersionTrend {
    kind: TrendKind,
}

#[derive(Clone, Debug)]
enum TrendKind {
    Constant { value: f64 },
    Local { fit: LocalFit },
}

impl Deseq2LocalDispersionTrend {
    /// Fit a DESeq2-style local dispersion trend.
    ///
    /// If no dispersions pass the `min_disp * 10` filter, this crate returns a
    /// constant trend at `min_disp`. This is an ergonomic Rust behavior for the
    /// initial implementation rather than a claim of exact DESeq2 behavior.
    pub fn fit(means: &[f64], disps: &[f64], min_disp: f64) -> Result<Self, LocfitError> {
        if means.len() != disps.len() {
            return Err(LocfitError::LengthMismatch {
                x: means.len(),
                y: disps.len(),
                weights: None,
            });
        }
        if means.is_empty() {
            return Err(LocfitError::EmptyInput);
        }
        if !min_disp.is_finite() || min_disp <= 0.0 {
            return Err(LocfitError::InvalidInput(
                "min_disp must be finite and greater than zero".to_string(),
            ));
        }

        let mut x = Vec::new();
        let mut y = Vec::new();
        let mut weights = Vec::new();
        let threshold = min_disp * 10.0;

        for (index, (&mean, &disp)) in means.iter().zip(disps).enumerate() {
            if !mean.is_finite() || mean <= 0.0 {
                return Err(LocfitError::InvalidInput(format!(
                    "mean at index {index} must be finite and greater than zero"
                )));
            }
            if !disp.is_finite() || disp <= 0.0 {
                return Err(LocfitError::InvalidInput(format!(
                    "dispersion at index {index} must be finite and greater than zero"
                )));
            }
            if disp >= threshold {
                x.push(mean.ln());
                y.push(disp.ln());
                weights.push(mean);
            }
        }

        if x.is_empty() {
            return Ok(Self {
                kind: TrendKind::Constant { value: min_disp },
            });
        }

        let fit = LocalFit::fit(
            &x,
            &y,
            Some(&weights),
            LocalRegressionConfig {
                prediction_method: PredictionMethod::LocfitHermiteApprox,
                ..LocalRegressionConfig::default()
            },
        )?;
        Ok(Self {
            kind: TrendKind::Local { fit },
        })
    }

    /// Predict one normal-scale dispersion for a normal-scale mean.
    pub fn predict_one(&self, mean: f64) -> Result<f64, LocfitError> {
        validate_prediction_mean(mean)?;
        match &self.kind {
            TrendKind::Constant { value } => Ok(*value),
            TrendKind::Local { fit } => fit.predict_one(mean.ln()).map(f64::exp),
        }
    }

    /// Predict normal-scale dispersions for normal-scale means.
    pub fn predict(&self, means: &[f64]) -> Result<Vec<f64>, LocfitError> {
        means.iter().map(|&mean| self.predict_one(mean)).collect()
    }
}

/// Convenience function for [`Deseq2LocalDispersionTrend::fit`].
pub fn fit_deseq2_local_dispersion_trend(
    means: &[f64],
    disps: &[f64],
    min_disp: f64,
) -> Result<Deseq2LocalDispersionTrend, LocfitError> {
    Deseq2LocalDispersionTrend::fit(means, disps, min_disp)
}

fn validate_prediction_mean(mean: f64) -> Result<(), LocfitError> {
    if !mean.is_finite() || mean <= 0.0 {
        return Err(LocfitError::InvalidInput(
            "prediction mean must be finite and greater than zero".to_string(),
        ));
    }
    Ok(())
}

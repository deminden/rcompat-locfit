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
    /// constant trend at `min_disp`. This is a pragmatic Rust fallback for
    /// all-filtered inputs rather than a claim of exact DESeq2 behavior.
    ///
    /// ```
    /// use rcompat_locfit::Deseq2LocalDispersionTrend;
    ///
    /// let means = [1.0, 2.0, 5.0, 10.0, 100.0, 1000.0];
    /// let disps = [0.5, 0.3, 0.2, 0.12, 0.06, 0.03];
    ///
    /// let trend = Deseq2LocalDispersionTrend::fit(&means, &disps, 1e-8)?;
    /// assert!(trend.predict_one(30.0)? > 0.0);
    /// # Ok::<(), rcompat_locfit::LocfitError>(())
    /// ```
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

        Self::fit_validated_logs(&x, &y, &weights, min_disp)
    }

    /// Fit a DESeq2-style local dispersion trend from precomputed log columns.
    ///
    /// `log_means` and `log_disps` are natural-log transformed means and
    /// dispersions. `means` are the corresponding normal-scale mean values,
    /// matching DESeq2's local trend fit. This constructor is useful when input
    /// data already carries R-computed log columns and callers want to avoid
    /// reintroducing normal-scale formatting or parsing noise.
    ///
    /// ```
    /// use rcompat_locfit::Deseq2LocalDispersionTrend;
    ///
    /// let means = [1.0, 2.0, 5.0, 10.0, 100.0, 1000.0];
    /// let disps = [0.5, 0.3, 0.2, 0.12, 0.06, 0.03];
    /// let log_means: Vec<_> = means.iter().map(|mean| f64::ln(*mean)).collect();
    /// let log_disps: Vec<_> = disps.iter().map(|disp| f64::ln(*disp)).collect();
    ///
    /// let trend =
    ///     Deseq2LocalDispersionTrend::fit_from_logs(&log_means, &log_disps, &means, 1e-8)?;
    /// assert!(trend.predict_log_dispersion_one(30.0_f64.ln())?.is_finite());
    /// # Ok::<(), rcompat_locfit::LocfitError>(())
    /// ```
    pub fn fit_from_logs(
        log_means: &[f64],
        log_disps: &[f64],
        means: &[f64],
        min_disp: f64,
    ) -> Result<Self, LocfitError> {
        if log_means.len() != log_disps.len() || means.len() != log_means.len() {
            return Err(LocfitError::LengthMismatch {
                x: log_means.len(),
                y: log_disps.len(),
                weights: Some(means.len()),
            });
        }
        if log_means.is_empty() {
            return Err(LocfitError::EmptyInput);
        }
        if !min_disp.is_finite() || min_disp <= 0.0 {
            return Err(LocfitError::InvalidInput(
                "min_disp must be finite and greater than zero".to_string(),
            ));
        }

        let threshold = (min_disp * 10.0).ln();
        let mut x = Vec::new();
        let mut y = Vec::new();
        let mut weights = Vec::new();

        for (index, ((&log_mean, &log_disp), &mean)) in
            log_means.iter().zip(log_disps).zip(means).enumerate()
        {
            if !log_mean.is_finite() {
                return Err(LocfitError::InvalidInput(format!(
                    "log mean at index {index} must be finite"
                )));
            }
            if !log_disp.is_finite() {
                return Err(LocfitError::InvalidInput(format!(
                    "log dispersion at index {index} must be finite"
                )));
            }
            if !mean.is_finite() || mean <= 0.0 {
                return Err(LocfitError::InvalidInput(format!(
                    "mean at index {index} must be finite and greater than zero"
                )));
            }
            if log_disp >= threshold {
                x.push(log_mean);
                y.push(log_disp);
                weights.push(mean);
            }
        }

        Self::fit_validated_logs(&x, &y, &weights, min_disp)
    }

    /// Predict one normal-scale dispersion for a normal-scale mean.
    ///
    /// ```
    /// use rcompat_locfit::fit_deseq2_local_dispersion_trend;
    ///
    /// let means = [1.0, 2.0, 5.0, 10.0, 100.0, 1000.0];
    /// let disps = [0.5, 0.3, 0.2, 0.12, 0.06, 0.03];
    /// let trend = fit_deseq2_local_dispersion_trend(&means, &disps, 1e-8)?;
    ///
    /// assert!(trend.predict_one(30.0)? > 0.0);
    /// # Ok::<(), rcompat_locfit::LocfitError>(())
    /// ```
    pub fn predict_one(&self, mean: f64) -> Result<f64, LocfitError> {
        validate_prediction_mean(mean)?;
        match &self.kind {
            TrendKind::Constant { value } => Ok(*value),
            TrendKind::Local { .. } => self.predict_log_dispersion_one(mean.ln()).map(f64::exp),
        }
    }

    /// Predict one log-scale dispersion for a log-scale mean.
    ///
    /// ```
    /// use rcompat_locfit::fit_deseq2_local_dispersion_trend_from_logs;
    ///
    /// let means = [1.0, 2.0, 5.0, 10.0, 100.0, 1000.0];
    /// let disps = [0.5, 0.3, 0.2, 0.12, 0.06, 0.03];
    /// let log_means: Vec<_> = means.iter().map(|mean| f64::ln(*mean)).collect();
    /// let log_disps: Vec<_> = disps.iter().map(|disp| f64::ln(*disp)).collect();
    /// let trend =
    ///     fit_deseq2_local_dispersion_trend_from_logs(&log_means, &log_disps, &means, 1e-8)?;
    ///
    /// assert!(trend.predict_log_dispersion_one(30.0_f64.ln())?.is_finite());
    /// # Ok::<(), rcompat_locfit::LocfitError>(())
    /// ```
    pub fn predict_log_dispersion_one(&self, log_mean: f64) -> Result<f64, LocfitError> {
        validate_prediction_log_mean(log_mean)?;
        match &self.kind {
            TrendKind::Constant { value } => Ok(value.ln()),
            TrendKind::Local { fit } => fit.predict_one(log_mean),
        }
    }

    /// Predict normal-scale dispersions for normal-scale means.
    ///
    /// ```
    /// use rcompat_locfit::fit_deseq2_local_dispersion_trend;
    ///
    /// let means = [1.0, 2.0, 5.0, 10.0, 100.0, 1000.0];
    /// let disps = [0.5, 0.3, 0.2, 0.12, 0.06, 0.03];
    /// let trend = fit_deseq2_local_dispersion_trend(&means, &disps, 1e-8)?;
    ///
    /// let predictions = trend.predict(&[3.0, 30.0, 300.0])?;
    /// assert_eq!(predictions.len(), 3);
    /// # Ok::<(), rcompat_locfit::LocfitError>(())
    /// ```
    pub fn predict(&self, means: &[f64]) -> Result<Vec<f64>, LocfitError> {
        means.iter().map(|&mean| self.predict_one(mean)).collect()
    }

    /// Predict log-scale dispersions for log-scale means.
    ///
    /// ```
    /// use rcompat_locfit::fit_deseq2_local_dispersion_trend_from_logs;
    ///
    /// let means = [1.0, 2.0, 5.0, 10.0, 100.0, 1000.0];
    /// let disps = [0.5, 0.3, 0.2, 0.12, 0.06, 0.03];
    /// let log_means: Vec<_> = means.iter().map(|mean| f64::ln(*mean)).collect();
    /// let log_disps: Vec<_> = disps.iter().map(|disp| f64::ln(*disp)).collect();
    /// let trend =
    ///     fit_deseq2_local_dispersion_trend_from_logs(&log_means, &log_disps, &means, 1e-8)?;
    ///
    /// let predictions = trend.predict_log_dispersion(&[3.0_f64.ln(), 30.0_f64.ln()])?;
    /// assert_eq!(predictions.len(), 2);
    /// # Ok::<(), rcompat_locfit::LocfitError>(())
    /// ```
    pub fn predict_log_dispersion(&self, log_means: &[f64]) -> Result<Vec<f64>, LocfitError> {
        log_means
            .iter()
            .map(|&log_mean| self.predict_log_dispersion_one(log_mean))
            .collect()
    }

    fn fit_validated_logs(
        log_means: &[f64],
        log_disps: &[f64],
        weights: &[f64],
        min_disp: f64,
    ) -> Result<Self, LocfitError> {
        if log_means.is_empty() {
            return Ok(Self {
                kind: TrendKind::Constant { value: min_disp },
            });
        }

        let fit = LocalFit::fit(
            log_means,
            log_disps,
            Some(weights),
            LocalRegressionConfig {
                prediction_method: PredictionMethod::LocfitHermiteApprox,
                ..LocalRegressionConfig::default()
            },
        )?;
        Ok(Self {
            kind: TrendKind::Local { fit },
        })
    }
}

/// Convenience function for [`Deseq2LocalDispersionTrend::fit`].
///
/// ```
/// use rcompat_locfit::fit_deseq2_local_dispersion_trend;
///
/// let means = [1.0, 2.0, 5.0, 10.0, 100.0, 1000.0];
/// let disps = [0.5, 0.3, 0.2, 0.12, 0.06, 0.03];
///
/// let trend = fit_deseq2_local_dispersion_trend(&means, &disps, 1e-8)?;
/// assert!(trend.predict_one(30.0)? > 0.0);
/// # Ok::<(), rcompat_locfit::LocfitError>(())
/// ```
pub fn fit_deseq2_local_dispersion_trend(
    means: &[f64],
    disps: &[f64],
    min_disp: f64,
) -> Result<Deseq2LocalDispersionTrend, LocfitError> {
    Deseq2LocalDispersionTrend::fit(means, disps, min_disp)
}

/// Convenience function for [`Deseq2LocalDispersionTrend::fit_from_logs`].
///
/// ```
/// use rcompat_locfit::fit_deseq2_local_dispersion_trend_from_logs;
///
/// let means = [1.0, 2.0, 5.0, 10.0, 100.0, 1000.0];
/// let disps = [0.5, 0.3, 0.2, 0.12, 0.06, 0.03];
/// let log_means: Vec<_> = means.iter().map(|mean| f64::ln(*mean)).collect();
/// let log_disps: Vec<_> = disps.iter().map(|disp| f64::ln(*disp)).collect();
///
/// let trend = fit_deseq2_local_dispersion_trend_from_logs(
///     &log_means,
///     &log_disps,
///     &means,
///     1e-8,
/// )?;
/// assert!(trend.predict_log_dispersion_one(30.0_f64.ln())?.is_finite());
/// # Ok::<(), rcompat_locfit::LocfitError>(())
/// ```
pub fn fit_deseq2_local_dispersion_trend_from_logs(
    log_means: &[f64],
    log_disps: &[f64],
    means: &[f64],
    min_disp: f64,
) -> Result<Deseq2LocalDispersionTrend, LocfitError> {
    Deseq2LocalDispersionTrend::fit_from_logs(log_means, log_disps, means, min_disp)
}

fn validate_prediction_mean(mean: f64) -> Result<(), LocfitError> {
    if !mean.is_finite() || mean <= 0.0 {
        return Err(LocfitError::InvalidInput(
            "prediction mean must be finite and greater than zero".to_string(),
        ));
    }
    Ok(())
}

fn validate_prediction_log_mean(log_mean: f64) -> Result<(), LocfitError> {
    if !log_mean.is_finite() {
        return Err(LocfitError::InvalidInput(
            "prediction log mean must be finite".to_string(),
        ));
    }
    Ok(())
}

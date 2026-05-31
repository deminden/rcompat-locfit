use crate::error::LocfitError;

/// Kernel used to convert scaled distance into local regression weight.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Kernel {
    /// Tricube kernel: `(1 - |u|^3)^3` for `|u| < 1`, otherwise zero.
    Tricube,
}

/// Prediction strategy for a fitted local regression.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum PredictionMethod {
    /// Fit the local polynomial directly at each requested prediction point.
    Direct,
    /// Approximate R `locfit`'s two-step prediction behavior for larger
    /// one-dimensional fits: compute direct fits on a small evaluation grid,
    /// then use cubic Hermite interpolation between those grid points and
    /// weighted-quadratic extrapolation outside the grid.
    LocfitHermiteApprox,
}

/// Configuration for one-dimensional local polynomial regression.
///
/// `alpha` is the nearest-neighbor fraction used to choose the adaptive
/// bandwidth. Prior weights supplied to [`crate::LocalFit`] are multiplied by
/// kernel weights at prediction time.
#[derive(Clone, Debug)]
pub struct LocalRegressionConfig {
    /// Nearest-neighbor fraction. The current implementation supports
    /// `0 < alpha <= 1`.
    pub alpha: f64,
    /// Local polynomial degree. Supported values are 0, 1, and 2.
    pub degree: usize,
    /// Distance kernel.
    pub kernel: Kernel,
    /// Minimum number of finite, positively weighted points required at fit
    /// time.
    pub min_points: usize,
    /// If true, singular local fits are retried with lower polynomial degree.
    pub allow_degree_downgrade: bool,
    /// Prediction strategy used by [`crate::LocalFit::predict_one`] and
    /// [`crate::LocalFit::predict`].
    pub prediction_method: PredictionMethod,
}

impl Default for LocalRegressionConfig {
    fn default() -> Self {
        Self {
            alpha: 0.7,
            degree: 2,
            kernel: Kernel::Tricube,
            min_points: 3,
            allow_degree_downgrade: true,
            prediction_method: PredictionMethod::Direct,
        }
    }
}

impl LocalRegressionConfig {
    pub(crate) fn validate(&self) -> Result<(), LocfitError> {
        if !self.alpha.is_finite() || self.alpha <= 0.0 {
            return Err(LocfitError::InvalidConfig(
                "alpha must be finite and greater than zero".to_string(),
            ));
        }
        if self.alpha > 1.0 {
            return Err(LocfitError::InvalidConfig(
                "alpha values greater than 1 are not supported yet".to_string(),
            ));
        }
        if self.degree > 2 {
            return Err(LocfitError::InvalidConfig(
                "degree must be 0, 1, or 2".to_string(),
            ));
        }
        if self.min_points == 0 {
            return Err(LocfitError::InvalidConfig(
                "min_points must be greater than zero".to_string(),
            ));
        }
        Ok(())
    }
}

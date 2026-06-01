//! Clean-room Rust implementation of selected R `locfit`-compatible local
//! regression behavior.
//!
//! This crate currently focuses on the one-dimensional local polynomial path
//! needed for DESeq2-style local dispersion trend fitting. It does not bind to
//! R, does not call R at runtime, and does not contain R `locfit` source code.
//! The DESeq2-oriented path is fixture-validated against black-box R `locfit`;
//! broader exact R `locfit` parity is still in progress.
//!
//! # Example
//!
//! ```
//! use rcompat_locfit::fit_deseq2_local_dispersion_trend;
//!
//! let means = [1.0, 2.0, 5.0, 10.0, 100.0, 1000.0];
//! let disps = [0.5, 0.3, 0.2, 0.12, 0.06, 0.03];
//!
//! let trend = fit_deseq2_local_dispersion_trend(&means, &disps, 1e-8)?;
//! let prediction = trend.predict_one(30.0)?;
//! assert!(prediction.is_finite() && prediction > 0.0);
//! # Ok::<(), rcompat_locfit::LocfitError>(())
//! ```

pub mod config;
pub mod deseq2;
pub mod error;
pub mod kernel;
pub mod local_fit;
pub mod wls;

pub use config::{Kernel, LocalRegressionConfig, PredictionMethod};
pub use deseq2::{
    fit_deseq2_local_dispersion_trend, fit_deseq2_local_dispersion_trend_from_logs,
    Deseq2LocalDispersionTrend,
};
pub use error::LocfitError;
pub use local_fit::LocalFit;

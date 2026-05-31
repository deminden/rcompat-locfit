//! Clean-room Rust implementation of selected R `locfit`-compatible local
//! regression behavior.
//!
//! This crate currently focuses on the one-dimensional local polynomial path
//! needed for DESeq2-style local dispersion trend fitting. It does not bind to
//! R, does not call R at runtime, and does not contain R `locfit` source code.
//! Numeric parity with R `locfit` is fixture-driven and still in progress.

pub mod config;
pub mod deseq2;
pub mod error;
pub mod kernel;
pub mod local_fit;
pub mod wls;

pub use config::{Kernel, LocalRegressionConfig, PredictionMethod};
pub use deseq2::{fit_deseq2_local_dispersion_trend, Deseq2LocalDispersionTrend};
pub use error::LocfitError;
pub use local_fit::LocalFit;

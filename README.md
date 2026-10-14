# rcompat-locfit

`rcompat-locfit` is a pure Rust crate for selected R `locfit`-compatible
one-dimensional local regression. Its main target is the local dispersion
trend used by DESeq2-style workflows.

This is a clean-room implementation. It does not call R at runtime and contains
no code copied or translated from R `locfit`, DESeq2, or other copyleft projects.
R is used only as a black-box numeric oracle during development.

Version **0.1.5** improves evaluation-grid construction and quadratic fitting
without changing the public API. See [the changelog](CHANGELOG.md) for details.
Broader exact `locfit` parity remains work in progress, including ties,
singular local designs, and behavior across R builds and control settings.

## Scope

- Local polynomial regression in one predictor, with degrees 0, 1, and 2.
- Tricube kernels, nearest-neighbor bandwidths, and prior weights.
- Direct fitted values and slopes, plus an adaptive evaluation grid with cubic
  Hermite interpolation and weighted-quadratic boundary extrapolation.
- A DESeq2-oriented wrapper that fits `log(dispersion)` against `log(mean)`,
  uses normal-scale means as prior weights, and excludes dispersions below
  `min_disp * 10`.
- Normal-scale predictions and a log-column constructor with log-dispersion
  predictions for callers that already have precomputed natural logs.

This crate implements selected behavior rather than the full R `locfit` API.

## Current precision

Measured maximum errors for version **0.1.5** on macOS arm64:

| Validation set | Error measure | Maximum error |
| --- | --- | ---: |
| Synthetic `locfit` matrix | Relative dispersion error | `1.8970e-13` |
| Curated real DESeq2-derived subset | Relative dispersion error | `1.3319e-13` |
| Hard-real subset, stored-log API: 72 predictions | Absolute log-dispersion error | `1.9645e-11` |
| Hard-real subset, normal API with identical oracle inputs | Relative dispersion error | `1.3053e-11` |
| Actual DESeq2 local fit: 1,198 gene predictions | Absolute log-dispersion error | `1.6741e-11` |

The quadratic-fitting improvements achieved **39.04–39.78× lower log error**
on seven difficult predictions, **292.69×** and **221.85×** on two others,
and an exact oracle match on a tenth. Maximum absolute log error across these
ten predictions is **`1.1902e-13`**. The table reports accuracy across the full
validation sets.

The DESeq2 comparison uses R 4.6.1, DESeq2 1.52.0, and `locfit` 1.5-9.12,
running `DESeq(fitType="local")` on seeded simulated counts with 1,200 genes
and 12 samples. These measurements validate the local dispersion trend.
[The fixture guide](fixtures/README.md) documents oracle provenance and the
per-case gains.

## Usage

```rust
use rcompat_locfit::{
    fit_deseq2_local_dispersion_trend, fit_deseq2_local_dispersion_trend_from_logs,
};

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let means = vec![1.0, 2.0, 5.0, 10.0, 100.0, 1000.0];
    let disps = vec![0.5, 0.3, 0.2, 0.12, 0.06, 0.03];

    let trend = fit_deseq2_local_dispersion_trend(&means, &disps, 1e-8)?;
    let predicted = trend.predict(&[3.0, 30.0, 300.0])?;

    let log_means: Vec<_> = means.iter().map(|mean| f64::ln(*mean)).collect();
    let log_disps: Vec<_> = disps.iter().map(|disp| f64::ln(*disp)).collect();
    let log_trend =
        fit_deseq2_local_dispersion_trend_from_logs(&log_means, &log_disps, &means, 1e-8)?;
    let predicted_log_disp = log_trend.predict_log_dispersion(&[30.0_f64.ln()])?;

    println!("{predicted:?}");
    println!("{predicted_log_disp:?}");
    Ok(())
}
```

## Validation

Run the self-contained checks without R:

```bash
cargo test --locked
cargo test --locked --release
```

Tests cover synthetic edge cases, real-derived subsets, exact input transport,
precision regressions, and the DESeq2 local-fit comparison.

R is needed only to regenerate fixtures. See [the fixture guide](fixtures/README.md)
for commands, input formats, and dataset details.

## License

MIT. See [LICENSE](LICENSE).

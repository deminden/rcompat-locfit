# rcompat-locfit

**Version 0.2.0** — pure Rust one-dimensional local regression for the local dispersion trend used
by DESeq2. This clean-room implementation uses **R 4.6.2** as its numerical
reference. It has no R runtime dependency and contains no code copied or
translated from locfit or DESeq2.

Build with **Rust 1.85 or newer**, using the **Rust 2024 edition**.

## Behavior

- Local polynomial degrees 0, 1, and 2, tricube kernels, nearest-neighbor
  bandwidths, and prior weights.
- Adaptive evaluation grids with cubic Hermite interpolation and quadratic
  boundary extrapolation.
- DESeq2 fitting of `log(dispersion)` against `log(mean)`, weighted by the
  original means, with the `min_disp * 10` fitting filter.
- Normal-input and precomputed-log APIs, with normal and log predictions.

The crate implements selected locfit behavior; broader compatibility remains
work in progress.

## Current precision

Maximum absolute log-dispersion errors across all **23,480 usable gene
predictions**, measured through both input APIs against actual DESeq2 1.52.0
local fits under R 4.6.2 and locfit 1.5-9.12 on macOS arm64:

| Dataset | Gene predictions | Maximum absolute log error |
| --- | ---: | ---: |
| Real airway counts | 22,369 | `2.6473e-13` |
| Airway subset | 512 | `1.6701e-11` |
| Seeded simulated counts | 599 | `8.3134e-13` |

Each entry reports the larger maximum from the two input APIs. Self-contained
regression tests check **133 fixed queries**, including **96 hard predictions**,
through both APIs with dataset-specific absolute error limits.

All reference values use binary64 transport. These measurements concern the
local dispersion trend, not every output of the DESeq2 pipeline.

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

```bash
cargo test --locked
cargo test --locked --release
cargo run --release --example precision
```

Default tests are self-contained and require no R installation or external
private data. See [the fixture guide](fixtures/README.md) for reference
provenance, regeneration commands, and the fixed-query benchmark.

Tests separate numerical unit checks, public API behavior and validation, R
reference regressions, and documentation examples. CI runs debug and release
tests on Linux and macOS with both stable Rust and the minimum Rust 1.85.0.

## License

MIT. See [LICENSE](LICENSE).

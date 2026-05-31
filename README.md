# rcompat-locfit

`rcompat-locfit` is a fixture-driven Rust crate for selected R
`locfit`-compatible one-dimensional local regression behavior. The first target
is the local dispersion trend path used by DESeq2-style workflows.

This is a clean-room Rust implementation. It is not a binding to R `locfit`,
does not call R at runtime, and does not contain code copied or translated from
R `locfit`, DESeq2, or other copyleft projects. R `locfit` is used only as a
black-box oracle for generating numeric test fixtures.

Current status: validated for targeted DESeq2 local-dispersion compatibility
work. The DESeq2 wrapper matches black-box R `locfit` very closely on the real
DESeq2-derived fixtures used during development. Broader exact R `locfit`
parity is still work in progress, especially around adaptive evaluation grids,
ties, and singular local designs.

## Scope

- Pure Rust local polynomial regression in one predictor.
- Tricube kernel, adaptive nearest-neighbor bandwidth, and prior weights.
- Local fitted values and local slopes, so future R-style Hermite
  interpolation can be implemented on top of the direct fits.
- DESeq2-oriented wrapper:
  - normal-scale means and dispersions as input;
  - log(mean) and log(dispersion) fit internally;
  - original means as prior weights;
  - an approximate R `locfit`-style Hermite evaluation grid for larger fits;
  - floor-based nearest-neighbor bandwidth selection for R-style prediction;
  - R-style handling for tiny rank-deficient quadratic cells;
  - R-style weighted-quadratic extrapolation outside the evaluation range;
  - normal-scale dispersion predictions as output.

This crate is not a full port of R `locfit` and does not attempt to clone its
public API.

## Parity

The current implementation is validated against black-box R `locfit` fixtures.
Recent local measurements were run with R 4.6.0 and `locfit` 1.5-9.12:

| Fixture set | Max relative error |
| --- | ---: |
| Committed real DESeq2-derived subset | `1.10e-8` |
| Ignored full real DESeq2 hard rows | `3.19e-9` |
| Ignored full real DESeq2 all rows | `3.19e-9` |
| Synthetic R `locfit` matrix | `3.16e-8` |

The synthetic matrix intentionally includes edge cases that are outside the
main DESeq2-sized path. It remains useful as a compatibility target for future
changes even though the current committed cases are now close to the R oracle.

## Usage

```rust
use rcompat_locfit::fit_deseq2_local_dispersion_trend;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let means = vec![1.0, 2.0, 5.0, 10.0, 100.0, 1000.0];
    let disps = vec![0.5, 0.3, 0.2, 0.12, 0.06, 0.03];

    let trend = fit_deseq2_local_dispersion_trend(&means, &disps, 1e-8)?;
    let predicted = trend.predict(&[3.0, 30.0, 300.0])?;

    println!("{predicted:?}");
    Ok(())
}
```

## Fixture Generation

Fixtures are generated with R `locfit` as a black-box oracle. The script does
not require DESeq2.

```bash
Rscript fixtures/r/generate_locfit_fixtures.R
cargo test --test fixtures -- --ignored
```

The committed real DESeq2-derived subset can be regenerated from the ignored
local `/data/` debug tables with:

```bash
Rscript fixtures/r/generate_real_deseq_subset_fixture.R
```

The default `cargo test` run includes the committed synthetic fixture matrix
and the committed real DESeq2-derived subset parity check. It does not require
R. Full real-data diagnostics are ignored by default because they depend on
local `/data/` tables.

## License

MIT. See `LICENSE`.

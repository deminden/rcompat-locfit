# R 4.6.2 fixtures

All numerical oracles use R 4.6.2. R is a black-box reference; no locfit or
DESeq2 implementation source is copied or translated into this crate.

| File | Purpose |
| --- | --- |
| `deseq2_r462_hard_cases.csv` | Complete fitting inputs and difficult predictions from actual DESeq2 runs |
| `deseq2_r462_hard_queries.csv` | Fixed gene identities: 32 hard queries per dataset plus mean-count quantiles |
| `deseq2_r462_baseline.csv` | Frozen pre-refinement errors for both APIs and expected prediction bits |
| `locfit_r462_cases.csv` | Synthetic coverage of boundaries, ties, small fits, and filtering |

The two case fixtures share a schema with case, row kind, gene/query identity,
and little-endian IEEE754 binary64 hexadecimal values. The query manifest
records gene identities and selection categories; the baseline records expected
prediction bits and decimal absolute errors. Fitting inputs and expected
predictions pass from R to Rust without decimal rounding. Default tests require
neither R nor external data.

## DESeq2 datasets

The reference environment has R 4.6.2, DESeq2 1.52.0, locfit 1.5-9.12, and
airway 1.32.0. Runtime, package, platform, and query-manifest metadata are
recorded in the fixture.

- `airway_full`: 22,369 real genes with raw row sums at least 10; design
  `~ cell + dex`, eight samples, 8,942 active fitting observations.
- `airway_512`: 512 equally spaced raw-mean ranks from the same counts,
  with 213 active fitting observations.
- `simulated_600`: 600 genes and 12 samples, seed 462 and `betaSD=1`, yielding
  599 usable gene predictions and 579 active fitting observations.

Every observation entering the local trend is retained. Predictions include
low-dispersion genes excluded from fitting, and queries outside the fitting
range exercise extrapolation.
The query manifest preserves 32 largest observed errors per case, considering
both APIs, plus mean-count quantiles with duplicates removed. The 133 query
identities are not reselected when the solver changes.

## Regeneration

Run from the repository root with R 4.6.2 and the package versions listed above:

```bash
Rscript fixtures/r/generate_deseq2_r462_hard_fixtures.R
Rscript fixtures/r/generate_locfit_r462_fixtures.R
cargo run --release --example precision
```

Generators reject other R versions. Repeated native generation must reproduce
the same numeric bytes. Compiler/package builds can affect floating-point
results, so fixture changes must be reviewed independently of solver changes.

## Precision gate

The baseline was frozen before numerical refinement. The strict test requires
each of the 96 selected hard predictions to have at most one tenth of its
baseline error, through both the normal-input and stored-log APIs. Each
case's maximum across every query must also shrink tenfold. The original
expected prediction bits are checked against the baseline.

The same test also enforces absolute log-error limits for every committed
prediction through both APIs, including range queries:

| Dataset | Committed queries | Maximum allowed absolute log error |
| --- | ---: | ---: |
| `airway_512` | 44 | `1e-11` |
| `airway_full` | 44 | `3e-13` |
| `simulated_600` | 45 | `9e-13` |

These ceilings round the measured 0.2.0 committed-query maxima upward with a
small margin for platform rounding. They prevent regressions that could still
pass the historical tenfold gate. Both checks share one evaluation of each
dataset through each API.

`examples/precision.rs` reports errors for the committed fixture by default,
or for a fixture supplied as its only argument. The committed baseline is a
permanent regression reference. Range queries remain part of the overall
maximum even when their baseline error is zero.

## Test organization

- Unit tests live beside the numerical code.
- `tests/api.rs` groups public API checks into `local_regression`, `deseq2`, and
  `validation` modules under `tests/api/`, compiled as one integration target.
- `tests/fixtures.rs` checks fixture integrity, both hard-case precision gates,
  and synthetic shape parity; `tests/support/` shares fixture loading and
  evaluation with the reporting example.
- Rust documentation examples run as doctests.

Run `cargo test --locked --test api` for API checks or
`cargo test --locked --test fixtures` for numerical reference checks. The full
debug and release suites run in CI on Linux and macOS with stable Rust and
Rust 1.85.0; formatting, Clippy, documentation, and packaging checks run once
on Linux stable. Publishing depends on all four test jobs succeeding.

## Numerical implementation

The local fitter preserves observation order even when predictors repeat,
reuses the reciprocal bandwidth, and groups diagonal normalization factors.
The centered global fit enforces its zero first-moment invariant. These are
uniform numerical rules, with no sample-count or fixture-specific dispatch.

## Verified 0.2.0 results

The strict committed test passes for all 96 hard queries through both APIs;
minimum individual error reduction is 16.048×. The maximum over all 133 queries
also improves more than tenfold in every case.

A broader native-R verification evaluates every usable gene from the same
actual DESeq2 runs, including predictions outside the committed query list:

```bash
Rscript fixtures/r/generate_deseq2_r462_hard_fixtures.R target/all_gene_r462.csv --all-queries
cargo run --release --example precision -- target/all_gene_r462.csv
```

| Dataset | Predictions | Before max log error | After max log error | Reduction |
| --- | ---: | ---: | ---: | ---: |
| `airway_full` | 22,369 | `4.283629007062473e-12` | `2.647326802218686e-13` | 16.18× |
| `airway_512` | 512 | `1.112738101660682e-8` | `1.670130700404115e-11` | 666.26× |
| `simulated_600` | 599 | `2.386317810021410e-10` | `8.313350008393172e-13` | 287.05× |

The selected-query airway-subset maximum is smaller (`9.8703e-12`); the table
reports the broader all-gene maximum. The native runtime is R 4.6.2; the locfit
binary's original build metadata is recorded separately for reproducibility.

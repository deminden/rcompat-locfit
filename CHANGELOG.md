# Changelog

## 0.2.0 — 2026-10-29

- Use R 4.6.2 as the sole numerical oracle and remove obsolete reference
  tables, generators, hardcoded legacy values, and missing-data diagnostics.
- Adopt Rust 2024, requiring Rust 1.85 or newer.
- Preserve input-order reductions for repeated predictors, reuse reciprocal
  bandwidths, group local covariance normalization, and enforce centered
  global first moments. Remove sample-count-specific arithmetic dispatch and
  residual-based selection between competing solves.
- Improve all 96 fixed difficult DESeq2 predictions by at least 16.048× through
  both APIs. Maximum errors over all 23,480 usable gene predictions improve
  16.18× on the full airway dataset, 666.26× on its 512-gene subset, and 287.05×
  on the seed-462 simulation.
- Add self-contained binary64 fixtures, an immutable before-change baseline,
  strict precision regressions, and the `precision` reporting example.
- Protect current precision with dataset-specific absolute log-error ceilings
  alongside the frozen-baseline improvement checks, using the same evaluations.
- Group public API tests by regression behavior, DESeq2 behavior, and input
  validation in one integration target. Run debug and release tests on Linux
  and macOS with stable Rust and the minimum Rust 1.85.0 in CI.

The public Rust API is unchanged. R is required only to regenerate references.

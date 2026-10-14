# TODO

## Precision and portability

- Reduce the largest remaining frozen-fixture log error, currently
  `blood_full/range_03` at approximately `1.96e-11`, while preserving the
  achieved ten-case improvement and all existing parity bounds.
- Validate the quadratic solver and floating-point traversal across R builds
  and Rust platforms using identical binary inputs.
- Broaden actual DESeq2 local-fit validation beyond the seed-461 holdout, with
  additional sample sizes, count distributions, and real datasets.
- Rerun the five optional full-data diagnostics when their external tables are
  available; the 0.1.5 validation covers the self-contained fixtures.

## Broader locfit compatibility

- Validate adaptive `rbox()` subdivision and Hermite interpolation across more
  nearest-neighbor fractions, clustered predictors, and repeated values.
- Confirm rank-deficient quadratic cells, singular designs, and floor-based
  bandwidth selection across more `locfit` control settings.
- Confirm weighted-quadratic boundary extrapolation across more control settings.
- Expand generated and compact real-derived fixtures when new cases expose
  distinct behavior.
- Decide whether to expose only DESeq2 behavior or a broader locfit-compatible API.
- Add a cargo license audit before publishing.

# TODO

- Confirm R `locfit` parity over additional generated/randomized fixture matrices.
- Investigate exact R `locfit` evaluation grid / interpolation behavior.
- Refine the approximate one-dimensional cubic Hermite interpolation over fitted evaluation points.
- Refine the approximate default `rbox()` evaluation-point generation for the DESeq2 path.
- Confirm rank-deficient quadratic cell behavior across more `locfit` control settings.
- Confirm floor-based nearest-neighbor bandwidth behavior across more `locfit` control settings.
- Investigate whether the global weighted-quadratic boundary extrapolation generalizes across more `locfit` control settings.
- Investigate tie behavior and repeated x values.
- Investigate singular local design behavior.
- Decide whether to expose only DESeq2 behavior or a broader locfit-compatible API later.
- Add cargo license audit before publishing.

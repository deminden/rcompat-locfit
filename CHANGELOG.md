# Changelog

## 0.1.5 — 2026-10-14

Improves DESeq2 local dispersion trend precision while retaining the existing
public API and original parity bounds.

- Build the Hermite evaluation grid by adaptive subdivision using predictor
  geometry and nearest-neighbor bandwidths instead of sample-count templates.
- Accumulate one shared cubic covariance moment and select between equivalent
  quadratic solves using the smaller relative equation residual.
- Add exact binary-input transport for R oracles, with separate normal-input
  and stored-log expectations and a bit-identity regression check.
- Reduce log error on ten difficult predictions: seven improve by
  39.04–39.78×, two by 292.69× and 221.85×, and one matches the oracle exactly.
  Maximum absolute log error on these cases is `1.1902e-13`.
- Add a reproducible actual DESeq2 local-fit holdout on seeded simulated counts.
  Maximum log error over 1,198 usable genes falls from `2.39e-2` to `1.67e-11`.

Maximum absolute log error across all 72 hard-real predictions is `1.9645e-11`.
Measurements use macOS arm64; the fixture guide records oracle provenance and
per-case results.

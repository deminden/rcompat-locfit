# Fixtures

The committed fixtures use R `locfit` and DESeq2 only as black-box numeric
oracles. Default tests require neither R nor the external full-data tables.

| Fixture | Purpose |
| --- | --- |
| `locfit_deseq2_cases.csv` | Synthetic trends, boundaries, repeated means, tiny weighted fits, and filtering |
| `real_deseq2_subset_cases.csv` | Compact real-derived subset |
| `hard_real_2026_subset_cases.csv` | Three 192-point hard-real cases and 72 predictions; original frozen log oracle |
| `hard_real_2026_oracle_inputs.csv` | Exact binary inputs for all 648 input and prediction rows |
| `hard_real_2026_identical_input_oracles.csv` | Separate normal-input and stored-log outputs from identical doubles |
| `hard_real_2026_worst10_baseline.csv` | Original ten largest log errors at commit `e5b68e2`, used for individual tenfold targets |
| `deseq2_r461_local_trend.csv` | Actual DESeq2 local-fit holdout on seeded simulated counts |

Generate the synthetic matrix from the repository root:

```bash
Rscript fixtures/r/generate_locfit_fixtures.R
```

Run the committed fixture parity checks with:

```bash
cargo test --test fixtures
```

The default `cargo test` run already checks the committed real DESeq2-derived
subset, the committed 2026 hard-real subset, the synthetic matrix, the binary
input identity, the frozen ten-case targets, and the DESeq2 holdout. Full real-data diagnostics remain ignored by default because they use
local `/data/` tables. When those ignored local tables are available, the
ignored fixture tests run full real-data diagnostics, including the compressed
all-row table.

The newer ignored 2026 hard-real bundle can also be used for local diagnostics
when present at `data/locfit_hard_real_2026-06-01`. Those tests use the
public DESeq2 wrapper across all exported hard rows, use the bundle's stored
log columns through the public log-column API for the tight global-hardest
precision probe, and also run a broader prediction-grid guard against raw
locfit grid outputs.

## Real DESeq2-Derived Subset

`real_deseq2_subset_cases.csv` is a compact committed fixture generated from
the ignored local `/data/` debug tables:

- `data/local_dispersion_fit_points.tsv`
- `data/local_dispersion_ranked_hard_rows.tsv`
- `data/local_dispersion_all_rows.tsv.gz`

It intentionally does not commit the full data tables. Instead, it samples
real fit rows across the base-mean range and prediction rows from low-mean hard
cases plus range quantiles. Expected values are regenerated with black-box R
`locfit` on the committed subset so the fixture is self-contained and does not
pretend that a small subset can reproduce the original full-data DESeq2 fit.

Regenerate it from the repository root with:

```bash
Rscript fixtures/r/generate_real_deseq_subset_fixture.R
```

## 2026 Hard-Real Subset

`hard_real_2026_subset_cases.csv` is a compact committed fixture generated
from selected contrasts in the ignored 2026 hard-real bundle. It samples
real-shaped fit rows, anonymizes row labels, and regenerates expected values
with black-box R `locfit` on the committed subset. The default tests exercise
both the normal-scale public wrapper and the log-column API on this fixture.

Regenerate it from the repository root with:

```bash
Rscript fixtures/r/generate_hard_real_2026_subset_fixture.R
```

## Identical binary inputs

R and Rust can parse the same decimal text into different last bits. The
binary-input snapshot transports each double as 16 hexadecimal characters
representing its eight IEEE754 bytes in little-endian order. The R generator
decodes these bytes with `readBin` and verifies a `writeBin` roundtrip. A Rust
test checks the snapshot against the actual API inputs bit for bit.

Regenerate the transport and its separate oracle with an R environment
containing `locfit`:

```bash
cargo run --example export_hard_real_oracle_inputs -- fixtures/hard_real_2026_oracle_inputs.csv
Rscript fixtures/r/generate_identical_input_oracles.R fixtures/hard_real_2026_oracle_inputs.csv fixtures/hard_real_2026_identical_input_oracles.csv
```

The oracle contains `expected_normal`, `expected_normal_log`, and
`expected_stored_log`. The normal fit computes logs inside R from the decoded
normal inputs; the stored-log fit uses the decoded log columns directly. The
legacy normal-wrapper check still compares against its original fixture. A
separate identical-input normal check uses the new oracle with a relative-error
limit of `5e-11`.

The original log-oracle values and `hard_real_2026_worst10_baseline.csv` are
regression anchors. The ten-case test compares against those original values,
with each error limited to one tenth of its frozen baseline. Regenerating the
separate oracle does not change those anchors. R version and build differences
can still affect numerical results even with identical input bytes; snapshot
metadata records the oracle version, platform, and input digest.

## Measured per-case precision gains

The following log-API results compare 0.1.5 against 0.1.4 on the ten predictions
with the largest original error. Both versions are compared against the same
original log-oracle values. Gains are the original absolute error divided by
the current absolute error, rounded to two decimal places.

| Case / prediction | Current absolute log error | Error reduction |
| --- | ---: | ---: |
| `testis_full/range_01` | `6.7779e-14` | 292.69× |
| `testis_full/range_02` | `3.2252e-14` | 221.85× |
| `testis_full/range_11` | `0` | Exact oracle match |
| `testis_full/hard_00` | `1.1902e-13` | 39.04× |
| `testis_full/hard_01` | `1.1546e-13` | 39.17× |
| `testis_full/hard_02` | `1.1369e-13` | 39.61× |
| `testis_full/hard_03` | `1.0836e-13` | 39.16× |
| `testis_full/hard_04` | `1.0658e-13` | 39.78× |
| `testis_full/hard_05` | `1.0658e-13` | 39.58× |
| `testis_full/hard_06` | `1.0658e-13` | 39.50× |

Across the entire 72-prediction hard-real set, maximum absolute log error is
`1.9645e-11` at `blood_full/range_03`. The legacy public-wrapper check against
its original stored-log-derived oracle also passes, with maximum relative error
`1.6101e-8`; its inputs include formatting and log-transformation differences.
The separate identical-input normal comparison reaches `1.3053e-11`.

## DESeq2 R 4.6.1 holdout

`deseq2_r461_local_trend.csv` was generated with R 4.6.1, DESeq2 1.52.0, and
`locfit` 1.5-9.12. The generator calls `makeExampleDESeqDataSet` with 1,200 genes,
12 samples, `betaSD=1`, and seed 461, then runs `DESeq(fitType="local")`. The
1,198 usable gene rows contain binary encodings of means, gene-wise dispersion
estimates, precomputed logs, and expected log fitted dispersions.

```bash
Rscript fixtures/r/generate_deseq2_r461_fixture.R
cargo test --locked --test fixtures real_deseq2_r461_local_trend_parity
```

The test checks both the normal-input and stored-log APIs at an absolute
log-error limit of `2e-11`. It validates the local trend on simulated counts
through the actual DESeq2 package; it does not validate every DESeq2 output.

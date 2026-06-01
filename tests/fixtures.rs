use std::collections::{BTreeMap, HashSet};
use std::fs;
use std::path::Path;
use std::process::Command;

use rcompat_locfit::{
    fit_deseq2_local_dispersion_trend, fit_deseq2_local_dispersion_trend_from_logs,
};

const FIXTURE_PATH: &str = "fixtures/locfit_deseq2_cases.csv";
const REAL_DESEQ2_FIXTURE_PATH: &str = "fixtures/real_deseq2_subset_cases.csv";
const HARD_REAL_2026_COMMITTED_FIXTURE_PATH: &str = "fixtures/hard_real_2026_subset_cases.csv";
const REAL_DESEQ2_FULL_FIT_POINTS_PATH: &str = "data/local_dispersion_fit_points.tsv";
const REAL_DESEQ2_FULL_HARD_ROWS_PATH: &str = "data/local_dispersion_ranked_hard_rows.tsv";
const REAL_DESEQ2_FULL_ALL_ROWS_PATH: &str = "data/local_dispersion_all_rows.tsv.gz";
const HARD_REAL_2026_ROOT: &str = "data/locfit_hard_real_2026-06-01";
const HARD_REAL_2026_GLOBAL_HARDEST_PATH: &str =
    "data/locfit_hard_real_2026-06-01/global_hardest_2048.tsv";
const HARD_REAL_2026_GLOBAL_HARD_ROWS_PATH: &str =
    "data/locfit_hard_real_2026-06-01/global_hard_rows.tsv";
const MIN_DISP: f64 = 1e-8;
const DEFAULT_SYNTHETIC_REL_TOLERANCE: f64 = 5e-13;
const FEW_POINTS_REL_TOLERANCE: f64 = 2e-9;
const FILTERED_MIN_DISP_REL_TOLERANCE: f64 = 2e-9;
const REAL_DESEQ2_COMMITTED_REL_TOLERANCE: f64 = 2e-13;
const FULL_REAL_HARD_ROWS_REL_TOLERANCE: f64 = 2e-12;
const FULL_REAL_ALL_ROWS_REL_TOLERANCE: f64 = 2e-12;
const HARD_REAL_2026_COMMITTED_PUBLIC_REL_TOLERANCE: f64 = 5e-8;
const HARD_REAL_2026_COMMITTED_LOG_TOLERANCE: f64 = 5e-11;
const HARD_REAL_2026_PUBLIC_REL_TOLERANCE: f64 = 5e-8;
const HARD_REAL_2026_LOG_REL_TOLERANCE: f64 = 5e-10;
const HARD_REAL_2026_GRID_LOG_REL_TOLERANCE: f64 = 3e-2;

#[derive(Debug, Default)]
struct Case {
    means: Vec<f64>,
    disps: Vec<f64>,
    predictions: Vec<Prediction>,
}

#[derive(Debug)]
struct Prediction {
    mean: f64,
    expected: f64,
}

#[derive(Debug, Default)]
struct RealCase {
    means: Vec<f64>,
    disps: Vec<f64>,
    predictions: Vec<RealPrediction>,
}

#[derive(Debug)]
struct RealPrediction {
    gene: String,
    mean: f64,
    expected: f64,
    source_expected: Option<f64>,
}

#[derive(Debug, Default)]
struct HardRealCase {
    means: Vec<f64>,
    disps: Vec<f64>,
    log_means: Vec<f64>,
    log_disps: Vec<f64>,
    predictions: Vec<HardRealPrediction>,
}

#[derive(Debug)]
struct HardRealPrediction {
    label: String,
    mean: f64,
    expected: f64,
    log_mean: f64,
    expected_log: f64,
}

#[test]
fn generated_fixture_file_is_optional_for_default_tests() {
    if !Path::new(FIXTURE_PATH).exists() {
        eprintln!(
            "No fixture CSV found. Generate it with: Rscript fixtures/r/generate_locfit_fixtures.R"
        );
    }
}

#[test]
fn generated_fixture_matrix_covers_targeted_shapes_when_present() {
    if !Path::new(FIXTURE_PATH).exists() {
        return;
    }

    let cases = load_cases();
    for name in [
        "smooth_monotone",
        "wide_dynamic_range",
        "wide_extrapolation",
        "edge_predictions",
        "five_weighted_points",
        "seven_weighted_points",
        "repeated_means",
        "clustered_repeated_means",
        "few_points",
        "filtered_min_disp",
        "constant_min_disp",
    ] {
        assert!(cases.contains_key(name), "missing fixture case {name}");
    }
}

#[test]
fn constant_min_disp_behavior_is_strict() {
    let means = [1.0, 2.0, 5.0, 10.0, 50.0, 100.0];
    let disps = [1e-9, 2e-9, 3e-9, 4e-9, 5e-9, 8e-9];
    let trend = fit_deseq2_local_dispersion_trend(&means, &disps, MIN_DISP).unwrap();
    for mean in [1.0, 10.0, 100.0, 1000.0] {
        assert_eq!(trend.predict_one(mean).unwrap(), MIN_DISP);
    }
}

#[test]
fn real_deseq2_committed_subset_is_present_and_compact() {
    let cases = load_real_cases(REAL_DESEQ2_FIXTURE_PATH);
    let case = cases
        .get("real_deseq2_curated")
        .expect("missing committed real DESeq2 curated fixture case");

    assert_eq!(case.means.len(), 384);
    assert!(case.predictions.len() >= 30);
    assert!(case
        .means
        .iter()
        .all(|value| value.is_finite() && *value > 0.0));
    assert!(case
        .disps
        .iter()
        .all(|value| value.is_finite() && *value > 0.0));
    assert!(case
        .predictions
        .iter()
        .all(|prediction| prediction.expected.is_finite() && prediction.expected > 0.0));
    assert!(case.predictions.iter().any(|prediction| {
        prediction.source_expected.is_some()
            && relative_error(prediction.expected, prediction.source_expected.unwrap()) > 1e-2
    }));
}

#[test]
fn hard_real_2026_committed_subset_is_present_and_compact() {
    let cases = load_hard_real_2026_committed_cases();

    assert_eq!(cases.len(), 3);
    for (name, case) in cases {
        assert_eq!(case.means.len(), 192, "{name}");
        assert_eq!(case.log_means.len(), case.means.len(), "{name}");
        assert_eq!(case.log_disps.len(), case.disps.len(), "{name}");
        assert!(case.predictions.len() >= 20, "{name}");
        assert!(case
            .means
            .iter()
            .all(|value| value.is_finite() && *value > 0.0));
        assert!(case.log_means.iter().all(|value| value.is_finite()));
        assert!(case.log_disps.iter().all(|value| value.is_finite()));
        assert!(case
            .predictions
            .iter()
            .all(|prediction| prediction.expected_log.is_finite()));
    }
}

#[test]
fn r_locfit_fixture_parity() {
    let cases = load_cases();
    if cases.is_empty() {
        eprintln!(
            "No fixture cases found. Generate them with: Rscript fixtures/r/generate_locfit_fixtures.R"
        );
        return;
    }

    let mut failures = Vec::new();
    let mut max_rel = 0.0_f64;
    let mut max_label = String::new();
    for (name, case) in cases {
        let trend = fit_deseq2_local_dispersion_trend(&case.means, &case.disps, MIN_DISP).unwrap();
        for prediction in case.predictions {
            let actual = trend.predict_one(prediction.mean).unwrap();
            let rel = relative_error(actual, prediction.expected);
            if rel > max_rel {
                max_rel = rel;
                max_label = format!("{name}/mean={}", prediction.mean);
            }
            // Case-specific tolerances make remaining compatibility gaps
            // explicit while keeping the usual generated fixtures near-exact.
            let tolerance = synthetic_case_tolerance(&name);
            if rel > tolerance {
                failures.push(format!(
                    "{name}: mean={} actual={} expected={} rel={} tolerance={}",
                    prediction.mean, actual, prediction.expected, rel, tolerance
                ));
            }
        }
    }

    eprintln!("synthetic R locfit fixture max relative error: {max_rel} at {max_label}");
    assert!(
        failures.is_empty(),
        "fixture parity failures:\n{}",
        failures.join("\n")
    );
}

#[test]
fn real_deseq2_curated_subset_parity() {
    let cases = load_real_cases(REAL_DESEQ2_FIXTURE_PATH);
    let mut failures = Vec::new();
    let mut max_rel = 0.0_f64;

    for (name, case) in cases {
        let trend = fit_deseq2_local_dispersion_trend(&case.means, &case.disps, MIN_DISP).unwrap();
        for prediction in case.predictions {
            let actual = trend.predict_one(prediction.mean).unwrap();
            let rel = relative_error(actual, prediction.expected);
            max_rel = max_rel.max(rel);
            if rel > REAL_DESEQ2_COMMITTED_REL_TOLERANCE {
                failures.push(format!(
                    "{name}/{}: mean={} actual={} expected={} rel={}",
                    prediction.gene, prediction.mean, actual, prediction.expected, rel
                ));
            }
        }
    }

    eprintln!("real DESeq2 committed subset max relative error: {max_rel}");
    assert!(
        failures.is_empty(),
        "real DESeq2 curated fixture parity failures:\n{}",
        failures.join("\n")
    );
}

#[test]
fn hard_real_2026_committed_subset_log_column_parity() {
    let cases = load_hard_real_2026_committed_cases();
    let mut failures = Vec::new();
    let mut max_abs_log = 0.0_f64;
    let mut max_label = String::new();

    for (name, case) in cases {
        let trend = fit_deseq2_local_dispersion_trend_from_logs(
            &case.log_means,
            &case.log_disps,
            &case.means,
            MIN_DISP,
        )
        .unwrap_or_else(|error| panic!("failed to fit committed 2026 case {name}: {error}"));

        for prediction in case.predictions {
            let actual_log = trend
                .predict_log_dispersion_one(prediction.log_mean)
                .unwrap();
            let abs_log = log_relative_error(actual_log, prediction.expected_log);
            if abs_log > max_abs_log {
                max_abs_log = abs_log;
                max_label = format!("{name}/{}", prediction.label);
            }
            if abs_log > HARD_REAL_2026_COMMITTED_LOG_TOLERANCE {
                failures.push(format!(
                    "{name}/{}: mean={} actual_log={} expected_log={} abs_log={} tolerance={}",
                    prediction.label,
                    prediction.mean,
                    actual_log,
                    prediction.expected_log,
                    abs_log,
                    HARD_REAL_2026_COMMITTED_LOG_TOLERANCE
                ));
            }
        }
    }

    eprintln!("committed 2026 hard-real subset max log error: {max_abs_log} at {max_label}");
    assert!(
        failures.is_empty(),
        "committed 2026 hard-real subset parity failures:\n{}",
        failures.join("\n")
    );
}

#[test]
fn hard_real_2026_committed_subset_public_wrapper_parity() {
    let cases = load_hard_real_2026_committed_cases();
    let mut failures = Vec::new();
    let mut max_rel = 0.0_f64;
    let mut max_label = String::new();

    for (name, case) in cases {
        let trend = fit_deseq2_local_dispersion_trend(&case.means, &case.disps, MIN_DISP)
            .unwrap_or_else(|error| panic!("failed to fit committed 2026 case {name}: {error}"));

        for prediction in case.predictions {
            let actual = trend.predict_one(prediction.mean).unwrap();
            let rel = relative_error(actual, prediction.expected);
            if rel > max_rel {
                max_rel = rel;
                max_label = format!("{name}/{}", prediction.label);
            }
            if rel > HARD_REAL_2026_COMMITTED_PUBLIC_REL_TOLERANCE {
                failures.push(format!(
                    "{name}/{}: mean={} actual={} expected={} rel={} tolerance={}",
                    prediction.label,
                    prediction.mean,
                    actual,
                    prediction.expected,
                    rel,
                    HARD_REAL_2026_COMMITTED_PUBLIC_REL_TOLERANCE
                ));
            }
        }
    }

    eprintln!("committed 2026 hard-real subset max relative error: {max_rel} at {max_label}");
    assert!(
        failures.is_empty(),
        "committed 2026 hard-real subset public-wrapper parity failures:\n{}",
        failures.join("\n")
    );
}

#[test]
#[ignore = "uses ignored /data real DESeq2 debug tables as a local compatibility diagnostic"]
fn full_real_deseq2_hard_rows_diagnostic() {
    if !Path::new(REAL_DESEQ2_FULL_FIT_POINTS_PATH).exists()
        || !Path::new(REAL_DESEQ2_FULL_HARD_ROWS_PATH).exists()
    {
        eprintln!(
            "Skipping full real DESeq2 diagnostic; expected {REAL_DESEQ2_FULL_FIT_POINTS_PATH} \
             and {REAL_DESEQ2_FULL_HARD_ROWS_PATH}"
        );
        return;
    }

    let (means, disps) = load_full_real_fit_points();
    let predictions = load_full_real_hard_rows();
    let trend = fit_deseq2_local_dispersion_trend(&means, &disps, MIN_DISP).unwrap();

    let mut failures = Vec::new();
    let mut max_rel = 0.0_f64;
    for prediction in predictions {
        let actual = trend.predict_one(prediction.mean).unwrap();
        let rel = relative_error(actual, prediction.expected);
        max_rel = max_rel.max(rel);
        if rel > FULL_REAL_HARD_ROWS_REL_TOLERANCE {
            failures.push(format!(
                "{}: mean={} actual={} source_disp_fit={} rel={} tolerance={}",
                prediction.gene,
                prediction.mean,
                actual,
                prediction.expected,
                rel,
                FULL_REAL_HARD_ROWS_REL_TOLERANCE
            ));
        }
    }

    eprintln!("full real hard-row max relative error: {max_rel}");
    assert!(
        failures.is_empty(),
        "full real DESeq2 hard-row diagnostic failures:\n{}",
        failures.join("\n")
    );
}

#[test]
#[ignore = "uses ignored /data full real DESeq2 table as a local precision diagnostic"]
fn full_real_deseq2_all_rows_diagnostic() {
    if !Path::new(REAL_DESEQ2_FULL_FIT_POINTS_PATH).exists()
        || !Path::new(REAL_DESEQ2_FULL_ALL_ROWS_PATH).exists()
    {
        eprintln!(
            "Skipping full real DESeq2 all-row diagnostic; expected \
             {REAL_DESEQ2_FULL_FIT_POINTS_PATH} and {REAL_DESEQ2_FULL_ALL_ROWS_PATH}"
        );
        return;
    }

    let (means, disps) = load_full_real_fit_points();
    let predictions = match load_full_real_all_rows_predictions() {
        Ok(predictions) => predictions,
        Err(message) => {
            eprintln!("Skipping full real DESeq2 all-row diagnostic: {message}");
            return;
        }
    };
    let trend = fit_deseq2_local_dispersion_trend(&means, &disps, MIN_DISP).unwrap();

    let mut failures = Vec::new();
    let mut max_rel = 0.0_f64;
    for prediction in predictions {
        let actual = trend.predict_one(prediction.mean).unwrap();
        let rel = relative_error(actual, prediction.expected);
        max_rel = max_rel.max(rel);
        if rel > FULL_REAL_ALL_ROWS_REL_TOLERANCE {
            failures.push(format!(
                "{}: mean={} actual={} source_disp_fit={} rel={} tolerance={}",
                prediction.gene,
                prediction.mean,
                actual,
                prediction.expected,
                rel,
                FULL_REAL_ALL_ROWS_REL_TOLERANCE
            ));
        }
    }

    eprintln!("full real all-row max relative error: {max_rel}");
    assert!(
        failures.is_empty(),
        "full real DESeq2 all-row diagnostic failures:\n{}",
        failures.into_iter().take(20).collect::<Vec<_>>().join("\n")
    );
}

#[test]
#[ignore = "uses ignored /data 2026 real DESeq2 bundle as a local precision diagnostic"]
fn hard_real_2026_global_hardest_diagnostic() {
    run_hard_real_2026_log_hard_rows_diagnostic(
        HARD_REAL_2026_GLOBAL_HARDEST_PATH,
        "global hardest rows",
    );
}

#[test]
#[ignore = "uses ignored /data 2026 real DESeq2 bundle as a local precision diagnostic"]
fn hard_real_2026_global_hard_rows_diagnostic() {
    run_hard_real_2026_public_hard_rows_diagnostic(
        HARD_REAL_2026_GLOBAL_HARD_ROWS_PATH,
        "global hard rows",
    );
}

#[test]
#[ignore = "uses ignored /data 2026 real DESeq2 prediction grids as a local precision diagnostic"]
fn hard_real_2026_prediction_grid_diagnostic() {
    if !Path::new(HARD_REAL_2026_ROOT).exists() {
        eprintln!("Skipping 2026 prediction-grid diagnostic; expected {HARD_REAL_2026_ROOT}");
        return;
    }

    // The grid files contain raw locfit prediction-grid values rather than the
    // final DESeq2 `dispFit` boundary semantics used by the public wrapper.
    // Keep this as a broad local regression guard instead of a tight parity
    // contract.
    let mut failures = Vec::new();
    let mut top_errors = Vec::<(f64, String)>::new();
    let mut max_rel = 0.0_f64;
    let mut max_label = String::new();

    for contrast in load_hard_real_2026_contrasts() {
        let (means, _, log_means, log_disps) = load_hard_real_2026_fit_points(&contrast);
        let trend =
            fit_deseq2_local_dispersion_trend_from_logs(&log_means, &log_disps, &means, MIN_DISP)
                .unwrap_or_else(|error| panic!("failed to fit contrast {contrast}: {error}"));
        let predictions = load_hard_real_2026_prediction_grid(&contrast);

        let mut contrast_max_rel = 0.0_f64;
        for prediction in predictions {
            let actual_log = trend
                .predict_log_dispersion_one(prediction.log_mean)
                .unwrap();
            let actual = actual_log.exp();
            let expected = prediction.expected_log.exp();
            let rel = log_relative_error(actual_log, prediction.expected_log);
            contrast_max_rel = contrast_max_rel.max(rel);
            let label = format!(
                "{contrast}/{}: mean={} actual={} expected={} rel={}",
                prediction.label, prediction.mean, actual, expected, rel
            );
            top_errors.push((rel, label.clone()));
            if rel > max_rel {
                max_rel = rel;
                max_label = format!("{contrast}/{}", prediction.label);
            }
            if rel > HARD_REAL_2026_GRID_LOG_REL_TOLERANCE {
                failures.push(format!(
                    "{label} tolerance={}",
                    HARD_REAL_2026_GRID_LOG_REL_TOLERANCE
                ));
            }
        }
        eprintln!("{contrast}: 2026 prediction-grid max relative error {contrast_max_rel}");
    }

    top_errors.sort_by(|a, b| b.0.total_cmp(&a.0));
    eprintln!("2026 prediction-grid max relative error: {max_rel} at {max_label}");
    eprintln!(
        "2026 prediction-grid top errors:\n{}",
        top_errors
            .iter()
            .take(20)
            .map(|(_, label)| label.as_str())
            .collect::<Vec<_>>()
            .join("\n")
    );
    assert!(
        failures.is_empty(),
        "2026 prediction-grid diagnostic failures:\n{}",
        failures.into_iter().take(20).collect::<Vec<_>>().join("\n")
    );
}

fn run_hard_real_2026_log_hard_rows_diagnostic(path: &str, diagnostic_name: &str) {
    if !Path::new(path).exists() {
        eprintln!("Skipping 2026 hard real {diagnostic_name} diagnostic; expected {path}");
        return;
    }

    let cases = load_hard_real_2026_cases(path);
    let mut failures = Vec::new();
    let mut top_errors = Vec::<(f64, String)>::new();
    let mut checked_rows = 0_usize;
    let mut max_rel = 0.0_f64;
    let mut max_label = String::new();

    for (contrast, case) in cases {
        let trend = fit_deseq2_local_dispersion_trend_from_logs(
            &case.log_means,
            &case.log_disps,
            &case.means,
            MIN_DISP,
        )
        .unwrap_or_else(|error| panic!("failed to fit contrast {contrast}: {error}"));
        let mut contrast_max_rel = 0.0_f64;
        for prediction in case.predictions {
            checked_rows += 1;
            let actual_log = trend
                .predict_log_dispersion_one(prediction.log_mean)
                .unwrap();
            let actual = actual_log.exp();
            let expected = prediction.expected_log.exp();
            let rel = log_relative_error(actual_log, prediction.expected_log);
            contrast_max_rel = contrast_max_rel.max(rel);
            let label = format!(
                "{contrast}/{}: mean={} actual={} expected={} rel={}",
                prediction.label, prediction.mean, actual, expected, rel
            );
            top_errors.push((rel, label.clone()));
            if rel > max_rel {
                max_rel = rel;
                max_label = format!("{contrast}/{}", prediction.label);
            }
            if rel > HARD_REAL_2026_LOG_REL_TOLERANCE {
                failures.push(format!(
                    "{label} tolerance={}",
                    HARD_REAL_2026_LOG_REL_TOLERANCE
                ));
            }
        }
        eprintln!("{contrast}: 2026 {diagnostic_name} max relative error {contrast_max_rel}");
    }

    top_errors.sort_by(|a, b| b.0.total_cmp(&a.0));
    eprintln!("2026 hard-real {diagnostic_name} checked rows: {checked_rows}");
    eprintln!("2026 hard-real {diagnostic_name} max relative error: {max_rel} at {max_label}");
    eprintln!(
        "2026 hard-real {diagnostic_name} top errors:\n{}",
        top_errors
            .iter()
            .take(20)
            .map(|(_, label)| label.as_str())
            .collect::<Vec<_>>()
            .join("\n")
    );
    assert!(
        failures.is_empty(),
        "2026 hard real {diagnostic_name} diagnostic failures:\n{}",
        failures.into_iter().take(20).collect::<Vec<_>>().join("\n")
    );
}

fn run_hard_real_2026_public_hard_rows_diagnostic(path: &str, diagnostic_name: &str) {
    if !Path::new(path).exists() {
        eprintln!("Skipping 2026 hard real {diagnostic_name} diagnostic; expected {path}");
        return;
    }

    let cases = load_hard_real_2026_cases(path);
    let mut failures = Vec::new();
    let mut top_errors = Vec::<(f64, String)>::new();
    let mut checked_rows = 0_usize;
    let mut max_rel = 0.0_f64;
    let mut max_label = String::new();

    for (contrast, case) in cases {
        let trend = fit_deseq2_local_dispersion_trend(&case.means, &case.disps, MIN_DISP)
            .unwrap_or_else(|error| panic!("failed to fit contrast {contrast}: {error}"));
        let mut contrast_max_rel = 0.0_f64;
        for prediction in case.predictions {
            checked_rows += 1;
            let actual = trend.predict_one(prediction.mean).unwrap();
            let rel = relative_error(actual, prediction.expected);
            contrast_max_rel = contrast_max_rel.max(rel);
            let label = format!(
                "{contrast}/{}: mean={} actual={} expected={} rel={}",
                prediction.label, prediction.mean, actual, prediction.expected, rel
            );
            top_errors.push((rel, label.clone()));
            if rel > max_rel {
                max_rel = rel;
                max_label = format!("{contrast}/{}", prediction.label);
            }
            if rel > HARD_REAL_2026_PUBLIC_REL_TOLERANCE {
                failures.push(format!(
                    "{label} tolerance={}",
                    HARD_REAL_2026_PUBLIC_REL_TOLERANCE
                ));
            }
        }
        eprintln!(
            "{contrast}: 2026 public-wrapper {diagnostic_name} max relative error {contrast_max_rel}"
        );
    }

    top_errors.sort_by(|a, b| b.0.total_cmp(&a.0));
    eprintln!("2026 public-wrapper {diagnostic_name} checked rows: {checked_rows}");
    eprintln!("2026 public-wrapper {diagnostic_name} max relative error: {max_rel} at {max_label}");
    eprintln!(
        "2026 public-wrapper {diagnostic_name} top errors:\n{}",
        top_errors
            .iter()
            .take(20)
            .map(|(_, label)| label.as_str())
            .collect::<Vec<_>>()
            .join("\n")
    );
    assert!(
        failures.is_empty(),
        "2026 public-wrapper {diagnostic_name} diagnostic failures:\n{}",
        failures.into_iter().take(20).collect::<Vec<_>>().join("\n")
    );
}

fn load_cases() -> BTreeMap<String, Case> {
    let Ok(contents) = fs::read_to_string(FIXTURE_PATH) else {
        return BTreeMap::new();
    };

    let mut cases = BTreeMap::<String, Case>::new();
    for (line_index, line) in contents.lines().enumerate() {
        let line = line.trim();
        if line.is_empty() || line.starts_with('#') {
            continue;
        }
        let fields: Vec<_> = line.split(',').collect();
        assert_eq!(
            fields.len(),
            6,
            "invalid fixture line {}: {}",
            line_index + 1,
            line
        );
        let case_name = fields[0].to_string();
        let kind = fields[1];
        let mean = parse_required(fields[3], line_index, "mean");
        let case = cases.entry(case_name).or_default();

        match kind {
            "input" => {
                case.means.push(mean);
                case.disps
                    .push(parse_required(fields[4], line_index, "disp"));
            }
            "predict" => {
                case.predictions.push(Prediction {
                    mean,
                    expected: parse_required(fields[5], line_index, "expected"),
                });
            }
            other => panic!("invalid fixture kind on line {}: {other}", line_index + 1),
        }
    }

    cases
}

fn load_real_cases(path: &str) -> BTreeMap<String, RealCase> {
    let contents = fs::read_to_string(path)
        .unwrap_or_else(|error| panic!("failed to read real fixture {path}: {error}"));

    let mut cases = BTreeMap::<String, RealCase>::new();
    for (line_index, line) in contents.lines().enumerate() {
        let line = line.trim();
        if line.is_empty() || line.starts_with('#') {
            continue;
        }
        let fields = split_csv_line(line);
        assert_eq!(
            fields.len(),
            8,
            "invalid real fixture line {}: {}",
            line_index + 1,
            line
        );
        let case_name = fields[0].to_string();
        let kind = fields[1].as_str();
        let gene = fields[3].to_string();
        let mean = parse_required(&fields[4], line_index, "mean");
        let case = cases.entry(case_name).or_default();

        match kind {
            "input" => {
                case.means.push(mean);
                case.disps
                    .push(parse_required(&fields[5], line_index, "disp"));
            }
            "predict" => {
                case.predictions.push(RealPrediction {
                    gene,
                    mean,
                    expected: parse_required(&fields[6], line_index, "expected"),
                    source_expected: parse_optional(&fields[7], line_index, "source_expected"),
                });
            }
            other => panic!(
                "invalid real fixture kind on line {}: {other}",
                line_index + 1
            ),
        }
    }

    cases
}

fn load_hard_real_2026_committed_cases() -> BTreeMap<String, HardRealCase> {
    let contents =
        fs::read_to_string(HARD_REAL_2026_COMMITTED_FIXTURE_PATH).unwrap_or_else(|error| {
            panic!(
                "failed to read committed 2026 hard-real fixture \
                 {HARD_REAL_2026_COMMITTED_FIXTURE_PATH}: {error}"
            )
        });

    let mut cases = BTreeMap::<String, HardRealCase>::new();
    for (line_index, line) in contents.lines().enumerate() {
        let line = line.trim();
        if line.is_empty() || line.starts_with('#') {
            continue;
        }
        let fields = split_csv_line(line);
        assert_eq!(
            fields.len(),
            10,
            "invalid committed 2026 fixture line {}: {}",
            line_index + 1,
            line
        );
        let case_name = fields[0].to_string();
        let kind = fields[1].as_str();
        let label = fields[3].to_string();
        let mean = parse_required(&fields[4], line_index, "mean");
        let log_mean = parse_required(&fields[6], line_index, "log_mean");
        let case = cases.entry(case_name).or_default();

        match kind {
            "input" => {
                case.means.push(mean);
                case.disps
                    .push(parse_required(&fields[5], line_index, "disp"));
                case.log_means.push(log_mean);
                case.log_disps
                    .push(parse_required(&fields[7], line_index, "log_disp"));
            }
            "predict" => {
                case.predictions.push(HardRealPrediction {
                    label,
                    mean,
                    expected: parse_required(&fields[8], line_index, "expected"),
                    log_mean,
                    expected_log: parse_required(&fields[9], line_index, "expected_log"),
                });
            }
            other => panic!(
                "invalid committed 2026 fixture kind on line {}: {other}",
                line_index + 1
            ),
        }
    }

    cases
}

fn load_full_real_fit_points() -> (Vec<f64>, Vec<f64>) {
    let contents = fs::read_to_string(REAL_DESEQ2_FULL_FIT_POINTS_PATH).unwrap_or_else(|error| {
        panic!("failed to read {REAL_DESEQ2_FULL_FIT_POINTS_PATH}: {error}")
    });
    let mut lines = contents.lines();
    let header = lines
        .next()
        .unwrap_or_else(|| panic!("{REAL_DESEQ2_FULL_FIT_POINTS_PATH} is empty"));
    let columns: Vec<_> = header.split('\t').collect();
    let mean_index = column_index(&columns, "baseMean");
    let disp_index = column_index(&columns, "dispGeneEst");

    let mut means = Vec::new();
    let mut disps = Vec::new();
    for (line_index, line) in lines.enumerate() {
        let fields: Vec<_> = line.split('\t').collect();
        let mean = parse_required(fields[mean_index], line_index + 1, "baseMean");
        let disp = parse_required(fields[disp_index], line_index + 1, "dispGeneEst");
        if mean.is_finite() && mean > 0.0 && disp.is_finite() && disp > 0.0 {
            means.push(mean);
            disps.push(disp);
        }
    }
    (means, disps)
}

fn load_full_real_hard_rows() -> Vec<RealPrediction> {
    let contents = fs::read_to_string(REAL_DESEQ2_FULL_HARD_ROWS_PATH).unwrap_or_else(|error| {
        panic!("failed to read {REAL_DESEQ2_FULL_HARD_ROWS_PATH}: {error}")
    });
    let mut lines = contents.lines();
    let header = lines
        .next()
        .unwrap_or_else(|| panic!("{REAL_DESEQ2_FULL_HARD_ROWS_PATH} is empty"));
    let columns: Vec<_> = header.split('\t').collect();
    let gene_index = column_index(&columns, "gene");
    let mean_index = column_index(&columns, "baseMean");
    let expected_index = column_index(&columns, "dispFit");

    let mut predictions = Vec::new();
    let mut seen_means = HashSet::new();
    for (line_index, line) in lines.enumerate() {
        let fields: Vec<_> = line.split('\t').collect();
        let mean = parse_required(fields[mean_index], line_index + 1, "baseMean");
        let expected = parse_required(fields[expected_index], line_index + 1, "dispFit");
        if !seen_means.insert(mean.to_bits()) {
            continue;
        }
        predictions.push(RealPrediction {
            gene: fields[gene_index].to_string(),
            mean,
            expected,
            source_expected: Some(expected),
        });
    }
    predictions
}

fn load_full_real_all_rows_predictions() -> Result<Vec<RealPrediction>, String> {
    let output = Command::new("gzip")
        .args(["-cd", REAL_DESEQ2_FULL_ALL_ROWS_PATH])
        .output()
        .map_err(|error| format!("failed to run gzip: {error}"))?;
    if !output.status.success() {
        return Err(format!("gzip exited with status {}", output.status));
    }
    let contents = String::from_utf8(output.stdout)
        .map_err(|error| format!("gzip output was not valid UTF-8: {error}"))?;

    let mut lines = contents.lines();
    let header = lines
        .next()
        .ok_or_else(|| format!("{REAL_DESEQ2_FULL_ALL_ROWS_PATH} is empty"))?;
    let columns: Vec<_> = header.split('\t').collect();
    let gene_index = column_index(&columns, "gene");
    let mean_index = column_index(&columns, "baseMean");
    let expected_index = column_index(&columns, "dispFit");

    let mut predictions = Vec::new();
    let mut seen_means = HashSet::new();
    for (line_index, line) in lines.enumerate() {
        let fields: Vec<_> = line.split('\t').collect();
        if fields.get(mean_index) == Some(&"NA") || fields.get(expected_index) == Some(&"NA") {
            continue;
        }
        let Some(mean_field) = fields.get(mean_index) else {
            continue;
        };
        let Some(expected_field) = fields.get(expected_index) else {
            continue;
        };
        let Ok(mean) = mean_field.parse::<f64>() else {
            continue;
        };
        let Ok(expected) = expected_field.parse::<f64>() else {
            continue;
        };
        if !mean.is_finite() || mean <= 0.0 || !expected.is_finite() || expected <= 0.0 {
            continue;
        }
        if !seen_means.insert(mean.to_bits()) {
            continue;
        }
        let gene = fields.get(gene_index).map_or_else(
            || format!("line_{}", line_index + 2),
            |value| value.to_string(),
        );
        predictions.push(RealPrediction {
            gene,
            mean,
            expected,
            source_expected: Some(expected),
        });
    }
    Ok(predictions)
}

fn load_hard_real_2026_cases(prediction_path: &str) -> BTreeMap<String, HardRealCase> {
    let predictions_by_contrast = load_hard_real_2026_predictions(prediction_path);
    let mut cases = BTreeMap::new();
    for (contrast, predictions) in predictions_by_contrast {
        let (means, disps, log_means, log_disps) = load_hard_real_2026_fit_points(&contrast);
        cases.insert(
            contrast,
            HardRealCase {
                means,
                disps,
                log_means,
                log_disps,
                predictions,
            },
        );
    }
    cases
}

fn load_hard_real_2026_contrasts() -> Vec<String> {
    let path = format!("{HARD_REAL_2026_ROOT}/completed_contrasts.tsv");
    let contents =
        fs::read_to_string(&path).unwrap_or_else(|error| panic!("failed to read {path}: {error}"));
    let mut lines = contents.lines();
    let header = lines.next().unwrap_or_else(|| panic!("{path} is empty"));
    let columns: Vec<_> = header.split('\t').collect();
    let contrast_index = column_index(&columns, "contrast");

    lines
        .filter_map(|line| {
            let fields: Vec<_> = line.split('\t').collect();
            fields
                .get(contrast_index)
                .map(|contrast| contrast.to_string())
        })
        .collect()
}

fn load_hard_real_2026_fit_points(contrast: &str) -> (Vec<f64>, Vec<f64>, Vec<f64>, Vec<f64>) {
    let path = format!("{HARD_REAL_2026_ROOT}/contrasts/{contrast}/locfit/fit_points.tsv");
    let contents =
        fs::read_to_string(&path).unwrap_or_else(|error| panic!("failed to read {path}: {error}"));
    let mut lines = contents.lines();
    let header = lines.next().unwrap_or_else(|| panic!("{path} is empty"));
    let columns: Vec<_> = header.split('\t').collect();
    let mean_index = column_index(&columns, "baseMean");
    let disp_index = column_index(&columns, "dispGeneEst");
    let log_mean_index = column_index(&columns, "logBaseMean");
    let log_disp_index = column_index(&columns, "logDispGeneEst");

    let mut means = Vec::new();
    let mut disps = Vec::new();
    let mut log_means = Vec::new();
    let mut log_disps = Vec::new();
    for (line_index, line) in lines.enumerate() {
        let fields: Vec<_> = line.split('\t').collect();
        let mean = parse_required(fields[mean_index], line_index + 1, "baseMean");
        let disp = parse_required(fields[disp_index], line_index + 1, "dispGeneEst");
        let log_mean = parse_required(fields[log_mean_index], line_index + 1, "logBaseMean");
        let log_disp = parse_required(fields[log_disp_index], line_index + 1, "logDispGeneEst");
        if mean.is_finite()
            && mean > 0.0
            && disp.is_finite()
            && disp > 0.0
            && log_mean.is_finite()
            && log_disp.is_finite()
        {
            means.push(mean);
            disps.push(disp);
            log_means.push(log_mean);
            log_disps.push(log_disp);
        }
    }
    (means, disps, log_means, log_disps)
}

fn load_hard_real_2026_predictions(path: &str) -> BTreeMap<String, Vec<HardRealPrediction>> {
    let contents =
        fs::read_to_string(path).unwrap_or_else(|error| panic!("failed to read {path}: {error}"));
    let mut lines = contents.lines();
    let header = lines.next().unwrap_or_else(|| panic!("{path} is empty"));
    let columns: Vec<_> = header.split('\t').collect();
    let contrast_index = column_index(&columns, "contrast");
    let gene_index = column_index(&columns, "gene");
    let mean_index = column_index(&columns, "baseMean");
    let expected_index = column_index(&columns, "dispFit");
    let log_mean_index = column_index(&columns, "logBaseMean");
    let expected_log_index = column_index(&columns, "logDispFit");

    let mut predictions_by_contrast = BTreeMap::<String, Vec<HardRealPrediction>>::new();
    for (line_index, line) in lines.enumerate() {
        let fields: Vec<_> = line.split('\t').collect();
        let contrast = fields[contrast_index].to_string();
        let mean = parse_required(fields[mean_index], line_index + 1, "baseMean");
        let expected = parse_required(fields[expected_index], line_index + 1, "dispFit");
        let log_mean = parse_required(fields[log_mean_index], line_index + 1, "logBaseMean");
        let expected_log = parse_required(fields[expected_log_index], line_index + 1, "logDispFit");
        if !mean.is_finite()
            || mean <= 0.0
            || !expected.is_finite()
            || expected <= 0.0
            || !log_mean.is_finite()
            || !expected_log.is_finite()
        {
            continue;
        }
        predictions_by_contrast
            .entry(contrast)
            .or_default()
            .push(HardRealPrediction {
                label: fields[gene_index].to_string(),
                mean,
                expected,
                log_mean,
                expected_log,
            });
    }
    predictions_by_contrast
}

fn load_hard_real_2026_prediction_grid(contrast: &str) -> Vec<HardRealPrediction> {
    let path = format!("{HARD_REAL_2026_ROOT}/contrasts/{contrast}/locfit/prediction_grid.tsv");
    let contents =
        fs::read_to_string(&path).unwrap_or_else(|error| panic!("failed to read {path}: {error}"));
    let mut lines = contents.lines();
    let header = lines.next().unwrap_or_else(|| panic!("{path} is empty"));
    let columns: Vec<_> = header.split('\t').collect();
    let kind_index = column_index(&columns, "gridKind");
    let mean_index = column_index(&columns, "baseMean");
    let log_mean_index = column_index(&columns, "logBaseMean");
    let expected_log_index = column_index(&columns, "logDirectLocalDispFit");

    let mut predictions = Vec::new();
    let mut seen = HashSet::new();
    for (line_index, line) in lines.enumerate() {
        let fields: Vec<_> = line.split('\t').collect();
        let mean = parse_required(fields[mean_index], line_index + 1, "baseMean");
        let log_mean = parse_required(fields[log_mean_index], line_index + 1, "logBaseMean");
        let expected_log = parse_required(
            fields[expected_log_index],
            line_index + 1,
            "logDirectLocalDispFit",
        );
        if !mean.is_finite() || mean <= 0.0 || !log_mean.is_finite() || !expected_log.is_finite() {
            continue;
        }
        if !seen.insert((fields[kind_index].to_string(), mean.to_bits())) {
            continue;
        }
        predictions.push(HardRealPrediction {
            label: fields[kind_index].to_string(),
            mean,
            expected: expected_log.exp(),
            log_mean,
            expected_log,
        });
    }
    predictions
}

fn column_index(columns: &[&str], name: &str) -> usize {
    columns
        .iter()
        .position(|column| *column == name)
        .unwrap_or_else(|| panic!("missing required column {name}"))
}

fn split_csv_line(line: &str) -> Vec<String> {
    let mut fields = Vec::new();
    let mut field = String::new();
    let mut chars = line.chars().peekable();
    let mut in_quotes = false;

    while let Some(ch) = chars.next() {
        match ch {
            '"' if in_quotes && chars.peek() == Some(&'"') => {
                field.push('"');
                chars.next();
            }
            '"' => in_quotes = !in_quotes,
            ',' if !in_quotes => {
                fields.push(std::mem::take(&mut field));
            }
            _ => field.push(ch),
        }
    }
    fields.push(field);
    fields
}

fn parse_required(value: &str, line_index: usize, field_name: &str) -> f64 {
    value.parse::<f64>().unwrap_or_else(|error| {
        panic!(
            "invalid {field_name} field on line {}: {} ({error})",
            line_index + 1,
            value
        )
    })
}

fn parse_optional(value: &str, line_index: usize, field_name: &str) -> Option<f64> {
    if value.is_empty() {
        return None;
    }
    Some(parse_required(value, line_index, field_name))
}

fn relative_error(actual: f64, expected: f64) -> f64 {
    if expected == 0.0 {
        actual.abs()
    } else {
        ((actual - expected) / expected).abs()
    }
}

fn log_relative_error(actual_log: f64, expected_log: f64) -> f64 {
    (actual_log - expected_log).abs()
}

fn synthetic_case_tolerance(case_name: &str) -> f64 {
    match case_name {
        "few_points" => FEW_POINTS_REL_TOLERANCE,
        "filtered_min_disp" => FILTERED_MIN_DISP_REL_TOLERANCE,
        _ => DEFAULT_SYNTHETIC_REL_TOLERANCE,
    }
}

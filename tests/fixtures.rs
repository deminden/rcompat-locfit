use std::collections::{BTreeMap, HashSet};
use std::fs;
use std::path::Path;
use std::process::Command;

use rcompat_locfit::fit_deseq2_local_dispersion_trend;

const FIXTURE_PATH: &str = "fixtures/locfit_deseq2_cases.csv";
const REAL_DESEQ2_FIXTURE_PATH: &str = "fixtures/real_deseq2_subset_cases.csv";
const REAL_DESEQ2_FULL_FIT_POINTS_PATH: &str = "data/local_dispersion_fit_points.tsv";
const REAL_DESEQ2_FULL_HARD_ROWS_PATH: &str = "data/local_dispersion_ranked_hard_rows.tsv";
const REAL_DESEQ2_FULL_ALL_ROWS_PATH: &str = "data/local_dispersion_all_rows.tsv.gz";
const MIN_DISP: f64 = 1e-8;
const DEFAULT_SYNTHETIC_REL_TOLERANCE: f64 = 2e-9;
const FEW_POINTS_REL_TOLERANCE: f64 = 2e-9;
const FILTERED_MIN_DISP_REL_TOLERANCE: f64 = 2e-9;
const REAL_DESEQ2_COMMITTED_REL_TOLERANCE: f64 = 1e-9;
const FULL_REAL_HARD_ROWS_REL_TOLERANCE: f64 = 1e-10;
const FULL_REAL_ALL_ROWS_REL_TOLERANCE: f64 = 1e-10;

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

fn synthetic_case_tolerance(case_name: &str) -> f64 {
    match case_name {
        "few_points" => FEW_POINTS_REL_TOLERANCE,
        "filtered_min_disp" => FILTERED_MIN_DISP_REL_TOLERANCE,
        _ => DEFAULT_SYNTHETIC_REL_TOLERANCE,
    }
}

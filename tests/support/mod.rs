#![allow(dead_code)]

use std::collections::BTreeMap;
use std::fs;
use std::path::Path;

use rcompat_locfit::{
    LocfitError, fit_deseq2_local_dispersion_trend, fit_deseq2_local_dispersion_trend_from_logs,
};

pub const HARD_FIXTURE: &str = "fixtures/deseq2_r462_hard_cases.csv";
pub const HARD_QUERIES: &str = "fixtures/deseq2_r462_hard_queries.csv";
pub const BASELINE: &str = "fixtures/deseq2_r462_baseline.csv";

#[derive(Default)]
pub struct Case {
    pub means: Vec<f64>,
    pub disps: Vec<f64>,
    pub xs: Vec<f64>,
    pub ys: Vec<f64>,
    pub predictions: Vec<Prediction>,
}

pub struct Prediction {
    pub gene: String,
    pub mean: f64,
    pub x: f64,
    pub expected_log: f64,
}

pub struct ErrorRow {
    pub case: String,
    pub gene: String,
    pub expected_log: f64,
    pub log_error: f64,
    pub normal_log_error: f64,
}

pub fn decode(hex: &str) -> f64 {
    assert_eq!(hex.len(), 16, "expected one binary64 value");
    let mut bytes = [0; 8];
    for (index, byte) in bytes.iter_mut().enumerate() {
        *byte = u8::from_str_radix(&hex[2 * index..2 * index + 2], 16).unwrap();
    }
    f64::from_le_bytes(bytes)
}

pub fn load(path: impl AsRef<Path>) -> BTreeMap<String, Case> {
    let source = fs::read_to_string(path).unwrap();
    assert!(
        source.contains("# r_version=4.6.2\n"),
        "R 4.6.2 reference required"
    );
    let mut cases: BTreeMap<String, Case> = BTreeMap::new();
    for line in source
        .lines()
        .filter(|line| !line.starts_with('#') && !line.starts_with("case,") && !line.is_empty())
    {
        let fields: Vec<_> = line.split(',').collect();
        assert_eq!(fields.len(), 8);
        let case = cases.entry(fields[0].into()).or_default();
        let mean = decode(fields[3]);
        let x = decode(fields[5]);
        assert!(mean.is_finite() && mean > 0.0 && x.is_finite());
        match fields[1] {
            "input" => {
                let disp = decode(fields[4]);
                let y = decode(fields[6]);
                assert!(disp.is_finite() && disp > 0.0 && y.is_finite());
                case.means.push(mean);
                case.disps.push(disp);
                case.xs.push(x);
                case.ys.push(y);
            }
            "predict" => {
                let expected_log = decode(fields[7]);
                assert!(expected_log.is_finite());
                case.predictions.push(Prediction {
                    gene: fields[2].into(),
                    mean,
                    x,
                    expected_log,
                });
            }
            kind => panic!("unknown fixture row kind: {kind}"),
        }
    }
    cases
}

pub fn evaluate(cases: &BTreeMap<String, Case>) -> Result<Vec<ErrorRow>, LocfitError> {
    let mut rows = Vec::new();
    for (name, case) in cases {
        let logs =
            fit_deseq2_local_dispersion_trend_from_logs(&case.xs, &case.ys, &case.means, 1e-8)?;
        let normal = fit_deseq2_local_dispersion_trend(&case.means, &case.disps, 1e-8)?;
        for prediction in &case.predictions {
            rows.push(ErrorRow {
                case: name.clone(),
                gene: prediction.gene.clone(),
                expected_log: prediction.expected_log,
                log_error: (logs.predict_log_dispersion_one(prediction.x)?
                    - prediction.expected_log)
                    .abs(),
                normal_log_error: (normal.predict_one(prediction.mean)?.ln()
                    - prediction.expected_log)
                    .abs(),
            });
        }
    }
    Ok(rows)
}

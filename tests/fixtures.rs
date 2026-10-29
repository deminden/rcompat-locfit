mod support;

use std::collections::{BTreeMap, BTreeSet};
use std::fs;
use support::{BASELINE, HARD_FIXTURE, HARD_QUERIES, decode, evaluate, load};

#[test]
fn r462_hard_fixture_integrity() {
    let cases = load(HARD_FIXTURE);
    let manifest = fs::read_to_string(HARD_QUERIES).unwrap();
    let mut queries = BTreeMap::new();
    let mut hard_counts: BTreeMap<&str, usize> = BTreeMap::new();
    for line in manifest
        .lines()
        .filter(|line| !line.starts_with('#') && !line.starts_with("case,"))
    {
        let fields: Vec<_> = line.split(',').collect();
        assert_eq!(fields.len(), 3);
        assert!(queries.insert((fields[0], fields[1]), fields[2]).is_none());
        if fields[2] == "hard" {
            *hard_counts.entry(fields[0]).or_default() += 1;
        }
    }
    assert_eq!(cases.len(), 3);
    for (name, inputs, predictions) in [
        ("airway_512", 213, 44),
        ("simulated_600", 579, 45),
        ("airway_full", 8942, 44),
    ] {
        let case = &cases[name];
        assert_eq!(case.means.len(), inputs);
        assert_eq!(case.predictions.len(), predictions);
        assert_eq!(hard_counts[name], 32);
        let mut seen = BTreeSet::new();
        for query in &case.predictions {
            assert!(seen.insert(&query.gene));
            assert!(queries.contains_key(&(name, query.gene.as_str())));
        }
    }
    assert_eq!(queries.len(), 133);
}

#[test]
fn r462_hard_predictions_meet_precision_gates() {
    let rows = evaluate(&load(HARD_FIXTURE)).unwrap();
    let manifest = fs::read_to_string(HARD_QUERIES).unwrap();
    let hard: BTreeSet<_> = manifest
        .lines()
        .filter(|line| line.ends_with(",hard"))
        .map(|line| {
            let fields: Vec<_> = line.split(',').collect();
            (fields[0].to_owned(), fields[1].to_owned())
        })
        .collect();
    let source = fs::read_to_string(BASELINE).unwrap();
    let mut baseline = BTreeMap::new();
    for line in source
        .lines()
        .filter(|line| !line.starts_with('#') && !line.starts_with("case,"))
    {
        let fields: Vec<_> = line.split(',').collect();
        assert_eq!(fields.len(), 5);
        assert!(
            baseline
                .insert(
                    (fields[0].to_owned(), fields[1].to_owned()),
                    (
                        decode(fields[2]),
                        fields[3].parse::<f64>().unwrap(),
                        fields[4].parse::<f64>().unwrap()
                    )
                )
                .is_none()
        );
    }
    assert_eq!(rows.len(), baseline.len());
    let mut failures = Vec::new();
    let mut maxima: BTreeMap<&str, (f64, f64, f64, f64)> = BTreeMap::new();
    let mut checked = 0;
    let mut worst_hard_ratio = 0.0_f64;
    for row in &rows {
        let key = (row.case.clone(), row.gene.clone());
        let &(expected, old_log, old_normal) = &baseline[&key];
        assert_eq!(
            row.expected_log.to_bits(),
            expected.to_bits(),
            "oracle changed for {key:?}"
        );
        assert!(row.log_error.is_finite() && row.normal_log_error.is_finite());
        // Round the measured 0.2.0 maxima upward, allowing a small margin for
        // platform rounding while protecting today's precision independently
        // of the historical tenfold-improvement gate below.
        let absolute_limit = match row.case.as_str() {
            "airway_512" => 1e-11,
            "airway_full" => 3e-13,
            "simulated_600" => 9e-13,
            name => panic!("missing absolute error limit for {name}"),
        };
        if row.log_error > absolute_limit || row.normal_log_error > absolute_limit {
            failures.push(format!(
                "{}/{}: log {:.17e}, normal {:.17e} exceed absolute log-error limit {:.17e}",
                row.case, row.gene, row.log_error, row.normal_log_error, absolute_limit
            ));
        }
        let m = maxima.entry(&row.case).or_default();
        m.0 = m.0.max(old_log);
        m.1 = m.1.max(old_normal);
        m.2 = m.2.max(row.log_error);
        m.3 = m.3.max(row.normal_log_error);
        if hard.contains(&key) {
            checked += 1;
            assert!(old_log > 0.0 && old_normal > 0.0);
            worst_hard_ratio = worst_hard_ratio
                .max((row.log_error / old_log).max(row.normal_log_error / old_normal));
            if row.log_error > old_log / 10.0 || row.normal_log_error > old_normal / 10.0 {
                failures.push(format!(
                    "{}/{}: log {:.17e} (limit {:.17e}), normal {:.17e} (limit {:.17e})",
                    row.case,
                    row.gene,
                    row.log_error,
                    old_log / 10.0,
                    row.normal_log_error,
                    old_normal / 10.0
                ));
            }
        }
    }
    assert_eq!(checked, 96);
    eprintln!(
        "Minimum improvement across 96 hard predictions and both APIs: {:.6}x",
        1.0 / worst_hard_ratio
    );
    for (name, (old_log, old_normal, log, normal)) in maxima {
        eprintln!(
            "R4.6.2 {name}: log={log:.17e}, normal={normal:.17e}, baseline log={old_log:.17e}"
        );
        if log > old_log / 10.0 || normal > old_normal / 10.0 {
            failures.push(format!(
                "{name}: full query-set maximum has not improved tenfold"
            ));
        }
    }
    assert!(
        failures.is_empty(),
        "R4.6.2 precision failures:\n{}",
        failures.join("\n")
    );
}

#[test]
fn r462_synthetic_shape_parity() {
    let cases = load("fixtures/locfit_r462_cases.csv");
    assert_eq!(cases.len(), 10);
    for row in evaluate(&cases).unwrap() {
        assert!(
            row.log_error <= 2e-8 && row.normal_log_error <= 2e-8,
            "{}/{}: log={}, normal={}",
            row.case,
            row.gene,
            row.log_error,
            row.normal_log_error
        );
    }
}

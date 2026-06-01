use rcompat_locfit::{
    fit_deseq2_local_dispersion_trend, fit_deseq2_local_dispersion_trend_from_logs, LocalFit,
    LocalRegressionConfig, LocfitError, PredictionMethod,
};

fn assert_close(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "actual={actual}, expected={expected}, tolerance={tolerance}"
    );
}

fn assert_relative_close(actual: f64, expected: f64, tolerance: f64) {
    let relative_error = ((actual - expected) / expected).abs();
    assert!(
        relative_error <= tolerance,
        "actual={actual}, expected={expected}, relative_error={relative_error}, tolerance={tolerance}"
    );
}

#[test]
fn generic_constant_fit_predicts_constant() {
    let x = [-2.0, -1.0, 0.0, 1.0, 2.0];
    let y = [7.25; 5];
    let fit = LocalFit::fit(&x, &y, None, LocalRegressionConfig::default()).unwrap();
    assert_close(fit.predict_one(0.3).unwrap(), 7.25, 1e-10);
}

#[test]
fn generic_linear_fit_predicts_line() {
    let x = [-3.0, -2.0, -1.0, 0.0, 1.0, 2.0, 3.0];
    let y: Vec<_> = x.iter().map(|&x| 2.0 * x + 1.5).collect();
    let config = LocalRegressionConfig {
        alpha: 1.0,
        degree: 1,
        ..LocalRegressionConfig::default()
    };
    let fit = LocalFit::fit(&x, &y, None, config).unwrap();
    assert_close(fit.predict_one(0.75).unwrap(), 3.0, 1e-10);
}

#[test]
fn generic_linear_fit_reports_local_slope() {
    let x = [-3.0, -2.0, -1.0, 0.0, 1.0, 2.0, 3.0];
    let y: Vec<_> = x.iter().map(|&x| 2.0 * x + 1.5).collect();
    let config = LocalRegressionConfig {
        alpha: 1.0,
        degree: 1,
        ..LocalRegressionConfig::default()
    };
    let fit = LocalFit::fit(&x, &y, None, config).unwrap();
    let (value, slope) = fit.predict_one_with_derivative(0.75).unwrap();
    assert_close(value, 3.0, 1e-10);
    assert_close(slope, 2.0, 1e-10);
}

#[test]
fn hermite_prediction_preserves_exact_line() {
    let x: Vec<_> = (0..120).map(|index| index as f64 / 10.0).collect();
    let y: Vec<_> = x.iter().map(|&x| -0.5 * x + 2.0).collect();
    let config = LocalRegressionConfig {
        alpha: 0.7,
        degree: 1,
        prediction_method: PredictionMethod::LocfitHermiteApprox,
        ..LocalRegressionConfig::default()
    };
    let fit = LocalFit::fit(&x, &y, None, config).unwrap();
    assert_close(fit.predict_one(-2.0).unwrap(), 3.0, 1e-10);
    assert_close(fit.predict_one(3.25).unwrap(), 0.375, 1e-10);
    assert_close(fit.predict_one(13.0).unwrap(), -4.5, 1e-10);
}

#[test]
fn generic_quadratic_fit_predicts_quadratic_near_training_range() {
    let x = [-3.0, -2.0, -1.0, 0.0, 1.0, 2.0, 3.0];
    let y: Vec<_> = x.iter().map(|&x| 0.5 * x * x - 2.0 * x + 4.0).collect();
    let config = LocalRegressionConfig {
        alpha: 1.0,
        degree: 2,
        ..LocalRegressionConfig::default()
    };
    let fit = LocalFit::fit(&x, &y, None, config).unwrap();
    assert_close(fit.predict_one(0.5).unwrap(), 3.125, 1e-10);
}

#[test]
fn deseq2_wrapper_returns_positive_predictions() {
    let means = [1.0, 2.0, 5.0, 10.0, 100.0, 1000.0];
    let disps = [0.5, 0.3, 0.2, 0.12, 0.06, 0.03];
    let trend = fit_deseq2_local_dispersion_trend(&means, &disps, 1e-8).unwrap();
    let predictions = trend.predict(&[3.0, 30.0, 300.0]).unwrap();
    assert!(predictions
        .iter()
        .all(|value| value.is_finite() && *value > 0.0));
}

#[test]
fn deseq2_log_space_constructor_matches_normal_scale_constructor() {
    let means = [1.0, 2.0, 5.0, 10.0, 100.0, 1000.0];
    let disps = [0.5, 0.3, 0.2, 0.12, 0.06, 0.03];
    let log_means: Vec<_> = means.iter().map(|mean| f64::ln(*mean)).collect();
    let log_disps: Vec<_> = disps.iter().map(|disp| f64::ln(*disp)).collect();

    let normal_trend = fit_deseq2_local_dispersion_trend(&means, &disps, 1e-8).unwrap();
    let log_trend =
        fit_deseq2_local_dispersion_trend_from_logs(&log_means, &log_disps, &means, 1e-8).unwrap();

    for mean in [3.0, 30.0, 300.0] {
        assert_close(
            log_trend.predict_one(mean).unwrap(),
            normal_trend.predict_one(mean).unwrap(),
            1e-14,
        );
        assert_close(
            log_trend.predict_log_dispersion_one(mean.ln()).unwrap(),
            normal_trend.predict_one(mean).unwrap().ln(),
            1e-14,
        );
    }
}

#[test]
fn deseq2_log_space_constructor_filters_min_disp_points() {
    let means = [1.0, 2.0, 5.0, 10.0];
    let disps = [1e-9, 2e-9, 5e-9, 8e-9];
    let log_means: Vec<_> = means.iter().map(|mean| f64::ln(*mean)).collect();
    let log_disps: Vec<_> = disps.iter().map(|disp| f64::ln(*disp)).collect();

    let trend =
        fit_deseq2_local_dispersion_trend_from_logs(&log_means, &log_disps, &means, 1e-8).unwrap();

    assert_close(trend.predict_one(100.0).unwrap(), 1e-8, 0.0);
    assert_close(
        trend.predict_log_dispersion_one(100.0_f64.ln()).unwrap(),
        1e-8_f64.ln(),
        0.0,
    );
}

#[test]
fn deseq2_log_space_batch_prediction_matches_scalar_prediction() {
    let means = [1.0, 2.0, 5.0, 10.0, 100.0, 1000.0];
    let disps = [0.5, 0.3, 0.2, 0.12, 0.06, 0.03];
    let log_means: Vec<_> = means.iter().map(|mean| f64::ln(*mean)).collect();
    let log_disps: Vec<_> = disps.iter().map(|disp| f64::ln(*disp)).collect();
    let trend =
        fit_deseq2_local_dispersion_trend_from_logs(&log_means, &log_disps, &means, 1e-8).unwrap();

    let query_log_means = [3.0_f64.ln(), 30.0_f64.ln(), 300.0_f64.ln()];
    let batch = trend.predict_log_dispersion(&query_log_means).unwrap();

    for (&log_mean, &batch_prediction) in query_log_means.iter().zip(&batch) {
        assert_close(
            batch_prediction,
            trend.predict_log_dispersion_one(log_mean).unwrap(),
            0.0,
        );
    }
}

#[test]
fn deseq2_log_space_constructor_validates_input_shapes() {
    let log_means = [0.0, 1.0];
    let log_disps = [-1.0];
    let means = [1.0, 2.0];

    let error = fit_deseq2_local_dispersion_trend_from_logs(&log_means, &log_disps, &means, 1e-8)
        .unwrap_err();
    assert_eq!(
        error,
        LocfitError::LengthMismatch {
            x: 2,
            y: 1,
            weights: Some(2)
        }
    );
}

#[test]
fn deseq2_log_space_constructor_rejects_invalid_values() {
    let error = fit_deseq2_local_dispersion_trend_from_logs(
        &[0.0, f64::NAN],
        &[-1.0, -2.0],
        &[1.0, 2.0],
        1e-8,
    )
    .unwrap_err();
    assert!(matches!(error, LocfitError::InvalidInput(_)));

    let error =
        fit_deseq2_local_dispersion_trend_from_logs(&[0.0, 1.0], &[-1.0, -2.0], &[1.0, 0.0], 1e-8)
            .unwrap_err();
    assert!(matches!(error, LocfitError::InvalidInput(_)));
}

#[test]
fn deseq2_log_space_prediction_rejects_non_finite_log_means() {
    let means = [1.0, 2.0, 5.0, 10.0, 100.0, 1000.0];
    let disps = [0.5, 0.3, 0.2, 0.12, 0.06, 0.03];
    let log_means: Vec<_> = means.iter().map(|mean| f64::ln(*mean)).collect();
    let log_disps: Vec<_> = disps.iter().map(|disp| f64::ln(*disp)).collect();
    let trend =
        fit_deseq2_local_dispersion_trend_from_logs(&log_means, &log_disps, &means, 1e-8).unwrap();

    let error = trend.predict_log_dispersion_one(f64::INFINITY).unwrap_err();
    assert!(matches!(error, LocfitError::InvalidInput(_)));
}

#[test]
fn deseq2_wrapper_tracks_generated_edge_extrapolation_fixture() {
    let means = [
        2.0,
        2.6222678748431285,
        3.438144403717148,
        4.507867809469581,
        5.910418470405774,
        7.74935024091226,
        10.160436093826037,
        13.321692581618315,
        17.466523247656863,
        22.900951398765624,
        30.026214578163447,
        39.36838894573223,
        51.61723080836147,
        67.67710306856463,
        88.73374661957222,
        116.34182658748706,
        152.53971718046878,
        199.99999999999991,
    ];
    let disps = [
        0.21946457081379978,
        0.20000793290109392,
        0.18254935065251637,
        0.16688363823612287,
        0.15282668091742993,
        0.14021327121045887,
        0.12889516724038114,
        0.11873935049824415,
        0.10962646251165668,
        0.10144940205806127,
        0.0941120664340345,
        0.08752822198710482,
        0.08162049063576482,
        0.07631944046653105,
        0.07156276972009307,
        0.0672945745761695,
        0.06346469213155108,
        0.06002811084953578,
    ];
    let trend = fit_deseq2_local_dispersion_trend(&means, &disps, 1e-8).unwrap();

    assert_relative_close(trend.predict_one(0.5).unwrap(), 0.36846083780576766, 5e-8);
    assert_relative_close(trend.predict_one(800.0).unwrap(), 0.04697201772137147, 5e-8);
}

#[test]
fn deseq2_wrapper_tracks_generated_few_point_fixture() {
    let means = [2.0, 5.0, 20.0, 100.0];
    let disps = [0.4, 0.25, 0.1, 0.06];
    let trend = fit_deseq2_local_dispersion_trend(&means, &disps, 1e-8).unwrap();

    assert_relative_close(trend.predict_one(2.0).unwrap(), 0.4, 5e-8);
    assert_relative_close(trend.predict_one(10.0).unwrap(), 0.20717049640480428, 5e-8);
    assert_relative_close(trend.predict_one(50.0).unwrap(), 0.07275476727180408, 5e-8);
    assert_relative_close(trend.predict_one(100.0).unwrap(), 0.06, 5e-8);
}

#[test]
fn deseq2_wrapper_tracks_generated_filtered_min_disp_fixture() {
    let means = [1.0, 2.0, 5.0, 10.0, 20.0, 50.0, 100.0, 200.0, 500.0, 1000.0];
    let disps = [2e-9, 5e-9, 8e-9, 0.11, 0.08, 0.05, 0.03, 0.02, 0.015, 0.01];
    let trend = fit_deseq2_local_dispersion_trend(&means, &disps, 1e-8).unwrap();

    assert_relative_close(trend.predict_one(10.0).unwrap(), 0.1100000014273344, 5e-8);
    assert_relative_close(trend.predict_one(50.0).unwrap(), 0.049572198654013806, 5e-8);
    assert_relative_close(
        trend.predict_one(200.0).unwrap(),
        0.020009763405300963,
        5e-8,
    );
    assert_relative_close(
        trend.predict_one(1000.0).unwrap(),
        0.010000000000160977,
        5e-8,
    );
}

#[test]
fn deseq2_wrapper_filters_min_disp_points() {
    let means = [1.0, 2.0, 5.0, 10.0];
    let disps = [1e-9, 2e-9, 5e-9, 8e-9];
    let trend = fit_deseq2_local_dispersion_trend(&means, &disps, 1e-8).unwrap();
    assert_close(trend.predict_one(100.0).unwrap(), 1e-8, 0.0);
}

#[test]
fn deseq2_wrapper_rejects_non_positive_means() {
    let means = [1.0, 0.0, 3.0];
    let disps = [0.5, 0.3, 0.2];
    let error = fit_deseq2_local_dispersion_trend(&means, &disps, 1e-8).unwrap_err();
    assert!(matches!(error, LocfitError::InvalidInput(_)));
}

#[test]
fn length_mismatch_is_error() {
    let x = [1.0, 2.0];
    let y = [1.0];
    let error = LocalFit::fit(&x, &y, None, LocalRegressionConfig::default()).unwrap_err();
    assert_eq!(
        error,
        LocfitError::LengthMismatch {
            x: 2,
            y: 1,
            weights: None
        }
    );
}

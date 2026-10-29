use super::assert_close;
use rcompat_locfit::{
    fit_deseq2_local_dispersion_trend, fit_deseq2_local_dispersion_trend_from_logs,
};

#[test]
fn deseq2_wrapper_returns_positive_predictions() {
    let means = [1.0, 2.0, 5.0, 10.0, 100.0, 1000.0];
    let disps = [0.5, 0.3, 0.2, 0.12, 0.06, 0.03];
    let trend = fit_deseq2_local_dispersion_trend(&means, &disps, 1e-8).unwrap();
    let predictions = trend.predict(&[3.0, 30.0, 300.0]).unwrap();
    assert!(
        predictions
            .iter()
            .all(|value| value.is_finite() && *value > 0.0)
    );
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
fn deseq2_wrapper_filters_min_disp_points() {
    let means = [1.0, 2.0, 5.0, 10.0];
    let disps = [1e-9, 2e-9, 5e-9, 8e-9];
    let trend = fit_deseq2_local_dispersion_trend(&means, &disps, 1e-8).unwrap();
    assert_close(trend.predict_one(100.0).unwrap(), 1e-8, 0.0);
}

use rcompat_locfit::{
    LocalFit, LocalRegressionConfig, LocfitError, fit_deseq2_local_dispersion_trend,
    fit_deseq2_local_dispersion_trend_from_logs,
};

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

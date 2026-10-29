use super::assert_close;
use rcompat_locfit::{LocalFit, LocalRegressionConfig, PredictionMethod};

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

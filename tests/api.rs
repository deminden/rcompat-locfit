//! Public API behavior and validation, grouped in a single integration target.

#[path = "api/deseq2.rs"]
mod deseq2;
#[path = "api/local_regression.rs"]
mod local_regression;
#[path = "api/validation.rs"]
mod validation;

fn assert_close(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "actual={actual}, expected={expected}, tolerance={tolerance}"
    );
}

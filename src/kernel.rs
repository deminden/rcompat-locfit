use crate::config::Kernel;

/// Evaluate a configured kernel at scaled distance `u`.
pub fn evaluate(kernel: Kernel, u: f64) -> f64 {
    match kernel {
        Kernel::Tricube => tricube(u),
    }
}

/// Tricube kernel.
///
/// Returns `(1 - |u|^3)^3` for `|u| < 1` and zero otherwise.
pub fn tricube(u: f64) -> f64 {
    let abs_u = u.abs();
    if !abs_u.is_finite() || abs_u >= 1.0 {
        return 0.0;
    }
    let inner = 1.0 - abs_u * abs_u * abs_u;
    inner * inner * inner
}

#[cfg(test)]
mod tests {
    use super::tricube;

    #[test]
    fn tricube_basics() {
        assert_eq!(tricube(0.0), 1.0);
        assert_eq!(tricube(1.0), 0.0);
        assert_eq!(tricube(-1.0), 0.0);
        assert!(tricube(0.5) > 0.0);
        assert_eq!(tricube(0.25), tricube(-0.25));
    }
}

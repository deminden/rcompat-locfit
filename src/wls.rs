use crate::error::LocfitError;

/// Pivot threshold for tiny dense normal-equation systems.
///
/// The local polynomial degree is capped at 2, so this solver prioritizes
/// deterministic behavior and clear singular-fit handling over generality.
pub(crate) const PIVOT_THRESHOLD: f64 = 1e-12;

/// Fit a weighted polynomial through centered predictor values and return
/// coefficients `[beta_0, beta_1, beta_2]`.
///
/// Coefficients above the requested degree are returned as zero. `beta_0` is
/// the fitted value at the center and `beta_1` is the local slope, which is the
/// quantity needed for Hermite interpolation compatibility work.
pub(crate) fn weighted_polynomial_coefficients(
    z: &[f64],
    y: &[f64],
    weights: &[f64],
    degree: usize,
) -> Result<[f64; 3], LocfitError> {
    if z.len() != y.len() || z.len() != weights.len() {
        return Err(LocfitError::LengthMismatch {
            x: z.len(),
            y: y.len(),
            weights: Some(weights.len()),
        });
    }
    if degree > 2 {
        return Err(LocfitError::InvalidConfig(
            "degree must be 0, 1, or 2".to_string(),
        ));
    }

    let n_params = degree + 1;
    let mut matrix = [[0.0_f64; 3]; 3];
    let mut rhs = [0.0_f64; 3];

    for ((&zi, &yi), &wi) in z.iter().zip(y).zip(weights) {
        if wi <= 0.0 || !wi.is_finite() || !zi.is_finite() || !yi.is_finite() {
            continue;
        }
        let basis = [1.0, zi, zi * zi];
        for row in 0..n_params {
            rhs[row] += wi * basis[row] * yi;
            for col in 0..n_params {
                matrix[row][col] += wi * basis[row] * basis[col];
            }
        }
    }

    gaussian_solve(&mut matrix, &mut rhs, n_params)
}

fn gaussian_solve(
    matrix: &mut [[f64; 3]; 3],
    rhs: &mut [f64; 3],
    n_params: usize,
) -> Result<[f64; 3], LocfitError> {
    for col in 0..n_params {
        let mut pivot = col;
        let mut pivot_abs = matrix[col][col].abs();
        for (row, values) in matrix.iter().enumerate().take(n_params).skip(col + 1) {
            let candidate = values[col].abs();
            if candidate > pivot_abs {
                pivot = row;
                pivot_abs = candidate;
            }
        }

        if pivot_abs <= PIVOT_THRESHOLD {
            return Err(LocfitError::SingularFit);
        }

        if pivot != col {
            matrix.swap(col, pivot);
            rhs.swap(col, pivot);
        }

        let pivot_value = matrix[col][col];
        let pivot_row = matrix[col];
        for row in (col + 1)..n_params {
            let factor = matrix[row][col] / pivot_value;
            matrix[row][col] = 0.0;
            for (entry, &pivot_entry) in matrix[row]
                .iter_mut()
                .zip(pivot_row.iter())
                .take(n_params)
                .skip(col + 1)
            {
                *entry -= factor * pivot_entry;
            }
            rhs[row] -= factor * rhs[col];
        }
    }

    let mut solution = [0.0_f64; 3];
    for row in (0..n_params).rev() {
        let tail: f64 = matrix[row]
            .iter()
            .zip(solution.iter())
            .enumerate()
            .take(n_params)
            .skip(row + 1)
            .map(|(_, (&a, &b))| a * b)
            .sum();
        let diagonal = matrix[row][row];
        if diagonal.abs() <= PIVOT_THRESHOLD {
            return Err(LocfitError::SingularFit);
        }
        solution[row] = (rhs[row] - tail) / diagonal;
    }

    Ok(solution)
}

#[cfg(test)]
mod tests {
    use super::weighted_polynomial_coefficients;
    use crate::LocfitError;

    fn close(actual: f64, expected: f64) {
        assert!(
            (actual - expected).abs() < 1e-10,
            "actual={actual}, expected={expected}"
        );
    }

    #[test]
    fn exact_constant_data() {
        let z = [-2.0, -1.0, 0.0, 1.0, 2.0];
        let y = [3.5; 5];
        let weights = [1.0; 5];
        close(
            weighted_polynomial_coefficients(&z, &y, &weights, 0).unwrap()[0],
            3.5,
        );
    }

    #[test]
    fn exact_line_degree_one() {
        let z = [-2.0, -1.0, 0.0, 1.0, 2.0];
        let y = [-1.0, 1.0, 3.0, 5.0, 7.0];
        let weights = [1.0; 5];
        close(
            weighted_polynomial_coefficients(&z, &y, &weights, 1).unwrap()[0],
            3.0,
        );
    }

    #[test]
    fn exact_line_returns_slope() {
        let z = [-2.0, -1.0, 0.0, 1.0, 2.0];
        let y = [-1.0, 1.0, 3.0, 5.0, 7.0];
        let weights = [1.0; 5];
        let coefficients = weighted_polynomial_coefficients(&z, &y, &weights, 1).unwrap();
        close(coefficients[0], 3.0);
        close(coefficients[1], 2.0);
        close(coefficients[2], 0.0);
    }

    #[test]
    fn exact_quadratic_degree_two() {
        let z = [-2.0, -1.0, 0.0, 1.0, 2.0];
        let y = [9.0, 4.0, 3.0, 6.0, 13.0];
        let weights = [1.0; 5];
        close(
            weighted_polynomial_coefficients(&z, &y, &weights, 2).unwrap()[0],
            3.0,
        );
    }

    #[test]
    fn exact_quadratic_returns_local_slope() {
        let z = [-2.0, -1.0, 0.0, 1.0, 2.0];
        let y = [9.0, 4.0, 3.0, 6.0, 13.0];
        let weights = [1.0; 5];
        let coefficients = weighted_polynomial_coefficients(&z, &y, &weights, 2).unwrap();
        close(coefficients[0], 3.0);
        close(coefficients[1], 1.0);
        close(coefficients[2], 2.0);
    }

    #[test]
    fn weighted_constant_data() {
        let z = [-2.0, -1.0, 0.0, 1.0, 2.0];
        let y = [4.25; 5];
        let weights = [1.0, 10.0, 3.0, 5.0, 2.0];
        close(
            weighted_polynomial_coefficients(&z, &y, &weights, 0).unwrap()[0],
            4.25,
        );
    }

    #[test]
    fn repeated_x_is_singular_for_line() {
        let z = [0.0, 0.0, 0.0];
        let y = [1.0, 2.0, 3.0];
        let weights = [1.0, 1.0, 1.0];
        assert_eq!(
            weighted_polynomial_coefficients(&z, &y, &weights, 1),
            Err(LocfitError::SingularFit)
        );
    }
}

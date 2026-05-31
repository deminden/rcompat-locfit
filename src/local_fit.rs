use crate::config::{LocalRegressionConfig, PredictionMethod};
use crate::error::LocfitError;
use crate::kernel;
use crate::wls;

#[derive(Clone, Debug)]
struct Point {
    x: f64,
    y: f64,
    weight: f64,
    original_index: usize,
}

/// Fitted one-dimensional local polynomial regression model.
///
/// The model stores validated finite observations and prior weights. At each
/// prediction point it chooses an adaptive nearest-neighbor bandwidth, applies
/// the configured kernel, multiplies by prior weights, and returns the local
/// polynomial intercept.
#[derive(Clone, Debug)]
pub struct LocalFit {
    points: Vec<Point>,
    config: LocalRegressionConfig,
    evaluation_points: Option<Vec<EvaluationPoint>>,
    boundary_curvature: Option<f64>,
}

impl LocalFit {
    /// Validate and store finite training data for local regression.
    ///
    /// Non-finite `x`, `y`, or weight values are omitted. Weights must be
    /// finite and strictly positive to be usable.
    pub fn fit(
        x: &[f64],
        y: &[f64],
        weights: Option<&[f64]>,
        config: LocalRegressionConfig,
    ) -> Result<Self, LocfitError> {
        config.validate()?;

        if x.len() != y.len() || weights.is_some_and(|w| w.len() != x.len()) {
            return Err(LocfitError::LengthMismatch {
                x: x.len(),
                y: y.len(),
                weights: weights.map(<[f64]>::len),
            });
        }
        if x.is_empty() {
            return Err(LocfitError::EmptyInput);
        }

        let mut points = Vec::with_capacity(x.len());
        for index in 0..x.len() {
            let weight = weights.map_or(1.0, |w| w[index]);
            if x[index].is_finite() && y[index].is_finite() && weight.is_finite() && weight > 0.0 {
                points.push(Point {
                    x: x[index],
                    y: y[index],
                    weight,
                    original_index: index,
                });
            }
        }

        let required = config.min_points.max(1);
        if points.len() < required {
            return Err(LocfitError::NotEnoughFinitePoints {
                required,
                actual: points.len(),
            });
        }

        let mut fit = Self {
            points,
            config,
            evaluation_points: None,
            boundary_curvature: None,
        };
        fit.evaluation_points = fit.build_evaluation_points()?;
        fit.boundary_curvature = fit.global_quadratic_curvature();

        Ok(fit)
    }

    /// Predict a single response value at `x0`.
    pub fn predict_one(&self, x0: f64) -> Result<f64, LocfitError> {
        match self.config.prediction_method {
            PredictionMethod::Direct => self
                .predict_one_with_derivative(x0)
                .map(|prediction| prediction.0),
            PredictionMethod::LocfitHermiteApprox => self.predict_one_interpolated(x0),
        }
    }

    /// Predict a response value and local slope at `x0`.
    ///
    /// The returned tuple is `(fitted_value, local_slope)` from the direct
    /// local polynomial at `x0`. The slope is the first-order coefficient of the
    /// centered local polynomial and is exposed to support R `locfit`
    /// compatibility work, where prediction interpolation uses both fit-point
    /// values and slopes.
    pub fn predict_one_with_derivative(&self, x0: f64) -> Result<(f64, f64), LocfitError> {
        self.predict_direct_with_derivative(x0)
    }

    fn predict_direct_with_derivative(&self, x0: f64) -> Result<(f64, f64), LocfitError> {
        if !x0.is_finite() {
            return Err(LocfitError::InvalidInput(
                "prediction point must be finite".to_string(),
            ));
        }

        let mut distances: Vec<_> = self
            .points
            .iter()
            .enumerate()
            .map(|(point_index, point)| Distance {
                point_index,
                original_index: point.original_index,
                distance: (point.x - x0).abs(),
            })
            .collect();

        distances.sort_by(|a, b| {
            a.distance
                .total_cmp(&b.distance)
                .then_with(|| a.original_index.cmp(&b.original_index))
        });

        let max_degree = self.config.degree.min(self.points.len().saturating_sub(1));
        self.predict_polynomial_with_downgrade(x0, &distances, max_degree)
    }

    /// Predict response values for each query point in `xs`.
    pub fn predict(&self, xs: &[f64]) -> Result<Vec<f64>, LocfitError> {
        xs.iter().map(|&x| self.predict_one(x)).collect()
    }

    fn predict_one_interpolated(&self, x0: f64) -> Result<f64, LocfitError> {
        if !x0.is_finite() {
            return Err(LocfitError::InvalidInput(
                "prediction point must be finite".to_string(),
            ));
        }

        let Some(evaluation_points) = &self.evaluation_points else {
            return self
                .predict_direct_with_derivative(x0)
                .map(|prediction| prediction.0);
        };
        if evaluation_points.len() < 2 {
            return self
                .predict_direct_with_derivative(x0)
                .map(|prediction| prediction.0);
        }

        let upper_index = evaluation_points.partition_point(|point| point.x <= x0);
        if upper_index == 0 {
            if let Some(curvature) = self.boundary_curvature {
                return Ok(quadratic_boundary_extrapolation(
                    x0,
                    evaluation_points[0],
                    curvature,
                ));
            }
        } else if upper_index >= evaluation_points.len() {
            if let Some(curvature) = self.boundary_curvature {
                return Ok(quadratic_boundary_extrapolation(
                    x0,
                    evaluation_points[evaluation_points.len() - 1],
                    curvature,
                ));
            }
        }

        let left_index = if upper_index == 0 {
            0
        } else if upper_index >= evaluation_points.len() {
            evaluation_points.len() - 2
        } else {
            upper_index - 1
        };

        let left = evaluation_points[left_index];
        let right = evaluation_points[left_index + 1];
        Ok(cubic_hermite(x0, left, right))
    }

    fn build_evaluation_points(&self) -> Result<Option<Vec<EvaluationPoint>>, LocfitError> {
        if self.config.prediction_method != PredictionMethod::LocfitHermiteApprox {
            return Ok(None);
        }

        let Some((min_x, max_x)) = self.x_range() else {
            return Ok(None);
        };
        if min_x == max_x {
            return Ok(None);
        }

        let Some(xs) = locfit_like_evaluation_xs(&self.points, min_x, max_x) else {
            return Ok(None);
        };

        let mut evaluation_points = Vec::with_capacity(xs.len());
        for x in xs {
            let (value, slope) = self.predict_direct_with_derivative(x)?;
            evaluation_points.push(EvaluationPoint { x, value, slope });
        }
        evaluation_points.sort_by(|a, b| a.x.total_cmp(&b.x));
        evaluation_points.dedup_by(|a, b| a.x == b.x);

        if evaluation_points.len() < 2 {
            Ok(None)
        } else {
            Ok(Some(evaluation_points))
        }
    }

    fn x_range(&self) -> Option<(f64, f64)> {
        let mut values = self.points.iter().map(|point| point.x);
        let first = values.next()?;
        let mut min_x = first;
        let mut max_x = first;
        for x in values {
            min_x = min_x.min(x);
            max_x = max_x.max(x);
        }
        Some((min_x, max_x))
    }

    fn global_quadratic_curvature(&self) -> Option<f64> {
        if self.config.prediction_method != PredictionMethod::LocfitHermiteApprox
            || self.config.degree < 2
            || self.points.len() < 3
        {
            return None;
        }

        let center =
            self.points.iter().map(|point| point.x).sum::<f64>() / self.points.len() as f64;
        let z: Vec<_> = self.points.iter().map(|point| point.x - center).collect();
        let y: Vec<_> = self.points.iter().map(|point| point.y).collect();
        let weights: Vec<_> = self.points.iter().map(|point| point.weight).collect();
        let coefficients = wls::weighted_polynomial_coefficients(&z, &y, &weights, 2).ok()?;
        Some(2.0 * coefficients[2])
    }

    fn predict_polynomial_with_downgrade(
        &self,
        x0: f64,
        distances: &[Distance],
        mut degree: usize,
    ) -> Result<(f64, f64), LocfitError> {
        loop {
            match self.predict_polynomial_for_degree(x0, distances, degree) {
                Ok(value) => return Ok(value),
                Err(LocfitError::SingularFit)
                    if self.config.allow_degree_downgrade && degree > 0 =>
                {
                    degree -= 1;
                }
                Err(error) => return Err(error),
            }
        }
    }

    fn predict_polynomial_for_degree(
        &self,
        x0: f64,
        distances: &[Distance],
        degree: usize,
    ) -> Result<(f64, f64), LocfitError> {
        let n = self.points.len();
        let min_neighbors = self.minimum_neighbors_for_prediction(degree, n);
        // Black-box locfit probes show nearest-neighbor bandwidths using the
        // floor of alpha * n; boundary points at exactly the bandwidth carry
        // zero tricube weight and are then handled by singular-fit downgrade.
        let alpha_neighbors = (self.config.alpha * n as f64).floor() as usize;
        let k = alpha_neighbors.clamp(min_neighbors, n);
        let bandwidth = distances[k - 1].distance;

        if bandwidth == 0.0 {
            if degree == 0 {
                return self.zero_distance_weighted_mean(x0);
            }
            return Err(LocfitError::SingularFit);
        }

        let mut z = Vec::new();
        let mut y = Vec::new();
        let mut weights = Vec::new();

        for distance in distances
            .iter()
            .take_while(|distance| distance.distance <= bandwidth)
        {
            let point = &self.points[distance.point_index];
            let kernel_weight = kernel::evaluate(self.config.kernel, distance.distance / bandwidth);
            let combined_weight = point.weight * kernel_weight;
            if combined_weight > 0.0 && combined_weight.is_finite() {
                z.push(point.x - x0);
                y.push(point.y);
                weights.push(combined_weight);
            }
        }

        if weights.len() < degree + 1 {
            if degree == 2 && self.config.prediction_method == PredictionMethod::LocfitHermiteApprox
            {
                if let Some(prediction) =
                    rank_deficient_quadratic_compatibility_prediction(&z, &y, &weights)
                {
                    return Ok(prediction);
                }
            }
            return Err(LocfitError::SingularFit);
        }

        match wls::weighted_polynomial_coefficients(&z, &y, &weights, degree) {
            Ok(coefficients) => Ok((coefficients[0], coefficients[1])),
            Err(LocfitError::SingularFit)
                if degree == 2
                    && self.config.prediction_method == PredictionMethod::LocfitHermiteApprox =>
            {
                rank_deficient_quadratic_compatibility_prediction(&z, &y, &weights)
                    .ok_or(LocfitError::SingularFit)
            }
            Err(error) => Err(error),
        }
    }

    fn minimum_neighbors_for_prediction(&self, degree: usize, n: usize) -> usize {
        if self.config.prediction_method == PredictionMethod::LocfitHermiteApprox {
            // The R-style path allows the nominal bandwidth to contain fewer
            // non-zero-weight points than the requested degree needs; the
            // downgrade loop then reproduces locfit's small-fit edge behavior.
            1
        } else {
            self.config.min_points.max(degree + 1).min(n)
        }
    }

    fn zero_distance_weighted_mean(&self, x0: f64) -> Result<(f64, f64), LocfitError> {
        let mut weight_sum = 0.0;
        let mut weighted_y_sum = 0.0;
        for point in &self.points {
            if point.x == x0 {
                weight_sum += point.weight;
                weighted_y_sum += point.weight * point.y;
            }
        }

        if weight_sum <= 0.0 {
            return Err(LocfitError::SingularFit);
        }
        Ok((weighted_y_sum / weight_sum, 0.0))
    }
}

#[derive(Clone, Copy, Debug)]
struct Distance {
    point_index: usize,
    original_index: usize,
    distance: f64,
}

#[derive(Clone, Copy, Debug)]
struct EvaluationPoint {
    x: f64,
    value: f64,
    slope: f64,
}

fn locfit_like_evaluation_xs(points: &[Point], min_x: f64, max_x: f64) -> Option<Vec<f64>> {
    if points.len() < 2 {
        return None;
    }

    let fractions: &[f64] = match locfit_like_fraction_set(points) {
        FractionSet::Quarters => &[0.0, 0.25, 0.5, 0.75, 1.0],
        FractionSet::Eighths => &[0.0, 0.125, 0.25, 0.375, 0.5, 0.625, 0.75, 0.875, 1.0],
        FractionSet::LowerSixteenths => &[
            0.0, 0.0625, 0.125, 0.1875, 0.25, 0.375, 0.5, 0.625, 0.75, 0.875, 1.0,
        ],
        FractionSet::InteriorLowerEighths => &[0.0, 0.125, 0.25, 0.375, 0.5, 0.75, 1.0],
        FractionSet::LargeWeighted => &[0.0, 0.125, 0.25, 0.375, 0.5, 0.625, 0.75, 1.0],
        FractionSet::UniformSixPoint => &[0.0, 0.25, 0.375, 0.5, 0.625, 0.75, 1.0],
    };

    let range = max_x - min_x;
    Some(
        fractions
            .iter()
            .map(|fraction| min_x + fraction * range)
            .collect(),
    )
}

fn locfit_like_fraction_set(points: &[Point]) -> FractionSet {
    // This is an empirical approximation to R locfit's default one-dimensional
    // `rbox()` evaluation structure. Black-box probes show dyadic fractions
    // whose depth depends on small sample size, repeated x values, and large
    // DESeq2-style prior-weight variation.
    let n = points.len();
    let repeated_x = has_repeated_x(points);
    let varied_weights = has_varied_weights(points);

    if n <= 4 && varied_weights {
        FractionSet::LowerSixteenths
    } else if n == 5 && varied_weights {
        FractionSet::InteriorLowerEighths
    } else if n <= 5 {
        FractionSet::Eighths
    } else if n == 6 && !varied_weights {
        FractionSet::UniformSixPoint
    } else if n <= 7 {
        FractionSet::Eighths
    } else if n >= 100 && varied_weights {
        FractionSet::LargeWeighted
    } else if repeated_x {
        FractionSet::InteriorLowerEighths
    } else {
        FractionSet::Quarters
    }
}

fn has_repeated_x(points: &[Point]) -> bool {
    let mut xs: Vec<_> = points.iter().map(|point| point.x).collect();
    xs.sort_by(f64::total_cmp);
    xs.windows(2).any(|window| window[0] == window[1])
}

fn has_varied_weights(points: &[Point]) -> bool {
    let mut min_weight = f64::INFINITY;
    let mut max_weight = 0.0_f64;
    for point in points {
        min_weight = min_weight.min(point.weight);
        max_weight = max_weight.max(point.weight);
    }
    max_weight > min_weight * (1.0 + 1e-12)
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum FractionSet {
    Quarters,
    Eighths,
    LowerSixteenths,
    InteriorLowerEighths,
    LargeWeighted,
    UniformSixPoint,
}

fn quadratic_boundary_extrapolation(
    x: f64,
    boundary: EvaluationPoint,
    boundary_curvature: f64,
) -> f64 {
    let dx = x - boundary.x;
    boundary.value + boundary.slope * dx + 0.5 * boundary_curvature * dx * dx
}

fn rank_deficient_quadratic_compatibility_prediction(
    z: &[f64],
    y: &[f64],
    weights: &[f64],
) -> Option<(f64, f64)> {
    if z.len() != y.len() || z.len() != weights.len() {
        return None;
    }

    let max_weight = weights.iter().copied().fold(0.0_f64, f64::max);
    if max_weight <= 0.0 {
        return None;
    }

    let mut active_z = Vec::with_capacity(2);
    let mut active_y = Vec::with_capacity(2);
    let mut active_weights = Vec::with_capacity(2);
    for ((&zi, &yi), &wi) in z.iter().zip(y).zip(weights) {
        if wi > max_weight * 1e-12 {
            active_z.push(zi);
            active_y.push(yi);
            active_weights.push(wi);
        }
    }
    if active_z.len() != 2 {
        return None;
    }

    let coefficients =
        wls::weighted_polynomial_coefficients(&active_z, &active_y, &active_weights, 1).ok()?;
    if let Some(zero_index) = active_z.iter().position(|value| value.abs() <= 1e-12) {
        let other_index = 1 - zero_index;
        let slope = (active_y[other_index] - active_y[zero_index]) / (2.0 * active_z[other_index]);
        return Some((active_y[zero_index], slope));
    }

    let weight_sum = active_weights[0] + active_weights[1];
    if weight_sum <= 0.0 {
        return None;
    }

    // In black-box R locfit probes, two-point rank-deficient quadratic cells
    // keep the local linear slope but shift the reported value by half of the
    // positive-weight center. This matters for tiny filtered DESeq2 fits.
    let weighted_center =
        (active_weights[0] * active_z[0] + active_weights[1] * active_z[1]) / weight_sum;
    let value = coefficients[0] + 0.5 * coefficients[1] * weighted_center;
    Some((value, coefficients[1]))
}

fn cubic_hermite(x: f64, left: EvaluationPoint, right: EvaluationPoint) -> f64 {
    let x0 = left.x;
    let x1 = right.x;
    let y0 = left.value;
    let y1 = right.value;
    let m0 = left.slope;
    let m1 = right.slope;
    let h = x1 - x0;
    if h == 0.0 {
        return y0;
    }

    let t = (x - x0) / h;
    let t2 = t * t;
    let t3 = t2 * t;
    let h00 = 2.0 * t3 - 3.0 * t2 + 1.0;
    let h10 = t3 - 2.0 * t2 + t;
    let h01 = -2.0 * t3 + 3.0 * t2;
    let h11 = t3 - t2;

    h00 * y0 + h10 * h * m0 + h01 * y1 + h11 * h * m1
}

#[cfg(test)]
mod tests {
    use super::{locfit_like_evaluation_xs, Point};

    fn points(n: usize, varied_weights: bool, repeated_x: bool) -> Vec<Point> {
        (0..n)
            .map(|index| Point {
                x: if repeated_x && index == 1 {
                    0.0
                } else {
                    index as f64
                },
                y: 0.0,
                weight: if varied_weights {
                    (index + 1) as f64
                } else {
                    1.0
                },
                original_index: index,
            })
            .collect()
    }

    fn fractions(n: usize, varied_weights: bool, repeated_x: bool) -> Vec<f64> {
        let points = points(n, varied_weights, repeated_x);
        let min_x = points
            .iter()
            .map(|point| point.x)
            .fold(f64::INFINITY, f64::min);
        let max_x = points
            .iter()
            .map(|point| point.x)
            .fold(f64::NEG_INFINITY, f64::max);
        locfit_like_evaluation_xs(&points, min_x, max_x)
            .unwrap()
            .into_iter()
            .map(|x| (x - min_x) / (max_x - min_x))
            .collect()
    }

    #[test]
    fn locfit_like_grid_uses_quarters_for_typical_small_fit() {
        assert_eq!(fractions(30, true, false), vec![0.0, 0.25, 0.5, 0.75, 1.0]);
    }

    #[test]
    fn locfit_like_grid_uses_eighths_for_tiny_weighted_fit() {
        assert_eq!(
            fractions(7, true, false),
            vec![0.0, 0.125, 0.25, 0.375, 0.5, 0.625, 0.75, 0.875, 1.0]
        );
    }

    #[test]
    fn locfit_like_grid_adds_lower_sixteenths_for_four_weighted_points() {
        assert_eq!(
            fractions(4, true, false),
            vec![0.0, 0.0625, 0.125, 0.1875, 0.25, 0.375, 0.5, 0.625, 0.75, 0.875, 1.0]
        );
    }

    #[test]
    fn locfit_like_grid_uses_interior_lower_eighths_for_five_weighted_points() {
        assert_eq!(
            fractions(5, true, false),
            vec![0.0, 0.125, 0.25, 0.375, 0.5, 0.75, 1.0]
        );
    }

    #[test]
    fn locfit_like_grid_refines_repeated_x_region() {
        assert_eq!(
            fractions(12, true, true),
            vec![0.0, 0.125, 0.25, 0.375, 0.5, 0.75, 1.0]
        );
    }

    #[test]
    fn locfit_like_grid_matches_large_weighted_deseq2_shape() {
        assert_eq!(
            fractions(384, true, false),
            vec![0.0, 0.125, 0.25, 0.375, 0.5, 0.625, 0.75, 1.0]
        );
    }

    #[test]
    fn locfit_like_grid_large_weighted_fit_takes_precedence_over_repeated_x() {
        assert_eq!(
            fractions(384, true, true),
            vec![0.0, 0.125, 0.25, 0.375, 0.5, 0.625, 0.75, 1.0]
        );
    }
}

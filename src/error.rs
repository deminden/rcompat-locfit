use std::fmt;

/// Errors returned by local regression fitting and prediction.
#[derive(Debug, Clone, PartialEq)]
pub enum LocfitError {
    /// Input slices have incompatible lengths.
    LengthMismatch {
        /// Length of the predictor slice.
        x: usize,
        /// Length of the response slice.
        y: usize,
        /// Length of the optional weight slice.
        weights: Option<usize>,
    },
    /// No observations were provided.
    EmptyInput,
    /// Too few finite, positively weighted points are available.
    NotEnoughFinitePoints {
        /// Number of points required.
        required: usize,
        /// Number of usable points found.
        actual: usize,
    },
    /// Invalid local regression configuration.
    InvalidConfig(String),
    /// Invalid user input.
    InvalidInput(String),
    /// Local weighted least-squares system was singular.
    SingularFit,
}

impl fmt::Display for LocfitError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::LengthMismatch { x, y, weights } => {
                write!(f, "input length mismatch: x={x}, y={y}")?;
                if let Some(weights) = weights {
                    write!(f, ", weights={weights}")?;
                }
                Ok(())
            }
            Self::EmptyInput => write!(f, "input is empty"),
            Self::NotEnoughFinitePoints { required, actual } => write!(
                f,
                "not enough finite positively weighted points: required {required}, got {actual}"
            ),
            Self::InvalidConfig(message) => write!(f, "invalid config: {message}"),
            Self::InvalidInput(message) => write!(f, "invalid input: {message}"),
            Self::SingularFit => write!(f, "local fit is singular"),
        }
    }
}

impl std::error::Error for LocfitError {}

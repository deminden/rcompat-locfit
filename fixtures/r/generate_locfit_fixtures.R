suppressPackageStartupMessages(library(locfit))

out_path <- file.path("fixtures", "locfit_deseq2_cases.csv")
min_disp_default <- 1e-8

write_header <- function(con) {
  writeLines("# generated_by=fixtures/r/generate_locfit_fixtures.R", con)
  writeLines(paste0("# r_version=", paste(R.version$major, R.version$minor, sep = ".")), con)
  writeLines(paste0("# locfit_version=", as.character(packageVersion("locfit"))), con)
  writeLines("# columns:case,kind,index,mean,disp,expected", con)
}

compute_expected <- function(means, disps, pred_means, min_disp) {
  keep <- is.finite(means) &
    is.finite(disps) &
    means > 0 &
    disps > 0 &
    disps >= min_disp * 10

  if (!any(keep)) {
    return(rep(min_disp, length(pred_means)))
  }

  d <- data.frame(
    logDisps = log(disps[keep]),
    logMeans = log(means[keep])
  )

  fit <- locfit(
    logDisps ~ lp(logMeans, nn = 0.7, deg = 2),
    data = d,
    weights = means[keep]
  )

  as.numeric(exp(predict(fit, newdata = data.frame(logMeans = log(pred_means)))))
}

write_case <- function(con, name, means, disps, pred_means, min_disp = min_disp_default) {
  for (i in seq_along(means)) {
    writeLines(
      sprintf("%s,input,%d,%.17g,%.17g,", name, i - 1L, means[[i]], disps[[i]]),
      con
    )
  }

  expected <- compute_expected(means, disps, pred_means, min_disp)
  for (i in seq_along(pred_means)) {
    writeLines(
      sprintf("%s,predict,%d,%.17g,,%.17g", name, i - 1L, pred_means[[i]], expected[[i]]),
      con
    )
  }
}

set.seed(1)
dir.create(dirname(out_path), recursive = TRUE, showWarnings = FALSE)

con <- file(out_path, open = "w")
on.exit(close(con), add = TRUE)

write_header(con)

means <- exp(seq(log(1), log(1000), length.out = 30))
disps <- 0.02 + 0.6 / sqrt(means) + exp(rnorm(length(means), sd = 0.03)) * 0.005
write_case(
  con,
  "smooth_monotone",
  means,
  disps,
  exp(seq(log(1.5), log(800), length.out = 12))
)

means <- exp(seq(log(0.1), log(1e6), length.out = 36))
disps <- 1e-4 + 0.8 / (1 + means ^ 0.35)
disps <- disps * exp(rnorm(length(means), sd = 0.04))
write_case(
  con,
  "wide_dynamic_range",
  means,
  disps,
  c(0.1, 1, 10, 1e3, 1e5, 1e6)
)
write_case(
  con,
  "wide_extrapolation",
  means,
  disps,
  c(0.01, 0.1, 1e3, 1e6, 1e7, 1e8)
)

means <- exp(seq(log(2), log(200), length.out = 18))
disps <- 0.03 + 0.25 / means ^ 0.4
write_case(
  con,
  "edge_predictions",
  means,
  disps,
  c(0.5, min(means), 15, max(means), 800)
)

means <- c(1, 2, 6, 30, 150)
disps <- c(0.62, 0.41, 0.21, 0.085, 0.045)
write_case(
  con,
  "five_weighted_points",
  means,
  disps,
  c(0.5, 1, 3, 15, 80, 300)
)

means <- c(1, 1.8, 3.5, 7, 20, 60, 180)
disps <- c(0.58, 0.44, 0.31, 0.21, 0.105, 0.062, 0.038)
write_case(
  con,
  "seven_weighted_points",
  means,
  disps,
  c(0.7, 1, 2.5, 10, 40, 180, 500)
)

means <- c(1, 1, 1, 3, 3, 10, 30, 30, 100, 300, 300, 1000)
disps <- c(0.55, 0.5, 0.52, 0.34, 0.31, 0.22, 0.12, 0.14, 0.07, 0.05, 0.052, 0.03)
write_case(
  con,
  "repeated_means",
  means,
  disps,
  c(1, 2, 3, 30, 300, 1000)
)

means <- c(1, 1, 2, 2, 4, 8, 8, 16, 64, 256)
disps <- c(0.6, 0.55, 0.42, 0.39, 0.27, 0.18, 0.19, 0.12, 0.06, 0.035)
write_case(
  con,
  "clustered_repeated_means",
  means,
  disps,
  c(0.5, 1, 2, 5, 8, 32, 256, 512)
)

means <- c(2, 5, 20, 100)
disps <- c(0.4, 0.25, 0.1, 0.06)
write_case(
  con,
  "few_points",
  means,
  disps,
  c(2, 10, 50, 100)
)

means <- c(1, 2, 5, 10, 20, 50, 100, 200, 500, 1000)
disps <- c(2e-9, 5e-9, 8e-9, 1.1e-1, 8e-2, 5e-2, 3e-2, 2e-2, 1.5e-2, 1e-2)
write_case(
  con,
  "filtered_min_disp",
  means,
  disps,
  c(10, 50, 200, 1000)
)

means <- c(1, 2, 5, 10, 50, 100)
disps <- c(1e-9, 2e-9, 3e-9, 4e-9, 5e-9, 8e-9)
write_case(
  con,
  "constant_min_disp",
  means,
  disps,
  c(1, 10, 100, 1000)
)

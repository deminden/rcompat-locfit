suppressPackageStartupMessages(library(locfit))

root_path <- file.path("data", "locfit_hard_real_2026-06-01")
out_path <- file.path("fixtures", "hard_real_2026_subset_cases.csv")
min_disp_default <- 1e-8

selected_contrasts <- c(
  "blood_full_blocked_permutation_rep01",
  "testis_full_blocked_permutation_rep01",
  "thyroid_full_blocked_permutation_rep01"
)
target_fit_count <- 192L
target_hard_predictions <- 12L
target_range_predictions <- 12L

finite_positive <- function(x) is.finite(x) & x > 0

read_real_rows <- function(path) {
  read.delim(
    path,
    stringsAsFactors = FALSE,
    check.names = FALSE,
    na.strings = c("NA", "")
  )
}

csv_field <- function(value) {
  value <- as.character(value)
  value[is.na(value)] <- ""
  needs_quote <- grepl("[,\"\n\r]", value)
  value <- gsub("\"", "\"\"", value, fixed = TRUE)
  ifelse(needs_quote, paste0("\"", value, "\""), value)
}

write_row <- function(con, fields) {
  writeLines(paste(csv_field(fields), collapse = ","), con)
}

select_evenly <- function(rows, count) {
  if (nrow(rows) <= count) {
    return(rows)
  }
  indices <- unique(as.integer(round(seq(1L, nrow(rows), length.out = count))))
  rows[indices, , drop = FALSE]
}

dedupe_by_log_mean <- function(rows) {
  rows[!duplicated(sprintf("%.17g", rows$logBaseMean)), , drop = FALSE]
}

write_header <- function(con, input_count, prediction_count) {
  writeLines("# generated_by=fixtures/r/generate_hard_real_2026_subset_fixture.R", con)
  writeLines("# source=ignored_2026_hard_real_bundle", con)
  writeLines("# expected=black_box_R_locfit_on_committed_subset", con)
  writeLines(paste0("# r_version=", paste(R.version$major, R.version$minor, sep = ".")), con)
  writeLines(paste0("# locfit_version=", as.character(packageVersion("locfit"))), con)
  writeLines(paste0("# input_count=", input_count), con)
  writeLines(paste0("# prediction_count=", prediction_count), con)
  writeLines("# columns:case,kind,index,label,mean,disp,log_mean,log_disp,expected,expected_log", con)
}

if (!dir.exists(root_path)) {
  stop("missing ignored 2026 hard-real bundle")
}

case_rows <- list()
total_inputs <- 0L
total_predictions <- 0L

for (contrast in selected_contrasts) {
  fit_points_path <- file.path(root_path, "contrasts", contrast, "locfit", "fit_points.tsv")
  hard_rows_path <- file.path(root_path, "contrasts", contrast, "locfit", "hard_rows.tsv")
  if (!file.exists(fit_points_path)) {
    stop("missing fit-points file for contrast: ", contrast)
  }
  if (!file.exists(hard_rows_path)) {
    stop("missing hard-rows file for contrast: ", contrast)
  }

  fit_points <- read_real_rows(fit_points_path)
  fit_points <- fit_points[
    finite_positive(fit_points$baseMean) &
      finite_positive(fit_points$dispGeneEst) &
      fit_points$dispGeneEst >= min_disp_default * 10 &
      is.finite(fit_points$logBaseMean) &
      is.finite(fit_points$logDispGeneEst),
    ,
    drop = FALSE
  ]
  fit_points <- fit_points[order(fit_points$baseMean, fit_points$gene), , drop = FALSE]
  selected_fit <- select_evenly(fit_points, target_fit_count)

  hard_rows <- read_real_rows(hard_rows_path)
  hard_rows <- hard_rows[
    finite_positive(hard_rows$baseMean) &
      is.finite(hard_rows$logBaseMean),
    ,
    drop = FALSE
  ]
  hard_rows <- hard_rows[order(-hard_rows$hardScore, hard_rows$baseMean, hard_rows$gene), , drop = FALSE]
  hard_predictions <- dedupe_by_log_mean(hard_rows)
  hard_predictions <- head(hard_predictions, target_hard_predictions)
  hard_predictions$label <- sprintf("hard_%02d", seq_len(nrow(hard_predictions)) - 1L)

  range_predictions <- select_evenly(fit_points, target_range_predictions)
  range_predictions$label <- sprintf("range_%02d", seq_len(nrow(range_predictions)) - 1L)

  prediction_rows <- rbind(
    hard_predictions[, c("label", "baseMean", "logBaseMean"), drop = FALSE],
    range_predictions[, c("label", "baseMean", "logBaseMean"), drop = FALSE]
  )
  prediction_rows <- dedupe_by_log_mean(prediction_rows[order(prediction_rows$baseMean), , drop = FALSE])

  locfit_data <- data.frame(
    logDisps = selected_fit$logDispGeneEst,
    logMeans = selected_fit$logBaseMean
  )
  fit <- locfit(
    logDisps ~ lp(logMeans, nn = 0.7, deg = 2),
    data = locfit_data,
    weights = selected_fit$baseMean
  )
  expected_log <- as.numeric(predict(
    fit,
    newdata = data.frame(logMeans = prediction_rows$logBaseMean)
  ))
  prediction_rows$expectedLog <- expected_log
  prediction_rows$expected <- exp(expected_log)

  case_name <- sub("_blocked_permutation_rep01$", "", contrast)
  case_rows[[case_name]] <- list(
    input = selected_fit,
    predictions = prediction_rows
  )
  total_inputs <- total_inputs + nrow(selected_fit)
  total_predictions <- total_predictions + nrow(prediction_rows)
}

dir.create(dirname(out_path), recursive = TRUE, showWarnings = FALSE)
con <- file(out_path, open = "w")
on.exit(close(con), add = TRUE)

write_header(con, total_inputs, total_predictions)

for (case_name in names(case_rows)) {
  selected_fit <- case_rows[[case_name]]$input
  prediction_rows <- case_rows[[case_name]]$predictions

  for (i in seq_len(nrow(selected_fit))) {
    write_row(con, c(
      case_name,
      "input",
      i - 1L,
      sprintf("input_%03d", i - 1L),
      sprintf("%.17g", selected_fit$baseMean[[i]]),
      sprintf("%.17g", selected_fit$dispGeneEst[[i]]),
      sprintf("%.17g", selected_fit$logBaseMean[[i]]),
      sprintf("%.17g", selected_fit$logDispGeneEst[[i]]),
      "",
      ""
    ))
  }

  for (i in seq_len(nrow(prediction_rows))) {
    write_row(con, c(
      case_name,
      "predict",
      i - 1L,
      prediction_rows$label[[i]],
      sprintf("%.17g", prediction_rows$baseMean[[i]]),
      "",
      sprintf("%.17g", prediction_rows$logBaseMean[[i]]),
      "",
      sprintf("%.17g", prediction_rows$expected[[i]]),
      sprintf("%.17g", prediction_rows$expectedLog[[i]])
    ))
  }
}

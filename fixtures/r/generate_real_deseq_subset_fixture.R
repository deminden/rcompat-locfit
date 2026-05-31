suppressPackageStartupMessages(library(locfit))

fit_points_path <- file.path("data", "local_dispersion_fit_points.tsv")
hard_rows_path <- file.path("data", "local_dispersion_ranked_hard_rows.tsv")
all_rows_path <- file.path("data", "local_dispersion_all_rows.tsv.gz")
out_path <- file.path("fixtures", "real_deseq2_subset_cases.csv")

if (!file.exists(fit_points_path)) {
  stop("missing ignored data file: ", fit_points_path)
}
if (!file.exists(hard_rows_path)) {
  stop("missing ignored data file: ", hard_rows_path)
}
if (!file.exists(all_rows_path)) {
  stop("missing ignored data file: ", all_rows_path)
}

read_real_rows <- function(path) {
  read.delim(
    path,
    stringsAsFactors = FALSE,
    check.names = FALSE,
    na.strings = c("NA", "")
  )
}

finite_positive <- function(x) is.finite(x) & x > 0

write_header <- function(con, input_count, prediction_count) {
  writeLines("# generated_by=fixtures/r/generate_real_deseq_subset_fixture.R", con)
  writeLines(paste0("# selected_from=", fit_points_path, ";", hard_rows_path, ";", all_rows_path), con)
  writeLines("# source=real_deseq2_local_dispersion_debug_tables", con)
  writeLines("# expected=black_box_R_locfit_on_committed_subset", con)
  writeLines("# source_expected=original_full_data_dispFit_for_context_only", con)
  writeLines(paste0("# r_version=", paste(R.version$major, R.version$minor, sep = ".")), con)
  writeLines(paste0("# locfit_version=", as.character(packageVersion("locfit"))), con)
  writeLines(paste0("# input_count=", input_count), con)
  writeLines(paste0("# prediction_count=", prediction_count), con)
  writeLines("# columns:case,kind,index,gene,mean,disp,expected,source_expected", con)
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

dedupe_by_mean <- function(rows) {
  rows[!duplicated(sprintf("%.17g", rows$baseMean)), , drop = FALSE]
}

fit_points <- read_real_rows(fit_points_path)
fit_points <- fit_points[
  finite_positive(fit_points$baseMean) &
    finite_positive(fit_points$dispGeneEst) &
    is.finite(fit_points$log_baseMean) &
    is.finite(fit_points$log_dispGeneEst),
  ,
  drop = FALSE
]
fit_points <- fit_points[order(fit_points$baseMean, fit_points$gene), , drop = FALSE]

# The committed subset is quantile-spaced over the real fit points. This keeps
# the fixture compact while preserving low, middle, high, and extreme means.
target_fit_count <- 384L
fit_indices <- unique(as.integer(round(seq(1L, nrow(fit_points), length.out = target_fit_count))))
selected_fit <- fit_points[fit_indices, , drop = FALSE]

hard_rows <- read_real_rows(hard_rows_path)
hard_rows <- hard_rows[
  finite_positive(hard_rows$baseMean) &
    finite_positive(hard_rows$dispFit),
  ,
  drop = FALSE
]
hard_rows <- dedupe_by_mean(hard_rows[order(hard_rows$baseMean, hard_rows$gene), , drop = FALSE])

all_rows <- read_real_rows(gzfile(all_rows_path))
all_rows <- all_rows[
  finite_positive(all_rows$baseMean) &
    finite_positive(all_rows$dispFit),
  ,
  drop = FALSE
]
all_rows <- dedupe_by_mean(all_rows[order(all_rows$baseMean, all_rows$gene), , drop = FALSE])
all_indices <- unique(as.integer(round(seq(1L, nrow(all_rows), length.out = 24L))))
range_predictions <- all_rows[all_indices, , drop = FALSE]

prediction_rows <- rbind(
  hard_rows[, names(range_predictions), drop = FALSE],
  range_predictions
)
prediction_rows <- dedupe_by_mean(prediction_rows[order(prediction_rows$baseMean, prediction_rows$gene), , drop = FALSE])

locfit_data <- data.frame(
  logDisps = log(selected_fit$dispGeneEst),
  logMeans = log(selected_fit$baseMean)
)
fit <- locfit(
  logDisps ~ lp(logMeans, nn = 0.7, deg = 2),
  data = locfit_data,
  weights = selected_fit$baseMean
)
expected <- as.numeric(exp(predict(
  fit,
  newdata = data.frame(logMeans = log(prediction_rows$baseMean))
)))

dir.create(dirname(out_path), recursive = TRUE, showWarnings = FALSE)
con <- file(out_path, open = "w")
on.exit(close(con), add = TRUE)

write_header(con, nrow(selected_fit), nrow(prediction_rows))

case_name <- "real_deseq2_curated"
for (i in seq_len(nrow(selected_fit))) {
  write_row(con, c(
    case_name,
    "input",
    i - 1L,
    selected_fit$gene[[i]],
    sprintf("%.17g", selected_fit$baseMean[[i]]),
    sprintf("%.17g", selected_fit$dispGeneEst[[i]]),
    "",
    ""
  ))
}

for (i in seq_len(nrow(prediction_rows))) {
  write_row(con, c(
    case_name,
    "predict",
    i - 1L,
    prediction_rows$gene[[i]],
    sprintf("%.17g", prediction_rows$baseMean[[i]]),
    "",
    sprintf("%.17g", expected[[i]]),
    sprintf("%.17g", prediction_rows$dispFit[[i]])
  ))
}

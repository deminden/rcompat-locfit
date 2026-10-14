# Black-box locfit oracle: decode Rust-parsed input bits instead of parsing
# decimal doubles in R. The normal API and log API get separate expectations.
suppressPackageStartupMessages(library(locfit))
args <- commandArgs(trailingOnly = TRUE)
input_path <- if (length(args) >= 1L) args[[1L]] else "target/hard_real_oracle_inputs.csv"
output_path <- if (length(args) >= 2L) args[[2L]] else "fixtures/hard_real_2026_identical_input_oracles.csv"
input <- read.csv(input_path, comment.char = "#", colClasses = "character", check.names = FALSE)

decode_double <- function(hex) {
  if (hex == "") return(NA_real_)
  stopifnot(nchar(hex) == 16L)
  bytes <- as.raw(strtoi(substring(hex, seq(1L, 15L, 2L), seq(2L, 16L, 2L)), base = 16L))
  value <- readBin(bytes, what = "double", n = 1L, size = 8L, endian = "little")
  stopifnot(identical(bytes, writeBin(value, raw(), size = 8L, endian = "little")))
  value
}
for (name in c("mean", "disp", "log_mean", "log_disp")) {
  input[[name]] <- vapply(input[[paste0(name, "_bits")]], decode_double, numeric(1L))
}

out <- file(output_path, open = "w")
writeLines(c(
  "# generated_by=fixtures/r/generate_identical_input_oracles.R",
  "# input_encoding=IEEE754_binary64_little_endian_hex",
  paste0("# input_md5=", unname(tools::md5sum(input_path))),
  paste0("# r_version=", getRversion()),
  paste0("# locfit_version=", packageVersion("locfit")),
  paste0("# platform=", R.version$platform),
  "case,label,expected_normal,expected_normal_log,expected_stored_log"
), out)
for (name in unique(input$case)) {
  d <- input[input$case == name & input$kind == "input", ]
  p <- input[input$case == name & input$kind == "predict", ]
  keep_normal <- d$disp >= 1e-7
  keep_log <- d$log_disp >= log(1e-7)
  normal_data <- data.frame(x = log(d$mean[keep_normal]), y = log(d$disp[keep_normal]))
  log_data <- data.frame(x = d$log_mean[keep_log], y = d$log_disp[keep_log])
  normal_fit <- locfit(y ~ lp(x, nn = .7, deg = 2), data = normal_data, weights = d$mean[keep_normal])
  log_fit <- locfit(y ~ lp(x, nn = .7, deg = 2), data = log_data, weights = d$mean[keep_log])
  expected_normal_log <- as.numeric(predict(normal_fit, newdata = data.frame(x = log(p$mean))))
  expected_stored_log <- as.numeric(predict(log_fit, newdata = data.frame(x = p$log_mean)))
  writeLines(sprintf("%s,%s,%.17g,%.17g,%.17g", name, p$label, exp(expected_normal_log), expected_normal_log, expected_stored_log), out)
}
close(out)
cat("Wrote", output_path, "\n")

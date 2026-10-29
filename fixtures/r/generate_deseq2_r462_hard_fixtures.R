# Generate numeric oracles through actual DESeq2 calls. No implementation
# source from DESeq2 or locfit is inspected or copied.
suppressPackageStartupMessages(library(DESeq2))
suppressPackageStartupMessages(library(airway))
stopifnot(getRversion() == "4.6.2")

args <- commandArgs(trailingOnly = TRUE)
output_path <- if (length(args)) args[[1L]] else "fixtures/deseq2_r462_hard_cases.csv"
all_queries <- length(args) >= 2L && args[[2L]] == "--all-queries"
query_path <- "fixtures/deseq2_r462_hard_queries.csv"
queries <- read.csv(query_path, comment.char = "#", colClasses = "character")
data(airway)

hex <- function(x) {
  paste(sprintf("%02x", as.integer(writeBin(as.double(x), raw(), size = 8L,
                                          endian = "little"))), collapse = "")
}
out <- file(output_path, "w")
writeLines(c(
  "# generated_by=fixtures/r/generate_deseq2_r462_hard_fixtures.R",
  paste0("# r_version=", getRversion()),
  paste0("# DESeq2_version=", packageVersion("DESeq2")),
  paste0("# locfit_version=", packageVersion("locfit")),
  paste0("# locfit_build=", packageDescription("locfit")$Built),
  paste0("# airway_version=", packageVersion("airway")),
  paste0("# platform=", R.version$platform),
  paste0("# query_md5=", unname(tools::md5sum(query_path))),
  "# input_encoding=IEEE754_binary64_little_endian_hex",
  "# simulated_seed=462; simulated_genes=600; simulated_samples=12; betaSD=1",
  "# airway_subset=512 equally spaced ranks of raw mean counts; rowSums>=10",
  "case,kind,gene,mean_bits,disp_bits,log_mean_bits,log_disp_bits,expected_log_bits"
), out)

run_case <- function(name, dds) {
  dds <- DESeq(dds, fitType = "local", quiet = TRUE)
  stopifnot(identical(attr(dispersionFunction(dds), "fitType"), "local"))
  m <- as.data.frame(mcols(dds))
  usable <- is.finite(m$baseMean) & m$baseMean > 0 &
    is.finite(m$dispGeneEst) & m$dispGeneEst > 0 &
    is.finite(m$dispFit) & m$dispFit > 0
  m <- m[usable, , drop = FALSE]
  # Keep every observation used by DESeq2's local trend. Predictions also
  # include low-dispersion genes excluded from fitting, exercising extrapolation.
  fit_rows <- which(m$dispGeneEst >= 1e-7)
  for (i in fit_rows) {
    writeLines(paste(c(name, "input", rownames(m)[i], hex(m$baseMean[i]),
      hex(m$dispGeneEst[i]), hex(log(m$baseMean[i])),
      hex(log(m$dispGeneEst[i])), ""), collapse = ","), out)
  }
  query_genes <- if (all_queries) rownames(m) else queries$gene[queries$case == name]
  indices <- match(query_genes, rownames(m))
  stopifnot(length(indices) > 0L, !anyNA(indices))
  for (i in indices) {
    writeLines(paste(c(name, "predict", rownames(m)[i], hex(m$baseMean[i]),
      "", hex(log(m$baseMean[i])), "", hex(log(m$dispFit[i]))),
      collapse = ","), out)
  }
  cat(name, ":", length(fit_rows), "fit rows,", length(indices), "queries\n")
}

counts <- assay(airway)
keep <- which(rowSums(counts) >= 10)
ordered <- keep[order(rowMeans(counts[keep, , drop = FALSE]), rownames(counts)[keep])]
selected <- ordered[unique(round(seq(1, length(ordered), length.out = 512)))]
set.seed(462)
run_case("airway_512", DESeqDataSet(airway[selected, ], design = ~ cell + dex))
set.seed(462)
run_case("simulated_600", makeExampleDESeqDataSet(n = 600L, m = 12L, betaSD = 1))
run_case("airway_full", DESeqDataSet(airway[keep, ], design = ~ cell + dex))
close(out)

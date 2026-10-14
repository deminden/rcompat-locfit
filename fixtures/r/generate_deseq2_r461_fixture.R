suppressPackageStartupMessages(library(DESeq2))
set.seed(461)
dds <- makeExampleDESeqDataSet(n=1200L,m=12L,betaSD=1)
dds <- DESeq(dds,fitType="local",quiet=TRUE)
meta <- as.data.frame(mcols(dds))
keep <- is.finite(meta$baseMean) & meta$baseMean > 0 & is.finite(meta$dispGeneEst) & meta$dispGeneEst > 0 & is.finite(meta$dispFit) & meta$dispFit > 0
meta<-meta[keep,,drop=FALSE]
hex <- function(x) paste(sprintf("%02x",as.integer(writeBin(as.double(x),raw(),size=8L,endian="little"))),collapse="")
f <- file("fixtures/deseq2_r461_local_trend.csv","w")
writeLines(c(paste0("# r=",getRversion()),paste0("# DESeq2=",packageVersion("DESeq2")),paste0("# locfit=",packageVersion("locfit")),"gene,mean_bits,disp_bits,log_mean_bits,log_disp_bits,expected_log_bits"),f)
for (i in seq_len(nrow(meta))) writeLines(paste(c(rownames(meta)[[i]],hex(meta$baseMean[[i]]),hex(meta$dispGeneEst[[i]]),hex(log(meta$baseMean[[i]])),hex(log(meta$dispGeneEst[[i]])),hex(log(meta$dispFit[[i]]))),collapse=","),f)
close(f)
cat("Real DESeq2 local fit: ",nrow(meta),"genes; fitType",attr(dispersionFunction(dds),"fitType"),"\n")


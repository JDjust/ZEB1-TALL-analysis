# Install ALLCatchR2 in user lib and classify GSE280250. No new science.
userlib <- "/data-b/liangfuhua/R/library"
dir.create(userlib, recursive = TRUE, showWarnings = FALSE)
.libPaths(c(userlib, .libPaths()))
need <- c("rlang","caret","glmnet","LiblineaR","kknn","randomForest",
          "ranger","caTools","elasticnet","ggplot2","ggrepel","umap")
miss <- need[!vapply(need, requireNamespace, quietly = TRUE, FUN.VALUE = TRUE)]
if (length(miss)) {
  install.packages(miss, repos = "https://cloud.r-project.org")
}
if (!requireNamespace("singscore", quietly = TRUE)) {
  if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager", repos = "https://cloud.r-project.org")
  }
  BiocManager::install("singscore", update = FALSE, ask = FALSE)
}
if (!requireNamespace("ALLCatchR2", quietly = TRUE)) {
  install.packages("/data-b/liangfuhua/zeb1_a8/ALLCatchR2", repos = NULL, type = "source")
}
library(ALLCatchR2)
cnt <- read.delim("/data-b/liangfuhua/zeb1_a8/a8_counts.tsv", check.names = FALSE)
dir.create("/data-b/liangfuhua/zeb1_a8/allcatchr2_plots", showWarnings = FALSE)
out <- allcatchr_1.1(
  Lineage = "T-ALL",
  Counts.file = cnt,
  ID_class = "symbol",
  sep = "\t",
  out.file = "/data-b/liangfuhua/zeb1_a8/allcatchr2_predictions.tsv",
  plot.path = "/data-b/liangfuhua/zeb1_a8/allcatchr2_plots"
)
cat("done class", class(out), "len", length(out), "\n")
if (is.list(out)) cat("names", paste(names(out), collapse = ","), "\n")

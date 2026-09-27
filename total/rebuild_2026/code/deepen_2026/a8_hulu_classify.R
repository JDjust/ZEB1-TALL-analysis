userlib <- "/data-b/liangfuhua/R/library"
.libPaths(c(userlib, .libPaths()))
library(ALLCatchR2)
dir.create("/data-b/liangfuhua/zeb1_a8/allcatchr2_plots", showWarnings = FALSE)
out <- allcatchr_1.1(
  Lineage = "T-ALL",
  Counts.file = "/data-b/liangfuhua/zeb1_a8/a8_counts.tsv",
  ID_class = "symbol",
  sep = "\t",
  out.file = "/data-b/liangfuhua/zeb1_a8/allcatchr2_predictions.tsv",
  plot.path = "/data-b/liangfuhua/zeb1_a8/allcatchr2_plots"
)
cat("done\n")
if (is.list(out) || is.data.frame(out)) {
  cat("class", paste(class(out), collapse = ","), "\n")
  if (is.data.frame(out)) cat("nrow", nrow(out), "ncol", ncol(out), "\n")
}

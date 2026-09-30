userlib <- "/data-b/liangfuhua/R/library"
dir.create(userlib, recursive = TRUE, showWarnings = FALSE)
.libPaths(c(userlib, .libPaths()))
Sys.setenv(HOME = "/data-b/liangfuhua")
dir.create("/data-b/liangfuhua/.R", showWarnings = FALSE)
writeLines(c(
  "CC = /usr/bin/gcc",
  "CXX = /usr/bin/g++",
  "FC = /data-a/bioinfo/public-envs/r-bio-2026.08/bin/x86_64-conda-linux-gnu-gfortran",
  "F77 = /data-a/bioinfo/public-envs/r-bio-2026.08/bin/x86_64-conda-linux-gnu-gfortran"
), "/data-b/liangfuhua/.R/Makevars")
install.packages(c("randomForest", "lars", "elasticnet"),
                 repos = "https://cloud.r-project.org")
if (!requireNamespace("singscore", quietly = TRUE)) {
  # Bioconductor index is blocked; GitHub source is R-only.
  if (!requireNamespace("remotes", quietly = TRUE)) {
    install.packages("remotes", repos = "https://cloud.r-project.org")
  }
  remotes::install_github("DavisLaboratory/singscore", upgrade = "never", dependencies = FALSE)
}
if (!requireNamespace("ALLCatchR2", quietly = TRUE)) {
  install.packages("/data-b/liangfuhua/zeb1_a8/ALLCatchR2", repos = NULL, type = "source")
}
for (p in c("randomForest", "elasticnet", "singscore", "ALLCatchR2")) {
  cat(p, requireNamespace(p, quietly = TRUE), "\n")
}

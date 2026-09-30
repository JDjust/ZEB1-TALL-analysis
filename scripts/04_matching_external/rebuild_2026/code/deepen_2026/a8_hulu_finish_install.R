userlib <- "/data-b/liangfuhua/R/library"
dir.create(userlib, recursive = TRUE, showWarnings = FALSE)
.libPaths(c(userlib, .libPaths()))
options(timeout = 600)
Sys.setenv(MAKEFLAGS = "-j4")

need <- c("caTools", "ggrepel", "ggplot2")
have <- rownames(installed.packages())
miss <- setdiff(need, have)
if (length(miss)) {
  cat("CRAN missing:", paste(miss, collapse = ", "), "\n")
  install.packages(miss, lib = userlib, repos = "https://cloud.r-project.org")
} else {
  cat("CRAN extras already present\n")
}

if (!requireNamespace("singscore", quietly = TRUE)) {
  cat("install singscore stub\n")
  install.packages("/data-b/liangfuhua/zeb1_a8/singscore_stub",
                   repos = NULL, type = "source", lib = userlib)
} else {
  cat("singscore already installed\n")
}

if (!requireNamespace("ALLCatchR2", quietly = TRUE)) {
  cat("install ALLCatchR2 from local clone\n")
  install.packages("/data-b/liangfuhua/zeb1_a8/ALLCatchR2",
                   repos = NULL, type = "source", lib = userlib)
} else {
  cat("ALLCatchR2 already installed\n")
}

stopifnot(requireNamespace("singscore", quietly = TRUE))
stopifnot(requireNamespace("ALLCatchR2", quietly = TRUE))
cat("packages ready\n")

source("/data-b/liangfuhua/zeb1_a8/a8_hulu_classify.R")

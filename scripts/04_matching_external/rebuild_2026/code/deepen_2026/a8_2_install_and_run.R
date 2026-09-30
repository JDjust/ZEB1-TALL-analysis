# A8-2: ALLCatchR2 annotation of GSE280250. Subtype context only.
options(timeout = 600)
if (!requireNamespace("ALLCatchR2", quietly = TRUE)) {
  if (!requireNamespace("devtools", quietly = TRUE)) {
    install.packages("devtools", repos = "https://cloud.r-project.org")
  }
  devtools::install_github("ThomasBeder/ALLCatchR2", upgrade = "never", quiet = FALSE)
}
suppressPackageStartupMessages(library(ALLCatchR2))
cat("ALLCatchR2 loaded\n"); flush.console()
print(ls("package:ALLCatchR2"))

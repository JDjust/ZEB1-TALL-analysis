# rebuild_all_figures_nature.R
# Master script: rebuild all main + supplementary figures with theme_nature
# Run from: D:\_bioinformation\ZEB1\total\rebuild_2026\code\
#
# Usage:
#   Rscript rebuild_all_figures_nature.R          # all figures
#   Rscript rebuild_all_figures_nature.R F1 F3    # specific figures only

suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
  library(ragg)
})

args <- commandArgs(trailingOnly = TRUE)
script_dir <- normalizePath(dirname(sub("^--file=", "", grep("^--file=",
  commandArgs(FALSE), value = TRUE))), winslash = "/")
source(file.path(script_dir, "theme_nature.R"))

root <- normalizePath(file.path(script_dir, "../../.."), winslash = "/")
cat("Project root:", root, "\n")
cat("Theme: theme_nature (7pt base, 600dpi, Nature/Blood format)\n\n")

# Build list
all_figures <- c("F1","F2","F3","F4","F5","F6","F7",
                 "S1","S2","S3","S4","S5","S6","S7")
targets <- if (length(args) > 0) args else all_figures

for (fig in targets) {
  script_name <- switch(fig,
    F1 = "build_f1_nature.R",
    F2 = "build_f2_ggplot.R",
    F3 = "build_f3_ggplot.R",
    F4 = "build_f4_ggplot.R",
    F5 = "build_f5_ggplot.R",
    F6 = "build_f6_ggplot.R",
    S1 = "build_s1_ggplot.R",
    S2 = "build_s2_ggplot.R",
    S3 = "build_s3_ggplot.R",
    S4 = "build_s4_ggplot.R",
    S5 = "build_s5_ggplot.R",
    S6 = "build_s6_ggplot.R",
    S7 = "build_s7_ggplot.R",
    NULL
  )
  
  if (is.null(script_name)) {
    cat(sprintf("[SKIP] %s: no build script mapped\n", fig))
    next
  }
  
  script_path <- file.path(script_dir, script_name)
  if (!file.exists(script_path)) {
    cat(sprintf("[SKIP] %s: %s not found\n", fig, script_name))
    next
  }
  
  cat(sprintf("[BUILD] %s from %s ... ", fig, script_name))
  t0 <- Sys.time()
  tryCatch({
    source(script_path, local = new.env(parent = globalenv()))
    dt <- round(as.numeric(difftime(Sys.time(), t0, units = "secs")), 1)
    cat(sprintf("OK (%.1fs)\n", dt))
  }, error = function(e) {
    cat(sprintf("FAILED: %s\n", conditionMessage(e)))
  })
}

cat("\nDone. Check edition_20260925/figures/ for outputs.\n")

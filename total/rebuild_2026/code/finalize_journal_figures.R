# Export builder-designed figures at their native final physical dimensions.
# Builders own typography, geometry, layout and page size. This exporter never
# rescales plot objects or changes the aspect ratio.
suppressPackageStartupMessages({library(ggplot2); library(ragg)})

argv <- grep("^--file=", commandArgs(FALSE), value = TRUE)
script <- normalizePath(sub("^--file=", "", argv), winslash = "/")
code <- dirname(script)
root <- normalizePath(file.path(code, "../../.."), winslash = "/")
sub_root <- file.path(root, "submission/figures")
dir.create(file.path(sub_root, "main"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(sub_root, "supplementary"), recursive = TRUE, showWarnings = FALSE)

run_capture <- function(builder, keep) {
  e <- new.env(parent = .GlobalEnv)
  captured <- new.env(parent = emptyenv())
  e$ggsave <- function(filename, plot = last_plot(), width = NA, height = NA, ...) {
    key <- sub("\\.[A-Za-z]+$", "", basename(filename))
    if (key %in% keep && !exists(key, envir = captured, inherits = FALSE))
      assign(key, list(plot = plot, width = width, height = height), envir = captured)
    invisible(filename)
  }
  e$save_nature <- function(plot, basename, width, height, ...) {
    key <- basename(basename)
    if (key %in% keep && !exists(key, envir = captured, inherits = FALSE))
      assign(key, list(plot = plot, width = width, height = height), envir = captured)
    invisible(basename)
  }
  sys.source(file.path(code, builder), envir = e)
  as.list(captured)
}

export_native <- function(obj, outbase) {
  p <- obj$plot
  w <- obj$width; h <- obj$height
  stopifnot(is.numeric(w), is.numeric(h), length(w) == 1, length(h) == 1,
            is.finite(w), is.finite(h), w > 0, h > 0)
  ggsave(paste0(outbase, ".pdf"), p, width = w, height = h, units = "in",
         device = cairo_pdf, bg = "white", limitsize = FALSE)
  ggsave(paste0(outbase, ".tiff"), p, width = w, height = h, units = "in",
         device = agg_tiff, dpi = 600, compression = "lzw", bg = "white",
         limitsize = FALSE)
  ggsave(paste0(outbase, ".png"), p, width = w, height = h, units = "in",
         device = agg_png, dpi = 300, bg = "white", limitsize = FALSE)
  cat(sprintf("  %-16s %.0f x %.0f mm\n", basename(outbase),
              w * 25.4, h * 25.4))
}

jobs <- list(
  list(b = "build_f1_nature.R", keys = "F1", out = "Figure1", dir = "main"),
  list(b = "build_f2_ggplot.R", keys = "F2", out = "Figure2", dir = "main"),
  list(b = "build_f3_revision2.R", keys = "F3", out = "Figure3", dir = "main"),
  list(b = "build_f4_ggplot.R", keys = "F4", out = "Figure4", dir = "main"),
  list(b = "build_main_validation_relayout.R", keys = "F5",
       out = "Figure5", dir = "main"),
  list(b = "build_f6_gateA.R", keys = "F6",
       out = "Figure6", dir = "main"),
  list(b = "build_f7_revision3.R", keys = "F7",
       out = "Figure7", dir = "main"),
  list(b = "build_main_validation_relayout.R", keys = "F7",
       out = "FigureS9", dir = "supplementary"),
  list(b = "build_s3_audit.R", keys = c("S3_p1", "S3_p2"),
       out = c("FigureS3_page1", "FigureS3_page2"), dir = "supplementary"),
  list(b = "build_s8_clinical_audit.R", keys = "S8",
       out = "FigureS8", dir = "supplementary"),
  list(b = "build_supplement_consolidated.R",
       keys = c(paste0(rep(paste0("S", 1:2), each = 2), "_p", 1:2),
                paste0("S", 4:7)),
       out = c(paste0(rep(paste0("FigureS", 1:2), each = 2),
                      "_page", 1:2), paste0("FigureS", 4:7)),
       dir = "supplementary")
)

want <- commandArgs(trailingOnly = TRUE)
for (j in jobs) {
  sel <- if (length(want)) intersect(j$keys, want) else j$keys
  if (!length(sel)) next
  cat("[", j$b, "]\n", sep = "")
  got <- run_capture(j$b, sel)
  for (k in sel) {
    if (is.null(got[[k]])) stop("Missing plot: ", k, " in ", j$b)
    idx <- match(k, j$keys)
    export_native(got[[k]], file.path(sub_root, j$dir, j$out[idx]))
  }
}
writeLines(capture.output(sessionInfo()), file.path(sub_root, "R_sessionInfo.txt"))
cat("Done ->", sub_root, "\n")

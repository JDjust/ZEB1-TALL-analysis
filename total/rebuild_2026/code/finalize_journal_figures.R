# Journal-width finalization of all main and supplementary figures.
# Each existing figure builder is executed unchanged; its final ggplot/patchwork
# object is captured, fonts and mark sizes are rescaled for a 180-mm canvas,
# and the result is exported to ZEB1/submission/figures.
#
# Usage: Rscript finalize_journal_figures.R            (all figures)
#        Rscript finalize_journal_figures.R F2 S3      (selected figures)
suppressPackageStartupMessages({
  library(ggplot2); library(patchwork); library(ragg)
})

argv <- grep("^--file=", commandArgs(FALSE), value = TRUE)
script <- normalizePath(sub("^--file=", "", argv), winslash = "/")
code <- dirname(script)
root <- normalizePath(file.path(code, "../../.."), winslash = "/")
sub_root <- file.path(root, "submission/figures")
dir.create(file.path(sub_root, "main"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(sub_root, "supplementary"), recursive = TRUE, showWarnings = FALSE)

WIDTH_IN <- 180 / 25.4
MAX_H_IN <- 240 / 25.4

getp <- function(el, nm) {
  v <- tryCatch(el[[nm]], error = function(e) NULL)
  if (is.null(v)) v <- tryCatch(S7::prop(el, nm), error = function(e) NULL)
  v
}
is_abs <- function(x) is.numeric(x) && length(x) == 1 && !inherits(x, "rel") && is.finite(x)

rescale_theme <- function(p, ft, fg) {
  th <- p$theme
  if (is.null(th) || !length(th)) return(p)
  adds <- list()
  for (nm in names(th)) {
    el <- th[[nm]]
    if (inherits(el, "element_text") || inherits(el, "ggplot2::element_text")) {
      s <- getp(el, "size"); if (is_abs(s)) adds[[nm]] <- element_text(size = s * ft)
    } else if (inherits(el, "element_line") || inherits(el, "ggplot2::element_line")) {
      s <- getp(el, "linewidth"); if (is_abs(s)) adds[[nm]] <- element_line(linewidth = s * fg)
    } else if (inherits(el, "element_rect") || inherits(el, "ggplot2::element_rect")) {
      s <- getp(el, "linewidth"); if (is_abs(s)) adds[[nm]] <- element_rect(linewidth = s * fg)
    } else if (inherits(el, "element_geom") || inherits(el, "ggplot2::element_geom")) {
      args <- list()
      for (k in c("linewidth", "borderwidth", "pointsize")) {
        s <- getp(el, k); if (is_abs(s)) args[[k]] <- s * fg
      }
      s <- getp(el, "fontsize"); if (is_abs(s)) args$fontsize <- s * ft
      if (length(args)) adds[[nm]] <- do.call(element_geom, args)
    }
  }
  if (length(adds)) p <- p + do.call(theme, adds)
  p
}

rescale_layers <- function(p, ft, fg) {
  for (i in seq_along(p$layers)) {
    ly <- p$layers[[i]]
    is_text <- inherits(ly$geom, c("GeomText", "GeomLabel", "GeomTextRepel", "GeomLabelRepel"))
    for (k in c("size", "linewidth", "stroke")) {
      v <- ly$aes_params[[k]]
      if (is.numeric(v)) {
        f <- if (k == "size" && is_text) ft else fg
        ly$aes_params[[k]] <- v * f
      }
    }
    p$layers[[i]] <- ly
  }
  p
}

rescale_plot <- function(p, ft, fg) {
  if (inherits(p, "patchwork")) {
    pl <- p$patches$plots
    for (i in seq_along(pl)) pl[[i]] <- rescale_plot(pl[[i]], ft, fg)
    p$patches$plots <- pl
  }
  if (inherits(p, "ggplot")) {
    p <- rescale_layers(p, ft, fg)
    p <- rescale_theme(p, ft, fg)
  }
  p
}

# ── Capture: run a builder with ggsave intercepted ─────────────────────────
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
      assign(key, list(plot = plot, width = width, height = height, final = TRUE),
             envir = captured)
    invisible(basename)
  }
  sys.source(file.path(code, builder), envir = e)
  as.list(captured)
}

export <- function(obj, outbase, ft, fg, hscale) {
  p <- obj$plot
  if (isTRUE(obj$final)) {
    w <- obj$width; h <- obj$height
  } else {
    p <- rescale_plot(p, ft, fg)
    w <- WIDTH_IN
    h <- min(MAX_H_IN, obj$height * (WIDTH_IN / obj$width) * hscale)
  }
  ggsave(paste0(outbase, ".pdf"), p, width = w, height = h, units = "in",
         device = cairo_pdf, bg = "white", limitsize = FALSE)
  ggsave(paste0(outbase, ".tiff"), p, width = w, height = h, units = "in",
         device = agg_tiff, dpi = 600, compression = "lzw", bg = "white",
         limitsize = FALSE)
  ggsave(paste0(outbase, ".png"), p, width = w, height = h, units = "in",
         device = agg_png, dpi = 300, bg = "white", limitsize = FALSE)
  cat(sprintf("  %-14s %.2f x %.2f in (%.0f x %.0f mm)\n", basename(outbase),
              w, h, w * 25.4, h * 25.4))
}

# builder, captured keys, output names, text factor, geom factor, height scale
jobs <- list(
  list(b = "build_f1_nature.R", keys = "F1", out = "Figure1", dir = "main",
       ft = 1, fg = 1, hs = 1),
  list(b = "build_f2_ggplot.R", keys = "F2", out = "Figure2", dir = "main",
       ft = .64, fg = .55, hs = 1.30),
  list(b = "build_f3_revision2.R", keys = "F3", out = "Figure3", dir = "main",
       ft = .64, fg = .55, hs = 1.25),
  list(b = "build_f4_ggplot.R", keys = "F4", out = "Figure4", dir = "main",
       ft = .64, fg = .55, hs = 1.25),
  list(b = "build_main_validation_relayout.R", keys = c("F5", "F6", "F7"),
       out = c("Figure5", "Figure6", "Figure7"), dir = "main",
       ft = .60, fg = .50, hs = c(1.60, 1.25, 1.25)),
  list(b = "build_supplement_consolidated.R", keys = paste0("S", 1:7),
       out = paste0("FigureS", 1:7), dir = "supplementary",
       ft = .56, fg = .45, hs = c(1.40, 1.20, 1.20, 1.25, 1.25, 1.25, 1.70))
)

want <- commandArgs(trailingOnly = TRUE)
for (j in jobs) {
  sel <- if (length(want)) intersect(j$keys, want) else j$keys
  if (!length(sel)) next
  cat("[", j$b, "]\n", sep = "")
  got <- run_capture(j$b, sel)
  for (k in sel) {
    if (is.null(got[[k]])) { cat("  MISSING ", k, "\n"); next }
    idx <- match(k, j$keys)
    hs <- if (length(j$hs) > 1) j$hs[idx] else j$hs
    export(got[[k]], file.path(sub_root, j$dir, j$out[idx]), j$ft, j$fg, hs)
  }
}
writeLines(capture.output(sessionInfo()), file.path(sub_root, "R_sessionInfo.txt"))
cat("Done ->", sub_root, "\n")

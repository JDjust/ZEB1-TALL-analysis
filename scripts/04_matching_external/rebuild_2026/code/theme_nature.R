# theme_nature.R - Publication-grade ggplot2 theme for top-tier journals
# Designed for Blood / Leukemia / Nature Communications figure standards
# Consistent typography, minimalist axes, high information density
#
# Usage: source("theme_nature.R") before building any figure script.

suppressPackageStartupMessages({
  library(ggplot2)
  library(grid)
})

# ── Master palette ──────────────────────────────────────────────────────────
pal <- list(
  ink       = "#1a1a2e",       # near-black for text and primary elements
  axis      = "#4a4a5a",       # dark grey for axis text
  grid      = "#e8e8ee",       # very light grey for optional gridlines
  border    = "#c5c5d2",       # panel border (when needed)
  bg        = "white",         # background
  annot     = "#6b6b7b",       # annotation text

  # Subtype palette (17 T-ALL molecular subtypes, colorblind-optimized)
  BCL11B           = "#6a3d9a",
  SPI1             = "#8968b0",
  LMO2_gd          = "#a67fb8",
  TME              = "#3182bd",
  NUP98            = "#6baed6",
  MLLT10           = "#41ab5d",
  ETP              = "#00796b",
  HOXA9_TCR        = "#7c9e5a",
  STAG2_LMO2       = "#b8a042",
  TLX1             = "#d4722c",
  KMT2A            = "#5c7ea0",
  NKX2_5           = "#8c9da6",
  TAL1_ab          = "#c94c5a",
  NKX2_1           = "#8c6b6b",
  TAL1_DP          = "#d94452",
  NUP214           = "#7d8e6e",
  TLX3             = "#d4a017",

  # Functional colors
  zeb1      = "#c0392b",       # deep red for ZEB1
  zeb2      = "#2471a3",       # steel blue for ZEB2
  lmo2      = "#b8860b",       # dark goldenrod for LMO2
  highlight = "#e74c3c",       # bright red for highlights
  positive  = "#c0392b",       # positive direction
  negative  = "#2471a3",       # negative direction
  neutral   = "#95a5a6",       # neutral/non-significant

  # ETP/non-ETP
  etp       = "#00796b",       # teal for ETP
  nonetp    = "#d4722c",       # burnt orange for non-ETP

  # Diverging scale endpoints
  div_low   = "#315f82",       # blue end
  div_mid   = "#f6f3ec",       # light neutral midpoint
  div_high  = "#b84922",       # warm orange-red end

  # Sequential scale
  seq_low   = "#f0f0f4",       # light
  seq_high  = "#1a1a2e"        # dark
)

subtype_palette <- c(
  "BCL11B"               = pal$BCL11B,
  "SPI1"                 = pal$SPI1,
  "LMO2 gamma-delta-like"= pal$LMO2_gd,
  "TME-enriched"         = pal$TME,
  "NUP98"                = pal$NUP98,
  "MLLT10"               = pal$MLLT10,
  "ETP-like"             = pal$ETP,
  "HOXA9 TCR"            = pal$HOXA9_TCR,
  "STAG2&LMO2"           = pal$STAG2_LMO2,
  "TLX1"                 = pal$TLX1,
  "KMT2A"                = pal$KMT2A,
  "NKX2-5"               = pal$NKX2_5,
  "TAL1 alpha-beta-like" = pal$TAL1_ab,
  "NKX2-1"               = pal$NKX2_1,
  "TAL1 DP-like"         = pal$TAL1_DP,
  "NUP214"               = pal$NUP214,
  "TLX3"                 = pal$TLX3
)

# Immature-to-mature ordered subtype levels
subtype_order <- c(
  "BCL11B", "LMO2 gamma-delta-like", "TME-enriched", "SPI1",
  "ETP-like", "STAG2&LMO2", "NUP98", "MLLT10",
  "KMT2A", "HOXA9 TCR", "NKX2-5", "NKX2-1",
  "TAL1 alpha-beta-like", "TLX1", "TAL1 DP-like", "NUP214", "TLX3"
)

# ── Core theme ──────────────────────────────────────────────────────────────
theme_nature <- function(base_size = 7, base_family = "Arial") {
  theme_classic(base_size = base_size, base_family = base_family) %+replace%
    theme(
      # Text hierarchy
      plot.title       = element_text(size = rel(1.25), face = "bold",
                                      colour = pal$ink, hjust = 0,
                                      margin = margin(b = 3)),
      plot.subtitle    = element_text(size = rel(0.95), colour = pal$axis,
                                      hjust = 0, margin = margin(b = 2)),
      plot.caption     = element_text(size = rel(0.82), colour = pal$annot,
                                      hjust = 1),

      # Axes
      axis.title       = element_text(size = rel(1.0), colour = pal$ink),
      axis.title.x     = element_text(margin = margin(t = 4)),
      axis.title.y     = element_text(angle = 90, margin = margin(r = 4)),
      axis.text        = element_text(size = rel(0.88), colour = pal$axis),
      axis.line        = element_line(colour = pal$axis, linewidth = 0.3),
      axis.ticks       = element_line(colour = pal$axis, linewidth = 0.25),
      axis.ticks.length = unit(1.5, "pt"),

      # Panel
      panel.background = element_rect(fill = pal$bg, colour = NA),
      panel.grid       = element_blank(),
      plot.background  = element_rect(fill = pal$bg, colour = NA),

      # Legend
      legend.title     = element_text(size = rel(0.88), face = "bold",
                                      colour = pal$ink),
      legend.text      = element_text(size = rel(0.82), colour = pal$axis),
      legend.key.size  = unit(8, "pt"),
      legend.key       = element_rect(fill = NA, colour = NA),
      legend.background = element_rect(fill = NA, colour = NA),
      legend.margin    = margin(2, 2, 2, 2),

      # Strip (facets)
      strip.background = element_blank(),
      strip.text       = element_text(size = rel(0.95), face = "bold",
                                      colour = pal$ink, hjust = 0),

      # Margins
      plot.margin      = margin(4, 4, 4, 4, "pt")
    )
}

# ── Variant themes ──────────────────────────────────────────────────────────
theme_nature_minimal <- function(base_size = 7, base_family = "Arial") {
  theme_nature(base_size, base_family) %+replace%
    theme(
      axis.line  = element_blank(),
      axis.ticks = element_blank(),
      panel.border = element_rect(colour = pal$border, fill = NA,
                                  linewidth = 0.3)
    )
}

theme_nature_heatmap <- function(base_size = 7, base_family = "Arial") {
  theme_nature(base_size, base_family) %+replace%
    theme(
      axis.line  = element_blank(),
      axis.ticks = element_blank(),
      panel.grid = element_blank()
    )
}

theme_nature_umap <- function(base_size = 7, base_family = "Arial") {
  theme_nature(base_size, base_family) %+replace%
    theme(
      axis.line  = element_blank(),
      axis.ticks = element_blank(),
      axis.text  = element_blank(),
      panel.border = element_rect(colour = pal$border, fill = NA,
                                  linewidth = 0.25)
    )
}

# ── Convenience scale functions ─────────────────────────────────────────────
scale_fill_diverging <- function(midpoint = 0, limits = NULL, ...) {
  scale_fill_gradient2(
    low = pal$div_low, mid = pal$div_mid, high = pal$div_high,
    midpoint = midpoint, limits = limits,
    oob = scales::squish, ...
  )
}

scale_colour_subtype <- function(...) {
  scale_colour_manual(values = subtype_palette, ...)
}

scale_fill_subtype <- function(...) {
  scale_fill_manual(values = subtype_palette, ...)
}

# ── Panel label: bold letter drawn inside the panel title (always visible) ──
tag_plot <- function(p, letter, subtitle = NULL) {
  p + labs(title = letter, subtitle = subtitle) +
    theme(
      plot.title = element_text(size = 9, face = "bold", colour = pal$ink,
                                hjust = 0, margin = margin(b = 1)),
      plot.subtitle = element_text(size = 6.5, colour = pal$axis,
                                   hjust = 0, margin = margin(b = 2))
    )
}

# ── Standard figure dimensions (inches, Nature style) ──────────────────────
fig_dims <- list(
  single_col  = c(w = 3.35,  h = 2.8),   # ~85 mm
  one_half    = c(w = 5.51,  h = 4.5),   # ~140 mm
  two_col     = c(w = 7.09,  h = 5.5),   # ~180 mm
  full_page   = c(w = 7.09,  h = 9.45)   # ~180 x 240 mm
)

# ── Export helper ───────────────────────────────────────────────────────────
save_nature <- function(plot, basename, width, height,
                        dpi = 600, formats = c("pdf", "png", "tiff", "svg")) {
  if ("pdf" %in% formats)
    ggsave(paste0(basename, ".pdf"), plot, width = width, height = height,
           units = "in", device = cairo_pdf, bg = "white", limitsize = FALSE)
  if ("png" %in% formats)
    ggsave(paste0(basename, ".png"), plot, width = width, height = height,
           units = "in", device = ragg::agg_png, dpi = dpi, bg = "white",
           limitsize = FALSE)
  if ("tiff" %in% formats)
    ggsave(paste0(basename, ".tiff"), plot, width = width, height = height,
           units = "in", device = ragg::agg_tiff, dpi = dpi, bg = "white",
           limitsize = FALSE,
           compression = "lzw")
  if ("svg" %in% formats)
    ggsave(paste0(basename, ".svg"), plot, width = width, height = height,
           units = "in", device = svg, bg = "white", limitsize = FALSE)
  invisible(basename)
}

cat("theme_nature.R loaded: base_size=7pt, Nature/Blood/Leukemia format\n")

# Figure 3: developmental inheritance versus subtype-specific ZEB balance.
# All estimates come from frozen tables. No model is refitted here.
suppressPackageStartupMessages({
  library(ggplot2); library(patchwork); library(ggridges); library(ggrepel); library(ragg)
})
argv <- grep("^--file=", commandArgs(FALSE), value = TRUE)
script <- normalizePath(sub("^--file=", "", argv), winslash = "/")
root <- normalizePath(file.path(dirname(script), "../../.."), winslash = "/")
vdir <- file.path(root, "total/data/validation/zeb_developmental_residual")
out <- file.path(root, "total/rebuild_2026/edition_20260925/figures/main")
sdir <- file.path(root, "total/rebuild_2026/data/source_data_rebuilt/F3")
dir.create(out, recursive = TRUE, showWarnings = FALSE)
dir.create(sdir, recursive = TRUE, showWarnings = FALSE)
rd <- function(x) as.data.frame(readr::read_tsv(file.path(vdir, x),
  locale = readr::locale(encoding = "UTF-8"), show_col_types = FALSE,
  progress = FALSE, name_repair = "minimal"))
normal <- rd("normal_stage_units_scores.tsv")
patient <- rd("polonen_patient_residual.tsv")
summary <- rd("subtype_residual_summary.tsv")
variance <- rd("variance_partition.tsv")
stopifnot(nrow(normal) == 20L, nrow(patient) == 1309L,
          nrow(summary) == 17L, nrow(variance) == 4L)
normalize_label <- function(x) {
  x <- as.character(x)
  x[grepl("^LMO2.*-like$", x)] <- "LMO2 gamma-delta-like"
  x[grepl("^TAL1.*-like$", x) & !grepl("DP", x)] <- "TAL1 \u03b1\u03b2-like"
  x
}
patient$label <- normalize_label(patient$subtype)
summary$label <- normalize_label(summary$subtype)
normal_range <- range(normal$dev)
curve <- unique(patient[, c("dev", "expected_balance")])
curve <- curve[order(curve$dev), ]
curve$in_range <- curve$dev >= normal_range[1] & curve$dev <= normal_range[2]

ink <- "#152b39"; muted <- "#70838e"; red <- "#bc3446"; blue <- "#246b96"
teal <- "#168b82"; gold <- "#c28830"; purple <- "#776298"
subcols <- c("BCL11B" = purple, "ETP-like" = teal, "TLX3" = gold)
font <- "Arial"
theme_pub <- theme_classic(base_family = font, base_size = 11.5) +
  theme(plot.title = element_text(colour = ink, face = "bold", size = 12.8,
                                  margin = margin(b = 5)),
        plot.subtitle = element_text(colour = ink, size = 9.5),
        plot.caption = element_text(colour = ink, size = 9),
        axis.title = element_text(colour = ink, size = 10.2),
        axis.text = element_text(colour = ink, size = 9.3),
        legend.title = element_text(colour = ink, size = 9),
        legend.text = element_text(colour = ink, size = 9),
        plot.margin = margin(5, 5, 5, 5),
        panel.grid.major.y = element_line(colour = "#e6edef", linewidth = .3))

normal$window <- factor(normal$window, levels = c("early", "mid", "late"))
pA <- ggplot() +
  geom_vline(xintercept = normal_range, colour = "#a9bbc3", linetype = "dashed", linewidth = .35) +
  geom_line(data = curve[curve$in_range, ],
            aes(dev, expected_balance), colour = ink, linewidth = 1.05) +
  geom_point(data = normal, aes(dev, balance, fill = window), shape = 21,
             size = 3.2, colour = "white", stroke = .45) +
  scale_fill_manual(values = c(early = blue, mid = teal, late = gold),
                    labels = c("Early", "Cortical / DP", "Later")) +
  annotate("text", x = 1.9, y = min(normal$balance) + .2,
           label = "R2 = 0.73 / 20 stage units", colour = ink, size = 3.1,
           hjust = 1, family = font) +
  labs(title = "A  Normal developmental reference", subtitle = "Frozen three-gene coordinate and expected balance",
       x = "Developmental coordinate", y = "ZEB1-ZEB2 balance", fill = NULL) +
  theme_pub + theme(legend.position = "bottom", legend.key.width = grid::unit(8, "pt"))

key <- patient[patient$label %in% c("BCL11B", "ETP-like", "TLX3"), ]
pB <- ggplot(patient, aes(dev, balance)) +
  annotate("rect", xmin = normal_range[1], xmax = normal_range[2],
           ymin = -Inf, ymax = Inf, fill = "#edf5f4", alpha = .7) +
  geom_point(colour = "#c8d2d5", alpha = .48, size = .55) +
  geom_line(data = curve[curve$in_range, ], aes(dev, expected_balance),
            inherit.aes = FALSE, colour = ink, linewidth = 1.0) +
  geom_line(data = curve[!curve$in_range & curve$dev < normal_range[1], ],
            aes(dev, expected_balance), inherit.aes = FALSE,
            colour = ink, linetype = "dashed", linewidth = .7) +
  geom_line(data = curve[!curve$in_range & curve$dev > normal_range[2], ],
            aes(dev, expected_balance), inherit.aes = FALSE,
            colour = ink, linetype = "dashed", linewidth = .7) +
  geom_point(data = key, aes(colour = label), alpha = .8, size = 1.2) +
  scale_colour_manual(values = subcols) +
  labs(title = "B  Diagnostic T-ALL against the reference", subtitle = "1,309 patients / shaded area is the normal observed range",
       x = "Developmental coordinate", y = "ZEB1-ZEB2 balance", colour = NULL) +
  theme_pub + theme(legend.position = "bottom", legend.key.width = grid::unit(7, "pt"))

summary <- summary[order(summary$residual_normal_median), ]
levels_order <- rev(summary$label)
patient$label <- factor(patient$label, levels = levels_order)
summary$label <- factor(summary$label, levels = levels_order)
summary$state <- ifelse(as.character(summary$label) == "BCL11B", "BCL11B",
  ifelse(summary$residual_normal_median < -.55, "ZEB2 shift",
  ifelse(summary$residual_normal_median > .5, "ZEB1 shift", "Near reference")))
patient$state <- summary$state[match(as.character(patient$label), as.character(summary$label))]
patient$subtype_n <- summary$n[match(as.character(patient$label), as.character(summary$label))]
ridge_pat <- patient[patient$subtype_n >= 10L, ]
fill_state <- c("BCL11B" = "#b1a7cf", "ZEB2 shift" = "#a8d6d0",
                "ZEB1 shift" = "#e8ce93", "Near reference" = "#dce6e9")
set.seed(2026)
pC <- ggplot(patient, aes(residual_normal, label)) +
  geom_vline(xintercept = 0, colour = "#8fa5af", linewidth = .55) +
  geom_density_ridges(data = ridge_pat, aes(fill = state), scale = .72, alpha = .76,
                      colour = "#8599a3", linewidth = .3,
                      rel_min_height = .025, bandwidth = .45) +
  geom_point(position = position_jitter(height = .095, width = 0),
             colour = ink, alpha = .24, size = .39) +
  geom_point(data = summary, aes(residual_normal_median, label),
             inherit.aes = FALSE, shape = 18, size = 2.5, colour = ink) +
  scale_fill_manual(values = fill_state, guide = "none") +
  labs(title = "C  Residual landscape by molecular subtype",
       subtitle = "Raw patients and median; density only for subtypes with n >= 10",
       x = "Actual minus normal-expected balance", y = NULL) +
  theme_pub + theme(axis.text.y = element_text(size = 10.2, colour = ink),
                    panel.grid.major.y = element_line(colour = "#edf2f3", linewidth = .25))

focus <- c("BCL11B", "SPI1", "LMO2 gamma-delta-like", "TME-enriched", "MLLT10",
           "ETP-like", "TLX1", "TAL1 DP-like", "TLX3")
sens <- summary[as.character(summary$label) %in% focus, ]
sens$label <- factor(as.character(sens$label), levels = rev(focus))
pD <- ggplot(sens, aes(y = label)) +
  geom_vline(xintercept = 0, colour = "#98abb3", linewidth = .4) +
  geom_segment(aes(x = residual_normal_median, xend = residual_within_median,
                   yend = label), colour = "#879ba4", linewidth = .7) +
  geom_point(aes(x = residual_normal_median), colour = purple, size = 2.4) +
  geom_point(aes(x = residual_within_median), colour = teal, size = 2.4, shape = 17) +
  labs(title = "D  Residual definition sensitivity",
       subtitle = "Circle: normal curve / triangle: within-cohort spline",
       caption = "BCL11B: 14/18 coordinates lie beyond the normal observed range",
       x = "Subtype median residual", y = NULL) +
  theme_pub + theme(axis.text.y = element_text(size = 8.7, colour = ink))

vp <- variance[grepl("Polonen", variance$model), ]
vp$label <- c("Development only", "Subtype only", "Combined")
vp$label <- factor(vp$label, levels = rev(c("Development only", "Subtype only", "Combined")))
vp$colour <- c(teal, purple, ink)
pE <- ggplot(vp, aes(r_squared, label)) +
  geom_segment(aes(x = 0, xend = r_squared, yend = label),
               colour = "#b9c8ce", linewidth = 1.0) +
  geom_point(aes(colour = label), size = 3.3) +
  geom_text(aes(label = sprintf("%.2f", r_squared)), hjust = -0.4,
            colour = ink, size = 3.3, family = font) +
  scale_colour_manual(values = setNames(vp$colour, as.character(vp$label)), guide = "none") +
  scale_x_continuous(limits = c(0, .51), breaks = c(0, .2, .4)) +
  labs(title = "E  Subtype adds explanatory variance",
       subtitle = "Additional R2 after development = 0.14",
       x = "Explained variance (R2)", y = NULL) +
  theme_pub + theme(axis.text.y = element_text(size = 8.7, colour = ink))

summary$label_chr <- as.character(summary$label)
summary$highlight <- ifelse(summary$label_chr %in% c("BCL11B", "SPI1", "LMO2 gamma-delta-like",
 "TME-enriched", "MLLT10", "ETP-like", "TLX3"), "Key", "Other")
pF <- ggplot(summary, aes(dev_median, residual_normal_median)) +
  geom_hline(yintercept = 0, colour = "#9aadb5", linewidth = .45) +
  geom_vline(xintercept = 0, colour = "#d4dfe3", linewidth = .35) +
  geom_point(aes(colour = highlight), size = 2.6) +
  geom_text_repel(data = subset(summary, highlight == "Key"), aes(label = label_chr),
                  colour = ink, size = 3.0, family = font, max.overlaps = Inf,
                  min.segment.length = 0, segment.color = "#9bacb3") +
  scale_colour_manual(values = c(Key = purple, Other = "#9cb0b8"), guide = "none") +
  labs(title = "F  Inheritance and reconfiguration",
       subtitle = "Subtype medians in developmental-residual space",
       x = "Median developmental coordinate", y = "Median residual") +
  theme_pub

top <- wrap_plots(pA, pB, ncol = 2, widths = c(1, 1.3))
right <- wrap_plots(pD, pE, pF, ncol = 1, heights = c(1.2, .8, 1.15))
bottom <- wrap_plots(pC, right, ncol = 2, widths = c(1.3, 1))
full <- wrap_plots(top, bottom, ncol = 1, heights = c(.8, 1.55))

wt <- function(x, name) write.table(x, file.path(sdir, name), sep = "\t",
                                    quote = FALSE, row.names = FALSE, na = "NA")
wt(normal, "F3A_normal_stage_units.tsv")
wt(curve, "F3A-B_frozen_expected_curve.tsv")
wt(patient[, c("sample_id", "subtype", "dev", "balance", "expected_balance",
               "residual_normal", "residual_within", "dev_outside_normal_range")],
   "F3B-C_patient_level.tsv")
wt(summary[, c("subtype", "n", "dev_median", "balance_median",
               "residual_normal_median", "residual_within_median",
               "n_outside_normal_dev_range", "wilcox_fdr")],
   "F3C-D-F_subtype_summary.tsv")
wt(variance, "F3E_variance_partition.tsv")

base <- file.path(out, "F3")
ggsave(paste0(base, ".pdf"), full, width = 14.4, height = 11.3,
       units = "in", device = cairo_pdf, bg = "white", limitsize = FALSE)
ggsave(paste0(base, ".svg"), full, width = 14.4, height = 11.3,
       units = "in", device = svg, bg = "white", limitsize = FALSE)
ggsave(paste0(base, ".png"), full, width = 14.4, height = 11.3,
       units = "in", device = agg_png, dpi = 300, bg = "white", limitsize = FALSE)
ggsave(paste0(base, ".tiff"), full, width = 14.4, height = 11.3,
       units = "in", device = agg_tiff, dpi = 300, bg = "white", limitsize = FALSE)
writeLines(capture.output(sessionInfo()), paste0(base, "_R_sessionInfo.txt"))
manifest <- data.frame(
  figure = "F3", panel = LETTERS[1:6],
  source = c("normal_stage_units_scores.tsv + polonen_patient_residual.tsv",
             "polonen_patient_residual.tsv", "polonen_patient_residual.tsv + subtype_residual_summary.tsv",
             "subtype_residual_summary.tsv", "variance_partition.tsv",
             "subtype_residual_summary.tsv"),
  displayed_data = c("F3A_normal_stage_units.tsv + F3A-B_frozen_expected_curve.tsv",
    "F3B-C_patient_level.tsv", "F3B-C_patient_level.tsv + F3C-D-F_subtype_summary.tsv",
    "F3C-D-F_subtype_summary.tsv", "F3E_variance_partition.tsv",
    "F3C-D-F_subtype_summary.tsv"),
  display_transform = c("Frozen spline predictions plotted without refitting",
    "All 1309 patients; normal observed range shaded; extrapolation dashed",
    "Subtype order by frozen median; density only when n >= 10; raw points retained",
    "Paired frozen subtype medians under normal and within-cohort references",
    "Frozen model R2 values shown as lollipops",
    "Frozen subtype median coordinates; labels positioned for readability"),
  script = "total/rebuild_2026/code/build_f3_ggplot.R")
write.table(manifest, paste0(base, "_panel_manifest.tsv"), sep = "\t",
            quote = FALSE, row.names = FALSE)
cat(base, "\n")

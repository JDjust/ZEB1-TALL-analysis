# Revision 2 Figure 3: within-cohort inference, normal thymus as shape anchor.
suppressPackageStartupMessages({
  library(ggplot2); library(patchwork); library(ggridges); library(ggrepel);
  library(ragg); library(readr)
})
argv <- grep("^--file=", commandArgs(FALSE), value = TRUE)
script <- normalizePath(sub("^--file=", "", argv), winslash = "/")
base <- normalizePath(file.path(dirname(script), ".."), winslash = "/")
valid <- file.path(base, "../data/validation/zeb_developmental_residual")
stats <- file.path(base, "data/source_data_rebuilt/revision2_statistics")
out <- file.path(base, "edition_20260925/figures/main")
sdir <- file.path(base, "data/source_data_rebuilt/F3")
dir.create(out, recursive = TRUE, showWarnings = FALSE)
dir.create(sdir, recursive = TRUE, showWarnings = FALSE)
rd <- function(path) as.data.frame(read_tsv(path, locale = locale(encoding = "UTF-8"),
  show_col_types = FALSE, progress = FALSE, name_repair = "minimal"))
normal <- rd(file.path(valid, "normal_stage_units_scores.tsv"))
patient <- rd(file.path(valid, "polonen_patient_residual.tsv"))
sub <- rd(file.path(stats, "subtype_effects_bootstrap.tsv"))
old <- rd(file.path(valid, "subtype_residual_summary.tsv"))
stopifnot(nrow(normal) == 20L, nrow(patient) == 1309L,
          nrow(sub) == 17L, nrow(old) == 17L)
label <- function(x) {
  x <- as.character(x)
  x[grepl("^LMO2.*-like$", x)] <- "LMO2 gamma-delta-like"
  x[grepl("^TAL1.*-like$", x) & !grepl("DP", x)] <- "TAL1 alpha-beta-like"
  x
}
sub$label <- label(sub$subtype)
patient$label <- label(patient$subtype)
old$label <- label(old$subtype)
sub <- sub[order(sub$residual_within_median), ]
ord <- rev(sub$label)
sub$label <- factor(sub$label, levels = ord)
patient$label <- factor(patient$label, levels = ord)
sub$normal_reference_median <- old$residual_normal_median[match(as.character(sub$label), old$label)]
sub$dev_median <- old$dev_median[match(as.character(sub$label), old$label)]
sub$outside_normal <- old$n_outside_normal_dev_range[match(as.character(sub$label), old$label)]
stopifnot(!anyNA(sub$normal_reference_median))

ink <- "#152b39"; muted <- "#78919c"; zeb1 <- "#bc3446"
zeb2 <- "#246b96"; teal <- "#168b82"; gold <- "#c28830"
purple <- "#776298"; font <- "Arial"
theme_pub <- theme_classic(base_family = font, base_size = 11.4) +
  theme(plot.title = element_text(colour = ink, face = "bold", size = 12.3,
                                  margin = margin(b = 4)),
        plot.subtitle = element_text(colour = ink, size = 9.2),
        plot.caption = element_text(colour = ink, size = 8.7),
        axis.title = element_text(colour = ink, size = 10),
        axis.text = element_text(colour = ink, size = 9.2),
        legend.text = element_text(colour = ink, size = 8.9),
        legend.title = element_text(colour = ink, size = 9),
        panel.grid.major.y = element_line(colour = "#e8eef0", linewidth = .28),
        plot.margin = margin(6, 7, 6, 7))
state <- function(x) ifelse(x == "BCL11B", "BCL11B",
  ifelse(x == "ETP-like", "ETP-like", ifelse(x == "TLX3", "TLX3", "Other")))
cols <- c(BCL11B = purple, `ETP-like` = teal, TLX3 = gold, Other = muted)

# Normal values are shown as within-source standardized *shape* only.
normal$window <- factor(normal$window, levels = c("early", "mid", "late"))
normal_curve <- unique(patient[, c("dev", "expected_balance")])
normal_curve <- normal_curve[order(normal_curve$dev), ]
normal_range <- range(normal$dev)
normal_curve <- subset(normal_curve, dev >= normal_range[1] & dev <= normal_range[2])
pA <- ggplot(normal, aes(dev, balance)) +
  geom_line(data = normal_curve, aes(dev, expected_balance), inherit.aes = FALSE,
            colour = ink, linewidth = .85) +
  geom_point(aes(fill = window), shape = 21, colour = "white", stroke = .4,
             size = 3.0) +
  scale_fill_manual(values = c(early = zeb2, mid = teal, late = gold),
                    labels = c("Early", "Cortical / DP", "Later")) +
  labs(title = "A  Normal thymus defines the trajectory shape",
       subtitle = "20 stage units / normal-source standardization / R² = 0.73",
       x = "Developmental coordinate", y = "ZEB1 - ZEB2 balance", fill = NULL) +
  theme_pub + theme(legend.position = "bottom")

within_curve <- unique(patient[, c("dev", "expected_within")])
within_curve <- within_curve[order(within_curve$dev), ]
patient$state <- state(as.character(patient$label))
key <- subset(patient, state != "Other")
pB <- ggplot(patient, aes(dev, balance)) +
  geom_point(colour = "#b8c8cd", alpha = .45, size = .57) +
  geom_line(data = within_curve, aes(dev, expected_within), inherit.aes = FALSE,
            colour = ink, linewidth = 1.05) +
  geom_point(data = key, aes(colour = state), alpha = .78, size = 1.05) +
  scale_colour_manual(values = cols[c("BCL11B", "ETP-like", "TLX3")]) +
  labs(title = "B  Developmental expectation fitted in T-ALL",
       subtitle = "1,309 patients / curve fitted within the same RNA cohort",
       x = "Developmental coordinate", y = "ZEB1 - ZEB2 balance", colour = NULL) +
  theme_pub + theme(legend.position = "bottom")

patient$subtype_n <- sub$n[match(as.character(patient$label), as.character(sub$label))]
ridge <- subset(patient, subtype_n >= 10)
set.seed(2026)
pC <- ggplot(patient, aes(residual_within, label)) +
  geom_vline(xintercept = 0, colour = ink, linewidth = .55) +
  geom_density_ridges(data = ridge, fill = "#c4dcd9", colour = "#8baeb2",
                      alpha = .68, linewidth = .26, scale = .68,
                      rel_min_height = .025, bandwidth = .43) +
  geom_point(position = position_jitter(width = 0, height = .09),
             colour = ink, alpha = .23, size = .37) +
  geom_point(data = sub, aes(residual_within_median, label), inherit.aes = FALSE,
             shape = 18, colour = purple, size = 2.2) +
  labs(title = "C  Within-cohort residual landscape",
       subtitle = "All patients retained; diamonds show subtype medians",
       x = "Observed minus cohort-expected balance", y = NULL) +
  scale_y_discrete(limits = ord, drop = FALSE) +
  theme_pub + theme(axis.text.y = element_text(size = 9.2, colour = ink))

sub$state <- state(as.character(sub$label))
pD <- ggplot(sub, aes(residual_within_median, label)) +
  geom_vline(xintercept = 0, colour = ink, linewidth = .55) +
  geom_segment(aes(x = residual_within_ci_low, xend = residual_within_ci_high,
                   yend = label), colour = "#9db4ba", linewidth = .8) +
  geom_point(aes(colour = state), size = 2.6) +
  scale_colour_manual(values = cols, guide = "none") +
  labs(title = "D  Residual medians with refit uncertainty",
       subtitle = "3,000 stratified patient bootstraps / percentile 95% CI",
       x = "Median within-cohort residual", y = NULL) +
  theme_pub + theme(axis.text.y = element_text(size = 9.2, colour = ink))

sub$sig <- ifelse(sub$adjusted_effect_bh_fdr < .05, "BH FDR < 0.05", "BH FDR >= 0.05")
pE <- ggplot(sub, aes(adjusted_effect, label)) +
  geom_vline(xintercept = 0, colour = ink, linewidth = .55) +
  geom_segment(aes(x = adjusted_effect_ci_low, xend = adjusted_effect_ci_high,
                   yend = label), colour = "#9db4ba", linewidth = .8) +
  geom_point(aes(fill = sig), size = 2.8, shape = 21, colour = ink, stroke = .35) +
  scale_fill_manual(values = c("BH FDR < 0.05" = purple,
                               "BH FDR >= 0.05" = "white")) +
  labs(title = "E  Joint-model subtype effects",
       subtitle = "Spline + 17 subtypes / sum contrasts / bootstrap 95% CI",
       x = "Adjusted effect vs subtype-average intercept", y = NULL,
       fill = NULL) +
  theme_pub + theme(axis.text.y = element_text(size = 9.2, colour = ink),
                    legend.position = "bottom")

focus <- c("BCL11B", "SPI1", "LMO2 gamma-delta-like", "TME-enriched",
           "MLLT10", "ETP-like", "TAL1 DP-like", "TLX3")
sen <- subset(sub, as.character(label) %in% focus)
sen$label <- factor(as.character(sen$label), levels = rev(focus))
pF <- ggplot(sen, aes(y = label)) +
  geom_vline(xintercept = 0, colour = ink, linewidth = .55) +
  geom_segment(aes(x = residual_within_median, xend = normal_reference_median,
                   yend = label), colour = "#a9bdc2", linewidth = .8) +
  geom_point(aes(x = residual_within_median), colour = teal, size = 2.65) +
  geom_point(aes(x = normal_reference_median), colour = purple,
             shape = 17, size = 2.55) +
  labs(title = "F  Cross-cohort projection is sensitivity only",
       subtitle = "Circle: within cohort  /  triangle: normal-stage projection",
       caption = "Separately standardized cohorts; BCL11B 14/18 outside normal domain",
       x = "Subtype median residual (different reference scales)", y = NULL) +
  theme_pub + theme(axis.text.y = element_text(size = 9.2, colour = ink))

top <- wrap_plots(pA, pB, ncol = 2, widths = c(1, 1.25))
mid <- wrap_plots(pC, pD, ncol = 2, widths = c(1.15, 1))
bottom <- wrap_plots(pE, pF, ncol = 2, widths = c(1.12, 1))
full <- wrap_plots(top, mid, bottom, ncol = 1, heights = c(.8, 1.43, 1.36))

wt <- function(x, name) write.table(x, file.path(sdir, name), sep = "\t",
                                    quote = FALSE, row.names = FALSE, na = "NA")
wt(normal, "F3A_normal_stage_units.tsv")
wt(normal_curve, "F3A_normal_shape_curve.tsv")
wt(patient[, c("sample_id", "subtype", "dev", "balance", "expected_within",
               "residual_within", "residual_normal", "dev_outside_normal_range")],
   "F3B-C_patient_level.tsv")
wt(within_curve, "F3B_within_cohort_curve.tsv")
wt(sub[, c("subtype", "n", "dev_median", "residual_within_median",
           "residual_within_ci_low", "residual_within_ci_high",
           "adjusted_effect", "adjusted_effect_ci_low", "adjusted_effect_ci_high",
           "adjusted_effect_hc3_p", "adjusted_effect_bh_fdr",
           "normal_reference_median", "outside_normal")],
   "F3D-F_subtype_statistics.tsv")

base_file <- file.path(out, "F3")
ggsave(paste0(base_file, ".pdf"), full, width = 14.4, height = 13.0,
       units = "in", device = cairo_pdf, bg = "white", limitsize = FALSE)
ggsave(paste0(base_file, ".svg"), full, width = 14.4, height = 13.0,
       units = "in", device = svg, bg = "white", limitsize = FALSE)
ggsave(paste0(base_file, ".png"), full, width = 14.4, height = 13.0,
       units = "in", device = agg_png, dpi = 300, bg = "white", limitsize = FALSE)
ggsave(paste0(base_file, ".tiff"), full, width = 14.4, height = 13.0,
       units = "in", device = agg_tiff, dpi = 300, bg = "white", limitsize = FALSE)
writeLines(capture.output(sessionInfo()), paste0(base_file, "_R_sessionInfo.txt"))
manifest <- data.frame(figure = "F3", panel = LETTERS[1:6],
  source = c("normal_stage_units_scores.tsv", "polonen_patient_residual.tsv",
    "polonen_patient_residual.tsv", "subtype_effects_bootstrap.tsv",
    "subtype_effects_bootstrap.tsv", "subtype_residual_summary.tsv + subtype_effects_bootstrap.tsv"),
  displayed_data = c("F3A_normal_stage_units.tsv + F3A_normal_shape_curve.tsv",
    "F3B-C_patient_level.tsv + F3B_within_cohort_curve.tsv",
    "F3B-C_patient_level.tsv", "F3D-F_subtype_statistics.tsv",
    "F3D-F_subtype_statistics.tsv", "F3D-F_subtype_statistics.tsv"),
  display_transform = c("Normal-stage standardized shape only",
    "Frozen within-cohort fitted spline and all 1309 patients",
    "Frozen within-cohort residuals; density only if n >= 10",
    "Median residual and 3000-stratified-patient bootstrap percentile 95% CI",
    "Sum-contrast adjusted effects; bootstrap CI; HC3 P values BH-corrected",
    "Within-cohort primary compared with separately standardized normal projection sensitivity"),
  script = "total/rebuild_2026/code/build_f3_revision2.R")
write.table(manifest, paste0(base_file, "_panel_manifest.tsv"), sep = "\t",
            quote = FALSE, row.names = FALSE)
cat(base_file, "\n")

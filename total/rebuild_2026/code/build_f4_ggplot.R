# Figure 4: BCL11B-rearranged convergence from frozen patient and case tables.
suppressPackageStartupMessages({
  library(ggplot2); library(patchwork); library(ragg)
})
argv <- grep("^--file=", commandArgs(FALSE), value = TRUE)
script <- normalizePath(sub("^--file=", "", argv), winslash = "/")
root <- normalizePath(file.path(dirname(script), "../../.."), winslash = "/")
valid <- file.path(root, "total/data/validation")
out <- file.path(root, "total/rebuild_2026/edition_20260925/figures/main")
sdir <- file.path(root, "total/rebuild_2026/data/source_data_rebuilt/F4")
dir.create(out, recursive = TRUE, showWarnings = FALSE)
dir.create(sdir, recursive = TRUE, showWarnings = FALSE)
rd <- function(x) as.data.frame(readr::read_tsv(file.path(valid, x),
  locale = readr::locale(encoding = "UTF-8"), show_col_types = FALSE,
  progress = FALSE, name_repair = "minimal"))
sub <- rd("zeb_developmental_residual/subtype_residual_summary.tsv")
pol <- rd("zeb_developmental_axis/bcl11b_patients.tsv")
cases <- rd("gse162280/case_manifest_and_expression.tsv")
groups <- rd("gse162280/partner_group_summary.tsv")
stopifnot(nrow(sub) == 17L, nrow(pol) == 18L, nrow(cases) == 12L)

ink <- "#152b39"; purple <- "#776298"; blue <- "#246b96"
teal <- "#168b82"; gold <- "#c28830"; zeb1 <- "#bc3446"
zeb2 <- blue; font <- "Arial"
theme_pub <- theme_classic(base_family = font, base_size = 8) +
  theme(plot.title = element_text(colour = ink, face = "bold", size = 9,
                                  margin = margin(b = 5)),
        plot.subtitle = element_text(colour = ink, size = 7.5),
        plot.caption = element_text(colour = ink, size = 7),
        axis.title = element_text(colour = ink, size = 8),
        axis.text = element_text(colour = ink, size = 7.4),
        legend.title = element_text(colour = ink, size = 7.3),
        legend.text = element_text(colour = ink, size = 7.2),
        plot.margin = margin(5, 5, 5, 5),
        panel.grid.major.y = element_line(colour = "#e7eef0", linewidth = .3))
fix_label <- function(x) {
  x <- as.character(x)
  x[grepl("^LMO2.*-like$", x)] <- "LMO2 gamma-delta-like"
  x[grepl("^TAL1.*-like$", x) & !grepl("DP", x)] <- "TAL1 alpha-beta-like"
  x
}

sub$label <- fix_label(sub$subtype)
sub <- sub[order(sub$residual_within_median), ]
sub$label <- factor(sub$label, levels = rev(sub$label))
sub$highlight <- ifelse(as.character(sub$label) == "BCL11B", "BCL11B", "Other")
pA <- ggplot(sub, aes(residual_within_median, label)) +
  geom_vline(xintercept = 0, colour = "#9caeb6", linewidth = .5) +
  geom_segment(aes(x = 0, xend = residual_within_median, yend = label),
               colour = "#c2d0d5", linewidth = .75) +
  geom_point(aes(colour = highlight), size = 2.7) +
  geom_text(data = subset(sub, highlight == "BCL11B"),
            aes(label = sprintf("%.2f  /  n=%d", residual_within_median, n)),
            hjust = -.15, colour = purple, size = 3.2, fontface = "bold", family = font) +
  scale_colour_manual(values = c(BCL11B = purple, Other = "#8fa7b0"), guide = "none") +
  scale_y_discrete(labels = function(x) ifelse(x %in%
                     c("BCL11B", "ETP-like", "TLX3"), x, "")) +
  labs(title = "A  BCL11B residual extreme",
       subtitle = "Within-cohort adjustment; 17 subtypes",
       x = "Within-cohort residual", y = NULL) +
  theme_pub + theme(axis.text.y = element_text(size = 7.3, colour = ink))

pol$case_no <- seq_len(nrow(pol))
paired <- rbind(
  data.frame(sample_id = pol$sample_id, gene = "ZEB1", value = pol$ZEB1),
  data.frame(sample_id = pol$sample_id, gene = "ZEB2", value = pol$ZEB2))
paired$gene <- factor(paired$gene, levels = c("ZEB1", "ZEB2"))
pB <- ggplot(paired, aes(gene, value, group = sample_id)) +
  geom_line(colour = "#aabcc3", linewidth = .55) +
  geom_point(aes(colour = gene), size = 2.5) +
  scale_colour_manual(values = c(ZEB1 = zeb1, ZEB2 = zeb2), guide = "none") +
  labs(title = "B  Paired ZEB expression",
       subtitle = "18 BCL11B subtype samples; TMM log2 CPM",
       x = NULL, y = "Expression") +
  theme_pub + theme(panel.grid.major.y = element_line(colour = "#e5ecee", linewidth = .3))

cases$lesion <- ifelse(cases$zeb2_partner == "yes", "ZEB2 fusion",
  ifelse(grepl("ARID1B", cases$partner), "ARID1B enhancer", "CDK6 enhancer"))
cases$lesion <- factor(cases$lesion,
  levels = c("ZEB2 fusion", "ARID1B enhancer", "CDK6 enhancer"))
lesion_cols <- c("ZEB2 fusion" = blue, "ARID1B enhancer" = teal,
                 "CDK6 enhancer" = gold)
cases$phenotype_class <- ifelse(grepl("ETP-ALL", cases$phenotype), "ETP-ALL",
  ifelse(grepl("MPAL", cases$phenotype), "MPAL",
  ifelse(grepl("AML", cases$phenotype), "AML", "Unspecified")))
cases <- cases[order(cases$table_no), ]
cases$case <- factor(cases$case, levels = rev(cases$case))
case_expr <- rbind(
  data.frame(case = cases$case, gene = "ZEB1", cpm = cases$ZEB1_cpm),
  data.frame(case = cases$case, gene = "ZEB2", cpm = cases$ZEB2_cpm))
case_expr$gene <- factor(case_expr$gene, levels = c("ZEB1", "ZEB2"))
case_labels <- setNames(paste0(as.character(cases$case), "  |  ",
                               cases$phenotype_class), as.character(cases$case))
pC <- ggplot(cases, aes(y = case)) +
  geom_vline(xintercept = c(10, 100), colour = "#eef2f3", linewidth = .3) +
  geom_segment(aes(x = ZEB1_cpm, xend = ZEB2_cpm, yend = case),
               colour = "#a5b7be", linewidth = .75) +
  geom_point(data = case_expr, aes(x = cpm, colour = gene), size = 2.35) +
  geom_point(aes(x = 9, shape = lesion), colour = ink, size = 2.0) +
  geom_text(aes(x = 680, label = sprintf("%.2f", log2_ZEB1_over_ZEB2)),
            colour = ink, size = 2.65, fontface = "bold", family = font) +
  scale_colour_manual(values = c(ZEB1 = zeb1, ZEB2 = zeb2), name = NULL) +
  scale_shape_manual(values = c("ZEB2 fusion" = 15,
                                "ARID1B enhancer" = 16,
                                "CDK6 enhancer" = 17), name = NULL) +
  scale_x_log10(limits = c(8, 900), breaks = c(10, 30, 100, 300),
                labels = c("10", "30", "100", "300")) +
  scale_y_discrete(labels = case_labels) +
  labs(title = "C  Case-level expression and lesion class",
       subtitle = "12 cases; each segment connects ZEB1 and ZEB2 in one case; right label = log2 ratio",
       x = "Expression (CPM; log scale)", y = NULL) +
  theme_pub + theme(axis.text.y = element_text(size = 7.4, colour = ink),
                    legend.position = "bottom", legend.box = "horizontal",
                    legend.margin = margin(0, 0, 0, 0),
                    legend.key.width = grid::unit(4, "mm"))

mix <- as.data.frame(table(cases$lesion, cases$phenotype_class))
names(mix) <- c("lesion", "phenotype", "n")
mix$phenotype <- factor(mix$phenotype,
  levels = c("AML", "MPAL", "ETP-ALL", "Unspecified"))
pD <- ggplot(mix, aes(phenotype, lesion, fill = n)) +
  geom_tile(colour = "white", linewidth = .9) +
  geom_text(aes(label = ifelse(n > 0, n, "")), colour = ink,
            size = 4.0, fontface = "bold", family = font) +
  scale_fill_gradient(low = "#eff4f5", high = "#82aeb8", limits = c(0, 5),
                      name = "Cases") +
  labs(title = "D  Lesion \u00d7 phenotype",
       subtitle = "Cases per cell; series includes AML and MPAL",
       x = NULL, y = NULL) +
  theme_pub + theme(panel.grid = element_blank(), axis.line = element_blank(),
                    axis.ticks = element_blank(), axis.text.x = element_text(angle = 20, hjust = 1),
                    axis.text.y = element_text(size = 7.3, colour = ink),
                    legend.position = "none")

rng <- groups[groups$group %in% c("ZEB2-BCL11B fusion",
                                   "6q25 ARID1B enhancer", "7q21 CDK6 enhancer"), ]
rng$lesion <- c("ZEB2 fusion", "ARID1B enhancer", "CDK6 enhancer")[
  match(rng$group, c("ZEB2-BCL11B fusion", "6q25 ARID1B enhancer", "7q21 CDK6 enhancer"))]
rng$lesion <- factor(rng$lesion,
  levels = rev(c("ZEB2 fusion", "ARID1B enhancer", "CDK6 enhancer")))
pE <- ggplot(rng, aes(y = lesion)) +
  geom_vline(xintercept = 0, colour = ink, linewidth = .7) +
  geom_segment(aes(x = min_log2_ratio, xend = max_log2_ratio, yend = lesion),
               linewidth = 2.0, colour = "#b7c9ce") +
  geom_point(aes(x = median_log2_ZEB1_over_ZEB2, colour = lesion), size = 3.2) +
  geom_text(aes(x = min_log2_ratio, label = paste0(n_ZEB2_counts_above_ZEB1, "/", n)),
            hjust = 1.25, colour = ink, size = 3.1, family = font) +
  scale_colour_manual(values = lesion_cols, guide = "none") +
  scale_x_continuous(limits = c(-6.2, .2), breaks = c(-6, -4, -2, 0)) +
  labs(title = "E  ZEB2 dominance by partner",
       subtitle = "Dot, median; bar, range; label, ZEB2-dominant/n",
       x = "log2(ZEB1 / ZEB2)", y = NULL) +
  theme_pub + theme(axis.text.y = element_text(size = 7.3, colour = ink))

pol$sample_id <- factor(pol$sample_id,
  levels = pol$sample_id[order(pol$balance)])
flag <- rbind(
  data.frame(sample_id = pol$sample_id, feature = "ARID1B enhancer flag", present = pol$arid),
  data.frame(sample_id = pol$sample_id, feature = "BCL11B translocation flag", present = pol$tra),
  data.frame(sample_id = pol$sample_id, feature = "MYC enhancer flag", present = pol$myc))
flag$feature <- factor(flag$feature,
  levels = c("Balance", "MYC enhancer flag", "BCL11B translocation flag",
             "ARID1B enhancer flag"))
balance_label <- data.frame(sample_id = pol$sample_id,
                            feature = factor("Balance", levels = levels(flag$feature)),
                            label = sprintf("%.1f", pol$balance))
pF <- ggplot(flag, aes(sample_id, feature)) +
  geom_tile(fill = "#f0f5f6", colour = "white", linewidth = .9) +
  geom_point(data = subset(flag, present == 1), shape = 15,
             size = 3.6, colour = purple) +
  geom_text(data = balance_label, aes(label = label),
            colour = ink, size = 3.0, fontface = "bold", family = font) +
  labs(title = "F  Recorded BCL11B subtype lesion annotations",
       subtitle = "Main cohort / 18 diagnostic cases; flagged events may overlap",
       caption = "Numbers are patient ZEB balance; filled cells are recorded lesion flags",
       x = NULL, y = NULL) +
  theme_pub + theme(panel.grid = element_blank(), axis.line = element_blank(),
                    axis.ticks = element_blank(),
                    axis.text.x = element_text(angle = 42, hjust = 1, size = 7.2),
                    axis.text.y = element_text(colour = ink, size = 7.3))

top <- wrap_plots(pA, pB, ncol = 2, widths = c(1.3, 1))
full <- top / pC / (pD | pE) / pF +
  plot_layout(heights = c(.82, 1.47, .75, .58))

wt <- function(x, name) write.table(x, file.path(sdir, name), sep = "\t",
  quote = FALSE, row.names = FALSE, na = "NA")
wt(sub[, c("subtype", "n", "residual_normal_median", "residual_within_median")],
   "F4A_subtype_residual_medians.tsv")
wt(pol, "F4B_polonen_BCL11B_patients.tsv")
wt(cases[, c("table_no", "upn", "case", "partner", "zeb2_partner", "lesion",
             "phenotype", "phenotype_class", "ZEB1", "ZEB2", "ZEB1_cpm", "ZEB2_cpm",
             "log2_ZEB1_over_ZEB2")], "F4C_GSE162280_cases.tsv")
wt(mix, "F4D_lesion_phenotype_counts.tsv")
wt(rng, "F4E_lesion_group_ranges.tsv")
wt(pol[, c("sample_id", "balance", "arid", "tra", "myc")],
   "F4F_polonen_recorded_lesion_flags.tsv")

base <- file.path(out, "F4")
ggsave(paste0(base, ".pdf"), full, width = 7.09, height = 8.3,
       units = "in", device = cairo_pdf, bg = "white", limitsize = FALSE)
ggsave(paste0(base, ".svg"), full, width = 7.09, height = 8.3,
       units = "in", device = svg, bg = "white", limitsize = FALSE)
ggsave(paste0(base, ".png"), full, width = 7.09, height = 8.3,
       units = "in", device = agg_png, dpi = 300, bg = "white", limitsize = FALSE)
ggsave(paste0(base, ".tiff"), full, width = 7.09, height = 8.3,
       units = "in", device = agg_tiff, dpi = 300, bg = "white", limitsize = FALSE)
writeLines(capture.output(sessionInfo()), paste0(base, "_R_sessionInfo.txt"))
manifest <- data.frame(
  figure = "F4", panel = LETTERS[1:6],
  source = c("zeb_developmental_residual/subtype_residual_summary.tsv",
             "zeb_developmental_axis/bcl11b_patients.tsv",
             "gse162280/case_manifest_and_expression.tsv",
             "gse162280/case_manifest_and_expression.tsv",
             "gse162280/partner_group_summary.tsv",
             "zeb_developmental_axis/bcl11b_patients.tsv"),
  displayed_data = c("F4A_subtype_residual_medians.tsv",
                     "F4B_polonen_BCL11B_patients.tsv", "F4C_GSE162280_cases.tsv",
                     "F4D_lesion_phenotype_counts.tsv", "F4E_lesion_group_ranges.tsv",
                     "F4F_polonen_recorded_lesion_flags.tsv"),
  display_transform = c("Within-cohort subtype medians ordered; BCL11B highlighted",
    "Paired frozen ZEB1/ZEB2 TMM log2 CPM per diagnostic patient",
    "Frozen case-level CPM ratio; partner and phenotype labels arranged",
    "Case counts grouped by partner and reported phenotype",
    "Frozen observed range, median and dominance counts",
    "Recorded lesion flags and frozen balance per P\u00f6l\u00f6nen BCL11B patient"),
  script = "total/rebuild_2026/code/build_f4_ggplot.R")
write.table(manifest, paste0(base, "_panel_manifest.tsv"), sep = "\t",
            quote = FALSE, row.names = FALSE)
cat(base, "\n")

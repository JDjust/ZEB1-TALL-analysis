# Figure 5: independent bulk expression polarity and Lim patient-level state localization.
suppressPackageStartupMessages({library(ggplot2); library(patchwork); library(ragg)})
argv <- grep("^--file=", commandArgs(FALSE), value = TRUE)
script <- normalizePath(sub("^--file=", "", argv), winslash = "/")
root <- normalizePath(file.path(dirname(script), "../../.."), winslash = "/")
valid <- file.path(root, "total/data/validation")
out <- file.path(root, "total/rebuild_2026/edition_20260925/figures/main")
sdir <- file.path(root, "total/rebuild_2026/data/source_data_rebuilt/F5")
dir.create(out, recursive = TRUE, showWarnings = FALSE)
dir.create(sdir, recursive = TRUE, showWarnings = FALSE)
rd <- function(x) as.data.frame(readr::read_tsv(file.path(valid, x),
  locale = readr::locale(encoding = "UTF-8"), show_col_types = FALSE,
  progress = FALSE, name_repair = "minimal"))
bulk <- rd("gse146901/zeb_expression_by_sample.tsv")
bulk_stats <- rd("gse146901/zeb_expression_etp_vs_nonetp.tsv")
pairs <- rd("gateA_lim2025/q3_state_day0_pairs.tsv")
state_stats <- rd("gateA_lim2025/q3_state.tsv")
lim <- rd("gateA_lim2025/eligible_patient_timepoint.tsv")
baseline_stats <- rd("gateA_lim2025/q1_baseline.tsv")
stopifnot(nrow(bulk) == 22L, sum(bulk$group == "ETP") == 8L,
          sum(bulk$group == "non-ETP") == 10L,
          nrow(pairs) == 41L)
day0 <- lim[lim$timepoint == "Day0", ]
stopifnot(nrow(day0) == 54L, sum(day0$induction == "IF") == 19L,
          sum(day0$induction == "responsive") == 35L)

ink <- "#152b39"; zeb1 <- "#bc3446"; zeb2 <- "#246b96"
teal <- "#168b82"; gold <- "#c28830"; purple <- "#776298"
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

bulk$group <- factor(bulk$group, levels = c("ETP", "non-ETP", "normal T"))
group_cols <- c("ETP" = teal, "non-ETP" = gold, "normal T" = "#9baeb7")
set.seed(2026)
pA <- ggplot(bulk, aes(group, balance)) +
  geom_hline(yintercept = 0, colour = "#a9b9c0", linewidth = .45) +
  geom_boxplot(width = .42, outlier.shape = NA, fill = "white",
               colour = "#8299a3", linewidth = .5) +
  geom_point(aes(colour = group), position = position_jitter(width = .12),
             size = 2.4, alpha = .9) +
  scale_colour_manual(values = group_cols, guide = "none") +
  scale_x_discrete(labels = c("ETP", "non-ETP", "healthy\nperipheral T")) +
  labs(title = "A  ETP/non-ETP balance",
       subtitle = "GSE146901 / 8 ETP, 10 non-ETP; MW P = 0.0117",
       x = NULL, y = "z(ZEB1) - z(ZEB2)") +
  theme_pub

pB <- ggplot(bulk, aes(group, log2_ratio)) +
  geom_hline(yintercept = 0, colour = "#a9b9c0", linewidth = .45) +
  geom_boxplot(width = .42, outlier.shape = NA, fill = "white",
               colour = "#8299a3", linewidth = .5) +
  geom_point(aes(colour = group), position = position_jitter(width = .12),
             size = 2.4, alpha = .9) +
  scale_colour_manual(values = group_cols, guide = "none") +
  scale_x_discrete(labels = c("ETP", "non-ETP", "healthy\nperipheral T")) +
  labs(title = "B  Expression ratio",
       subtitle = "log2(ZEB1/ZEB2); MW P = 0.0155",
       x = NULL, y = "log2 expression ratio") +
  theme_pub

bulk_leuk <- bulk[bulk$group %in% c("ETP", "non-ETP"), ]
bulk_long <- rbind(data.frame(sample = bulk_leuk$sample, group = bulk_leuk$group,
                              gene = "ZEB1", value = bulk_leuk$log2_ZEB1),
                   data.frame(sample = bulk_leuk$sample, group = bulk_leuk$group,
                              gene = "ZEB2", value = bulk_leuk$log2_ZEB2))
bulk_long$gene <- factor(bulk_long$gene, levels = c("ZEB1", "ZEB2"))
pC <- ggplot(bulk_long, aes(gene, value, group = sample)) +
  geom_line(aes(colour = group), alpha = .42, linewidth = .7) +
  geom_point(aes(fill = gene), shape = 21, colour = "white", size = 2.5,
             stroke = .25) +
  facet_wrap(~group, nrow = 1) +
  scale_colour_manual(values = group_cols, guide = "none") +
  scale_fill_manual(values = c(ZEB1 = zeb1, ZEB2 = zeb2), guide = "none") +
  labs(title = "C  Paired genes per sample",
       subtitle = "Within-sample pairs; leukemia samples only",
       x = NULL, y = "log2 expression") +
  theme_pub + theme(strip.background = element_blank(),
                    strip.text = element_text(face = "bold", colour = ink))

pairs$induction <- factor(pairs$induction_pos, levels = c("IF", "responsive"))
pair_long <- rbind(data.frame(sample = pairs$biological_sample_id, induction = pairs$induction,
                              state = "ZBTB16+", balance = pairs$balance_pos),
                   data.frame(sample = pairs$biological_sample_id, induction = pairs$induction,
                              state = "ZBTB16-", balance = pairs$balance_neg))
pair_long$state <- factor(pair_long$state, levels = c("ZBTB16-", "ZBTB16+"))
pD <- ggplot(pair_long, aes(state, balance, group = sample)) +
  geom_line(colour = "#9badb5", linewidth = .5, alpha = .75) +
  geom_point(aes(colour = state), size = 2.1, alpha = .87) +
  scale_colour_manual(values = c("ZBTB16-" = "#8fa5ad", "ZBTB16+" = purple), guide = "none") +
  labs(title = "D  ZBTB16 marks lower-balance blast states",
       subtitle = "Lim / 41 paired Day0 malignant pseudobulks; paired P = 0.00395",
       x = NULL, y = "ZEB balance") +
  theme_pub

pE <- ggplot(pairs, aes(induction, d_balance)) +
  geom_hline(yintercept = 0, colour = "#a9b9c0", linewidth = .45) +
  geom_boxplot(width = .38, outlier.shape = NA, fill = "white",
               colour = "#8299a3", linewidth = .5) +
  geom_point(aes(colour = induction), position = position_jitter(width = .10),
             size = 2, alpha = .84) +
  scale_colour_manual(values = c("IF" = purple, "responsive" = teal), guide = "none") +
  labs(title = "E  State contrast by induction group",
       subtitle = "ZBTB16+ minus ZBTB16-; 18 IF and 23 responsive pairs",
       x = NULL, y = "Paired balance difference") +
  theme_pub

day0$induction <- factor(day0$induction, levels = c("IF", "responsive"))
pF <- ggplot(day0, aes(induction, balance)) +
  geom_hline(yintercept = 0, colour = "#a9b9c0", linewidth = .45) +
  geom_boxplot(width = .42, outlier.shape = NA, fill = "white",
               colour = "#8299a3", linewidth = .5) +
  geom_point(aes(colour = etp_status), position = position_jitter(width = .12),
             size = 2.2, alpha = .88) +
  scale_colour_manual(values = c("ETP" = teal, "nonETP" = gold), name = "ETP status") +
  labs(title = "F  Baseline balance and ETP composition",
       subtitle = "Day0 patient pseudobulk / 19 IF, 35 responsive",
       x = NULL, y = "ZEB balance") +
  theme_pub + theme(legend.position = "bottom", legend.key.width = grid::unit(9, "pt"))

comp <- as.data.frame(table(day0$induction, day0$etp_status))
names(comp) <- c("induction", "etp_status", "n")
comp$etp_status <- factor(comp$etp_status, levels = c("ETP", "nonETP"))
pG <- ggplot(comp, aes(induction, n, fill = etp_status)) +
  geom_col(width = .55, colour = "white", linewidth = .5) +
  geom_text(aes(label = n), position = position_stack(vjust = .5),
            colour = ink, fontface = "bold", size = 3.4, family = font) +
  scale_fill_manual(values = c("ETP" = teal, "nonETP" = gold), name = "ETP status") +
  labs(title = "G  ETP-state composition at diagnosis",
       subtitle = "Frozen patient-level counts, not cell counts",
       x = NULL, y = "Patients") +
  theme_pub + theme(legend.position = "bottom", legend.key.width = grid::unit(9, "pt"))

top <- wrap_plots(pA, pB, pC, ncol = 3, widths = c(1, 1, 1.15))
middle <- wrap_plots(pD, pE, ncol = 2, widths = c(1.35, 1))
bottom <- wrap_plots(pF, pG, ncol = 2, widths = c(1.25, 1))
full <- wrap_plots(top, middle, bottom, ncol = 1, heights = c(1.0, 1.0, .86))
wt <- function(x, name) write.table(x, file.path(sdir, name), sep = "\t",
  quote = FALSE, row.names = FALSE, na = "NA")
wt(bulk, "F5A-C_GSE146901_samples.tsv")
wt(bulk_stats, "F5A-B_GSE146901_frozen_tests.tsv")
wt(pairs[, c("biological_sample_id", "patient_id_pos", "cohort_pos", "induction_pos",
             "etp_status_pos", "n_cells_pos", "n_cells_neg", "balance_pos", "balance_neg",
             "d_balance")], "F5D-E_Lim_Day0_state_pairs.tsv")
wt(day0[, c("biological_sample_id", "patient_id", "cohort", "induction",
             "etp_status", "subtype_l1", "n_cells", "balance")],
   "F5F_Lim_Day0_patient_pseudobulk.tsv")
wt(comp, "F5G_Lim_induction_ETP_counts.tsv")
wt(state_stats[state_stats$metric == "balance", ], "F5D-E_Lim_frozen_state_tests.tsv")
wt(baseline_stats[baseline_stats$metric == "balance", ], "F5F_Lim_frozen_baseline_tests.tsv")

base <- file.path(out, "F5")
ggsave(paste0(base, ".pdf"), full, width = 14.4, height = 11.2,
       units = "in", device = cairo_pdf, bg = "white", limitsize = FALSE)
ggsave(paste0(base, ".svg"), full, width = 14.4, height = 11.2,
       units = "in", device = svg, bg = "white", limitsize = FALSE)
ggsave(paste0(base, ".png"), full, width = 14.4, height = 11.2,
       units = "in", device = agg_png, dpi = 300, bg = "white", limitsize = FALSE)
ggsave(paste0(base, ".tiff"), full, width = 14.4, height = 11.2,
       units = "in", device = agg_tiff, dpi = 300, bg = "white", limitsize = FALSE)
writeLines(capture.output(sessionInfo()), paste0(base, "_R_sessionInfo.txt"))
manifest <- data.frame(
  figure = "F5", panel = LETTERS[1:7],
  source = c(rep("gse146901/zeb_expression_by_sample.tsv", 3),
    rep("gateA_lim2025/q3_state_day0_pairs.tsv", 2),
    rep("gateA_lim2025/eligible_patient_timepoint.tsv", 2)),
  displayed_data = c("F5A-C_GSE146901_samples.tsv", "F5A-C_GSE146901_samples.tsv",
    "F5A-C_GSE146901_samples.tsv", "F5D-E_Lim_Day0_state_pairs.tsv",
    "F5D-E_Lim_Day0_state_pairs.tsv", "F5F_Lim_Day0_patient_pseudobulk.tsv",
    "F5G_Lim_induction_ETP_counts.tsv"),
  display_transform = c("Raw patient values; frozen MW comparison; healthy controls are peripheral T",
    "Raw patient values; frozen MW comparison; healthy controls are peripheral T",
    "Frozen within-sample ZEB1/ZEB2 values paired",
    "41 frozen Day0 ZBTB16 positive/negative malignant-pseudobulk pairs",
    "Frozen paired balance differences split by induction group",
    "54 frozen Day0 patient pseudobulks; ETP status marked",
    "Patient-level ETP/nonETP counts by induction response"),
  script = "total/rebuild_2026/code/build_f5_ggplot.R")
write.table(manifest, paste0(base, "_panel_manifest.tsv"), sep = "\t",
            quote = FALSE, row.names = FALSE)
cat(base, "\n")

# Figure 2: frozen 1,309-patient subtype configurations and phenotype mapping.
suppressPackageStartupMessages({
  library(ggplot2); library(patchwork); library(ggalluvial); library(ragg)
})
argv <- grep("^--file=", commandArgs(FALSE), value = TRUE)
script <- normalizePath(sub("^--file=", "", argv), winslash = "/")
root <- normalizePath(file.path(dirname(script), "../../.."), winslash = "/")
valid <- file.path(root, "total/data/validation")
out <- file.path(root, "total/rebuild_2026/edition_20260925/figures/main")
sdir <- file.path(root, "total/rebuild_2026/data/source_data_rebuilt/F2")
dir.create(out, recursive = TRUE, showWarnings = FALSE)
dir.create(sdir, recursive = TRUE, showWarnings = FALSE)
rd <- function(x) as.data.frame(readr::read_tsv(file.path(valid, x),
  locale = readr::locale(encoding = "UTF-8"), show_col_types = FALSE,
  progress = FALSE, name_repair = "minimal"))
patients <- rd("polonen_round1/patient_level_frozen_balance.tsv")
resid <- rd("zeb_developmental_residual/polonen_patient_residual.tsv")
land <- rd("polonen_round1/subtype_landscape.tsv")
sub <- rd("zeb_developmental_residual/subtype_residual_summary.tsv")
ip <- rd("zeb_developmental_axis/ip_by_subtype_counts.tsv")
variance <- rd("zeb_developmental_residual/variance_partition.tsv")
ip_variance <- rd("zeb_developmental_axis/polonen_balance_variance.tsv")
stopifnot(nrow(patients) == 1309L, nrow(resid) == 1309L,
          nrow(sub) == 17L, sum(ip$N) == 1309L)
idx <- match(patients$sample_id, resid$sample_id)
stopifnot(!anyNA(idx), all(patients$subtype == resid$subtype[idx]))
patients$dev <- resid$dev[idx]

clean <- function(x) {
  x <- as.character(x)
  x[grepl("^LMO2.*-like$", x)] <- "LMO2 gamma-delta-like"
  x[grepl("^TAL1.*-like$", x) & !grepl("DP", x)] <- "TAL1 alpha-beta-like"
  x
}
patients$label <- clean(patients$subtype)
sub$label <- clean(sub$subtype)
land$label <- clean(land$subtype)
ip$label <- clean(ip$subtype)
ip$ip_label <- ifelse(ip$ip == "DP-like", "DP-like IP",
  ifelse(ip$ip == "ETP-like", "ETP-like IP",
  ifelse(ip$ip == "Myeloid-like", "Myeloid-like IP",
  ifelse(ip$ip == "Other", "Other / unknown IP", "alpha-beta-like IP"))))

ink <- "#152b39"; purple <- "#776298"; teal <- "#168b82"
gold <- "#c28830"; zeb1 <- "#bc3446"; zeb2 <- "#246b96"
font <- "Arial"
sub_palette <- c(
 "BCL11B" = purple, "SPI1" = "#8565a9", "LMO2 gamma-delta-like" = "#a17d9b",
 "TME-enriched" = "#477f9e", "NUP98" = "#8498aa", "MLLT10" = "#5a9a9a",
 "ETP-like" = teal, "HOXA9 TCR" = "#8a9b77", "STAG2&LMO2" = "#b0a279",
 "TLX1" = "#b97b65", "KMT2A" = "#6d8799", "NKX2-5" = "#95a3a7",
 "TAL1 alpha-beta-like" = "#bd777b", "NKX2-1" = "#927d7d",
 "TAL1 DP-like" = "#c5626e", "NUP214" = "#929c87", "TLX3" = gold)
stopifnot(all(unique(patients$label) %in% names(sub_palette)))
theme_pub <- theme_classic(base_family = font, base_size = 11.5) +
  theme(plot.title = element_text(colour = ink, face = "bold", size = 12.8,
                                  margin = margin(b = 5)),
        plot.subtitle = element_text(colour = ink, size = 9.5),
        plot.caption = element_text(colour = ink, size = 9),
        axis.title = element_text(colour = ink, size = 10.2),
        axis.text = element_text(colour = ink, size = 9.3),
        legend.text = element_text(colour = ink, size = 9),
        plot.margin = margin(5, 5, 5, 5),
        panel.grid.major.y = element_line(colour = "#e8edef", linewidth = .3))

sub <- sub[order(sub$balance_median), ]
sub_order <- sub$label
ip$label <- factor(ip$label, levels = sub_order)
ip$ip_label <- factor(ip$ip_label,
  levels = c("DP-like IP", "ETP-like IP", "Myeloid-like IP",
             "alpha-beta-like IP", "Other / unknown IP"))
stopifnot(!anyNA(ip$ip_label))
# Adapted from D:/_bioinformation/code1/27桑基图/bioR27.ggalluvial.R:
# its to_lodes_form grammar is applied to
# frozen patient counts. The tutorial's setwd, sample data and palette are
# deliberately replaced by project provenance and the manuscript palette.
ip$flow_id <- seq_len(nrow(ip))
lodes <- ggalluvial::to_lodes_form(
  ip[, c('ip_label', 'label', 'N', 'flow_id')], axes = 1:2,
  id = 'flow_id', key = 'axis', value = 'stratum')
lodes$flow_subtype <- ip$label[match(lodes$flow_id, ip$flow_id)]
pA <- ggplot(lodes, aes(x = axis, stratum = stratum,
                       alluvium = flow_id, y = N)) +
  geom_alluvium(aes(fill = flow_subtype), width = .09, alpha = .54,
                colour = NA, knot.pos = .38) +
  geom_stratum(width = .09, fill = '#ffffff', colour = '#708c99', linewidth = .48) +
  geom_text(stat = "stratum", aes(label = after_stat(ifelse(
                       (x == 1 & ymax - ymin >= 250) |
                       (x == 2 & ymax - ymin >= 100),
                       as.character(stratum), ""))), size = 2.9,
            colour = ink, family = font) +
  scale_fill_manual(values = sub_palette, guide = "none") +
  scale_x_discrete(limits = c('ip_label','label'),
                   labels = c('Immunophenotype','Molecular subtype'),
                   expand = c(.18, .12)) +
  labs(title = "A  Immunophenotypes distribute across molecular subtypes",
       subtitle = "Observed patient counts; all 47 nonzero flows retained",
       x = NULL, y = "Patients") +
  theme_pub + theme(axis.line.x = element_blank(), axis.ticks.x = element_blank(),
                    panel.grid = element_blank())

heat <- land[, c("metric", "label", "median")]
heat$metric <- ifelse(heat$metric == "balance", "Balance", heat$metric)
extra <- data.frame(metric = "Development", label = sub$label,
                    median = sub$dev_median)
heat <- rbind(heat, extra)
heat$z <- ave(heat$median, heat$metric, FUN = function(x) as.numeric(scale(x)))
heat$label <- factor(heat$label, levels = rev(sub_order))
heat$metric <- factor(heat$metric,
  levels = c("ZEB1", "ZEB2", "LMO2", "Balance", "Development"))
pB <- ggplot(heat, aes(metric, label, fill = z)) +
  geom_tile(colour = "white", linewidth = .75) +
  scale_fill_gradient2(low = "#5482a2", mid = "#f7f8f5", high = "#bf6170",
                       midpoint = 0, limits = c(-2.5, 2.5),
                       oob = scales::squish, name = "Across-subtype z") +
  labs(title = "B  Configuration matrix",
       subtitle = "Subtype medians, z-scored per metric",
       x = NULL, y = NULL, fill = "z") +
  theme_pub + theme(panel.grid = element_blank(), axis.line = element_blank(),
                    axis.ticks = element_blank(),
                    axis.text.x = element_text(angle = 30, hjust = 1, size = 9.2),
                    axis.text.y = element_text(size = 9.2, colour = ink),
                    legend.position = "right", legend.key.width = grid::unit(7, "pt"),
                    legend.key.height = grid::unit(22, "pt"))

patients$label <- factor(patients$label, levels = rev(sub_order))
q <- do.call(rbind, lapply(split(patients, patients$label), function(d) {
  data.frame(label = as.character(d$label[1]), n = nrow(d),
    q25 = as.numeric(quantile(d$balance, .25)),
    median = median(d$balance), q75 = as.numeric(quantile(d$balance, .75)))
}))
q$label <- factor(q$label, levels = rev(sub_order))
set.seed(2026)
pC <- ggplot(patients, aes(balance, label)) +
  geom_vline(xintercept = 0, colour = "#9eb0b8", linewidth = .5) +
  geom_point(position = position_jitter(height = .15, width = 0),
             colour = "#667f8b", alpha = .32, size = .5) +
  geom_segment(data = q, aes(x = q25, xend = q75, y = label, yend = label),
               inherit.aes = FALSE, colour = ink, linewidth = 2.6, alpha = .84) +
  geom_point(data = q, aes(x = median, y = label), inherit.aes = FALSE,
             shape = 18, size = 2.5, colour = purple) +
  labs(title = "C  Balance by subtype",
       subtitle = "Dots, patients; bar, IQR; diamond, median",
       x = "z(ZEB1) - z(ZEB2)", y = NULL) +
  theme_pub + theme(axis.text.y = element_text(size = 9.2, colour = ink))

patients$key <- ifelse(patients$label %in% c("BCL11B", "ETP-like", "TLX3"),
                       as.character(patients$label), "Other")
pD <- ggplot(patients, aes(z_ZEB1, z_ZEB2)) +
  geom_abline(slope = 1, intercept = 0, colour = "#a8b9c1", linetype = "dashed") +
  geom_point(data = subset(patients, key == "Other"), colour = "#abbcc2",
             size = .55, alpha = .35) +
  geom_point(data = subset(patients, key != "Other"), aes(colour = key),
             size = 1.25, alpha = .8) +
  scale_colour_manual(values = c("BCL11B" = purple, "ETP-like" = teal,
                                 "TLX3" = gold)) +
  coord_equal() +
  labs(title = "E  ZEB1 versus ZEB2",
       subtitle = "One point per patient; equality dashed",
       x = "ZEB1, within-cohort z-score", y = "ZEB2, within-cohort z-score",
       colour = NULL) +
  theme_pub + theme(legend.position = "bottom", legend.key.width = grid::unit(9, "pt"))

dev <- sub[, c("label", "n", "dev_median")]
dev$label <- factor(dev$label, levels = rev(sub_order))
dev$highlight <- ifelse(as.character(dev$label) %in% c("BCL11B", "ETP-like", "TLX3"),
                        as.character(dev$label), "Other")
pE <- ggplot(dev, aes(dev_median, label)) +
  geom_vline(xintercept = 0, colour = "#a9bbc2", linewidth = .5) +
  geom_segment(aes(x = 0, xend = dev_median, yend = label),
               colour = "#bfd0d5", linewidth = .7) +
  geom_point(aes(colour = highlight), size = 2.3) +
  scale_colour_manual(values = c("BCL11B" = purple, "ETP-like" = teal,
                                 "TLX3" = gold, "Other" = "#869fa9"), guide = "none") +
  labs(title = "D  Development",
       subtitle = "Median coordinate",
       x = "Three-gene coordinate", y = NULL) +
  theme_pub + theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(),
                    panel.grid.major.y = element_line(colour = "#e8edef", linewidth = .3))

vp <- variance[grepl("Polonen", variance$model), ]
model <- data.frame(label = c("Immunophenotype", "Development", "Molecular subtype",
                              "Development + subtype"),
                    r2 = c(ip_variance$r_squared[ip_variance$model == "author immunophenotype (4)"],
                           vp$r_squared[grepl("ns\\(dev", vp$model) & !grepl("subtype", vp$model)],
                           vp$r_squared[grepl("balance ~ subtype", vp$model)],
                           vp$r_squared[grepl("subtype", vp$model) & grepl("ns\\(dev", vp$model)]))
stopifnot(nrow(model) == 4L, all(is.finite(model$r2)))
model$label <- factor(model$label, levels = rev(model$label))
pF <- ggplot(model, aes(r2, label)) +
  geom_segment(aes(x = 0, xend = r2, yend = label), colour = "#bfd0d5", linewidth = .9) +
  geom_point(colour = purple, size = 3.0) +
  geom_text(aes(label = sprintf("%.2f", r2)), hjust = -.3, colour = ink,
            size = 3.3, family = font) +
  scale_x_continuous(limits = c(0, .5), breaks = c(0, .2, .4)) +
  labs(title = "F  Variance explained",
       subtitle = "Subtype adds \u0394R\u00b2 = 0.14 after development",
       x = expression("Explained variance ("*R^2*")"), y = NULL) +
  theme_pub

top <- pA
middle <- wrap_plots(pB, pC, pE, ncol = 3, widths = c(1.18, 1.4, .82))
bottom <- wrap_plots(pD, pF, ncol = 2, widths = c(1, 1))
full <- wrap_plots(top, middle, bottom, ncol = 1, heights = c(.8, 1.23, .82))
wt <- function(x, name) write.table(x, file.path(sdir, name), sep = "\t",
  quote = FALSE, row.names = FALSE, na = "NA")
wt(ip[, c("ip", "subtype", "N")], "F2A_immunophenotype_subtype_counts.tsv")
wt(heat[, c("metric", "label", "median", "z")], "F2B_subtype_matrix.tsv")
wt(patients[, c("sample_id", "subtype", "ZEB1", "ZEB2", "LMO2",
                  "z_ZEB1", "z_ZEB2", "balance", "dev")], "F2C-D_patients.tsv")
wt(q, "F2C_subtype_IQR.tsv")
wt(dev[, c("label", "n", "dev_median")], "F2D_developmental_medians.tsv")
wt(model, "F2F_model_variance.tsv")
wt(data.frame(subtype = names(sub_palette), colour = unname(sub_palette)),
   "F2_subtype_palette.tsv")

base <- file.path(out, "F2")
ggsave(paste0(base, ".pdf"), full, width = 14.4, height = 12.6,
       units = "in", device = cairo_pdf, bg = "white", limitsize = FALSE)
ggsave(paste0(base, ".svg"), full, width = 14.4, height = 12.6,
       units = "in", device = svg, bg = "white", limitsize = FALSE)
ggsave(paste0(base, ".png"), full, width = 14.4, height = 12.6,
       units = "in", device = agg_png, dpi = 300, bg = "white", limitsize = FALSE)
ggsave(paste0(base, ".tiff"), full, width = 14.4, height = 12.6,
       units = "in", device = agg_tiff, dpi = 300, bg = "white", limitsize = FALSE)
writeLines(capture.output(sessionInfo()), paste0(base, "_R_sessionInfo.txt"))
manifest <- data.frame(
  figure = "F2", panel = LETTERS[1:6],
  source = c("zeb_developmental_axis/ip_by_subtype_counts.tsv",
    "polonen_round1/subtype_landscape.tsv + zeb_developmental_residual/subtype_residual_summary.tsv",
    "polonen_round1/patient_level_frozen_balance.tsv",
    "zeb_developmental_residual/subtype_residual_summary.tsv",
    "polonen_round1/patient_level_frozen_balance.tsv",
    "zeb_developmental_residual/variance_partition.tsv + zeb_developmental_axis/polonen_balance_variance.tsv"),
  displayed_data = c("F2A_immunophenotype_subtype_counts.tsv",
    "F2B_subtype_matrix.tsv", "F2C-D_patients.tsv + F2C_subtype_IQR.tsv",
    "F2D_developmental_medians.tsv", "F2C-D_patients.tsv", "F2F_model_variance.tsv"),
  display_transform = c("47 frozen counts via code1 to_lodes_form grammar; left labels >=250, right labels >=100",
    "Frozen metric medians standardized across subtypes for heatmap colour",
    "All 1309 raw balances; display-only quartiles and medians",
    "Frozen subtype median developmental coordinate, aligned to B-C order",
    "All 1309 within-cohort z(ZEB1), z(ZEB2) pairs; three subtypes highlighted",
    "Frozen R2 models shown as lollipops"),
  script = "total/rebuild_2026/code/build_f2_ggplot.R")
write.table(manifest, paste0(base, "_panel_manifest.tsv"), sep = "\t",
            quote = FALSE, row.names = FALSE)
cat(base, "\n")

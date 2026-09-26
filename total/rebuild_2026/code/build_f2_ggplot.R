# Figure 2: frozen 1,309-patient subtype configurations and phenotype mapping.
suppressPackageStartupMessages({
  library(ggplot2); library(patchwork); library(ragg)
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
 "TAL1 alpha-beta-like" = "#b45e42", "NKX2-1" = "#927d7d",
 "TAL1 DP-like" = "#d4743e", "NUP214" = "#929c87", "TLX3" = gold)
stopifnot(all(unique(patients$label) %in% names(sub_palette)))
theme_pub <- theme_classic(base_family = font, base_size = 8) +
  theme(plot.title = element_text(colour = ink, face = "bold", size = 9,
                                  margin = margin(b = 5)),
        plot.subtitle = element_text(colour = ink, size = 7.5),
        plot.caption = element_text(colour = ink, size = 7),
        axis.title = element_text(colour = ink, size = 8),
        axis.text = element_text(colour = ink, size = 7.4),
        legend.text = element_text(colour = ink, size = 7.2),
        plot.margin = margin(5, 5, 5, 5),
        panel.grid.major.y = element_line(colour = "#e8edef", linewidth = .3))

sub <- sub[order(sub$balance_median), ]
sub_order <- sub$label
ip$label <- factor(ip$label, levels = sub_order)
ip$ip_label <- factor(ip$ip_label,
  levels = c("DP-like IP", "ETP-like IP", "Myeloid-like IP",
             "alpha-beta-like IP", "Other / unknown IP"))
stopifnot(!anyNA(ip$ip_label))
# A count matrix is legible at final print size and preserves the 47 observed
# transitions. At 17 target subtypes the previous alluvial hid target labels.
ip_grid <- expand.grid(label = sub_order,
                       ip_label = levels(ip$ip_label), stringsAsFactors = FALSE)
ip_grid <- merge(ip_grid, ip[, c("label", "ip_label", "N")],
                 by = c("label", "ip_label"), all.x = TRUE, sort = FALSE)
ip_grid$N[is.na(ip_grid$N)] <- 0L
ip_grid$total <- sub$n[match(as.character(ip_grid$label), as.character(sub$label))]
ip_grid$fraction <- ip_grid$N / ip_grid$total
ip_grid$label <- factor(ip_grid$label, levels = rev(sub_order))
ip_grid$ip_label <- factor(ip_grid$ip_label, levels = levels(ip$ip_label),
  labels = c("DP", "ETP", "Myeloid", "alpha-beta", "Other"))
stopifnot(sum(ip_grid$N) == 1309L, sum(ip_grid$N > 0) == 47L)
pA <- ggplot(ip_grid, aes(ip_label, label, fill = fraction)) +
  geom_tile(colour = "#ffffff", linewidth = .35) +
  geom_text(aes(label = ifelse(N > 0, N, "")), size = 2.05,
            colour = ink, family = font) +
  scale_fill_gradientn(colours = c("#f3f7f8", "#b9d7d7", "#5ba3a2", "#1a6c78"),
                       limits = c(0, 1), guide = "none") +
  labs(title = "A  Phenotype counts",
       subtitle = "47 nonzero cells / n printed",
       x = NULL, y = NULL) +
  theme_pub + theme(panel.grid = element_blank(), axis.line = element_blank(),
                    axis.ticks = element_blank(),
                    axis.text.x = element_text(angle = 35, hjust = 1, size = 6.8),
                    axis.text.y = element_text(size = 6.9, colour = ink))

heat <- land[, c("metric", "label", "median")]
heat$metric <- ifelse(heat$metric == "balance", "Balance", heat$metric)
extra <- data.frame(metric = "Development", label = sub$label,
                    median = sub$dev_median)
heat <- rbind(heat, extra)
heat$z <- ave(heat$median, heat$metric, FUN = function(x) as.numeric(scale(x)))
heat$label <- factor(heat$label, levels = rev(sub_order))
heat$metric <- factor(heat$metric,
  levels = c("ZEB1", "ZEB2", "LMO2", "Balance", "Development"))
heat_ann <- data.frame(label = factor(sub$label, levels = rev(sub_order)),
                       n = sub$n)
stopifnot(sum(heat_ann$n) == 1309L, nrow(heat_ann) == 17L)
pB <- ggplot(heat, aes(metric, label, fill = z)) +
  geom_tile(colour = "#e4e9e9", linewidth = .28) +
  geom_text(data = heat_ann, aes(x = 5.65, y = label, label = n),
            inherit.aes = FALSE, colour = ink, size = 2.55,
            family = font) +
  scale_fill_gradientn(colours = c("#315f82", "#a7c0ce", "#f6f3ec",
                                     "#e4a47d", "#b84922"),
                       values = scales::rescale(c(-2, -1, 0, 1, 2)),
                       limits = c(-2, 2),
                       oob = scales::squish, name = "Across-subtype z") +
  labs(title = "B  Expression matrix",
       subtitle = "Subtype medians; right: n",
       x = NULL, y = NULL, fill = "z") +
  theme_pub + theme(panel.grid = element_blank(), axis.line = element_blank(),
                    axis.ticks = element_blank(),
                    axis.text.x = element_text(angle = 30, hjust = 1, size = 7.4),
                    axis.text.y = element_blank(),
                    axis.ticks.y = element_blank(),
                    legend.position = "right", legend.key.width = grid::unit(7, "pt"),
                    legend.key.height = grid::unit(22, "pt")) +
  scale_x_discrete(expand = expansion(add = c(.45, 1.10)))

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
               inherit.aes = FALSE, colour = ink, linewidth = 1.15, alpha = .84) +
  geom_point(data = q, aes(x = median, y = label), inherit.aes = FALSE,
             shape = 18, size = 2.5, colour = purple) +
  labs(title = "C  ZEB balance",
       subtitle = "1,309 patients",
       x = "z(ZEB1) - z(ZEB2)", y = NULL) +
  theme_pub + theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())

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
       subtitle = "Patients; within-cohort z-scores",
       x = "ZEB1 z-score", y = "ZEB2 z-score",
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
  labs(title = "D  Coordinate",
       subtitle = "Subtype medians",
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
model$fill <- as.character(model$label)
pF <- ggplot(model, aes(r2, label)) +
  geom_col(aes(fill = fill), width = .58, show.legend = FALSE) +
  geom_text(aes(label = sprintf("%.1f%%", 100 * r2)), hjust = -.15,
            colour = ink, size = 3.5, family = font) +
  scale_fill_manual(values = c("Development + subtype" = purple, "Development" = teal,
                               "Molecular subtype" = "#7a9aac", "Immunophenotype" = "#bbcbd0")) +
  scale_x_continuous(limits = c(0, .48), breaks = c(0, .2, .4),
                     labels = function(x) paste0(round(x * 100), "%")) +
  labs(title = "F  Explained variance",
       subtitle = "+13.9 points beyond development",
       x = NULL, y = NULL) +
  theme_pub + theme(panel.grid.major.y = element_blank(),
                    panel.grid.major.x = element_line(colour = "#e8edef", linewidth = .3))

row1 <- pA + pB + pC + plot_layout(widths = c(1.35, 1.25, 1.0))
full <- row1 /
  (pE | pD) /
  pF +
  plot_layout(heights = c(2.0, 1.15, .56))
wt <- function(x, name) write.table(x, file.path(sdir, name), sep = "\t",
  quote = FALSE, row.names = FALSE, na = "NA")
wt(ip[, c("ip", "subtype", "N")], "F2A_immunophenotype_subtype_counts.tsv")
wt(heat[, c("metric", "label", "median", "z")], "F2B_subtype_matrix.tsv")
wt(heat_ann, "F2B_subtype_counts.tsv")
wt(patients[, c("sample_id", "subtype", "ZEB1", "ZEB2", "LMO2",
                  "z_ZEB1", "z_ZEB2", "balance", "dev")], "F2C-D_patients.tsv")
wt(q, "F2C_subtype_IQR.tsv")
wt(dev[, c("label", "n", "dev_median")], "F2D_developmental_medians.tsv")
wt(model, "F2F_model_variance.tsv")
wt(data.frame(subtype = names(sub_palette), colour = unname(sub_palette)),
   "F2_subtype_palette.tsv")

base <- file.path(out, "F2")
ggsave(paste0(base, ".pdf"), full, width = 7.09, height = 8.15,
       units = "in", device = cairo_pdf, bg = "white", limitsize = FALSE)
ggsave(paste0(base, ".svg"), full, width = 7.09, height = 8.15,
       units = "in", device = svg, bg = "white", limitsize = FALSE)
ggsave(paste0(base, ".png"), full, width = 7.09, height = 8.15,
       units = "in", device = agg_png, dpi = 300, bg = "white", limitsize = FALSE)
ggsave(paste0(base, ".tiff"), full, width = 7.09, height = 8.15,
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
    "F2B_subtype_matrix.tsv + F2B_subtype_counts.tsv", "F2C-D_patients.tsv + F2C_subtype_IQR.tsv",
    "F2D_developmental_medians.tsv", "F2C-D_patients.tsv", "F2F_model_variance.tsv"),
  display_transform = c("47 frozen nonzero phenotype-by-subtype counts in a 17x5 matrix; shade is within-subtype fraction; zero cells remain blank",
    "Frozen metric medians standardized across subtypes; colour display capped at z +/-2; frozen subtype n printed at right",
    "All 1309 raw balances; display-only quartiles and medians",
    "Frozen subtype median developmental coordinate, aligned to B-C order",
    "All 1309 within-cohort z(ZEB1), z(ZEB2) pairs; three subtypes highlighted",
    "Frozen R2 models shown as direct-labelled bars"),
  script = "total/rebuild_2026/code/build_f2_ggplot.R")
write.table(manifest, paste0(base, "_panel_manifest.tsv"), sep = "\t",
            quote = FALSE, row.names = FALSE)
cat(base, "\n")

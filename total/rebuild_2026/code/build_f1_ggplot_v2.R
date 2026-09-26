# Figure 1, rebuilt from frozen normal-thymus tables and saved cell embeddings.
# Display-only transformations: within-atlas stage z-scores for expression tracks.
suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
  library(ragg)
})

argv <- grep("^--file=", commandArgs(FALSE), value = TRUE)
stopifnot(length(argv) == 1L)
script <- normalizePath(sub("^--file=", "", argv), winslash = "/")
root <- normalizePath(file.path(dirname(script), "../../.."), winslash = "/")
vdir <- file.path(root, "total/data/validation")
out <- file.path(root, "total/rebuild_2026/edition_20260925/figures/main")
dir.create(out, recursive = TRUE, showWarnings = FALSE)
rd <- function(x) read.delim(file.path(vdir, x), check.names = FALSE)
traj <- rd("zeb_developmental_axis/normal_thymus_balance_trajectory.tsv")
pattern <- rd("zeb_developmental_axis/trajectory_pattern_calls.tsv")
marker <- rd("zeb_developmental_residual/gene_direction_by_dataset.tsv")
park <- read.delim(file.path(root, "total/data/SourceData_ParkHTA_donor_celltype.tsv"))
cells <- read.delim(gzfile(file.path(root,
  "modules/module2/tables/M2_r2_GSE195812_plot_cells.tsv.gz")),
  check.names = FALSE)
cells <- cells[cells$lineage == "T" & is.finite(cells$UMAP1) & is.finite(cells$UMAP2), ]

ink <- "#152b39"; zeb1 <- "#bc3446"; zeb2 <- "#246b96"
teal <- "#168b82"; gold <- "#ca8a32"; purple <- "#776298"
font <- "Arial"
theme_pub <- theme_classic(base_family = font, base_size = 11.5) +
  theme(plot.title = element_text(size = 12.5, face = "bold", colour = ink,
                                  margin = margin(b = 5)),
        plot.subtitle = element_text(size = 9.2, colour = ink),
        plot.caption = element_text(size = 9, colour = ink),
        axis.title = element_text(size = 10.2, colour = ink),
        axis.text = element_text(size = 9.5, colour = ink),
        legend.title = element_text(size = 9.3, colour = ink),
        legend.text = element_text(size = 9.0, colour = ink),
        plot.margin = margin(4, 5, 4, 5))

stage_levels <- c("DN1", "DN2", "DN3", "ISP", "DP_CD3neg", "DP_CD3pos", "CD4SP", "CD8SP")
stage_cols <- c("#315d79", "#4e7b91", "#72a0aa", "#a57caa",
                "#c17b76", "#d69b66", "#b35460", "#8c4c69")
cells$facs_stage <- factor(cells$facs_stage, levels = stage_levels)
set.seed(2026)
cells <- cells[sample.int(nrow(cells)), ]
umap_base <- theme_pub + theme(axis.line = element_blank(), axis.ticks = element_blank(),
  axis.text = element_blank(), panel.border = element_rect(colour = "#d8e3e7", fill = NA, linewidth = .35))

pA <- ggplot(cells, aes(UMAP1, UMAP2)) +
  geom_point(aes(colour = facs_stage), size = .56, alpha = .90, stroke = 0) +
  scale_colour_manual(values = setNames(stage_cols, stage_levels),
    labels = c("DN1", "DN2", "DN3", "ISP", "DP CD3-", "DP CD3+", "CD4 SP", "CD8 SP")) +
  guides(colour = guide_legend(override.aes = list(size = 2.2, alpha = 1), ncol = 4)) +
  labs(title = "A  Stage-labeled thymocyte embedding", subtitle = "GSE195812 / T-lineage cells; original embedding",
       x = "UMAP 1", y = "UMAP 2", colour = NULL) +
  umap_base + theme(legend.position = "bottom", legend.spacing.x = grid::unit(2, "pt"),
                    legend.key.width = grid::unit(9, "pt"),
                    legend.margin = margin(t = -2, b = 0))

make_feature <- function(gene, letter, colour) {
  vals <- cells[[paste0("g_", gene)]]
  d <- cells[order(vals), ]; d$expr <- vals[order(vals)]
  ggplot(d, aes(UMAP1, UMAP2)) +
    geom_point(colour = "#dce4e7", size = .33, alpha = .50, stroke = 0) +
    geom_point(data = d[d$expr > 0, ], aes(colour = expr), size = .72,
               alpha = .94, stroke = 0) +
    scale_colour_gradient(low = "#efd8d0", high = colour, name = "log1p expression") +
    labs(title = paste0(letter, "  ", gene, " cellular expression"),
         subtitle = "Same embedding / detected cells highlighted",
         x = "UMAP 1", y = NULL) +
    umap_base + theme(legend.position = "bottom",
                      legend.key.width = grid::unit(18, "pt"),
                      legend.key.height = grid::unit(3, "pt")) +
    guides(colour = guide_colourbar(title.position = "top", title.hjust = 0))
}
pB <- make_feature("ZEB1", "B", zeb1)
pC <- make_feature("ZEB2", "C", zeb2)

zstage <- function(x) {
  if (all(is.na(x)) || sd(x, na.rm = TRUE) == 0) return(rep(0, length(x)))
  as.numeric(scale(x))
}
traj <- traj[order(traj$dataset, traj$order), ]
traj$z1 <- ave(traj$ZEB1, traj$dataset, FUN = zstage)
traj$z2 <- ave(traj$ZEB2, traj$dataset, FUN = zstage)
tracks <- rbind(
  data.frame(dataset = traj$dataset, order = traj$order, stage = traj$stage,
             track = "ZEB1", value = traj$z1),
  data.frame(dataset = traj$dataset, order = traj$order, stage = traj$stage,
             track = "ZEB2", value = traj$z2),
  data.frame(dataset = traj$dataset, order = traj$order, stage = traj$stage,
             track = "Balance", value = traj$balance)
)
tracks$track <- factor(tracks$track, levels = c("Balance", "ZEB2", "ZEB1"))
tracks$dataset <- factor(tracks$dataset,
  levels = c("GSE142522", "GSE195812", "GSE206710", "Park/HTA"))
# This is each atlas's recorded stage order, not inferred pseudotime.
tracks$stage_key <- paste(tracks$dataset, sprintf("%02d", tracks$order), tracks$stage, sep = "|")
tracks$stage_key <- factor(tracks$stage_key,
  levels = unique(tracks$stage_key[order(tracks$dataset, tracks$order)]))
pD <- ggplot(tracks, aes(stage_key, track, fill = value)) +
  geom_tile(width = .93, height = .86, colour = "white", linewidth = .6) +
  geom_text(data = subset(tracks, track == "Balance"),
            aes(label = sprintf("%+.1f", value)), size = 3.0,
            colour = ink, fontface = "bold", family = font) +
  facet_wrap(~dataset, ncol = 2, scales = "free_x") +
  scale_fill_gradient2(low = "#4c789b", mid = "#f7f7f5", high = "#bf5260",
                       midpoint = 0, limits = c(-2.4, 2.4),
                       oob = scales::squish, name = "Within-atlas\nstandardized value") +
  scale_x_discrete(labels = function(x) gsub("_", " ", sub("^.*\\|", "", x)),
                   expand = expansion(add = .5)) +
  labs(title = "D  Four atlases recover the stage-patterned axis",
       subtitle = "Stage labels are cohort-specific; printed values are frozen ZEB1-ZEB2 balance",
       x = NULL, y = NULL) +
  theme_pub + theme(panel.grid = element_blank(), axis.line = element_blank(),
                    axis.ticks = element_blank(),
                    strip.background = element_blank(),
                    strip.text = element_text(face = "bold", colour = ink, size = 10.5),
                    axis.text.x = element_text(angle = 28, hjust = 1, size = 9.2),
                    legend.position = "right")

pat <- rbind(
  data.frame(dataset = pattern$dataset, window = "Early", balance = pattern$balance_early),
  data.frame(dataset = pattern$dataset, window = "Cortical / DP", balance = pattern$balance_mid),
  data.frame(dataset = pattern$dataset, window = "Later", balance = pattern$balance_late)
)
pat$window <- factor(pat$window, levels = c("Early", "Cortical / DP", "Later"))
pE <- ggplot(pat, aes(window, balance, group = dataset, colour = dataset)) +
  geom_hline(yintercept = 0, colour = "#b5c6cd", linewidth = .4) +
  geom_line(linewidth = .7) + geom_point(size = 2.2) +
  scale_colour_manual(values = c("GSE142522" = zeb1, "GSE195812" = teal,
                                 "GSE206710" = gold, "Park/HTA" = purple)) +
  labs(title = "E  Cross-atlas cortical / DP elevation", x = NULL,
       y = "Frozen balance", colour = NULL) +
  theme_pub + theme(legend.position = "bottom", legend.key.width = grid::unit(7, "pt"),
                    axis.text.x = element_text(angle = 20, hjust = 1),
                    panel.grid.major.y = element_line(colour = "#e4ebee", linewidth = .3)) +
  guides(colour = guide_legend(nrow = 2))

don <- park[park$celltype %in% c("double negative thymocyte",
                   "CD4-positive, alpha-beta T cell") & park$n_cells >= 20, ]
don$stage <- ifelse(don$celltype == "double negative thymocyte", "DN", "CD4 SP")
wide <- reshape(don[, c("donor", "stage", "ZEB1")], idvar = "donor",
                timevar = "stage", direction = "wide")
wide <- wide[complete.cases(wide), ]
don <- don[don$donor %in% wide$donor, ]
don$stage <- factor(don$stage, levels = c("DN", "CD4 SP"))
pF <- ggplot(don, aes(stage, ZEB1, group = donor)) +
  geom_line(colour = "#90a4ad", linewidth = .44) +
  geom_point(aes(colour = stage), size = 2.1) +
  scale_colour_manual(values = c("DN" = zeb1, "CD4 SP" = teal), guide = "none") +
  labs(title = "E  Matched-donor stage contrast", subtitle = paste0("Park/HTA / ", nrow(wide), " donors"),
       x = NULL, y = "Mean ZEB1 / donor") +
  theme_pub + theme(panel.grid.major.y = element_line(colour = "#e4ebee", linewidth = .3))

mk <- expand.grid(gene = c("CD34", "LYL1", "CD1A"),
                  dataset = c("GSE195812", "GSE206710", "Park/HTA", "GSE142522"))
key <- paste(marker$gene, marker$dataset, sep = "|")
mk$delta <- marker$mid_minus_early[match(paste(mk$gene, mk$dataset, sep = "|"), key)]
mk$gene <- factor(mk$gene, levels = c("CD1A", "LYL1", "CD34"))
mk$dataset <- factor(mk$dataset,
  levels = c("GSE195812", "GSE206710", "Park/HTA", "GSE142522"))
mk$direction <- ifelse(is.na(mk$delta), "NA", ifelse(mk$delta > 0, "+", "-"))
lim <- max(abs(mk$delta), na.rm = TRUE)
pG <- ggplot(mk, aes(dataset, gene, fill = delta)) +
  geom_tile(colour = "white", linewidth = .8) +
  geom_text(aes(label = direction), size = 3.5, colour = ink,
            family = font, fontface = "bold") +
  scale_fill_gradient2(low = "#688aa5", mid = "#f7f7f5", high = "#cb7b73",
                       midpoint = 0, limits = c(-lim, lim), na.value = "#e9eef0",
                       name = "Mid - early") +
  labs(title = "F  Three-gene developmental coordinate",
       subtitle = "Signs across the four normal references",
       caption = "z(CD1A) - mean[z(CD34), z(LYL1)]; ZEB genes excluded",
       x = NULL, y = NULL) +
  theme_pub + theme(panel.grid = element_blank(), axis.line = element_blank(),
    axis.ticks = element_blank(), axis.text.x = element_text(angle = 40, hjust = 1, size = 7.4),
    legend.position = "none")

top <- wrap_plots(pA, pB, pC, ncol = 3, widths = c(1.4, 1, 1))
bottom <- wrap_plots(pF, pG, ncol = 2, widths = c(1, 1.35))
full <- wrap_plots(top, pD, bottom, ncol = 1, heights = c(1.2, 1.15, .8))

source_dir <- file.path(root, "total/rebuild_2026/data/source_data_rebuilt/F1")
dir.create(source_dir, recursive = TRUE, showWarnings = FALSE)
write_tsv <- function(x, name) write.table(x, file.path(source_dir, name),
  sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
write_tsv(cells[, c("sample", "gsm", "facs_stage", "lineage", "UMAP1", "UMAP2",
                     "g_ZEB1", "g_ZEB2")], "F1A-C_displayed_cells.tsv")
write_tsv(tracks[, c("dataset", "stage", "order", "track", "value")],
          "F1D_stage_tracks.tsv")
write_tsv(don[, c("donor", "stage", "ZEB1", "n_cells")],
          "F1E_matched_donor_ZEB1.tsv")
write_tsv(mk[, c("gene", "dataset", "delta", "direction")],
          "F1F_score_marker_directions.tsv")
write_tsv(pattern, "F1_cross_atlas_window_summary.tsv")

base <- file.path(out, "F1")
ggsave(paste0(base, ".pdf"), full, width = 14.4, height = 10.2,
       units = "in", device = cairo_pdf, bg = "white", limitsize = FALSE)
ggsave(paste0(base, ".svg"), full, width = 14.4, height = 10.2,
       units = "in", device = svg, bg = "white", limitsize = FALSE)
ggsave(paste0(base, ".png"), full, width = 14.4, height = 10.2,
       units = "in", device = agg_png, dpi = 300, bg = "white", limitsize = FALSE)
ggsave(paste0(base, ".tiff"), full, width = 14.4, height = 10.2,
       units = "in", device = agg_tiff, dpi = 300, bg = "white", limitsize = FALSE)
writeLines(capture.output(sessionInfo()), paste0(base, "_R_sessionInfo.txt"))
manifest <- data.frame(
  figure = "F1", panel = LETTERS[1:6],
  source = c("modules/module2/tables/M2_r2_GSE195812_plot_cells.tsv.gz",
             "modules/module2/tables/M2_r2_GSE195812_plot_cells.tsv.gz",
             "modules/module2/tables/M2_r2_GSE195812_plot_cells.tsv.gz",
             "total/data/validation/zeb_developmental_axis/normal_thymus_balance_trajectory.tsv",
             "total/data/SourceData_ParkHTA_donor_celltype.tsv",
             "total/data/validation/zeb_developmental_residual/gene_direction_by_dataset.tsv"),
  displayed_data = c(rep("F1A-C_displayed_cells.tsv", 3),
                     "F1D_stage_tracks.tsv", "F1E_matched_donor_ZEB1.tsv",
                     "F1F_score_marker_directions.tsv"),
  display_transform = c(rep("Existing UMAP; T-lineage cells only; no re-embedding", 3),
    "Stage-order heatmap; within-atlas ZEB1/ZEB2 z; frozen balance unchanged",
    "Matched donors with >=20 cells per DN and CD4 SP stage",
    "Frozen middle-minus-early direction; NA means unavailable"),
  script = "total/rebuild_2026/code/build_f1_ggplot_v2.R")
write.table(manifest, paste0(base, "_panel_manifest.tsv"), sep = "\t",
            quote = FALSE, row.names = FALSE)
cat(base, "\n")

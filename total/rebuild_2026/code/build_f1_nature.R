# Figure 1 — Nature/Blood-grade rebuild using theme_nature.R
# ZEB1-ZEB2 axis across normal human thymopoiesis
# Outputs: F1.pdf/png/tiff/svg at 180mm (full two-column width)

suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
  library(ragg)
})

argv <- grep("^--file=", commandArgs(FALSE), value = TRUE)
stopifnot(length(argv) == 1L)
script <- normalizePath(sub("^--file=", "", argv), winslash = "/")
root   <- normalizePath(file.path(dirname(script), "../../.."), winslash = "/")
source(file.path(dirname(script), "theme_nature.R"))

vdir <- file.path(root, "total/data/validation")
out  <- file.path(root, "total/rebuild_2026/edition_20260925/figures/main")
dir.create(out, recursive = TRUE, showWarnings = FALSE)

rd <- function(x) read.delim(file.path(vdir, x), check.names = FALSE)
traj    <- rd("zeb_developmental_axis/normal_thymus_balance_trajectory.tsv")
pattern <- rd("zeb_developmental_axis/trajectory_pattern_calls.tsv")
marker  <- rd("zeb_developmental_residual/gene_direction_by_dataset.tsv")
park    <- read.delim(file.path(root, "total/data/SourceData_ParkHTA_donor_celltype.tsv"))
cells   <- read.delim(gzfile(file.path(root,
  "modules/module2/tables/M2_r2_GSE195812_plot_cells.tsv.gz")),
  check.names = FALSE)
cells <- cells[cells$lineage == "T" & is.finite(cells$UMAP1) & is.finite(cells$UMAP2), ]

# ── Panel A: Stage-labeled UMAP ─────────────────────────────────────────────
stage_levels <- c("DN1","DN2","DN3","ISP","DP_CD3neg","DP_CD3pos","CD4SP","CD8SP")
stage_labels <- c("DN1","DN2","DN3","ISP","DP CD3\u2212","DP CD3\u207a","CD4 SP","CD8 SP")
stage_cols   <- c("#3b5998","#5b7bb3","#7da5b8","#9e7bb0",
                  "#c47570","#d99c5e","#b34e5e","#8c4c69")
cells$facs_stage <- factor(cells$facs_stage, levels = stage_levels)
set.seed(2026); cells <- cells[sample.int(nrow(cells)), ]

pA <- ggplot(cells, aes(UMAP1, UMAP2)) +
  geom_point(aes(colour = facs_stage), size = 0.28, alpha = 0.9, stroke = 0,
             shape = 16) +
  scale_colour_manual(values = setNames(stage_cols, stage_levels),
                      labels = stage_labels,
                      name = NULL) +
  guides(colour = guide_legend(override.aes = list(size = 1.8, alpha = 1),
                               ncol = 4, keywidth = unit(6, "pt"))) +
  labs(x = "UMAP 1", y = "UMAP 2") +
  theme_nature_umap() +
  theme(legend.position = "bottom",
        legend.text = element_text(size = 5.5),
        legend.spacing.x = unit(1, "pt"),
        legend.margin = margin(t = -2))

# ── Panels B,C: Feature UMAPs ──────────────────────────────────────────────
make_feature <- function(gene, colour_high, colour_low) {
  vals <- cells[[paste0("g_", gene)]]
  d <- cells[order(vals), ]
  d$expr <- vals[order(vals)]
  bg_pts <- d[d$expr == 0, ]
  fg_pts <- d[d$expr > 0, ]
  
  ggplot(mapping = aes(UMAP1, UMAP2)) +
    geom_point(data = bg_pts, colour = "#e4e4ea", size = 0.12,
               alpha = 0.45, stroke = 0, shape = 16) +
    geom_point(data = fg_pts, aes(colour = expr), size = 0.38,
               alpha = 0.92, stroke = 0, shape = 16) +
    scale_colour_gradient(low = colour_low, high = colour_high,
                          name = paste0(gene, "\nlog1p")) +
    labs(x = "UMAP 1", y = NULL) +
    theme_nature_umap() +
    theme(legend.position = "bottom",
          legend.key.width = unit(14, "pt"),
          legend.key.height = unit(3, "pt"),
          legend.title = element_text(size = 5.5),
          legend.text = element_text(size = 5)) +
    guides(colour = guide_colourbar(title.position = "left",
                                    title.hjust = 0.5))
}
pB <- make_feature("ZEB1", pal$zeb1, "#f0bfb7")
pC <- make_feature("ZEB2", pal$zeb2, "#b9d0e6")

# ── Panel D: Four-atlas stage heatmap ──────────────────────────────────────
zstage <- function(x) {
  if (all(is.na(x)) || sd(x, na.rm = TRUE) == 0) return(rep(0, length(x)))
  as.numeric(scale(x))
}
traj <- traj[order(traj$dataset, traj$order), ]
traj$z1 <- ave(traj$ZEB1, traj$dataset, FUN = zstage)
traj$z2 <- ave(traj$ZEB2, traj$dataset, FUN = zstage)

tracks <- rbind(
  data.frame(dataset = traj$dataset, order = traj$order, stage = traj$stage,
             track = "ZEB1", value = traj$z1, stringsAsFactors = FALSE),
  data.frame(dataset = traj$dataset, order = traj$order, stage = traj$stage,
             track = "ZEB2", value = traj$z2, stringsAsFactors = FALSE),
  data.frame(dataset = traj$dataset, order = traj$order, stage = traj$stage,
             track = "Balance", value = traj$balance, stringsAsFactors = FALSE)
)
tracks$track <- factor(tracks$track, levels = c("Balance","ZEB2","ZEB1"))
tracks$dataset <- factor(tracks$dataset,
  levels = c("GSE142522","GSE195812","GSE206710","Park/HTA"))
tracks$stage_key <- paste(tracks$dataset, sprintf("%02d", tracks$order),
                          tracks$stage, sep = "|")
tracks$stage_key <- factor(tracks$stage_key,
  levels = unique(tracks$stage_key[order(tracks$dataset, tracks$order)]))

pD <- ggplot(tracks, aes(stage_key, track, fill = value)) +
  geom_tile(width = 0.90, height = 0.82, colour = "white", linewidth = 0.4) +
  geom_text(data = subset(tracks, track == "Balance"),
            aes(label = sprintf("%+.1f", value)), size = 1.9,
            colour = pal$ink, fontface = "bold", family = "Arial") +
  facet_wrap(~dataset, ncol = 2, scales = "free_x") +
  scale_fill_diverging(limits = c(-2.5, 2.5),
                       name = "Within-atlas\nstandardized") +
  scale_x_discrete(labels = function(x) gsub("_"," ", sub("^.*\\|","",x)),
                   expand = expansion(add = 0.4)) +
  labs(x = NULL, y = NULL) +
  theme_nature_heatmap() +
  theme(strip.text = element_text(face = "bold", size = 6.5, colour = pal$ink),
        axis.text.x = element_text(angle = 30, hjust = 1, size = 5.5),
        axis.text.y = element_text(size = 6),
        legend.key.width = unit(6, "pt"),
        legend.key.height = unit(16, "pt"),
        legend.position = "right")

# ── Panel E: Matched-donor ZEB1 ───────────────────────────────────────────
don <- park[park$celltype %in% c("double negative thymocyte",
                                  "CD4-positive, alpha-beta T cell") &
              park$n_cells >= 20, ]
don$stage <- ifelse(don$celltype == "double negative thymocyte", "DN", "CD4 SP")
wide <- reshape(don[, c("donor","stage","ZEB1")], idvar = "donor",
                timevar = "stage", direction = "wide")
wide <- wide[complete.cases(wide), ]
don <- don[don$donor %in% wide$donor, ]
don$stage <- factor(don$stage, levels = c("DN","CD4 SP"))

pE <- ggplot(don, aes(stage, ZEB1, group = donor)) +
  geom_line(colour = "#a0aab0", linewidth = 0.35) +
  geom_point(aes(colour = stage), size = 1.4, shape = 16) +
  scale_colour_manual(values = c("DN" = pal$zeb1, "CD4 SP" = pal$etp),
                      guide = "none") +
  labs(x = NULL, y = "Mean ZEB1 / donor",
       subtitle = paste0("Park/HTA, n = ", nrow(wide), " donors")) +
  theme_nature() +
  theme(panel.grid.major.y = element_line(colour = pal$grid, linewidth = 0.2))

# ── Panel F: Coordinate marker directions ──────────────────────────────────
mk <- expand.grid(
  gene = c("CD34","LYL1","CD1A"),
  dataset = c("GSE195812","GSE206710","Park/HTA","GSE142522"),
  stringsAsFactors = FALSE)
key <- paste(marker$gene, marker$dataset, sep = "|")
mk$delta <- marker$mid_minus_early[match(paste(mk$gene, mk$dataset, sep="|"), key)]
mk$gene    <- factor(mk$gene, levels = c("CD1A","LYL1","CD34"))
mk$dataset <- factor(mk$dataset,
  levels = c("GSE195812","GSE206710","Park/HTA","GSE142522"))
mk$direction <- ifelse(is.na(mk$delta), "NA",
                       ifelse(mk$delta > 0, "+", "\u2212"))
lim <- max(abs(mk$delta), na.rm = TRUE)

pF <- ggplot(mk, aes(dataset, gene, fill = delta)) +
  geom_tile(colour = "white", linewidth = 0.6) +
  geom_text(aes(label = direction), size = 2.5, colour = pal$ink,
            family = "Arial", fontface = "bold") +
  scale_fill_gradient2(low = pal$div_low, mid = pal$div_mid, high = pal$div_high,
                       midpoint = 0, limits = c(-lim, lim),
                       na.value = "#e9eef0", guide = "none") +
  labs(x = NULL, y = NULL,
       caption = "z(CD1A) \u2212 mean[z(CD34), z(LYL1)]; ZEB genes excluded") +
  theme_nature_heatmap() +
  theme(axis.text.x = element_text(angle = 35, hjust = 1, size = 5.5),
        axis.text.y = element_text(size = 6),
        plot.caption = element_text(size = 5))

# ── Assembly ───────────────────────────────────────────────────────────────
pA <- tag_plot(pA, "A", "GSE195812 T-lineage cells")
pB <- tag_plot(pB, "B", "Detected cells highlighted")
pC <- tag_plot(pC, "C", "Detected cells highlighted")
pD <- tag_plot(pD, "D", "Printed values are ZEB1\u2013ZEB2 balance")
pE <- tag_plot(pE, "E", paste0("Park/HTA matched donors, n = ", nrow(wide)))
pF <- tag_plot(pF, "F", "Coordinate marker direction")

top <- wrap_elements(full = wrap_plots(pA, pB, pC, ncol = 3, widths = c(1.35, 1, 1)))
bot <- wrap_elements(full = wrap_plots(pE, pF, ncol = 2, widths = c(1, 1.25)))
full <- wrap_plots(top, wrap_elements(full = pD), bot, ncol = 1,
                   heights = c(1.15, 1.15, 0.78))

# ── Export ─────────────────────────────────────────────────────────────────
save_nature(full, file.path(out, "F1"),
            width = 7.09, height = 8.5, dpi = 600)

writeLines(capture.output(sessionInfo()),
           file.path(out, "F1_R_sessionInfo.txt"))
cat("F1 exported at 180mm / 600dpi Nature format\n")

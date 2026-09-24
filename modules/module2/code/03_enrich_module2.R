# Module 2 figure enrichment from code2 (not another violin/lollipop pass).
# Recipes:
#   【189】 bump / rank trajectory
#   【66】  heatmap + line combo
#   【40】  PCA + density + ellipse
#   【53】  correlation tile
#   【67】  dual-axis bar + line
#   【78】  alluvial / sankey
#   【150】 radial bar
#   【191】 dumbbell
#   【193】 scatter + lm equation
#   【206】/【209】 volcano
#   【243】 ternary (ggtern if present)
#   【244】 polar radar
#   【88】  coefficient forest
#   【346】 bidirectional bar
#   【210】 annotated heatmap

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})
has_bump <- requireNamespace("ggbump", quietly = TRUE)
has_alluvial <- requireNamespace("ggalluvial", quietly = TRUE)
has_tern <- requireNamespace("ggtern", quietly = TRUE)
has_aplot <- requireNamespace("aplot", quietly = TRUE)
has_patch <- requireNamespace("patchwork", quietly = TRUE)
has_cow <- requireNamespace("cowplot", quietly = TRUE)
has_cht <- requireNamespace("ComplexHeatmap", quietly = TRUE)
has_circlize <- requireNamespace("circlize", quietly = TRUE)
has_ggrepel <- requireNamespace("ggrepel", quietly = TRUE)
has_beeswarm <- requireNamespace("ggbeeswarm", quietly = TRUE)
if (has_beeswarm) library(ggbeeswarm)
if (has_patch) library(patchwork)
if (has_aplot) library(aplot)
if (has_cow) library(cowplot)

root <- "/data-b/liangfuhua/projects/ZEB1_TALL_analysis/module2"
fig <- file.path(root, "figures")
tab <- file.path(root, "tables")
proc <- file.path(root, "processed")
dir.create(fig, FALSE, TRUE)
dir.create(tab, FALSE, TRUE)

pal_gene <- c(ZEB1 = "#E64B35", ZEB2 = "#3C5488", LMO2 = "#F39B7F", IL7R = "#00A087",
              CD34 = "#7E6148", CD1A = "#4DBBD5", TAL1 = "#8491B4", LYL1 = "#91D1C2",
              KIT = "#B09C85", SPI1 = "#DC0000", RAG1 = "#91D1C2", DNTT = "#3C5488")
pal_stage <- c("CD34+" = "#E64B35", ISP = "#F39B7F", "Early cortical" = "#00A087",
               "Late cortical" = "#3C5488", SP4 = "#4DBBD5", SP8 = "#8491B4")
pal_ip <- c(Immature = "#E64B35", Pre_ab_Cortical = "#00A087", TCR_Mature = "#3C5488")
pal_group <- c(sorted_thymocyte = "#4DBBD5", total_thymus = "#91D1C2", "adult_T-ALL" = "#E64B35")
pal_sub <- c(TAL = "#E64B35", TLX = "#4DBBD5", HOXA = "#00A087",
             "LMO2/LYL1" = "#3C5488", NKX2_1 = "#F39B7F", Unknown = "#8491B4")

theme_sci <- function(base_size = 11) {
  theme_classic(base_size = base_size) +
    theme(
      axis.text = element_text(color = "black"),
      axis.title = element_text(color = "black"),
      axis.line = element_line(linewidth = 0.45, color = "black"),
      axis.ticks = element_line(linewidth = 0.45, color = "black"),
      plot.title = element_text(face = "bold", size = base_size + 1),
      plot.subtitle = element_text(size = base_size - 1, color = "grey25"),
      legend.title = element_text(size = base_size - 1),
      strip.background = element_blank()
    )
}
save_both <- function(plot, file, width, height) {
  tryCatch({
    pdf(file, width = width, height = height, useDingbats = FALSE); print(plot); dev.off()
    png(sub("\\.pdf$", ".png", file), width = width * 180, height = height * 180, res = 180)
    print(plot); dev.off()
    message("saved ", basename(file))
  }, error = function(e) message("FAIL save ", basename(file), ": ", conditionMessage(e)))
}
run <- function(label, fn) {
  tryCatch({ fn(); message("OK ", label) },
           error = function(e) message("FAIL ", label, ": ", conditionMessage(e)))
}
fmt_p <- function(p) ifelse(is.na(p), "NA", ifelse(p < 0.001, "P < 0.001", sprintf("P = %.3f", p)))

sc <- fread(file.path(proc, "sample_scores.tsv"))
kg <- fread(file.path(proc, "keygenes_samples.tsv"))
st <- fread(file.path(proc, "stage_gene_spearman.tsv"))
pca <- fread(file.path(proc, "pca_samples.tsv"))
near <- fread(file.path(proc, "nearest_stage.tsv"))
reg <- fread(file.path(proc, "regressions.tsv"))
de <- fread(file.path(proc, "deseq2_tall_vs_normal.tsv"))
heat <- fread(file.path(proc, "heatmap_stage_genes.tsv"))
tgt <- fread(file.path(proc, "target_maturity_panel.tsv"))
phf <- file.path(proc, "pharmacotype_maturity_panel.tsv")
ph <- if (file.exists(phf)) fread(phf) else data.table()
pvar <- fread(file.path(proc, "pca_variance.txt"), header = FALSE)
setnames(pvar, c("pc", "val"))
pc1 <- pvar[pc == "PC1", as.numeric(val)]
pc2 <- pvar[pc == "PC2", as.numeric(val)]

stg <- kg[group == "sorted_thymocyte"]
stg[, stage_label := factor(stage_label, levels = names(pal_stage))]
stage_lv <- names(pal_stage)
genes_tr <- intersect(c("ZEB1", "ZEB2", "LMO2", "IL7R", "CD34", "CD1A", "TAL1", "LYL1"), names(stg))
hm_genes <- intersect(c("CD34", "KIT", "SPI1", "MEF2C", "LMO2", "LYL1", "IL7R",
                        "CD1A", "RAG1", "DNTT", "CD4", "CD8A", "CD3E", "BCL11B",
                        "ZEB1", "ZEB2", "TAL1"), names(stg))
tall <- sc[is_tall == TRUE]

tryCatch({

# ---------------------------------------------------------------------
# 【189】 bump chart: gene rank across stages
# ---------------------------------------------------------------------
long <- melt(stg, id.vars = c("stage_label", "stage_index"), measure.vars = genes_tr,
             variable.name = "gene", value.name = "expr")
long[, rk := frank(-expr, ties.method = "first"), by = stage_index]
if (has_bump) {
  p_bump <- ggplot(long, aes(x = stage_index, y = rk, color = gene)) +
    ggbump::geom_bump(linewidth = 1.15, smooth = 8) +
    geom_point(size = 3.4) +
    scale_y_reverse(breaks = seq_along(genes_tr)) +
    scale_x_continuous(breaks = 1:6, labels = stage_lv) +
    scale_color_manual(values = pal_gene) +
    labs(title = "Marker rank along thymocyte stages (code2 189 bump)",
         subtitle = "Rank 1 = highest expression at that stage; n = 1 / stage",
         x = NULL, y = "Rank (1 = highest)", color = NULL) +
    theme_sci() + theme(axis.text.x = element_text(angle = 25, hjust = 1), legend.position = "right")
} else {
  p_bump <- ggplot(long, aes(x = stage_index, y = rk, color = gene, group = gene)) +
    geom_line(linewidth = 1.1) + geom_point(size = 3.2) +
    scale_y_reverse(breaks = seq_along(genes_tr)) +
    scale_x_continuous(breaks = 1:6, labels = stage_lv) +
    scale_color_manual(values = pal_gene) +
    labs(title = "Marker rank along thymocyte stages (code2 189)",
         x = NULL, y = "Rank (1 = highest)", color = NULL) +
    theme_sci() + theme(axis.text.x = element_text(angle = 25, hjust = 1))
}
save_both(p_bump, file.path(fig, "M2_bump_marker_ranks.pdf"), 8.2, 5.4)

# scaled area / ribbon trends (code2 change-area idea)
long[, scaled := (expr - min(expr)) / (max(expr) - min(expr) + 1e-8), by = gene]
p_area <- ggplot(long, aes(x = stage_index, y = scaled, fill = gene, color = gene)) +
  geom_area(alpha = 0.18, position = "identity") +
  geom_line(linewidth = 1.05) +
  geom_point(size = 2.2) +
  scale_x_continuous(breaks = 1:6, labels = stage_lv) +
  scale_color_manual(values = pal_gene) + scale_fill_manual(values = pal_gene) +
  facet_wrap(~gene, ncol = 4) +
  labs(title = "Min-max scaled developmental trajectories",
       subtitle = "Each gene scaled 0-1 within the 6 stages",
       x = NULL, y = "Scaled expression") +
  theme_sci() + theme(legend.position = "none", axis.text.x = element_text(angle = 30, hjust = 1))
save_both(p_area, file.path(fig, "M2_scaled_area_trajectories.pdf"), 10.2, 5.6)

# ---------------------------------------------------------------------
# 【66】 heatmap tiles of markers + ZEB1 line
# ---------------------------------------------------------------------
hm_genes <- intersect(c("CD34", "KIT", "SPI1", "MEF2C", "LMO2", "LYL1", "IL7R",
                        "CD1A", "RAG1", "DNTT", "CD4", "CD8A", "CD3E", "BCL11B",
                        "ZEB1", "ZEB2", "TAL1"), names(stg))
hm_long <- melt(stg, id.vars = "stage_label", measure.vars = hm_genes,
                variable.name = "gene", value.name = "expr")
hm_long[, z := as.numeric(scale(expr)), by = gene]
hm_long[, gene := factor(gene, levels = rev(hm_genes))]
p_tile <- ggplot(hm_long, aes(x = stage_label, y = gene, fill = z)) +
  geom_tile(color = "white", linewidth = 0.35) +
  scale_fill_gradient2(low = "#3C5488", mid = "white", high = "#E64B35", midpoint = 0) +
  labs(title = "Stage marker heatmap (code2 66/210)", x = NULL, y = NULL, fill = "row z") +
  theme_sci() + theme(axis.text.x = element_text(angle = 30, hjust = 1), axis.line = element_blank(),
                      axis.ticks = element_blank())
p_zeb_line <- ggplot(stg, aes(x = stage_label, y = ZEB1, group = 1)) +
  geom_line(color = "#E64B35", linewidth = 1.1) +
  geom_point(size = 2.8, fill = "#E64B35", shape = 21, color = "black") +
  labs(x = NULL, y = "ZEB1") + theme_sci() +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
if (has_patch) {
  save_both(p_zeb_line / p_tile + patchwork::plot_layout(heights = c(0.28, 1)),
            file.path(fig, "M2_heatmap_line_combo.pdf"), 7.6, 8.0)
} else {
  save_both(p_tile, file.path(fig, "M2_heatmap_line_combo.pdf"), 7.6, 6.5)
}

# ---------------------------------------------------------------------
# 【40】 PCA + density + ellipse
# ---------------------------------------------------------------------
pca2 <- copy(pca)
pca2[, grp := fifelse(group == "adult_T-ALL", "Adult T-ALL", "Normal thymus")]
p_pca <- ggplot(pca2, aes(PC1, PC2, color = grp, fill = grp)) +
  stat_ellipse(geom = "polygon", alpha = 0.12, level = 0.8, color = NA) +
  geom_point(aes(shape = grp), size = 2.6, stroke = 0.35) +
  scale_color_manual(values = c("Adult T-ALL" = "#E64B35", "Normal thymus" = "#4DBBD5")) +
  scale_fill_manual(values = c("Adult T-ALL" = "#E64B35", "Normal thymus" = "#4DBBD5")) +
  scale_shape_manual(values = c("Adult T-ALL" = 17, "Normal thymus" = 16)) +
  labs(title = "PCA with 80% ellipses (code2 40/95)",
       x = sprintf("PC1 (%.1f%%)", 100 * pc1), y = sprintf("PC2 (%.1f%%)", 100 * pc2),
       color = NULL, fill = NULL, shape = NULL) + theme_sci() + theme(legend.position = "top")
d1 <- ggplot(pca2, aes(PC1, fill = grp)) +
  geom_density(alpha = 0.55, color = "black", linewidth = 0.25) +
  scale_fill_manual(values = c("Adult T-ALL" = "#E64B35", "Normal thymus" = "#4DBBD5")) +
  theme_void() + theme(legend.position = "none")
d2 <- ggplot(pca2, aes(PC2, fill = grp)) +
  geom_density(alpha = 0.55, color = "black", linewidth = 0.25) +
  scale_fill_manual(values = c("Adult T-ALL" = "#E64B35", "Normal thymus" = "#4DBBD5")) +
  coord_flip() + theme_void() + theme(legend.position = "none")
if (has_aplot) {
  g <- aplot::insert_top(p_pca, d1, height = 0.28)
  g <- aplot::insert_right(g, d2, width = 0.28)
  pdf(file.path(fig, "M2_PCA_density.pdf"), width = 7.6, height = 6.4, useDingbats = FALSE)
  print(g); dev.off()
  png(file.path(fig, "M2_PCA_density.png"), width = 7.6 * 180, height = 6.4 * 180, res = 180)
  print(g); dev.off()
} else if (has_patch) {
  save_both(d1 + patchwork::plot_spacer() + p_pca + d2 + patchwork::plot_layout(ncol = 2, widths = c(4, 1), heights = c(1, 4)),
            file.path(fig, "M2_PCA_density.pdf"), 7.8, 6.6)
} else save_both(p_pca, file.path(fig, "M2_PCA_density.pdf"), 7.2, 5.8)

# ---------------------------------------------------------------------
# 【53】 T-ALL gene-gene correlation tile
# ---------------------------------------------------------------------
tall <- sc[is_tall == TRUE]
guse <- intersect(c("ZEB1", "ZEB2", "LMO2", "IL7R", "maturity_literature", "maturity_empirical"), names(tall))
cm <- cor(as.matrix(tall[, ..guse]), method = "spearman", use = "pairwise.complete.obs")
cdt <- as.data.table(as.table(cm))
setnames(cdt, c("a", "b", "r"))
cdt[, a := factor(a, levels = guse)]
cdt[, b := factor(b, levels = guse)]
p_cor <- ggplot(cdt, aes(a, b, fill = r)) +
  geom_tile(color = "white", linewidth = 0.4) +
  geom_point(aes(size = abs(r)), shape = 21, color = "grey20") +
  scale_fill_gradient2(low = "#3C5488", mid = "white", high = "#E64B35", midpoint = 0, limits = c(-1, 1)) +
  scale_size(range = c(1, 8), guide = "none") +
  labs(title = "Adult T-ALL Spearman correlations (code2 53)",
       subtitle = sprintf("n = %d", nrow(tall)), x = NULL, y = NULL, fill = "rho") +
  theme_sci() + theme(axis.text.x = element_text(angle = 35, hjust = 1), axis.line = element_blank(), axis.ticks = element_blank())
save_both(p_cor, file.path(fig, "M2_corr_tile_TALL.pdf"), 6.6, 5.8)

# ---------------------------------------------------------------------
# 【67】 dual axis: LMO2 bars + ZEB1 line
# ---------------------------------------------------------------------
ratio <- mean(stg$ZEB1) / mean(stg$LMO2)
p_dual <- ggplot(stg, aes(x = stage_label)) +
  geom_col(aes(y = LMO2), fill = "#F39B7F", alpha = 0.85, width = 0.62) +
  geom_line(aes(y = ZEB1 / ratio, group = 1), color = "#E64B35", linewidth = 1.15) +
  geom_point(aes(y = ZEB1 / ratio), fill = "#E64B35", shape = 21, size = 3, color = "black") +
  scale_y_continuous(name = "LMO2 (bars)",
                     sec.axis = sec_axis(~ . * ratio, name = "ZEB1 (line)")) +
  labs(title = "Dual-axis LMO2 vs ZEB1 (code2 67)", x = NULL) +
  theme_sci() +
  theme(axis.text.x = element_text(angle = 30, hjust = 1),
        axis.line.y.left = element_line(color = "#F39B7F"),
        axis.line.y.right = element_line(color = "#E64B35"),
        axis.text.y.left = element_text(color = "#B86A3A"),
        axis.text.y.right = element_text(color = "#E64B35"),
        axis.title.y.left = element_text(color = "#B86A3A"),
        axis.title.y.right = element_text(color = "#E64B35"))
save_both(p_dual, file.path(fig, "M2_dualaxis_LMO2_ZEB1.pdf"), 7.4, 5.2)

# ---------------------------------------------------------------------
# 【78】 alluvial: immunophenotype -> nearest stage
# ---------------------------------------------------------------------
if (has_alluvial && nrow(near)) {
  flow <- near[, .N, by = .(immunophenotype, nearest_label)]
  p_al <- ggplot(flow, aes(axis1 = immunophenotype, axis2 = nearest_label, y = N)) +
    ggalluvial::geom_alluvium(aes(fill = immunophenotype), width = 1 / 8, alpha = 0.8, knot.pos = 0.25) +
    ggalluvial::geom_stratum(width = 1 / 8, fill = "grey97", color = "grey30", linewidth = 0.3) +
    geom_text(stat = "stratum", aes(label = after_stat(stratum)), size = 3) +
    scale_fill_manual(values = pal_ip) +
    scale_x_discrete(limits = c("Immunophenotype", "Nearest stage"), expand = c(0.08, 0.08)) +
    labs(title = "Adult T-ALL mapped onto nearest thymocyte stage (code2 78)",
         y = "n patients") +
    theme_sci() + theme(legend.position = "none", axis.line.y = element_blank(), axis.ticks.y = element_blank())
  save_both(p_al, file.path(fig, "M2_alluvial_nearest_stage.pdf"), 8.4, 5.6)
}

# 【150】 radial bar of nearest stage
nn <- near[, .N, by = nearest_label]
nn[, nearest_label := factor(nearest_label, levels = stage_lv)]
p_rad <- ggplot(nn, aes(x = nearest_label, y = N, fill = nearest_label)) +
  geom_col(width = 0.9, color = "white") +
  coord_polar(theta = "y", start = 0) +
  geom_text(aes(y = 0, label = sprintf("%s  %d", nearest_label, N)), hjust = 1.05, size = 3.2) +
  scale_fill_manual(values = pal_stage) +
  labs(title = "Nearest-stage composition (code2 150 radial)") +
  theme_void() + theme(legend.position = "none", plot.title = element_text(face = "bold", hjust = 0.5))
save_both(p_rad, file.path(fig, "M2_radial_nearest_stage.pdf"), 6.8, 6.4)

# ---------------------------------------------------------------------
# 【191】 dumbbell: CD34+ vs mean(SP4,SP8)
# ---------------------------------------------------------------------
db <- st[gene %in% hm_genes, .(gene, cd34, sp = (sp4 + sp8) / 2)]
db <- db[is.finite(cd34) & is.finite(sp)]
setorder(db, cd34)
db[, gene := factor(gene, levels = gene)]
p_db <- ggplot(db) +
  geom_segment(aes(x = cd34, xend = sp, y = gene, yend = gene), color = "#C4B48A", linewidth = 1.35) +
  geom_point(aes(x = cd34, y = gene), color = "#E64B35", size = 3.3) +
  geom_point(aes(x = sp, y = gene), color = "#4DBBD5", size = 3.3) +
  labs(title = "CD34+ vs SP (mean SP4/SP8) (code2 191 dumbbell)",
       subtitle = "Red = CD34+; blue = mature SP", x = "log2 count", y = NULL) +
  theme_sci()
save_both(p_db, file.path(fig, "M2_dumbbell_CD34_vs_SP.pdf"), 7.0, 6.4)

# ---------------------------------------------------------------------
# 【346】 bidirectional: T-ALL mean minus stage mean
# ---------------------------------------------------------------------
stage_mean <- kg[group == "sorted_thymocyte", .(ZEB1_s = mean(ZEB1), ZEB2_s = mean(ZEB2), LMO2_s = mean(LMO2))]
delta <- data.table(
  gene = c("ZEB1", "ZEB2", "LMO2"),
  delta = c(mean(tall$ZEB1) - stage_mean$ZEB1_s,
            mean(tall$ZEB2) - stage_mean$ZEB2_s,
            mean(tall$LMO2) - stage_mean$LMO2_s)
)
delta[, side := fifelse(delta >= 0, "T-ALL higher", "T-ALL lower")]
p_bi <- ggplot(delta, aes(x = reorder(gene, delta), y = delta, fill = side)) +
  geom_hline(yintercept = 0, linetype = 2, color = "grey40") +
  geom_col(width = 0.62, color = "white") +
  geom_text(aes(label = sprintf("%+.2f", delta), hjust = fifelse(delta >= 0, -0.15, 1.15)), size = 3.3) +
  coord_flip() +
  scale_fill_manual(values = c("T-ALL higher" = "#E64B35", "T-ALL lower" = "#3C5488")) +
  labs(title = "Adult T-ALL minus mean sorted thymocyte (code2 346)",
       x = NULL, y = "Delta log2 count", fill = NULL) +
  theme_sci() + theme(legend.position = "top")
save_both(p_bi, file.path(fig, "M2_bidirectional_TALL_vs_stages.pdf"), 7.2, 4.6)

# ---------------------------------------------------------------------
# 【244】 polar radar of stage marker z-scores
# ---------------------------------------------------------------------
rad <- copy(hm_long)
rad_mean <- dcast(rad, stage_label ~ gene, value.var = "z")
rad_long <- melt(rad_mean, id.vars = "stage_label", variable.name = "gene", value.name = "z")
p_radar <- ggplot(rad_long, aes(x = gene, y = z, group = stage_label, color = stage_label, fill = stage_label)) +
  geom_polygon(alpha = 0.12, linewidth = 0.7) +
  geom_point(size = 1.6) +
  coord_polar() +
  scale_color_manual(values = pal_stage) + scale_fill_manual(values = pal_stage) +
  labs(title = "Stage marker radar (code2 244)", y = "row z", x = NULL, color = NULL, fill = NULL) +
  theme_sci() + theme(axis.line = element_blank(), axis.ticks = element_blank(),
                      axis.text.x = element_text(size = 8), legend.position = "right")
save_both(p_radar, file.path(fig, "M2_radar_stage_markers.pdf"), 8.0, 6.2)

# ---------------------------------------------------------------------
# 【243】 ternary ZEB1 / ZEB2 / LMO2
# ---------------------------------------------------------------------
run("ternary 243", function() {
  if (!has_tern) stop("ggtern not installed")
  td <- sc[group %in% c("sorted_thymocyte", "adult_T-ALL") & is.finite(ZEB1) & is.finite(ZEB2) & is.finite(LMO2)]
  td[, a := 2^ZEB1]; td[, b := 2^ZEB2]; td[, c := 2^LMO2]
  p_ter <- ggtern::ggtern(td, ggplot2::aes(x = a, y = b, z = c, color = group)) +
    ggplot2::geom_point(size = 2.4, alpha = 0.85) +
    ggplot2::scale_color_manual(values = pal_group) +
    ggplot2::labs(title = "ZEB1-ZEB2-LMO2 ternary (code2 243)", color = NULL) +
    ggplot2::theme(legend.position = "top")
  pdf(file.path(fig, "M2_ternary_ZEB1_ZEB2_LMO2.pdf"), width = 6.8, height = 6.2, useDingbats = FALSE)
  print(p_ter); dev.off()
  png(file.path(fig, "M2_ternary_ZEB1_ZEB2_LMO2.png"), width = 6.8 * 180, height = 6.2 * 180, res = 180)
  print(p_ter); dev.off()
  if ("package:ggtern" %in% search()) detach("package:ggtern", unload = FALSE)
})

# ---------------------------------------------------------------------
# 【206】 volcano: author DESeq2 T-ALL vs normal
# ---------------------------------------------------------------------
if (nrow(de) && "log2FC_tall_vs_normal" %in% names(de)) {
  dv <- copy(de)
  if ("gene_u" %in% names(dv)) dv[, gname := gene_u] else dv[, gname := as.character(gene)]
  dv[, log2FC := as.numeric(log2FC_tall_vs_normal)]
  dv[, fdr := as.numeric(fdr_tall_vs_normal)]
  dv <- dv[is.finite(log2FC) & is.finite(fdr) & fdr > 0]
  dv[, neglog := -log10(fdr)]
  dv[, cls := fcase(fdr < 0.05 & log2FC > 1, "up in T-ALL",
                    fdr < 0.05 & log2FC < -1, "down in T-ALL",
                    default = "NS")]
  highlight <- dv[gname %in% c("ZEB1", "ZEB2", "LMO2", "IL7R", "CD34", "CD1A", "TAL1", "LYL1")]
  p_vol <- ggplot(dv, aes(log2FC, neglog, color = cls)) +
    geom_point(size = 0.35, alpha = 0.45) +
    geom_point(data = highlight, size = 2.4, shape = 21, fill = "white", stroke = 0.9) +
    geom_vline(xintercept = c(-1, 1), linetype = 2, linewidth = 0.35, color = "grey50") +
    geom_hline(yintercept = -log10(0.05), linetype = 2, linewidth = 0.35, color = "grey50") +
    scale_color_manual(values = c("up in T-ALL" = "#E64B35", "down in T-ALL" = "#3C5488", NS = "grey75")) +
    labs(title = "T-ALL vs normal thymus (authors' DESeq2, code2 206)",
         subtitle = "Labels: key genes. Sorted stages n=1 so this is disease vs bulk normal, not stage DEG.",
         x = "log2 fold-change (T-ALL vs normal)", y = "-log10 FDR", color = NULL) +
    theme_sci() + theme(legend.position = "top")
  if (has_ggrepel && nrow(highlight)) {
    p_vol <- p_vol + ggrepel::geom_text_repel(data = highlight, aes(label = gname),
                                              size = 3, color = "black", max.overlaps = 20)
  }
  save_both(p_vol, file.path(fig, "M2_volcano_TALL_vs_normal.pdf"), 7.6, 6.0)
}

# ---------------------------------------------------------------------
# 【193】 scatter + equation
# ---------------------------------------------------------------------
ct12 <- suppressWarnings(cor.test(tall$maturity_literature, tall$ZEB1, method = "spearman", exact = FALSE))
fit <- lm(ZEB1 ~ maturity_literature, data = tall)
eq <- sprintf("ZEB1 = %.2f %+.2f * maturity\nSpearman rho = %.2f  %s",
              coef(fit)[1], coef(fit)[2], unname(ct12$estimate), fmt_p(ct12$p.value))
p_sc <- ggplot(tall, aes(maturity_literature, ZEB1, fill = immunophenotype)) +
  geom_smooth(method = "lm", formula = y ~ x, color = "grey25", fill = "#F3E5C8", linewidth = 0.6,
              inherit.aes = FALSE, aes(x = maturity_literature, y = ZEB1)) +
  geom_point(size = 2.6, shape = 21, color = "black", stroke = 0.3) +
  annotate("text", x = min(tall$maturity_literature, na.rm = TRUE),
           y = max(tall$ZEB1, na.rm = TRUE), label = eq, hjust = 0, vjust = 1, size = 3.2) +
  scale_fill_manual(values = pal_ip) +
  labs(title = "ZEB1 vs maturity with fitted line (code2 193)",
       x = "Literature maturity score", y = "log2 count ZEB1", fill = NULL) +
  theme_sci() + theme(legend.position = "top")
save_both(p_sc, file.path(fig, "M2_scatter_equation_ZEB1_maturity.pdf"), 7.4, 5.6)

# ---------------------------------------------------------------------
# 【88】 regression forest (restyle)
# ---------------------------------------------------------------------
if (nrow(reg)) {
  use <- reg[term != "intercept"]
  use[, lab := paste0(model, " | ", term)]
  setorder(use, coef)
  use[, lab := factor(lab, levels = unique(lab))]
  p_f <- ggplot(use, aes(x = coef, y = lab)) +
    geom_vline(xintercept = 0, linetype = 2, color = "grey40") +
    geom_errorbar(aes(xmin = ci95_lo, xmax = ci95_hi), width = 0.2, linewidth = 0.5) +
    geom_point(aes(fill = coef < 0), shape = 22, size = 3.4, color = "black") +
    scale_fill_manual(values = c("TRUE" = "#3C5488", "FALSE" = "#E64B35"), guide = "none") +
    labs(title = "OLS coefficients, 95% CI (code2 88 forest)",
         subtitle = "After maturity, the T-ALL indicator for ZEB1 is no longer significant",
         x = "Coefficient", y = NULL) +
    theme_sci()
  save_both(p_f, file.path(fig, "M2_forest_OLS.pdf"), 9.6, max(4.8, 0.36 * nrow(use) + 1.8))
}

# ---------------------------------------------------------------------
# 2.18 TARGET + Pharmacotype
# ---------------------------------------------------------------------
if (nrow(tgt) && "maturity_literature" %in% names(tgt)) {
  tgt[is.na(subtype) | subtype == "", subtype := "Unknown"]
  ct_t <- suppressWarnings(cor.test(tgt$ZEB1_log, tgt$maturity_literature, method = "spearman", exact = FALSE))
  p18 <- ggplot(tgt, aes(maturity_literature, ZEB1_log, fill = subtype)) +
    geom_smooth(method = "lm", formula = y ~ x, color = "grey20", fill = "#F3E5C8", linewidth = 0.55,
                inherit.aes = FALSE, aes(x = maturity_literature, y = ZEB1_log)) +
    geom_point(size = 2.0, shape = 21, color = "black", stroke = 0.25, alpha = 0.85) +
    scale_fill_manual(values = pal_sub, drop = FALSE) +
    labs(title = "TARGET replication: ZEB1 vs maturity",
         subtitle = sprintf("Pediatric RNA-seq; n = %d; rho = %.2f  %s",
                            nrow(tgt), unname(ct_t$estimate), fmt_p(ct_t$p.value)),
         x = "Maturity score", y = "log2(TPM+1) ZEB1", fill = "Subtype") +
    theme_sci()
  save_both(p18, file.path(fig, "M2_2.18_TARGET_ZEB1_maturity.pdf"), 7.6, 5.5)
  # subtype rain is OK as ONE panel, plus bump-like box of maturity by subtype
  p_mat <- ggplot(tgt[subtype != "Unknown"], aes(x = subtype, y = maturity_literature, fill = subtype)) +
    geom_boxplot(width = 0.55, outlier.shape = NA, alpha = 0.9) +
    geom_jitter(width = 0.12, size = 0.7, alpha = 0.45, color = "grey20") +
    scale_fill_manual(values = pal_sub) +
    labs(title = "TARGET maturity score by Liu subtype",
         subtitle = "Literature thymocyte signature applied to pediatric T-ALL",
         x = NULL, y = "Maturity score") +
    theme_sci() + theme(legend.position = "none", axis.text.x = element_text(angle = 25, hjust = 1))
  save_both(p_mat, file.path(fig, "M2_2.18_TARGET_maturity_by_subtype.pdf"), 7.2, 5.2)
}
if (nrow(ph) && "maturity_literature" %in% names(ph) && "ZEB1_log" %in% names(ph) && nrow(ph) >= 10) {
  if (!"etp" %in% names(ph)) ph[, etp := "notETP"]
  ph[, etp_lab := fifelse(etp == "ETP", "ETP", "notETP")]
  ct_p <- suppressWarnings(cor.test(ph$ZEB1_log, ph$maturity_literature, method = "spearman", exact = FALSE))
  p18b <- ggplot(ph, aes(maturity_literature, ZEB1_log, fill = etp_lab)) +
    geom_smooth(method = "lm", formula = y ~ x, color = "grey20", fill = "#F3E5C8", linewidth = 0.55,
                inherit.aes = FALSE, aes(x = maturity_literature, y = ZEB1_log)) +
    geom_point(size = 2.4, shape = 21, color = "black", stroke = 0.3) +
    scale_fill_manual(values = c(ETP = "#E64B35", notETP = "#4DBBD5")) +
    labs(title = "Pharmacotype replication: ZEB1 vs maturity",
         subtitle = sprintf("n = %d; rho = %.2f  %s", nrow(ph), unname(ct_p$estimate), fmt_p(ct_p$p.value)),
         x = "Maturity score", y = "log2(FPKM+1) ZEB1", fill = NULL) +
    theme_sci()
  save_both(p18b, file.path(fig, "M2_2.18_Pharmacotype_ZEB1_maturity.pdf"), 7.4, 5.5)
  fwrite(data.table(
    cohort = c("TARGET", "Pharmacotype"),
    rho = c(if (exists("ct_t")) unname(ct_t$estimate) else NA_real_, unname(ct_p$estimate)),
    p = c(if (exists("ct_t")) ct_t$p.value else NA_real_, ct_p$p.value),
    n = c(nrow(tgt), nrow(ph))
  ), file.path(tab, "M2_2.18_replication.tsv"), sep = "\t")
}

}, error = function(e) message("FAIL main enrich block: ", conditionMessage(e)))

# ---------------------------------------------------------------------
# Extra code2 recipes that are NOT violin / lollipop
# 【118】 nested donut / sunburst sample map
# 【141】 arrow (CD34+ -> SP change)
# 【238】 bar + beeswarm (immunophenotype)
# 【136】 stacked area of nearest-stage composition by immunophenotype
# 【34】  circular heatmap of stage markers
# 【142】 signed stage Spearman of key genes
# 【235】 ZEB1 vs LMO2 interaction scatter
# 【60】  waffle of nearest stage
# ---------------------------------------------------------------------
run("sunburst 118", function() {
  inner <- sc[, .N, by = group]
  inner[, frac := N / sum(N)]
  inner[, ymax := cumsum(frac)]
  inner[, ymin := c(0, head(ymax, -1))]
  outer <- sc[, .N, by = .(group, immunophenotype)]
  outer[, frac := N / sum(N)]
  outer[, ymax := cumsum(frac)]
  outer[, ymin := c(0, head(ymax, -1))]
  p <- ggplot() +
    geom_rect(data = outer, aes(xmin = 3.2, xmax = 4.2, ymin = ymin, ymax = ymax, fill = immunophenotype),
              color = "white", linewidth = 0.35) +
    geom_rect(data = inner, aes(xmin = 2.0, xmax = 3.15, ymin = ymin, ymax = ymax, fill = group),
              color = "white", linewidth = 0.35) +
    coord_polar(theta = "y") +
    xlim(0.6, 4.3) +
    scale_fill_manual(values = c(pal_group, pal_ip), drop = FALSE) +
    labs(title = "Sample composition (code2 118 sunburst)", fill = NULL) +
    theme_void() + theme(legend.position = "right", plot.title = element_text(face = "bold", hjust = 0.5))
  save_both(p, file.path(fig, "M2_sunburst_samples.pdf"), 7.2, 6.2)
})

run("arrow 141", function() {
  ar <- st[gene %in% hm_genes]
  ar[, dlt := sp_mean - cd34]
  ar[, dir := fifelse(dlt >= 0, "up toward SP", "down toward SP")]
  setorder(ar, dlt)
  ar[, gene := factor(gene, levels = gene)]
  p <- ggplot(ar, aes(x = gene, y = dlt, color = dir)) +
    geom_hline(yintercept = 0, linewidth = 0.3, color = "grey40") +
    geom_segment(aes(xend = gene, y = 0, yend = dlt),
                 arrow = arrow(type = "closed", length = unit(1.6, "mm")), linewidth = 0.9) +
    scale_color_manual(values = c("up toward SP" = "#00A087", "down toward SP" = "#E64B35")) +
    labs(title = "CD34+ to SP expression change (code2 141 arrows)",
         subtitle = "Positive = higher in SP4/SP8 than CD34+",
         x = NULL, y = "Delta log2 (SP mean - CD34+)", color = NULL) +
    theme_sci() + theme(axis.text.x = element_text(angle = 40, hjust = 1), legend.position = "top")
  save_both(p, file.path(fig, "M2_arrow_CD34_to_SP.pdf"), 8.4, 5.2)
})

run("bar beeswarm 238", function() {
  d <- copy(tall)
  d[, immunophenotype := factor(immunophenotype, levels = names(pal_ip))]
  p <- ggplot(d, aes(x = immunophenotype, y = ZEB1, fill = immunophenotype)) +
    stat_summary(fun = mean, geom = "col", width = 0.62, color = "black", linewidth = 0.35, alpha = 0.85) +
    stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.18, linewidth = 0.45) +
    {
      if (has_beeswarm) ggbeeswarm::geom_beeswarm(shape = 21, color = "black", size = 2.2, stroke = 0.3, cex = 2.2)
      else geom_jitter(shape = 21, color = "black", size = 2.2, width = 0.12, height = 0)
    } +
    scale_fill_manual(values = pal_ip) +
    labs(title = "Adult T-ALL ZEB1 by immunophenotype (code2 238 bar+beeswarm)",
         subtitle = sprintf("n = %d; not a violin", nrow(d)),
         x = NULL, y = "log2 count ZEB1") +
    theme_sci() + theme(legend.position = "none", axis.text.x = element_text(angle = 20, hjust = 1))
  save_both(p, file.path(fig, "M2_bar_beeswarm_immunophenotype.pdf"), 6.8, 5.4)
})

run("stacked nearest 136", function() {
  flow <- near[, .N, by = .(immunophenotype, nearest_label)]
  flow[, nearest_label := factor(nearest_label, levels = stage_lv)]
  p <- ggplot(flow, aes(x = immunophenotype, y = N, fill = nearest_label)) +
    geom_col(width = 0.72, color = "white") +
    scale_fill_manual(values = pal_stage, drop = FALSE) +
    labs(title = "Nearest thymocyte stage by immunophenotype (code2 136 stacked)",
         x = NULL, y = "n patients", fill = "Nearest stage") +
    theme_sci() + theme(axis.text.x = element_text(angle = 20, hjust = 1))
  save_both(p, file.path(fig, "M2_stacked_nearest_by_IP.pdf"), 7.4, 5.4)
})

run("circular heatmap 34", function() {
  genes_c <- intersect(hm_genes, names(stg))
  mat <- as.matrix(stg[, ..genes_c])
  mat <- scale(mat)
  rownames(mat) <- as.character(stg$stage_label)
  md <- melt(as.data.table(mat, keep.rownames = "stage"), id.vars = "stage",
             variable.name = "gene", value.name = "z")
  md[, stage := factor(stage, levels = stage_lv)]
  p <- ggplot(md, aes(x = gene, y = stage, fill = z)) +
    geom_tile(color = "white", linewidth = 0.25) +
    coord_polar() +
    scale_fill_gradient2(low = "#3C5488", mid = "white", high = "#E64B35", midpoint = 0) +
    labs(title = "Circular stage-marker heatmap (code2 34)", fill = "z") +
    theme_sci() + theme(axis.line = element_blank(), axis.ticks = element_blank(),
                        axis.text.x = element_text(size = 7), axis.title = element_blank())
  save_both(p, file.path(fig, "M2_circular_heatmap_markers.pdf"), 8.0, 7.2)
})

run("signed spearman 142", function() {
  ks <- st[gene %in% c("ZEB1", "ZEB2", "LMO2", "IL7R", "CD34", "CD1A", "TAL1", "LYL1",
                       "KIT", "RAG1", "DNTT", "BCL11B")]
  setorder(ks, rho_stage)
  ks[, gene := factor(gene, levels = gene)]
  ks[, side := fifelse(rho_stage >= 0, "increases with stage", "decreases with stage")]
  p <- ggplot(ks, aes(x = rho_stage, y = gene, fill = side)) +
    geom_vline(xintercept = 0, linetype = 2, color = "grey40") +
    geom_col(width = 0.7, color = "white") +
    geom_text(aes(label = sprintf("%.2f", rho_stage),
                  hjust = fifelse(rho_stage >= 0, -0.15, 1.15)), size = 3) +
    scale_fill_manual(values = c("increases with stage" = "#E64B35", "decreases with stage" = "#3C5488")) +
    labs(title = "Stage Spearman of markers (code2 142 signed bars)",
         subtitle = "n = 1 per stage; descriptive only",
         x = "Spearman rho vs stage order", y = NULL, fill = NULL) +
    theme_sci() + theme(legend.position = "top")
  save_both(p, file.path(fig, "M2_signed_stage_spearman.pdf"), 7.6, 6.0)
})

run("interaction scatter 235", function() {
  d <- sc[group %in% c("sorted_thymocyte", "adult_T-ALL")]
  p <- ggplot(d, aes(LMO2, ZEB1, fill = group, shape = group)) +
    geom_smooth(method = "lm", formula = y ~ x, se = TRUE, linewidth = 0.55,
                aes(color = group, fill = group), alpha = 0.12) +
    geom_point(size = 2.8, color = "black", stroke = 0.3) +
    scale_fill_manual(values = pal_group) +
    scale_color_manual(values = pal_group) +
    scale_shape_manual(values = c(sorted_thymocyte = 21, "adult_T-ALL" = 24)) +
    labs(title = "ZEB1 vs LMO2 in stages vs adult T-ALL (code2 235)",
         x = "log2 count LMO2", y = "log2 count ZEB1", fill = NULL, color = NULL, shape = NULL) +
    theme_sci() + theme(legend.position = "top")
  save_both(p, file.path(fig, "M2_scatter_ZEB1_vs_LMO2.pdf"), 7.2, 5.6)
})

run("waffle nearest 60", function() {
  nn <- near[, .N, by = nearest_label]
  nn[, nearest_label := factor(nearest_label, levels = stage_lv)]
  cells <- nn[, .(id = seq_len(N)), by = nearest_label]
  cells[, i := .I]
  cells[, x := ((i - 1) %% 10) + 1]
  cells[, y := floor((i - 1) / 10) + 1]
  p <- ggplot(cells, aes(x, y, fill = nearest_label)) +
    geom_tile(color = "white", linewidth = 0.6, width = 0.92, height = 0.92) +
    coord_equal() +
    scale_fill_manual(values = pal_stage, drop = FALSE) +
    labs(title = "Nearest-stage waffle (code2 60; 1 tile = 1 patient)",
         fill = NULL) +
    theme_void() + theme(legend.position = "right", plot.title = element_text(face = "bold"))
  save_both(p, file.path(fig, "M2_waffle_nearest_stage.pdf"), 7.0, 4.6)
})

run("complexheatmap 210", function() {
  if (!has_cht || !has_circlize) stop("ComplexHeatmap/circlize not installed")
  genes_c <- intersect(hm_genes, names(stg))
  mat <- t(scale(as.matrix(stg[, ..genes_c])))
  colnames(mat) <- as.character(stg$stage_label)
  pdf(file.path(fig, "M2_complexheatmap_stages.pdf"), width = 6.8, height = 7.2)
  ComplexHeatmap::Heatmap(
    mat, name = "row z",
    col = circlize::colorRamp2(c(-2, 0, 2), c("#3C5488", "white", "#E64B35")),
    cluster_columns = FALSE, cluster_rows = TRUE,
    column_title = "Stage marker heatmap (code2 210)",
    rect_gp = grid::gpar(col = "white", lwd = 0.4),
    heatmap_legend_param = list(direction = "horizontal")
  )
  dev.off()
  png(file.path(fig, "M2_complexheatmap_stages.png"), width = 6.8 * 180, height = 7.2 * 180, res = 180)
  ComplexHeatmap::Heatmap(
    mat, name = "row z",
    col = circlize::colorRamp2(c(-2, 0, 2), c("#3C5488", "white", "#E64B35")),
    cluster_columns = FALSE, cluster_rows = TRUE,
    column_title = "Stage marker heatmap (code2 210)",
    rect_gp = grid::gpar(col = "white", lwd = 0.4)
  )
  dev.off()
})

run("LMO2 ZEB2 vs maturity 193 extra", function() {
  ct_l <- suppressWarnings(cor.test(tall$maturity_literature, tall$LMO2, method = "spearman", exact = FALSE))
  ct_z <- suppressWarnings(cor.test(tall$maturity_literature, tall$ZEB2, method = "spearman", exact = FALSE))
  long <- melt(tall, id.vars = c("sample", "maturity_literature", "immunophenotype"),
               measure.vars = c("LMO2", "ZEB2"), variable.name = "gene", value.name = "expr")
  p <- ggplot(long, aes(maturity_literature, expr, fill = immunophenotype)) +
    geom_smooth(method = "lm", formula = y ~ x, color = "grey25", fill = "#F3E5C8",
                linewidth = 0.55, inherit.aes = FALSE, aes(x = maturity_literature, y = expr)) +
    geom_point(size = 2.3, shape = 21, color = "black", stroke = 0.3) +
    facet_wrap(~gene, scales = "free_y") +
    scale_fill_manual(values = pal_ip) +
    labs(title = "LMO2 and ZEB2 vs maturity (code2 193)",
         subtitle = sprintf("LMO2 rho = %.2f %s; ZEB2 rho = %.2f %s",
                            unname(ct_l$estimate), fmt_p(ct_l$p.value),
                            unname(ct_z$estimate), fmt_p(ct_z$p.value)),
         x = "Literature maturity score", y = "log2 count", fill = NULL) +
    theme_sci() + theme(legend.position = "top")
  save_both(p, file.path(fig, "M2_scatter_LMO2_ZEB2_maturity.pdf"), 9.2, 5.2)
})

cat("MODULE2 ENRICH DONE\n")
cat("FIG", length(list.files(fig, pattern = "\\.(pdf|png)$")), "\n")


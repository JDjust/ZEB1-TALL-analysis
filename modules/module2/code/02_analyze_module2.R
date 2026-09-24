# Module 2: ZEB1 vs normal T-cell development (GSE142521/GSE142522)
# code2 recipes: raincloud, lollipop, scatter+lm, heatmap, PCA.
# Sorted stages n=1: trajectory plots are descriptive. T-ALL n=41 is the statistical unit.

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})
has_beeswarm <- requireNamespace("ggbeeswarm", quietly = TRUE)
has_cht <- requireNamespace("ComplexHeatmap", quietly = TRUE)
has_circlize <- requireNamespace("circlize", quietly = TRUE)
has_ggrepel <- requireNamespace("ggrepel", quietly = TRUE)
if (has_beeswarm) library(ggbeeswarm)

root <- "/data-b/liangfuhua/projects/ZEB1_TALL_analysis/module2"
fig <- file.path(root, "figures")
tab <- file.path(root, "tables")
proc <- file.path(root, "processed")
dir.create(fig, FALSE, TRUE)
dir.create(tab, FALSE, TRUE)

pal_stage <- c(
  "CD34+" = "#E64B35", ISP = "#F39B7F", "Early cortical" = "#00A087",
  "Late cortical" = "#3C5488", SP4 = "#4DBBD5", SP8 = "#8491B4",
  "total_thymus" = "#91D1C2", "T-ALL" = "#7E6148"
)
pal_group <- c(sorted_thymocyte = "#4DBBD5", total_thymus = "#91D1C2", "adult_T-ALL" = "#E64B35")
pal_ip <- c(Normal = "#91D1C2", Immature = "#E64B35", Pre_ab_Cortical = "#00A087", TCR_Mature = "#3C5488")
pal_onco <- c(Normal = "#91D1C2", TAL1 = "#E64B35", TLX1 = "#4DBBD5", TLX3 = "#00A087", unknown = "grey70")
pal_gene <- c(ZEB1 = "#E64B35", ZEB2 = "#3C5488", LMO2 = "#F39B7F", IL7R = "#00A087")
pal_sub <- c(
  TAL = "#E64B35", TLX = "#4DBBD5", HOXA = "#00A087",
  "LMO2/LYL1" = "#3C5488", NKX2_1 = "#F39B7F", Unknown = "#8491B4"
)

theme_sci <- function(base_size = 11) {
  theme_classic(base_size = base_size) +
    theme(
      axis.text = element_text(color = "black"),
      axis.title = element_text(color = "black"),
      axis.line = element_line(linewidth = 0.45, color = "black"),
      axis.ticks = element_line(linewidth = 0.45, color = "black"),
      plot.title = element_text(face = "bold", size = base_size + 1, color = "black"),
      plot.subtitle = element_text(size = base_size - 1, color = "grey25"),
      legend.title = element_text(size = base_size - 1),
      strip.background = element_blank(),
      strip.text = element_text(face = "bold")
    )
}
save_both <- function(plot, file, width, height) {
  pdf(file, width = width, height = height, useDingbats = FALSE); print(plot); dev.off()
  png(sub("\\.pdf$", ".png", file), width = width * 180, height = height * 180, res = 180)
  print(plot); dev.off()
}
add_points <- function(size = 0.5, alpha = 0.5, width = 0.16) {
  if (has_beeswarm) geom_quasirandom(size = size, alpha = alpha, width = width, color = "grey20")
  else geom_jitter(size = size, alpha = alpha, width = width, color = "grey20")
}
fmt_p <- function(p) ifelse(is.na(p), "NA", ifelse(p < 0.001, "P < 0.001", sprintf("P = %.3f", p)))
p_star <- function(p) ifelse(is.na(p), "ns", ifelse(p < 0.001, "***", ifelse(p < 0.01, "**", ifelse(p < 0.05, "*", "ns"))))

plot_rain <- function(df, x, y, fill_var, pal, ylab, title, subtitle = NULL, angle = 35) {
  ggplot(df, aes(x = .data[[x]], y = .data[[y]], fill = .data[[fill_var]])) +
    geom_violin(trim = FALSE, scale = "width", color = NA, alpha = 0.72, width = 0.88) +
    add_points() +
    geom_boxplot(width = 0.16, fill = "white", outlier.shape = NA, linewidth = 0.4) +
    scale_fill_manual(values = pal, drop = FALSE) +
    labs(x = NULL, y = ylab, title = title, subtitle = subtitle) +
    theme_sci() +
    theme(legend.position = "none",
          axis.text.x = element_text(angle = angle, hjust = if (angle > 0) 1 else 0.5))
}

file.copy(file.path(proc, "M2_data_note.txt"), file.path(tab, "M2_data_note.txt"), overwrite = TRUE)

sc <- fread(file.path(proc, "sample_scores.tsv"))
kg <- fread(file.path(proc, "keygenes_samples.tsv"))
meta <- fread(file.path(proc, "sample_metadata.tsv"))
st <- fread(file.path(proc, "stage_gene_spearman.tsv"))
pca <- fread(file.path(proc, "pca_samples.tsv"))
near <- fread(file.path(proc, "nearest_stage.tsv"))
reg <- fread(file.path(proc, "regressions.tsv"))
ct <- fread(file.path(proc, "zeb1_contrasts.tsv"))
heat <- fread(file.path(proc, "heatmap_stage_genes.tsv"))
tgt <- fread(file.path(proc, "target_maturity_panel.tsv"))
ph <- fread(file.path(proc, "pharmacotype_maturity_panel.tsv"))
pvar <- fread(file.path(proc, "pca_variance.txt"), header = FALSE)
setnames(pvar, c("pc", "val"))

fwrite(meta[, .N, by = .(group, immunophenotype, oncogene)], file.path(tab, "M2_2.1_sample_counts.tsv"), sep = "\t")
fwrite(st[gene %in% c("ZEB1", "ZEB2", "LMO2", "IL7R", "TAL1", "LYL1")], file.path(tab, "M2_stage_keygene_rho.tsv"), sep = "\t")
fwrite(ct, file.path(tab, "M2_ZEB1_contrasts.tsv"), sep = "\t")
fwrite(reg, file.path(tab, "M2_2.16_2.17_regressions.tsv"), sep = "\t")
fwrite(near[, .N, by = .(nearest_label, immunophenotype)], file.path(tab, "M2_2.11_nearest_stage.tsv"), sep = "\t")

# 2.1 sample map
sc[, group_lab := factor(group, levels = c("sorted_thymocyte", "total_thymus", "adult_T-ALL"),
                         labels = c("Sorted thymocytes", "Total thymus", "Adult T-ALL"))]
cnt <- sc[, .N, by = group_lab]
p1 <- ggplot(cnt, aes(x = group_lab, y = N, fill = group_lab)) +
  geom_col(width = 0.7, color = "white") +
  geom_text(aes(label = N), vjust = -0.25, size = 4) +
  scale_fill_manual(values = c("Sorted thymocytes" = "#4DBBD5", "Total thymus" = "#91D1C2", "Adult T-ALL" = "#E64B35")) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.18))) +
  labs(title = "GSE142522 samples used in Module 2",
       subtitle = "Cieslak et al. SOLiD RNA-seq; processed log2 counts (Table S1)",
       x = NULL, y = "n samples") +
  theme_sci() + theme(legend.position = "none")
save_both(p1, file.path(fig, "M2_2.1_sample_map.pdf"), 6.4, 5.0)

sc[group == "adult_T-ALL", immunophenotype := factor(immunophenotype, levels = c("Immature", "Pre_ab_Cortical", "TCR_Mature"))]
cnt2 <- sc[group == "adult_T-ALL", .N, by = immunophenotype]
p1b <- ggplot(cnt2, aes(x = immunophenotype, y = N, fill = immunophenotype)) +
  geom_col(width = 0.7, color = "white") +
  geom_text(aes(label = N), vjust = -0.25, size = 4) +
  scale_fill_manual(values = pal_ip) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.18))) +
  labs(title = "Adult T-ALL immunophenotype (GSE142520)", x = NULL, y = "n patients") +
  theme_sci() + theme(legend.position = "none", axis.text.x = element_text(angle = 20, hjust = 1))
save_both(p1b, file.path(fig, "M2_2.1_TALL_immunophenotype.pdf"), 6.0, 5.0)

# 2.2-2.5 trajectories
stg <- kg[group == "sorted_thymocyte"]
stg[, stage_label := factor(stage_label, levels = c("CD34+", "ISP", "Early cortical", "Late cortical", "SP4", "SP8"))]
plot_traj <- function(gname, ylab, title) {
  d <- stg[, .(stage_label, stage_index, val = get(gname))]
  rho <- st[gene == gname, rho_stage][1]
  pp <- st[gene == gname, p_stage][1]
  ggplot(d, aes(x = stage_label, y = val, group = 1)) +
    geom_line(color = pal_gene[[gname]], linewidth = 1.05) +
    geom_point(size = 3.2, fill = pal_gene[[gname]], color = "black", shape = 21, stroke = 0.35) +
    labs(title = title,
         subtitle = sprintf("n = 1 / stage; Spearman vs order rho = %.2f  %s", rho, fmt_p(pp)),
         x = NULL, y = ylab) +
    theme_sci() + theme(axis.text.x = element_text(angle = 30, hjust = 1))
}
save_both(plot_traj("ZEB1", "log2 count ZEB1", "2.2  Normal thymocytes: ZEB1 along T-cell development"),
          file.path(fig, "M2_2.2_ZEB1_trajectory.pdf"), 7.2, 5.0)
save_both(plot_traj("ZEB2", "log2 count ZEB2", "2.3  Normal thymocytes: ZEB2"),
          file.path(fig, "M2_2.3_ZEB2_trajectory.pdf"), 7.2, 5.0)
save_both(plot_traj("LMO2", "log2 count LMO2", "2.4  Normal thymocytes: LMO2"),
          file.path(fig, "M2_2.4_LMO2_trajectory.pdf"), 7.2, 5.0)
save_both(plot_traj("IL7R", "log2 count IL7R", "2.5  Normal thymocytes: IL7R"),
          file.path(fig, "M2_2.5_IL7R_trajectory.pdf"), 7.2, 5.0)

long4 <- melt(stg, id.vars = c("stage_label", "stage_index"),
              measure.vars = c("ZEB1", "ZEB2", "LMO2", "IL7R"),
              variable.name = "gene", value.name = "expr")
pcomb <- ggplot(long4, aes(x = stage_label, y = expr, color = gene, group = gene)) +
  geom_line(linewidth = 1.0) +
  geom_point(size = 2.6) +
  scale_color_manual(values = pal_gene) +
  labs(title = "ZEB1 / ZEB2 / LMO2 / IL7R along thymocyte stages",
       subtitle = "Same y-scale; n = 1 per stage", x = NULL, y = "log2 count", color = NULL) +
  theme_sci() + theme(axis.text.x = element_text(angle = 30, hjust = 1), legend.position = "top")
save_both(pcomb, file.path(fig, "M2_2.2_2.5_four_gene_trajectories.pdf"), 8.0, 5.2)

# 2.7 GAM/loess
p27 <- ggplot(stg, aes(x = stage_index, y = ZEB1)) +
  geom_smooth(method = "loess", span = 1, color = "#B07D3A", fill = "#F3E5C8", linewidth = 0.9) +
  geom_point(size = 3.2, fill = "#E64B35", color = "black", shape = 21) +
  geom_text(aes(label = stage_label), vjust = -0.9, size = 3) +
  scale_x_continuous(breaks = 1:6, labels = c("CD34+", "ISP", "EC", "LC", "SP4", "SP8")) +
  labs(title = "2.7  ZEB1 smoother across ordered stages",
       subtitle = "loess on n = 6; not a biological-replicate GAM",
       x = "Developmental order", y = "log2 count ZEB1") +
  theme_sci()
save_both(p27, file.path(fig, "M2_2.7_ZEB1_loess.pdf"), 7.2, 5.2)

# 2.6 / 2.8 heatmap
hm <- as.matrix(heat[, -1, with = FALSE])
rownames(hm) <- heat[[1]]
stage_cols <- intersect(c("CD34plus", "ISP", "EC", "LC", "SP4", "SP8"), colnames(hm))
hm_s <- hm[, stage_cols, drop = FALSE]
hm_z <- t(scale(t(hm_s)))
if (has_cht) {
  if (has_circlize) {
    col_fun <- circlize::colorRamp2(c(-2, 0, 2), c("#3C5488", "white", "#E64B35"))
  } else col_fun <- NULL
  pdf(file.path(fig, "M2_2.6_stage_marker_heatmap.pdf"), width = 6.4, height = 7.2)
  ComplexHeatmap::draw(ComplexHeatmap::Heatmap(
    hm_z, name = "row z", col = col_fun, cluster_columns = FALSE, cluster_rows = FALSE,
    column_labels = c("CD34+", "ISP", "EC", "LC", "SP4", "SP8"),
    column_names_rot = 45, row_names_gp = grid::gpar(fontsize = 9),
    heatmap_legend_param = list(title = "row z"),
    column_title = "2.6  Marker genes across thymocyte stages"
  ))
  dev.off()
  png(file.path(fig, "M2_2.6_stage_marker_heatmap.png"), width = 6.4 * 180, height = 7.2 * 180, res = 180)
  ComplexHeatmap::draw(ComplexHeatmap::Heatmap(
    hm_z, name = "row z", col = col_fun, cluster_columns = FALSE, cluster_rows = FALSE,
    column_labels = c("CD34+", "ISP", "EC", "LC", "SP4", "SP8"),
    column_names_rot = 45, row_names_gp = grid::gpar(fontsize = 9),
    column_title = "2.6  Marker genes across thymocyte stages"
  ))
  dev.off()
  top_up <- st[is.finite(rho_stage)][order(-rho_stage)][1:25]
  top_dn <- st[is.finite(rho_stage)][order(rho_stage)][1:25]
  topg <- unique(c(top_up$gene, top_dn$gene, "ZEB1", "ZEB2", "LMO2"))
  expr <- fread(file.path(proc, "expression_log2.tsv"))
  if ("gene_u" %in% names(expr)) {
    gkeep <- expr$gene_u %in% topg
    em <- as.matrix(expr[gkeep, ..stage_cols])
    rownames(em) <- expr$gene_u[gkeep]
  } else {
    em <- as.matrix(expr[expr[[1]] %in% topg, ..stage_cols])
    rownames(em) <- expr[[1]][expr[[1]] %in% topg]
  }
  emz <- t(scale(t(em)))
  pdf(file.path(fig, "M2_2.8_stage_correlated_heatmap.pdf"), width = 6.6, height = 9.0)
  ComplexHeatmap::draw(ComplexHeatmap::Heatmap(
    emz, name = "row z", col = col_fun, cluster_columns = FALSE,
    column_labels = c("CD34+", "ISP", "EC", "LC", "SP4", "SP8"),
    column_names_rot = 45, row_names_gp = grid::gpar(fontsize = 7.5),
    column_title = "2.8  Top stage-correlated genes (descriptive, n=1/stage)"
  ))
  dev.off()
}

stg_sc <- sc[group == "sorted_thymocyte"]
stg_sc[, stage_label := factor(stage_label, levels = c("CD34+", "ISP", "Early cortical", "Late cortical", "SP4", "SP8"))]
p29 <- ggplot(stg_sc, aes(x = stage_label, y = maturity_literature, group = 1)) +
  geom_line(color = "#3C5488", linewidth = 1.05) +
  geom_point(size = 3.2, fill = "#3C5488", color = "black", shape = 21) +
  labs(title = "2.9  Literature thymocyte maturation score",
       subtitle = "late (CD1A/CD3/CD4/CD8/RAG1/BCL11B/TCF7) minus early (CD34/KIT/IL7R/LMO2/LYL1/MEF2C)",
       x = NULL, y = "Maturity score") +
  theme_sci() + theme(axis.text.x = element_text(angle = 30, hjust = 1))
save_both(p29, file.path(fig, "M2_2.9_maturity_score_stages.pdf"), 7.2, 5.0)

# 2.10 PCA
pc1 <- pvar[pc == "PC1", as.numeric(val)]
pc2 <- pvar[pc == "PC2", as.numeric(val)]
pca[, stage_label := factor(stage_label, levels = names(pal_stage))]
p10 <- ggplot(pca, aes(PC1, PC2, fill = group, shape = group)) +
  geom_point(size = 3.0, color = "black", stroke = 0.3, alpha = 0.9) +
  scale_fill_manual(values = pal_group) +
  scale_shape_manual(values = c(sorted_thymocyte = 21, total_thymus = 22, "adult_T-ALL" = 24)) +
  labs(title = "2.10  PCA of GSE142522 (normal + adult T-ALL)",
       subtitle = sprintf("Top variable coding genes; PC1 %.1f%%  PC2 %.1f%%", 100 * pc1, 100 * pc2),
       x = sprintf("PC1 (%.1f%%)", 100 * pc1), y = sprintf("PC2 (%.1f%%)", 100 * pc2),
       fill = NULL, shape = NULL) +
  theme_sci() + theme(legend.position = "top")
save_both(p10, file.path(fig, "M2_2.10_PCA.pdf"), 7.4, 5.6)

p10b <- ggplot(pca[group != "total_thymus"], aes(PC1, PC2, fill = immunophenotype)) +
  geom_point(size = 3.0, color = "black", shape = 21, stroke = 0.3, alpha = 0.9) +
  scale_fill_manual(values = pal_ip) +
  labs(title = "PCA colored by immunophenotype",
       x = sprintf("PC1 (%.1f%%)", 100 * pc1), y = sprintf("PC2 (%.1f%%)", 100 * pc2), fill = NULL) +
  theme_sci() + theme(legend.position = "top")
save_both(p10b, file.path(fig, "M2_2.10_PCA_immunophenotype.pdf"), 7.4, 5.6)

# 2.11 nearest stage
near[, nearest_label := factor(nearest_label, levels = c("CD34+", "ISP", "Early cortical", "Late cortical", "SP4", "SP8"))]
nn <- near[, .N, by = nearest_label]
p11 <- ggplot(nn, aes(x = nearest_label, y = N, fill = nearest_label)) +
  geom_col(width = 0.72, color = "white") +
  geom_text(aes(label = N), vjust = -0.25, size = 3.6) +
  scale_fill_manual(values = pal_stage) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.18))) +
  labs(title = "2.11  Adult T-ALL nearest normal thymocyte stage",
       subtitle = "Pearson correlation to each stage using 1500 most variable stage genes",
       x = NULL, y = "n T-ALL") +
  theme_sci() + theme(legend.position = "none", axis.text.x = element_text(angle = 30, hjust = 1))
save_both(p11, file.path(fig, "M2_2.11_nearest_stage.pdf"), 7.2, 5.2)

p11b <- ggplot(near, aes(x = immunophenotype, fill = nearest_label)) +
  geom_bar(position = "fill", width = 0.72, color = "white", linewidth = 0.2) +
  scale_fill_manual(values = pal_stage) +
  scale_y_continuous(labels = function(x) paste0(100 * x, "%")) +
  labs(title = "Nearest stage by T-ALL immunophenotype", x = NULL, y = "Proportion", fill = "Nearest stage") +
  theme_sci() + theme(axis.text.x = element_text(angle = 20, hjust = 1))
save_both(p11b, file.path(fig, "M2_2.11_nearest_by_immunophenotype.pdf"), 7.4, 5.2)

# 2.12-2.15 T-ALL scatter vs maturity
tall <- sc[is_tall == TRUE]
ct12 <- suppressWarnings(cor.test(tall$ZEB1, tall$maturity_literature, method = "spearman", exact = FALSE))
scatter_gene <- function(y, ylab, title, file, color_var = "immunophenotype", pal = pal_ip) {
  ct0 <- suppressWarnings(cor.test(tall$maturity_literature, tall[[y]], method = "spearman", exact = FALSE))
  p <- ggplot(tall, aes(x = maturity_literature, y = .data[[y]], fill = .data[[color_var]])) +
    geom_point(size = 2.6, shape = 21, color = "black", stroke = 0.3, alpha = 0.9) +
    geom_smooth(method = "lm", formula = y ~ x, color = "grey20", se = TRUE, linewidth = 0.55,
                inherit.aes = FALSE, aes(x = maturity_literature, y = .data[[y]])) +
    scale_fill_manual(values = pal) +
    labs(title = title,
         subtitle = sprintf("n = %d adult T-ALL; Spearman rho = %.2f  %s",
                            nrow(tall), unname(ct0$estimate), fmt_p(ct0$p.value)),
         x = "Maturity score (literature)", y = ylab, fill = NULL) +
    theme_sci() + theme(legend.position = "top")
  save_both(p, file, 7.2, 5.4)
}
scatter_gene("ZEB1", "log2 count ZEB1", "2.12  Adult T-ALL: ZEB1 vs maturity",
             file.path(fig, "M2_2.12_ZEB1_vs_maturity.pdf"))
scatter_gene("LMO2", "log2 count LMO2", "2.13  Adult T-ALL: LMO2 vs maturity",
             file.path(fig, "M2_2.13_LMO2_vs_maturity.pdf"))
scatter_gene("ZEB2", "log2 count ZEB2", "2.14  Adult T-ALL: ZEB2 vs maturity",
             file.path(fig, "M2_2.14_ZEB2_vs_maturity.pdf"))
scatter_gene("ratio_lmo", "LMO2 − ZEB1 (log2)", "2.15  Adult T-ALL: LMO2/ZEB1 vs maturity",
             file.path(fig, "M2_2.15_LMO2_ZEB1_vs_maturity.pdf"))

p_ip <- plot_rain(tall, "immunophenotype", "ZEB1", "immunophenotype", pal_ip,
                  "log2 count ZEB1", "Adult T-ALL: ZEB1 by immunophenotype")
save_both(p_ip, file.path(fig, "M2_ZEB1_by_immunophenotype.pdf"), 6.6, 5.3)

# combined normal vs T-ALL ZEB1
sc[, cls := fifelse(is_tall, as.character(immunophenotype),
                    fifelse(group == "sorted_thymocyte", as.character(stage_label), "Total thymus"))]
ord <- c("CD34+", "ISP", "Early cortical", "Late cortical", "SP4", "SP8", "Total thymus",
         "Immature", "Pre_ab_Cortical", "TCR_Mature")
sc[, cls := factor(cls, levels = intersect(ord, unique(cls)))]
pbox <- ggplot(sc, aes(x = cls, y = ZEB1, fill = group)) +
  geom_boxplot(width = 0.55, outlier.shape = NA, alpha = 0.85) +
  add_points(size = 0.7, alpha = 0.65, width = 0.12) +
  scale_fill_manual(values = pal_group) +
  labs(title = "ZEB1 in normal thymocyte stages vs adult T-ALL",
       subtitle = "Sorted stages are single samples; T-ALL n = 41",
       x = NULL, y = "log2 count ZEB1", fill = NULL) +
  theme_sci() + theme(axis.text.x = element_text(angle = 40, hjust = 1), legend.position = "top")
save_both(pbox, file.path(fig, "M2_ZEB1_stages_vs_TALL.pdf"), 9.2, 5.6)

# 2.16-2.17 coefficient forest
if (nrow(reg)) {
  use <- reg[term != "intercept"]
  use[, lab := paste0(model, " | ", term)]
  setorder(use, coef)
  use[, lab := factor(lab, levels = unique(lab))]
  pcoef <- ggplot(use, aes(x = coef, y = lab)) +
    geom_vline(xintercept = 0, linetype = 2, color = "grey40") +
    geom_errorbar(aes(xmin = ci95_lo, xmax = ci95_hi), width = 0.18, linewidth = 0.45) +
    geom_point(aes(fill = coef < 0), shape = 21, size = 3.6, color = "black", stroke = 0.35) +
    scale_fill_manual(values = c("TRUE" = "#3C5488", "FALSE" = "#E64B35"), guide = "none") +
    geom_text(aes(label = sprintf("%.2f  %s", coef, p_star(p))),
              hjust = ifelse(use$coef >= 0, -0.15, 1.15), size = 2.7) +
    labs(title = "2.16–2.17  OLS coefficients (95% CI)",
         subtitle = "Sorted stages + adult T-ALL; maturity = literature score",
         x = "Coefficient", y = NULL) +
    theme_sci()
  save_both(pcoef, file.path(fig, "M2_2.16_2.17_coefficients.pdf"), 9.5, max(4.5, 0.38 * nrow(use) + 1.8))
}

# 2.18 TARGET / Pharmacotype
if (nrow(tgt) && "maturity_literature" %in% names(tgt)) {
  tgt[is.na(subtype) | subtype == "", subtype := "Unknown"]
  ct_t <- suppressWarnings(cor.test(tgt$ZEB1_log, tgt$maturity_literature, method = "spearman", exact = FALSE))
  p18 <- ggplot(tgt, aes(maturity_literature, ZEB1_log, fill = subtype)) +
    geom_point(size = 2.0, shape = 21, color = "black", stroke = 0.25, alpha = 0.85) +
    geom_smooth(method = "lm", formula = y ~ x, color = "grey20", se = TRUE, linewidth = 0.55,
                inherit.aes = FALSE, aes(x = maturity_literature, y = ZEB1_log)) +
    scale_fill_manual(values = pal_sub, drop = FALSE) +
    labs(title = "2.18  TARGET: ZEB1 vs thymocyte maturity score",
         subtitle = sprintf("Independent pediatric RNA-seq; n = %d; rho = %.2f  %s",
                            nrow(tgt), unname(ct_t$estimate), fmt_p(ct_t$p.value)),
         x = "Maturity score (same literature gene set)", y = "log2(TPM+1) ZEB1", fill = "Subtype") +
    theme_sci()
  save_both(p18, file.path(fig, "M2_2.18_TARGET_ZEB1_maturity.pdf"), 7.6, 5.5)
  fwrite(data.table(cohort = "TARGET", rho = unname(ct_t$estimate), p = ct_t$p.value, n = nrow(tgt)),
         file.path(tab, "M2_2.18_replication.tsv"), sep = "\t")
}
if (nrow(ph) && "maturity_literature" %in% names(ph) && "ZEB1_log" %in% names(ph)) {
  if (!"etp" %in% names(ph)) ph[, etp := "notETP"]
  ph[, etp_lab := fifelse(etp == "ETP", "ETP", "notETP")]
  ct_p <- suppressWarnings(cor.test(ph$ZEB1_log, ph$maturity_literature, method = "spearman", exact = FALSE))
  p18b <- ggplot(ph, aes(maturity_literature, ZEB1_log, fill = etp_lab)) +
    geom_point(size = 2.4, shape = 21, color = "black", stroke = 0.3, alpha = 0.9) +
    geom_smooth(method = "lm", formula = y ~ x, color = "grey20", se = TRUE, linewidth = 0.55,
                inherit.aes = FALSE, aes(x = maturity_literature, y = ZEB1_log)) +
    scale_fill_manual(values = c(ETP = "#E64B35", notETP = "#4DBBD5"), name = NULL) +
    labs(title = "2.18  Pharmacotype: ZEB1 vs thymocyte maturity score",
         subtitle = sprintf("Independent pediatric RNA-seq; n = %d; rho = %.2f  %s",
                            nrow(ph), unname(ct_p$estimate), fmt_p(ct_p$p.value)),
         x = "Maturity score", y = "log2(FPKM+1) ZEB1") +
    theme_sci()
  save_both(p18b, file.path(fig, "M2_2.18_Pharmacotype_ZEB1_maturity.pdf"), 7.4, 5.5)
  fwrite(data.table(cohort = c("TARGET", "Pharmacotype"),
                    rho = c(if (exists("ct_t")) unname(ct_t$estimate) else NA_real_, unname(ct_p$estimate)),
                    p = c(if (exists("ct_t")) ct_t$p.value else NA_real_, ct_p$p.value),
                    n = c(nrow(tgt), nrow(ph))),
         file.path(tab, "M2_2.18_replication.tsv"), sep = "\t")
}

# lollipop of key contrasts
if ("hedges_g" %in% names(ct)) {
  cc <- ct[is.finite(hedges_g)]
  if (nrow(cc)) {
    setorder(cc, hedges_g)
    cc[, contrast := factor(contrast, levels = unique(contrast))]
    pl <- ggplot(cc, aes(x = hedges_g, y = contrast)) +
      geom_vline(xintercept = 0, linetype = 2, color = "grey40") +
      geom_segment(aes(x = 0, xend = hedges_g, yend = contrast), color = "#B07D3A", linewidth = 1) +
      geom_errorbar(aes(xmin = ci95_lo, xmax = ci95_hi), width = 0.18, linewidth = 0.45) +
      geom_point(aes(fill = hedges_g < 0), shape = 21, size = 3.6, color = "black") +
      scale_fill_manual(values = c("TRUE" = "#3C5488", "FALSE" = "#E64B35"), guide = "none") +
      labs(title = "ZEB1 effect sizes in GSE142522",
           subtitle = "Negative g = lower in group A; sorted stages n=6 are not i.i.d. replicates",
           x = "Hedges' g", y = NULL) +
      theme_sci()
    save_both(pl, file.path(fig, "M2_ZEB1_effect_lollipop.pdf"), 8.4, max(4.2, 0.4 * nrow(cc) + 1.6))
  }
}

sink(file.path(tab, "sessionInfo.txt")); print(sessionInfo()); sink()
cat("MODULE2 ANALYZE DONE\n")
cat("FIG", length(list.files(fig)), "TAB", length(list.files(tab)), "\n")

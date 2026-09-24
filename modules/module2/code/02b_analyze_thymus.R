# Module 2 rebuild figures. FACS thymus GSE195812 + pediatric GSE206710.
# Recipes: code2 229 UMAP, 238 bar+beeswarm, 210 heatmap, 88 forest-like.
suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})
has_beeswarm <- requireNamespace("ggbeeswarm", quietly = TRUE)
has_cht <- requireNamespace("ComplexHeatmap", quietly = TRUE)
has_circlize <- requireNamespace("circlize", quietly = TRUE)
if (has_beeswarm) library(ggbeeswarm)

root <- "/data-b/liangfuhua/projects/ZEB1_TALL_analysis/module2"
fig <- file.path(root, "figures"); tab <- file.path(root, "tables")
dir.create(fig, FALSE, TRUE)

theme_sci <- function(base_size = 11) {
  theme_classic(base_size = base_size) +
    theme(axis.text = element_text(color = "black"), axis.title = element_text(color = "black"),
          axis.line = element_line(linewidth = 0.45, color = "black"),
          plot.title = element_text(face = "bold", size = base_size + 1),
          plot.subtitle = element_text(size = base_size - 1, color = "grey25"),
          strip.background = element_blank())
}
save_both <- function(plot, file, width, height) {
  tryCatch({
    pdf(file, width = width, height = height, useDingbats = FALSE); print(plot); dev.off()
    png(sub("\\.pdf$", ".png", file), width = width * 160, height = height * 160, res = 160)
    print(plot); dev.off(); message("saved ", basename(file))
  }, error = function(e) message("FAIL ", basename(file), ": ", conditionMessage(e)))
}
run <- function(label, fn) tryCatch({ fn(); message("OK ", label) },
                                    error = function(e) message("FAIL ", label, ": ", conditionMessage(e)))
rd <- function(f) {
  p <- file.path(tab, f)
  if (!file.exists(p)) return(data.table())
  fread(p)
}
add_pts <- function(size = 0.35, alpha = 0.35, width = 0.18) {
  if (has_beeswarm) geom_quasirandom(size = size, alpha = alpha, width = width, color = "grey20")
  else geom_jitter(size = size, alpha = alpha, width = width, color = "grey20")
}

pal_stage <- c(
  DN1 = "#3C5488", DN2 = "#4DBBD5", DN3 = "#00A087", ISP = "#91D1C2",
  DP_CD3neg = "#F39B7F", DP_CD3pos = "#E64B35", CD4SP = "#7E6148", CD8SP = "#B09C85",
  CD3neg = "#4DBBD5", CD3pos = "#E64B35",
  DN_early = "#3C5488", DN_late_ISP = "#00A087", DP_immature = "#F39B7F",
  DP_mature = "#E64B35", other = "#8491B4"
)
stage_ord <- c("DN1", "DN2", "DN3", "ISP", "DP_CD3neg", "DP_CD3pos", "CD4SP", "CD8SP")

c195 <- rd("M2_r2_GSE195812_plot_cells.tsv.gz")
m195 <- rd("M2_r2_GSE195812_sample_means.tsv")
tr195 <- rd("M2_r2_GSE195812_gene_vs_stage.tsv")
c206 <- rd("M2_r2_GSE206710_plot_cells.tsv.gz")
d206 <- rd("M2_r2_GSE206710_donor_fraction.tsv")
pat <- rd("M2_r2_map_GSE227122_patient.tsv")
comp <- rd("M2_r2_map_GSE227122_stage_comp.tsv")
cent_z <- rd("M2_r2_GSE195812_ZEB1_by_stage.tsv")
stg_g <- rd("M2_r2_GSE195812_stage_gene_means.tsv")

if (nrow(c195) && "facs_stage" %in% names(c195)) {
  c195[, facs_stage := factor(facs_stage, levels = intersect(stage_ord, unique(facs_stage)))]
}
if (nrow(m195) && "stage" %in% names(m195)) {
  m195[, stage := factor(stage, levels = intersect(stage_ord, unique(stage)))]
}

run("2.1 UMAP FACS stage", function() {
  if (!nrow(c195)) stop("no 195 cells")
  p <- ggplot(c195, aes(UMAP1, UMAP2, color = facs_stage)) +
    geom_point(size = 0.15, alpha = 0.75) +
    scale_color_manual(values = pal_stage, na.value = "grey80") +
    labs(title = "GSE195812 postnatal thymus", subtitle = "FACS-sorted libraries; 6 donors pooled",
         color = "FACS") + theme_sci() +
    theme(axis.ticks = element_blank(), axis.text = element_blank())
  save_both(p, file.path(fig, "M2_r2_2.1_umap_facs.pdf"), 7.2, 5.8)
})

run("2.2 UMAP ZEB1", function() {
  if (!nrow(c195) || !"g_ZEB1" %in% names(c195)) stop("no ZEB1")
  p <- ggplot(c195, aes(UMAP1, UMAP2, color = g_ZEB1)) +
    geom_point(size = 0.15, alpha = 0.8) +
    scale_color_gradientn(colors = c("#3C5488", "#EEEEEE", "#E64B35")) +
    labs(title = "ZEB1 in normal thymus", color = "ZEB1") + theme_sci() +
    theme(axis.ticks = element_blank(), axis.text = element_blank())
  save_both(p, file.path(fig, "M2_r2_2.2_umap_ZEB1.pdf"), 7.0, 5.8)
})

run("2.3 UMAP LMO2", function() {
  if (!"g_LMO2" %in% names(c195)) stop("no LMO2")
  p <- ggplot(c195, aes(UMAP1, UMAP2, color = g_LMO2)) +
    geom_point(size = 0.15, alpha = 0.8) +
    scale_color_gradientn(colors = c("#3C5488", "#EEEEEE", "#E64B35")) +
    labs(title = "LMO2 in normal thymus", color = "LMO2") + theme_sci() +
    theme(axis.ticks = element_blank(), axis.text = element_blank())
  save_both(p, file.path(fig, "M2_r2_2.3_umap_LMO2.pdf"), 7.0, 5.8)
})

run("2.4 ZEB1 by FACS T-lineage", function() {
  d <- c195[lineage == "T"]
  p <- ggplot(d, aes(facs_stage, g_ZEB1, fill = facs_stage)) +
    geom_violin(scale = "width", linewidth = 0.3, alpha = 0.85) +
    add_pts(size = 0.12, alpha = 0.12, width = 0.22) +
    scale_fill_manual(values = pal_stage, guide = "none") +
    labs(title = "ZEB1 along FACS stages (T-lineage cells)",
         subtitle = "Descriptive cell distributions; biological n = 1 pooled-donor library/stage",
         x = NULL, y = "ZEB1 (log1p)") + theme_sci() +
    theme(axis.text.x = element_text(angle = 35, hjust = 1))
  save_both(p, file.path(fig, "M2_r2_2.4_ZEB1_facs_violin.pdf"), 8.4, 5.4)
})

run("2.5 key genes sample-mean trajectories", function() {
  if (!nrow(m195)) stop("no means")
  genes <- intersect(c("g_ZEB1", "g_ZEB2", "g_LMO2", "g_LYL1", "g_TCF7", "g_GATA3",
                       "g_BCL11B", "g_IL7R", "g_CD34", "g_MEF2C"), names(m195))
  long <- melt(m195, id.vars = c("sample", "stage"), measure.vars = genes,
               variable.name = "gene", value.name = "expr")
  long[, gene := sub("^g_", "", gene)]
  p <- ggplot(long, aes(stage, expr, group = gene, color = gene)) +
    geom_line(linewidth = 0.7) + geom_point(size = 2) +
    facet_wrap(~gene, scales = "free_y", ncol = 5) +
    labs(title = "Library-mean expression along FACS order",
         subtitle = "One pooled-donor library per stage", x = NULL, y = "mean log1p") +
    theme_sci() + theme(axis.text.x = element_text(angle = 40, hjust = 1), legend.position = "none")
  save_both(p, file.path(fig, "M2_r2_2.5_gene_trajectories.pdf"), 11.5, 6.5)
})

run("2.6 stage Spearman", function() {
  if (!nrow(tr195)) stop("no trends")
  d <- copy(tr195)
  d[, lab := ifelse(gene %in% c("ZEB1", "ZEB2", "LMO2", "LYL1", "TCF7", "GATA3", "BCL11B", "CD34", "IL7R"), gene, "")]
  p <- ggplot(d, aes(rho_vs_stage, -log10(pmax(p, 1e-12)), color = abs(rho_vs_stage))) +
    geom_point(size = 2.2, alpha = 0.85) +
    geom_text(aes(label = lab), nudge_y = 0.15, size = 3, color = "black") +
    scale_color_gradient(low = "#8491B4", high = "#E64B35") +
    labs(title = "Gene vs FACS stage order (library means)",
         x = "Spearman rho (stage index)", y = "-log10 P") + theme_sci()
  save_both(p, file.path(fig, "M2_r2_2.6_gene_vs_stage.pdf"), 6.8, 5.6)
})

run("2.7 stage heatmap", function() {
  if (!nrow(stg_g) || !has_cht || !has_circlize) stop("no heatmap")
  mat <- as.matrix(stg_g[, setdiff(names(stg_g), names(stg_g)[1]), with = FALSE])
  if (names(stg_g)[1] != "") {
    rn <- stg_g[[1]]
    cn <- names(stg_g)[-1]
    mat <- as.matrix(stg_g[, -1, with = FALSE])
    rownames(mat) <- rn
    colnames(mat) <- sub("^g_", "", cn)
  }
  keep <- intersect(c("ZEB1", "ZEB2", "LMO2", "LYL1", "MEF2C", "SPI1", "CD34", "KIT",
                      "IL7R", "CD1A", "RAG1", "DNTT", "CD4", "CD8A", "CD3D",
                      "BCL11B", "TCF7", "GATA3", "LEF1"), colnames(mat))
  if (!length(keep)) stop("no genes")
  mat <- t(mat[, keep, drop = FALSE])
  mat <- mat[, intersect(stage_ord, colnames(mat)), drop = FALSE]
  mat <- t(scale(t(mat)))
  pdf(file.path(fig, "M2_r2_2.7_stage_heatmap.pdf"), width = 7.2, height = 6.4)
  ComplexHeatmap::Heatmap(mat, name = "z", cluster_columns = FALSE,
                          column_title = "T-lineage mean (z)")
  dev.off()
})

run("2.8 GSE206710 UMAP fraction", function() {
  if (!nrow(c206)) stop("no 206")
  p <- ggplot(c206, aes(UMAP1, UMAP2, color = fraction)) +
    geom_point(size = 0.15, alpha = 0.75) +
    scale_color_manual(values = pal_stage) +
    facet_wrap(~donor) +
    labs(title = "GSE206710 pediatric thymus", subtitle = "3 donors; CD3− vs CD3+") +
    theme_sci() + theme(axis.ticks = element_blank(), axis.text = element_blank())
  save_both(p, file.path(fig, "M2_r2_2.8_umap_206_donor.pdf"), 9.5, 4.2)
})

run("2.9 GSE206710 donor ZEB1", function() {
  if (!nrow(d206)) stop("no donor tab")
  p <- ggplot(d206, aes(fraction, ZEB1, color = donor, group = donor)) +
    geom_point(size = 3) + geom_line(linewidth = 0.8) +
    labs(title = "Donor-mean ZEB1: CD3− vs CD3+ T-lineage",
         subtitle = "Biological n = 3 pediatric donors", x = NULL, y = "mean ZEB1") +
    theme_sci()
  save_both(p, file.path(fig, "M2_r2_2.9_donor_ZEB1_fraction.pdf"), 5.8, 5.2)
})

run("2.10 map patient majority", function() {
  if (!nrow(pat)) stop("no map")
  p <- ggplot(pat, aes(majority_stage)) +
    geom_bar(fill = "#3C5488", width = 0.7) +
    labs(title = "GSE227122 Dx malignant cells: majority FACS stage",
         subtitle = "Mapped by cosine to GSE195812 T-lineage centroids",
         x = NULL, y = "patients") + theme_sci() +
    theme(axis.text.x = element_text(angle = 35, hjust = 1))
  save_both(p, file.path(fig, "M2_r2_2.10_patient_majority_stage.pdf"), 6.8, 5.0)
})

run("2.11 map composition stacked", function() {
  if (!nrow(comp)) stop("no comp")
  comp[, nearest_stage := factor(nearest_stage, levels = intersect(stage_ord, unique(nearest_stage)))]
  p <- ggplot(comp, aes(patient, frac, fill = nearest_stage)) +
    geom_col(width = 0.85) +
    scale_fill_manual(values = pal_stage) +
    labs(title = "T-ALL malignant cells mapped to FACS stages",
         x = NULL, y = "fraction of malignant cells", fill = "nearest") +
    theme_sci()
  save_both(p, file.path(fig, "M2_r2_2.11_patient_stage_stack.pdf"), 8.5, 5.2)
})

run("2.12 residual vs entropy", function() {
  if (!nrow(pat) || !"ZEB1_residual" %in% names(pat)) stop("no residual")
  p <- ggplot(pat, aes(map_entropy, ZEB1_residual, size = n_cells)) +
    geom_hline(yintercept = 0, linetype = 2, color = "grey50") +
    geom_point(shape = 21, fill = "#E64B35", alpha = 0.85) +
    labs(title = "ZEB1 residual after FACS-stage expectation",
         subtitle = "Residual = observed − mean(ZEB1 | nearest normal stage)",
         x = "mapping entropy (normalized)", y = "ZEB1 residual") + theme_sci()
  save_both(p, file.path(fig, "M2_r2_2.12_residual_entropy.pdf"), 6.6, 5.4)
})

run("2.13 ZEB1 vs expected", function() {
  if (!nrow(pat) || !"g_ZEB1" %in% names(pat)) stop("no zeb")
  # expected from majority stage if present in centroid table
  p <- ggplot(pat, aes(majority_stage, g_ZEB1)) +
    geom_boxplot(width = 0.5, outlier.shape = NA, fill = "#4DBBD5", alpha = 0.8) +
    geom_point(size = 2.4, alpha = 0.85) +
    labs(title = "T-ALL malignant-mean ZEB1 by mapped stage",
         x = NULL, y = "patient mean ZEB1") + theme_sci() +
    theme(axis.text.x = element_text(angle = 35, hjust = 1))
  save_both(p, file.path(fig, "M2_r2_2.13_TALL_ZEB1_by_mapped_stage.pdf"), 7.0, 5.2)
})

run("2.14 DPT density", function() {
  if (!"dpt_pseudotime" %in% names(c195)) stop("no dpt")
  d <- c195[lineage == "T" & is.finite(dpt_pseudotime) & dpt_pseudotime < Inf]
  if (!nrow(d)) stop("DPT non-finite")
  p <- ggplot(d, aes(dpt_pseudotime, fill = facs_stage, color = facs_stage)) +
    geom_density(alpha = 0.25, linewidth = 0.6) +
    scale_fill_manual(values = pal_stage) + scale_color_manual(values = pal_stage) +
    labs(title = "Diffusion pseudotime by FACS stage", x = "DPT", y = "density") + theme_sci()
  save_both(p, file.path(fig, "M2_r2_2.14_dpt_density.pdf"), 7.4, 5.2)
})

message("M2 r2 figures done")
message("FIG ", length(list.files(fig, pattern = "^M2_r2_")))

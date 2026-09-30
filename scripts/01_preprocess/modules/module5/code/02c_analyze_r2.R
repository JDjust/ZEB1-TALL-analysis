# Module 5 round-2 figures: thymus mapping, cores, composition, CITE ADT.
# Recipes: code2 78 stacked, 88 forest, 193 scatter, 238 bar.
# Biological n is patient.
suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})
has_beeswarm <- requireNamespace("ggbeeswarm", quietly = TRUE)
if (has_beeswarm) library(ggbeeswarm)

root <- "/data-b/liangfuhua/projects/ZEB1_TALL_analysis/module5"
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
fmt_p <- function(p) ifelse(is.na(p), "NA", ifelse(p < 0.001, "P < 0.001", sprintf("P = %.3f", p)))

pal_stage <- c(
  DN1 = "#3C5488", DN2 = "#4DBBD5", DN3 = "#00A087", ISP = "#91D1C2",
  DP_CD3neg = "#F39B7F", DP_CD3pos = "#E64B35", CD4SP = "#7E6148", CD8SP = "#B09C85",
  `DN-like` = "#3C5488", `DP-like` = "#E64B35", `SP-like` = "#7E6148"
)
stage_ord <- c("DN1", "DN2", "DN3", "ISP", "DP_CD3neg", "DP_CD3pos", "CD4SP", "CD8SP")

pat <- rd("M5_r2_patient_merged.tsv")
comp <- rd("M5_r2_map_stage_comp.tsv")
assoc <- rd("M5_r2_ZEB1_vs_mapped.tsv")
core <- rd("M5_r2_core_patient.tsv")
cite <- rd("M5_r2_cite_ZEB1_vs_ADT.tsv")
cite_m <- rd("M5_r2_cite_ADT_RNA_merged.tsv")
compbin <- rd("M5_r2_composition_by_ZEB1bin.tsv")

run("5r2.1 mapped stage stack", function() {
  if (!nrow(comp)) stop("no composition")
  d <- copy(comp)
  d[, nearest_stage := factor(nearest_stage, levels = intersect(stage_ord, unique(nearest_stage)))]
  p <- ggplot(d, aes(patient, frac, fill = nearest_stage)) +
    geom_col(width = 0.85, color = "white", linewidth = 0.2) +
    scale_fill_manual(values = pal_stage, na.value = "grey80") +
    scale_y_continuous(expand = expansion(mult = c(0, 0.02))) +
    labs(title = "GSE227122 malignant cells mapped to FACS thymus stages",
         subtitle = "Reference = GSE195812; patient is the biological n",
         x = NULL, y = "Fraction of cells", fill = "Nearest stage") +
    theme_sci() + theme(axis.text.x = element_text(angle = 45, hjust = 1))
  save_both(p, file.path(fig, "M5_r2_5.7_stage_stack.pdf"), 8.4, 5.4)
})

run("5r2.2 entropy vs ZEB1", function() {
  d <- copy(pat)
  if (!nrow(d) || !"g_ZEB1" %in% names(d) || !"map_entropy" %in% names(d)) stop("no merge")
  d <- d[is.finite(g_ZEB1) & is.finite(map_entropy)]
  ct <- cor.test(d$g_ZEB1, d$map_entropy, method = "spearman", exact = FALSE)
  p <- ggplot(d, aes(g_ZEB1, map_entropy)) +
    geom_point(size = 2.6, color = "#E64B35") +
    geom_smooth(method = "lm", se = TRUE, linewidth = 0.5, color = "grey20") +
    labs(title = "Mapped developmental entropy vs ZEB1",
         subtitle = sprintf("n = %d patients; Spearman rho = %.2f, %s",
                            nrow(d), unname(ct$estimate), fmt_p(ct$p.value)),
         x = "Malignant-cell mean ZEB1", y = "Mapping entropy") + theme_sci()
  save_both(p, file.path(fig, "M5_r2_5.25_entropy_vs_ZEB1.pdf"), 6.2, 5.4)
})

run("5r2.3 residual vs ZEB1", function() {
  d <- copy(pat)
  if (!"ZEB1_residual" %in% names(d)) stop("no residual")
  d <- d[is.finite(g_ZEB1) & is.finite(ZEB1_residual)]
  ct <- cor.test(d$g_ZEB1, d$ZEB1_residual, method = "spearman", exact = FALSE)
  p <- ggplot(d, aes(g_ZEB1, ZEB1_residual, fill = majority_stage)) +
    geom_hline(yintercept = 0, linetype = 2, color = "grey50") +
    geom_point(size = 3, shape = 21, color = "grey20") +
    scale_fill_manual(values = pal_stage, na.value = "grey80") +
    labs(title = "ZEB1 residual after expected normal-stage mean",
         subtitle = sprintf("n = %d; rho = %.2f, %s", nrow(d), unname(ct$estimate), fmt_p(ct$p.value)),
         x = "Observed ZEB1", y = "Residual ZEB1", fill = "Majority stage") + theme_sci()
  save_both(p, file.path(fig, "M5_r2_5.16_residual.pdf"), 6.8, 5.4)
})

run("5r2.4 mapping forest", function() {
  if (!nrow(assoc)) stop("no assoc")
  d <- copy(assoc)
  keep <- d[feature %in% c("map_entropy", "ZEB1_residual", "core_pos", "core_neg",
                           "p_DN1", "p_DN2", "p_DN3", "p_ISP", "p_DP_CD3neg",
                           "p_DP_CD3pos", "p_CD4SP", "p_CD8SP", "sig_GATA3prog",
                           "sig_LMO2prog", "sig_ETP", "sig_Tcell_diff")]
  if (!nrow(keep)) keep <- d[order(p)][1:min(14, nrow(d))]
  setorder(keep, rho)
  keep[, feature := factor(feature, levels = feature)]
  p <- ggplot(keep, aes(rho, feature)) +
    geom_vline(xintercept = 0, linetype = 2, color = "grey40") +
    geom_point(aes(fill = p < 0.05), shape = 21, size = 3.2, color = "grey20") +
    scale_fill_manual(values = c("TRUE" = "#E64B35", "FALSE" = "#8491B4"), guide = "none") +
    labs(title = "Patient-level Spearman: ZEB1 vs mapped / core features",
         subtitle = "GSE227122 Dx; n = 10. Filled = nominal P < 0.05",
         x = "rho vs ZEB1", y = NULL) + theme_sci()
  save_both(p, file.path(fig, "M5_r2_5.12_mapping_forest.pdf"), 7.2, 5.8)
})

run("5r2.5 core vs ZEB1", function() {
  d <- if (nrow(core)) copy(core) else copy(pat)
  if (!"core_pos" %in% names(d) || !"g_ZEB1" %in% names(d)) stop("no core")
  d <- d[is.finite(g_ZEB1) & is.finite(core_pos)]
  ct <- cor.test(d$g_ZEB1, d$core_pos, method = "spearman", exact = FALSE)
  p <- ggplot(d, aes(g_ZEB1, core_pos)) +
    geom_point(size = 2.8, color = "#E64B35") +
    geom_smooth(method = "lm", se = TRUE, linewidth = 0.5, color = "grey20") +
    labs(title = "M4 subtype-adjusted ZEB1-positive core in GSE227122",
         subtitle = sprintf("Patient malignant means, n = %d; rho = %.2f, %s",
                            nrow(d), unname(ct$estimate), fmt_p(ct$p.value)),
         x = "ZEB1", y = "Positive core score") + theme_sci()
  save_both(p, file.path(fig, "M5_r2_5.9_core_vs_ZEB1.pdf"), 6.2, 5.4)
})

run("5r2.6 CITE ADT forest", function() {
  if (!nrow(cite)) stop("no CITE assoc")
  d <- copy(cite)
  setorder(d, rho)
  d[, feature := factor(feature, levels = feature)]
  p <- ggplot(d, aes(rho, feature)) +
    geom_vline(xintercept = 0, linetype = 2) +
    geom_point(shape = 21, size = 3.2, fill = "#4DBBD5", color = "grey20") +
    labs(title = "GSE248287: surface ADT vs ZEB1 RNA",
         subtitle = "Patient means; WNN skipped (no muon). Combined RNA mtx truncated.",
         x = "Spearman rho vs ZEB1", y = NULL) + theme_sci()
  save_both(p, file.path(fig, "M5_r2_5.18_CITE_ADT_forest.pdf"), 7.0, 5.2)
})

run("5r2.7 CITE CD1A scatter", function() {
  if (!nrow(cite_m)) stop("no CITE merge")
  cd1 <- grep("^ADT_CD1", names(cite_m), value = TRUE)
  if (!length(cd1) || !"g_ZEB1" %in% names(cite_m)) stop("no CD1A")
  d <- copy(cite_m)
  d[, y := d[[cd1[1]]]]
  d <- d[is.finite(g_ZEB1) & is.finite(y)]
  ct <- cor.test(d$g_ZEB1, d$y, method = "spearman", exact = FALSE)
  p <- ggplot(d, aes(g_ZEB1, y)) +
    geom_point(size = 2.8, color = "#3C5488") +
    geom_smooth(method = "lm", se = TRUE, linewidth = 0.5, color = "grey20") +
    labs(title = "CITE-seq CD1a protein vs ZEB1 RNA",
         subtitle = sprintf("n = %d Dx patients; rho = %.2f, %s", nrow(d), unname(ct$estimate), fmt_p(ct$p.value)),
         x = "ZEB1 RNA", y = cd1[1]) + theme_sci()
  save_both(p, file.path(fig, "M5_r2_5.18_CITE_CD1A.pdf"), 6.2, 5.4)
})

message("M5 R2 FIGURES DONE")

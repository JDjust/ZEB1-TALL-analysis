# Module 5 figures. Checklist IDs 5.1-5.42. Recipes: code2 229 UMAP,
# 238 bar+beeswarm, 78 stacked composition, 233 paired, 88 forest,
# 193 scatter, 346 bidirectional, 40 density, 210 heatmap.
# Biological n is patient, not cell.

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})
has_ggrepel <- requireNamespace("ggrepel", quietly = TRUE)
has_cht <- requireNamespace("ComplexHeatmap", quietly = TRUE)
has_circlize <- requireNamespace("circlize", quietly = TRUE)
has_patch <- requireNamespace("patchwork", quietly = TRUE)
if (has_patch) library(patchwork)

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
    print(plot); dev.off()
    message("saved ", basename(file))
  }, error = function(e) message("FAIL ", basename(file), ": ", conditionMessage(e)))
}
run <- function(label, fn) tryCatch({ fn(); message("OK ", label) },
                                    error = function(e) message("FAIL ", label, ": ", conditionMessage(e)))
rd <- function(f) {
  p <- file.path(tab, f)
  if (!file.exists(p)) return(data.table())
  fread(p)
}

pal_tp <- c(Dx = "#E64B35", EOI = "#4DBBD5", Rel = "#00A087")
pal_lin <- c(T = "#E64B35", NK = "#4DBBD5", B = "#00A087", Myeloid = "#3C5488",
             Ery = "#F39B7F", Unknown = "#8491B4")
pal_mal <- c(malignant = "#E64B35", nonmalignant = "#8491B4")

qc <- rd("M5_5.1_sample_qc.tsv")
qcell <- rd("M5_qc_cells_prefilter.tsv.gz")
cells <- rd("M5_plot_cells.tsv.gz")
comp <- rd("M5_5.13_composition.tsv")
malc <- rd("M5_malignant_composition.tsv")
pat <- rd("M5_5.36_patient_pseudobulk.tsv")
pairs <- rd("M5_5.37_Dx_EOI_paired.tsv")
assoc <- rd("M5_5.41_patient_ZEB1_associations.tsv")
clsum <- rd("M5_cluster_summary.tsv")
v248 <- rd("M5_5.42_ZEB1_patient_assoc.tsv")
pb248 <- rd("M5_5.42_patient_pseudobulk.tsv")

if (nrow(cells) && "UMAP1" %in% names(cells)) {
  cells[, UMAP1 := as.numeric(UMAP1)]
  cells[, UMAP2 := as.numeric(UMAP2)]
}

umap_base <- function(d, color_col, title, pal = NULL, discrete = TRUE) {
  dd <- copy(d)
  dd[, colorv := dd[[color_col]]]
  p <- ggplot(dd, aes(UMAP1, UMAP2, color = colorv)) +
    geom_point(size = 0.15, alpha = 0.7) +
    labs(title = title, x = "UMAP1", y = "UMAP2", color = NULL) +
    theme_sci() +
    theme(axis.ticks = element_blank(), axis.text = element_blank())
  if (discrete && !is.null(pal)) p <- p + scale_color_manual(values = pal, na.value = "grey80")
  if (!discrete) p <- p + scale_color_gradientn(colors = c("#3C5488", "#EEEEEE", "#E64B35"))
  p + theme(legend.position = "right")
}

run("5.1 cells per sample", function() {
  if (!nrow(qc)) stop("no qc")
  setorder(qc, patient, timepoint)
  qc[, sample := factor(sample, levels = sample)]
  p <- ggplot(qc, aes(n_cells, sample, fill = timepoint)) +
    geom_col(width = 0.75, color = "white") +
    geom_text(aes(label = n_cells), hjust = -0.08, size = 3) +
    scale_fill_manual(values = pal_tp) +
    labs(title = "5.1  Cells per sample (capped at 8000 before QC)",
         x = "Cells", y = NULL, fill = NULL) + theme_sci()
  save_both(p, file.path(fig, "M5_5.1_cells_per_sample.pdf"), 7.6, 6.4)
})

run("5.2 nFeature", function() {
  if (!nrow(qcell)) stop("no qcell")
  p <- ggplot(qcell, aes(sample, n_genes_by_counts, fill = timepoint)) +
    geom_boxplot(outlier.size = 0.2, width = 0.7, alpha = 0.9) +
    scale_fill_manual(values = pal_tp) +
    labs(title = "5.2  nFeature (genes / cell)", x = NULL, y = "nFeature", fill = NULL) +
    theme_sci() + theme(axis.text.x = element_text(angle = 45, hjust = 1))
  save_both(p, file.path(fig, "M5_5.2_nFeature.pdf"), 9.2, 5.4)
})

run("5.3 nCount", function() {
  if (!nrow(qcell)) stop("no qcell")
  p <- ggplot(qcell, aes(sample, total_counts, fill = timepoint)) +
    geom_boxplot(outlier.size = 0.2, width = 0.7, alpha = 0.9) +
    scale_y_log10() + scale_fill_manual(values = pal_tp) +
    labs(title = "5.3  nCount (UMIs / cell)", x = NULL, y = "nCount", fill = NULL) +
    theme_sci() + theme(axis.text.x = element_text(angle = 45, hjust = 1))
  save_both(p, file.path(fig, "M5_5.3_nCount.pdf"), 9.2, 5.4)
})

run("5.4 mito", function() {
  if (!nrow(qcell)) stop("no qcell")
  set.seed(1)
  d <- qcell[sample(.N, min(.N, 25000))]
  p <- ggplot(d, aes(total_counts, pct_counts_mt, color = timepoint)) +
    geom_point(size = 0.25, alpha = 0.35) +
    geom_hline(yintercept = 20, linetype = 2, color = "grey40") +
    scale_x_log10() + scale_color_manual(values = pal_tp) +
    labs(title = "5.4  Mitochondrial percent vs nCount",
         subtitle = "Dashed line = 20% filter", x = "nCount", y = "mito %", color = NULL) +
    theme_sci()
  save_both(p, file.path(fig, "M5_5.4_mito.pdf"), 7.4, 5.6)
})

run("5.5 doublet", function() {
  if (!nrow(cells) || !"doublet_score" %in% names(cells)) stop("no doublet")
  d <- cells[is.finite(as.numeric(doublet_score))]
  if (!nrow(d)) stop("doublet all NA (scrublet unavailable)")
  p <- ggplot(d, aes(as.numeric(doublet_score), fill = timepoint)) +
    geom_histogram(bins = 50, color = "white", linewidth = 0.15, position = "identity", alpha = 0.7) +
    scale_fill_manual(values = pal_tp) +
    labs(title = "5.5  Scrublet doublet score", x = "doublet score", y = "Cells", fill = NULL) +
    theme_sci()
  save_both(p, file.path(fig, "M5_5.5_doublet.pdf"), 7.0, 5.0)
})

run("5.6 umap patient", function() {
  if (!nrow(cells)) stop("no cells")
  p <- umap_base(cells, "patient", "5.6  UMAP colored by patient (batch)")
  save_both(p, file.path(fig, "M5_5.6_umap_patient.pdf"), 8.2, 6.4)
})

run("5.7 harmony vs unint", function() {
  if (!nrow(cells) || !"UMAP1_unint" %in% names(cells)) stop("no unint umap")
  a <- ggplot(cells, aes(UMAP1_unint, UMAP2_unint, color = patient)) +
    geom_point(size = 0.12, alpha = 0.65) +
    labs(title = "Before Harmony", x = "UMAP1", y = "UMAP2") + theme_sci() +
    theme(legend.position = "none", axis.text = element_blank(), axis.ticks = element_blank())
  b <- ggplot(cells, aes(UMAP1, UMAP2, color = patient)) +
    geom_point(size = 0.12, alpha = 0.65) +
    labs(title = "After Harmony (sample batch)", x = "UMAP1", y = "UMAP2") + theme_sci() +
    theme(axis.text = element_blank(), axis.ticks = element_blank())
  if (has_patch) save_both(a + b, file.path(fig, "M5_5.7_umap_harmony_vs_unint.pdf"), 12.5, 5.8)
  else {
    save_both(a, file.path(fig, "M5_5.7_umap_unint.pdf"), 6.2, 5.6)
    save_both(b, file.path(fig, "M5_5.7_umap_harmony.pdf"), 7.4, 5.6)
  }
})

run("5.8 leiden UMAP", function() {
  if (!nrow(cells) || !"leiden" %in% names(cells)) stop("no leiden")
  lab <- cells[, .(UMAP1 = median(UMAP1), UMAP2 = median(UMAP2)), by = leiden]
  p <- ggplot(cells, aes(UMAP1, UMAP2, color = factor(leiden))) +
    geom_point(size = 0.15, alpha = 0.7) +
    geom_point(data = lab, size = 5, color = "white") +
    geom_text(data = lab, aes(label = leiden), size = 3, color = "black") +
    labs(title = "5.8  Global UMAP (leiden 0.6, code2 229)", x = "UMAP1", y = "UMAP2", color = NULL) +
    theme_sci() + theme(legend.position = "none", axis.text = element_blank(), axis.ticks = element_blank())
  save_both(p, file.path(fig, "M5_5.8_umap_leiden.pdf"), 7.2, 6.4)
})

run("5.9 marker heatmap", function() {
  mm <- rd("M5_5.9_marker_means.tsv")
  if (!nrow(mm)) stop("no marker means")
  nm <- names(mm)[1]
  long <- melt(mm, id.vars = nm, variable.name = "lineage", value.name = "mean")
  setnames(long, nm, "gene")
  long[, gene := gsub("^g_", "", gene)]
  p <- ggplot(long, aes(lineage, gene, fill = mean)) +
    geom_tile(color = "white") +
    scale_fill_gradient2(low = "#3C5488", mid = "white", high = "#E64B35", midpoint = median(long$mean, na.rm = TRUE)) +
    labs(title = "5.9  Canonical marker means by lineage", x = NULL, y = NULL, fill = "mean") +
    theme_sci() + theme(axis.line = element_blank())
  save_both(p, file.path(fig, "M5_5.9_marker_means.pdf"), 6.4, 6.8)
})

run("5.10 lineage UMAP", function() {
  if (!nrow(cells) || !"lineage" %in% names(cells)) stop("no lineage")
  p <- umap_base(cells, "lineage", "5.10  Lineage scores (T/NK/B/myeloid/Ery)", pal_lin)
  save_both(p, file.path(fig, "M5_5.10_umap_lineage.pdf"), 8.0, 6.4)
})

run("5.11 malignant UMAP", function() {
  if (!nrow(cells) || !"malignant" %in% names(cells)) stop("no mal")
  p <- umap_base(cells, "malignant", "5.11  Malignant call (T-high, patient-private clusters)", pal_mal)
  save_both(p, file.path(fig, "M5_5.11_umap_malignant.pdf"), 8.0, 6.4)
})

run("5.12 cnv proxy", function() {
  if (!nrow(cells) || !"cnv_proxy" %in% names(cells)) stop("no cnv")
  cells[, cnv_proxy := as.numeric(cnv_proxy)]
  p <- umap_base(cells, "cnv_proxy", "5.12  CNV proxy = T-score minus myeloid-score (inferCNV not installed)", discrete = FALSE)
  save_both(p, file.path(fig, "M5_5.12_umap_cnv_proxy.pdf"), 8.0, 6.4)
})

run("5.13 composition", function() {
  if (!nrow(comp)) stop("no comp")
  p <- ggplot(comp, aes(sample, frac, fill = lineage)) +
    geom_col(width = 0.88, color = "white", linewidth = 0.15) +
    scale_fill_manual(values = pal_lin) +
    labs(title = "5.13  Lineage composition by sample", x = NULL, y = "Fraction", fill = NULL) +
    theme_sci() + theme(axis.text.x = element_text(angle = 45, hjust = 1))
  save_both(p, file.path(fig, "M5_5.13_composition.pdf"), 9.4, 5.6)
})

gene_umap <- function(col, title, outfile) {
  if (!nrow(cells) || !col %in% names(cells)) stop(col)
  d <- copy(cells)
  d[, v := as.numeric(d[[col]])]
  p <- ggplot(d, aes(UMAP1, UMAP2, color = v)) +
    geom_point(size = 0.15, alpha = 0.75) +
    scale_color_gradientn(colors = c("#3C5488", "#EEEEEE", "#E64B35"), na.value = "grey90") +
    labs(title = title, x = "UMAP1", y = "UMAP2", color = NULL) +
    theme_sci() + theme(axis.text = element_blank(), axis.ticks = element_blank())
  save_both(p, outfile, 7.6, 6.2)
}

run("5.14 ZEB1 feature", function() gene_umap("g_ZEB1", "5.14  ZEB1", file.path(fig, "M5_5.14_umap_ZEB1.pdf")))
run("5.17 ZEB2 feature", function() gene_umap("g_ZEB2", "5.17  ZEB2", file.path(fig, "M5_5.17_umap_ZEB2.pdf")))
run("5.18 LMO2 feature", function() gene_umap("g_LMO2", "5.18  LMO2", file.path(fig, "M5_5.18_umap_LMO2.pdf")))
run("5.19 IL7R feature", function() gene_umap("g_IL7R", "5.19  IL7R", file.path(fig, "M5_5.19_umap_IL7R.pdf")))

run("5.15 ZEB1 by patient", function() {
  if (!nrow(cells) || !"g_ZEB1" %in% names(cells)) stop("no zeb1")
  d <- cells[malignant == "malignant"]
  if (!nrow(d)) d <- cells
  p <- ggplot(d, aes(patient, as.numeric(g_ZEB1), fill = timepoint)) +
    geom_boxplot(outlier.size = 0.15, width = 0.7, position = position_dodge(0.8)) +
    scale_fill_manual(values = pal_tp) +
    labs(title = "5.15  ZEB1 in malignant cells by patient",
         subtitle = "Boxes are cells; statistics use patient means",
         x = NULL, y = "ZEB1 log1p", fill = NULL) + theme_sci()
  save_both(p, file.path(fig, "M5_5.15_ZEB1_by_patient.pdf"), 8.6, 5.6)
})

run("5.16 malignant vs other", function() {
  if (!nrow(pat)) {
    # cell-level display only
    d <- cells[, .(mean_z = mean(as.numeric(g_ZEB1), na.rm = TRUE)), by = .(patient, malignant)]
  } else {
    # still show cells but subtitle warns
    d <- cells[, .(mean_z = mean(as.numeric(g_ZEB1), na.rm = TRUE), n = .N), by = .(patient, malignant)]
  }
  p <- ggplot(d, aes(malignant, mean_z)) +
    geom_boxplot(width = 0.5, outlier.shape = NA, fill = "#D9D9D9") +
    geom_point(aes(color = patient), size = 2.6, position = position_jitter(width = 0.08, seed = 1)) +
    labs(title = "5.16  Patient-mean ZEB1, malignant vs other",
         x = NULL, y = "Mean ZEB1", color = "patient") + theme_sci()
  save_both(p, file.path(fig, "M5_5.16_ZEB1_malignant_vs_other.pdf"), 7.2, 5.6)
})

run("5.20 GATA3 program", function() gene_umap("sig_GATA3prog", "5.20  GATA3/TCF7/BCL11B program (ZEB1 regulon substitute)", file.path(fig, "M5_5.20_umap_GATA3prog.pdf")))
run("5.21 ETP state", function() gene_umap("sig_ETP", "5.21  ETP / immature score (not assumed to be ZEB1-low)", file.path(fig, "M5_5.21_umap_ETP.pdf")))
run("5.22 cytotrace proxy", function() gene_umap("cytotrace_proxy", "5.22  CytoTRACE substitute (scaled nFeature in malignant)", file.path(fig, "M5_5.22_umap_cytotrace_proxy.pdf")))
run("5.24 DPT UMAP", function() gene_umap("dpt_pseudotime", "5.24  Diffusion pseudotime (Slingshot/Monocle3 not installed)", file.path(fig, "M5_5.24_umap_dpt.pdf")))

dpt_gene <- function(gcol, title, outfile) {
  if (!nrow(cells) || !"dpt_pseudotime" %in% names(cells) || !gcol %in% names(cells)) stop(gcol)
  d <- cells[malignant == "malignant" & is.finite(as.numeric(dpt_pseudotime))]
  if (!nrow(d)) d <- cells[is.finite(as.numeric(dpt_pseudotime))]
  d[, x := as.numeric(dpt_pseudotime)]
  d[, y := as.numeric(d[[gcol]])]
  p <- ggplot(d, aes(x, y)) +
    geom_point(size = 0.2, alpha = 0.25, color = "#8491B4") +
    geom_smooth(method = "loess", se = TRUE, color = "#E64B35", linewidth = 0.9, fill = "#F39B7F") +
    labs(title = title, x = "Diffusion pseudotime", y = gcol) + theme_sci()
  save_both(p, outfile, 7.2, 5.4)
}
run("5.25 ZEB1 DPT", function() dpt_gene("g_ZEB1", "5.25  ZEB1 along DPT (malignant cells)", file.path(fig, "M5_5.25_ZEB1_dpt.pdf")))
run("5.26 ZEB2 DPT", function() dpt_gene("g_ZEB2", "5.26  ZEB2 along DPT", file.path(fig, "M5_5.26_ZEB2_dpt.pdf")))
run("5.27 LMO2 DPT", function() dpt_gene("g_LMO2", "5.27  LMO2 along DPT", file.path(fig, "M5_5.27_LMO2_dpt.pdf")))
run("5.28 stemness DPT", function() dpt_gene("sig_Stemness", "5.28  Stemness score along DPT", file.path(fig, "M5_5.28_stemness_dpt.pdf")))

run("5.29 no thymocyte atlas", function() {
  p <- ggplot() + annotate("text", x = 0, y = 0, label = "5.29  No sorted thymocyte reference in the download inventory.\nGSE227122 is T-ALL only. Mapping skipped rather than mixing platforms.") +
    theme_void() + xlim(-1, 1) + ylim(-1, 1)
  save_both(p, file.path(fig, "M5_5.29_no_thymocyte_reference.pdf"), 8.0, 3.2)
})

score_umap <- function(col, title, outfile) gene_umap(col, title, outfile)
run("5.30 Tdiff", function() score_umap("sig_Tcell_diff", "5.30  T-cell differentiation score", file.path(fig, "M5_5.30_umap_Tdiff.pdf")))
run("5.31 stemness", function() score_umap("sig_Stemness", "5.31  Stemness score", file.path(fig, "M5_5.31_umap_stemness.pdf")))
run("5.32 prolif", function() score_umap("sig_Prolif", "5.32  Cell-cycle / proliferation score", file.path(fig, "M5_5.32_umap_prolif.pdf")))
run("5.33 JAK", function() score_umap("sig_JAKSTAT", "5.33  IL7/JAK/STAT score", file.path(fig, "M5_5.33_umap_JAKSTAT.pdf")))
run("5.34 NOTCH", function() score_umap("sig_NOTCH", "5.34  NOTCH score", file.path(fig, "M5_5.34_umap_NOTCH.pdf")))
run("5.35 TGFb", function() score_umap("sig_TGFb", "5.35  TGF-beta score", file.path(fig, "M5_5.35_umap_TGFb.pdf")))

run("5.36 patient heatmap", function() {
  if (!nrow(pat)) stop("no patient pb")
  cols <- intersect(c("g_ZEB1", "g_ZEB2", "g_LMO2", "g_IL7R", "sig_ETP", "sig_Stemness",
                      "sig_Tcell_diff", "sig_JAKSTAT", "sig_NOTCH", "sig_TGFb", "sig_Prolif"), names(pat))
  if (length(cols) < 3) stop("few columns")
  m <- as.matrix(pat[, ..cols])
  rownames(m) <- pat$patient
  m <- scale(m)
  m[!is.finite(m)] <- 0
  if (has_cht && has_circlize) {
    pdf(file.path(fig, "M5_5.36_patient_pseudobulk.pdf"), width = 8.4, height = 6.2)
    print(ComplexHeatmap::Heatmap(t(m), name = "z",
                                  column_title = "5.36  Dx malignant pseudobulk (patient n)",
                                  col = circlize::colorRamp2(c(-2, 0, 2), c("#3C5488", "white", "#E64B35"))))
    dev.off()
  } else {
    d <- data.table(patient = rownames(m)); d <- cbind(d, as.data.table(m))
    long <- melt(d, id.vars = "patient")
    p <- ggplot(long, aes(patient, variable, fill = value)) + geom_tile(color = "white") +
      scale_fill_gradient2(low = "#3C5488", mid = "white", high = "#E64B35") +
      labs(title = "5.36  Dx malignant pseudobulk", x = NULL, y = NULL) + theme_sci()
    save_both(p, file.path(fig, "M5_5.36_patient_pseudobulk.pdf"), 8.4, 6.2)
  }
})

run("5.37 paired Dx EOI", function() {
  if (!nrow(pairs) || !"g_ZEB1_Dx" %in% names(pairs)) stop("no pairs")
  long <- rbind(
    pairs[, .(patient, timepoint = "Dx", ZEB1 = g_ZEB1_Dx, ETP = sig_ETP_Dx, LMO2 = g_LMO2_Dx)],
    pairs[, .(patient, timepoint = "EOI", ZEB1 = g_ZEB1_EOI, ETP = sig_ETP_EOI, LMO2 = g_LMO2_EOI)]
  )
  p1 <- ggplot(long, aes(timepoint, ZEB1, group = patient, color = patient)) +
    geom_point(size = 2.8) + geom_line(linewidth = 0.5) +
    labs(title = "5.37  Paired Dx vs EOI ZEB1 (patient means)", y = "Mean ZEB1", x = NULL) + theme_sci()
  save_both(p1, file.path(fig, "M5_5.37_paired_ZEB1.pdf"), 6.6, 5.4)
})

run("5.38 residual blast", function() {
  if (!nrow(malc)) stop("no malc")
  d <- malc[malignant == "malignant"]
  p <- ggplot(d, aes(timepoint, frac, color = patient, group = patient)) +
    geom_point(size = 2.8) + geom_line(linewidth = 0.45) +
    labs(title = "5.38  Malignant-call fraction (EOI residual proxy)",
         y = "Malignant fraction", x = NULL) + theme_sci()
  save_both(p, file.path(fig, "M5_5.38_malignant_fraction.pdf"), 6.8, 5.4)
})

run("5.39 relapse", function() {
  rel <- rd("M5_5.39_relapse_counts.tsv")
  if (!nrow(rel)) stop("no relapse")
  p <- ggplot(rel, aes(lineage, n, fill = malignant)) +
    geom_col(position = "stack", width = 0.7, color = "white") +
    scale_fill_manual(values = pal_mal) +
    labs(title = "5.39  T11 relapse composition (n = 1 patient)", x = NULL, y = "Cells") + theme_sci()
  save_both(p, file.path(fig, "M5_5.39_relapse.pdf"), 6.8, 5.2)
})

run("5.40 ETP-high fraction", function() {
  if (!nrow(pat) || !"frac_ETP_high" %in% names(pat)) stop("no frac")
  setorder(pat, frac_ETP_high)
  pat[, patient := factor(patient, levels = patient)]
  p <- ggplot(pat, aes(frac_ETP_high, patient, fill = as.numeric(g_ZEB1))) +
    geom_col(width = 0.7, color = "white") +
    scale_fill_gradientn(colors = c("#3C5488", "#EEEEEE", "#E64B35")) +
    labs(title = "5.40  Fraction ETP-high among Dx malignant cells",
         subtitle = "Not a hard ZEB1-low class; ETP score > cohort median",
         x = "Fraction", y = NULL, fill = "mean ZEB1") + theme_sci()
  save_both(p, file.path(fig, "M5_5.40_ETP_high_fraction.pdf"), 7.2, 5.8)
})

run("5.41 patient forest", function() {
  if (!nrow(assoc)) stop("no assoc")
  d <- copy(assoc)
  setorder(d, rho)
  d[, feature := factor(feature, levels = feature)]
  p <- ggplot(d, aes(rho, feature, fill = rho > 0)) +
    geom_vline(xintercept = 0, linetype = 2) +
    geom_col(width = 0.7, color = "white") +
    geom_text(aes(label = sprintf("rho=%.2f P=%.3f n=%d", rho, p, n)),
              hjust = ifelse(d$rho >= 0, -0.05, 1.05), size = 2.7) +
    scale_fill_manual(values = c("TRUE" = "#E64B35", "FALSE" = "#3C5488"), guide = "none") +
    labs(title = "5.41  Patient-level Spearman vs mean ZEB1 (Dx malignant)",
         x = "rho", y = NULL) + theme_sci()
  save_both(p, file.path(fig, "M5_5.41_patient_forest.pdf"), 8.4, 6.6)
})

run("5.42 GSE248287", function() {
  if (nrow(v248)) {
    setorder(v248, rho)
    v248[, feature := factor(feature, levels = feature)]
    p <- ggplot(v248, aes(rho, feature, fill = rho > 0)) +
      geom_vline(xintercept = 0, linetype = 2) +
      geom_col(width = 0.7, color = "white") +
      scale_fill_manual(values = c("TRUE" = "#E64B35", "FALSE" = "#3C5488"), guide = "none") +
      labs(title = "5.42  GSE248287 CITE-seq patient-level vs ZEB1",
           x = "rho", y = NULL) + theme_sci()
    save_both(p, file.path(fig, "M5_5.42_GSE248287_patient_forest.pdf"), 8.0, 6.4)
  } else if (nrow(pb248)) {
    p <- ggplot(pb248, aes(g_ZEB1, sig_ETP)) +
      geom_point(size = 3, shape = 21, fill = "#4DBBD5") +
      geom_smooth(method = "lm", se = TRUE, color = "#E64B35", linewidth = 0.8) +
      labs(title = "5.42  GSE248287 patient means", x = "ZEB1", y = "ETP score") + theme_sci()
    save_both(p, file.path(fig, "M5_5.42_GSE248287_scatter.pdf"), 6.6, 5.8)
  } else stop("no 248287 table")
})

cat("MODULE5 ANALYZE DONE\n")
cat("FIG", length(list.files(fig, pattern = "\\.(pdf|png)$")), "\n")

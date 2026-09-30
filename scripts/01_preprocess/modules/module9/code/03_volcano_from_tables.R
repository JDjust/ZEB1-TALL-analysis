# Rebuild Module 9 volcanoes / concordance from already-written DE tables.
suppressPackageStartupMessages({ library(data.table); library(ggplot2) })
has_ggrepel <- requireNamespace("ggrepel", quietly = TRUE)
root <- "/data-b/liangfuhua/projects/ZEB1_TALL_analysis/module9"
fig <- file.path(root, "figures"); tab <- file.path(root, "tables")
highlight <- c("ZEB1","ZEB2","LMO2","TLX1","TLX3","TAL1","LYL1","GATA3","TCF7","BCL11B",
               "MEF2C","RUNX1","MYC","IL7R","CD34","NOTCH1","SPI1","RAG1","DNTT","LEF1")
theme_sci <- function(base_size = 11) {
  theme_classic(base_size = base_size) +
    theme(axis.text = element_text(color = "black"), axis.title = element_text(color = "black"),
          axis.line = element_line(linewidth = 0.45, color = "black"),
          plot.title = element_text(face = "bold", size = base_size + 1),
          plot.subtitle = element_text(size = base_size - 1, color = "grey25"),
          strip.background = element_blank())
}
save_both <- function(plot, file, width, height) {
  pdf(file, width = width, height = height, useDingbats = FALSE); print(plot); dev.off()
  png(sub("\\.pdf$", ".png", file), width = width * 180, height = height * 180, res = 180)
  print(plot); dev.off()
  message("saved ", basename(file))
}
volcano_dt <- function(dt, title, file, lfc = "log2FoldChange", padj = "padj", lfc_cut = 0.5, sub = NULL) {
  dv <- copy(dt)
  if (!"gene" %in% names(dv) && "Gene.Name" %in% names(dv)) dv[, gene := as.character(dv[["Gene.Name"]])]
  if (!"gene" %in% names(dv)) stop("no gene column: ", paste(names(dv), collapse = ","))
  lfc_name <- lfc; padj_name <- padj
  if (!lfc_name %in% names(dv)) stop("missing LFC column ", lfc_name)
  if (!padj_name %in% names(dv)) stop("missing padj column ", padj_name)
  lfc_vals <- as.numeric(dv[[lfc_name]])
  padj_vals <- as.numeric(dv[[padj_name]])
  dv[, lfc_plot := lfc_vals]
  dv[, ypadj := padj_vals]
  dv <- dv[is.finite(lfc_plot) & is.finite(ypadj) & ypadj >= 0]
  if (!nrow(dv)) {
    # limma n=2 may leave padj NA; fall back to pvalue
    if ("pvalue" %in% names(dt)) {
      dv <- copy(dt)
      if (!"gene" %in% names(dv) && "Gene.Name" %in% names(dv)) dv[, gene := as.character(dv[["Gene.Name"]])]
      lfc_vals <- as.numeric(dv[[lfc_name]])
      padj_vals <- as.numeric(dv[["pvalue"]])
      dv[, lfc_plot := lfc_vals]
      dv[, ypadj := padj_vals]
      dv <- dv[is.finite(lfc_plot) & is.finite(ypadj) & ypadj >= 0]
      sub <- paste(sub, "; y-axis uses nominal P (FDR all NA)")
    }
  }
  if (!nrow(dv)) stop("volcano: no finite LFC/P rows")
  dv[, mlogp := -log10(pmin(pmax(ypadj, 1e-300), 1))]
  dv[, sig_class := fifelse(ypadj < 0.05 & lfc_plot > lfc_cut, "up",
                      fifelse(ypadj < 0.05 & lfc_plot < -lfc_cut, "down", "NS"))]
  dv[, gene_up := toupper(gsub("-", "", as.character(gene)))]
  hit <- dv[gene_up %in% highlight]
  p <- ggplot(as.data.frame(dv), aes(lfc_plot, mlogp, color = sig_class)) +
    geom_point(size = 0.35, alpha = 0.4) +
    geom_point(data = as.data.frame(hit), size = 2.3, shape = 21, fill = "white", stroke = 0.8, color = "black") +
    geom_vline(xintercept = c(-lfc_cut, lfc_cut), linetype = 2, linewidth = 0.35, color = "grey50") +
    geom_hline(yintercept = -log10(0.05), linetype = 2, linewidth = 0.35, color = "grey50") +
    scale_color_manual(values = c(up = "#E64B35", down = "#3C5488", NS = "grey75"), drop = FALSE) +
    labs(title = title, subtitle = sub, x = "log2 fold change", y = "-log10 P/FDR", color = NULL) +
    theme_sci() + theme(legend.position = "top")
  if (has_ggrepel && nrow(hit))
    p <- p + ggrepel::geom_text_repel(data = as.data.frame(hit), aes(label = gene_up),
                                      size = 3, color = "black", max.overlaps = 40)
  save_both(p, file, 7.8, 6.2)
}
concord_scatter <- function(a, b, lab_a, lab_b, file, title) {
  aa <- copy(a)[, gene := toupper(as.character(gene))]
  bb <- copy(b)[, gene := toupper(as.character(gene))]
  m <- merge(aa[, .(gene, log2FoldChange, padj)], bb[, .(gene, log2FoldChange, padj)], by = "gene")
  setnames(m, c("gene", "lfc_a", "padj_a", "lfc_b", "padj_b"))
  m <- m[is.finite(lfc_a) & is.finite(lfc_b)]
  rho <- suppressWarnings(cor(m$lfc_a, m$lfc_b, method = "spearman"))
  hit <- m[gene %in% highlight]
  p <- ggplot(as.data.frame(m), aes(lfc_a, lfc_b)) +
    geom_hline(yintercept = 0, linetype = 2, color = "grey70") +
    geom_vline(xintercept = 0, linetype = 2, color = "grey70") +
    geom_abline(slope = 1, intercept = 0, linetype = 3, color = "grey50") +
    geom_point(size = 0.3, alpha = 0.25, color = "grey50") +
    geom_point(data = as.data.frame(hit), shape = 21, size = 2.4, fill = "white", color = "black") +
    labs(title = title, subtitle = sprintf("Spearman rho = %.2f; n = %d genes", rho, nrow(m)),
         x = lab_a, y = lab_b) + theme_sci()
  if (has_ggrepel && nrow(hit))
    p <- p + ggrepel::geom_text_repel(data = as.data.frame(hit), aes(label = gene), size = 3, max.overlaps = 40)
  save_both(p, file, 7.2, 6.4)
  invisible(rho)
}

specs <- list(
  list(f = "M9_gse110635_combined_KD.tsv", title = "TLX1 KD vs control (total RNA, combined siRNAs)",
       out = "M9_9.3_tlx1_kd_total_volcano.pdf", sub = "GSE110635 ALL-SIL; DESeq2; n = 6 KD vs 3 NTC"),
  list(f = "M9_gse110632_combined_KD.tsv", title = "TLX1 KD vs control (polyA, combined siRNAs)",
       out = "M9_9.3b_tlx1_kd_polya_volcano.pdf", sub = "GSE110632 ALL-SIL; DESeq2"),
  list(f = "M9_gse62143_TLX1_OE.tsv", title = "TLX1 overexpression vs empty (CD34+ thymocytes, Agilent)",
       out = "M9_9.6_tlx1_oe_volcano.pdf", sub = "GSE62143; limma-trend; n=2 exploratory"),
  list(f = "M9_gse62141_KD.tsv", title = "TLX1 KD vs scrambled (GSE62141 Agilent, 24h)",
       out = "M9_9.3c_gse62141_volcano.pdf", sub = "Independent ALL-SIL KD; limma-trend"),
  list(f = "M9_gse110636_TLX_vs_rest.tsv", title = "Primary T-ALL: TLX vs rest (GSE110636)",
       out = "M9_9.19_patient_TLX_volcano.pdf", sub = "DESeq2"),
  list(f = "M9_gse186943_TG_vs_WT.tsv", title = "Mouse LMO2 transgenic vs WT",
       out = "M9_9.16_lmo2_oe_volcano.pdf", sub = "GSE186943 FPKM; limma-trend"),
  list(f = "M9_gse144035_precomputed.tsv", title = "Reversible Lmo2 model (author DE)",
       out = "M9_9.17_lmo2_OFF_volcano.pdf", sub = "GSE144035")
)
for (s in specs) {
  fp <- file.path(tab, s$f)
  if (!file.exists(fp)) { message("skip ", s$f); next }
  d <- fread(fp)
  message("TRY ", s$f, " nrow=", nrow(d), " names=", paste(names(d), collapse = "|"))
  tryCatch(volcano_dt(d, s$title, file.path(fig, s$out), sub = s$sub),
           error = function(e) message("FAIL ", s$f, ": ", conditionMessage(e),
                                       " call=", paste(deparse(conditionCall(e)), collapse=" ")))
}

a <- fread(file.path(tab, "M9_gse110635_combined_KD.tsv"))
b <- fread(file.path(tab, "M9_gse110632_combined_KD.tsv"))
tryCatch(concord_scatter(a, b, "total RNA log2FC", "polyA log2FC",
                         file.path(fig, "M9_9.5_total_polya_concordance.pdf"),
                         "TLX1 KD: total RNA vs polyA"),
         error = function(e) message("FAIL concord: ", conditionMessage(e)))
s1 <- fread(file.path(tab, "M9_gse110635_siTLX1_1.tsv"))
s2 <- fread(file.path(tab, "M9_gse110635_siTLX1_2.tsv"))
tryCatch(concord_scatter(s1, s2, "siTLX1-1 log2FC", "siTLX1-2 log2FC",
                         file.path(fig, "M9_9.4_two_sirna_concordance.pdf"),
                         "Two independent TLX1 siRNAs (total RNA)"),
         error = function(e) message("FAIL sirna: ", conditionMessage(e)))
oe <- fread(file.path(tab, "M9_gse62143_TLX1_OE.tsv"))
tryCatch(concord_scatter(a, oe, "TLX1 KD log2FC", "TLX1 OE log2FC",
                         file.path(fig, "M9_9.7_kd_oe_inverse.pdf"),
                         "KD vs OE concordance"),
         error = function(e) message("FAIL kdoe: ", conditionMessage(e)))
pt <- fread(file.path(tab, "M9_gse110636_TLX_vs_rest.tsv"))
tryCatch(concord_scatter(pt, a, "patient TLX vs rest log2FC", "TLX1 KD log2FC",
                         file.path(fig, "M9_9.11_patient_vs_kd_scatter.pdf"),
                         "Patient TLX association vs TLX1 KD"),
         error = function(e) message("FAIL patient: ", conditionMessage(e)))
message("VOLCANO SCRIPT DONE")
message("ZEB1 110635: ", paste(a[toupper(gene)=="ZEB1", .(log2FoldChange, padj, pvalue)], collapse=" "))
message("TLX1 110635: ", paste(a[toupper(gene)=="TLX1", .(log2FoldChange, padj, pvalue)], collapse=" "))

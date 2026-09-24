# Module 4 revision: subtype-adjusted ZEB1 transcriptome (limma)
# β1 = association of each gene with ZEB1 after subtype (+ age).
# ZEB1 stays continuous. Do not treat Q1/Q4 as biology.

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(limma)
})
has_fgsea <- requireNamespace("fgsea", quietly = TRUE)
has_ggrepel <- requireNamespace("ggrepel", quietly = TRUE)

root <- "/data-b/liangfuhua/projects/ZEB1_TALL_analysis/module4"
m3 <- "/data-b/liangfuhua/projects/ZEB1_TALL_analysis/module3/processed"
fig <- file.path(root, "figures"); tab <- file.path(root, "tables"); proc <- file.path(root, "processed")
gmt_dir <- "/data-b/liangfuhua/projects/TALL_dataset_download/data/MSigDB"
dir.create(fig, FALSE, TRUE); dir.create(tab, FALSE, TRUE)

theme_sci <- function(base_size = 11) {
  theme_classic(base_size = base_size) +
    theme(axis.text = element_text(color = "black"),
          axis.title = element_text(color = "black"),
          axis.line = element_line(linewidth = 0.45, color = "black"),
          plot.title = element_text(face = "bold", size = base_size + 1),
          plot.subtitle = element_text(size = base_size - 1, color = "grey25"),
          strip.background = element_blank())
}
save_both <- function(plot, file, width, height) {
  pdf(file, width = width, height = height, useDingbats = FALSE); print(plot); dev.off()
  png(sub("\\.pdf$", ".png", file), width = width * 180, height = height * 180, res = 180)
  print(plot); dev.off()
}

read_mat <- function(path) {
  d <- as.data.frame(fread(path), stringsAsFactors = FALSE)
  rownames(d) <- d[[1]]; d[[1]] <- NULL
  as.matrix(d)
}

expr <- read_mat(file.path(proc, "target_log2tpm_filtered.tsv.gz"))
lab <- fread(file.path(m3, "target_tall_subtype.tsv"))
lab <- lab[usi %in% colnames(expr)]
lab <- lab[subtype != "Unknown" & !is.na(subtype) & subtype != ""]
expr <- expr[, lab$usi, drop = FALSE]
lab[, ZEB1_z := as.numeric(scale(log2(ZEB1 + 1)))]
lab[, subtype := factor(subtype)]
if ("age_years" %in% names(lab)) lab[, age_years := as.numeric(age_years)]

ok <- is.finite(lab$ZEB1_z)
if ("age_years" %in% names(lab) && sum(is.finite(lab$age_years)) >= 50) {
  ok <- ok & is.finite(lab$age_years)
  form <- ~ ZEB1_z + subtype + age_years
  form0 <- ~ ZEB1_z
} else {
  form <- ~ ZEB1_z + subtype
  form0 <- ~ ZEB1_z
}
lab <- lab[ok]
expr <- expr[, lab$usi, drop = FALSE]
message("n samples ", ncol(expr), " genes ", nrow(expr), " formula ", deparse(form))
message(paste(capture.output(table(lab$subtype)), collapse = "\n"))

design0 <- model.matrix(form0, data = lab)
design <- model.matrix(form, data = lab)
stopifnot(nrow(design) == ncol(expr))

fit0 <- eBayes(lmFit(expr, design0), trend = TRUE, robust = TRUE)
fit  <- eBayes(lmFit(expr, design),  trend = TRUE, robust = TRUE)

tt0 <- as.data.table(topTable(fit0, coef = "ZEB1_z", number = Inf, sort.by = "none"), keep.rownames = "gene")
tt  <- as.data.table(topTable(fit,  coef = "ZEB1_z", number = Inf, sort.by = "none"), keep.rownames = "gene")
setnames(tt0, c("logFC", "P.Value", "adj.P.Val", "t"),
         c("logFC_unadj", "p_unadj", "fdr_unadj", "t_unadj"), skip_absent = TRUE)
setnames(tt,  c("logFC", "P.Value", "adj.P.Val", "t"),
         c("logFC_adj", "p_adj", "fdr_adj", "t_adj"), skip_absent = TRUE)
mrg <- merge(tt0[, .(gene, logFC_unadj, t_unadj, p_unadj, fdr_unadj)],
             tt[,  .(gene, logFC_adj, t_adj, p_adj, fdr_adj)], by = "gene")
fwrite(mrg, file.path(tab, "M4_r2_subtype_adjusted_ZEB1.tsv"), sep = "\t")

n_un <- mrg[fdr_unadj < 0.05, .N]
n_ad <- mrg[fdr_adj < 0.05, .N]
n_both <- mrg[fdr_unadj < 0.05 & fdr_adj < 0.05, .N]
n_lost <- mrg[fdr_unadj < 0.05 & fdr_adj >= 0.05, .N]
n_gain <- mrg[fdr_unadj >= 0.05 & fdr_adj < 0.05, .N]
sumdt <- data.table(
  n_samples = ncol(expr), n_genes = nrow(mrg),
  n_FDR05_unadj = n_un, n_FDR05_adj = n_ad,
  n_both = n_both, n_lost_after_subtype = n_lost, n_gain_after_subtype = n_gain,
  spearman_t_unadj_vs_adj = cor(mrg$t_unadj, mrg$t_adj, method = "spearman", use = "complete")
)
fwrite(sumdt, file.path(tab, "M4_r2_adjust_summary.tsv"), sep = "\t")
print(sumdt)

keys <- c("LMO2", "LYL1", "ZEB2", "GATA3", "TCF7", "BCL11B", "TLX1", "TLX3",
          "TAL1", "MEF2C", "MYC", "IL7R", "CD34", "NOTCH1", "RUNX1")
kg <- mrg[gene %in% keys]
fwrite(kg, file.path(tab, "M4_r2_keygenes_adjusted.tsv"), sep = "\t")

# volcano adjusted
mrg[, sig := fifelse(fdr_adj < 0.05 & logFC_adj > 0, "up",
              fifelse(fdr_adj < 0.05 & logFC_adj < 0, "down", "NS"))]
labg <- mrg[gene %in% keys]
p_v <- ggplot(mrg, aes(logFC_adj, -log10(pmax(p_adj, 1e-300)))) +
  geom_point(aes(color = sig), size = 0.35, alpha = 0.55) +
  geom_hline(yintercept = -log10(0.05), linetype = 2, color = "grey50", linewidth = 0.35) +
  geom_vline(xintercept = 0, linetype = 2, color = "grey50", linewidth = 0.35) +
  scale_color_manual(values = c(up = "#E64B35", down = "#4DBBD5", NS = "grey75")) +
  labs(title = "TARGET: ZEB1 effect after subtype (+ age)",
       subtitle = sprintf("limma trend+robust; FDR<0.05: %d / %d genes (unadjusted %d)",
                          n_ad, nrow(mrg), n_un),
       x = "logFC per SD of ZEB1 (subtype-adjusted)", y = "-log10 P", color = NULL) +
  theme_sci() + theme(legend.position = "top")
if (has_ggrepel && nrow(labg)) {
  p_v <- p_v + ggrepel::geom_text_repel(data = labg, aes(label = gene), size = 2.7, max.overlaps = 40)
} else if (nrow(labg)) {
  p_v <- p_v + geom_text(data = labg, aes(label = gene), size = 2.5, vjust = -0.4)
}
save_both(p_v, file.path(fig, "M4_r2_4.1_adjusted_volcano.pdf"), 7.6, 6.4)

# unadj vs adj scatter
p_s <- ggplot(mrg, aes(t_unadj, t_adj)) +
  geom_point(size = 0.3, alpha = 0.35, color = "grey40") +
  geom_abline(slope = 1, intercept = 0, linetype = 2, color = "grey30") +
  geom_hline(yintercept = 0, color = "grey80") + geom_vline(xintercept = 0, color = "grey80") +
  labs(title = "Unadjusted vs subtype-adjusted ZEB1 t-statistic",
       subtitle = sprintf("Spearman rho = %.2f; lost FDR %d; gained %d",
                          sumdt$spearman_t_unadj_vs_adj, n_lost, n_gain),
       x = "t (ZEB1 only)", y = "t (ZEB1 + subtype + age)") +
  theme_sci()
if (nrow(labg)) {
  p_s <- p_s + geom_point(data = labg, color = "#E64B35", size = 1.6) +
    geom_text(data = labg, aes(label = gene), size = 2.6, vjust = -0.5)
}
save_both(p_s, file.path(fig, "M4_r2_4.1_unadj_vs_adj_scatter.pdf"), 6.8, 6.4)

# key gene forest of adjusted logFC
kg2 <- copy(kg)
setorder(kg2, logFC_adj)
kg2[, gene := factor(gene, levels = gene)]
p_f <- ggplot(kg2, aes(logFC_adj, gene, fill = logFC_adj < 0)) +
  geom_vline(xintercept = 0, linetype = 2, color = "grey40") +
  geom_segment(aes(x = 0, xend = logFC_adj, yend = gene), color = "#B07D3A", linewidth = 0.9) +
  geom_point(shape = 21, size = 3.4, color = "black") +
  geom_text(aes(label = sprintf("FDR=%s", ifelse(fdr_adj < 0.001, "<0.001", sprintf("%.3f", fdr_adj)))),
            hjust = ifelse(kg2$logFC_adj >= 0, -0.1, 1.1), size = 2.6) +
  scale_fill_manual(values = c("TRUE" = "#3C5488", "FALSE" = "#E64B35"), guide = "none") +
  labs(title = "Pre-specified genes, subtype-adjusted ZEB1 effect",
       x = "logFC per SD ZEB1", y = NULL) + theme_sci()
save_both(p_f, file.path(fig, "M4_r2_4.1_keygene_forest.pdf"), 7.4, 6.2)

# GSEA on adjusted t
if (has_fgsea) {
  ranks <- mrg[!is.na(t_adj), setNames(t_adj, gene)]
  ranks <- sort(ranks[!duplicated(names(ranks))], decreasing = TRUE)
  gmt <- file.path(gmt_dir, "h.all.v2023.2.Hs.symbols.gmt")
  if (file.exists(gmt)) {
    gs <- fgsea::gmtPathways(gmt)
    res <- as.data.table(fgsea::fgsea(gs, ranks, minSize = 10, maxSize = 500, nPermSimple = 4000))
    res[, leadingEdge := vapply(leadingEdge, function(x) paste(head(x, 10), collapse = ","), "")]
    res[, absNES := abs(NES)]
    setorder(res, padj, -absNES)
    fwrite(res, file.path(tab, "M4_r2_4.6_hallmark_adjusted.tsv"), sep = "\t")
    d <- res[!is.na(NES)][1:min(18, nrow(res))]
    d[, pathway := factor(pathway, levels = rev(pathway))]
    p_g <- ggplot(d, aes(NES, pathway)) +
      geom_point(aes(size = size, fill = -log10(pmax(padj, 1e-8))), shape = 21, color = "grey20") +
      geom_vline(xintercept = 0, linetype = 2) +
      scale_fill_gradient(low = "#4DBBD5", high = "#E64B35") +
      labs(title = "Hallmark GSEA on subtype-adjusted ZEB1 t",
           x = "NES", y = NULL, fill = "-log10 FDR") + theme_sci()
    save_both(p_g, file.path(fig, "M4_r2_4.6_hallmark_adj_bubble.pdf"), 8.4, 6.8)
  }
}

# Pharmacotype scheme B: ZEB1 + ETP + age
pexpr <- read_mat(file.path(proc, "pharmacotype_log2fpkm_filtered.tsv.gz"))
plab <- fread(file.path(m3, "pharmacotype_tall_subtype.tsv"))
idcol <- if ("sample" %in% names(plab)) "sample" else if ("k_patient" %in% names(plab)) "k_patient" else names(plab)[1]
if (!idcol %in% names(plab)) {
  message("Pharmacotype labels have no sample id; skip scheme B")
} else {
  # align columns
  common <- intersect(colnames(pexpr), as.character(plab[[idcol]]))
  if (length(common) < 30 && "usi" %in% names(plab)) {
    common <- intersect(colnames(pexpr), as.character(plab$usi))
    idcol <- "usi"
  }
  message("Pharmacotype overlap ", length(common), " idcol ", idcol)
  if (length(common) >= 30) {
    plab <- plab[as.character(get(idcol)) %in% common]
    pexpr <- pexpr[, as.character(plab[[idcol]]), drop = FALSE]
    plab[, ZEB1_z := as.numeric(scale(log2(ZEB1 + 1)))]
    if ("etp" %in% names(plab)) {
      plab[, etp := factor(etp)]
      formp <- if ("age_years" %in% names(plab) && sum(is.finite(plab$age_years)) >= 30)
        ~ ZEB1_z + etp + age_years else ~ ZEB1_z + etp
      okp <- is.finite(plab$ZEB1_z)
      if (grepl("age", deparse(formp))) okp <- okp & is.finite(plab$age_years)
      plab <- plab[okp]; pexpr <- pexpr[, as.character(plab[[idcol]]), drop = FALSE]
      des <- model.matrix(formp, data = plab)
      fitp <- eBayes(lmFit(pexpr, des), trend = TRUE, robust = TRUE)
      ttp <- as.data.table(topTable(fitp, coef = "ZEB1_z", number = Inf, sort.by = "none"),
                           keep.rownames = "gene")
      fwrite(ttp, file.path(tab, "M4_r2_pharmacotype_ETP_adjusted.tsv"), sep = "\t")
      message("Pharmacotype FDR<0.05 ", ttp[adj.P.Val < 0.05, .N], " / ", nrow(ttp))
    }
  }
}

writeLines(capture.output(sessionInfo()), file.path(tab, "M4_r2_sessionInfo.txt"))
message("M4_R2_DONE")

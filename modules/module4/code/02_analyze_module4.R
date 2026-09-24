# Module 4 figures from code2: volcano 206, scatter 193, heatmap 210,
# signed bars 142/346, forest 88. Not violin/lollipop.

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})
has_fgsea <- requireNamespace("fgsea", quietly = TRUE)
has_ggrepel <- requireNamespace("ggrepel", quietly = TRUE)
has_cht <- requireNamespace("ComplexHeatmap", quietly = TRUE)
has_circlize <- requireNamespace("circlize", quietly = TRUE)

root <- "/data-b/liangfuhua/projects/ZEB1_TALL_analysis/module4"
fig <- file.path(root, "figures"); tab <- file.path(root, "tables"); proc <- file.path(root, "processed")
dir.create(fig, FALSE, TRUE); dir.create(tab, FALSE, TRUE)
gmt_dir <- "/data-b/liangfuhua/projects/TALL_dataset_download/data/MSigDB"

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
    png(sub("\\.pdf$", ".png", file), width = width * 180, height = height * 180, res = 180)
    print(plot); dev.off()
    message("saved ", basename(file))
  }, error = function(e) message("FAIL ", basename(file), ": ", conditionMessage(e)))
}
run <- function(label, fn) tryCatch({ fn(); message("OK ", label) },
                                    error = function(e) message("FAIL ", label, ": ", conditionMessage(e)))
fmt_p <- function(p) ifelse(is.na(p), "NA", ifelse(p < 0.001, "P < 0.001", sprintf("P = %.3f", p)))

st <- fread(file.path(proc, "target_zeb1_gene_stats.tsv"))
ph <- fread(file.path(proc, "pharmacotype_zeb1_gene_stats.tsv"))
both <- fread(file.path(proc, "cross_cohort_gene_stats.tsv"))
lab <- fread(file.path(proc, "target_samples.tsv"))
highlight <- c("ZEB1", "ZEB2", "LMO2", "IL7R", "CD34", "CD1A", "TAL1", "LYL1", "NOTCH1",
               "MYC", "SOX4", "BCL11B", "TCF7", "GATA3", "MEF2C", "SPI1", "RAG1", "DNTT")

run("volcano TARGET 206", function() {
  dv <- copy(st)
  dv[, cls := fcase(fdr_wilcox < 0.05 & lfc_Q4_minus_Q1 > 0.5, "higher in ZEB1-Q4",
                    fdr_wilcox < 0.05 & lfc_Q4_minus_Q1 < -0.5, "higher in ZEB1-Q1",
                    default = "NS")]
  hit <- dv[gene %in% highlight]
  p <- ggplot(dv, aes(lfc_Q4_minus_Q1, -log10(pmax(fdr_wilcox, 1e-300)), color = cls)) +
    geom_point(size = 0.35, alpha = 0.4) +
    geom_point(data = hit, size = 2.3, shape = 21, fill = "white", stroke = 0.8, color = "black") +
    geom_vline(xintercept = c(-0.5, 0.5), linetype = 2, linewidth = 0.35, color = "grey50") +
    geom_hline(yintercept = -log10(0.05), linetype = 2, linewidth = 0.35, color = "grey50") +
    scale_color_manual(values = c("higher in ZEB1-Q4" = "#E64B35", "higher in ZEB1-Q1" = "#3C5488", NS = "grey75")) +
    labs(title = "TARGET Q4 vs Q1 by ZEB1 (code2 206 volcano)",
         subtitle = sprintf("n = %d T-ALL; ZEB1 continuous quartiles are display-only", unique(st$n)[1]),
         x = "log2 mean difference (Q4 - Q1)", y = "-log10 FDR", color = NULL) +
    theme_sci() + theme(legend.position = "top")
  if (has_ggrepel) p <- p + ggrepel::geom_text_repel(data = hit, aes(label = gene), size = 3, color = "black", max.overlaps = 30)
  save_both(p, file.path(fig, "M4_4.3_volcano_Q4_vs_Q1.pdf"), 7.8, 6.2)
})

run("rho volcano", function() {
  dv <- copy(st)
  dv[, cls := fcase(fdr_rho < 0.05 & rho_zeb1 > 0.3, "pos rho>0.3",
                    fdr_rho < 0.05 & rho_zeb1 < -0.3, "neg rho<-0.3",
                    default = "NS/weak")]
  hit <- dv[gene %in% highlight]
  p <- ggplot(dv, aes(rho_zeb1, -log10(pmax(fdr_rho, 1e-300)), color = cls)) +
    geom_point(size = 0.35, alpha = 0.4) +
    geom_point(data = hit, size = 2.3, shape = 21, fill = "white", stroke = 0.8, color = "black") +
    geom_vline(xintercept = c(-0.3, 0.3), linetype = 2, linewidth = 0.35, color = "grey50") +
    scale_color_manual(values = c("pos rho>0.3" = "#E64B35", "neg rho<-0.3" = "#3C5488", "NS/weak" = "grey75")) +
    labs(title = "TARGET genome-wide Spearman vs ZEB1 (code2 206)",
         x = "Spearman rho vs ZEB1", y = "-log10 FDR", color = NULL) +
    theme_sci() + theme(legend.position = "top")
  if (has_ggrepel) p <- p + ggrepel::geom_text_repel(data = hit, aes(label = gene), size = 3, color = "black", max.overlaps = 30)
  save_both(p, file.path(fig, "M4_4.1_volcano_rho.pdf"), 7.8, 6.2)
})

run("cross-cohort scatter 193", function() {
  d <- copy(both)
  d[, cls := fcase(fdr_rho_TARGET < 0.05 & fdr_rho_Pharmacotype < 0.05 &
                     sign(rho_zeb1_TARGET) == sign(rho_zeb1_Pharmacotype), "shared FDR<0.05",
                   default = "other")]
  hit <- d[gene %in% highlight]
  p <- ggplot(d, aes(rho_zeb1_TARGET, rho_zeb1_Pharmacotype, color = cls)) +
    geom_hline(yintercept = 0, linetype = 2, color = "grey60", linewidth = 0.3) +
    geom_vline(xintercept = 0, linetype = 2, color = "grey60", linewidth = 0.3) +
    geom_point(size = 0.45, alpha = 0.35) +
    geom_point(data = hit, size = 2.4, shape = 21, fill = "white", color = "black", stroke = 0.8) +
    scale_color_manual(values = c("shared FDR<0.05" = "#E64B35", other = "grey75")) +
    labs(title = "ZEB1 correlation replicated in Pharmacotype (code2 193)",
         subtitle = sprintf("shared same-direction FDR<0.05: %d / %d genes",
                            sum(d$cls == "shared FDR<0.05"), nrow(d)),
         x = "rho vs ZEB1 (TARGET)", y = "rho vs ZEB1 (Pharmacotype)", color = NULL) +
    theme_sci() + theme(legend.position = "top")
  if (has_ggrepel) p <- p + ggrepel::geom_text_repel(data = hit, aes(label = gene), size = 3, color = "black", max.overlaps = 25)
  save_both(p, file.path(fig, "M4_4.5_cross_cohort_rho.pdf"), 7.6, 6.4)
})

run("signed top rho 346", function() {
  top <- rbind(
    st[fdr_rho < 0.05][order(-rho_zeb1)][1:15][, side := "positive"],
    st[fdr_rho < 0.05][order(rho_zeb1)][1:15][, side := "negative"]
  )
  setorder(top, rho_zeb1)
  top[, gene := factor(gene, levels = gene)]
  p <- ggplot(top, aes(rho_zeb1, gene, fill = side)) +
    geom_vline(xintercept = 0, linetype = 2, color = "grey40") +
    geom_col(width = 0.75, color = "white") +
    scale_fill_manual(values = c(positive = "#E64B35", negative = "#3C5488")) +
    labs(title = "Strongest TARGET Spearman vs ZEB1 (code2 346)",
         x = "rho", y = NULL, fill = NULL) +
    theme_sci() + theme(legend.position = "top")
  save_both(p, file.path(fig, "M4_4.1_top_rho_bars.pdf"), 7.4, 7.2)
})

run("key gene forest 88", function() {
  k <- st[gene %in% highlight]
  setorder(k, rho_zeb1)
  k[, gene := factor(gene, levels = gene)]
  p <- ggplot(k, aes(rho_zeb1, gene)) +
    geom_vline(xintercept = 0, linetype = 2, color = "grey40") +
    geom_point(aes(fill = rho_zeb1 > 0), shape = 22, size = 3.6, color = "black") +
    geom_text(aes(label = sprintf("rho=%.2f  FDR=%s", rho_zeb1,
                                  ifelse(fdr_rho < 0.001, "<0.001", sprintf("%.3f", fdr_rho)))),
              hjust = ifelse(k$rho_zeb1 >= 0, -0.12, 1.12), size = 2.8) +
    scale_fill_manual(values = c("TRUE" = "#E64B35", "FALSE" = "#3C5488"), guide = "none") +
    labs(title = "Pre-specified genes vs ZEB1 (code2 88)",
         subtitle = "TARGET STAR TPM, n = 265", x = "Spearman rho", y = NULL) +
    theme_sci()
  save_both(p, file.path(fig, "M4_4.1_keygene_forest.pdf"), 8.2, 6.0)
})

run("heatmap top 210", function() {
  if (!has_cht || !has_circlize) stop("ComplexHeatmap missing")
  topg <- unique(c(
    st[fdr_rho < 0.05][order(-abs(rho_zeb1)), gene][1:40],
    intersect(highlight, st$gene)
  ))
  matf <- file.path(proc, "target_log2tpm_filtered.tsv.gz")
  expr <- fread(matf)
  gn <- names(expr)[1]
  expr <- expr[get(gn) %in% topg]
  m <- as.matrix(expr[, -1, with = FALSE])
  rownames(m) <- expr[[gn]]
  ord <- lab[order(ZEB1_log), usi]
  ord <- intersect(ord, colnames(m))
  m <- m[, ord, drop = FALSE]
  m <- t(scale(t(m)))
  m[!is.finite(m)] <- 0
  ann <- lab[match(ord, usi)]
  ha <- ComplexHeatmap::HeatmapAnnotation(
    ZEB1 = ann$ZEB1_log,
    subtype = ann$subtype,
    col = list(
      ZEB1 = circlize::colorRamp2(range(ann$ZEB1_log, na.rm = TRUE), c("#3C5488", "#E64B35")),
      subtype = c(TAL = "#E64B35", TLX = "#4DBBD5", HOXA = "#00A087",
                  "LMO2/LYL1" = "#3C5488", NKX2_1 = "#F39B7F", Unknown = "#8491B4")
    ),
    show_annotation_name = TRUE
  )
  pdf(file.path(fig, "M4_4.4_top_corr_heatmap.pdf"), width = 10, height = 8.2)
  print(ComplexHeatmap::Heatmap(
    m, name = "row z", top_annotation = ha, cluster_columns = FALSE,
    show_column_names = FALSE, column_title = "Samples ordered by ZEB1 (code2 210)",
    col = circlize::colorRamp2(c(-2, 0, 2), c("#3C5488", "white", "#E64B35")),
    row_names_gp = grid::gpar(fontsize = 7)
  ))
  dev.off()
})

run("fgsea 4.6-4.8", function() {
  if (!has_fgsea) stop("fgsea missing")
  gmt_files <- list.files(gmt_dir, pattern = "\\.gmt$", full.names = TRUE)
  if (!length(gmt_files)) stop("no gmt in MSigDB dir")
  ranks <- st[!is.na(rho_zeb1), setNames(rho_zeb1, gene)]
  ranks <- ranks[!duplicated(names(ranks))]
  ranks <- sort(ranks, decreasing = TRUE)
  all_res <- list()
  for (f in gmt_files) {
    gs <- fgsea::gmtPathways(f)
    res <- fgsea::fgsea(pathways = gs, stats = ranks, minSize = 10, maxSize = 500, nPermSimple = 5000)
    res <- as.data.table(res)
    res[, leadingEdge := vapply(leadingEdge, function(x) paste(head(x, 15), collapse = ","), character(1))]
    res[, collection := basename(f)]
    all_res[[basename(f)]] <- res
    fwrite(res[order(padj, -abs(NES))], file.path(tab, paste0("M4_fgsea_", gsub("[^A-Za-z0-9]+", "_", basename(f)), ".tsv")), sep = "\t")
    top <- rbind(res[NES > 0][order(padj)][1:12], res[NES < 0][order(padj)][1:12])
    top <- top[!is.na(pathway)]
    setorder(top, NES)
    top[, pathway := factor(pathway, levels = pathway)]
    p <- ggplot(top, aes(NES, pathway, fill = NES > 0)) +
      geom_vline(xintercept = 0, linetype = 2, color = "grey40") +
      geom_col(width = 0.75, color = "white") +
      scale_fill_manual(values = c("TRUE" = "#E64B35", "FALSE" = "#3C5488"), guide = "none") +
      labs(title = paste("fgsea", basename(f)),
           subtitle = "Ranked by Spearman rho vs ZEB1 (TARGET)",
           x = "NES", y = NULL) +
      theme_sci()
    save_both(p, file.path(fig, paste0("M4_fgsea_", gsub("[^A-Za-z0-9]+", "_", tools::file_path_sans_ext(basename(f))), ".pdf")),
              9.5, 6.8)
  }
  if (length(all_res)) fwrite(rbindlist(all_res, fill = TRUE), file.path(tab, "M4_fgsea_all.tsv"), sep = "\t")
})

fwrite(st[gene %in% highlight], file.path(tab, "M4_keygene_stats.tsv"), sep = "\t")
cat("MODULE4 ANALYZE DONE\n")
cat("FIG", length(list.files(fig, pattern = "\\.(pdf|png)$")), "\n")

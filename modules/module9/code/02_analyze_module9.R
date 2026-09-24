# Module 9: perturbation public data. Recipes: code2 206 volcano, 193 scatter,
# 88 forest, 210 heatmap, 237 bubble, 346 signed bars. ZEB1 kept continuous in
# patient ranks. No causal LMO2→ZEB1 wording in titles.

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})
has_deseq <- requireNamespace("DESeq2", quietly = TRUE)
has_limma <- requireNamespace("limma", quietly = TRUE)
has_fgsea <- requireNamespace("fgsea", quietly = TRUE)
has_ggrepel <- requireNamespace("ggrepel", quietly = TRUE)
has_cht <- requireNamespace("ComplexHeatmap", quietly = TRUE)
has_circlize <- requireNamespace("circlize", quietly = TRUE)
has_patch <- requireNamespace("patchwork", quietly = TRUE)
if (has_patch) library(patchwork)

root <- "/data-b/liangfuhua/projects/ZEB1_TALL_analysis/module9"
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
                                    error = function(e) {
                                      message("FAIL ", label, ": ", conditionMessage(e),
                                              " | call=", paste(deparse(conditionCall(e)), collapse = " "))
                                    })
fmt_p <- function(p) ifelse(is.na(p), "NA", ifelse(p < 1e-3, "P < 0.001", sprintf("P = %.3f", p)))
p_star <- function(p) ifelse(is.na(p), "ns", ifelse(p < 0.001, "***", ifelse(p < 0.01, "**", ifelse(p < 0.05, "*", "ns"))))

highlight <- c("ZEB1","ZEB2","LMO2","TLX1","TLX3","TAL1","LYL1","GATA3","TCF7","BCL11B",
               "MEF2C","RUNX1","MYC","IL7R","CD34","NOTCH1","SPI1","RAG1","DNTT","LEF1")
pal_grp <- c(control="#4DBBD5", empty="#4DBBD5", siTLX1_1="#E64B35", siTLX1_2="#F39B7F",
             KD="#E64B35", TLX1_OE="#3C5488", WT="#4DBBD5", TG="#E64B35")

read_mat <- function(path) {
  d <- fread(path)
  g <- d[[1]]
  m <- as.matrix(d[, -1, with = FALSE])
  storage.mode(m) <- "double"
  rownames(m) <- as.character(g)
  m
}

limma_two_group <- function(expr, coldata, treat, ctrl) {
  if (!has_limma) stop("limma missing")
  coldata <- as.data.frame(coldata)
  coldata$sample <- as.character(coldata$sample)
  expr <- expr[, coldata$sample, drop = FALSE]
  grp <- factor(ifelse(coldata$group == treat, treat, ifelse(coldata$group == ctrl, ctrl, NA)),
                levels = c(ctrl, treat))
  keep <- !is.na(grp)
  expr <- expr[, keep, drop = FALSE]
  grp <- droplevels(grp[keep])
  design <- model.matrix(~ grp)
  fit <- limma::lmFit(expr, design)
  fit <- limma::eBayes(fit, trend = TRUE)
  tt <- as.data.table(limma::topTable(fit, coef = 2, number = Inf), keep.rownames = "gene")
  setnames(tt, c("logFC", "adj.P.Val", "P.Value"), c("log2FoldChange", "padj", "pvalue"), skip_absent = TRUE)
  if (!"log2FoldChange" %in% names(tt)) tt[, log2FoldChange := logFC]
  if (!"padj" %in% names(tt)) tt[, padj := adj.P.Val]
  if (!"pvalue" %in% names(tt)) tt[, pvalue := P.Value]
  tt[, signed_stat := sign(log2FoldChange) * (-log10(pmax(pvalue, 1e-300)))]
  tt[, lfcSE := NA_real_]
  tt[]
}

deseq_contrast <- function(counts, coldata, group_col, treat, ctrl, min_count = 10) {
  if (!has_deseq) stop("DESeq2 missing")
  coldata <- as.data.frame(coldata)
  coldata$sample <- as.character(coldata$sample)
  counts <- counts[, coldata$sample, drop = FALSE]
  keep <- rowSums(counts >= min_count) >= max(2, floor(ncol(counts) / 3))
  counts <- round(counts[keep, , drop = FALSE])
  coldata[[group_col]] <- factor(as.character(coldata[[group_col]]))
  if (!treat %in% levels(coldata[[group_col]]) || !ctrl %in% levels(coldata[[group_col]])) {
    stop("contrast levels missing: ", paste(levels(coldata[[group_col]]), collapse = ","))
  }
  rownames(coldata) <- coldata$sample
  dds <- DESeq2::DESeqDataSetFromMatrix(countData = counts, colData = coldata,
                                        design = as.formula(paste("~", group_col)))
  dds <- DESeq2::DESeq(dds, quiet = TRUE, fitType = "local")
  res <- as.data.frame(DESeq2::results(dds, contrast = c(group_col, treat, ctrl), independentFiltering = TRUE))
  setDT(res, keep.rownames = "gene")
  res[, signed_stat := sign(log2FoldChange) * (-log10(pmax(pvalue, 1e-300)))]
  res[]
}

volcano_dt <- function(dt, title, file, lfc = "log2FoldChange", padj = "padj",
                       lfc_cut = 0.5, sub = NULL) {
  dv <- copy(dt)
  if (!"gene" %in% names(dv) && "Gene.Name" %in% names(dv)) dv[, gene := as.character(`Gene.Name`)]
  if (!"gene" %in% names(dv)) stop("no gene column: ", paste(names(dv), collapse = ","))
  lfc_name <- lfc; padj_name <- padj
  if (!lfc_name %in% names(dv)) stop("missing LFC column ", lfc_name)
  if (!padj_name %in% names(dv)) stop("missing padj column ", padj_name)
  lfc_vals <- as.numeric(dv[[lfc_name]])
  padj_vals <- as.numeric(dv[[padj_name]])
  dv[, lfc_plot := lfc_vals]
  dv[, ypadj := padj_vals]
  dv <- dv[is.finite(lfc_plot) & is.finite(ypadj) & ypadj >= 0]
  if (!nrow(dv)) stop("volcano: no finite LFC/FDR rows")
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
    labs(title = title, subtitle = sub, x = "log2 fold change", y = "-log10 FDR", color = NULL) +
    theme_sci() + theme(legend.position = "top")
  if (has_ggrepel && nrow(hit))
    p <- p + ggrepel::geom_text_repel(data = as.data.frame(hit), aes(label = gene_up), size = 3, color = "black", max.overlaps = 40, inherit.aes = TRUE)
  save_both(p, file, 7.8, 6.2)
}

keygene_forest <- function(dt_list, file, title) {
  rows <- rbindlist(lapply(names(dt_list), function(nm) {
    d <- copy(dt_list[[nm]])
    d[, gene_up := toupper(as.character(gene))]
    d[gene_up %in% highlight, .(experiment = nm, gene = gene_up, log2FoldChange, padj, lfcSE)]
  }), fill = TRUE)
  if (!nrow(rows)) return(invisible(NULL))
  fwrite(rows, file.path(tab, "M9_keygene_effects.tsv"), sep = "\t")
  rows[, gene := factor(gene, levels = intersect(highlight, unique(gene)))]
  rows[, lo := log2FoldChange - 1.96 * ifelse(is.finite(lfcSE), lfcSE, 0)]
  rows[, hi := log2FoldChange + 1.96 * ifelse(is.finite(lfcSE), lfcSE, 0)]
  p <- ggplot(rows, aes(log2FoldChange, gene, fill = experiment)) +
    geom_vline(xintercept = 0, linetype = 2, color = "grey40") +
    geom_errorbar(aes(xmin = lo, xmax = hi), width = 0.2, position = position_dodge(0.6), color = "grey20") +
    geom_point(shape = 21, size = 2.6, position = position_dodge(0.6), color = "black") +
    labs(title = title, x = "log2 fold change", y = NULL, fill = NULL,
         subtitle = "Whiskers 95% CI when DESeq2 SE available") +
    theme_sci() + theme(legend.position = "bottom")
  save_both(p, file, 8.2, 7.2)
}

pca_plot <- function(counts, coldata, file, title) {
  keep <- rowSums(counts >= 10) >= 2
  vs <- log2(counts[keep, , drop = FALSE] + 1)
  vs <- vs[order(apply(vs, 1, var), decreasing = TRUE)[seq_len(min(2000, nrow(vs)))], , drop = FALSE]
  pc <- prcomp(t(vs), scale. = TRUE)
  dd <- data.table(sample = rownames(pc$x), PC1 = pc$x[, 1], PC2 = pc$x[, 2])
  dd <- merge(dd, as.data.table(coldata), by = "sample")
  ve <- summary(pc)$importance[2, 1:2] * 100
  p <- ggplot(dd, aes(PC1, PC2, fill = group)) +
    geom_point(shape = 21, size = 3.4, color = "black", stroke = 0.35) +
    scale_fill_manual(values = pal_grp, na.value = "grey70") +
    labs(title = title, x = sprintf("PC1 (%.1f%%)", ve[1]), y = sprintf("PC2 (%.1f%%)", ve[2])) +
    theme_sci()
  save_both(p, file, 6.4, 5.4)
}

concord_scatter <- function(a, b, lab_a, lab_b, file, title) {
  aa <- copy(a)[, gene := toupper(as.character(gene))]
  bb <- copy(b)[, gene := toupper(as.character(gene))]
  m <- merge(aa[, .(gene, log2FoldChange, padj)], bb[, .(gene, log2FoldChange, padj)], by = "gene")
  setnames(m, c("gene", "lfc_a", "padj_a", "lfc_b", "padj_b"))
  m <- m[is.finite(lfc_a) & is.finite(lfc_b)]
  rho <- suppressWarnings(cor(m$lfc_a, m$lfc_b, method = "spearman"))
  hit <- m[gene %in% highlight]
  p <- ggplot(m, aes(lfc_a, lfc_b)) +
    geom_hline(yintercept = 0, linetype = 2, color = "grey70") +
    geom_vline(xintercept = 0, linetype = 2, color = "grey70") +
    geom_abline(slope = 1, intercept = 0, linetype = 3, color = "grey50") +
    geom_point(size = 0.3, alpha = 0.25, color = "grey50") +
    geom_point(data = hit, shape = 21, size = 2.4, fill = "white", color = "black") +
    labs(title = title, subtitle = sprintf("Spearman rho = %.2f; n = %d genes", rho, nrow(m)),
         x = lab_a, y = lab_b) +
    theme_sci()
  if (has_ggrepel && nrow(hit))
    p <- p + ggrepel::geom_text_repel(data = hit, aes(label = gene), size = 3, max.overlaps = 40)
  save_both(p, file, 7.2, 6.4)
  invisible(rho)
}

simple_rrho <- function(stat_a, stat_b, file, title, nstep = 20) {
  genes <- intersect(names(stat_a), names(stat_b))
  genes <- genes[!is.na(genes) & genes != ""]
  a <- sort(stat_a[genes], decreasing = TRUE)
  b <- sort(stat_b[genes], decreasing = TRUE)
  n <- length(genes)
  if (n < 200) stop("too few shared genes for RRHO")
  steps <- unique(pmax(50L, as.integer(round(seq(50, n, length.out = nstep)))))
  grid <- rbindlist(lapply(steps, function(i) {
    data.table(i = i, j = steps)
  }))
  na <- names(a); nb <- names(b)
  grid[, p := mapply(function(ii, jj) {
    k <- length(intersect(na[seq_len(ii)], nb[seq_len(jj)]))
    phyper(k - 1, ii, n - ii, jj, lower.tail = FALSE)
  }, i, j)]
  grid[, mlogp := -log10(pmax(p, 1e-16))]
  p <- ggplot(grid, aes(i / n, j / n, fill = mlogp)) +
    geom_tile() +
    scale_fill_gradient(low = "white", high = "#E64B35") +
    labs(title = title, x = "top fraction list A", y = "top fraction list B",
         fill = "-log10 P", subtitle = "Hypergeometric overlap of top ranks (RRHO-style)") +
    theme_sci() + coord_equal()
  save_both(p, file, 6.6, 5.8)
}

mean_rank_rra <- function(stat_list) {
  ranks <- lapply(stat_list, function(s) {
    s <- s[is.finite(s)]
    r <- rank(-s, ties.method = "average") / length(s)
    setNames(r, names(s))
  })
  genes <- Reduce(intersect, lapply(ranks, names))
  mat <- sapply(ranks, `[`, genes)
  data.table(gene = genes, mean_rank = rowMeans(mat), n_exp = ncol(mat))[order(mean_rank)]
}

gsea_cores <- function(stats, file_prefix, title) {
  if (!has_fgsea) return(invisible(NULL))
  stats <- stats[is.finite(stats)]
  stats <- sort(stats[!duplicated(names(stats))], decreasing = TRUE)
  pos <- scan(file.path(proc, "m4_zeb1_pos_core.txt"), what = "", quiet = TRUE)
  neg <- scan(file.path(proc, "m4_zeb1_neg_core.txt"), what = "", quiet = TRUE)
  gs <- list(M4_ZEB1_positive_core = pos, M4_ZEB1_negative_core = neg)
  hm <- file.path(gmt_dir, "h.all.v2023.2.Hs.symbols.gmt")
  if (file.exists(hm)) gs <- c(gs, fgsea::gmtPathways(hm))
  res <- as.data.table(fgsea::fgsea(pathways = gs, stats = stats, minSize = 8, maxSize = 800, nPermSimple = 4000))
  res[, leadingEdge := vapply(leadingEdge, function(x) paste(head(x, 10), collapse = ","), character(1))]
  fwrite(res[order(padj, -abs(NES))], paste0(file_prefix, ".tsv"), sep = "\t")
  d <- res[pathway %in% c("M4_ZEB1_positive_core", "M4_ZEB1_negative_core") |
             grepl("^HALLMARK", pathway)][order(padj)][1:min(18, .N)]
  if (!nrow(d)) return(invisible(res))
  d[, pathway := factor(pathway, levels = rev(pathway))]
  p <- ggplot(d, aes(NES, pathway)) +
    geom_point(aes(size = size, fill = -log10(pmax(padj, 1e-8))), shape = 21, color = "grey20") +
    geom_vline(xintercept = 0, linetype = 2, color = "grey40") +
    scale_fill_gradient(low = "#4DBBD5", high = "#E64B35") +
    labs(title = title, x = "NES", y = NULL, fill = "-log10 FDR") +
    theme_sci()
  save_both(p, paste0(file_prefix, ".pdf"), 8.4, 6.4)
  invisible(res)
}

# ---------- load KD total / polyA ----------
load_kd <- function(prefix) {
  counts <- read_mat(file.path(proc, paste0(prefix, "_symbol_counts.tsv.gz")))
  cd <- fread(file.path(proc, paste0(prefix, "_symbol_coldata.tsv")))
  if (!"sample" %in% names(cd)) cd[, sample := gsm]
  list(counts = counts, cd = cd)
}

res_list <- list()

run("GSE110635 DESeq2", function() {
  x <- load_kd("gse110635")
  cd <- copy(x$cd)
  cd[, sample := gsm]
  cd[, kd := ifelse(group == "control", "control", "KD")]
  pca_plot(x$counts, cd, file.path(fig, "M9_9.2_tlx1_kd_total_PCA.pdf"),
           "ALL-SIL TLX1 KD total RNA-seq PCA")
  r1 <- deseq_contrast(x$counts, cd[group %in% c("control", "siTLX1_1")], "group", "siTLX1_1", "control")
  r2 <- deseq_contrast(x$counts, cd[group %in% c("control", "siTLX1_2")], "group", "siTLX1_2", "control")
  rc <- deseq_contrast(x$counts, cd, "kd", "KD", "control")
  fwrite(r1, file.path(tab, "M9_gse110635_siTLX1_1.tsv"), sep = "\t")
  fwrite(r2, file.path(tab, "M9_gse110635_siTLX1_2.tsv"), sep = "\t")
  fwrite(rc, file.path(tab, "M9_gse110635_combined_KD.tsv"), sep = "\t")
  res_list$gse110635_combined <<- rc
  res_list$gse110635_s1 <<- r1
  res_list$gse110635_s2 <<- r2
  volcano_dt(rc, "TLX1 KD vs control (total RNA, combined siRNAs)",
             file.path(fig, "M9_9.3_tlx1_kd_total_volcano.pdf"),
             sub = "GSE110635 ALL-SIL; DESeq2; n = 6 KD vs 3 NTC")
  concord_scatter(r1, r2, "siTLX1-1 log2FC", "siTLX1-2 log2FC",
                  file.path(fig, "M9_9.4_two_sirna_concordance.pdf"),
                  "Two independent TLX1 siRNAs (total RNA)")
  z <- rc[gene == "ZEB1"]
  message(sprintf("GSE110635 ZEB1 LFC=%.3f padj=%s", z$log2FoldChange[1], z$padj[1]))
})

run("GSE110632 DESeq2", function() {
  x <- load_kd("gse110632")
  cd <- copy(x$cd); cd[, sample := gsm]
  cd[, kd := ifelse(group == "control", "control", "KD")]
  pca_plot(x$counts, cd, file.path(fig, "M9_9.2b_tlx1_kd_polya_PCA.pdf"),
           "ALL-SIL TLX1 KD polyA RNA-seq PCA")
  r1 <- deseq_contrast(x$counts, cd[group %in% c("control", "siTLX1_1")], "group", "siTLX1_1", "control")
  r2 <- deseq_contrast(x$counts, cd[group %in% c("control", "siTLX1_2")], "group", "siTLX1_2", "control")
  rc <- deseq_contrast(x$counts, cd, "kd", "KD", "control")
  fwrite(rc, file.path(tab, "M9_gse110632_combined_KD.tsv"), sep = "\t")
  res_list$gse110632_combined <<- rc
  volcano_dt(rc, "TLX1 KD vs control (polyA, combined siRNAs)",
             file.path(fig, "M9_9.3b_tlx1_kd_polya_volcano.pdf"),
             sub = "GSE110632 ALL-SIL; DESeq2")
  res_list$gse110632_combined <<- rc
  if (!is.null(res_list$gse110635_combined)) {
    concord_scatter(res_list$gse110635_combined, rc,
                    "total RNA log2FC", "polyA log2FC",
                    file.path(fig, "M9_9.5_total_polya_concordance.pdf"),
                    "TLX1 KD: total RNA vs polyA technical designs")
  }
})

run("GSE62143 OE limma", function() {
  expr <- read_mat(file.path(proc, "gse62143_expr.tsv.gz"))
  cd <- fread(file.path(proc, "gse62143_coldata.tsv"))
  if (!"sample" %in% names(cd)) cd[, sample := gsm]
  ro <- limma_two_group(expr, cd, "TLX1_OE", "empty")
  fwrite(ro, file.path(tab, "M9_gse62143_TLX1_OE.tsv"), sep = "\t")
  volcano_dt(ro, "TLX1 overexpression vs empty (CD34+ thymocytes, Agilent)",
             file.path(fig, "M9_9.6_tlx1_oe_volcano.pdf"),
             sub = "GSE62143; two donors; limma-trend on processed signals; n=2 is exploratory")
  res_list$gse62143_oe <<- ro
  if (!is.null(res_list$gse110635_combined)) {
    concord_scatter(res_list$gse110635_combined, ro,
                    "TLX1 KD log2FC (ALL-SIL RNA-seq)", "TLX1 OE log2FC (CD34+ array)",
                    file.path(fig, "M9_9.7_kd_oe_inverse.pdf"),
                    "KD vs OE: inverse concordance expected if TLX1 drives the program")
  }
})

run("GSE62141 KD limma", function() {
  expr <- read_mat(file.path(proc, "gse62141_expr.tsv.gz"))
  cd <- fread(file.path(proc, "gse62141_coldata.tsv"))
  if (!"sample" %in% names(cd)) cd[, sample := gsm]
  cd[, group2 := ifelse(group == "control", "control", "KD")]
  # temporarily map group for limma
  cd2 <- copy(cd)
  cd2[, group := group2]
  rc <- limma_two_group(expr, cd2, "KD", "control")
  fwrite(rc, file.path(tab, "M9_gse62141_KD.tsv"), sep = "\t")
  volcano_dt(rc, "TLX1 KD vs scrambled (GSE62141 Agilent, 24h)",
             file.path(fig, "M9_9.3c_gse62141_volcano.pdf"),
             sub = "Independent ALL-SIL KD platform; limma-trend")
  res_list$gse62141 <<- rc
  if (!is.null(res_list$gse110635_combined)) {
    concord_scatter(res_list$gse110635_combined, rc,
                    "GSE110635 RNA-seq log2FC", "GSE62141 array log2FC",
                    file.path(fig, "M9_9.5b_cross_platform_kd.pdf"),
                    "TLX1 KD RNA-seq vs Agilent concordance")
  }
})

run("keygene forest + RRA + GSEA", function() {
  use <- res_list[c("gse110635_combined", "gse110632_combined", "gse62143_oe", "gse62141")]
  use <- use[!sapply(use, is.null)]
  keygene_forest(use, file.path(fig, "M9_9.8_zeb1_bidirectional_forest.pdf"),
                 "Pre-specified genes across TLX1 perturbations")
  stats <- lapply(use, function(d) setNames(d$signed_stat, d$gene))
  rra <- mean_rank_rra(stats)
  fwrite(rra, file.path(tab, "M9_9.10_RRA_mean_rank.tsv"), sep = "\t")
  top <- rra[1:30]
  top[, gene := factor(gene, levels = rev(gene))]
  p <- ggplot(top, aes(mean_rank, gene)) +
    geom_segment(aes(x = 0, xend = mean_rank, yend = gene), color = "#B07D3A", linewidth = 1) +
    geom_point(shape = 21, size = 3, fill = "#E64B35", color = "black") +
    labs(title = "Consensus TLX1 perturbation ranking (mean fractional rank)",
         subtitle = "Lower rank = more consistently up after ranking by signed statistic",
         x = "mean fractional rank", y = NULL) +
    theme_sci()
  save_both(p, file.path(fig, "M9_9.10_RRA_consensus.pdf"), 7.4, 7.0)
  if (!is.null(res_list$gse110635_combined)) {
    st <- setNames(res_list$gse110635_combined$signed_stat, res_list$gse110635_combined$gene)
    gsea_cores(st, file.path(fig, "M9_9.13_kd_core_gsea"),
               "GSEA of TLX1 KD ranks on M4 ZEB1 cores + Hallmark")
    # also write tsv next to pdf: gsea_cores writes fig prefix tsv; copy
  }
})

run("patient GSE110636 TLX vs rest", function() {
  counts <- read_mat(file.path(proc, "gse110636_symbol_counts.tsv.gz"))
  cd <- fread(file.path(proc, "gse110636_symbol_coldata.tsv"))
  cd[, sample := gsm]
  cd[, group := ifelse(toupper(subgroup) == "TLX", "TLX", "rest")]
  # sample composition
  cc <- cd[, .N, by = subgroup]
  p <- ggplot(cc, aes(reorder(subgroup, -N), N, fill = subgroup)) +
    geom_col(color = "black", width = 0.7) +
    geom_text(aes(label = N), vjust = -0.3, size = 3.2) +
    labs(title = "GSE110636 primary T-ALL subgroups", x = NULL, y = "n patients") +
    theme_sci() + theme(legend.position = "none")
  save_both(p, file.path(fig, "M9_9.1_gse110636_subgroups.pdf"), 6.2, 4.8)
  rc <- deseq_contrast(counts, cd, "group", "TLX", "rest", min_count = 10)
  fwrite(rc, file.path(tab, "M9_gse110636_TLX_vs_rest.tsv"), sep = "\t")
  volcano_dt(rc, "Primary T-ALL: TLX vs rest (GSE110636)",
             file.path(fig, "M9_9.19_patient_TLX_volcano.pdf"),
             sub = sprintf("n TLX = %d, rest = %d", sum(cd$group == "TLX"), sum(cd$group == "rest")))
  res_list$gse110636_tlx <<- rc
  if (!is.null(res_list$gse110635_combined)) {
    concord_scatter(rc, res_list$gse110635_combined,
                    "patient TLX vs rest log2FC", "TLX1 KD log2FC",
                    file.path(fig, "M9_9.11_patient_vs_kd_scatter.pdf"),
                    "Patient TLX association vs TLX1 KD")
    sa <- setNames(rc$signed_stat, rc$gene)
    sb <- setNames(res_list$gse110635_combined$signed_stat, res_list$gse110635_combined$gene)
    simple_rrho(sa, sb, file.path(fig, "M9_9.11_patient_vs_kd_RRHO.pdf"),
                "RRHO: patient TLX vs rest  vs  TLX1 KD")
  }
})

run("TARGET TLX RRHO", function() {
  f <- file.path(proc, "target_tlx_vs_rest.tsv.gz")
  if (!file.exists(f) || is.null(res_list$gse110635_combined)) return()
  pt <- fread(f)
  sa <- setNames(pt$signed_stat, pt$gene)
  sb <- setNames(res_list$gse110635_combined$signed_stat, res_list$gse110635_combined$gene)
  simple_rrho(sa, sb, file.path(fig, "M9_9.11b_TARGET_vs_kd_RRHO.pdf"),
              "RRHO: TARGET TLX vs rest  vs  TLX1 KD")
  concord_scatter(
    pt[, .(gene, log2FoldChange = lfc_TLX_minus_rest, padj = p)],
    res_list$gse110635_combined,
    "TARGET TLX vs rest LFC", "TLX1 KD log2FC",
    file.path(fig, "M9_9.11b_TARGET_vs_kd_scatter.pdf"),
    "TARGET TLX patients vs TLX1 KD"
  )
})

run("GSE186943 LMO2 TG vs WT", function() {
  fpkm <- read_mat(file.path(proc, "gse186943_fpkm.tsv.gz"))
  cd <- fread(file.path(proc, "gse186943_coldata.tsv"))
  # keep RNA samples only (exp_ columns that are WT/TG, not chip)
  cd <- cd[genotype %in% c("WT", "TG")]
  common <- intersect(colnames(fpkm), cd$sample)
  fpkm <- fpkm[, common, drop = FALSE]
  cd <- cd[sample %in% common]
  loge <- log2(fpkm + 1)
  # gene-wise limma ~ genotype + stage if stages vary
  design <- model.matrix(~ genotype + stage, data = as.data.frame(cd))
  if (!has_limma) stop("limma missing")
  fit <- limma::lmFit(loge, design)
  fit <- limma::eBayes(fit, trend = TRUE)
  coef <- grep("genotypeTG|genotypeWT", colnames(design), value = TRUE)
  # design genotype is factor; TG vs WT depends on reference
  cn <- colnames(fit$coefficients)
  gcoef <- cn[grepl("genotype", cn)][1]
  tt <- as.data.table(limma::topTable(fit, coef = gcoef, number = Inf), keep.rownames = "gene")
  setnames(tt, c("logFC", "adj.P.Val"), c("log2FoldChange", "padj"), skip_absent = TRUE)
  if (!"log2FoldChange" %in% names(tt) && "logFC" %in% names(tt)) tt[, log2FoldChange := logFC]
  if (!"padj" %in% names(tt) && "adj.P.Val" %in% names(tt)) tt[, padj := `adj.P.Val`]
  fwrite(tt, file.path(tab, "M9_gse186943_TG_vs_WT.tsv"), sep = "\t")
  volcano_dt(tt, "Mouse LMO2 transgenic vs WT (all stages, limma-trend)",
             file.path(fig, "M9_9.16_lmo2_oe_volcano.pdf"),
             sub = "GSE186943 FPKM; genotype + stage")
  res_list$gse186943 <<- tt
  # Zeb1 by stage
  if ("Zeb1" %in% rownames(loge) || "ZEB1" %in% rownames(loge)) {
    g <- intersect(c("Zeb1", "ZEB1"), rownames(loge))[1]
    dd <- data.table(sample = colnames(loge), expr = as.numeric(loge[g, ]))
    dd <- merge(dd, cd, by = "sample")
    p <- ggplot(dd, aes(stage, expr, fill = genotype)) +
      geom_boxplot(width = 0.6, outlier.shape = NA, position = position_dodge(0.75)) +
      geom_point(position = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.75),
                 shape = 21, size = 1.8, color = "black") +
      scale_fill_manual(values = c(WT = "#4DBBD5", TG = "#E64B35")) +
      labs(title = paste(g, "in LMO2 transgenic vs WT thymocytes"),
           y = "log2(FPKM+1)", x = NULL) +
      theme_sci() + theme(axis.text.x = element_text(angle = 35, hjust = 1))
    save_both(p, file.path(fig, "M9_9.16b_Zeb1_by_stage.pdf"), 8.0, 5.2)
  }
})

run("GSE186943 ChIP peaks", function() {
  pk <- fread(file.path(proc, "gse186943_chip_zeb1.tsv"))
  fwrite(pk, file.path(tab, "M9_9.19_lmo2_chip_zeb1.tsv"), sep = "\t")
  p <- ggplot(pk, aes(reorder(file, n_Zeb1_window), n_Zeb1_window, fill = n_Zeb1_window > 0)) +
    geom_col(color = "black", width = 0.7) +
    coord_flip() +
    scale_fill_manual(values = c("TRUE" = "#E64B35", "FALSE" = "grey75"), guide = "none") +
    labs(title = "LMO2 ChIP peaks overlapping the mouse Zeb1 window",
         subtitle = "chr18:5.60–5.85 Mb (mm10/mm39-tolerant)",
         x = NULL, y = "n peaks in Zeb1 window") +
    theme_sci()
  save_both(p, file.path(fig, "M9_9.19_lmo2_chip_Zeb1.pdf"), 7.6, 4.8)
})

run("GSE144035 Lmo2 ON/OFF", function() {
  f <- file.path(proc, "gse144035_precomputed_de.tsv.gz")
  if (!file.exists(f)) stop("no precomputed DE")
  de <- fread(f)
  if ("Gene.Name" %in% names(de)) de[, gene := `Gene.Name`]
  lfc_col <- names(de)[grepl("Independent|Dependent|log", names(de), ignore.case = TRUE)][1]
  fdr_col <- names(de)[grepl("^FDR$|adj", names(de), ignore.case = TRUE)][1]
  if (is.na(lfc_col) || is.na(fdr_col)) stop(paste("cols", paste(names(de), collapse = ",")))
  de[, log2FoldChange := as.numeric(de[[lfc_col]])]
  de[, padj := as.numeric(de[[fdr_col]])]
  fwrite(de, file.path(tab, "M9_gse144035_precomputed.tsv"), sep = "\t")
  volcano_dt(de, paste0("Reversible Lmo2 model: ", lfc_col),
             file.path(fig, "M9_9.17_lmo2_OFF_volcano.pdf"),
             sub = "GSE144035 author-computed table; check ZEB1/Zeb1 direction")
  res_list$gse144035 <<- de
})

run("GSE188225 Lmo2 TG", function() {
  counts <- read_mat(file.path(proc, "gse188225_counts.tsv.gz"))
  cd <- fread(file.path(proc, "gse188225_coldata.tsv"))
  cd <- cd[genotype %in% c("WT", "TG")]
  # Entrez IDs: map a few key genes only for forest; genome-wide keep Entrez
  rc <- deseq_contrast(counts, cd, "genotype", "TG", "WT", min_count = 5)
  fwrite(rc, file.path(tab, "M9_gse188225_TG_vs_WT.tsv"), sep = "\t")
  # map Entrez of mouse key genes
  entrez <- c(Zeb1 = "21417", Zeb2 = "24136", Lmo2 = "16909", Gata3 = "14462",
              Tcf7 = "21414", Bcl11b = "58208", Tal1 = "21349", Lyl1 = "17095",
              Il7r = "16197", Cd34 = "12490", Myc = "17869", Runx1 = "12394")
  rc[, gene_symbol := names(entrez)[match(gene, entrez)]]
  hit <- rc[!is.na(gene_symbol)]
  fwrite(hit, file.path(tab, "M9_gse188225_keygenes.tsv"), sep = "\t")
  if (nrow(hit)) {
    p <- ggplot(hit, aes(log2FoldChange, gene_symbol, fill = log2FoldChange > 0)) +
      geom_vline(xintercept = 0, linetype = 2) +
      geom_col(width = 0.7, color = "black") +
      scale_fill_manual(values = c("TRUE" = "#E64B35", "FALSE" = "#3C5488"), guide = "none") +
      labs(title = "Lmo2-transgenic DN2/DN3 (GSE188225)", y = NULL, x = "log2FC TG vs WT") +
      theme_sci()
    save_both(p, file.path(fig, "M9_9.16c_gse188225_keygenes.pdf"), 6.4, 5.0)
  }
})

run("human LMO2 ChIP occupancy", function() {
  f <- file.path(proc, "gse154675_chip_zeb1.tsv")
  if (!file.exists(f)) stop("no chip table")
  d <- fread(f)
  fwrite(d, file.path(tab, "M9_9.19b_human_chip_zeb1.tsv"), sep = "\t")
  if (!"file" %in% names(d) || nrow(d) == 0 || !"hg38_gene" %in% names(d) && !"status" %in% names(d)) {
    message("human ChIP skipped or empty")
    return()
  }
  if ("status" %in% names(d)) return()
  dm <- melt(d, id.vars = "file", measure.vars = intersect(c("hg38_gene", "hg19_gene", "hg38p_prom", "hg19p_prom"), names(d)),
             variable.name = "window", value.name = "mean_signal")
  p <- ggplot(dm, aes(window, mean_signal, fill = grepl("LMO2", file))) +
    geom_boxplot(outlier.shape = NA, width = 0.6) +
    geom_point(position = position_jitter(width = 0.12), shape = 21, size = 2) +
    labs(title = "GSE154675 bigWig mean signal at ZEB1",
         subtitle = "Human T-ALL LMO2/TAL1 ChIP; occupancy is not proof of repression",
         y = "mean bigWig signal", x = NULL, fill = "LMO2 track") +
    theme_sci()
  save_both(p, file.path(fig, "M9_9.19b_human_LMO2_ChIP_ZEB1.pdf"), 7.4, 5.2)
})

run("evidence matrix", function() {
  rows <- list()
  for (nm in names(res_list)) {
    d <- res_list[[nm]]
    if (is.null(d) || !nrow(d) || !"log2FoldChange" %in% names(d)) next
    tmp <- copy(d)
    gcol <- if ("gene_symbol" %in% names(tmp) && any(!is.na(tmp$gene_symbol))) "gene_symbol" else "gene"
    if (!gcol %in% names(tmp)) next
    tmp[, g := toupper(as.character(get(gcol)))]
    tmp <- tmp[g %in% highlight]
    if (!nrow(tmp)) next
    if (!"padj" %in% names(tmp)) tmp[, padj := NA_real_]
    rows[[nm]] <- tmp[, .(experiment = nm, gene = g[1], log2FoldChange = log2FoldChange[1], padj = padj[1]), by = g]
  }
  ev <- rbindlist(rows, fill = TRUE)
  if (!nrow(ev)) stop("empty evidence")
  fwrite(ev, file.path(tab, "M9_9.24_evidence_matrix.tsv"), sep = "\t")
  mat <- dcast(ev, gene ~ experiment, value.var = "log2FoldChange")
  fwrite(mat, file.path(tab, "M9_9.24_evidence_LFC.tsv"), sep = "\t")
  mm <- as.matrix(mat[, -1, with = FALSE])
  rownames(mm) <- mat$gene
  mm[!is.finite(mm)] <- 0
  if (has_cht && has_circlize) {
    col_fun <- circlize::colorRamp2(c(-1.5, 0, 1.5), c("#3C5488", "white", "#E64B35"))
    ht <- ComplexHeatmap::Heatmap(mm, name = "log2FC", col = col_fun,
                                  cluster_rows = FALSE, cluster_columns = FALSE,
                                  column_title = "Perturbation evidence matrix (log2FC)",
                                  row_names_gp = grid::gpar(fontsize = 10),
                                  column_names_rot = 55,
                                  cell_fun = function(j, i, x, y, w, h, fill) {
                                    grid::grid.text(sprintf("%.2f", mm[i, j]), x, y, gp = grid::gpar(fontsize = 7))
                                  })
    pdf(file.path(fig, "M9_9.24_evidence_heatmap.pdf"), width = 9.5, height = 6.5)
    print(ht); dev.off()
    png(file.path(fig, "M9_9.24_evidence_heatmap.png"), width = 9.5 * 180, height = 6.5 * 180, res = 180)
    print(ht); dev.off()
  }
})

run("KD heatmap key genes", function() {
  x <- load_kd("gse110635")
  cd <- copy(x$cd); cd[, sample := gsm]
  genes <- intersect(highlight, rownames(x$counts))
  mat <- log2(x$counts[genes, cd$sample, drop = FALSE] + 1)
  keep <- apply(mat, 1, function(z) sd(z) > 1e-8)
  mat <- mat[keep, , drop = FALSE]
  mat <- t(scale(t(mat)))
  ha <- NULL
  if (has_cht) {
    if (has_circlize) {
      col_fun <- circlize::colorRamp2(c(-2, 0, 2), c("#3C5488", "white", "#E64B35"))
    } else col_fun <- NULL
    ha <- ComplexHeatmap::HeatmapAnnotation(group = cd$group, col = list(group = pal_grp))
    pdf(file.path(fig, "M9_9.14_kd_keygene_heatmap.pdf"), width = 7.2, height = 6.2)
    print(ComplexHeatmap::Heatmap(mat, name = "z", top_annotation = ha, cluster_columns = FALSE,
                                  column_split = cd$group, column_title = "GSE110635 key genes"))
    dev.off()
  }
})

run("sample design counts", function() {
  d635 <- fread(file.path(proc, "gse110635_symbol_coldata.tsv"))
  d632 <- fread(file.path(proc, "gse110632_symbol_coldata.tsv"))
  d643 <- fread(file.path(proc, "gse62143_coldata.tsv"))
  a <- rbind(
    d635[, .(dataset = "GSE110635 total KD", group, n = 1)],
    d632[, .(dataset = "GSE110632 polyA KD", group, n = 1)],
    d643[, .(dataset = "GSE62143 TLX1 OE", group, n = 1)]
  )
  cc <- a[, .N, by = .(dataset, group)]
  p <- ggplot(cc, aes(group, N, fill = group)) +
    geom_col(color = "black", width = 0.7) +
    facet_wrap(~dataset, scales = "free_x") +
    scale_fill_manual(values = pal_grp, na.value = "grey70") +
    labs(title = "Module 9 perturbation sample sizes", y = "n samples", x = NULL) +
    theme_sci() + theme(legend.position = "none", axis.text.x = element_text(angle = 30, hjust = 1))
  save_both(p, file.path(fig, "M9_9.1_design_samples.pdf"), 9.0, 4.6)
})

run("rebuild volcanoes from written tables", function() {
  specs <- list(
    list(f = "M9_gse110635_combined_KD.tsv", title = "TLX1 KD vs control (total RNA, combined siRNAs)",
         out = "M9_9.3_tlx1_kd_total_volcano.pdf", sub = "GSE110635 ALL-SIL; DESeq2; n = 6 KD vs 3 NTC"),
    list(f = "M9_gse110632_combined_KD.tsv", title = "TLX1 KD vs control (polyA, combined siRNAs)",
         out = "M9_9.3b_tlx1_kd_polya_volcano.pdf", sub = "GSE110632 ALL-SIL; DESeq2"),
    list(f = "M9_gse62143_TLX1_OE.tsv", title = "TLX1 overexpression vs empty (CD34+ thymocytes, Agilent)",
         out = "M9_9.6_tlx1_oe_volcano.pdf", sub = "GSE62143; limma-trend; n=2 is exploratory"),
    list(f = "M9_gse62141_KD.tsv", title = "TLX1 KD vs scrambled (GSE62141 Agilent, 24h)",
         out = "M9_9.3c_gse62141_volcano.pdf", sub = "Independent ALL-SIL KD platform; limma-trend"),
    list(f = "M9_gse110636_TLX_vs_rest.tsv", title = "Primary T-ALL: TLX vs rest (GSE110636)",
         out = "M9_9.19_patient_TLX_volcano.pdf", sub = "DESeq2"),
    list(f = "M9_gse186943_TG_vs_WT.tsv", title = "Mouse LMO2 transgenic vs WT (limma-trend)",
         out = "M9_9.16_lmo2_oe_volcano.pdf", sub = "GSE186943 FPKM; genotype + stage"),
    list(f = "M9_gse144035_precomputed.tsv", title = "Reversible Lmo2 model (author DE table)",
         out = "M9_9.17_lmo2_OFF_volcano.pdf", sub = "GSE144035")
  )
  for (s in specs) {
    fp <- file.path(tab, s$f)
    if (!file.exists(fp)) { message("skip missing ", s$f); next }
    d <- fread(fp)
    message("volcano ", s$f, " nrow=", nrow(d), " cols=", paste(names(d), collapse = ","))
    tryCatch(
      volcano_dt(d, s$title, file.path(fig, s$out), sub = s$sub),
      error = function(e) message("VOLCANO FAIL ", s$f, ": ", conditionMessage(e), " call=", deparse(conditionCall(e)))
    )
  }
})

message("ANALYZE DONE")
message("FIG_N ", length(list.files(fig)))
writeLines(capture.output(sessionInfo()), file.path(tab, "sessionInfo.txt"))

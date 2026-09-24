# Module 4 high-order analyses, numbered M4_4.x
# Recipes: code1 33 multi-GSEA, 25 cor-network; code2 206 volcano, 230 enrich-bar,
# 237 bubble, 210 heatmap, 88 forest, 312 ssGSEA, 84 GSEA running-ES.

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})
has_fgsea <- requireNamespace("fgsea", quietly = TRUE)
has_enrichplot <- requireNamespace("enrichplot", quietly = TRUE)
has_cp <- requireNamespace("clusterProfiler", quietly = TRUE)
has_gsva <- requireNamespace("GSVA", quietly = TRUE)
has_igraph <- requireNamespace("igraph", quietly = TRUE)
has_ggrepel <- requireNamespace("ggrepel", quietly = TRUE)
has_cht <- requireNamespace("ComplexHeatmap", quietly = TRUE)
has_circlize <- requireNamespace("circlize", quietly = TRUE)
has_patch <- requireNamespace("patchwork", quietly = TRUE)
if (has_patch) library(patchwork)

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

st <- fread(file.path(proc, "target_zeb1_gene_stats.tsv"))
ph <- fread(file.path(proc, "pharmacotype_zeb1_gene_stats.tsv"))
both <- fread(file.path(proc, "cross_cohort_gene_stats.tsv"))
lab <- fread(file.path(proc, "target_samples.tsv"))
labp <- fread(file.path(proc, "pharmacotype_samples.tsv"))
expr <- as.data.frame(fread(file.path(proc, "target_log2tpm_filtered.tsv.gz")), stringsAsFactors = FALSE)
rownames(expr) <- expr[[1]]; expr[[1]] <- NULL
expr <- as.matrix(expr)
pexpr <- as.data.frame(fread(file.path(proc, "pharmacotype_log2fpkm_filtered.tsv.gz")), stringsAsFactors = FALSE)
rownames(pexpr) <- pexpr[[1]]; pexpr[[1]] <- NULL
pexpr <- as.matrix(pexpr)

ranks <- st[!is.na(rho_zeb1), setNames(rho_zeb1, gene)]
ranks <- ranks[!duplicated(names(ranks))]
ranks <- sort(ranks, decreasing = TRUE)

read_gmt <- function(path) {
  if (!file.exists(path) || !has_fgsea) return(NULL)
  fgsea::gmtPathways(path)
}

gsea_one <- function(gs, ranks, minSize = 10, maxSize = 500) {
  if (is.null(gs) || !has_fgsea) return(data.table())
  res <- fgsea::fgsea(pathways = gs, stats = ranks, minSize = minSize, maxSize = maxSize, nPermSimple = 4000)
  res <- as.data.table(res)
  res[, leadingEdge := vapply(leadingEdge, function(x) paste(head(x, 12), collapse = ","), character(1))]
  res[order(padj, -abs(NES))]
}

bubble_gsea <- function(res, title, outfile, n = 18) {
  d <- res[!is.na(NES) & !is.na(padj)][1:min(n, nrow(res))]
  d[, pathway := factor(pathway, levels = rev(pathway))]
  p <- ggplot(d, aes(NES, pathway)) +
    geom_point(aes(size = size, fill = -log10(pmax(padj, 1e-8))), shape = 21, color = "grey20") +
    geom_vline(xintercept = 0, linetype = 2, color = "grey40") +
    scale_fill_gradient(low = "#4DBBD5", high = "#E64B35") +
    scale_size(range = c(2.5, 8)) +
    labs(title = title, x = "NES (ranked by ZEB1 Spearman)", y = NULL,
         fill = "-log10 padj", size = "n genes") +
    theme_sci()
  save_both(p, outfile, 9.4, 6.8)
}

bar_onbar <- function(res, title, outfile, n = 12) {
  d <- res[!is.na(NES)][1:min(n, nrow(res))]
  setorder(d, NES)
  d[, pathway := factor(pathway, levels = pathway)]
  d[, lab := gsub("^HALLMARK_|^REACTOME_|^GOBP_", "", pathway)]
  p <- ggplot(d, aes(NES, pathway, fill = NES > 0)) +
    geom_col(width = 0.78, color = "white", alpha = 0.55) +
    geom_text(aes(x = 0, label = lab), hjust = 0, size = 2.7, nudge_x = 0.02) +
    scale_fill_manual(values = c("TRUE" = "#E64B35", "FALSE" = "#3C5488"), guide = "none") +
    labs(title = title, x = "NES", y = NULL) +
    theme_sci() + theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(), axis.line.y = element_blank())
  save_both(p, outfile, 8.8, 6.4)
}

# ---------------------------------------------------------------------
# 4.2 Q1 vs Q4 MA-style (mean vs difference)
# ---------------------------------------------------------------------
run("4.2 MA", function() {
  d <- copy(st)
  d[, mean_ab := (mean_Q1 + mean_Q4) / 2]
  d[, cls := fcase(fdr_wilcox < 0.05 & lfc_Q4_minus_Q1 > 0.5, "Q4 high",
                   fdr_wilcox < 0.05 & lfc_Q4_minus_Q1 < -0.5, "Q1 high", default = "NS")]
  hit <- d[gene %in% c("ZEB1", "ZEB2", "LMO2", "LYL1", "GATA3", "TCF7", "BCL11B", "IL7R", "RAG1")]
  p <- ggplot(d, aes(mean_ab, lfc_Q4_minus_Q1, color = cls)) +
    geom_hline(yintercept = 0, linetype = 2, color = "grey50") +
    geom_point(size = 0.35, alpha = 0.4) +
    geom_point(data = hit, size = 2.2, shape = 21, fill = "white", color = "black") +
    scale_color_manual(values = c("Q4 high" = "#E64B35", "Q1 high" = "#3C5488", NS = "grey75")) +
    labs(title = "4.2  TARGET Q4 vs Q1 MA plot",
         subtitle = "Quartiles are display-only; statistics also reported as continuous rho",
         x = "Mean log2(TPM+1)", y = "Q4 - Q1", color = NULL) +
    theme_sci() + theme(legend.position = "top")
  if (has_ggrepel) p <- p + ggrepel::geom_text_repel(data = hit, aes(label = gene), size = 3, color = "black")
  save_both(p, file.path(fig, "M4_4.2_Q1Q4_MA.pdf"), 7.6, 6.0)
})

# ---------------------------------------------------------------------
# 4.6-4.8 GSEA
# ---------------------------------------------------------------------
h_gmt <- file.path(gmt_dir, "h.all.v2023.2.Hs.symbols.gmt")
r_gmt <- file.path(gmt_dir, "c2.cp.reactome.v2023.2.Hs.symbols.gmt")
g_gmt <- file.path(gmt_dir, "c5.go.bp.v2023.2.Hs.symbols.gmt")
t_gmt <- file.path(gmt_dir, "c3.tft.v2023.2.Hs.symbols.gmt")
hs <- read_gmt(h_gmt); rs <- read_gmt(r_gmt); gs <- read_gmt(g_gmt); ts <- read_gmt(t_gmt)

res_h <- gsea_one(hs, ranks)
res_r <- gsea_one(rs, ranks, maxSize = 400)
res_g <- gsea_one(gs, ranks, minSize = 20, maxSize = 250)
res_t <- gsea_one(ts, ranks, minSize = 10, maxSize = 400)
if (nrow(res_h)) fwrite(res_h, file.path(tab, "M4_4.6_hallmark_fgsea.tsv"), sep = "\t")
if (nrow(res_r)) fwrite(res_r, file.path(tab, "M4_4.7_reactome_fgsea.tsv"), sep = "\t")
if (nrow(res_g)) fwrite(res_g, file.path(tab, "M4_4.8_gobp_fgsea.tsv"), sep = "\t")
if (nrow(res_t)) fwrite(res_t, file.path(tab, "M4_4.20_tft_fgsea.tsv"), sep = "\t")

run("4.6 hallmark vis", function() {
  if (!nrow(res_h)) stop("no hallmark")
  bubble_gsea(res_h, "4.6  Hallmark GSEA (code2 237 bubble)", file.path(fig, "M4_4.6_hallmark_bubble.pdf"))
  bar_onbar(res_h, "4.6  Hallmark NES (code2 230)", file.path(fig, "M4_4.6_hallmark_bar.pdf"))
  if (has_cp && has_enrichplot && length(hs)) {
    t2g <- rbindlist(lapply(names(hs), function(p) data.table(term = p, gene = hs[[p]])))
    eg <- clusterProfiler::GSEA(ranks, TERM2GENE = t2g, pvalueCutoff = 1, minGSSize = 10, maxGSSize = 500, verbose = FALSE)
    ids <- head(eg@result$ID[order(eg@result$p.adjust)], 3)
    if (length(ids)) {
      p <- enrichplot::gseaplot2(eg, geneSetID = ids, pvalue_table = TRUE, base_size = 10)
      pdf(file.path(fig, "M4_4.6_hallmark_runningES.pdf"), width = 8.6, height = 7.2)
      print(p); dev.off()
      png(file.path(fig, "M4_4.6_hallmark_runningES.png"), width = 8.6 * 180, height = 7.2 * 180, res = 180)
      print(p); dev.off()
    }
  }
})
run("4.7 reactome vis", function() {
  if (!nrow(res_r)) stop("no reactome gmt")
  bubble_gsea(res_r, "4.7  Reactome GSEA (code2 237)", file.path(fig, "M4_4.7_reactome_bubble.pdf"), n = 20)
  bar_onbar(res_r, "4.7  Reactome NES (code2 230)", file.path(fig, "M4_4.7_reactome_bar.pdf"), n = 14)
})
run("4.8 GO-BP vis", function() {
  if (!nrow(res_g)) stop("no gobp gmt")
  bubble_gsea(res_g, "4.8  GO-BP GSEA (code2 237)", file.path(fig, "M4_4.8_gobp_bubble.pdf"), n = 20)
  bar_onbar(res_g, "4.8  GO-BP NES (code2 230)", file.path(fig, "M4_4.8_gobp_bar.pdf"), n = 14)
})

# ---------------------------------------------------------------------
# 4.9-4.19 signature scores (ssGSEA if GSVA present)
# ---------------------------------------------------------------------
sigs <- list(
  Tcell_differentiation = c("CD1A","CD1B","CD3E","CD3D","CD4","CD8A","RAG1","RAG2","DNTT","BCL11B","TCF7","LEF1","CD5","CD2","LCK","TRAC"),
  ETP = c("CD34","KIT","IL7R","LMO2","LYL1","HHEX","WT1","BAALC","IGLL1","SPINK2","IGFBP7","MEF2C","SPI1","CD44","HES1"),
  Stemness = c("PROM1","BMI1","KIT","ABCG2","SOX4","MYCN","ALDH1A1","THY1","NES","MSI2","HOXA9","MEIS1"),
  Proliferation = c("MKI67","PCNA","TOP2A","CDK1","CCNB1","CCNA2","MCM2","MCM7","BIRC5","UBE2C"),
  Apoptosis = c("BAX","BAK1","BCL2","BCL2L1","MCL1","CASP3","CASP8","FAS","PMAIP1","BBC3"),
  TGFb = c("TGFB1","TGFBR1","TGFBR2","SMAD2","SMAD3","SMAD4","SERPINE1","ID1","SKIL"),
  IL7_JAK_STAT = c("IL7R","JAK1","JAK3","STAT5A","STAT5B","STAT1","IL2RG","SOCS1","CISH"),
  NOTCH = c("NOTCH1","NOTCH3","RBPJ","HES1","HEY1","DTX1","MYC","NRARP"),
  MYC = c("MYC","MAX","CAD","LDHA","NCL","NPM1","EIF4E","ODC1"),
  PI3K_AKT_mTOR = c("AKT1","MTOR","RPS6KB1","EIF4EBP1","PIK3CA","PTEN","GSK3B","TSC2"),
  DNA_repair = c("BRCA1","BRCA2","RAD51","XRCC6","PRKDC","ATM","ATR","CHEK1","PCNA")
)
# overlay hallmark gene sets when present
if (!is.null(hs)) {
  pick <- c(
    Proliferation = "HALLMARK_E2F_TARGETS",
    Apoptosis = "HALLMARK_APOPTOSIS",
    TGFb = "HALLMARK_TGF_BETA_SIGNALING",
    IL7_JAK_STAT = "HALLMARK_IL2_STAT5_SIGNALING",
    NOTCH = "HALLMARK_NOTCH_SIGNALING",
    MYC = "HALLMARK_MYC_TARGETS_V1",
    PI3K_AKT_mTOR = "HALLMARK_PI3K_AKT_MTOR_SIGNALING",
    DNA_repair = "HALLMARK_DNA_REPAIR"
  )
  for (nm in names(pick)) if (pick[[nm]] %in% names(hs)) sigs[[nm]] <- hs[[pick[[nm]]]]
}

score_sets <- function(mat, sets) {
  out <- sapply(sets, function(g) {
    gg <- intersect(g, rownames(mat))
    if (length(gg) < 3) return(rep(NA_real_, ncol(mat)))
    z <- t(scale(t(mat[gg, , drop = FALSE])))
    colMeans(z, na.rm = TRUE)
  })
  as.data.frame(out, check.names = FALSE)
}

run("4.9-4.19 signatures", function() {
  sc <- score_sets(expr, sigs)
  sc$usi <- if (!is.null(rownames(sc)) && !identical(rownames(sc), as.character(seq_len(nrow(sc))))) rownames(sc) else colnames(expr)
  sc <- merge(lab[, .(usi, ZEB1_log, subtype, etp)], sc, by = "usi")
  if (has_gsva) {
    ss <- tryCatch(GSVA::gsva(expr, sigs, method = "ssgsea", ssgsea.norm = TRUE, verbose = FALSE), error = function(e) NULL)
    if (!is.null(ss)) {
      ss <- as.data.frame(t(ss)); ss$usi <- rownames(ss)
      sc <- merge(lab[, .(usi, ZEB1_log, subtype, etp)], ss, by = "usi")
    }
  }
  fwrite(sc, file.path(tab, "M4_4.9_4.19_signature_scores.tsv"), sep = "\t")
  rows <- lapply(names(sigs), function(nm) {
    ct <- suppressWarnings(cor.test(sc$ZEB1_log, sc[[nm]], method = "spearman", exact = FALSE))
    data.table(signature = nm, rho = unname(ct$estimate), p = ct$p.value, n = sum(is.finite(sc[[nm]])))
  })
  sigtab <- rbindlist(rows)
  sigtab[, fdr := p.adjust(p, "BH")]
  setorder(sigtab, rho)
  sigtab[, signature := factor(signature, levels = signature)]
  fwrite(sigtab, file.path(tab, "M4_4.9_4.19_signature_vs_ZEB1.tsv"), sep = "\t")
  p <- ggplot(sigtab, aes(rho, signature, fill = rho > 0)) +
    geom_vline(xintercept = 0, linetype = 2, color = "grey40") +
    geom_col(width = 0.7, color = "white") +
    geom_text(aes(label = sprintf("rho=%.2f  FDR=%s", rho, ifelse(fdr < 0.001, "<0.001", sprintf("%.3f", fdr)))),
              hjust = ifelse(sigtab$rho >= 0, -0.08, 1.08), size = 2.8) +
    scale_fill_manual(values = c("TRUE" = "#E64B35", "FALSE" = "#3C5488"), guide = "none") +
    labs(title = "4.9-4.19  Signature scores vs ZEB1 (ssGSEA/z-mean)",
         subtitle = "TARGET n = 265; Hallmark sets used when the GMT is present",
         x = "Spearman rho vs ZEB1", y = NULL) +
    theme_sci()
  save_both(p, file.path(fig, "M4_4.9_4.19_signature_forest.pdf"), 8.6, 6.4)

  long <- melt(as.data.table(sc), id.vars = c("usi", "ZEB1_log", "subtype", "etp"),
               measure.vars = names(sigs), variable.name = "signature", value.name = "score")
  p2 <- ggplot(long[signature %in% c("Tcell_differentiation", "ETP", "Stemness", "Proliferation")],
               aes(ZEB1_log, score, fill = subtype)) +
    geom_smooth(method = "lm", formula = y ~ x, color = "grey25", fill = "#F3E5C8",
                inherit.aes = FALSE, aes(x = ZEB1_log, y = score), linewidth = 0.5) +
    geom_point(size = 1.6, shape = 21, color = "black", stroke = 0.2, alpha = 0.85) +
    facet_wrap(~signature, scales = "free_y") +
    scale_fill_manual(values = c(TAL = "#E64B35", TLX = "#4DBBD5", HOXA = "#00A087",
                                 "LMO2/LYL1" = "#3C5488", NKX2_1 = "#F39B7F", Unknown = "#8491B4")) +
    labs(title = "Key signatures vs continuous ZEB1", x = "log2(TPM+1) ZEB1", y = "Signature score", fill = NULL) +
    theme_sci() + theme(legend.position = "top")
  save_both(p2, file.path(fig, "M4_4.9_4.11_signature_scatter.pdf"), 9.2, 6.8)

  tryCatch({
    sdf <- as.data.frame(sc)
    keep <- intersect(names(sigs), names(sdf))
    mm <- sapply(keep, function(nm) as.numeric(as.character(sdf[[nm]])))
    if (is.null(dim(mm))) mm <- matrix(mm, ncol = 1, dimnames = list(sdf$usi, keep))
    rownames(mm) <- as.character(sdf$usi)
    mm <- mm[is.finite(rowSums(mm)), , drop = FALSE]
    if (has_cht && has_circlize && nrow(mm) > 10) {
      mat <- t(scale(mm))
      mat[!is.finite(mat)] <- 0
      ord <- as.character(sc$usi)[order(as.numeric(sc$ZEB1_log))]
      ord <- intersect(ord, colnames(mat))
      mat <- mat[, ord, drop = FALSE]
      ha <- ComplexHeatmap::HeatmapAnnotation(
        ZEB1 = as.numeric(sc$ZEB1_log[match(ord, as.character(sc$usi))]),
        subtype = as.character(sc$subtype[match(ord, as.character(sc$usi))]),
        col = list(
          ZEB1 = circlize::colorRamp2(range(as.numeric(sc$ZEB1_log), na.rm = TRUE), c("#3C5488", "#E64B35")),
          subtype = c(TAL = "#E64B35", TLX = "#4DBBD5", HOXA = "#00A087",
                      "LMO2/LYL1" = "#3C5488", NKX2_1 = "#F39B7F", Unknown = "#8491B4")
        )
      )
      pdf(file.path(fig, "M4_4.9_4.19_signature_heatmap.pdf"), width = 10.2, height = 4.8)
      print(ComplexHeatmap::Heatmap(mat, name = "z", top_annotation = ha, cluster_columns = FALSE,
                                    show_column_names = FALSE, cluster_rows = TRUE,
                                    column_title = "4.9-4.19 signatures, samples ordered by ZEB1",
                                    col = circlize::colorRamp2(c(-2, 0, 2), c("#3C5488", "white", "#E64B35"))))
      dev.off()
    }
  }, error = function(e) message("FAIL signature heatmap: ", conditionMessage(e)))
  # Pharmacotype replication of the same signatures
  scp <- score_sets(pexpr, sigs)
  scp$sample <- if (!is.null(rownames(scp)) && !identical(rownames(scp), as.character(seq_len(nrow(scp))))) rownames(scp) else colnames(pexpr)
  scp <- merge(labp[, .(sample, ZEB1_log, etp)], scp, by = "sample")
  rows2 <- lapply(intersect(names(sigs), names(scp)), function(nm) {
    x <- as.numeric(scp$ZEB1_log); y <- as.numeric(scp[[nm]])
    ok <- is.finite(x) & is.finite(y)
    if (sum(ok) < 20) return(NULL)
    ct <- suppressWarnings(cor.test(x[ok], y[ok], method = "spearman", exact = FALSE))
    data.table(signature = nm, rho = unname(ct$estimate), p = ct$p.value)
  })
  sigp <- rbindlist(rows2); sigp[, fdr := p.adjust(p, "BH")]
  fwrite(sigp, file.path(tab, "M4_4.32_pharmacotype_signature_vs_ZEB1.tsv"), sep = "\t")
  cmp <- merge(sigtab[, .(signature, rho_T = rho, fdr_T = fdr)],
               sigp[, .(signature, rho_P = rho, fdr_P = fdr)], by = "signature")
  p3 <- ggplot(cmp, aes(rho_T, rho_P, label = signature)) +
    geom_hline(yintercept = 0, linetype = 2, color = "grey50") +
    geom_vline(xintercept = 0, linetype = 2, color = "grey50") +
    geom_point(size = 3.2, shape = 21, fill = "#E64B35") +
    { if (has_ggrepel) ggrepel::geom_text_repel(size = 3) else geom_text(size = 3, vjust = -0.6) } +
    labs(title = "4.32  Signature–ZEB1 rho: TARGET vs Pharmacotype",
         x = "rho TARGET", y = "rho Pharmacotype") + theme_sci()
  save_both(p3, file.path(fig, "M4_4.32_signature_replication.pdf"), 7.2, 6.2)
})

# ---------------------------------------------------------------------
# 4.20-4.25 TF / regulon via C3 TFT GSEA + target-mean activity
# ---------------------------------------------------------------------
run("4.20-4.25 TF", function() {
  if (!nrow(res_t)) stop("no TFT gmt")
  bubble_gsea(res_t, "4.20  TF motif / TFT GSEA (C3)", file.path(fig, "M4_4.20_tft_bubble.pdf"), n = 18)
  # activity = mean z of leading-edge / set genes for selected TFs
  want_tf <- c("ZEB1", "ZEB2", "LMO2", "TAL1", "GATA3", "TCF7", "MYC", "NOTCH1", "SPI1", "MEF2C", "BCL11B")
  act <- list()
  for (tf in want_tf) {
    hits <- names(ts)[grepl(paste0("(^|_)", tf, "(_|$)"), names(ts), ignore.case = TRUE)]
    genes <- unique(unlist(ts[hits], use.names = FALSE))
    genes <- intersect(genes, rownames(expr))
    if (length(genes) < 8) next
    z <- t(scale(t(expr[genes, , drop = FALSE])))
    act[[tf]] <- colMeans(z, na.rm = TRUE)
  }
  if (!length(act)) stop("no TF sets matched")
  ad <- as.data.frame(act)
  ad$usi <- colnames(expr)
  ad <- merge(lab[, .(usi, ZEB1_log, subtype)], ad, by = "usi")
  fwrite(ad, file.path(tab, "M4_4.20_tf_activity.tsv"), sep = "\t")
  rows <- lapply(setdiff(names(ad), c("usi", "ZEB1_log", "subtype")), function(nm) {
    ct <- suppressWarnings(cor.test(ad$ZEB1_log, ad[[nm]], method = "spearman", exact = FALSE))
    data.table(tf = nm, rho = unname(ct$estimate), p = ct$p.value)
  })
  tftab <- rbindlist(rows); tftab[, fdr := p.adjust(p, "BH")]
  setorder(tftab, rho)
  tftab[, tf := factor(tf, levels = tf)]
  fwrite(tftab, file.path(tab, "M4_4.20_tf_activity_vs_ZEB1.tsv"), sep = "\t")
  p <- ggplot(tftab, aes(rho, tf, fill = rho > 0)) +
    geom_vline(xintercept = 0, linetype = 2) +
    geom_col(width = 0.7, color = "white") +
    scale_fill_manual(values = c("TRUE" = "#E64B35", "FALSE" = "#3C5488"), guide = "none") +
    labs(title = "4.20-4.24  Regulon activity vs ZEB1 expression",
         subtitle = "Target-mean of C3 TFT gene sets; not VIPER (package unavailable)",
         x = "Spearman rho vs ZEB1 mRNA", y = NULL) + theme_sci()
  save_both(p, file.path(fig, "M4_4.20_4.24_tf_activity_forest.pdf"), 7.8, 5.6)
  if ("ZEB1" %in% names(ad)) {
    p2 <- ggplot(ad, aes(ZEB1_log, ZEB1, fill = subtype)) +
      geom_smooth(method = "lm", formula = y ~ x, color = "grey25", fill = "#F3E5C8",
                  inherit.aes = FALSE, aes(x = ZEB1_log, y = ZEB1)) +
      geom_point(shape = 21, size = 2.2, color = "black", stroke = 0.25) +
      labs(title = "4.24  ZEB1 mRNA vs ZEB1 TFT activity",
           x = "ZEB1 expression", y = "ZEB1 regulon score", fill = NULL) +
      theme_sci() + theme(legend.position = "top")
    save_both(p2, file.path(fig, "M4_4.24_ZEB1_expr_vs_activity.pdf"), 7.2, 5.4)
  }
})

# ---------------------------------------------------------------------
# 4.26-4.32 correlation network + hierarchical modules
# ---------------------------------------------------------------------
run("4.26-4.32 network", function() {
  if (!has_igraph) stop("igraph missing")
  top <- st[fdr_rho < 0.05][order(-abs(rho_zeb1)), gene][1:80]
  top <- intersect(top, rownames(expr))
  sub <- expr[top, , drop = FALSE]
  cm <- cor(t(sub), method = "spearman", use = "pairwise.complete.obs")
  ut <- which(upper.tri(cm), arr.ind = TRUE)
  ed <- data.table(g1 = rownames(cm)[ut[, 1]], g2 = colnames(cm)[ut[, 2]], r = cm[ut])
  ed <- ed[is.finite(r) & abs(r) >= 0.55]
  if (nrow(ed) < 8) stop("too few coexpression edges")
  g <- igraph::graph_from_data_frame(ed[, .(g1, g2, r)], directed = FALSE)
  igraph::E(g)$color <- ifelse(ed$r > 0, "#E64B35", "#3C5488")
  igraph::E(g)$width <- abs(ed$r) * 2.2
  igraph::V(g)$color <- ifelse(st$rho_zeb1[match(igraph::V(g)$name, st$gene)] > 0, "#F39B7F", "#4DBBD5")
  pdf(file.path(fig, "M4_4.26_correlation_network.pdf"), width = 8.4, height = 8.0)
  set.seed(1)
  plot(g, vertex.size = 8, vertex.label.cex = 0.55, vertex.label.color = "black",
       layout = igraph::layout_with_fr(g), main = "4.26  Spearman network |rho|>=0.55 among top 80 ZEB1 genes")
  dev.off()
  png(file.path(fig, "M4_4.26_correlation_network.png"), width = 8.4 * 180, height = 8.0 * 180, res = 180)
  set.seed(1)
  plot(g, vertex.size = 8, vertex.label.cex = 0.55, vertex.label.color = "black",
       layout = igraph::layout_with_fr(g), main = "4.26  Spearman network |rho|>=0.55 among top 80 ZEB1 genes")
  dev.off()

  # hierarchical modules on top 400 genes (WGCNA substitute: signed cor + hclust)
  top400 <- st[fdr_rho < 0.05][order(-abs(rho_zeb1)), gene][1:400]
  top400 <- intersect(top400, rownames(expr))
  sub4 <- t(scale(t(expr[top400, , drop = FALSE])))
  d <- as.dist(1 - cor(t(sub4), method = "spearman", use = "pairwise.complete.obs"))
  # cluster genes (rows of sub4), not samples
  pcs <- stats::prcomp(sub4, rank. = 5, scale. = FALSE)$x
  set.seed(1)
  cl <- stats::kmeans(pcs, centers = 4, nstart = 25)$cluster
  mods <- data.table(gene = names(cl), module = paste0("M", cl))
  fwrite(mods, file.path(tab, "M4_4.27_gene_modules.tsv"), sep = "\t")
  me <- sapply(split(mods$gene, mods$module), function(gg) {
    colMeans(sub4[intersect(gg, rownames(sub4)), , drop = FALSE], na.rm = TRUE)
  })
  me <- as.data.frame(me)
  me$usi <- rownames(me)
  me <- merge(lab[, .(usi, ZEB1_log, subtype)], me, by = "usi")
  mnames <- grep("^M[0-9]", names(me), value = TRUE)
  rows <- lapply(mnames, function(nm) {
    x <- as.numeric(me$ZEB1_log); y <- as.numeric(me[[nm]])
    ok <- is.finite(x) & is.finite(y)
    if (sum(ok) < 20 || stats::sd(y[ok]) == 0) return(NULL)
    ct <- suppressWarnings(cor.test(x[ok], y[ok], method = "spearman", exact = FALSE))
    data.table(module = nm, n_genes = sum(mods$module == nm), rho = unname(ct$estimate), p = ct$p.value)
  })
  mt <- rbindlist(rows)
  if (!nrow(mt)) stop("no module vs ZEB1 tests")
  mt[, fdr := p.adjust(p, "BH")]
  fwrite(mt, file.path(tab, "M4_4.28_module_vs_ZEB1.tsv"), sep = "\t")
  p <- ggplot(mt, aes(rho, reorder(module, rho), fill = rho > 0)) +
    geom_vline(xintercept = 0, linetype = 2) +
    geom_col(width = 0.65, color = "white") +
    geom_text(aes(label = sprintf("n=%d  rho=%.2f", n_genes, rho)),
              hjust = ifelse(mt$rho >= 0, -0.1, 1.1), size = 3) +
    scale_fill_manual(values = c("TRUE" = "#E64B35", "FALSE" = "#3C5488"), guide = "none") +
    labs(title = "4.28  Co-expression modules vs ZEB1",
         subtitle = "k-means on gene PCs of top 400 ZEB1-correlated genes; WGCNA package not installed",
         x = "Spearman rho vs ZEB1", y = NULL) + theme_sci()
  save_both(p, file.path(fig, "M4_4.28_module_vs_ZEB1.pdf"), 7.6, 5.2)

  long <- melt(as.data.table(me), id.vars = c("usi", "ZEB1_log", "subtype"),
               measure.vars = mnames, variable.name = "module", value.name = "ME")
  p2 <- ggplot(long, aes(subtype, ME, fill = subtype)) +
    geom_boxplot(width = 0.62, outlier.size = 0.4, alpha = 0.9) +
    facet_wrap(~module, scales = "free_y", nrow = 2) +
    scale_fill_manual(values = c(TAL = "#E64B35", TLX = "#4DBBD5", HOXA = "#00A087",
                                 "LMO2/LYL1" = "#3C5488", NKX2_1 = "#F39B7F", Unknown = "#8491B4")) +
    labs(title = "4.29  Module eigengenes by Liu subtype", x = NULL, y = "Module mean z", fill = NULL) +
    theme_sci() + theme(legend.position = "none", axis.text.x = element_text(angle = 30, hjust = 1))
  save_both(p2, file.path(fig, "M4_4.29_module_by_subtype.pdf"), 9.0, 6.2)

  # hubs = highest intramodular connectivity in the 400-gene matrix
  cm4 <- cor(t(sub4), method = "spearman", use = "pairwise.complete.obs")
  hubs <- rbindlist(lapply(split(mods, by = "module"), function(dd) {
    gg <- intersect(dd$gene, rownames(cm4))
    if (length(gg) < 5) return(NULL)
    k <- rowSums(abs(cm4[gg, gg, drop = FALSE]), na.rm = TRUE)
    data.table(module = dd$module[1], gene = names(sort(k, decreasing = TRUE))[1:min(8, length(k))],
               k = as.numeric(sort(k, decreasing = TRUE)[1:min(8, length(k))]))
  }))
  fwrite(hubs, file.path(tab, "M4_4.30_hub_genes.tsv"), sep = "\t")
  p3 <- ggplot(hubs, aes(k, reorder(gene, k), fill = module)) +
    geom_col(width = 0.75, color = "white") +
    labs(title = "4.30  Intramodular hubs", x = "Sum |Spearman| within module", y = NULL) +
    theme_sci() + theme(legend.position = "right")
  save_both(p3, file.path(fig, "M4_4.30_hub_genes.pdf"), 8.0, 6.6)

  # 4.31 coexpression of hubs (PPI substitute)
  hg <- unique(hubs$gene)
  hg <- intersect(hg, rownames(expr))
  cmh <- cor(t(expr[hg, , drop = FALSE]), method = "spearman")
  dh <- melt(as.data.table(cmh, keep.rownames = "a"), id.vars = "a", variable.name = "b", value.name = "r")
  dh <- dh[a < as.character(b)]
  p4 <- ggplot(dh, aes(a, b, fill = r)) +
    geom_tile(color = "white") +
    scale_fill_gradient2(low = "#3C5488", mid = "white", high = "#E64B35", midpoint = 0, limits = c(-1, 1)) +
    labs(title = "4.31  Hub-hub coexpression (PPI database not used)", x = NULL, y = NULL, fill = "rho") +
    theme_sci() + theme(axis.text.x = element_text(angle = 45, hjust = 1), axis.line = element_blank())
  save_both(p4, file.path(fig, "M4_4.31_hub_coexpression.pdf"), 7.4, 6.6)

  # 4.32 module scores on Pharmacotype
  mep <- sapply(split(mods$gene, mods$module), function(gg) {
    gg <- intersect(gg, rownames(pexpr))
    if (length(gg) < 5) return(rep(NA_real_, ncol(pexpr)))
    z <- t(scale(t(pexpr[gg, , drop = FALSE])))
    colMeans(z, na.rm = TRUE)
  })
  mep <- as.data.frame(mep); mep$sample <- colnames(pexpr)
  mep <- merge(labp[, .(sample, ZEB1_log)], mep, by = "sample")
  rows <- lapply(intersect(mnames, names(mep)), function(nm) {
    x <- as.numeric(mep$ZEB1_log); y <- as.numeric(mep[[nm]])
    ok <- is.finite(x) & is.finite(y)
    if (sum(ok) < 20) return(NULL)
    ct <- suppressWarnings(cor.test(x[ok], y[ok], method = "spearman", exact = FALSE))
    data.table(module = nm, rho = unname(ct$estimate), p = ct$p.value)
  })
  mtp <- rbindlist(rows); mtp[, fdr := p.adjust(p, "BH")]
  cmp <- merge(mt[, .(module, rho_T = rho)], mtp[, .(module, rho_P = rho)], by = "module")
  fwrite(cmp, file.path(tab, "M4_4.32_module_replication.tsv"), sep = "\t")
  p5 <- ggplot(cmp, aes(rho_T, rho_P, label = module)) +
    geom_hline(yintercept = 0, linetype = 2, color = "grey50") +
    geom_vline(xintercept = 0, linetype = 2, color = "grey50") +
    geom_abline(slope = 1, intercept = 0, linetype = 3, color = "grey70") +
    geom_point(size = 4, shape = 21, fill = "#4DBBD5") +
    geom_text(vjust = -0.8, size = 3.3) +
    labs(title = "4.32  Module–ZEB1 rho replicated in Pharmacotype",
         x = "rho TARGET", y = "rho Pharmacotype") + theme_sci()
  save_both(p5, file.path(fig, "M4_4.32_module_replication.pdf"), 6.8, 6.0)
})

# drop the unnumbered leftover NES bar if still present
unlink(file.path(fig, c("M4_fgsea_h_all_v2023_2_Hs_symbols.pdf",
                        "M4_fgsea_h_all_v2023_2_Hs_symbols.png",
                        "M4_4.6_hallmark_gsea_nes.pdf",
                        "M4_4.6_hallmark_gsea_nes.png")))

cat("MODULE4 HIGHORDER DONE\n")
cat("FIG", length(list.files(fig, pattern = "\\.(pdf|png)$")), "\n")

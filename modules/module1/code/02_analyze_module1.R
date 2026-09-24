# Module 1: ZEB1 disease specificity in T-ALL
# Plotting style adapted from code1 (bioR11 violin, bioR07 boxplot,
# bioR40 ROC, bioR38 forest, bioR17 heatmap, bioR06 sorted boxplot).

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(pheatmap)
  library(reshape2)
})

root <- "/data-b/liangfuhua/projects/ZEB1_TALL_analysis/module1"
fig <- file.path(root, "figures")
tab <- file.path(root, "tables")
proc <- file.path(root, "processed")
dir.create(fig, showWarnings = FALSE, recursive = TRUE)
dir.create(tab, showWarnings = FALSE, recursive = TRUE)

has_pROC <- requireNamespace("pROC", quietly = TRUE)
has_ggpubr <- requireNamespace("ggpubr", quietly = TRUE)
has_metafor <- requireNamespace("metafor", quietly = TRUE)

theme_mod <- theme_bw(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, color = "black"),
    axis.text.y = element_text(color = "black"),
    plot.subtitle = element_text(size = 10),
    legend.position = "right",
    panel.grid.minor = element_blank()
  )

save_pdf <- function(plot, file, width, height) {
  pdf(file, width = width, height = height)
  print(plot)
  dev.off()
  png(sub("\\.pdf$", ".png", file), width = width * 150, height = height * 150, res = 150)
  print(plot)
  dev.off()
}

annot_two_group <- function(a, b) {
  hg <- hedges_g(a, b)
  p <- wilcox_p(a, b)
  sprintf("n=%d vs %d; P=%s; g=%.2f [%.2f, %.2f]",
          sum(is.finite(a)), sum(is.finite(b)), fmt_p(p), hg$g, hg$lo, hg$hi)
}

wilcox_p <- function(a, b) {
  a <- a[is.finite(a)]
  b <- b[is.finite(b)]
  if (length(a) < 3 || length(b) < 3) return(NA_real_)
  wilcox.test(a, b, alternative = "two.sided", exact = FALSE)$p.value
}

hedges_g <- function(a, b) {
  a <- a[is.finite(a)]
  b <- b[is.finite(b)]
  n1 <- length(a); n2 <- length(b)
  if (n1 < 3 || n2 < 3) return(list(g = NA_real_, se = NA_real_, lo = NA_real_, hi = NA_real_))
  s1 <- var(a); s2 <- var(b)
  sp <- sqrt(((n1 - 1) * s1 + (n2 - 1) * s2) / (n1 + n2 - 2))
  if (!is.finite(sp) || sp == 0) return(list(g = 0, se = 0, lo = 0, hi = 0))
  d <- (mean(a) - mean(b)) / sp
  J <- 1 - 3 / (4 * (n1 + n2) - 9)
  g <- J * d
  se <- sqrt((n1 + n2) / (n1 * n2) + d^2 / (2 * (n1 + n2))) * J
  list(g = g, se = se, lo = g - 1.96 * se, hi = g + 1.96 * se)
}

fmt_p <- function(p) {
  ifelse(is.na(p), "NA", ifelse(p < 0.001, "<0.001", sprintf("%.3f", p)))
}

bh <- function(p) p.adjust(p, method = "BH")

compare_row <- function(cohort, contrast, gene, a, b, n_a, n_b, group_a, group_b) {
  hg <- hedges_g(a, b)
  p <- wilcox_p(a, b)
  data.table(
    cohort = cohort, contrast = contrast, gene = gene,
    group_a = group_a, n_a = n_a, mean_a = mean(a, na.rm = TRUE),
    group_b = group_b, n_b = n_b, mean_b = mean(b, na.rm = TRUE),
    hedges_g = hg$g, se = hg$se, ci95_lo = hg$lo, ci95_hi = hg$hi,
    wilcox_p = p
  )
}

roc_auc <- function(label, score) {
  ok <- is.finite(score) & !is.na(label)
  label <- label[ok]; score <- score[ok]
  if (length(unique(label)) < 2) return(list(auc = NA_real_, lo = NA_real_, hi = NA_real_, spec = NA, sens = NA))
  if (has_pROC) {
    obj <- pROC::roc(label, score, quiet = TRUE, direction = "<")
    ci <- pROC::ci.auc(obj)
    return(list(auc = as.numeric(obj$auc), lo = ci[1], hi = ci[3], spec = obj$specificities, sens = obj$sensitivities, obj = obj))
  }
  # Mann-Whitney AUC
  pos <- score[label == 1]; neg <- score[label == 0]
  auc <- (sum(outer(pos, neg, ">")) + 0.5 * sum(outer(pos, neg, "=="))) / (length(pos) * length(neg))
  # direction: higher score should predict class 1. For ZEB1-low in T-ALL we may invert.
  list(auc = auc, lo = NA_real_, hi = NA_real_, spec = NA, sens = NA)
}

mw_auc <- function(pos, neg) {
  n1 <- length(pos); n0 <- length(neg)
  (sum(outer(pos, neg, ">")) + 0.5 * sum(outer(pos, neg, "=="))) / (n1 * n0)
}

plot_roc <- function(label, score, title, file) {
  ok <- is.finite(score) & !is.na(label)
  label <- as.integer(label[ok]); score <- as.numeric(score[ok])
  pos <- score[label == 1]; neg <- score[label == 0]
  auc_hi <- mw_auc(pos, neg)
  auc_lo <- mw_auc(-pos, -neg)
  use_high <- auc_hi >= auc_lo
  auc <- if (use_high) auc_hi else auc_lo
  sc <- if (use_high) score else -score
  set.seed(1)
  boots <- replicate(1000, {
    i <- sample.int(length(label), replace = TRUE)
    lab <- label[i]; s <- sc[i]
    if (length(unique(lab)) < 2) return(NA_real_)
    mw_auc(s[lab == 1], s[lab == 0])
  })
  ci <- quantile(boots, c(0.025, 0.975), na.rm = TRUE)
  o <- order(sc, decreasing = TRUE)
  y <- label[o]
  tp <- cumsum(y == 1); fp <- cumsum(y == 0)
  tpr <- c(0, tp / max(sum(y == 1), 1), 1)
  fpr <- c(0, fp / max(sum(y == 0), 1), 1)
  pdf(file, width = 5, height = 5)
  plot(fpr, tpr, type = "l", col = "red", lwd = 2, xlab = "1 - Specificity", ylab = "Sensitivity",
       main = sprintf("%s\nAUC=%.3f (95%% CI %.3f-%.3f)", title, auc, ci[1], ci[2]),
       xlim = c(0, 1), ylim = c(0, 1))
  abline(0, 1, lty = 2, col = "grey50")
  dev.off()
  png(sub("\\.pdf$", ".png", file), width = 750, height = 750, res = 150)
  plot(fpr, tpr, type = "l", col = "red", lwd = 2, xlab = "1 - Specificity", ylab = "Sensitivity",
       main = sprintf("%s\nAUC=%.3f (95%% CI %.3f-%.3f)", title, auc, ci[1], ci[2]),
       xlim = c(0, 1), ylim = c(0, 1))
  abline(0, 1, lty = 2, col = "grey50")
  dev.off()
  data.table(title = title, auc = auc, ci_lo = unname(ci[1]), ci_hi = unname(ci[2]),
             direction = if (use_high) "score_high_predicts_TALL" else "score_low_predicts_TALL",
             n_pos = sum(label == 1), n_neg = sum(label == 0))
}

# forest adapted from bioR38, but for Hedges g (vline at 0)
plot_forest_g <- function(dt, file, xlab = "Hedges' g (T-ALL minus comparator)") {
  dt <- dt[is.finite(hedges_g)]
  if (nrow(dt) == 0) return(invisible(NULL))
  dt <- dt[order(hedges_g)]
  n <- nrow(dt)
  nRow <- n + 1
  ylim <- c(1, nRow)
  pdf(file, width = 9, height = max(4, 0.42 * n + 2))
  layout(matrix(c(1, 2), nc = 2), width = c(4, 2.2))
  par(mar = c(4, 1, 1.5, 1))
  plot(1, xlim = c(0, 3), ylim = ylim, type = "n", axes = FALSE, xlab = "", ylab = "")
  lab <- paste0(dt$cohort, " | ", dt$contrast, " | ", dt$gene)
  text(0, n:1, lab, adj = 0, cex = 0.7)
  gtxt <- sprintf("%.2f (%.2f, %.2f)", dt$hedges_g, dt$ci95_lo, dt$ci95_hi)
  ptxt <- fmt_p(dt$wilcox_p)
  text(2.2, n:1, ptxt, adj = 1, cex = 0.7)
  text(2.2, n + 1, "P", cex = 0.75, font = 2, adj = 1)
  text(3, n:1, gtxt, adj = 1, cex = 0.65)
  text(3, n + 1, "g (95% CI)", cex = 0.75, font = 2, adj = 1)
  par(mar = c(4, 1, 1.5, 1), mgp = c(2, 0.5, 0))
  xr <- range(c(dt$ci95_lo, dt$ci95_hi, -0.2, 0.2), na.rm = TRUE)
  plot(1, xlim = xr, ylim = ylim, type = "n", axes = FALSE, ylab = "", xaxs = "i", xlab = xlab)
  arrows(dt$ci95_lo, n:1, dt$ci95_hi, n:1, angle = 90, code = 3, length = 0.03, col = "darkblue", lwd = 2)
  abline(v = 0, col = "black", lty = 2, lwd = 2)
  boxcolor <- ifelse(dt$hedges_g > 0, "red", "blue")
  points(dt$hedges_g, n:1, pch = 15, col = boxcolor, cex = 1.2)
  axis(1)
  dev.off()
}

dl_meta <- function(g, se) {
  ok <- is.finite(g) & is.finite(se) & se > 0
  g <- g[ok]; se <- se[ok]
  k <- length(g)
  if (k < 2) return(list(g = NA, se = NA, lo = NA, hi = NA, Q = NA, I2 = NA, k = k))
  w <- 1 / se^2
  gfix <- sum(w * g) / sum(w)
  Q <- sum(w * (g - gfix)^2)
  df <- k - 1
  C <- sum(w) - sum(w^2) / sum(w)
  tau2 <- max(0, (Q - df) / C)
  w2 <- 1 / (se^2 + tau2)
  gRE <- sum(w2 * g) / sum(w2)
  seRE <- sqrt(1 / sum(w2))
  I2 <- max(0, (Q - df) / Q) * 100
  list(g = gRE, se = seRE, lo = gRE - 1.96 * seRE, hi = gRE + 1.96 * seRE, Q = Q, I2 = I2, k = k, p_Q = pchisq(Q, df, lower.tail = FALSE))
}

# ---------------- load ----------------
mile <- fread(file.path(proc, "gse13159_keygenes.tsv"))
mile[, disease := factor(disease, levels = c("T-ALL", "B-ALL", "AML", "CML", "CLL", "MDS", "Normal BM", "Other"))]
mile <- mile[disease != "Other" | is.na(disease)]

# 1.1 sample flowchart / counts
cnt <- mile[, .N, by = .(leukemia_class, disease, sample_type)][order(-N)]
fwrite(cnt, file.path(tab, "M1_1.1_sample_counts.tsv"), sep = "\t")
p_flow <- ggplot(mile[!is.na(disease)], aes(x = disease, fill = disease)) +
  geom_bar() +
  geom_text(stat = "count", aes(label = after_stat(count)), vjust = -0.3, size = 3) +
  labs(title = "GSE13159 / MILE sample composition", x = "Disease class", y = "n samples") +
  theme_mod + theme(legend.position = "none")
save_pdf(p_flow, file.path(fig, "M1_1.1_sample_flowchart.pdf"), 8, 5)

# helper long for 3 genes
genes <- c("ZEB1", "ZEB2", "LMO2")
stopifnot(all(genes %in% names(mile)))

# 1.2 T-ALL vs B-ALL violin (bioR11)
vb <- mile[disease %in% c("T-ALL", "B-ALL")]
vb$Type <- factor(vb$disease, levels = c("B-ALL", "T-ALL"))
p_v <- ggplot(vb, aes(x = Type, y = ZEB1, fill = Type)) +
  geom_violin(trim = FALSE, alpha = 0.85) +
  geom_boxplot(width = 0.12, fill = "white", outlier.shape = NA) +
  labs(x = "Disease", y = "ZEB1 (GSE13159 series matrix, scaled intensity)",
       title = "GSE13159: ZEB1 in T-ALL vs B-ALL",
       subtitle = annot_two_group(vb[Type == "T-ALL", ZEB1], vb[Type == "B-ALL", ZEB1])) +
  theme_mod + theme(legend.position = "none", axis.text.x = element_text(angle = 0, hjust = 0.5))
save_pdf(p_v, file.path(fig, "M1_1.2_ZEB1_TALL_vs_BALL_violin.pdf"), 6.8, 5.3)

# 1.3 T-ALL vs AML
va <- mile[disease %in% c("T-ALL", "AML")]
va$Type <- factor(va$disease, levels = c("AML", "T-ALL"))
p_a <- ggplot(va, aes(x = Type, y = ZEB1, fill = Type)) +
  geom_violin(trim = FALSE, alpha = 0.85) +
  geom_boxplot(width = 0.12, fill = "white", outlier.shape = NA) +
  labs(x = "Disease", y = "ZEB1 (GSE13159 series matrix, scaled intensity)",
       title = "GSE13159: ZEB1 in T-ALL vs AML",
       subtitle = annot_two_group(va[Type == "T-ALL", ZEB1], va[Type == "AML", ZEB1])) +
  theme_mod + theme(legend.position = "none", axis.text.x = element_text(angle = 0, hjust = 0.5))
save_pdf(p_a, file.path(fig, "M1_1.3_ZEB1_TALL_vs_AML_violin.pdf"), 6.8, 5.3)

# 1.4 T-ALL vs CML/CLL/MDS boxplot (bioR07/09)
v4 <- mile[disease %in% c("T-ALL", "CML", "CLL", "MDS")]
v4$Type <- factor(v4$disease, levels = c("T-ALL", "CML", "CLL", "MDS"))
p4 <- ggplot(v4, aes(x = Type, y = ZEB1, color = Type)) +
  geom_boxplot(outlier.shape = NA, width = 0.6) +
  geom_jitter(width = 0.15, size = 0.4, alpha = 0.45) +
  labs(x = "Disease", y = "ZEB1 (GSE13159 series matrix, scaled intensity)",
       title = "GSE13159: ZEB1 in T-ALL vs CML/CLL/MDS") +
  theme_mod + theme(legend.position = "none")
save_pdf(p4, file.path(fig, "M1_1.4_ZEB1_TALL_vs_CML_CLL_MDS_box.pdf"), 6.5, 5)

# 1.5 T-ALL vs normal BM
vn <- mile[disease %in% c("T-ALL", "Normal BM")]
vn$Type <- factor(vn$disease, levels = c("Normal BM", "T-ALL"))
p5 <- ggplot(vn, aes(x = Type, y = ZEB1, color = Type)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.15, size = 0.6, alpha = 0.5) +
  labs(x = "Disease", y = "ZEB1 (GSE13159 series matrix, scaled intensity)",
       title = "GSE13159: ZEB1 in T-ALL vs normal bone marrow",
       subtitle = annot_two_group(vn[Type == "T-ALL", ZEB1], vn[Type == "Normal BM", ZEB1])) +
  theme_mod + theme(legend.position = "none", axis.text.x = element_text(angle = 0, hjust = 0.5))
save_pdf(p5, file.path(fig, "M1_1.5_ZEB1_TALL_vs_normalBM_box.pdf"), 5.5, 5)

# 1.6 ranking across all leukemias (bioR06 sorted boxplot)
vr <- mile[!is.na(disease)]
med <- vr[, .(med = median(ZEB1, na.rm = TRUE)), by = disease]
vr$disease <- factor(vr$disease, levels = med[order(-med)]$disease)
p6 <- ggplot(vr, aes(x = disease, y = ZEB1, color = disease)) +
  geom_boxplot(outlier.shape = NA) +
  labs(x = "Disease class", y = "ZEB1 (GSE13159 series matrix, scaled intensity)", title = "GSE13159: ZEB1 ranked by median across leukemias") +
  theme_mod + theme(legend.position = "none")
save_pdf(p6, file.path(fig, "M1_1.6_ZEB1_disease_ranking_box.pdf"), 8, 5.2)

# standardized Z ranking table
vr[, ZEB1_z := as.numeric(scale(ZEB1))]
rank_tab <- vr[, .(
  n = .N,
  mean = mean(ZEB1, na.rm = TRUE),
  sd = sd(ZEB1, na.rm = TRUE),
  median = median(ZEB1, na.rm = TRUE),
  mean_z = mean(ZEB1_z, na.rm = TRUE)
), by = disease][order(-mean_z)]
fwrite(rank_tab, file.path(tab, "M1_1.6_ZEB1_disease_ranking.tsv"), sep = "\t")

# pairwise stats vs T-ALL
tall <- mile[disease == "T-ALL", ZEB1]
stats <- rbindlist(lapply(setdiff(unique(as.character(mile$disease)), c("T-ALL", "Other", NA)), function(g) {
  b <- mile[disease == g, ZEB1]
  compare_row("GSE13159", paste0("T-ALL vs ", g), "ZEB1", tall, b, length(tall), length(b), "T-ALL", g)
}))

# 1.7 ZEB1 vs ZEB2 parallel
stats27 <- rbindlist(lapply(genes, function(gn) {
  rbindlist(lapply(c("B-ALL", "AML", "Normal BM"), function(g) {
    a <- mile[disease == "T-ALL"][[gn]]
    b <- mile[disease == g][[gn]]
    compare_row("GSE13159", paste0("T-ALL vs ", g), gn, a, b, length(a), length(b), "T-ALL", g)
  }))
}))
stats27[, fdr := bh(wilcox_p)]
fwrite(stats27, file.path(tab, "M1_1.7_ZEB1_ZEB2_LMO2_effects.tsv"), sep = "\t")
plot_forest_g(stats27, file.path(fig, "M1_1.7_ZEB1_ZEB2_LMO2_paired_forest.pdf"))

# 1.8 heatmap of 3 genes (class medians + sample heatmap)
mat_cls <- as.matrix(mile[!is.na(disease), lapply(.SD, median, na.rm = TRUE), by = disease, .SDcols = genes][, -1])
rownames(mat_cls) <- mile[!is.na(disease), lapply(.SD, median, na.rm = TRUE), by = disease, .SDcols = genes]$disease
pdf(file.path(fig, "M1_1.8_ZEB1_ZEB2_LMO2_class_heatmap.pdf"), width = 5, height = 5)
pheatmap(t(mat_cls), scale = "row", color = colorRampPalette(c("blue", "white", "red"))(50),
         fontsize = 10, main = "Median expression by disease class")
dev.off()

set.seed(1)
# downsample to <= 40 per class for sample heatmap
ss <- mile[!is.na(disease)][, .SD[sample(.N, min(.N, 40))], by = disease]
hm <- t(as.matrix(ss[, ..genes]))
colnames(hm) <- ss$gsm
ann <- data.frame(Disease = ss$disease, row.names = ss$gsm)
pdf(file.path(fig, "M1_1.8_ZEB1_ZEB2_LMO2_sample_heatmap.pdf"), width = 10, height = 4)
pheatmap(hm, annotation_col = ann, scale = "row", show_colnames = FALSE,
         color = colorRampPalette(c("blue", "white", "red"))(50), fontsize_row = 10)
dev.off()

# 1.9 ROC T-ALL vs B-ALL  (higher/lower ZEB1)
roc19 <- plot_roc(as.integer(vb$Type == "T-ALL"), vb$ZEB1,
                  "ZEB1: T-ALL vs B-ALL", file.path(fig, "M1_1.9_ROC_TALL_vs_BALL.pdf"))
# 1.10 ROC T-ALL vs non-T-ALL
non <- mile[!is.na(disease)]
roc10 <- plot_roc(as.integer(non$disease == "T-ALL"), non$ZEB1,
                  "ZEB1: T-ALL vs non-T-ALL", file.path(fig, "M1_1.10_ROC_TALL_vs_nonTALL.pdf"))
roc_tab <- rbind(roc19, roc10)
fwrite(roc_tab, file.path(tab, "M1_1.9_1.10_ROC.tsv"), sep = "\t")

# ---------------- TARGET 1.11: diagnosis primary, one sample per patient ----------------
tgt_f <- file.path(proc, "target_keygenes.tsv")
tgt <- if (file.exists(tgt_f)) fread(tgt_f) else data.table()
if (nrow(tgt)) {
  tgt <- tgt[grepl("^Primary", sample_type) & disease %in% c("T-ALL", "B-ALL") & is.finite(ZEB1)]
  tgt[, pref := fifelse(grepl("Bone Marrow", sample_type), 1L, 2L)]
  setorder(tgt, usi, pref)
  tgt <- tgt[!duplicated(usi)]
  tgt$Type <- factor(tgt$disease, levels = c("B-ALL", "T-ALL"))
  tgt[, ZEB1_log := log2(ZEB1 + 1)]
  p11 <- ggplot(tgt, aes(x = Type, y = ZEB1_log, fill = Type)) +
    geom_violin(trim = FALSE, alpha = 0.85) +
    geom_boxplot(width = 0.12, fill = "white", outlier.shape = NA) +
    labs(x = "Disease", y = "log2(TPM+1) ZEB1",
         title = "TARGET-ALL-P2 primary samples: ZEB1 in T-ALL vs B-ALL",
         subtitle = annot_two_group(tgt[disease == "T-ALL", ZEB1_log], tgt[disease == "B-ALL", ZEB1_log])) +
    theme_mod + theme(legend.position = "none", axis.text.x = element_text(angle = 0, hjust = 0.5))
  save_pdf(p11, file.path(fig, "M1_1.11_TARGET_ZEB1_violin.pdf"), 6.8, 5.3)
  stats_t <- compare_row("TARGET-ALL-P2", "T-ALL vs B-ALL", "ZEB1",
                         tgt[disease == "T-ALL", ZEB1_log],
                         tgt[disease == "B-ALL", ZEB1_log],
                         sum(tgt$disease == "T-ALL"), sum(tgt$disease == "B-ALL"), "T-ALL", "B-ALL")
} else {
  stats_t <- data.table()
}

# ---------------- Pharmacotype as 1.12 substitute for PanALL ----------------
ph_f <- file.path(proc, "pharmacotype_keygenes.tsv")
ph <- if (file.exists(ph_f)) fread(ph_f) else data.table()
note_panall <- "St. Jude PanALL public ProteinPaint matrix is B-ALL-only; T-ALL vs B-ALL independent RNA-seq uses St. Jude Pharmacotype."
writeLines(note_panall, file.path(tab, "M1_1.12_PanALL_note.txt"))
if (nrow(ph) && "ZEB1" %in% names(ph)) {
  ph <- ph[disease %in% c("T-ALL", "B-ALL") & is.finite(ZEB1)]
  if (nrow(ph)) {
    ph$Type <- factor(ph$disease, levels = c("B-ALL", "T-ALL"))
    # FPKM or already log? inspect range
    y <- if (max(ph$ZEB1, na.rm = TRUE) > 50) log2(ph$ZEB1 + 1) else ph$ZEB1
    ph$ZEB1_plot <- y
    p12 <- ggplot(ph, aes(x = Type, y = ZEB1_plot, fill = Type)) +
      geom_violin(trim = FALSE, alpha = 0.85) +
      geom_boxplot(width = 0.12, fill = "white", outlier.shape = NA) +
      labs(x = "Disease", y = "log2(FPKM+1) ZEB1",
           title = "St. Jude Pharmacotype: ZEB1 in T-ALL vs B-ALL",
           subtitle = annot_two_group(ph[disease == "T-ALL", ZEB1_plot], ph[disease == "B-ALL", ZEB1_plot])) +
      theme_mod + theme(legend.position = "none", axis.text.x = element_text(angle = 0, hjust = 0.5))
    save_pdf(p12, file.path(fig, "M1_1.12_Pharmacotype_ZEB1_violin.pdf"), 6.8, 5.3)
    a <- ph[disease == "T-ALL", ZEB1_plot]
    b <- ph[disease == "B-ALL", ZEB1_plot]
    stats_p <- compare_row("StJude_Pharmacotype", "T-ALL vs B-ALL", "ZEB1", a, b, length(a), length(b), "T-ALL", "B-ALL")
  } else stats_p <- data.table()
} else stats_p <- data.table()

# ---------------- 1.13-1.16 meta ----------------
meta_in <- rbindlist(list(
  compare_row("GSE13159", "T-ALL vs B-ALL", "ZEB1",
              mile[disease == "T-ALL", ZEB1], mile[disease == "B-ALL", ZEB1],
              sum(mile$disease == "T-ALL"), sum(mile$disease == "B-ALL"), "T-ALL", "B-ALL"),
  stats_t, stats_p
), fill = TRUE)
meta_in <- meta_in[is.finite(hedges_g)]
meta_in[, fdr := bh(wilcox_p)]
fwrite(meta_in, file.path(tab, "M1_1.13_cohort_effect_sizes.tsv"), sep = "\t")
plot_forest_g(meta_in, file.path(fig, "M1_1.13_cohort_effect_forest.pdf"))

ms <- dl_meta(meta_in$hedges_g, meta_in$se)
meta_sum <- data.table(
  model = "random_effects_DL",
  k = ms$k, hedges_g = ms$g, se = ms$se, ci95_lo = ms$lo, ci95_hi = ms$hi,
  Q = ms$Q, I2 = ms$I2, p_heterogeneity = ms$p_Q
)
# leave one out
loo <- rbindlist(lapply(seq_len(nrow(meta_in)), function(i) {
  sub <- meta_in[-i]
  m <- dl_meta(sub$hedges_g, sub$se)
  data.table(left_out = meta_in$cohort[i], k = m$k, hedges_g = m$g, ci95_lo = m$lo, ci95_hi = m$hi, I2 = m$I2)
}))
fwrite(meta_sum, file.path(tab, "M1_1.14_random_effects_meta.tsv"), sep = "\t")
fwrite(loo, file.path(tab, "M1_1.16_leave_one_cohort_out.tsv"), sep = "\t")
plot_forest_g(rbind(
  meta_in[, .(cohort, contrast = "T-ALL vs B-ALL", gene = "ZEB1", hedges_g, se, ci95_lo, ci95_hi, wilcox_p)],
  data.table(cohort = "Pooled RE", contrast = "T-ALL vs B-ALL", gene = "ZEB1",
             hedges_g = ms$g, se = ms$se, ci95_lo = ms$lo, ci95_hi = ms$hi, wilcox_p = NA_real_)
), file.path(fig, "M1_1.14_1.15_meta_forest.pdf"))
if (nrow(loo)) {
  loo2 <- data.table(cohort = paste0("LOO-", loo$left_out), contrast = "meta", gene = "ZEB1",
                     hedges_g = loo$hedges_g, se = (loo$ci95_hi - loo$ci95_lo) / 3.92,
                     ci95_lo = loo$ci95_lo, ci95_hi = loo$ci95_hi, wilcox_p = NA_real_)
  plot_forest_g(loo2, file.path(fig, "M1_1.16_leave_one_cohort_out_forest.pdf"))
}

# 1.17 child vs adult: TARGET age only (pediatric cohort). Adult T-ALL pending GSE142522.
if (nrow(tgt) && "Age at Diagnosis in Days" %in% names(tgt)) {
  tgt[, age_years := as.numeric(`Age at Diagnosis in Days`) / 365.25]
  fwrite(tgt[, .(n = .N, min_age = min(age_years, na.rm = TRUE), max_age = max(age_years, na.rm = TRUE)), by = disease],
         file.path(tab, "M1_1.17_TARGET_age_range.tsv"), sep = "\t")
  writeLines("GSE13159 has no age annotation. TARGET and Pharmacotype are pediatric. Adult T-ALL comparison requires GSE142522 (still downloading).",
             file.path(tab, "M1_1.17_age_limitation.txt"))
}

# 1.18 platform sensitivity: array vs RNA-seq as subgroups of the same contrast
plat <- copy(meta_in)
plat[, platform := ifelse(grepl("GSE13159", cohort), "microarray", "RNA-seq")]
fwrite(plat, file.path(tab, "M1_1.18_platform_effects.tsv"), sep = "\t")
plot_forest_g(plat[, .(cohort = paste(platform, cohort), contrast, gene, hedges_g, se, ci95_lo, ci95_hi, wilcox_p)],
              file.path(fig, "M1_1.18_platform_sensitivity_forest.pdf"))

# combined pairwise (1.2-1.5) stats
stats[, fdr := bh(wilcox_p)]
fwrite(stats, file.path(tab, "M1_1.2_1.5_pairwise_vs_TALL.tsv"), sep = "\t")
plot_forest_g(stats, file.path(fig, "M1_1.2_1.5_pairwise_forest.pdf"))

# omnibus Kruskal-Wallis for 1.4-like
kw <- kruskal.test(ZEB1 ~ disease, data = mile[disease %in% c("T-ALL", "CML", "CLL", "MDS")])
writeLines(capture.output(kw), file.path(tab, "M1_1.4_kruskal_wallis.txt"))

# session
writeLines(capture.output(sessionInfo()), file.path(tab, "sessionInfo.txt"))

# probe-level sensitivity: canonical ZEB1 212764_at
prb <- file.path(proc, "gse13159_keygene_probes.tsv")
if (file.exists(prb) && "212764_at" %in% names(fread(prb, nrows = 1))) {
  pd <- fread(prb)
  if (!"gsm" %in% names(pd)) setnames(pd, names(pd)[1], "gsm")
  m2 <- merge(mile[, .(gsm, disease)], pd[, .(gsm, `212764_at`)], by = "gsm")
  a <- m2[disease == "T-ALL", `212764_at`]
  sens <- rbindlist(lapply(c("B-ALL", "AML", "Normal BM"), function(g) {
    b <- m2[disease == g, `212764_at`]
    compare_row("GSE13159_212764_at", paste0("T-ALL vs ", g), "ZEB1", a, b, length(a), length(b), "T-ALL", g)
  }))
  fwrite(sens, file.path(tab, "M1_sensitivity_ZEB1_probe_212764_at.tsv"), sep = "\t")
}

cat("ANALYZE DONE\n")
cat("meta k=", ms$k, " g=", ms$g, " I2=", ms$I2, "\n")

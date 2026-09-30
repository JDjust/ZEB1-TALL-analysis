# Module 3 revision: independent subtype validation
# Recipes from module3 / code2: raincloud, lollipop, forest, heatmap.
# ZEB1 kept continuous. Within-cohort z used only for the cross-cohort heatmap.

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})
has_beeswarm <- requireNamespace("ggbeeswarm", quietly = TRUE)
has_cht <- requireNamespace("ComplexHeatmap", quietly = TRUE) && requireNamespace("circlize", quietly = TRUE)
if (has_beeswarm) library(ggbeeswarm)

root <- "/data-b/liangfuhua/projects/ZEB1_TALL_analysis/module3"
fig <- file.path(root, "figures")
tab <- file.path(root, "tables")
proc <- file.path(root, "processed")
dir.create(fig, FALSE, TRUE)
dir.create(tab, FALSE, TRUE)

pal_l1 <- c(
  Immature = "#3C5488", TLX = "#4DBBD5", TAL = "#E64B35",
  HOXA = "#00A087", NKX = "#F39B7F", Other = "#8491B4"
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
  pdf(file, width = width, height = height, useDingbats = FALSE)
  print(plot); dev.off()
  png(sub("\\.pdf$", ".png", file), width = width * 180, height = height * 180, res = 180)
  print(plot); dev.off()
}
fmt_p <- function(p) ifelse(is.na(p), "NA", ifelse(p < 0.001, "P < 0.001", sprintf("P = %.3f", p)))
p_star <- function(p) ifelse(is.na(p), "ns", ifelse(p < 0.001, "***", ifelse(p < 0.01, "**", ifelse(p < 0.05, "*", "ns"))))
add_points <- function(size = 0.5, alpha = 0.5, width = 0.16) {
  if (has_beeswarm) geom_quasirandom(size = size, alpha = alpha, width = width, color = "grey20")
  else geom_jitter(size = size, alpha = alpha, width = width, color = "grey20")
}
wilcox_p <- function(a, b) {
  a <- a[is.finite(a)]; b <- b[is.finite(b)]
  if (length(a) < 3 || length(b) < 3) return(NA_real_)
  wilcox.test(a, b, exact = FALSE)$p.value
}
hedges_g <- function(a, b) {
  a <- a[is.finite(a)]; b <- b[is.finite(b)]
  n1 <- length(a); n2 <- length(b)
  if (n1 < 3 || n2 < 3) return(list(g = NA_real_, se = NA_real_, lo = NA_real_, hi = NA_real_))
  sp <- sqrt(((n1 - 1) * var(a) + (n2 - 1) * var(b)) / (n1 + n2 - 2))
  if (!is.finite(sp) || sp == 0) return(list(g = 0, se = 0, lo = 0, hi = 0))
  d <- (mean(a) - mean(b)) / sp
  J <- 1 - 3 / (4 * (n1 + n2) - 9)
  g <- J * d
  se <- sqrt((n1 + n2) / (n1 * n2) + d^2 / (2 * (n1 + n2))) * J
  list(g = g, se = se, lo = g - 1.96 * se, hi = g + 1.96 * se)
}
dl_meta <- function(g, se) {
  ok <- is.finite(g) & is.finite(se) & se > 0
  g <- g[ok]; se <- se[ok]
  k <- length(g)
  if (k < 2) return(list(g = NA, se = NA, lo = NA, hi = NA, Q = NA, I2 = NA, k = k, p_Q = NA))
  w <- 1 / se^2
  gfix <- sum(w * g) / sum(w)
  Q <- sum(w * (g - gfix)^2)
  df <- k - 1
  C <- sum(w) - sum(w^2) / sum(w)
  tau2 <- max(0, (Q - df) / C)
  w2 <- 1 / (se^2 + tau2)
  gRE <- sum(w2 * g) / sum(w2)
  seRE <- sqrt(1 / sum(w2))
  I2 <- if (Q > 0) max(0, (Q - df) / Q) * 100 else 0
  list(g = gRE, se = seRE, lo = gRE - 1.96 * seRE, hi = gRE + 1.96 * seRE,
       Q = Q, I2 = I2, k = k, p_Q = pchisq(Q, df, lower.tail = FALSE))
}

dt <- fread(file.path(proc, "independent_keygenes.tsv"))
tall <- dt[disease == "T-ALL" & is.finite(ZEB1)]
tall[, level1 := factor(level1, levels = intersect(names(pal_l1), unique(level1)))]
tall[, ZEB1_z := scale(ZEB1), by = cohort]

# 3.1 composition
cnt <- tall[, .N, by = .(cohort, level1, level2)]
fwrite(cnt, file.path(tab, "M3_r2_label_counts.tsv"), sep = "\t")
sc <- tall[, .N, by = .(cohort, level1)]
p1 <- ggplot(sc, aes(x = cohort, y = N, fill = level1)) +
  geom_col(width = 0.72, color = "black", linewidth = 0.25) +
  scale_fill_manual(values = pal_l1, drop = FALSE) +
  labs(title = "Independent T-ALL cohorts, unified Level-1 subtype",
       subtitle = "Immature = LMO2/LYL1 / IMM / immature; author Level-2 kept in tables",
       x = NULL, y = "n patients", fill = "Level 1") +
  theme_sci() + theme(axis.text.x = element_text(angle = 20, hjust = 1))
save_both(p1, file.path(fig, "M3_r2_3.1_level1_counts.pdf"), 7.4, 5.0)

p1b <- ggplot(cnt[N > 0], aes(x = level1, y = N, fill = level2)) +
  geom_col(width = 0.8, color = "grey20", linewidth = 0.15) +
  facet_wrap(~ cohort, scales = "free_y") +
  labs(title = "Author labels (Level 2) inside unified Level 1",
       x = NULL, y = "n", fill = "Level 2") +
  theme_sci() + theme(axis.text.x = element_text(angle = 30, hjust = 1),
                      legend.position = "bottom", legend.text = element_text(size = 7))
save_both(p1b, file.path(fig, "M3_r2_3.2_label_alluvial_proxy.pdf"), 10.5, 7.2)

# 3.2 ZEB1 raincloud by cohort
plot_rain <- function(df, title) {
  ggplot(df, aes(x = level1, y = ZEB1, fill = level1)) +
    geom_violin(trim = FALSE, scale = "width", color = NA, alpha = 0.72, width = 0.88) +
    add_points() +
    geom_boxplot(width = 0.16, fill = "white", outlier.shape = NA, linewidth = 0.4) +
    scale_fill_manual(values = pal_l1, drop = FALSE) +
    labs(x = NULL, y = "ZEB1 (cohort native scale)", title = title) +
    theme_sci() + theme(legend.position = "none",
                        axis.text.x = element_text(angle = 25, hjust = 1))
}
for (co in unique(as.character(tall$cohort))) {
  d <- tall[cohort == co]
  fn <- sprintf("M3_r2_3.2_ZEB1_%s.pdf", gsub("[^A-Za-z0-9]+", "_", co))
  save_both(plot_rain(d, sprintf("%s: ZEB1 by unified subtype", co)),
            file.path(fig, fn), 7.2, 5.2)
}
p_facet <- ggplot(tall, aes(x = level1, y = ZEB1_z, fill = level1)) +
  geom_violin(trim = FALSE, scale = "width", color = NA, alpha = 0.7, width = 0.9) +
  add_points(size = 0.4, alpha = 0.4) +
  geom_boxplot(width = 0.16, fill = "white", outlier.shape = NA, linewidth = 0.35) +
  geom_hline(yintercept = 0, linetype = 2, color = "grey40", linewidth = 0.35) +
  facet_wrap(~ cohort, nrow = 2) +
  scale_fill_manual(values = pal_l1, drop = FALSE) +
  labs(title = "Within-cohort ZEB1 z-score by unified subtype",
       subtitle = "z computed inside each cohort so array and RNA-seq are not mixed on the raw scale",
       x = NULL, y = "ZEB1 z") +
  theme_sci() + theme(legend.position = "none",
                      axis.text.x = element_text(angle = 30, hjust = 1))
save_both(p_facet, file.path(fig, "M3_r2_3.2_ZEB1_z_facet.pdf"), 10.2, 7.0)

# vs-rest effects among labeled subtypes only (exclude Other/unknown, matching TARGET round 1)
genes <- c("ZEB1", "ZEB2", "LMO2")
subs <- c("Immature", "TLX", "TAL", "HOXA", "NKX")
rows <- list()
for (co in unique(as.character(tall$cohort))) {
  d <- tall[cohort == co & as.character(level1) != "Other"]
  for (gn in genes) {
    if (!gn %in% names(d) || all(!is.finite(d[[gn]]))) next
    for (s in subs) {
      a <- d[as.character(level1) == s][[gn]]
      b <- d[as.character(level1) != s][[gn]]
      if (sum(is.finite(a)) < 3 || sum(is.finite(b)) < 3) next
      hg <- hedges_g(a, b)
      rows[[length(rows) + 1]] <- data.table(
        cohort = co, gene = gn, subtype = s, contrast = paste(s, "vs rest"),
        n_a = sum(is.finite(a)), n_b = sum(is.finite(b)),
        mean_a = mean(a, na.rm = TRUE), mean_b = mean(b, na.rm = TRUE),
        hedges_g = hg$g, se = hg$se, ci95_lo = hg$lo, ci95_hi = hg$hi,
        wilcox_p = wilcox_p(a, b)
      )
    }
  }
}
eff <- rbindlist(rows)
eff[, fdr := p.adjust(wilcox_p, method = "BH"), by = .(cohort, gene)]
fwrite(eff, file.path(tab, "M3_r2_subtype_vs_rest.tsv"), sep = "\t")

zeb <- eff[gene == "ZEB1"]
zeb[, lab := paste(cohort, subtype)]
setorder(zeb, hedges_g)
zeb[, lab := factor(lab, levels = unique(lab))]
p_lol <- ggplot(zeb, aes(x = hedges_g, y = lab, fill = hedges_g < 0)) +
  geom_vline(xintercept = 0, linetype = 2, color = "grey40", linewidth = 0.4) +
  geom_segment(aes(x = 0, xend = hedges_g, yend = lab),
               color = "#B07D3A", linewidth = 0.9) +
  geom_errorbar(aes(xmin = ci95_lo, xmax = ci95_hi), width = 0.2, linewidth = 0.4) +
  geom_point(shape = 21, size = 3.4, color = "black", stroke = 0.35) +
  geom_text(aes(label = sprintf("g=%.2f %s", hedges_g, p_star(wilcox_p))),
            hjust = ifelse(zeb$hedges_g >= 0, -0.1, 1.1), size = 2.7) +
  scale_fill_manual(values = c("TRUE" = "#3C5488", "FALSE" = "#E64B35"), guide = "none") +
  labs(title = "ZEB1 subtype vs rest across independent cohorts",
       subtitle = "Negative g = lower ZEB1 in the named subtype; whiskers 95% CI",
       x = "Hedges' g", y = NULL) +
  theme_sci()
save_both(p_lol, file.path(fig, "M3_r2_3.5_ZEB1_lollipop.pdf"), 8.8, 7.4)

# effect heatmap
hm <- dcast(zeb, subtype ~ cohort, value.var = "hedges_g")
fwrite(hm, file.path(tab, "M3_r2_3.5_effect_matrix.tsv"), sep = "\t")
mat <- as.matrix(hm[, -1, with = FALSE])
rownames(mat) <- hm$subtype
if (has_cht) {
  library(ComplexHeatmap)
  library(circlize)
  col_fun <- colorRamp2(c(-1.2, 0, 1.2), c("#3C5488", "white", "#E64B35"))
  ht <- Heatmap(mat, name = "g", col = col_fun, cluster_rows = FALSE, cluster_columns = FALSE,
                cell_fun = function(j, i, x, y, w, h, fill) {
                  grid.text(sprintf("%.2f", mat[i, j]), x, y, gp = gpar(fontsize = 9))
                },
                column_title = "ZEB1 Hedges g (subtype vs rest)",
                row_names_gp = gpar(fontsize = 11), column_names_gp = gpar(fontsize = 10))
  pdf(file.path(fig, "M3_r2_3.5_effect_heatmap.pdf"), width = 7.2, height = 4.6)
  draw(ht); dev.off()
  png(file.path(fig, "M3_r2_3.5_effect_heatmap.png"), width = 7.2 * 180, height = 4.6 * 180, res = 180)
  draw(ht); dev.off()
} else {
  long <- melt(hm, id.vars = "subtype", variable.name = "cohort", value.name = "g")
  p_hm <- ggplot(long, aes(cohort, subtype, fill = g)) +
    geom_tile(color = "white") +
    geom_text(aes(label = sprintf("%.2f", g)), size = 3) +
    scale_fill_gradient2(low = "#3C5488", mid = "white", high = "#E64B35", midpoint = 0) +
    labs(title = "ZEB1 Hedges g (subtype vs rest)", x = NULL, y = NULL) +
    theme_sci() + theme(axis.text.x = element_text(angle = 25, hjust = 1))
  save_both(p_hm, file.path(fig, "M3_r2_3.5_effect_heatmap.pdf"), 7.2, 4.8)
}

# meta forests for Immature and TLX
plot_meta <- function(sub, title, file) {
  sub <- copy(sub)[is.finite(hedges_g)]
  if (nrow(sub) < 1) return(invisible(NULL))
  ms <- dl_meta(sub$hedges_g, sub$se)
  pooled <- data.table(cohort = "Pooled RE", hedges_g = ms$g, se = ms$se,
                       ci95_lo = ms$lo, ci95_hi = ms$hi, wilcox_p = NA_real_,
                       n_a = NA, n_b = NA)
  plotd <- rbind(sub[, .(cohort, hedges_g, se, ci95_lo, ci95_hi, wilcox_p, n_a, n_b)],
                 pooled, fill = TRUE)
  plotd[, lab := factor(cohort, levels = rev(cohort))]
  p <- ggplot(plotd, aes(x = hedges_g, y = lab)) +
    geom_vline(xintercept = 0, linetype = 2, color = "grey40") +
    geom_errorbar(aes(xmin = ci95_lo, xmax = ci95_hi), width = 0.18, linewidth = 0.45) +
    geom_point(aes(shape = cohort == "Pooled RE", fill = hedges_g < 0),
               size = 3.6, color = "black", stroke = 0.35) +
    scale_shape_manual(values = c("TRUE" = 23, "FALSE" = 21), guide = "none") +
    scale_fill_manual(values = c("TRUE" = "#3C5488", "FALSE" = "#E64B35"), guide = "none") +
    labs(title = title,
         subtitle = sprintf("Random-effects DL; k = %d; I2 = %.0f%%; pooled g = %.2f [%.2f, %.2f]",
                            ms$k, ifelse(is.finite(ms$I2), ms$I2, NA), ms$g, ms$lo, ms$hi),
         x = "Hedges' g (subtype minus rest)", y = NULL) +
    theme_sci()
  save_both(p, file, 7.6, max(4.2, 0.55 * nrow(plotd) + 1.8))
  invisible(ms)
}
meta_imm <- plot_meta(zeb[subtype == "Immature"],
                      "Immature / LMO2-like vs rest: ZEB1",
                      file.path(fig, "M3_r2_3.5_Immature_meta_forest.pdf"))
meta_tlx <- plot_meta(zeb[subtype == "TLX"],
                     "TLX vs rest: ZEB1",
                     file.path(fig, "M3_r2_3.5_TLX_meta_forest.pdf"))
meta_tab <- rbindlist(list(
  data.table(contrast = "Immature vs rest", k = meta_imm$k, hedges_g = meta_imm$g,
             ci95_lo = meta_imm$lo, ci95_hi = meta_imm$hi, I2 = meta_imm$I2, p_Q = meta_imm$p_Q),
  data.table(contrast = "TLX vs rest", k = meta_tlx$k, hedges_g = meta_tlx$g,
             ci95_lo = meta_tlx$lo, ci95_hi = meta_tlx$hi, I2 = meta_tlx$I2, p_Q = meta_tlx$p_Q)
), fill = TRUE)
fwrite(meta_tab, file.path(tab, "M3_r2_3.5_meta.tsv"), sep = "\t")

# 3.8 LMO2 sanity: Immature should be LMO2-high
lmo <- eff[gene == "LMO2" & subtype == "Immature"]
if (nrow(lmo)) {
  p_lmo <- ggplot(lmo, aes(x = hedges_g, y = reorder(cohort, hedges_g))) +
    geom_vline(xintercept = 0, linetype = 2) +
    geom_errorbar(aes(xmin = ci95_lo, xmax = ci95_hi), width = 0.18) +
    geom_point(shape = 21, fill = "#E64B35", size = 3.5) +
    labs(title = "LMO2 in Immature vs rest (sanity check)",
         x = "Hedges' g", y = NULL) + theme_sci()
  save_both(p_lmo, file.path(fig, "M3_r2_3.8_LMO2_Immature_forest.pdf"), 6.8, 4.4)
}

# GSE26713 vs normal BM
bm <- dt[cohort == "GSE26713" & disease == "normal"]
ta <- dt[cohort == "GSE26713" & disease == "T-ALL"]
if (nrow(bm) >= 3 && nrow(ta) >= 3) {
  hg <- hedges_g(ta$ZEB1, bm$ZEB1)
  fwrite(data.table(contrast = "T-ALL vs BM", n_tall = nrow(ta), n_bm = nrow(bm),
                    hedges_g = hg$g, ci95_lo = hg$lo, ci95_hi = hg$hi,
                    wilcox_p = wilcox_p(ta$ZEB1, bm$ZEB1)),
         file.path(tab, "M3_r2_GSE26713_vs_BM.tsv"), sep = "\t")
}

# within-cohort z ~ level1 (simple OLS, Immature as reference)
fit_rows <- list()
for (co in unique(as.character(tall$cohort))) {
  d <- tall[cohort == co & as.character(level1) != "Other"]
  if (length(unique(as.character(d$level1))) < 2) next
  d[, level1 := relevel(factor(as.character(level1)), ref = "Immature")]
  fit <- try(lm(ZEB1_z ~ level1, data = d), silent = TRUE)
  if (inherits(fit, "try-error")) next
  sm <- summary(fit)$coefficients
  fit_rows[[length(fit_rows) + 1]] <- data.table(
    cohort = co, term = rownames(sm), beta = sm[, 1], se = sm[, 2], p = sm[, 4]
  )
}
if (length(fit_rows)) fwrite(rbindlist(fit_rows), file.path(tab, "M3_r2_3.6_level1_lm.tsv"), sep = "\t")

writeLines(capture.output(sessionInfo()), file.path(tab, "M3_r2_sessionInfo.txt"))
cat("INDEPENDENT_DONE figures", length(list.files(fig, pattern = "^M3_r2_")), "\n")

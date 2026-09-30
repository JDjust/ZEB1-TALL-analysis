# Module 1 figure enrichment
# Plot recipes adapted from code2:
#   【33】 nested violin + box
#   【202】 violin/box + jitter
#   【216】/【295】 split violin
#   【191】/【81】 lollipop + dumbbell
#   【12】 multi-gene grouped violin
#   【6】 sorted boxplot
#   【78】 alluvial
#   【88】/【257】 forest
#   【90】/【91】 ROC / multi-ROC
#   【238】 beeswarm
#   【210】 heatmap
# Statistics are not recomputed except display-only ROC curves.

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(reshape2)
})

has_beeswarm <- requireNamespace("ggbeeswarm", quietly = TRUE)
has_alluvial <- requireNamespace("ggalluvial", quietly = TRUE)
has_cowplot <- requireNamespace("cowplot", quietly = TRUE)
has_cht <- requireNamespace("ComplexHeatmap", quietly = TRUE)
has_circlize <- requireNamespace("circlize", quietly = TRUE)
if (has_beeswarm) library(ggbeeswarm)

root <- "/data-b/liangfuhua/projects/ZEB1_TALL_analysis/module1"
fig <- file.path(root, "figures")
tab <- file.path(root, "tables")
proc <- file.path(root, "processed")
dir.create(fig, showWarnings = FALSE, recursive = TRUE)

pal_dz <- c(
  "T-ALL" = "#E64B35",
  "B-ALL" = "#4DBBD5",
  "AML" = "#00A087",
  "CML" = "#3C5488",
  "CLL" = "#F39B7F",
  "MDS" = "#8491B4",
  "Normal BM" = "#91D1C2"
)
pal_gene <- c(ZEB1 = "#E64B35", ZEB2 = "#3C5488", LMO2 = "#F39B7F")
pal_tb <- c("B-ALL" = "#4DBBD5", "T-ALL" = "#E64B35")

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
      strip.text = element_text(face = "bold", color = "black")
    )
}

save_both <- function(plot, file, width, height) {
  pdf(file, width = width, height = height, useDingbats = FALSE)
  print(plot)
  dev.off()
  png(sub("\\.pdf$", ".png", file), width = width * 180, height = height * 180, res = 180)
  print(plot)
  dev.off()
}

fmt_p <- function(p) {
  ifelse(is.na(p), "NA", ifelse(p < 0.001, "P < 0.001", sprintf("P = %.3f", p)))
}

p_star <- function(p) {
  ifelse(is.na(p), "ns", ifelse(p < 0.001, "***", ifelse(p < 0.01, "**", ifelse(p < 0.05, "*", "ns"))))
}

add_points <- function(mapping = NULL, size = 0.45, alpha = 0.45, width = 0.18) {
  if (has_beeswarm) {
    geom_quasirandom(mapping = mapping, size = size, alpha = alpha, width = width, color = "grey20")
  } else {
    geom_jitter(mapping = mapping, size = size, alpha = alpha, width = width, color = "grey20")
  }
}

# code2 【295】 split violin, plyr-free
GeomSplitViolin <- ggproto("GeomSplitViolin", GeomViolin,
  draw_group = function(self, data, ..., draw_quantiles = NULL) {
    data <- transform(data,
                      xminv = x - violinwidth * (x - xmin),
                      xmaxv = x + violinwidth * (xmax - x))
    grp <- data[1, "group"]
    newdata <- transform(data, x = if (grp %% 2 == 1) xminv else xmaxv)
    ord <- if (grp %% 2 == 1) order(newdata$y) else order(-newdata$y)
    newdata <- newdata[ord, ]
    newdata <- rbind(newdata[1, ], newdata, newdata[nrow(newdata), ], newdata[1, ])
    newdata[c(1, nrow(newdata) - 1, nrow(newdata)), "x"] <- round(newdata[1, "x"])
    ggplot2:::ggname("geom_split_violin", GeomPolygon$draw_panel(newdata, ...))
  }
)

geom_split_violin <- function(mapping = NULL, data = NULL, stat = "ydensity",
                              position = "identity", ..., trim = FALSE, scale = "width",
                              na.rm = FALSE, show.legend = NA, inherit.aes = TRUE) {
  layer(data = data, mapping = mapping, stat = stat, geom = GeomSplitViolin,
        position = position, show.legend = show.legend, inherit.aes = inherit.aes,
        params = list(trim = trim, scale = scale, na.rm = na.rm, ...))
}

hedges_g <- function(a, b) {
  a <- a[is.finite(a)]; b <- b[is.finite(b)]
  n1 <- length(a); n2 <- length(b)
  sp <- sqrt(((n1 - 1) * var(a) + (n2 - 1) * var(b)) / (n1 + n2 - 2))
  d <- (mean(a) - mean(b)) / sp
  J <- 1 - 3 / (4 * (n1 + n2) - 9)
  g <- J * d
  se <- sqrt((n1 + n2) / (n1 * n2) + d^2 / (2 * (n1 + n2))) * J
  list(g = g, se = se, lo = g - 1.96 * se, hi = g + 1.96 * se)
}

wilcox_p <- function(a, b) {
  wilcox.test(a[is.finite(a)], b[is.finite(b)], alternative = "two.sided", exact = FALSE)$p.value
}

annot_two <- function(a, b) {
  hg <- hedges_g(a, b)
  sprintf("n = %d vs %d; %s; g = %.2f [%.2f, %.2f]",
          sum(is.finite(a)), sum(is.finite(b)), fmt_p(wilcox_p(a, b)), hg$g, hg$lo, hg$hi)
}

mw_auc <- function(pos, neg) {
  (sum(outer(pos, neg, ">")) + 0.5 * sum(outer(pos, neg, "=="))) / (length(pos) * length(neg))
}

roc_df <- function(label, score, name) {
  ok <- is.finite(score) & !is.na(label)
  label <- as.integer(label[ok]); score <- as.numeric(score[ok])
  pos <- score[label == 1]; neg <- score[label == 0]
  auc_hi <- mw_auc(pos, neg)
  auc_lo <- mw_auc(-pos, -neg)
  use_high <- auc_hi >= auc_lo
  auc <- if (use_high) auc_hi else auc_lo
  sc <- if (use_high) score else -score
  set.seed(1)
  boots <- replicate(800, {
    i <- sample.int(length(label), replace = TRUE)
    lab <- label[i]; s <- sc[i]
    if (length(unique(lab)) < 2) return(NA_real_)
    mw_auc(s[lab == 1], s[lab == 0])
  })
  ci <- quantile(boots, c(0.025, 0.975), na.rm = TRUE)
  o <- order(sc, decreasing = TRUE)
  y <- label[o]
  tp <- cumsum(y == 1); fp <- cumsum(y == 0)
  data.table(
    name = name,
    fpr = c(0, fp / max(sum(y == 0), 1), 1),
    tpr = c(0, tp / max(sum(y == 1), 1), 1),
    auc = auc, lo = unname(ci[1]), hi = unname(ci[2]),
    direction = if (use_high) "high" else "low"
  )
}

plot_rain <- function(df, x, y, fill_var, pal, ylab, title, subtitle = NULL, angle = 0) {
  p <- ggplot(df, aes(x = .data[[x]], y = .data[[y]], fill = .data[[fill_var]])) +
    geom_violin(trim = FALSE, scale = "width", color = NA, alpha = 0.72, width = 0.88) +
    add_points(size = 0.42, alpha = 0.38, width = 0.16) +
    geom_boxplot(width = 0.16, fill = "white", outlier.shape = NA, linewidth = 0.4) +
    scale_fill_manual(values = pal, drop = FALSE) +
    labs(x = NULL, y = ylab, title = title, subtitle = subtitle) +
    theme_sci() +
    theme(legend.position = "none",
          axis.text.x = element_text(angle = angle, hjust = if (angle > 0) 1 else 0.5))
  p
}

add_bracket <- function(p, x1, x2, y, label) {
  p +
    annotate("segment", x = x1, xend = x2, y = y, yend = y, linewidth = 0.4) +
    annotate("text", x = (x1 + x2) / 2, y = y, label = label, vjust = -0.35, size = 3.6)
}

# ---------------- load ----------------
mile <- fread(file.path(proc, "gse13159_keygenes.tsv"))
mile <- mile[disease %in% names(pal_dz)]
mile[, disease := factor(disease, levels = names(pal_dz))]
genes <- c("ZEB1", "ZEB2", "LMO2")

tgt <- if (file.exists(file.path(proc, "target_keygenes.tsv"))) fread(file.path(proc, "target_keygenes.tsv")) else data.table()
if (nrow(tgt)) {
  tgt <- tgt[grepl("^Primary", sample_type) & disease %in% c("T-ALL", "B-ALL") & is.finite(ZEB1)]
  tgt[, pref := fifelse(grepl("Bone Marrow", sample_type), 1L, 2L)]
  setorder(tgt, usi, pref)
  tgt <- tgt[!duplicated(usi)]
  tgt[, ZEB1_log := log2(ZEB1 + 1)]
  tgt[, ZEB2_log := log2(ZEB2 + 1)]
  tgt[, LMO2_log := log2(LMO2 + 1)]
  tgt[, Type := factor(disease, levels = c("B-ALL", "T-ALL"))]
}

ph <- if (file.exists(file.path(proc, "pharmacotype_keygenes.tsv"))) fread(file.path(proc, "pharmacotype_keygenes.tsv")) else data.table()
if (nrow(ph) && "ZEB1" %in% names(ph)) {
  ph <- ph[disease %in% c("T-ALL", "B-ALL") & is.finite(ZEB1)]
  log_needed <- max(ph$ZEB1, na.rm = TRUE) > 50
  for (gn in genes) {
    if (gn %in% names(ph)) ph[[paste0(gn, "_plot")]] <- if (log_needed) log2(ph[[gn]] + 1) else ph[[gn]]
  }
  ph[, Type := factor(disease, levels = c("B-ALL", "T-ALL"))]
}

stats_pair <- fread(file.path(tab, "M1_1.2_1.5_pairwise_vs_TALL.tsv"))
stats_gene <- fread(file.path(tab, "M1_1.7_ZEB1_ZEB2_LMO2_effects.tsv"))
meta_in <- fread(file.path(tab, "M1_1.13_cohort_effect_sizes.tsv"))
meta_sum <- fread(file.path(tab, "M1_1.14_random_effects_meta.tsv"))
loo <- fread(file.path(tab, "M1_1.16_leave_one_cohort_out.tsv"))

# =====================================================================
# 1.1 sample composition: donut + stacked alluvial-like bar (code2 78/79)
# =====================================================================
cnt <- mile[, .N, by = disease]
cnt[, pct := 100 * N / sum(N)]
cnt[, lab := sprintf("%s\n%d (%.0f%%)", disease, N, pct)]
p_bar <- ggplot(cnt, aes(x = disease, y = N, fill = disease)) +
  geom_col(width = 0.72, color = "white", linewidth = 0.3) +
  geom_text(aes(label = sprintf("%d\n%.0f%%", N, pct)), vjust = -0.15, size = 3.1, lineheight = 0.9) +
  scale_fill_manual(values = pal_dz) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.18))) +
  labs(title = "GSE13159 / MILE sample composition", x = NULL, y = "n samples") +
  theme_sci() + theme(legend.position = "none")
save_both(p_bar, file.path(fig, "M1_1.1_sample_flowchart.pdf"), 8.2, 5.2)

p_donut <- ggplot(cnt, aes(x = 2, y = N, fill = disease)) +
  geom_col(width = 0.85, color = "white", linewidth = 0.4) +
  coord_polar(theta = "y", start = 0) +
  xlim(0.4, 2.6) +
  scale_fill_manual(values = pal_dz) +
  annotate("text", x = 0.4, y = 0, label = sprintf("n = %d", sum(cnt$N)), size = 4, fontface = "bold") +
  labs(title = "MILE disease mix", fill = NULL) +
  theme_void(base_size = 11) +
  theme(plot.title = element_text(face = "bold", hjust = 0.5), legend.position = "right")
save_both(p_donut, file.path(fig, "M1_1.1_sample_donut.pdf"), 7.2, 5.2)

if (has_alluvial) {
  flow <- fread(file.path(tab, "M1_1.1_sample_counts.tsv"))
  flow <- flow[disease %in% names(pal_dz)]
  flow[, lineage := fcase(
    disease == "T-ALL", "T lymphoid",
    disease %in% c("B-ALL", "CLL"), "B lymphoid",
    disease %in% c("AML", "CML", "MDS"), "Myeloid",
    disease == "Normal BM", "Normal",
    default = "Other"
  )]
  p_al <- ggplot(flow, aes(y = N, axis1 = lineage, axis2 = disease, axis3 = sample_type)) +
    ggalluvial::geom_alluvium(aes(fill = disease), width = 1 / 8, alpha = 0.75, knot.pos = 0.2) +
    ggalluvial::geom_stratum(width = 1 / 8, fill = "grey97", color = "grey35", linewidth = 0.3) +
    geom_text(stat = "stratum", aes(label = after_stat(stratum)), size = 3) +
    scale_fill_manual(values = pal_dz) +
    scale_x_discrete(limits = c("Lineage", "Disease", "Tissue"), expand = c(0.05, 0.05)) +
    labs(title = "MILE sample flow: lineage to disease to tissue", y = "n samples") +
    theme_sci() + theme(legend.position = "none", axis.line.y = element_blank(), axis.ticks.y = element_blank())
  save_both(p_al, file.path(fig, "M1_1.1_sample_alluvial.pdf"), 9.5, 6)
}

# =====================================================================
# 1.2 / 1.3 / 1.5 two-group raincloud (code2 33 + 202)
# =====================================================================
vb <- mile[disease %in% c("T-ALL", "B-ALL")]
vb[, Type := factor(disease, levels = c("B-ALL", "T-ALL"))]
p12 <- plot_rain(vb, "Type", "ZEB1", "Type", pal_tb,
                 "ZEB1 (GSE13159 scaled intensity)",
                 "GSE13159: ZEB1 in T-ALL vs B-ALL",
                 annot_two(vb[Type == "T-ALL", ZEB1], vb[Type == "B-ALL", ZEB1]))
p12 <- add_bracket(p12, 1, 2, max(vb$ZEB1, na.rm = TRUE) * 1.04,
                   p_star(wilcox_p(vb[Type == "T-ALL", ZEB1], vb[Type == "B-ALL", ZEB1])))
save_both(p12, file.path(fig, "M1_1.2_ZEB1_TALL_vs_BALL_violin.pdf"), 5.6, 5.4)

va <- mile[disease %in% c("T-ALL", "AML")]
va[, Type := factor(disease, levels = c("AML", "T-ALL"))]
p13 <- plot_rain(va, "Type", "ZEB1", "Type", c("AML" = pal_dz[["AML"]], "T-ALL" = pal_dz[["T-ALL"]]),
                 "ZEB1 (GSE13159 scaled intensity)",
                 "GSE13159: ZEB1 in T-ALL vs AML",
                 annot_two(va[Type == "T-ALL", ZEB1], va[Type == "AML", ZEB1]))
p13 <- add_bracket(p13, 1, 2, max(va$ZEB1, na.rm = TRUE) * 1.04,
                   p_star(wilcox_p(va[Type == "T-ALL", ZEB1], va[Type == "AML", ZEB1])))
save_both(p13, file.path(fig, "M1_1.3_ZEB1_TALL_vs_AML_violin.pdf"), 5.6, 5.4)

vn <- mile[disease %in% c("T-ALL", "Normal BM")]
vn[, Type := factor(disease, levels = c("Normal BM", "T-ALL"))]
p15 <- plot_rain(vn, "Type", "ZEB1", "Type",
                 c("Normal BM" = pal_dz[["Normal BM"]], "T-ALL" = pal_dz[["T-ALL"]]),
                 "ZEB1 (GSE13159 scaled intensity)",
                 "GSE13159: ZEB1 in T-ALL vs normal bone marrow",
                 annot_two(vn[Type == "T-ALL", ZEB1], vn[Type == "Normal BM", ZEB1]))
p15 <- add_bracket(p15, 1, 2, max(vn$ZEB1, na.rm = TRUE) * 1.04,
                   p_star(wilcox_p(vn[Type == "T-ALL", ZEB1], vn[Type == "Normal BM", ZEB1])))
save_both(p15, file.path(fig, "M1_1.5_ZEB1_TALL_vs_normalBM_box.pdf"), 5.6, 5.4)

# 1.4 four-group raincloud
v4 <- mile[disease %in% c("T-ALL", "CML", "CLL", "MDS")]
v4[, Type := factor(disease, levels = c("T-ALL", "CML", "CLL", "MDS"))]
p14 <- plot_rain(v4, "Type", "ZEB1", "Type", pal_dz,
                 "ZEB1 (GSE13159 scaled intensity)",
                 "GSE13159: ZEB1 in T-ALL vs CML / CLL / MDS")
save_both(p14, file.path(fig, "M1_1.4_ZEB1_TALL_vs_CML_CLL_MDS_box.pdf"), 6.8, 5.3)

# =====================================================================
# 1.6 ranked raincloud (code2 6 + 33)
# =====================================================================
vr <- copy(mile)
med <- vr[, .(med = median(ZEB1, na.rm = TRUE), n = .N), by = disease]
setorder(med, -med)
vr[, disease := factor(disease, levels = as.character(med$disease))]
vr <- merge(vr, med[, .(disease, n)], by = "disease")
vr[, xlab := factor(sprintf("%s\n(n=%d)", as.character(disease), n),
                    levels = sprintf("%s\n(n=%d)", as.character(med$disease), med$n))]
p16 <- plot_rain(vr, "xlab", "ZEB1", "disease", pal_dz,
                 "ZEB1 (GSE13159 scaled intensity)",
                 "GSE13159: ZEB1 ranked by median across leukemias")
save_both(p16, file.path(fig, "M1_1.6_ZEB1_disease_ranking_box.pdf"), 8.4, 5.6)

# =====================================================================
# 1.2-1.5 lollipop of Hedges g (code2 81 / 191)
# =====================================================================
lp <- copy(stats_pair)
lp[, lab := sub("^T-ALL vs ", "", contrast)]
lp[, lab := factor(lab, levels = lab[order(hedges_g)])]
p_lp <- ggplot(lp, aes(x = hedges_g, y = lab)) +
  geom_vline(xintercept = 0, linetype = 2, color = "grey40", linewidth = 0.4) +
  geom_segment(aes(x = 0, xend = hedges_g, yend = lab), color = "#B07D3A", linewidth = 1.05) +
  geom_errorbar(aes(xmin = ci95_lo, xmax = ci95_hi), width = 0.18, color = "grey25", linewidth = 0.45) +
  geom_point(aes(fill = hedges_g > 0.5), shape = 21, size = 4.2, color = "black", stroke = 0.4) +
  geom_text(aes(label = sprintf("g = %.2f", hedges_g)), hjust = -0.15, size = 3.1, nudge_y = 0.22) +
  scale_fill_manual(values = c("TRUE" = "#E64B35", "FALSE" = "#4DBBD5"), guide = "none") +
  labs(title = "T-ALL vs each comparator: Hedges' g for ZEB1",
       subtitle = "Positive g = higher ZEB1 in T-ALL; whiskers = 95% CI",
       x = "Hedges' g (T-ALL minus comparator)", y = NULL) +
  theme_sci() +
  coord_cartesian(xlim = c(min(lp$ci95_lo) - 0.15, max(lp$ci95_hi) + 0.55))
save_both(p_lp, file.path(fig, "M1_1.2_1.5_effect_lollipop.pdf"), 7.2, 4.8)

# =====================================================================
# 1.7 split violin ZEB1 vs ZEB2 (code2 216 / 295)
# =====================================================================
long_zz <- melt(mile, id.vars = "disease", measure.vars = c("ZEB1", "ZEB2"),
                variable.name = "Gene", value.name = "Expression")
p_split <- ggplot(long_zz, aes(x = disease, y = Expression, fill = Gene)) +
  geom_split_violin(colour = NA, alpha = 0.9, scale = "width") +
  stat_summary(fun = median, geom = "point",
               position = position_dodge(width = 0.35),
               color = "white", size = 1.7) +
  scale_fill_manual(values = pal_gene[c("ZEB1", "ZEB2")]) +
  labs(title = "ZEB1 and ZEB2 occupy opposite disease ranks",
       subtitle = "Split violin: left/right halves are the two genes in the same samples",
       x = NULL, y = "Scaled intensity (GSE13159)", fill = NULL) +
  theme_sci() +
  theme(axis.text.x = element_text(angle = 30, hjust = 1), legend.position = "top")
save_both(p_split, file.path(fig, "M1_1.7_ZEB1_ZEB2_splitviolin.pdf"), 8.4, 5.4)

# 1.7 grouped multi-gene violin T vs B (code2 12)
long_tb <- melt(vb, id.vars = "Type", measure.vars = genes,
                variable.name = "Gene", value.name = "Expression")
p_mg <- ggplot(long_tb, aes(x = Gene, y = Expression, fill = Type)) +
  geom_violin(trim = FALSE, scale = "width", color = NA, alpha = 0.7,
              position = position_dodge(width = 0.82), width = 0.78) +
  geom_boxplot(width = 0.14, outlier.shape = NA, fill = "white",
               position = position_dodge(width = 0.82), linewidth = 0.35) +
  scale_fill_manual(values = pal_tb) +
  labs(title = "T-ALL vs B-ALL: ZEB1 / ZEB2 / LMO2",
       subtitle = "Same MILE samples; ZEB1 up in T-ALL, ZEB2 down, LMO2 modestly down",
       x = NULL, y = "Scaled intensity (GSE13159)", fill = NULL) +
  theme_sci() + theme(legend.position = "top")
save_both(p_mg, file.path(fig, "M1_1.7_multigene_violin_TvsB.pdf"), 7.2, 5.3)

# 1.7 dumbbell of median ZEB1 vs ZEB2 (code2 191)
medz <- mile[, .(ZEB1 = median(ZEB1, na.rm = TRUE), ZEB2 = median(ZEB2, na.rm = TRUE)), by = disease]
medz[, disease := factor(disease, levels = medz[order(ZEB1)]$disease)]
med_long <- melt(medz, id.vars = "disease", measure.vars = c("ZEB1", "ZEB2"),
                 variable.name = "Gene", value.name = "median")
p_db <- ggplot() +
  geom_segment(data = medz, aes(x = ZEB1, xend = ZEB2, y = disease, yend = disease),
               color = "#C4B48A", linewidth = 1.4) +
  geom_point(data = med_long, aes(x = median, y = disease, color = Gene), size = 3.6) +
  scale_color_manual(values = pal_gene[c("ZEB1", "ZEB2")]) +
  labs(title = "Median ZEB1 vs ZEB2 by disease",
       subtitle = "T-ALL is the only class with high ZEB1 and low ZEB2.",
       x = "Median scaled intensity", y = NULL, color = NULL) +
  theme_sci() + theme(legend.position = "top")
save_both(p_db, file.path(fig, "M1_1.7_dumbbell_ZEB1_ZEB2.pdf"), 7.2, 4.8)

# 1.7 forest ggplot
fg <- copy(stats_gene)
fg[, lab := paste0(gene, " | ", sub("^T-ALL vs ", "", contrast))]
fg[, lab := factor(lab, levels = lab[order(hedges_g)])]
p_f7 <- ggplot(fg, aes(x = hedges_g, y = lab, color = gene)) +
  geom_vline(xintercept = 0, linetype = 2, color = "grey40") +
  geom_errorbar(aes(xmin = ci95_lo, xmax = ci95_hi), width = 0.18, linewidth = 0.55) +
  geom_point(size = 2.8) +
  scale_color_manual(values = pal_gene) +
  labs(title = "T-ALL vs comparators: three genes",
       x = "Hedges' g (T-ALL minus comparator)", y = NULL, color = NULL) +
  theme_sci() + theme(legend.position = "top")
save_both(p_f7, file.path(fig, "M1_1.7_ZEB1_ZEB2_LMO2_paired_forest.pdf"), 7.6, 5.2)

# =====================================================================
# 1.8 heatmaps
# =====================================================================
mat_cls <- as.matrix(mile[, lapply(.SD, median, na.rm = TRUE), by = disease, .SDcols = genes][, -1, with = FALSE])
rownames(mat_cls) <- as.character(mile[, lapply(.SD, median, na.rm = TRUE), by = disease, .SDcols = genes]$disease)
heat_col <- colorRampPalette(c("#3C5488", "white", "#E64B35"))(60)

if (has_cht) {
  library(ComplexHeatmap)
  if (has_circlize) {
    col_fun <- circlize::colorRamp2(c(-1.5, 0, 1.5), c("#3C5488", "white", "#E64B35"))
  } else {
    col_fun <- heat_col
  }
  z_cls <- t(scale(mat_cls))
  pdf(file.path(fig, "M1_1.8_ZEB1_ZEB2_LMO2_class_heatmap.pdf"), width = 6.2, height = 3.6)
  ht1 <- Heatmap(z_cls, name = "row z",
                 col = col_fun,
                 cluster_rows = TRUE, cluster_columns = TRUE,
                 column_names_rot = 45,
                 heatmap_legend_param = list(direction = "horizontal"),
                 column_title = "Median expression by disease (row-scaled)")
  draw(ht1, heatmap_legend_side = "bottom")
  dev.off()
  png(file.path(fig, "M1_1.8_ZEB1_ZEB2_LMO2_class_heatmap.png"), width = 6.2 * 180, height = 3.6 * 180, res = 180)
  draw(ht1, heatmap_legend_side = "bottom")
  dev.off()

  set.seed(1)
  ss <- mile[, .SD[sample(.N, min(.N, 40))], by = disease]
  hm <- t(scale(as.matrix(ss[, ..genes])))
  colnames(hm) <- paste0(ss$disease, "_", seq_len(nrow(ss)))
  ha <- HeatmapAnnotation(
    Disease = ss$disease,
    col = list(Disease = pal_dz),
    show_annotation_name = TRUE,
    annotation_name_side = "left"
  )
  pdf(file.path(fig, "M1_1.8_ZEB1_ZEB2_LMO2_sample_heatmap.pdf"), width = 10.5, height = 3.8)
  ht2 <- Heatmap(hm, name = "row z", col = col_fun,
                 top_annotation = ha, show_column_names = FALSE,
                 cluster_columns = TRUE, cluster_rows = FALSE,
                 column_split = ss$disease, column_title = NULL,
                 column_gap = unit(1.2, "mm"))
  draw(ht2)
  dev.off()
  png(file.path(fig, "M1_1.8_ZEB1_ZEB2_LMO2_sample_heatmap.png"), width = 10.5 * 180, height = 3.8 * 180, res = 180)
  draw(ht2)
  dev.off()
} else {
  library(pheatmap)
  pdf(file.path(fig, "M1_1.8_ZEB1_ZEB2_LMO2_class_heatmap.pdf"), width = 5.5, height = 4)
  pheatmap(t(mat_cls), scale = "row", color = heat_col, border_color = NA,
           fontsize = 11, main = "Median expression by disease class")
  dev.off()
  set.seed(1)
  ss <- mile[, .SD[sample(.N, min(.N, 40))], by = disease]
  hm <- t(as.matrix(ss[, ..genes]))
  colnames(hm) <- ss$gsm
  ann <- data.frame(Disease = ss$disease, row.names = ss$gsm)
  pdf(file.path(fig, "M1_1.8_ZEB1_ZEB2_LMO2_sample_heatmap.pdf"), width = 10, height = 4)
  pheatmap(hm, annotation_col = ann, scale = "row", show_colnames = FALSE,
           color = heat_col, border_color = NA, fontsize_row = 11,
           annotation_colors = list(Disease = pal_dz))
  dev.off()
}

# =====================================================================
# 1.9 / 1.10 ROC ggplot + 1.9 multi-ROC (code2 90 / 91)
# =====================================================================
plot_roc_gg <- function(dfs, title, file, width = 5.4, height = 5.2) {
  leg <- dfs[, .(lab = sprintf("%s (%s->T-ALL)  AUC = %.3f (%.3f-%.3f)",
                               name[1], direction[1], auc[1], lo[1], hi[1])), by = name]
  dfs[, name := factor(name, levels = unique(name))]
  pal <- setNames(c("#E64B35", "#3C5488", "#00A087", "#F39B7F")[seq_along(levels(dfs$name))], levels(dfs$name))
  p <- ggplot(dfs, aes(fpr, tpr, color = name)) +
    geom_abline(slope = 1, intercept = 0, linetype = 2, color = "grey55", linewidth = 0.4) +
    geom_line(linewidth = 1.05) +
    scale_color_manual(values = pal, labels = leg$lab[match(levels(dfs$name), leg$name)]) +
    scale_x_continuous(limits = c(0, 1), expand = c(0, 0)) +
    scale_y_continuous(limits = c(0, 1), expand = c(0, 0)) +
    coord_fixed(ratio = 1) +
    labs(title = title, x = "1 - Specificity", y = "Sensitivity", color = NULL) +
    theme_sci() +
    theme(legend.position = c(0.62, 0.18), legend.background = element_blank(),
          legend.text = element_text(size = 8.5))
  save_both(p, file, width, height)
}

r19 <- roc_df(as.integer(vb$Type == "T-ALL"), vb$ZEB1, "ZEB1")
plot_roc_gg(r19, "ZEB1 distinguishes T-ALL from B-ALL", file.path(fig, "M1_1.9_ROC_TALL_vs_BALL.pdf"))

non <- copy(mile)
r10 <- roc_df(as.integer(non$disease == "T-ALL"), non$ZEB1, "ZEB1")
plot_roc_gg(r10, "ZEB1 distinguishes T-ALL from all other MILE classes",
            file.path(fig, "M1_1.10_ROC_TALL_vs_nonTALL.pdf"))

r_multi <- rbind(
  roc_df(as.integer(vb$Type == "T-ALL"), vb$ZEB1, "ZEB1"),
  roc_df(as.integer(vb$Type == "T-ALL"), vb$ZEB2, "ZEB2"),
  roc_df(as.integer(vb$Type == "T-ALL"), vb$LMO2, "LMO2")
)
plot_roc_gg(r_multi, "T-ALL vs B-ALL: three-gene ROC",
            file.path(fig, "M1_1.9_multiROC_ZEB1_ZEB2_LMO2.pdf"), 5.8, 5.3)
fwrite(r_multi[, .(gene = name[1], auc = auc[1], ci_lo = lo[1], ci_hi = hi[1], direction = direction[1]), by = name],
       file.path(tab, "M1_1.9_multiROC_three_genes.tsv"), sep = "\t")

# =====================================================================
# 1.11 / 1.12 independent cohorts + three-cohort panel
# =====================================================================
if (nrow(tgt)) {
  p11 <- plot_rain(tgt, "Type", "ZEB1_log", "Type", pal_tb,
                   "log2(TPM + 1) ZEB1",
                   "TARGET-ALL-P2 primary: ZEB1 in T-ALL vs B-ALL",
                   annot_two(tgt[disease == "T-ALL", ZEB1_log], tgt[disease == "B-ALL", ZEB1_log]))
  p11 <- add_bracket(p11, 1, 2, max(tgt$ZEB1_log, na.rm = TRUE) * 1.04,
                     p_star(wilcox_p(tgt[disease == "T-ALL", ZEB1_log], tgt[disease == "B-ALL", ZEB1_log])))
  save_both(p11, file.path(fig, "M1_1.11_TARGET_ZEB1_violin.pdf"), 5.6, 5.4)
}
if (nrow(ph)) {
  p12p <- plot_rain(ph, "Type", "ZEB1_plot", "Type", pal_tb,
                    "log2(FPKM + 1) ZEB1",
                    "St. Jude Pharmacotype: ZEB1 in T-ALL vs B-ALL",
                    annot_two(ph[disease == "T-ALL", ZEB1_plot], ph[disease == "B-ALL", ZEB1_plot]))
  p12p <- add_bracket(p12p, 1, 2, max(ph$ZEB1_plot, na.rm = TRUE) * 1.04,
                      p_star(wilcox_p(ph[disease == "T-ALL", ZEB1_plot], ph[disease == "B-ALL", ZEB1_plot])))
  save_both(p12p, file.path(fig, "M1_1.12_Pharmacotype_ZEB1_violin.pdf"), 5.6, 5.4)
}

three <- rbindlist(list(
  vb[, .(cohort = "GSE13159\nmicroarray", Type, expr = ZEB1)],
  if (nrow(tgt)) tgt[, .(cohort = "TARGET-ALL-P2\nRNA-seq", Type, expr = ZEB1_log)] else NULL,
  if (nrow(ph)) ph[, .(cohort = "Pharmacotype\nRNA-seq", Type, expr = ZEB1_plot)] else NULL
), fill = TRUE)
three[, cohort := factor(cohort, levels = unique(cohort))]
p3c <- ggplot(three, aes(x = Type, y = expr, fill = Type)) +
  geom_violin(trim = FALSE, scale = "width", color = NA, alpha = 0.72, width = 0.88) +
  add_points(size = 0.28, alpha = 0.28, width = 0.16) +
  geom_boxplot(width = 0.16, fill = "white", outlier.shape = NA, linewidth = 0.35) +
  scale_fill_manual(values = pal_tb) +
  facet_wrap(~ cohort, scales = "free_y", nrow = 1) +
  labs(title = "ZEB1 is higher in T-ALL than B-ALL in three independent cohorts",
       subtitle = "Y axes are platform-native and not pooled; compare within panel only",
       x = NULL, y = "ZEB1 (platform-native scale)") +
  theme_sci() + theme(legend.position = "none")
save_both(p3c, file.path(fig, "M1_1.11_1.12_three_cohort_violin.pdf"), 10.2, 4.8)

# =====================================================================
# 1.13-1.16 ggplot forests (code2 88 layout, ggplot)
# =====================================================================
plot_forest_gg <- function(dt, file, title, subtitle = NULL, width = 8.2, height = NULL,
                           xlab = "Hedges' g (T-ALL minus comparator)") {
  dt <- copy(dt)[is.finite(hedges_g)]
  if (!"label" %in% names(dt)) dt[, label := paste0(cohort, " | ", contrast)]
  dt[, pooled := grepl("Pooled|Leave out", label)]
  dt_np <- dt[pooled == FALSE][order(hedges_g)]
  dt_p <- dt[pooled == TRUE]
  dt <- rbind(dt_np, dt_p)
  dt[, lab := factor(label, levels = label)]
  if (is.null(height)) height <- max(3.8, 0.55 * nrow(dt) + 1.8)
  xr <- c(min(dt$ci95_lo, na.rm = TRUE) - 0.15, max(dt$ci95_hi, na.rm = TRUE) + 0.95)
  p <- ggplot(dt, aes(x = hedges_g, y = lab)) +
    geom_vline(xintercept = 0, linetype = 2, color = "grey40", linewidth = 0.45) +
    geom_errorbar(aes(xmin = ci95_lo, xmax = ci95_hi), width = 0.18, linewidth = 0.6, color = "#2B4C7E") +
    geom_point(aes(shape = pooled, fill = hedges_g > 0), size = 3.4, color = "black", stroke = 0.35) +
    geom_text(aes(x = ci95_hi, label = sprintf("%.2f (%.2f, %.2f)", hedges_g, ci95_lo, ci95_hi)),
              hjust = -0.08, size = 3.1, color = "grey20") +
    scale_shape_manual(values = c("FALSE" = 21, "TRUE" = 23), guide = "none") +
    scale_fill_manual(values = c("TRUE" = "#E64B35", "FALSE" = "#4DBBD5"), guide = "none") +
    labs(title = title, subtitle = subtitle, x = xlab, y = NULL) +
    theme_sci() +
    coord_cartesian(xlim = xr)
  save_both(p, file, width, height)
}

plot_forest_gg(
  meta_in[, .(cohort, contrast, hedges_g, ci95_lo, ci95_hi, label = cohort)],
  file.path(fig, "M1_1.13_cohort_effect_forest.pdf"),
  "Cohort-level effect: T-ALL vs B-ALL",
  xlab = "Hedges' g (T-ALL minus B-ALL)"
)

meta_f <- rbind(
  meta_in[, .(label = cohort, hedges_g, ci95_lo, ci95_hi, contrast)],
  data.table(label = sprintf("Pooled RE  (I2=%.0f%%)", meta_sum$I2[1]),
             hedges_g = meta_sum$hedges_g[1], ci95_lo = meta_sum$ci95_lo[1],
             ci95_hi = meta_sum$ci95_hi[1], contrast = "T-ALL vs B-ALL")
)
plot_forest_gg(meta_f, file.path(fig, "M1_1.14_1.15_meta_forest.pdf"),
               "Random-effects meta-analysis of ZEB1 (T-ALL vs B-ALL)",
               sprintf("DerSimonian-Laird; k = %d; pooled g = %.2f (%.2f-%.2f)",
                       meta_sum$k[1], meta_sum$hedges_g[1], meta_sum$ci95_lo[1], meta_sum$ci95_hi[1]),
               xlab = "Hedges' g (T-ALL minus B-ALL)")

if (nrow(loo)) {
  loo2 <- data.table(label = paste0("Leave out ", loo$left_out),
                     hedges_g = loo$hedges_g, ci95_lo = loo$ci95_lo, ci95_hi = loo$ci95_hi,
                     contrast = "meta")
  plot_forest_gg(loo2, file.path(fig, "M1_1.16_leave_one_cohort_out_forest.pdf"),
                 "Leave-one-cohort-out pooled g")
}

plat <- copy(meta_in)
plat[, label := paste0(ifelse(grepl("GSE13159", cohort), "microarray | ", "RNA-seq | "), cohort)]
plot_forest_gg(plat, file.path(fig, "M1_1.18_platform_sensitivity_forest.pdf"),
               "Platform sensitivity: array vs RNA-seq")

plot_forest_gg(
  stats_pair[, .(label = sub("^T-ALL vs ", "", contrast), hedges_g, ci95_lo, ci95_hi, contrast, cohort)],
  file.path(fig, "M1_1.2_1.5_pairwise_forest.pdf"),
  "Pairwise Hedges' g: T-ALL vs each MILE class"
)

cat("ENRICH DONE\n")
cat("beeswarm=", has_beeswarm, " alluvial=", has_alluvial,
    " ComplexHeatmap=", has_cht, " cowplot=", has_cowplot, "\n")

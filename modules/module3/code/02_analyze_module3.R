# Module 3: ZEB1 molecular context in T-ALL
# Plot recipes from code2 (raincloud violin+box+jitter, lollipop, forest).
# Statistics: two-sided Wilcoxon, Hedges' g + 95% CI, BH-FDR. ZEB1 kept continuous.

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

has_beeswarm <- requireNamespace("ggbeeswarm", quietly = TRUE)
has_cht <- requireNamespace("ComplexHeatmap", quietly = TRUE)
has_circlize <- requireNamespace("circlize", quietly = TRUE)
if (has_beeswarm) library(ggbeeswarm)

root <- "/data-b/liangfuhua/projects/ZEB1_TALL_analysis/module3"
fig <- file.path(root, "figures")
tab <- file.path(root, "tables")
proc <- file.path(root, "processed")
dir.create(fig, FALSE, TRUE)
dir.create(tab, FALSE, TRUE)

pal_sub <- c(
  "TAL" = "#E64B35", "TLX" = "#4DBBD5", "HOXA" = "#00A087",
  "LMO2/LYL1" = "#3C5488", "NKX2_1" = "#F39B7F", "Unknown" = "#8491B4"
)
pal_etp <- c("notETP" = "#4DBBD5", "nearETP" = "#F39B7F", "ETP" = "#E64B35")
pal_mut <- c(wt = "#91D1C2", mut = "#E64B35")
pal_cimp <- c("CIMP-low" = "#4DBBD5", "CIMP-high" = "#E64B35")
pal_fus <- c(
  "no_fusion" = "#91D1C2", TAL = "#E64B35", TLX = "#4DBBD5", HOXA = "#00A087",
  LMO = "#3C5488", KMT2A = "#7E6148", MLLT10 = "#B09C85",
  "NKX2-1" = "#F39B7F", other_fusion = "#8491B4", unknown = "grey70"
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
      strip.text = element_text(face = "bold", color = "black")
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

add_points <- function(size = 0.45, alpha = 0.45, width = 0.18) {
  if (has_beeswarm) geom_quasirandom(size = size, alpha = alpha, width = width, color = "grey20")
  else geom_jitter(size = size, alpha = alpha, width = width, color = "grey20")
}

plot_rain <- function(df, x, y, fill_var, pal, ylab, title, subtitle = NULL, angle = 35) {
  ggplot(df, aes(x = .data[[x]], y = .data[[y]], fill = .data[[fill_var]])) +
    geom_violin(trim = FALSE, scale = "width", color = NA, alpha = 0.72, width = 0.88) +
    add_points(size = 0.5, alpha = 0.45, width = 0.16) +
    geom_boxplot(width = 0.16, fill = "white", outlier.shape = NA, linewidth = 0.4) +
    scale_fill_manual(values = pal, drop = FALSE) +
    labs(x = NULL, y = ylab, title = title, subtitle = subtitle) +
    theme_sci() +
    theme(legend.position = "none",
          axis.text.x = element_text(angle = angle, hjust = if (angle > 0) 1 else 0.5))
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
compare_row <- function(cohort, contrast, gene, a, b, group_a, group_b) {
  hg <- hedges_g(a, b)
  data.table(cohort = cohort, contrast = contrast, gene = gene,
             group_a = group_a, n_a = sum(is.finite(a)), mean_a = mean(a, na.rm = TRUE),
             group_b = group_b, n_b = sum(is.finite(b)), mean_b = mean(b, na.rm = TRUE),
             hedges_g = hg$g, se = hg$se, ci95_lo = hg$lo, ci95_hi = hg$hi,
             wilcox_p = wilcox_p(a, b))
}
annot_two <- function(a, b) {
  hg <- hedges_g(a, b)
  sprintf("n = %d vs %d; %s; g = %.2f [%.2f, %.2f]",
          sum(is.finite(a)), sum(is.finite(b)), fmt_p(wilcox_p(a, b)), hg$g, hg$lo, hg$hi)
}

plot_lollipop <- function(dt, file, title, xlab = "Hedges' g (group A minus group B)", width = 8, height = NULL) {
  dt <- dt[is.finite(hedges_g)]
  if (!nrow(dt)) return(invisible(NULL))
  dt <- unique(copy(dt), by = c("cohort", "contrast", "gene"))
  dt[, lab := make.unique(paste0(cohort, " | ", contrast))]
  setorder(dt, hedges_g)
  dt[, lab := factor(lab, levels = unique(lab))]
  if (is.null(height)) height <- max(4.2, 0.38 * nrow(dt) + 1.8)
  hj <- ifelse(dt$hedges_g >= 0, -0.12, 1.12)
  p <- ggplot(dt, aes(x = hedges_g, y = lab)) +
    geom_vline(xintercept = 0, linetype = 2, color = "grey40", linewidth = 0.4) +
    geom_segment(aes(x = 0, xend = hedges_g, yend = lab), color = "#B07D3A", linewidth = 1.0) +
    geom_errorbar(aes(xmin = ci95_lo, xmax = ci95_hi), width = 0.18, color = "grey25", linewidth = 0.45) +
    geom_point(aes(fill = hedges_g < 0), shape = 21, size = 3.8, color = "black", stroke = 0.35) +
    geom_text(aes(label = sprintf("g=%.2f  %s", hedges_g, p_star(wilcox_p))),
              hjust = hj, size = 2.8) +
    scale_fill_manual(values = c("TRUE" = "#3C5488", "FALSE" = "#E64B35"), guide = "none") +
    labs(title = title, subtitle = "Negative g = lower ZEB1 in group A; whiskers = 95% CI",
         x = xlab, y = NULL) +
    theme_sci()
  save_both(p, file, width, height)
}

tgt <- fread(file.path(proc, "target_tall_subtype.tsv"))
tgt[, ZEB1_log := log2(ZEB1 + 1)]
tgt[, ZEB2_log := log2(ZEB2 + 1)]
tgt[, LMO2_log := log2(LMO2 + 1)]
tgt[, ratio_zeb := log2((ZEB1 + 1e-6) / (ZEB2 + 1e-6))]
tgt[, ratio_lmo := log2((LMO2 + 1e-6) / (ZEB1 + 1e-6))]
tgt[is.na(subtype) | subtype == "", subtype := "Unknown"]
tgt[, subtype := factor(subtype, levels = intersect(names(pal_sub), unique(as.character(subtype))))]
if ("etp" %in% names(tgt)) tgt[etp %in% c("nan", "NA", "Unevaluable", ""), etp := NA]

ph <- fread(file.path(proc, "pharmacotype_tall_subtype.tsv"))
if (nrow(ph) && max(ph$ZEB1, na.rm = TRUE) > 50) {
  ph[, ZEB1_log := log2(ZEB1 + 1)]
  ph[, ZEB2_log := log2(ZEB2 + 1)]
  ph[, LMO2_log := log2(LMO2 + 1)]
} else if (nrow(ph)) {
  ph[, ZEB1_log := ZEB1]; ph[, ZEB2_log := ZEB2]; ph[, LMO2_log := LMO2]
}

ci <- fread(file.path(proc, "gse272023_cimp.tsv"))
if (!"cimp_g" %in% names(ci)) {
  ci[, cimp_g := fifelse(grepl("high", cimp, ignore.case = TRUE), "CIMP-high",
                         fifelse(grepl("low", cimp, ignore.case = TRUE), "CIMP-low", NA_character_))]
}
if ("is_tumor" %in% names(ci)) ci <- ci[is_tumor == TRUE | is.na(is_tumor)]
ci2 <- ci[cimp_g %in% c("CIMP-high", "CIMP-low") & is.finite(ZEB1)]

stats <- list()

cnt <- tgt[, .N, by = .(subtype, etp)]
fwrite(cnt, file.path(tab, "M3_3.1_subtype_etp_counts.tsv"), sep = "\t")
fwrite(tgt[, .N, by = subtype], file.path(tab, "M3_3.1_subtype_counts.tsv"), sep = "\t")
sc <- tgt[, .N, by = subtype]
sc[, pct := 100 * N / sum(N)]
p1 <- ggplot(sc, aes(x = reorder(subtype, -N), y = N, fill = subtype)) +
  geom_col(width = 0.72, color = "white", linewidth = 0.3) +
  geom_text(aes(label = sprintf("%d\n%.0f%%", N, pct)), vjust = -0.12, size = 3.1, lineheight = 0.9) +
  scale_fill_manual(values = pal_sub, drop = FALSE) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.18))) +
  labs(title = "TARGET T-ALL molecular subtype composition",
       subtitle = sprintf("n = %d primary unique patients", nrow(tgt)),
       x = NULL, y = "n patients") +
  theme_sci() + theme(legend.position = "none")
save_both(p1, file.path(fig, "M3_3.1_subtype_counts.pdf"), 7.2, 5.0)

use <- tgt[!is.na(subtype) & subtype != "Unknown"]
long <- melt(use, id.vars = "subtype",
             measure.vars = c("ZEB1_log", "ZEB2_log", "LMO2_log"),
             variable.name = "gene", value.name = "expr")
setDT(long)
long[, gene := gsub("_log", "", gene)]
p234 <- ggplot(long, aes(x = subtype, y = expr, fill = subtype)) +
  geom_violin(trim = FALSE, scale = "width", color = NA, alpha = 0.72, width = 0.88) +
  add_points(size = 0.4, alpha = 0.4, width = 0.15) +
  geom_boxplot(width = 0.16, fill = "white", outlier.shape = NA, linewidth = 0.35) +
  facet_wrap(~gene, scales = "free_y") +
  scale_fill_manual(values = pal_sub) +
  labs(title = "TARGET: ZEB1 / ZEB2 / LMO2 by T-ALL subtype", y = "log2(TPM+1)") +
  theme_sci() + theme(legend.position = "none", axis.text.x = element_text(angle = 35, hjust = 1))
save_both(p234, file.path(fig, "M3_3.2_3.4_genes_by_subtype.pdf"), 10.5, 5.6)

p32 <- plot_rain(use, "subtype", "ZEB1_log", "subtype", pal_sub,
                 "log2(TPM+1) ZEB1", "TARGET: ZEB1 by molecular subtype")
save_both(p32, file.path(fig, "M3_3.2_ZEB1_by_subtype.pdf"), 7.4, 5.4)

sub_stats <- rbindlist(lapply(c("ZEB1_log", "ZEB2_log", "LMO2_log"), function(gn) {
  rbindlist(lapply(as.character(unique(use$subtype)), function(s) {
    a <- use[as.character(subtype) == s][[gn]]
    b <- use[as.character(subtype) != s][[gn]]
    compare_row("TARGET-ALL-P2", paste0(s, " vs rest"), gsub("_log", "", gn), a, b, s, "rest")
  }))
}))
sub_stats[, fdr := p.adjust(wilcox_p, "BH")]
fwrite(sub_stats, file.path(tab, "M3_3.2_3.4_subtype_effects.tsv"), sep = "\t")
plot_lollipop(sub_stats[gene == "ZEB1"], file.path(fig, "M3_3.2_ZEB1_subtype_lollipop.pdf"),
              "TARGET: ZEB1 subtype vs rest")
stats[["subtype"]] <- sub_stats

p5 <- plot_rain(use, "subtype", "ratio_zeb", "subtype", pal_sub,
                "log2((ZEB1+eps)/(ZEB2+eps))", "TARGET: ZEB1/ZEB2 ratio by subtype")
save_both(p5, file.path(fig, "M3_3.5_ZEB1_ZEB2_ratio.pdf"), 7.4, 5.4)
p6 <- plot_rain(use, "subtype", "ratio_lmo", "subtype", pal_sub,
                "log2((LMO2+eps)/(ZEB1+eps))", "TARGET: LMO2/ZEB1 ratio by subtype")
save_both(p6, file.path(fig, "M3_3.6_LMO2_ZEB1_ratio.pdf"), 7.4, 5.4)

et <- tgt[etp %in% c("ETP", "nearETP", "notETP")]
et[, etp := factor(etp, levels = c("notETP", "nearETP", "ETP"))]
p7 <- plot_rain(et, "etp", "ZEB1_log", "etp", pal_etp, "log2(TPM+1) ZEB1",
                "TARGET: ZEB1 by ETP status",
                annot_two(et[etp == "ETP", ZEB1_log], et[etp == "notETP", ZEB1_log]), angle = 0)
save_both(p7, file.path(fig, "M3_3.7_ZEB1_ETP.pdf"), 5.8, 5.4)
etp_stats <- rbind(
  compare_row("TARGET-ALL-P2", "ETP vs notETP", "ZEB1", et[etp == "ETP", ZEB1_log], et[etp == "notETP", ZEB1_log], "ETP", "notETP"),
  compare_row("TARGET-ALL-P2", "nearETP vs notETP", "ZEB1", et[etp == "nearETP", ZEB1_log], et[etp == "notETP", ZEB1_log], "nearETP", "notETP"),
  compare_row("TARGET-ALL-P2", "ETP vs notETP", "ZEB2", et[etp == "ETP", ZEB2_log], et[etp == "notETP", ZEB2_log], "ETP", "notETP"),
  compare_row("TARGET-ALL-P2", "ETP vs notETP", "LMO2", et[etp == "ETP", LMO2_log], et[etp == "notETP", LMO2_log], "ETP", "notETP")
)
stats[["etp"]] <- etp_stats

bin_stats <- rbind(
  compare_row("TARGET-ALL-P2", "LMO2/LYL1 vs rest", "ZEB1",
              tgt[subtype == "LMO2/LYL1", ZEB1_log], tgt[subtype != "LMO2/LYL1" & subtype != "Unknown", ZEB1_log], "LMO2/LYL1", "rest"),
  compare_row("TARGET-ALL-P2", "HOXA vs rest", "ZEB1",
              tgt[subtype == "HOXA", ZEB1_log], tgt[subtype != "HOXA" & subtype != "Unknown", ZEB1_log], "HOXA", "rest"),
  compare_row("TARGET-ALL-P2", "TAL vs rest", "ZEB1",
              tgt[subtype == "TAL", ZEB1_log], tgt[subtype != "TAL" & subtype != "Unknown", ZEB1_log], "TAL", "rest"),
  compare_row("TARGET-ALL-P2", "TLX vs rest", "ZEB1",
              tgt[subtype == "TLX", ZEB1_log], tgt[subtype != "TLX" & subtype != "Unknown", ZEB1_log], "TLX", "rest")
)
stats[["bin"]] <- bin_stats
for (s in c("LMO2/LYL1", "HOXA", "TAL")) {
  d <- copy(tgt[subtype != "Unknown"])
  d[, grp := factor(ifelse(as.character(subtype) == s, s, "rest"), levels = c("rest", s))]
  pal2 <- c("rest" = "#B8B8B8"); pal2[[s]] <- unname(pal_sub[[s]])
  pp <- plot_rain(d, "grp", "ZEB1_log", "grp", pal2, "log2(TPM+1) ZEB1",
                  sprintf("TARGET: ZEB1 in %s vs rest", s),
                  annot_two(d[grp == s, ZEB1_log], d[grp == "rest", ZEB1_log]), angle = 0)
  tag <- gsub("[^A-Za-z0-9]+", "_", s)
  save_both(pp, file.path(fig, sprintf("M3_%s_vs_rest.pdf", tag)), 5.6, 5.3)
}

mut_stats <- data.table()
for (g in c("NOTCH1_mut", "PTEN_mut", "RAS_JAK", "FBXW7_mut", "PHF6_mut")) {
  if (!g %in% names(tgt)) next
  lab <- sub("_mut$", "", g)
  if (sum(tgt[[g]] == 1, na.rm = TRUE) < 5) next
  mut_stats <- rbind(mut_stats, compare_row("TARGET-ALL-P2", paste0(lab, " mut vs wt"), "ZEB1",
                                            tgt[get(g) == 1, ZEB1_log], tgt[get(g) == 0, ZEB1_log], "mut", "wt"))
  pd <- copy(tgt[!is.na(get(g))])
  pd$Status <- factor(ifelse(pd[[g]] == 1, "mut", "wt"), levels = c("wt", "mut"))
  pp <- plot_rain(pd, "Status", "ZEB1_log", "Status", pal_mut, "log2(TPM+1) ZEB1",
                  sprintf("TARGET: ZEB1 vs %s mutation", lab),
                  annot_two(pd[Status == "mut", ZEB1_log], pd[Status == "wt", ZEB1_log]), angle = 0)
  save_both(pp, file.path(fig, sprintf("M3_mut_%s.pdf", lab)), 5.4, 5.3)
}
stats[["mut"]] <- mut_stats
fwrite(mut_stats, file.path(tab, "M3_3.11_3.13_mutation_effects.tsv"), sep = "\t")

if ("cnv_burden" %in% names(tgt) && sum(is.finite(tgt$cnv_burden)) >= 20) {
  ct <- cor.test(tgt$ZEB1_log, tgt$cnv_burden, method = "spearman", exact = FALSE)
  p14 <- ggplot(tgt[is.finite(cnv_burden)], aes(cnv_burden, ZEB1_log, color = subtype)) +
    geom_point(alpha = 0.7, size = 1.6) +
    geom_smooth(method = "lm", formula = y ~ x, color = "grey20", se = TRUE, linewidth = 0.6,
                inherit.aes = FALSE, aes(x = cnv_burden, y = ZEB1_log)) +
    scale_color_manual(values = pal_sub, drop = FALSE) +
    labs(title = sprintf("TARGET: ZEB1 vs CNV burden  rho = %.2f  %s", unname(ct$estimate), fmt_p(ct$p.value)),
         x = "CNV burden (fraction of genes with |CN-2| > 0.5)", y = "log2(TPM+1) ZEB1", color = "Subtype") +
    theme_sci()
  save_both(p14, file.path(fig, "M3_3.14_ZEB1_CNV_burden.pdf"), 7.2, 5.4)
  fwrite(data.table(rho = unname(ct$estimate), p = ct$p.value, n = sum(is.finite(tgt$cnv_burden))),
         file.path(tab, "M3_3.14_cnv_burden.tsv"), sep = "\t")
}

if ("fusion_class" %in% names(tgt)) {
  fus <- tgt[fusion_class %in% names(pal_fus) & fusion_class != "unknown"]
  fus[, fusion_class := factor(fusion_class, levels = intersect(names(pal_fus), unique(fusion_class)))]
  p15 <- plot_rain(fus, "fusion_class", "ZEB1_log", "fusion_class", pal_fus,
                   "log2(TPM+1) ZEB1", "TARGET: ZEB1 by fusion class")
  save_both(p15, file.path(fig, "M3_3.15_ZEB1_fusion_class.pdf"), 8.0, 5.4)
  fwrite(tgt[, .N, by = fusion_class], file.path(tab, "M3_3.15_fusion_counts.tsv"), sep = "\t")
  fus_stats <- rbindlist(lapply(setdiff(unique(as.character(fus$fusion_class)), "no_fusion"), function(s) {
    compare_row("TARGET-ALL-P2", paste0(s, " vs no_fusion"), "ZEB1",
                fus[as.character(fusion_class) == s, ZEB1_log],
                fus[fusion_class == "no_fusion", ZEB1_log], s, "no_fusion")
  }))
  stats[["fus"]] <- fus_stats
}

if ("ZEB1_cn" %in% names(tgt) && sum(is.finite(tgt$ZEB1_cn)) >= 20) {
  ct <- cor.test(tgt$ZEB1_log, tgt$ZEB1_cn, method = "spearman", exact = FALSE)
  p16 <- ggplot(tgt[is.finite(ZEB1_cn)], aes(ZEB1_cn, ZEB1_log, color = subtype)) +
    geom_point(alpha = 0.7, size = 1.6) +
    geom_smooth(method = "lm", formula = y ~ x, color = "grey20", se = TRUE, linewidth = 0.6,
                inherit.aes = FALSE, aes(x = ZEB1_cn, y = ZEB1_log)) +
    scale_color_manual(values = pal_sub, drop = FALSE) +
    labs(title = sprintf("TARGET: ZEB1 expression vs locus copy number  rho = %.2f  %s", unname(ct$estimate), fmt_p(ct$p.value)),
         x = "ZEB1 copy number", y = "log2(TPM+1) ZEB1", color = "Subtype") +
    theme_sci()
  save_both(p16, file.path(fig, "M3_3.16_ZEB1_locus_CNV.pdf"), 7.2, 5.4)
  fwrite(data.table(rho = unname(ct$estimate), p = ct$p.value, n = sum(is.finite(tgt$ZEB1_cn))),
         file.path(tab, "M3_3.16_ZEB1_locus_CNV.tsv"), sep = "\t")
}

prev <- data.table(
  gene = c("ZEB1", "ZEB2", "LMO2"),
  n_mut = c(sum(tgt$ZEB1_mut == 1, na.rm = TRUE),
            if ("ZEB2_mut" %in% names(tgt)) sum(tgt$ZEB2_mut == 1, na.rm = TRUE) else NA_real_,
            if ("LMO2_mut" %in% names(tgt)) sum(tgt$LMO2_mut == 1, na.rm = TRUE) else NA_real_),
  n_fusion = c(if ("zeb1_fusion" %in% names(tgt)) sum(tgt$zeb1_fusion, na.rm = TRUE) else 0, NA, NA),
  n = nrow(tgt)
)
fwrite(prev, file.path(tab, "M3_3.17_3.18_ZEB1_mutation_fusion.tsv"), sep = "\t")
writeLines(c(
  "ZEB1 / ZEB2 / LMO2 coding mutations and ZEB1 fusions in TARGET T-ALL.",
  sprintf("ZEB1 mut: %d / %d", prev$n_mut[1], prev$n[1]),
  sprintf("ZEB2 mut: %d / %d", prev$n_mut[2], prev$n[1]),
  sprintf("LMO2 mut: %d / %d", prev$n_mut[3], prev$n[1]),
  sprintf("ZEB1 fusion (Lesion string): %d / %d", prev$n_fusion[1], prev$n[1]),
  "Cis-coding alteration of ZEB1 is essentially absent; expression differences are not explained by ZEB1 mutation or fusion."
), file.path(tab, "M3_3.17_3.18_note.txt"))

fit <- copy(tgt[is.finite(ZEB1_log)])
fit[, etp_bin := fifelse(etp == "ETP", "ETP", fifelse(etp == "notETP", "notETP", NA_character_))]
coef_from <- function(model_name, fit_obj) {
  s <- summary(fit_obj)$coefficients
  data.table(model = model_name, term = rownames(s), estimate = s[, 1], se = s[, 2],
             t = s[, 3], p = s[, 4],
             ci95_lo = s[, 1] - 1.96 * s[, 2], ci95_hi = s[, 1] + 1.96 * s[, 2])
}
mod_rows <- list()
if (sum(!is.na(fit$etp_bin)) >= 20) mod_rows[[length(mod_rows) + 1]] <- coef_from("ETP", lm(ZEB1_log ~ etp_bin, data = fit))
mod_rows[[length(mod_rows) + 1]] <- coef_from("subtype", lm(ZEB1_log ~ subtype, data = fit[subtype != "Unknown"]))
if ("age_years" %in% names(fit) && sum(is.finite(fit$age_years)) >= 20) {
  mod_rows[[length(mod_rows) + 1]] <- coef_from("subtype+age", lm(ZEB1_log ~ subtype + age_years, data = fit[subtype != "Unknown" & is.finite(age_years)]))
  mod_rows[[length(mod_rows) + 1]] <- coef_from("ETP+age", lm(ZEB1_log ~ etp_bin + age_years, data = fit[is.finite(age_years)]))
}
if ("NOTCH1_mut" %in% names(fit)) {
  mod_rows[[length(mod_rows) + 1]] <- coef_from("subtype+NOTCH1", lm(ZEB1_log ~ subtype + NOTCH1_mut, data = fit[subtype != "Unknown"]))
}
mod_dt <- rbindlist(mod_rows)
fwrite(mod_dt, file.path(tab, "M3_3.21_3.23_adjusted_models.tsv"), sep = "\t")
md <- mod_dt[term != "(Intercept)" & is.finite(estimate)]
if (nrow(md)) {
  md[, lab := paste0(model, " | ", term)]
  setorder(md, estimate)
  md[, lab := factor(lab, levels = lab)]
  pmod <- ggplot(md, aes(x = estimate, y = lab)) +
    geom_vline(xintercept = 0, linetype = 2, color = "grey40") +
    geom_errorbar(aes(xmin = ci95_lo, xmax = ci95_hi), width = 0.18, linewidth = 0.45) +
    geom_point(aes(fill = estimate < 0), shape = 21, size = 3.4, color = "black") +
    scale_fill_manual(values = c("TRUE" = "#3C5488", "FALSE" = "#E64B35"), guide = "none") +
    labs(title = "TARGET linear models for log2(TPM+1) ZEB1",
         subtitle = "Coefficients vs reference level; whiskers = 95% CI",
         x = "Coefficient", y = NULL) +
    theme_sci()
  save_both(pmod, file.path(fig, "M3_3.21_3.23_adjusted_coefficients.pdf"), 9.2, max(5, 0.32 * nrow(md) + 1.6))
}

ph_stats <- data.table()
if (nrow(ph) && uniqueN(ph$etp) > 1) {
  ppe <- plot_rain(ph, "etp", "ZEB1_log", "etp", pal_etp[c("notETP", "ETP")],
                   "log2(FPKM+1) ZEB1", "Pharmacotype: ZEB1 in ETP vs conventional T-ALL",
                   annot_two(ph[etp == "ETP", ZEB1_log], ph[etp == "notETP", ZEB1_log]), angle = 0)
  save_both(ppe, file.path(fig, "M3_pharmacotype_ETP.pdf"), 5.6, 5.4)
  ph_stats <- rbind(
    compare_row("StJude_Pharmacotype", "ETP vs notETP", "ZEB1",
                ph[etp == "ETP", ZEB1_log], ph[etp == "notETP", ZEB1_log], "ETP", "notETP"),
    compare_row("StJude_Pharmacotype", "ETP vs notETP", "LMO2",
                ph[etp == "ETP", LMO2_log], ph[etp == "notETP", LMO2_log], "ETP", "notETP")
  )
  if ("age_years" %in% names(ph) && sum(is.finite(ph$age_years)) >= 20) {
    fwrite(coef_from("Pharmacotype ETP+age", lm(ZEB1_log ~ etp + age_years, data = ph)),
           file.path(tab, "M3_pharmacotype_ETP_age_model.tsv"), sep = "\t")
  }
}
if (nrow(ph) && "subtype" %in% names(ph)) {
  pph <- plot_rain(ph, "subtype", "ZEB1_log", "subtype",
                   c("ETP" = "#E64B35", "T-ALL (non-ETP)" = "#4DBBD5", "Unknown" = "grey70"),
                   "log2(FPKM+1) ZEB1", "Pharmacotype T-ALL: published molecular labels", angle = 20)
  save_both(pph, file.path(fig, "M3_pharmacotype_subtype.pdf"), 6.4, 5.4)
}
stats[["ph"]] <- ph_stats

cimp_stats <- data.table()
if (nrow(ci2) >= 10) {
  ci2[, cimp_g := factor(cimp_g, levels = c("CIMP-low", "CIMP-high"))]
  pc <- plot_rain(ci2, "cimp_g", "ZEB1", "cimp_g", pal_cimp, "ZEB1 VST",
                  "GSE272023: ZEB1 by CIMP class",
                  annot_two(ci2[cimp_g == "CIMP-high", ZEB1], ci2[cimp_g == "CIMP-low", ZEB1]), angle = 0)
  save_both(pc, file.path(fig, "M3_3.19_CIMP_ZEB1.pdf"), 5.6, 5.4)
  cimp_stats <- rbind(
    compare_row("GSE272023", "CIMP-high vs CIMP-low", "ZEB1",
                ci2[cimp_g == "CIMP-high", ZEB1], ci2[cimp_g == "CIMP-low", ZEB1], "CIMP-high", "CIMP-low"),
    if ("ZEB2" %in% names(ci2)) compare_row("GSE272023", "CIMP-high vs CIMP-low", "ZEB2",
                ci2[cimp_g == "CIMP-high", ZEB2], ci2[cimp_g == "CIMP-low", ZEB2], "CIMP-high", "CIMP-low") else NULL,
    if ("LMO2" %in% names(ci2)) compare_row("GSE272023", "CIMP-high vs CIMP-low", "LMO2",
                ci2[cimp_g == "CIMP-high", LMO2], ci2[cimp_g == "CIMP-low", LMO2], "CIMP-high", "CIMP-low") else NULL
  )
  fwrite(ci2[, .N, by = cimp_g], file.path(tab, "M3_3.19_cimp_counts.tsv"), sep = "\t")
}
stats[["cimp"]] <- cimp_stats
writeLines(c(
  "3.20 CIMP x MRD: GSE272023 GEO metadata contains CIMP class but not MRD.",
  "PMID 39841000 reports MRD in clinical tables not deposited with the series matrix.",
  "Cannot run CIMP x MRD interaction without that table. Left incomplete, not imputed."
), file.path(tab, "M3_3.20_CIMP_MRD_note.txt"))

if (has_cht && nrow(use) >= 10) {
  library(ComplexHeatmap)
  med_dt <- use[, lapply(.SD, median, na.rm = TRUE), by = subtype, .SDcols = c("ZEB1_log", "ZEB2_log", "LMO2_log")]
  med <- as.matrix(med_dt[, -1])
  rownames(med) <- as.character(med_dt$subtype)
  colnames(med) <- c("ZEB1", "ZEB2", "LMO2")
  z <- scale(med)
  pdf(file.path(fig, "M3_subtype_gene_heatmap.pdf"), width = 5.6, height = 4.2)
  draw(Heatmap(t(z), name = "col z", cluster_rows = TRUE, cluster_columns = TRUE,
               column_names_rot = 45, column_title = "Median log2(TPM+1), column-scaled"))
  dev.off()
  png(file.path(fig, "M3_subtype_gene_heatmap.png"), width = 5.6 * 180, height = 4.2 * 180, res = 180)
  draw(Heatmap(t(z), name = "col z", cluster_rows = TRUE, cluster_columns = TRUE, column_names_rot = 45))
  dev.off()
}

all_long <- rbindlist(stats, fill = TRUE)
all_long[, fdr := p.adjust(wilcox_p, "BH")]
fwrite(all_long, file.path(tab, "M3_all_contrasts.tsv"), sep = "\t")
zeb <- all_long[gene == "ZEB1" & is.finite(hedges_g)]
fwrite(zeb, file.path(tab, "M3_3.24_ZEB1_contrasts.tsv"), sep = "\t")
plot_lollipop(zeb, file.path(fig, "M3_3.24_interaction_forest.pdf"),
              "ZEB1 contrasts across subtype, ETP, mutation, CIMP",
              width = 10, height = max(6.5, 0.34 * nrow(zeb) + 1.8))

writeLines(capture.output(sessionInfo()), file.path(tab, "sessionInfo.txt"))
cat("MODULE3 ANALYZE DONE\n")
cat("n TARGET", nrow(tgt), "pharm", nrow(ph), "cimp", nrow(ci2), "\n")

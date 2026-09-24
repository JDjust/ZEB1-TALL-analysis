# Module 7: ZEB1 vs St. Jude Pharmacotype drug LC50
# code2 raincloud / scatter / lollipop. ZEB1 kept continuous. Two-sided Spearman, BH-FDR.
# MTX is not in this 18-drug panel (see processed/M7_MTX_note.txt).

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})
has_beeswarm <- requireNamespace("ggbeeswarm", quietly = TRUE)
has_cht <- requireNamespace("ComplexHeatmap", quietly = TRUE)
has_circlize <- requireNamespace("circlize", quietly = TRUE)
has_mass <- requireNamespace("MASS", quietly = TRUE)
if (has_beeswarm) library(ggbeeswarm)
colorRamp2_safe <- function(breaks, colors) {
  if (has_circlize) circlize::colorRamp2(breaks, colors) else colors
}

root <- "/data-b/liangfuhua/projects/ZEB1_TALL_analysis/module7"
fig <- file.path(root, "figures")
tab <- file.path(root, "tables")
proc <- file.path(root, "processed")
dir.create(fig, FALSE, TRUE)
dir.create(tab, FALSE, TRUE)
if (file.exists(file.path(proc, "M7_MTX_note.txt"))) {
  file.copy(file.path(proc, "M7_MTX_note.txt"), file.path(tab, "M7_MTX_note.txt"), overwrite = TRUE)
}

pal_etp <- c("notETP" = "#4DBBD5", "ETP" = "#E64B35")
pal_tb <- c("B-ALL" = "#4DBBD5", "T-ALL" = "#E64B35")

theme_sci <- function(base_size = 11) {
  theme_classic(base_size = base_size) +
    theme(
      axis.text = element_text(color = "black"),
      axis.title = element_text(color = "black"),
      axis.line = element_line(linewidth = 0.45, color = "black"),
      axis.ticks = element_line(linewidth = 0.45, color = "black"),
      plot.title = element_text(face = "bold", size = base_size + 1, color = "black"),
      plot.subtitle = element_text(size = base_size - 1, color = "grey25")
    )
}
save_both <- function(plot, file, width, height) {
  pdf(file, width = width, height = height, useDingbats = FALSE); print(plot); dev.off()
  png(sub("\\.pdf$", ".png", file), width = width * 180, height = height * 180, res = 180)
  print(plot); dev.off()
}
add_points <- function(size = 0.5, alpha = 0.5, width = 0.16) {
  if (has_beeswarm) geom_quasirandom(size = size, alpha = alpha, width = width, color = "grey20")
  else geom_jitter(size = size, alpha = alpha, width = width, color = "grey20")
}
fmt_p <- function(p) ifelse(is.na(p), "NA", ifelse(p < 0.001, "P < 0.001", sprintf("P = %.3f", p)))

short_drug <- function(x) {
  x <- gsub(" \\(.*\\)$", "", x)
  x <- gsub("_normalized$", "", x)
  x
}

ph <- fread(file.path(proc, "pharmacotype_tall_drugs.tsv"))
alld <- fread(file.path(proc, "pharmacotype_all_drugs.tsv"))
raw_drugs <- grep("\\((IU/ml|nM|µM|uM|μM)\\)$", names(ph), value = TRUE)
norm_drugs <- grep("_normalized$", names(ph), value = TRUE)
stopifnot(length(raw_drugs) > 0)

# 7.1 sample filter
cnt <- alld[, .N, by = disease]
fwrite(cnt, file.path(tab, "M7_7.1_sample_counts.tsv"), sep = "\t")
p1 <- ggplot(cnt[disease %in% c("T-ALL", "B-ALL")], aes(x = disease, y = N, fill = disease)) +
  geom_col(width = 0.65, color = "white") +
  geom_text(aes(label = N), vjust = -0.3, size = 4) +
  scale_fill_manual(values = pal_tb) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
  labs(title = "Pharmacotype RNA+drug cases used in Module 7", x = NULL, y = "n patients") +
  theme_sci() + theme(legend.position = "none")
save_both(p1, file.path(fig, "M7_7.1_sample_counts.pdf"), 5.2, 4.8)

# 7.2 ZEB1 distribution in T-ALL
p2 <- ggplot(ph, aes(x = ZEB1_log)) +
  geom_histogram(bins = 20, fill = "#E64B35", color = "white", linewidth = 0.3) +
  geom_vline(xintercept = median(ph$ZEB1_log, na.rm = TRUE), linetype = 2) +
  labs(title = "Pharmacotype T-ALL: ZEB1 distribution",
       subtitle = sprintf("n = %d; median log2(FPKM+1) = %.2f", nrow(ph), median(ph$ZEB1_log, na.rm = TRUE)),
       x = "log2(FPKM+1) ZEB1", y = "n") +
  theme_sci()
save_both(p2, file.path(fig, "M7_7.2_ZEB1_distribution.pdf"), 5.6, 4.6)

spearman_row <- function(gene, drug, x, y) {
  ok <- is.finite(x) & is.finite(y)
  x <- x[ok]; y <- y[ok]
  if (length(x) < 8) return(NULL)
  ct <- cor.test(x, y, method = "spearman", exact = FALSE)
  data.table(gene = gene, drug = drug, n = length(x), rho = unname(ct$estimate), p = ct$p.value)
}

# 7.7-7.8 18-drug screen (primary = raw LC50; higher LC50 = more resistant)
scr <- rbindlist(lapply(c("ZEB1_log", "ZEB2_log", "LMO2_log", "ratio_zeb"), function(gn) {
  rbindlist(lapply(raw_drugs, function(d) spearman_row(gn, d, ph[[gn]], ph[[d]])))
}))
scr[, fdr := p.adjust(p, "BH"), by = gene]
scr[, drug_short := short_drug(drug)]
setorder(scr, gene, rho)
fwrite(scr, file.path(tab, "M7_7.7_7.8_drug_spearman.tsv"), sep = "\t")

zeb <- scr[gene == "ZEB1_log"]
setorder(zeb, rho)
zeb[, drug_short := factor(drug_short, levels = drug_short)]
p_lp <- ggplot(zeb, aes(x = rho, y = drug_short)) +
  geom_vline(xintercept = 0, linetype = 2, color = "grey40") +
  geom_segment(aes(x = 0, xend = rho, yend = drug_short), color = "#B07D3A", linewidth = 1) +
  geom_point(aes(fill = fdr < 0.05), shape = 21, size = 3.6, color = "black") +
  geom_text(aes(label = sprintf("rho=%.2f %s", rho, ifelse(fdr < 0.05, "*", ""))),
            hjust = ifelse(zeb$rho >= 0, -0.12, 1.12), size = 2.8) +
  scale_fill_manual(values = c("TRUE" = "#E64B35", "FALSE" = "#4DBBD5"), guide = "none") +
  labs(title = "T-ALL: Spearman ZEB1 vs LC50 (18 drugs)",
       subtitle = "Negative rho = higher ZEB1, lower LC50 (more sensitive). * FDR < 0.05",
       x = "Spearman rho", y = NULL) +
  coord_cartesian(xlim = c(min(zeb$rho, na.rm = TRUE) - 0.12, max(zeb$rho, na.rm = TRUE) + 0.28)) +
  theme_sci()
save_both(p_lp, file.path(fig, "M7_7.7_7.10_ZEB1_drug_lollipop.pdf"), 8.2, 7.2)

# ranking table
rank <- copy(zeb)
setorder(rank, p)
fwrite(rank, file.path(tab, "M7_7.10_sensitivity_ranking.tsv"), sep = "\t")

# scatter for top 4 drugs by |rho| plus mercaptopurine/thioguanine/nelarabine if present
want <- unique(c(
  as.character(zeb[order(-abs(rho))][1:min(4, .N), drug]),
  grep("Mercaptopurine|Thioguanine|Nelarabine|Venetoclax|Dasatinib|Ruxolitinib", raw_drugs, value = TRUE)
))
for (d in want) {
  dd <- ph[is.finite(get(d)) & is.finite(ZEB1_log)]
  if (nrow(dd) < 8) next
  ct <- cor.test(dd$ZEB1_log, dd[[d]], method = "spearman", exact = FALSE)
  pp <- ggplot(dd, aes(x = ZEB1_log, y = .data[[d]], color = etp)) +
    geom_point(size = 2, alpha = 0.8) +
    geom_smooth(method = "lm", formula = y ~ x, color = "grey25", se = TRUE, linewidth = 0.55, inherit.aes = FALSE,
                aes(x = ZEB1_log, y = .data[[d]])) +
    scale_color_manual(values = pal_etp, na.value = "grey50") +
    labs(title = sprintf("ZEB1 vs %s", short_drug(d)),
         subtitle = sprintf("n = %d; Spearman rho = %.2f; %s", nrow(dd), unname(ct$estimate), fmt_p(ct$p.value)),
         x = "log2(FPKM+1) ZEB1", y = d, color = "ETP") +
    theme_sci()
  fn <- sprintf("M7_scatter_%s.pdf", gsub("[^A-Za-z0-9]+", "_", short_drug(d)))
  save_both(pp, file.path(fig, fn), 6.2, 5.2)
}

# 7.6 Q1 vs Q4 for each drug (sensitivity display only)
q <- quantile(ph$ZEB1_log, c(0.25, 0.75), na.rm = TRUE)
ph[, zeb_q := fifelse(ZEB1_log <= q[1], "Q1", fifelse(ZEB1_log >= q[2], "Q4", NA_character_))]
qrows <- rbindlist(lapply(raw_drugs, function(d) {
  a <- ph[zeb_q == "Q4"][[d]]; b <- ph[zeb_q == "Q1"][[d]]
  a <- a[is.finite(a)]; b <- b[is.finite(b)]
  if (length(a) < 4 || length(b) < 4) return(NULL)
  data.table(drug = d, n_q4 = length(a), n_q1 = length(b),
             median_q4 = median(a), median_q1 = median(b),
             wilcox_p = wilcox.test(a, b, exact = FALSE)$p.value)
}))
if (nrow(qrows)) {
  qrows[, fdr := p.adjust(wilcox_p, "BH")]
  fwrite(qrows, file.path(tab, "M7_7.6_Q1Q4.tsv"), sep = "\t")
}

# 7.5 robust / OLS for each drug ~ ZEB1 + age + etp
coef_rows <- list()
for (d in raw_drugs) {
  dat <- ph[is.finite(get(d)) & is.finite(ZEB1_log)]
  if (nrow(dat) < 12) next
  f1 <- as.formula(paste0("`", d, "` ~ ZEB1_log"))
  m1 <- lm(f1, data = dat)
  s1 <- summary(m1)$coefficients
  coef_rows[[length(coef_rows) + 1]] <- data.table(drug = d, model = "ZEB1",
    estimate = s1["ZEB1_log", 1], se = s1["ZEB1_log", 2], p = s1["ZEB1_log", 4], n = nrow(dat))
  if ("age_years" %in% names(dat) || "Age at diagnosis (years)" %in% names(dat)) {
    if (!"age_years" %in% names(dat)) dat[, age_years := `Age at diagnosis (years)`]
    f2 <- as.formula(paste0("`", d, "` ~ ZEB1_log + age_years"))
    m2 <- lm(f2, data = dat[is.finite(age_years)])
    s2 <- summary(m2)$coefficients
    if ("ZEB1_log" %in% rownames(s2)) {
      coef_rows[[length(coef_rows) + 1]] <- data.table(drug = d, model = "ZEB1+age",
        estimate = s2["ZEB1_log", 1], se = s2["ZEB1_log", 2], p = s2["ZEB1_log", 4], n = nobs(m2))
    }
  }
  if ("etp" %in% names(dat) && uniqueN(dat$etp) > 1) {
    f3 <- as.formula(paste0("`", d, "` ~ ZEB1_log + etp"))
    m3 <- lm(f3, data = dat[!is.na(etp)])
    s3 <- summary(m3)$coefficients
    if ("ZEB1_log" %in% rownames(s3)) {
      coef_rows[[length(coef_rows) + 1]] <- data.table(drug = d, model = "ZEB1+ETP",
        estimate = s3["ZEB1_log", 1], se = s3["ZEB1_log", 2], p = s3["ZEB1_log", 4], n = nobs(m3))
    }
    f4 <- as.formula(paste0("`", d, "` ~ ZEB1_log * etp"))
    m4 <- try(lm(f4, data = dat[!is.na(etp)]), silent = TRUE)
    if (!inherits(m4, "try-error")) {
      s4 <- summary(m4)$coefficients
      hit <- grep("ZEB1_log:etp", rownames(s4), value = TRUE)
      if (length(hit)) {
        coef_rows[[length(coef_rows) + 1]] <- data.table(drug = d, model = "ZEB1xETP",
          estimate = s4[hit[1], 1], se = s4[hit[1], 2], p = s4[hit[1], 4], n = nobs(m4))
      }
    }
  }
  if (has_mass) {
    m5 <- try(MASS::rlm(f1, data = dat), silent = TRUE)
    if (!inherits(m5, "try-error")) {
      sm <- summary(m5)$coefficients
      coef_rows[[length(coef_rows) + 1]] <- data.table(drug = d, model = "rlm_ZEB1",
        estimate = sm["ZEB1_log", 1], se = sm["ZEB1_log", 2],
        p = 2 * pnorm(-abs(sm["ZEB1_log", 3])), n = nrow(dat))
    }
  }
}
coef_dt <- rbindlist(coef_rows)
if (nrow(coef_dt)) {
  coef_dt[, fdr := p.adjust(p, "BH"), by = model]
  fwrite(coef_dt, file.path(tab, "M7_7.5_7.13_regressions.tsv"), sep = "\t")
}

# 7.9 heatmap of available LC50 (column-scaled), T-ALL
if (has_cht) {
  suppressPackageStartupMessages(library(ComplexHeatmap))
  mat <- as.matrix(ph[, ..raw_drugs])
  rownames(mat) <- ph$sample
  colnames(mat) <- short_drug(raw_drugs)
  keep_row <- rowSums(is.finite(mat)) >= 3
  mat <- mat[keep_row, , drop = FALSE]
  mat_z <- apply(mat, 2, function(x) {
    s <- sd(x, na.rm = TRUE); if (!is.finite(s) || s == 0) return(rep(NA_real_, length(x)))
    (x - mean(x, na.rm = TRUE)) / s
  })
  mat_z[!is.finite(mat_z)] <- 0
  ha <- HeatmapAnnotation(
    ZEB1 = ph$ZEB1_log[match(rownames(mat_z), ph$sample)],
    ETP = ph$etp[match(rownames(mat_z), ph$sample)],
    col = list(
      ZEB1 = colorRamp2_safe(range(ph$ZEB1_log, na.rm = TRUE), c("#3C5488", "#E64B35")),
      ETP = pal_etp
    ),
    annotation_name_side = "left"
  )
  ht <- Heatmap(t(mat_z), name = "LC50 z", na_col = "grey90",
                top_annotation = ha, cluster_columns = TRUE, cluster_rows = TRUE,
                column_title = "T-ALL patients (missing LC50 set to z=0 for clustering)",
                row_names_gp = gpar(fontsize = 9),
                show_column_names = FALSE)
  pdf(file.path(fig, "M7_7.9_drug_heatmap.pdf"), width = 9.5, height = 6.8)
  draw(ht); dev.off()
  png(file.path(fig, "M7_7.9_drug_heatmap.png"), width = 9.5 * 150, height = 6.8 * 150, res = 150)
  draw(ht); dev.off()
}

# 7.18 MRD
mrd_cols <- grep("MRD", names(ph), value = TRUE)
mrd_rows <- rbindlist(lapply(mrd_cols, function(mc) {
  spearman_row("ZEB1_log", mc, ph$ZEB1_log, as.numeric(ph[[mc]]))
}))
if (nrow(mrd_rows)) fwrite(mrd_rows, file.path(tab, "M7_7.18_ZEB1_MRD.tsv"), sep = "\t")
if (length(mrd_cols)) {
  mc <- mrd_cols[1]
  dd <- ph[is.finite(as.numeric(get(mc))) & is.finite(ZEB1_log)]
  if (nrow(dd) >= 8) {
    dd[, mrd := as.numeric(dd[[mc]])]
    ct <- cor.test(dd$ZEB1_log, dd$mrd, method = "spearman", exact = FALSE)
    p18 <- ggplot(dd, aes(x = ZEB1_log, y = mrd + 1e-3)) +
      geom_point(alpha = 0.7, size = 1.8, color = "#3C5488") +
      scale_y_log10() +
      geom_smooth(method = "lm", formula = y ~ x, color = "grey20", se = TRUE) +
      labs(title = sprintf("ZEB1 vs %s", mc),
           subtitle = sprintf("n = %d; rho = %.2f; %s", nrow(dd), unname(ct$estimate), fmt_p(ct$p.value)),
           x = "log2(FPKM+1) ZEB1", y = paste(mc, "(log10 scale)")) +
      theme_sci()
    save_both(p18, file.path(fig, "M7_7.18_ZEB1_MRD.pdf"), 6.0, 5.0)
  }
}

# 7.19 exploratory: ZEB1 -> MRD, drug -> MRD, ZEB1+drug -> MRD for top |rho| drug
if (nrow(zeb) && length(mrd_cols)) {
  dtop <- as.character(zeb[which.min(p), drug])
  mc <- mrd_cols[1]
  dat <- ph[is.finite(ZEB1_log) & is.finite(get(dtop)) & is.finite(as.numeric(get(mc)))]
  if (nrow(dat) >= 20) {
    dat[, mrd := as.numeric(dat[[mc]])]
    m_a <- lm(mrd ~ ZEB1_log, data = dat)
    m_b <- lm(as.formula(paste0("`", dtop, "` ~ ZEB1_log")), data = dat)
    m_c <- lm(as.formula(paste0("mrd ~ ZEB1_log + `", dtop, "`")), data = dat)
    med <- rbind(
      data.table(step = "MRD~ZEB1", t(summary(m_a)$coefficients["ZEB1_log", ])),
      data.table(step = paste0(short_drug(dtop), "~ZEB1"), t(summary(m_b)$coefficients["ZEB1_log", ])),
      data.table(step = "MRD~ZEB1+drug_ZEB1", t(summary(m_c)$coefficients["ZEB1_log", ]))
    )
    fwrite(med, file.path(tab, "M7_7.19_exploratory_mediation.tsv"), sep = "\t")
    writeLines(c(
      "Exploratory three-regression mediation only. Not a causal claim.",
      sprintf("Drug used: %s (smallest Spearman P vs ZEB1).", dtop),
      "Public observational LC50/MRD cannot establish that ZEB1 acts through this drug."
    ), file.path(tab, "M7_7.19_note.txt"))
  }
}

writeLines(c(
  "7.20 independent cell-line drug validation is deferred to Module 8 (DepMap/PRISM).",
  "Current DepMap download has CRISPR GeneEffect + expression, not a complete PRISM matrix in this project folder."
), file.path(tab, "M7_7.20_note.txt"))

writeLines(capture.output(sessionInfo()), file.path(tab, "sessionInfo.txt"))
cat("MODULE7 ANALYZE DONE\n")
cat("T-ALL n", nrow(ph), "drugs", length(raw_drugs), "\n")
print(zeb[order(p)][1:5, .(drug_short, n, rho, p, fdr)])

# Module 4 leftovers that are still statistically distinct:
# adjusted cores, GAM for key genes, TARGET vs Pharmacotype rank concordance.
# decoupleR / WGCNA / EnrichmentMap not installed; Hallmark GSEA already empty.
suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(mgcv)
})
has_ggrepel <- requireNamespace("ggrepel", quietly = TRUE)

root <- "/data-b/liangfuhua/projects/ZEB1_TALL_analysis/module4"
m3 <- "/data-b/liangfuhua/projects/ZEB1_TALL_analysis/module3/processed"
fig <- file.path(root, "figures"); tab <- file.path(root, "tables"); proc <- file.path(root, "processed")
dir.create(fig, FALSE, TRUE)

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
    png(sub("\\.pdf$", ".png", file), width = width * 160, height = height * 160, res = 160)
    print(plot); dev.off(); message("saved ", basename(file))
  }, error = function(e) message("FAIL ", basename(file), ": ", conditionMessage(e)))
}
run <- function(label, fn) tryCatch({ fn(); message("OK ", label) },
                                    error = function(e) message("FAIL ", label, ": ", conditionMessage(e)))
read_mat <- function(path) {
  d <- as.data.frame(fread(path), stringsAsFactors = FALSE)
  rownames(d) <- d[[1]]; d[[1]] <- NULL
  as.matrix(d)
}

st <- fread(file.path(tab, "M4_r2_subtype_adjusted_ZEB1.tsv"))
pos <- st[fdr_adj < 0.05 & t_adj > 0][order(-t_adj)][1:80]
neg <- st[fdr_adj < 0.05 & t_adj < 0][order(t_adj)][1:80]
fwrite(pos, file.path(tab, "M4_r2_pos_core80.tsv"), sep = "\t")
fwrite(neg, file.path(tab, "M4_r2_neg_core80.tsv"), sep = "\t")
writeLines(pos$gene, file.path(proc, "m4_adj_pos_core.txt"))
writeLines(neg$gene, file.path(proc, "m4_adj_neg_core.txt"))

# rank concordance TARGET unadjusted vs Pharmacotype (RRHO substitute)
tstat <- fread(file.path(proc, "target_zeb1_gene_stats.tsv"))
pstat <- fread(file.path(proc, "pharmacotype_zeb1_gene_stats.tsv"))
if ("gene" %in% names(tstat) && "gene" %in% names(pstat)) {
  m <- merge(tstat[, .(gene, rho_zeb1, fdr_rho)],
             pstat[, .(gene, rho_ph = rho_zeb1, fdr_ph = fdr_rho)], by = "gene")
  m <- m[is.finite(rho_zeb1) & is.finite(rho_ph)]
  ct <- cor.test(m$rho_zeb1, m$rho_ph, method = "spearman", exact = FALSE)
  fwrite(data.table(n = nrow(m), spearman = unname(ct$estimate), p = ct$p.value),
         file.path(tab, "M4_r2_TARGET_vs_Pharmacotype_rank.tsv"), sep = "\t")
  run("4r2 rank concordance", function() {
    set.seed(1)
    show <- m[sample(.N, min(.N, 4000))]
    p <- ggplot(show, aes(rho_zeb1, rho_ph)) +
      geom_point(size = 0.35, alpha = 0.35, color = "#3C5488") +
      geom_hline(yintercept = 0, linetype = 2, color = "grey70") +
      geom_vline(xintercept = 0, linetype = 2, color = "grey70") +
      geom_smooth(method = "lm", se = FALSE, linewidth = 0.5, color = "#E64B35") +
      labs(title = "TARGET vs Pharmacotype: gene-wise ZEB1 Spearman",
           subtitle = sprintf("n = %d genes; rho = %.2f, P = %.1e. RRHO2 not installed.",
                              nrow(m), unname(ct$estimate), ct$p.value),
           x = "TARGET rho", y = "Pharmacotype rho") + theme_sci()
    save_both(p, file.path(fig, "M4_r2_4.5_cohort_rank_concordance.pdf"), 6.4, 6.0)
  })
}

# GAM for key genes (nonlinear leftover)
expr <- read_mat(file.path(proc, "target_log2tpm_filtered.tsv.gz"))
lab <- fread(file.path(m3, "target_tall_subtype.tsv"))
lab <- lab[usi %in% colnames(expr) & subtype != "Unknown" & !is.na(subtype) & subtype != ""]
expr <- expr[, lab$usi, drop = FALSE]
lab[, ZEB1_z := as.numeric(scale(log2(ZEB1 + 1)))]
keys <- intersect(c("LMO2", "LYL1", "GATA3", "TCF7", "BCL11B", "MEF2C", "IL7R", "ZEB2"), rownames(expr))
gam_rows <- list()
pred_all <- list()
for (g in keys) {
  dat <- data.table(y = as.numeric(expr[g, ]), ZEB1_z = lab$ZEB1_z, subtype = factor(lab$subtype))
  dat <- dat[is.finite(y) & is.finite(ZEB1_z)]
  fit <- tryCatch(gam(y ~ s(ZEB1_z, k = 5) + subtype, data = dat, method = "REML"), error = function(e) NULL)
  if (is.null(fit)) next
  sm <- summary(fit)
  edf <- tryCatch(sm$edf[1], error = function(e) NA_real_)
  pval <- tryCatch(sm$s.pv[1], error = function(e) NA_real_)
  gam_rows[[g]] <- data.table(gene = g, n = nrow(dat), edf = edf, s_p = pval,
                              r2 = sm$r.sq, deviance_expl = sm$dev.expl)
  nd <- data.table(ZEB1_z = seq(min(dat$ZEB1_z), max(dat$ZEB1_z), length.out = 80),
                   subtype = names(sort(table(dat$subtype), decreasing = TRUE))[1])
  pr <- predict(fit, nd, se.fit = TRUE)
  pred_all[[g]] <- data.table(gene = g, ZEB1_z = nd$ZEB1_z, fit = as.numeric(pr$fit),
                              se = as.numeric(pr$se.fit), y = NA_real_, ZEB1_obs = NA_real_)
  pred_all[[g]] <- rbind(
    pred_all[[g]],
    data.table(gene = g, ZEB1_z = dat$ZEB1_z, fit = NA_real_, se = NA_real_,
               y = dat$y, ZEB1_obs = dat$ZEB1_z)
  )
}
gamt <- rbindlist(gam_rows)
fwrite(gamt, file.path(tab, "M4_r2_keygene_GAM.tsv"), sep = "\t")
pred <- rbindlist(pred_all)
fwrite(pred, file.path(tab, "M4_r2_keygene_GAM_curves.tsv"), sep = "\t")

run("4r2 GAM panels", function() {
  if (!nrow(pred)) stop("no GAM")
  pts <- pred[is.finite(y)]
  crv <- pred[is.finite(fit)]
  p <- ggplot() +
    geom_point(data = pts, aes(ZEB1_obs, y), size = 0.35, alpha = 0.25, color = "grey50") +
    geom_ribbon(data = crv, aes(ZEB1_z, ymin = fit - 1.96 * se, ymax = fit + 1.96 * se),
                fill = "#E64B35", alpha = 0.18) +
    geom_line(data = crv, aes(ZEB1_z, fit), color = "#E64B35", linewidth = 0.7) +
    facet_wrap(~gene, scales = "free_y", ncol = 4) +
    labs(title = "GAM: key genes vs ZEB1 after subtype offset",
         subtitle = "s(ZEB1_z, k=5) + subtype; TARGET. LMO2 remaining after subtype is weak/nonlinear check.",
         x = "ZEB1 z", y = "log2 TPM") + theme_sci()
  save_both(p, file.path(fig, "M4_r2_4.25_keygene_GAM.pdf"), 10.2, 6.2)
})

writeLines(c(
  "RRHO2 / decoupleR / WGCNA / EnrichmentMap not installed.",
  "Hallmark GSEA after subtype adjustment was already empty (padj none < 0.05).",
  "Rank concordance scatter is the cross-cohort substitute for RRHO.",
  "GAM tests residual nonlinearity after a subtype offset; it is not a second genome-wide discovery."
), file.path(tab, "M4_r2_leftover_note.txt"))
message("M4 leftover DONE")

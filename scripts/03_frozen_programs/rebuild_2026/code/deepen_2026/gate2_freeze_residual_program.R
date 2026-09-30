# Gate 2: freeze ZEB1-side50 / ZEB2-side50 from effect size, not FDR rank.
# Residual is z-scored so beta is voom log2 change per 1 SD of residual.
# ZEB1/ZEB2 t-statistics are not used as biological QC.
suppressPackageStartupMessages({
  library(data.table); library(edgeR); library(limma); library(splines)
})
argv <- grep("^--file=", commandArgs(FALSE), value = TRUE)
rebuild <- normalizePath(file.path(dirname(normalizePath(sub("^--file=", "", argv), winslash = "/")), "../.."), winslash = "/")
root <- normalizePath(file.path(rebuild, "../.."), winslash = "/")
out <- file.path(rebuild, "data/deepen_2026/gate2_program")
dir.create(out, recursive = TRUE, showWarnings = FALSE)

exclude <- toupper(c("ZEB1", "ZEB2", "CD34", "LYL1", "CD1A"))
fdr_cut <- 0.01
beta_cut <- 0.25
n_side <- 50L

pt <- fread(file.path(root, "total/data/validation/zeb_developmental_residual/polonen_patient_residual.tsv"),
            encoding = "UTF-8")
stopifnot(nrow(pt) == 1309L, uniqueN(pt$subtype) == 17L)
pt[, subtype := factor(subtype)]
pt[, residual_z := as.numeric(scale(residual_within))]
stopifnot(all(is.finite(pt$residual_z)))

cnt <- fread(file.path(root, "data/TALL_X01_counts.tsv"), nrows = 0)
ids <- pt$sample_id
stopifnot(all(ids %in% names(cnt)))
z <- fread(file.path(root, "data/TALL_X01_counts.tsv"), select = c(names(cnt)[1], ids))
gene_id <- z[[1]]
mat <- as.matrix(z[, ..ids]); storage.mode(mat) <- "integer"
rownames(mat) <- gene_id
e <- new.env(parent = emptyenv())
load(file.path(root, "data/polonen_syn54032669/Data_1309Samples.RData"), envir = e)
symbols <- rownames(e$gexp)
stopifnot(length(symbols) == nrow(mat))
rm(e, z, cnt); gc()

dev_ns <- ns(pt$dev, df = 3)
colnames(dev_ns) <- paste0("ns", seq_len(ncol(dev_ns)))
design <- model.matrix(~ dev_ns + subtype + residual_z, data = pt)
y <- DGEList(counts = mat)
keep <- filterByExpr(y, design = design)
y <- y[keep, , keep.lib.sizes = FALSE]
sym <- symbols[keep]
y <- calcNormFactors(y, method = "TMM")
v <- voom(y, design, plot = FALSE)
fit <- eBayes(lmFit(v, design), robust = TRUE)
tt <- topTable(fit, coef = "residual_z", number = Inf, sort.by = "none")
df_res <- fit$df.residual[1]
res <- as.data.table(tt)
res[, `:=`(gene_id = rownames(tt), symbol = sym,
           beta_R = logFC,
           partial_R2 = t^2 / (t^2 + df_res))]
res <- res[, .(gene_id, symbol, AveExpr, beta_R, t, P.Value, adj.P.Val, partial_R2)]
setorder(res, P.Value)
fwrite(res, file.path(out, "all_genes_residual_effects.tsv"), sep = "\t")

elig <- res[!(toupper(symbol) %in% exclude)]
pass <- elig[adj.P.Val < fdr_cut & abs(beta_R) >= beta_cut]
zeb1 <- pass[beta_R > 0][order(-partial_R2)]
zeb2 <- pass[beta_R < 0][order(-partial_R2)]
side1 <- zeb1[seq_len(min(n_side, .N))][, `:=`(side = "ZEB1-side50", rank = .I)]
side2 <- zeb2[seq_len(min(n_side, .N))][, `:=`(side = "ZEB2-side50", rank = .I)]
freeze <- rbind(side1, side2)
fwrite(pass, file.path(out, "effect_threshold_pass.tsv"), sep = "\t")
fwrite(freeze, file.path(out, "frozen_ZEB_side50.tsv"), sep = "\t")

# Technical axes: residual vs library size / TMM / expression PCs.
tech <- data.table(
  sample_id = ids,
  residual_z = pt$residual_z,
  log_lib_size = log10(y$samples$lib.size),
  tmm_factor = y$samples$norm.factors
)
pc <- prcomp(t(v$E), center = TRUE, scale. = TRUE, rank. = 5)
tech[, paste0("PC", 1:5) := as.data.table(pc$x[, 1:5])]
tech_cor <- rbindlist(lapply(c("log_lib_size", "tmm_factor", paste0("PC", 1:5)), function(nm) {
  ct <- suppressWarnings(cor.test(tech$residual_z, tech[[nm]], method = "spearman"))
  data.table(axis = nm, spearman = unname(ct$estimate), p = ct$p.value)
}))
fwrite(tech, file.path(out, "residual_vs_technical_axes.tsv"), sep = "\t")
fwrite(tech_cor, file.path(out, "residual_technical_spearman.tsv"), sep = "\t")

# Sensitivity: add PC1-5, compare residual ranking among eligible genes.
design_pc <- cbind(design, pc$x[, 1:5])
colnames(design_pc)[(ncol(design) + 1):ncol(design_pc)] <- paste0("PC", 1:5)
v2 <- voom(y, design_pc, plot = FALSE)
fit2 <- eBayes(lmFit(v2, design_pc), robust = TRUE)
tt2 <- topTable(fit2, coef = "residual_z", number = Inf, sort.by = "none")
cmp <- data.table(symbol = sym, t_primary = res$t, t_pc = tt2$t)
cmp <- cmp[!(toupper(symbol) %in% exclude)]
rho <- suppressWarnings(cor.test(cmp$t_primary, cmp$t_pc, method = "spearman"))
sens <- data.table(
  n_compared = nrow(cmp),
  spearman_t_primary_vs_pc = unname(rho$estimate),
  p = rho$p.value
)
fwrite(cmp, file.path(out, "residual_t_primary_vs_pc15.tsv"), sep = "\t")
fwrite(sens, file.path(out, "residual_pc_sensitivity.tsv"), sep = "\t")

summ <- data.table(
  n_patients = nrow(pt),
  n_tested = nrow(res),
  n_eligible = nrow(elig),
  fdr_cut = fdr_cut,
  beta_cut = beta_cut,
  n_pass = nrow(pass),
  n_pass_ZEB1_side = nrow(zeb1),
  n_pass_ZEB2_side = nrow(zeb2),
  n_frozen_ZEB1 = nrow(side1),
  n_frozen_ZEB2 = nrow(side2),
  residual_vs_PC1_rho = tech_cor[axis == "PC1", spearman],
  residual_vs_tmm_rho = tech_cor[axis == "tmm_factor", spearman],
  t_rank_stability_with_PC1to5 = sens$spearman_t_primary_vs_pc
)
fwrite(summ, file.path(out, "freeze_summary.tsv"), sep = "\t")
writeLines(c(
  "Primary model: voom/limma E ~ ns(dev,3) + subtype + residual_z.",
  "residual_z is the frozen development-adjusted ZEB balance, standardized to unit variance.",
  "beta_R is the voom log2 coefficient per 1 SD residual. partial R2 = t^2/(t^2+df).",
  "Excluded from program: ZEB1, ZEB2, CD34, LYL1, CD1A.",
  "ZEB1/ZEB2 coefficients are mathematical consequences of the residual definition and are not QC.",
  "Eligibility: BH FDR < 0.01 and |beta_R| >= 0.25. Freeze: top 50 partial R2 per side.",
  "Thresholds were locked before looking at the freeze lists.",
  "PC1-5 enter only a sensitivity model. They are not in the primary freeze.",
  "OXPHOS is not interpreted as a central program."
), file.path(out, "METHODS.txt"))
writeLines(capture.output(sessionInfo()), file.path(out, "R_sessionInfo.txt"))
cat("Gate 2 freeze complete\n")
print(summ)

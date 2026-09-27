# A7-0: specimen (BM vs blood) sensitivity on frozen 100 genes. No reselection.
suppressPackageStartupMessages({
  library(data.table); library(edgeR); library(limma); library(splines)
})
rebuild <- "D:/_bioinformation/ZEB1/total/rebuild_2026"
root <- "D:/_bioinformation/ZEB1"
out <- file.path(rebuild, "data/deepen_2026/a7_lineage")
dir.create(out, recursive = TRUE, showWarnings = FALSE)

fr <- fread(file.path(rebuild, "data/deepen_2026/gate2_program/frozen_ZEB_side50.tsv"))
prim <- fread(file.path(rebuild, "data/deepen_2026/gate2_program/all_genes_residual_effects.tsv"))
pt <- fread(file.path(root, "total/data/validation/zeb_developmental_residual/polonen_patient_residual.tsv"), encoding = "UTF-8")
batch <- fread(file.path(rebuild, "data/deepen_2026/gate2_program/tall_x01_sample_batch.tsv"))
clin <- fread(file.path(root, "data/polonen_syn54032669/audit_clinical_1309.tsv"), select = c("sample_id", "specimen"))
pt <- merge(pt, batch, by = "sample_id", all.x = TRUE, sort = FALSE)
pt <- merge(pt, clin, by = "sample_id", all.x = TRUE, sort = FALSE)
stopifnot(nrow(pt) == 1309L, !anyNA(pt$batch), !anyNA(pt$specimen))
pt[, `:=`(
  subtype = factor(subtype),
  batch = factor(batch),
  specimen = factor(specimen, levels = c("blood", "bone marrow")),
  residual_z = as.numeric(scale(residual_within))
)]
fwrite(pt[, .N, by = specimen], file.path(out, "a7_0_specimen_counts.tsv"), sep = "\t")

ids <- pt$sample_id
cnt <- fread(file.path(root, "data/TALL_X01_counts.tsv"), nrows = 0)
z <- fread(file.path(root, "data/TALL_X01_counts.tsv"), select = c(names(cnt)[1], ids))
mat <- as.matrix(z[, ..ids]); storage.mode(mat) <- "integer"
rownames(mat) <- z[[1]]
e <- new.env(parent = emptyenv())
load(file.path(root, "data/polonen_syn54032669/Data_1309Samples.RData"), envir = e)
symbols <- rownames(e$gexp)
rm(e); gc()

dev_ns <- ns(pt$dev, df = 3)
colnames(dev_ns) <- paste0("ns", seq_len(ncol(dev_ns)))
design <- model.matrix(~ dev_ns + subtype + batch + specimen + residual_z, data = pt)
y <- DGEList(counts = mat)
keep <- filterByExpr(y, design = design)
y <- y[keep, , keep.lib.sizes = FALSE]
sym <- symbols[keep]
y <- calcNormFactors(y, method = "TMM")
v <- voom(y, design, plot = FALSE)
fit <- eBayes(lmFit(v, design), robust = TRUE)
tt <- as.data.table(topTable(fit, coef = "residual_z", number = Inf, sort.by = "none"))
tt[, symbol := sym]
adj <- merge(fr[, .(symbol, side, rank, beta_primary = beta_R)],
             tt[, .(symbol, beta_specimen = logFC, t_specimen = t, p_specimen = P.Value)],
             by = "symbol")
adj[, `:=`(same_sign = sign(beta_primary) == sign(beta_specimen),
           attenuation = 1 - beta_specimen / beta_primary)]
fwrite(adj, file.path(out, "a7_0_frozen_genes_specimen_adjusted.tsv"), sep = "\t")

by_side <- adj[, .(
  n = .N,
  n_same_sign = sum(same_sign),
  spearman_beta = suppressWarnings(cor(beta_primary, beta_specimen, method = "spearman")),
  median_attenuation = median(attenuation)
), by = side]
overall <- data.table(
  n_frozen_recovered = nrow(adj),
  n_same_sign = sum(adj$same_sign),
  spearman_beta = suppressWarnings(cor(adj$beta_primary, adj$beta_specimen, method = "spearman")),
  median_attenuation = median(adj$attenuation)
)
fwrite(by_side, file.path(out, "a7_0_by_side.tsv"), sep = "\t")
fwrite(overall, file.path(out, "a7_0_overall.tsv"), sep = "\t")

# frozen side scores on this voom E
E <- v$E
rownames(E) <- sym
zscore <- function(x) {
  s <- sd(x)
  if (!is.finite(s) || s == 0) return(rep(NA_real_, length(x)))
  (x - mean(x)) / s
}
score_side <- function(genes) {
  hit <- intersect(genes, rownames(E))
  if (!length(hit)) return(rep(NA_real_, ncol(E)))
  Z <- t(apply(E[hit, , drop = FALSE], 1, zscore))
  colMeans(Z, na.rm = TRUE)
}
sc <- data.table(
  sample_id = ids,
  ZEB1_side50 = score_side(fr[side == "ZEB1-side50", symbol]),
  ZEB2_side50 = score_side(fr[side == "ZEB2-side50", symbol])
)
sc <- merge(sc, pt[, .(sample_id, dev, subtype, batch, specimen)], by = "sample_id")
fwrite(sc, file.path(out, "a7_0_patient_scores.tsv"), sep = "\t")
dev_ns2 <- ns(sc$dev, df = 3)
colnames(dev_ns2) <- paste0("ns", seq_len(ncol(dev_ns2)))
score_rows <- lapply(c("ZEB1_side50", "ZEB2_side50"), function(nm) {
  fit_s <- lm(sc[[nm]] ~ dev_ns2 + subtype + batch + specimen, data = sc)
  cf <- summary(fit_s)$coefficients
  rn <- grep("specimen", rownames(cf), value = TRUE)
  data.table(
    score = nm,
    specimen_coef = if (length(rn)) unname(cf[rn[1], "Estimate"]) else NA_real_,
    specimen_p = if (length(rn)) unname(cf[rn[1], "Pr(>|t|)"]) else NA_real_,
    specimen_term = if (length(rn)) rn[1] else NA_character_
  )
})
fwrite(rbindlist(score_rows), file.path(out, "a7_0_score_specimen.tsv"), sep = "\t")
print(overall)
print(by_side)
print(rbindlist(score_rows))
cat("A7-0 done\n")

# Two locked sensitivities. Neither retunes the Gate 2 thresholds.
# 1) Discover side50 after dropping every BCL11B and ETP-like patient.
#    Z-scores and the developmental spline are estimated in the 15-subtype
#    training set only. The 18-pair matched contrast is then a held-out test.
# 2) Add known library batch to the full-cohort model. Do not re-select genes.
suppressPackageStartupMessages({
  library(data.table); library(edgeR); library(limma); library(splines)
})
argv <- grep("^--file=", commandArgs(FALSE), value = TRUE)
rebuild <- normalizePath(file.path(dirname(normalizePath(sub("^--file=", "", argv), winslash = "/")), "../.."), winslash = "/")
root <- normalizePath(file.path(rebuild, "../.."), winslash = "/")
out <- file.path(rebuild, "data/deepen_2026/gate2_sensitivity")
dir.create(out, recursive = TRUE, showWarnings = FALSE)

exclude <- toupper(c("ZEB1", "ZEB2", "ZEB1-AS1", "ZEB2-AS1", "CD34", "LYL1", "CD1A"))
fdr_cut <- 0.01
beta_cut <- 0.25
n_side <- 50L
zsd <- function(x, m, s) (as.numeric(x) - m) / s

pt <- fread(file.path(root, "total/data/validation/zeb_developmental_residual/polonen_patient_residual.tsv"),
            encoding = "UTF-8")
stopifnot(nrow(pt) == 1309L)
batch <- fread(file.path(rebuild, "data/deepen_2026/gate2_program/tall_x01_sample_batch.tsv"))
pt <- merge(pt, batch, by = "sample_id", all.x = TRUE, sort = FALSE)
stopifnot(nrow(pt) == 1309L, !anyNA(pt$batch))
held_sub <- c("BCL11B", "ETP-like")
train <- pt[!(subtype %in% held_sub)]
stopifnot(uniqueN(train$subtype) == 15L, nrow(train) == 1309L - 18L - 235L)

# Training-only balance / coordinate / residual.
mu <- sapply(c("ZEB1", "ZEB2", "CD1A", "CD34", "LYL1"), function(g) mean(train[[g]]))
sg <- sapply(c("ZEB1", "ZEB2", "CD1A", "CD34", "LYL1"), function(g) sd(train[[g]]))
stopifnot(all(is.finite(sg)), all(sg > 0))
score <- function(d) {
  d[, `:=`(
    z_ZEB1_tr = zsd(ZEB1, mu[["ZEB1"]], sg[["ZEB1"]]),
    z_ZEB2_tr = zsd(ZEB2, mu[["ZEB2"]], sg[["ZEB2"]]),
    z_CD1A_tr = zsd(CD1A, mu[["CD1A"]], sg[["CD1A"]]),
    z_CD34_tr = zsd(CD34, mu[["CD34"]], sg[["CD34"]]),
    z_LYL1_tr = zsd(LYL1, mu[["LYL1"]], sg[["LYL1"]])
  )]
  d[, `:=`(
    balance_tr = z_ZEB1_tr - z_ZEB2_tr,
    dev_tr = z_CD1A_tr - (z_CD34_tr + z_LYL1_tr) / 2
  )]
}
train <- score(copy(train))
fit_tr <- lm(balance_tr ~ ns(dev_tr, df = 3), data = train)
train[, residual_tr := balance_tr - as.numeric(predict(fit_tr, newdata = train))]
train[, residual_z := as.numeric(scale(residual_tr))]
train[, subtype := droplevels(factor(subtype))]

orig <- fread(file.path(rebuild, "data/deepen_2026/gate2_program/frozen_ZEB_side50.tsv"))
deg <- fread(file.path(rebuild, "data/source_data_rebuilt/revision4_targeted/matched_BCL11B_ETP_gene_results.tsv"))

cnt <- fread(file.path(root, "data/TALL_X01_counts.tsv"), nrows = 0)
e <- new.env(parent = emptyenv())
load(file.path(root, "data/polonen_syn54032669/Data_1309Samples.RData"), envir = e)
symbols_all <- rownames(e$gexp)
rm(e); gc()

# ---- 1. held-out program discovery ----
ids_tr <- train$sample_id
stopifnot(all(ids_tr %in% names(cnt)))
z <- fread(file.path(root, "data/TALL_X01_counts.tsv"), select = c(names(cnt)[1], ids_tr))
mat <- as.matrix(z[, ..ids_tr]); storage.mode(mat) <- "integer"
rownames(mat) <- z[[1]]
stopifnot(length(symbols_all) == nrow(mat))
dev_ns <- ns(train$dev_tr, df = 3)
colnames(dev_ns) <- paste0("ns", seq_len(ncol(dev_ns)))
design_tr <- model.matrix(~ dev_ns + subtype + residual_z, data = train)
y <- DGEList(counts = mat)
keep <- filterByExpr(y, design = design_tr)
y <- y[keep, , keep.lib.sizes = FALSE]
sym <- symbols_all[keep]
y <- calcNormFactors(y, method = "TMM")
v <- voom(y, design_tr, plot = FALSE)
fit <- eBayes(lmFit(v, design_tr), robust = TRUE)
tt0 <- topTable(fit, coef = "residual_z", number = Inf, sort.by = "none")
tt <- as.data.table(tt0)
tt[, `:=`(gene_id = rownames(tt0), symbol = sym, beta_R = logFC,
          partial_R2 = t^2 / (t^2 + fit$df.residual[1]))]
setorder(tt, P.Value)
fwrite(tt[, .(gene_id, symbol, AveExpr, beta_R, t, P.Value, adj.P.Val, partial_R2)],
       file.path(out, "heldout_all_genes.tsv"), sep = "\t")
elig <- tt[!(toupper(symbol) %in% exclude)]
pass <- elig[adj.P.Val < fdr_cut & abs(beta_R) >= beta_cut]
side1 <- pass[beta_R > 0][order(-partial_R2)][seq_len(min(n_side, .N))][, `:=`(side = "ZEB1-side50", rank = .I)]
side2 <- pass[beta_R < 0][order(-partial_R2)][seq_len(min(n_side, .N))][, `:=`(side = "ZEB2-side50", rank = .I)]
held_fr <- rbind(side1, side2)
fwrite(held_fr, file.path(out, "heldout_ZEB_side50.tsv"), sep = "\t")

ov <- rbindlist(lapply(c("ZEB1-side50", "ZEB2-side50"), function(s) {
  a <- orig[side == s, symbol]
  b <- held_fr[side == s, symbol]
  data.table(side = s, n_orig = uniqueN(a), n_heldout = uniqueN(b),
             n_overlap = length(intersect(a, b)))
}))
fwrite(ov, file.path(out, "heldout_vs_original_overlap.tsv"), sep = "\t")

stat <- deg$t
names(stat) <- deg$symbol
stat <- stat[is.finite(stat) & !duplicated(names(stat))]
idx <- lapply(split(held_fr$symbol, held_fr$side), function(s) which(names(stat) %in% unique(s)))
cam <- as.data.table(cameraPR(stat, index = idx), keep.rownames = "program")
cam[, expected := fifelse(program == "ZEB2-side50", "Up", "Down")]
cam[, direction_matches := Direction == expected]
fwrite(cam[, .(program, NGenes, Direction, PValue, FDR, expected, direction_matches)],
       file.path(out, "heldout_matched_cameraPR.tsv"), sep = "\t")
rm(z, mat, y, v, fit, design_tr, dev_ns); gc()

# ---- 2. known-batch sensitivity on the original frozen 100 genes ----
ids <- pt$sample_id
z <- fread(file.path(root, "data/TALL_X01_counts.tsv"), select = c(names(cnt)[1], ids))
mat <- as.matrix(z[, ..ids]); storage.mode(mat) <- "integer"
rownames(mat) <- z[[1]]
pt[, `:=`(subtype = factor(subtype), batch = factor(batch),
          residual_z = as.numeric(scale(residual_within)))]
dev_ns <- ns(pt$dev, df = 3)
colnames(dev_ns) <- paste0("ns", seq_len(ncol(dev_ns)))
design_b <- model.matrix(~ dev_ns + subtype + batch + residual_z, data = pt)
y <- DGEList(counts = mat)
keep <- filterByExpr(y, design = design_b)
y <- y[keep, , keep.lib.sizes = FALSE]
sym <- symbols_all[keep]
y <- calcNormFactors(y, method = "TMM")
v <- voom(y, design_b, plot = FALSE)
fitb <- eBayes(lmFit(v, design_b), robust = TRUE)
ttb <- as.data.table(topTable(fitb, coef = "residual_z", number = Inf, sort.by = "none"))
ttb[, symbol := sym]
prim <- fread(file.path(rebuild, "data/deepen_2026/gate2_program/all_genes_residual_effects.tsv"))
fr <- merge(orig[, .(symbol, side, rank, beta_primary = beta_R)],
            ttb[, .(symbol, beta_batch = logFC, t_batch = t, p_batch = P.Value)],
            by = "symbol")
fr[, `:=`(same_sign = sign(beta_primary) == sign(beta_batch),
          attenuation = 1 - beta_batch / beta_primary)]
fwrite(fr, file.path(out, "frozen_genes_batch_adjusted.tsv"), sep = "\t")

kw <- kruskal.test(residual_z ~ batch, data = pt)
bt <- fr[, .(
  n = .N,
  n_same_sign = sum(same_sign),
  spearman_beta = suppressWarnings(cor(beta_primary, beta_batch, method = "spearman")),
  median_attenuation = median(attenuation),
  median_abs_primary = median(abs(beta_primary)),
  median_abs_batch = median(abs(beta_batch))
), by = side]
bt_all <- data.table(
  residual_vs_batch_kruskal_p = kw$p.value,
  n_frozen_recovered = nrow(fr),
  n_same_sign = sum(fr$same_sign),
  spearman_beta = suppressWarnings(cor(fr$beta_primary, fr$beta_batch, method = "spearman")),
  median_attenuation = median(fr$attenuation)
)
fwrite(as.data.table(pt[, .N, by = batch]), file.path(out, "batch_counts.tsv"), sep = "\t")
fwrite(bt, file.path(out, "batch_sensitivity_by_side.tsv"), sep = "\t")
fwrite(bt_all, file.path(out, "batch_sensitivity_overall.tsv"), sep = "\t")

summ <- data.table(
  train_n = nrow(train),
  heldout_n = nrow(pt) - nrow(train),
  heldout_n_pass = nrow(pass),
  heldout_n_ZEB1 = nrow(side1),
  heldout_n_ZEB2 = nrow(side2),
  overlap_ZEB1 = ov[side == "ZEB1-side50", n_overlap],
  overlap_ZEB2 = ov[side == "ZEB2-side50", n_overlap],
  heldout_ZEB2_FDR = cam[program == "ZEB2-side50", FDR],
  heldout_ZEB2_direction = cam[program == "ZEB2-side50", Direction],
  heldout_ZEB1_FDR = cam[program == "ZEB1-side50", FDR],
  heldout_ZEB1_direction = cam[program == "ZEB1-side50", Direction],
  batch_same_sign = bt_all$n_same_sign,
  batch_spearman = bt_all$spearman_beta,
  batch_median_attenuation = bt_all$median_attenuation,
  residual_batch_p = bt_all$residual_vs_batch_kruskal_p
)
fwrite(summ, file.path(out, "sensitivity_summary.tsv"), sep = "\t")
writeLines(c(
  "Held-out discovery drops every BCL11B and ETP-like patient before z-scores, spline and gene selection.",
  "Expression values remain the deposited TMM logCPM; only standardization and ns(D,3) are training-only.",
  "Thresholds remain FDR<0.01 and |beta_R|>=0.25, top 50 partial R2 per side.",
  "Matched CAMERA PR uses the existing 18-pair t-statistics. Those 36 patients did not enter gene selection.",
  "Batch sensitivity adds the deposited STRANDED_X01 / STRANDED_TARGET / UNSTRANDED_TARGET labels.",
  "Batch-adjusted models do not re-select genes. Unsupervised PCs are not covariates."
), file.path(out, "METHODS.txt"))
writeLines(capture.output(sessionInfo()), file.path(out, "R_sessionInfo.txt"))
cat("Held-out and batch sensitivity complete\n")
print(summ)
print(cam)
print(bt)

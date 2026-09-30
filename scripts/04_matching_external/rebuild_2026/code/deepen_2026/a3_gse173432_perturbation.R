# A3-1: GSE173432 BCL11B OE in human CD34+ cord blood.
# Do not refit developmental residual. Use frozen side50 only.
suppressPackageStartupMessages({
  library(data.table)
  library(limma)
})

argv <- grep("^--file=", commandArgs(FALSE), value = TRUE)
script <- normalizePath(sub("^--file=", "", argv), winslash = "/")
code <- dirname(script)
rebuild <- normalizePath(file.path(code, "../.."), winslash = "/")
out <- file.path(rebuild, "data/deepen_2026/a3_gse165209/gse173432")
dir.create(out, recursive = TRUE, showWarnings = FALSE)

tpm <- fread(file.path(out, "GSE173432_CD34_EV_BCL11B_TPM.txt.gz"))
stopifnot(all(c("geneID", "geneSymbol") %in% names(tpm)))
samples <- c("CD34_EV_1", "CD34_EV_2", "CD34_EV_3",
             "CD34_BCL11B_1", "CD34_BCL11B_2", "CD34_BCL11B_3")
stopifnot(all(samples %in% names(tpm)))

# collapse duplicate symbols by max mean TPM; keep first geneID
tpm[, mean_tpm := rowMeans(.SD), .SDcols = samples]
setorder(tpm, geneSymbol, -mean_tpm)
tpm <- unique(tpm, by = "geneSymbol")

pheno <- data.table(
  sample = samples,
  donor = factor(rep(c("D1", "D2", "D3"), 2)),
  condition = factor(rep(c("EV", "OE"), each = 3), levels = c("EV", "OE"))
)

grab <- function(sym) {
  hit <- tpm[toupper(geneSymbol) == toupper(sym)]
  if (!nrow(hit)) return(rep(NA_real_, length(samples)))
  as.numeric(hit[1, ..samples])
}

key <- data.table(
  sample = samples,
  donor = pheno$donor,
  condition = pheno$condition,
  BCL11B = grab("BCL11B"),
  ZEB1 = grab("ZEB1"),
  ZEB2 = grab("ZEB2"),
  CD34 = grab("CD34"),
  LYL1 = grab("LYL1"),
  CD1A = grab("CD1A")
)
key[, `:=`(
  log2_ZEB1 = log2(ZEB1 + 1),
  log2_ZEB2 = log2(ZEB2 + 1),
  log2_ratio = log2((ZEB1 + 1) / (ZEB2 + 1)),
  log2_CD34 = log2(CD34 + 1),
  log2_LYL1 = log2(LYL1 + 1),
  log2_CD1A = log2(CD1A + 1)
)]
# within-experiment z, then developmental score; not a residual
zsd <- function(x) {
  s <- sd(x, na.rm = TRUE)
  if (!is.finite(s) || s == 0) return(rep(NA_real_, length(x)))
  (x - mean(x, na.rm = TRUE)) / s
}
key[, `:=`(
  z_CD1A = zsd(log2_CD1A),
  z_CD34 = zsd(log2_CD34),
  z_LYL1 = zsd(log2_LYL1)
)]
key[, dev := z_CD1A - (z_CD34 + z_LYL1) / 2]
fwrite(key, file.path(out, "a3_key_gene_tpm.tsv"), sep = "\t")

paired <- dcast(key, donor ~ condition,
                value.var = c("BCL11B", "ZEB1", "ZEB2", "log2_ratio", "dev",
                              "CD34", "LYL1", "CD1A", "log2_CD34", "log2_LYL1", "log2_CD1A"))
paired[, `:=`(
  d_log2_ratio = log2_ratio_OE - log2_ratio_EV,
  d_dev = dev_OE - dev_EV,
  d_ZEB1 = ZEB1_OE - ZEB1_EV,
  d_ZEB2 = ZEB2_OE - ZEB2_EV,
  d_BCL11B = BCL11B_OE - BCL11B_EV,
  d_CD1A = CD1A_OE - CD1A_EV,
  d_CD34 = CD34_OE - CD34_EV,
  d_LYL1 = LYL1_OE - LYL1_EV
)]
fwrite(paired, file.path(out, "a3_paired_donor_deltas.tsv"), sep = "\t")

ttest_one <- function(item, x) {
  x <- x[is.finite(x)]
  if (length(x) < 2) {
    return(data.table(item = item, n = length(x), mean_delta = NA_real_, t = NA_real_, p = NA_real_))
  }
  h <- t.test(x)
  data.table(item = item, n = length(x), mean_delta = unname(h$estimate),
             t = unname(h$statistic), p = unname(h$p.value))
}
anchor <- rbindlist(list(
  ttest_one("paired_d_log2_ZEB1_over_ZEB2", paired$d_log2_ratio),
  ttest_one("paired_d_dev", paired$d_dev),
  ttest_one("paired_d_ZEB1_tpm", paired$d_ZEB1),
  ttest_one("paired_d_ZEB2_tpm", paired$d_ZEB2),
  ttest_one("paired_d_BCL11B_tpm", paired$d_BCL11B),
  ttest_one("paired_d_CD1A_tpm", paired$d_CD1A),
  ttest_one("paired_d_CD34_tpm", paired$d_CD34),
  ttest_one("paired_d_LYL1_tpm", paired$d_LYL1)
))
# expected for Gate A: BCL11B leukemia is ZEB2-side, so log2(ZEB1/ZEB2) should fall
anchor[, expected_for_ZEB2_side_state := fifelse(
  item == "paired_d_log2_ZEB1_over_ZEB2", "negative",
  fifelse(item == "paired_d_BCL11B_tpm", "positive_OE_check", "recorded")
)]
fwrite(anchor, file.path(out, "a3_paired_anchor_tests.tsv"), sep = "\t")

# limma on log2(TPM+1), paired by donor
mat <- log2(as.matrix(tpm[, ..samples]) + 1)
rownames(mat) <- tpm$geneSymbol
design <- model.matrix(~ donor + condition, data = pheno)
fit <- eBayes(lmFit(mat, design), robust = TRUE)
tt0 <- topTable(fit, coef = "conditionOE", number = Inf, sort.by = "none")
tt <- as.data.table(tt0)
tt[, symbol := rownames(tt0)]
fwrite(tt[, .(symbol, logFC, AveExpr, t, P.Value, adj.P.Val)],
       file.path(out, "a3_limma_OE_vs_EV.tsv"), sep = "\t")

frozen <- fread(file.path(rebuild, "data/deepen_2026/gate2_program/frozen_ZEB_side50.tsv"))
idx <- match(toupper(tt$symbol), toupper(frozen$symbol))
tt[, side := frozen$side[idx]]
row_idx <- match(toupper(frozen$symbol), toupper(rownames(mat)))
sets <- list(
  `ZEB2-side50` = which(tt$side == "ZEB2-side50"),
  `ZEB1-side50` = which(tt$side == "ZEB1-side50")
)
sets_mat <- list(
  `ZEB2-side50` = na.omit(row_idx[frozen$side == "ZEB2-side50"]),
  `ZEB1-side50` = na.omit(row_idx[frozen$side == "ZEB1-side50"])
)
n_in <- data.table(
  program = names(sets),
  n_in_matrix = lengths(sets),
  n_frozen = c(frozen[side == "ZEB2-side50", .N], frozen[side == "ZEB1-side50", .N])
)
fwrite(n_in, file.path(out, "a3_frozen_gene_recovery.tsv"), sep = "\t")

cam <- cameraPR(tt$t, sets)
cam <- as.data.table(cam, keep.rownames = "program")
cam[, expected := fifelse(program == "ZEB2-side50", "Up", "Down")]
cam[, direction_matches := Direction == expected]
fwrite(cam, file.path(out, "a3_cameraPR_frozen_side50.tsv"), sep = "\t")

rst <- tryCatch({
  as.data.table(mroast(mat, sets_mat, design, contrast = ncol(design), nrot = 9999),
                keep.rownames = "program")
}, error = function(e) data.table(program = NA_character_, error = conditionMessage(e)))
fwrite(rst, file.path(out, "a3_mroast_frozen_side50.tsv"), sep = "\t")

# Gate A/B lock file
gate_a_delta <- paired[, mean(d_log2_ratio)]
gate_a_p <- anchor[item == "paired_d_log2_ZEB1_over_ZEB2", p]
gate_a <- if (is.finite(gate_a_delta) && gate_a_delta < 0 && is.finite(gate_a_p) && gate_a_p < 0.05) {
  "PASS_toward_ZEB2_side"
} else if (is.finite(gate_a_delta) && gate_a_delta < 0) {
  "DIRECTION_ONLY_toward_ZEB2_side"
} else if (is.finite(gate_a_delta) && gate_a_delta > 0) {
  "FAIL_toward_ZEB1_side"
} else {
  "INDETERMINATE"
}
cam2 <- cam[program == "ZEB2-side50"]
gate_b <- if (nrow(cam2) && cam2$Direction == "Up" && cam2$FDR < 0.05) {
  "PASS"
} else if (nrow(cam2) && cam2$Direction == "Up") {
  "DIRECTION_ONLY"
} else {
  "FAIL"
}

summary <- data.table(
  item = c(
    "n_genes_tpm", "n_donors", "BCL11B_EV_mean", "BCL11B_OE_mean",
    "mean_d_log2_ZEB1_over_ZEB2", "paired_p_log2_ratio",
    "mean_d_dev", "paired_p_dev",
    "ZEB2_side50_recovered", "ZEB2_side50_Direction", "ZEB2_side50_P", "ZEB2_side50_FDR",
    "ZEB1_side50_recovered", "ZEB1_side50_Direction", "ZEB1_side50_P", "ZEB1_side50_FDR",
    "GateA_relative_ZEB", "GateB_frozen_ZEB2_side50",
    "interpretation_note"
  ),
  value = c(
    nrow(tpm), 3,
    key[condition == "EV", mean(BCL11B)], key[condition == "OE", mean(BCL11B)],
    gate_a_delta, gate_a_p,
    paired[, mean(d_dev)], anchor[item == "paired_d_dev", p],
    n_in[program == "ZEB2-side50", n_in_matrix], cam2$Direction, cam2$PValue, cam2$FDR,
    n_in[program == "ZEB1-side50", n_in_matrix],
    cam[program == "ZEB1-side50", Direction],
    cam[program == "ZEB1-side50", PValue],
    cam[program == "ZEB1-side50", FDR],
    gate_a, gate_b,
    "If ZEB2-side moves with little dev shift: reconfiguration. If both move: lineage transition. Do not hide the latter."
  )
)
fwrite(summary, file.path(out, "a3_gse173432_gateAB_summary.tsv"), sep = "\t")
print(key)
print(paired)
print(anchor)
print(cam)
print(summary)
cat("A3-1 GSE173432 complete\n")

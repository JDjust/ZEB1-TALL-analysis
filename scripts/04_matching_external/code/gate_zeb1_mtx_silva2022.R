# Gate: independent ZEB1-MTX test in Canevarolo et al. 2022 (GSE218348).
# Phenotype was locked in silva2022_mtx_ic50.tsv before this script was run.
# Primary question is T-ALL only (n = 7). BCP-ALL is reported separately.

suppressPackageStartupMessages({
  library(affy)
  library(hgu133plus2.db)
  library(AnnotationDbi)
})

val <- "d:/_bioinformation/ZEB1/total/data/validation"
cel_dir <- file.path(val, "GSE218348_CEL")
cels <- list.files(cel_dir, pattern = "CEL\\.gz$", full.names = TRUE)
stopifnot(length(cels) == 13)

raw <- ReadAffy(filenames = cels)
eset <- rma(raw)
mat <- exprs(eset)
colnames(mat) <- sub("\\.CEL\\.gz$", "", basename(colnames(mat)))
colnames(mat) <- sub("^GSM[0-9]+_", "", colnames(mat))

sym <- mapIds(hgu133plus2.db, keys = rownames(mat),
              column = "SYMBOL", keytype = "PROBEID", multiVals = "first")
genes <- c("ZEB1", "FPGS", "SLC19A1", "GGH", "DHFR", "TYMS", "SHMT1", "SHMT2")

# Frozen collapse: mean of every probeset annotated to the symbol.
# Sensitivity: the single probeset with the highest mean across all 13 lines
# (chosen without using IC50).
gene_mean <- sapply(genes, function(g) {
  idx <- which(sym == g)
  if (!length(idx)) return(rep(NA_real_, ncol(mat)))
  colMeans(mat[idx, , drop = FALSE])
})
gene_best <- sapply(genes, function(g) {
  idx <- which(sym == g)
  if (!length(idx)) return(rep(NA_real_, ncol(mat)))
  means <- rowMeans(mat[idx, , drop = FALSE])
  mat[idx[which.max(means)], ]
})
probe_note <- sapply(genes, function(g) {
  idx <- which(sym == g)
  if (!length(idx)) return("")
  means <- rowMeans(mat[idx, , drop = FALSE])
  paste(rownames(mat)[idx[order(-means)]], sprintf("%.2f", sort(means, decreasing = TRUE)),
        sep = ":", collapse = ";")
})

map_line <- c(
  ALL_SIL = "ALL-SIL", SIL = "ALL-SIL", CEM = "CCRF-CEM", HPB = "HPB-ALL",
  Jurkat = "Jurkat", Molt4 = "Molt-4", P12 = "P12-ICHIKAWA", TALL = "TALL-1",
  Nalm16 = "Nalm16", Nalm30 = "Nalm30", Nalm6 = "Nalm6", REH = "REH",
  RS4_11 = "RS4;11", `697` = "697"
)
sample_key <- sub("_Control$", "", colnames(mat))
cell_line <- unname(map_line[sample_key])
stopifnot(!any(is.na(cell_line)))

expr <- data.frame(cell_line = cell_line, gene_mean, check.names = FALSE)
best <- data.frame(cell_line = cell_line, gene_best, check.names = FALSE)
names(best)[-1] <- paste0(names(best)[-1], "_bestprobe")
write.table(expr, file.path(val, "silva2022_rma_genes.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
write.table(best, file.path(val, "silva2022_rma_bestprobe.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
writeLines(probe_note, file.path(val, "silva2022_probe_means.txt"))

pheno <- read.delim(file.path(val, "silva2022_mtx_ic50.tsv"), stringsAsFactors = FALSE)
gsh <- read.delim(file.path(val, "silva2022_basal_metabolites.tsv"), stringsAsFactors = FALSE)
dat <- merge(merge(expr, pheno, by = "cell_line"), gsh, by = "cell_line")
dat <- merge(dat, best, by = "cell_line")
write.table(dat, file.path(val, "silva2022_zeb1_mtx_joined.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

spearman <- function(x, y) {
  ok <- is.finite(x) & is.finite(y)
  ct <- suppressWarnings(cor.test(x[ok], y[ok], method = "spearman", exact = FALSE))
  c(n = sum(ok), rho = unname(ct$estimate), p = ct$p.value)
}

partial_spearman <- function(x, y, z) {
  ok <- is.finite(x) & is.finite(y) & is.finite(z)
  rx <- rank(x[ok]); ry <- rank(y[ok]); rz <- rank(z[ok])
  xr <- residuals(lm(rx ~ rz))
  yr <- residuals(lm(ry ~ rz))
  ct <- suppressWarnings(cor.test(xr, yr, method = "pearson"))
  c(n = sum(ok), rho = unname(ct$estimate), p = ct$p.value)
}

rows <- list()
add <- function(cohort, yname, xname, stat) {
  rows[[length(rows) + 1]] <<- data.frame(
    cohort = cohort, y = yname, x = xname,
    n = stat["n"], rho = stat["rho"], p = stat["p"],
    stringsAsFactors = FALSE
  )
}

for (cohort in c("T-ALL", "BCP-ALL")) {
  sub <- dat[dat$lineage == cohort, ]
  y <- sub$ic50_96h_nM
  add(cohort, "ic50_96h", "ZEB1", spearman(sub$ZEB1, y))
  add(cohort, "ic50_96h", "ZEB1_bestprobe", spearman(sub$ZEB1_bestprobe, y))
  add(cohort, "ic50_48h", "ZEB1", spearman(sub$ZEB1, sub$ic50_48h_nM))
  for (g in c("FPGS", "SLC19A1", "GGH", "DHFR", "TYMS", "SHMT1", "SHMT2", "Glutathione")) {
    add(cohort, "ic50_96h", g, spearman(sub[[g]], y))
    add(cohort, "ZEB1", g, spearman(sub$ZEB1, sub[[g]]))
    add(cohort, "ic50_96h_partial", paste0("ZEB1|", g),
        partial_spearman(sub$ZEB1, y, sub[[g]]))
  }
  # leave-one-out of the primary association
  for (i in seq_len(nrow(sub))) {
    keep <- sub[-i, ]
    st <- spearman(keep$ZEB1, keep$ic50_96h_nM)
    add(paste0(cohort, "_loo_", sub$cell_line[i]), "ic50_96h", "ZEB1", st)
  }
}

out <- do.call(rbind, rows)
write.table(out, file.path(val, "silva2022_zeb1_mtx_gate.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
print(out[out$cohort %in% c("T-ALL", "BCP-ALL"), ], digits = 3)
cat("probe note\n")
cat(paste(names(probe_note), probe_note, sep = "\t"), sep = "\n")

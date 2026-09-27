# A7-B: TMM logCPM on 15 Dx-malignant patient pseudobulks. Frozen side50 only.
suppressPackageStartupMessages(library(edgeR))

P <- "D:/_bioinformation/ZEB1/total/rebuild_2026/data/deepen_2026"
out <- file.path(P, "a7b_malignant")
fr <- read.delim(file.path(P, "gate2_program/frozen_ZEB_side50.tsv"), stringsAsFactors = FALSE)
cnt <- read.delim(file.path(out, "a7b_patient_counts.tsv.gz"), check.names = FALSE)
genes <- cnt[[1]]
mat <- as.matrix(cnt[, -1, drop = FALSE])
storage.mode(mat) <- "integer"
rownames(mat) <- genes
patients <- colnames(mat)
stopifnot(ncol(mat) == 15)

y <- DGEList(counts = mat)
y <- calcNormFactors(y, method = "TMM")
lcpm <- cpm(y, log = TRUE, prior.count = 1)
write.table(data.frame(symbol = rownames(lcpm), lcpm, check.names = FALSE),
            file.path(out, "a7b_patient_logcpm.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

axis <- c("ZEB1", "ZEB2", "CD1A", "CD34", "LYL1")
stopifnot(all(axis %in% rownames(lcpm)))

fr$symbol_u <- toupper(fr$symbol)
sym_u <- toupper(rownames(lcpm))
fr$idx <- match(fr$symbol_u, sym_u)
fr$recovered <- !is.na(fr$idx)
write.table(fr[, c("symbol", "gene_id", "side", "beta_R", "recovered")],
            file.path(out, "a7b_recovered.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

zrow <- function(v) {
  s <- sd(v)
  if (!is.finite(s) || s == 0) return(rep(NA_real_, length(v)))
  as.numeric((v - mean(v)) / s)
}

rec <- fr[fr$recovered, ]
E <- lcpm[rec$idx, , drop = FALSE]
Ez <- t(scale(t(E), center = TRUE, scale = TRUE))
z1 <- rec$side == "ZEB1-side50"
z2 <- rec$side == "ZEB2-side50"

sc <- data.frame(
  patient = patients,
  ZEB1 = zrow(lcpm["ZEB1", ]),
  ZEB2 = zrow(lcpm["ZEB2", ]),
  CD1A = zrow(lcpm["CD1A", ]),
  CD34 = zrow(lcpm["CD34", ]),
  LYL1 = zrow(lcpm["LYL1", ]),
  ZEB1_side50 = colMeans(Ez[z1, , drop = FALSE], na.rm = TRUE),
  ZEB2_side50 = colMeans(Ez[z2, , drop = FALSE], na.rm = TRUE),
  stringsAsFactors = FALSE
)
sc$B <- sc$ZEB1 - sc$ZEB2
sc$D <- sc$CD1A - (sc$CD34 + sc$LYL1) / 2
lib <- read.delim(file.path(out, "a7b_patient_lib.tsv"), stringsAsFactors = FALSE)
sc <- merge(lib, sc, by = "patient", sort = FALSE)
sc <- sc[match(patients, sc$patient), ]
write.table(sc, file.path(out, "a7b_patient_scores.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

# gene-level: Spearman(E_g, B) vs discovery beta_R
ge <- rec
ge$rho_E_vs_B <- apply(E, 1, function(v) suppressWarnings(cor(v, sc$B, method = "spearman")))
ge$expected <- sign(ge$rho_E_vs_B) == sign(ge$beta_R)
ge$expected[ge$rho_E_vs_B == 0 | is.na(ge$rho_E_vs_B)] <- FALSE
write.table(ge, file.path(out, "a7b_genes.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

rho_z1 <- suppressWarnings(cor(sc$ZEB1_side50, sc$B, method = "spearman"))
rho_z2 <- suppressWarnings(cor(sc$ZEB2_side50, sc$B, method = "spearman"))
rho_gene <- suppressWarnings(cor(ge$rho_E_vs_B, ge$beta_R, method = "spearman"))

# bootstrap patients: re-z, re-score, Spearman vs B
set.seed(248287)
nboot <- 5000
b1 <- b2 <- numeric(nboot)
n <- length(patients)
for (b in seq_len(nboot)) {
  ix <- sample.int(n, n, replace = TRUE)
  Eb <- E[, ix, drop = FALSE]
  Ezb <- t(scale(t(Eb), center = TRUE, scale = TRUE))
  s1 <- colMeans(Ezb[z1, , drop = FALSE], na.rm = TRUE)
  s2 <- colMeans(Ezb[z2, , drop = FALSE], na.rm = TRUE)
  Bb <- zrow(lcpm["ZEB1", ix]) - zrow(lcpm["ZEB2", ix])
  b1[b] <- suppressWarnings(cor(s1, Bb, method = "spearman"))
  b2[b] <- suppressWarnings(cor(s2, Bb, method = "spearman"))
}
ci <- function(x) as.numeric(quantile(x, c(0.025, 0.975), na.rm = TRUE))
c1 <- ci(b1)
c2 <- ci(b2)

sm <- data.frame(
  item = c(
    "n_patients", "n_cells", "recovered", "recovered_Z1", "recovered_Z2",
    "rho_Z1_vs_B", "rho_Z2_vs_B",
    "boot_Z1_lo", "boot_Z1_hi", "boot_Z2_lo", "boot_Z2_hi",
    "gene_expected_frac", "gene_expected_Z1", "gene_expected_Z2",
    "rho_gene_vs_betaR"
  ),
  value = c(
    15, sum(sc$n_cells), nrow(rec), sum(z1), sum(z2),
    rho_z1, rho_z2,
    c1[1], c1[2], c2[1], c2[2],
    mean(ge$expected),
    mean(ge$expected[z1]),
    mean(ge$expected[z2]),
    rho_gene
  ),
  stringsAsFactors = FALSE
)
write.table(sm, file.path(out, "a7b_summary.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
write.table(data.frame(boot_Z1 = b1, boot_Z2 = b2),
            file.path(out, "a7b_bootstrap.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
cat(capture.output(print(sc[, c("patient", "n_cells", "B", "ZEB1_side50", "ZEB2_side50")])), sep = "\n")
cat("\n")
print(sm)
cat("Z1_CI_excludes0", (c1[1] > 0 | c1[2] < 0), "\n")
cat("Z2_CI_excludes0", (c2[1] > 0 | c2[2] < 0), "\n")

# A8-1: frozen residual + side50 on all 79 GSE280250 samples. No subtype yet.
suppressPackageStartupMessages({
  library(edgeR)
  library(limma)
  library(splines)
})

P <- "D:/_bioinformation/ZEB1/data/analysis/total/rebuild_2026/data/deepen_2026"
out <- file.path(P, "a8_adult")
fr <- read.delim(file.path(P, "gate2_program/frozen_ZEB_side50.tsv"), stringsAsFactors = FALSE)
cnt <- read.delim(file.path(out, "a8_counts.tsv.gz"), check.names = FALSE)
ph <- read.delim(file.path(out, "a8_phenotype.tsv"), stringsAsFactors = FALSE)
genes <- cnt[[1]]
mat <- as.matrix(cnt[, -1, drop = FALSE])
storage.mode(mat) <- "integer"
rownames(mat) <- genes
sid <- colnames(mat)
stopifnot(ncol(mat) == 79L, setequal(sid, ph$geo_accession))
ph <- ph[match(sid, ph$geo_accession), ]

y <- DGEList(counts = mat)
y <- calcNormFactors(y, method = "TMM")
lcpm <- cpm(y, log = TRUE, prior.count = 1)
axis <- c(ZEB1 = "ZEB1", ZEB2 = "ZEB2", CD1A = "CD1A", CD34 = "CD34", LYL1 = "LYL1")
stopifnot(all(axis %in% rownames(lcpm)))

zsd <- function(x) {
  s <- sd(x)
  if (!is.finite(s) || s == 0) return(rep(NA_real_, length(x)))
  as.numeric((x - mean(x)) / s)
}

sc <- data.frame(
  geo_accession = sid,
  title = ph$title,
  age = ph$age,
  sex = ph$sex,
  disease = ph$disease,
  ZEB1 = as.numeric(lcpm["ZEB1", ]),
  ZEB2 = as.numeric(lcpm["ZEB2", ]),
  CD1A = as.numeric(lcpm["CD1A", ]),
  CD34 = as.numeric(lcpm["CD34", ]),
  LYL1 = as.numeric(lcpm["LYL1", ]),
  stringsAsFactors = FALSE
)
sc$B <- zsd(sc$ZEB1) - zsd(sc$ZEB2)
sc$D <- zsd(sc$CD1A) - (zsd(sc$CD34) + zsd(sc$LYL1)) / 2
fit <- lm(B ~ ns(D, df = 3), data = sc)
sc$R <- as.numeric(resid(fit))
sc$residual_z <- as.numeric(scale(sc$R))

fr$symbol_u <- toupper(fr$symbol)
fr$idx <- match(fr$symbol_u, toupper(rownames(lcpm)))
fr$recovered <- !is.na(fr$idx)
write.table(fr[, c("symbol", "gene_id", "side", "beta_R", "recovered")],
            file.path(out, "a8_recovered.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
rec <- fr[fr$recovered, ]
E <- lcpm[rec$idx, , drop = FALSE]
Ez <- t(scale(t(E), center = TRUE, scale = TRUE))
z1 <- rec$side == "ZEB1-side50"
z2 <- rec$side == "ZEB2-side50"
sc$ZEB1_side50 <- colMeans(Ez[z1, , drop = FALSE], na.rm = TRUE)
sc$ZEB2_side50 <- colMeans(Ez[z2, , drop = FALSE], na.rm = TRUE)
write.table(sc, file.path(out, "a8_patient_scores.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

rho_z1 <- suppressWarnings(cor(sc$ZEB1_side50, sc$R, method = "spearman"))
rho_z2 <- suppressWarnings(cor(sc$ZEB2_side50, sc$R, method = "spearman"))
p_z1 <- suppressWarnings(cor.test(sc$ZEB1_side50, sc$R, method = "spearman")$p.value)
p_z2 <- suppressWarnings(cor.test(sc$ZEB2_side50, sc$R, method = "spearman")$p.value)

# limma analog of discovery β_R without subtype: E ~ ns(D,3) + residual_z
keep <- filterByExpr(y)
yv <- y[keep, , keep.lib.sizes = FALSE]
yv <- calcNormFactors(yv, method = "TMM")
design <- model.matrix(~ ns(sc$D, df = 3) + sc$residual_z)
v <- voom(yv, design, plot = FALSE)
fitv <- eBayes(lmFit(v, design))
tt <- topTable(fitv, coef = "sc$residual_z", number = Inf, sort.by = "none")
tt$symbol_u <- toupper(rownames(tt))
ge <- rec
ge$beta_adult <- tt$logFC[match(ge$symbol_u, tt$symbol_u)]
ge$expected <- sign(ge$beta_adult) == sign(ge$beta_R)
ge$expected[is.na(ge$beta_adult) | ge$beta_adult == 0] <- FALSE
write.table(ge, file.path(out, "a8_genes.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
ok <- ge[!is.na(ge$beta_adult), ]
rho_g <- suppressWarnings(cor(ok$beta_adult, ok$beta_R, method = "spearman"))
p_g <- suppressWarnings(cor.test(ok$beta_adult, ok$beta_R, method = "spearman")$p.value)

set.seed(280250)
nboot <- 5000
b1 <- b2 <- numeric(nboot)
n <- nrow(sc)
for (b in seq_len(nboot)) {
  ix <- sample.int(n, n, replace = TRUE)
  b1[b] <- suppressWarnings(cor(sc$ZEB1_side50[ix], sc$R[ix], method = "spearman"))
  b2[b] <- suppressWarnings(cor(sc$ZEB2_side50[ix], sc$R[ix], method = "spearman"))
}
ci <- function(x) as.numeric(quantile(x, c(0.025, 0.975), na.rm = TRUE))
c1 <- ci(b1)
c2 <- ci(b2)

sm <- data.frame(
  item = c(
    "n_samples", "age_min", "age_max", "n_age_lt18",
    "recovered", "recovered_Z1", "recovered_Z2",
    "rho_Z1_vs_R", "p_Z1_vs_R", "boot_Z1_lo", "boot_Z1_hi",
    "rho_Z2_vs_R", "p_Z2_vs_R", "boot_Z2_lo", "boot_Z2_hi",
    "gene_n_with_beta", "gene_expected_frac", "gene_expected_Z1", "gene_expected_Z2",
    "rho_betaAdult_vs_betaR", "p_gene"
  ),
  value = c(
    79, min(sc$age), max(sc$age), sum(sc$age < 18),
    nrow(rec), sum(z1), sum(z2),
    rho_z1, p_z1, c1[1], c1[2],
    rho_z2, p_z2, c2[1], c2[2],
    nrow(ok), mean(ok$expected),
    mean(ok$expected[ok$side == "ZEB1-side50"]),
    mean(ok$expected[ok$side == "ZEB2-side50"]),
    rho_g, p_g
  ),
  stringsAsFactors = FALSE
)
write.table(sm, file.path(out, "a8_1_summary.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
print(sm)
cat("Z1_dir", rho_z1 > 0, "Z2_dir", rho_z2 < 0, "\n")

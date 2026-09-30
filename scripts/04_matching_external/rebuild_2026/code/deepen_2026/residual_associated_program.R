# Residual-associated transcriptome after developmental spline and subtype.
# Expression_g ~ ns(dev, df=3) + subtype + residual_within.
# The residual is the frozen development-adjusted ZEB balance, not raw high/low
# balance. ZEB1 and ZEB2 are reported but excluded from the frozen program.
suppressPackageStartupMessages({
  library(data.table)
  library(edgeR)
  library(limma)
  library(splines)
  library(msigdbr)
})

argv <- grep("^--file=", commandArgs(FALSE), value = TRUE)
script <- normalizePath(sub("^--file=", "", argv), winslash = "/")
rebuild <- normalizePath(file.path(dirname(script), "../.."), winslash = "/")
root <- normalizePath(file.path(rebuild, "../.."), winslash = "/")
out <- file.path(rebuild, "data/deepen_2026/residual_program")
dir.create(out, recursive = TRUE, showWarnings = FALSE)

pt <- fread(file.path(root, "total/data/validation/zeb_developmental_residual/polonen_patient_residual.tsv"),
            encoding = "UTF-8")
stopifnot(nrow(pt) == 1309L, uniqueN(pt$sample_id) == 1309L, uniqueN(pt$subtype) == 17L)
stopifnot(all(is.finite(pt$dev)), all(is.finite(pt$residual_within)))
pt[, subtype := factor(subtype)]

cnt <- fread(file.path(root, "data/TALL_X01_counts.tsv"), nrows = 0)
ids <- pt$sample_id
stopifnot(all(ids %in% names(cnt)))
z <- fread(file.path(root, "data/TALL_X01_counts.tsv"), select = c(names(cnt)[1], ids))
gene_id <- z[[1]]
mat <- as.matrix(z[, ..ids])
storage.mode(mat) <- "integer"
rownames(mat) <- gene_id
e <- new.env(parent = emptyenv())
load(file.path(root, "data/polonen_syn54032669/Data_1309Samples.RData"), envir = e)
symbols <- rownames(e$gexp)
stopifnot(length(symbols) == nrow(mat), all(nzchar(symbols)))
rm(e, z, cnt); gc()

# Design uses the already frozen residual; do not re-estimate it here.
dev_ns <- ns(pt$dev, df = 3)
colnames(dev_ns) <- paste0("ns", seq_len(ncol(dev_ns)))
design <- model.matrix(~ dev_ns + subtype + residual_within, data = pt)
stopifnot(nrow(design) == ncol(mat))

y <- DGEList(counts = mat)
keep <- filterByExpr(y, design = design)
y <- y[keep, , keep.lib.sizes = FALSE]
sym <- symbols[keep]
y <- calcNormFactors(y, method = "TMM")
v <- voom(y, design, plot = FALSE)
fit <- eBayes(lmFit(v, design), robust = TRUE)
t <- topTable(fit, coef = "residual_within", number = Inf, sort.by = "none")
t$gene_id <- rownames(t)
t$symbol <- sym
res <- as.data.table(t)[, .(gene_id, symbol, logFC, AveExpr, t, P.Value, adj.P.Val)]
setorder(res, P.Value)
fwrite(res, file.path(out, "residual_associated_genes.tsv"), sep = "\t")

exclude <- toupper(c("ZEB1", "ZEB2", "ZEB1-AS1"))
prog <- res[!(toupper(symbol) %in% exclude)]
sig <- prog[adj.P.Val < 0.05]
zeb1_side <- sig[t > 0][order(-t)]
zeb2_side <- sig[t < 0][order(t)]
n_freeze <- 100L
freeze <- rbind(
  zeb1_side[seq_len(min(n_freeze, .N))][, `:=`(side = "ZEB1-side", rank = .I)],
  zeb2_side[seq_len(min(n_freeze, .N))][, `:=`(side = "ZEB2-side", rank = .I)]
)
fwrite(freeze, file.path(out, "frozen_residual_program_top100.tsv"), sep = "\t")

gs <- msigdbr(species = "Homo sapiens", collection = "H")
sets <- split(gs$gene_symbol, gs$gs_name)
indices <- lapply(sets, function(s) which(sym %in% unique(s)))
indices <- indices[lengths(indices) >= 10L & lengths(indices) <= 500L]
cam <- camera(v, index = indices, design = design, contrast = "residual_within")
cam$pathway <- rownames(cam)
cam <- as.data.table(cam)[, .(pathway, NGenes, Direction, PValue, FDR)]
setorder(cam, FDR, PValue)
fwrite(cam, file.path(out, "residual_associated_hallmark_camera.tsv"), sep = "\t")

summ <- data.table(
  n_patients = nrow(pt),
  n_genes_tested = nrow(res),
  n_fdr05 = nrow(res[adj.P.Val < 0.05]),
  n_fdr05_excluding_ZEB = nrow(sig),
  n_ZEB1_side = nrow(zeb1_side),
  n_ZEB2_side = nrow(zeb2_side),
  n_frozen_ZEB1_side = sum(freeze$side == "ZEB1-side"),
  n_frozen_ZEB2_side = sum(freeze$side == "ZEB2-side"),
  n_hallmark_fdr05 = nrow(cam[FDR < 0.05]),
  seed_note = "No additional resampling; inference is limma/voom empirical Bayes on the frozen residual."
)
fwrite(summ, file.path(out, "program_summary.tsv"), sep = "\t")
writeLines(c(
  "Model: voom/limma Expression ~ ns(dev,3) + 17-level subtype + residual_within.",
  "residual_within is the frozen Pölönen development-adjusted ZEB balance.",
  "This is not high-versus-low raw balance, and it is not a causal ZEB knockout.",
  "ZEB1 and ZEB2 are omitted from the frozen side lists because they define the residual.",
  "Gene sets: MSigDB Hallmark via msigdbr; CAMERA competitive test.",
  "Freeze: up to 100 genes per side among FDR<0.05, ranked by |t|.",
  paste("edgeR", as.character(packageVersion("edgeR"))),
  paste("limma", as.character(packageVersion("limma"))),
  paste("msigdbr", as.character(packageVersion("msigdbr")))
), file.path(out, "METHODS.txt"))
writeLines(capture.output(sessionInfo()), file.path(out, "R_sessionInfo.txt"))
cat("Residual program complete\n")
print(summ)

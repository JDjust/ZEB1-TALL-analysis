# A7-A scores on frozen map. No reselection. No new top50.
suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(edgeR)
})

dir <- "/data-b/liangfuhua/zeb1_a7"
out <- file.path(dir, "scores")
dir.create(out, recursive = TRUE, showWarnings = FALSE)

fr <- read.delim(file.path(dir, "frozen_ZEB_side50.tsv"), stringsAsFactors = FALSE)
map <- read.delim(file.path(dir, "a7_label_map.tsv"), stringsAsFactors = FALSE)
l2_to_broad <- setNames(map$broad, map$cluster_anno_l2)

official <- c("HSC/MPP/HSPC", "lymphoid progenitor/CLP", "T lineage", "Myeloid")
myeloid_broads <- c("myeloid progenitor/GMP", "monocyte/macrophage", "granulocytic")

cat("load object\n"); flush.console()
obj <- readRDS(gzcon(gzfile(file.path(dir, "GSE253355_Normal_Bone_Marrow_Atlas_Seurat_SB_v2.rds.gz"), "rb")))
md <- obj[[]]
cnt <- tryCatch(
  GetAssayData(obj, assay = "RNA", layer = "counts"),
  error = function(e) LayerData(obj[["RNA"]], layer = "counts")
)
genes <- rownames(cnt)
cat("counts", paste(dim(cnt), collapse = "x"), "\n"); flush.console()
rm(obj); gc()

md$l2 <- as.character(md$cluster_anno_l2)
md$broad <- unname(l2_to_broad[md$l2])
md$donor <- as.character(md$orig.ident)
md$official <- md$broad
md$official[md$broad %in% myeloid_broads] <- "Myeloid"

# match frozen genes: symbol then ENSEMBL
sym_u <- toupper(genes)
ens <- sub("\\.\\d+$", "", genes)
fr$symbol_u <- toupper(fr$symbol)
fr$ens <- sub("\\.\\d+$", "", fr$gene_id)
hit <- match(fr$symbol_u, sym_u)
hit2 <- match(fr$ens, ens)
hit[is.na(hit)] <- hit2[is.na(hit)]
fr$idx <- hit
fr$recovered <- !is.na(fr$idx)
write.table(fr[, c("symbol", "gene_id", "side", "beta_R", "recovered")],
            file.path(out, "a7a_recovered.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
cat("recovered", sum(fr$recovered), "/100  ZEB1", sum(fr$recovered & fr$side == "ZEB1-side50"),
    "ZEB2", sum(fr$recovered & fr$side == "ZEB2-side50"), "\n"); flush.console()

# donor x official unit cell counts
cells <- colnames(cnt)
stopifnot(identical(cells, rownames(md)))
units <- unique(md[, c("donor", "official")])
units <- units[!is.na(units$official) & units$official %in% official, ]
n_tab <- aggregate(rep(1L, nrow(md)), by = list(donor = md$donor, official = md$official), FUN = length)
names(n_tab)[3] <- "n_cells"
n_tab <- n_tab[n_tab$official %in% official, ]
write.table(n_tab, file.path(out, "a7a_n_cells.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

make_pb <- function(min_n) {
  keep <- n_tab[n_tab$n_cells >= min_n, ]
  mat <- matrix(0, nrow = nrow(cnt), ncol = nrow(keep))
  colnames(mat) <- paste(keep$donor, keep$official, sep = "||")
  rownames(mat) <- genes
  for (i in seq_len(nrow(keep))) {
    ix <- which(md$donor == keep$donor[i] & md$official == keep$official[i])
    mat[, i] <- Matrix::rowSums(cnt[, ix, drop = FALSE])
  }
  list(keep = keep, mat = mat)
}

run_scores <- function(min_n, tag) {
  pb <- make_pb(min_n)
  keep <- pb$keep
  y <- DGEList(counts = pb$mat)
  y <- calcNormFactors(y, method = "TMM")
  lcpm <- cpm(y, log = TRUE, prior.count = 1)
  rec <- fr[fr$recovered, ]
  E <- lcpm[rec$idx, , drop = FALSE]
  rownames(E) <- rec$symbol
  # gene z across units
  Ez <- t(scale(t(E), center = TRUE, scale = TRUE))
  z1 <- rec$side == "ZEB1-side50"
  z2 <- rec$side == "ZEB2-side50"
  sc <- data.frame(
    donor = keep$donor,
    lineage = keep$official,
    n_cells = keep$n_cells,
    ZEB1_side50 = colMeans(Ez[z1, , drop = FALSE], na.rm = TRUE),
    ZEB2_side50 = colMeans(Ez[z2, , drop = FALSE], na.rm = TRUE),
    stringsAsFactors = FALSE
  )
  write.table(sc, file.path(out, paste0("a7a_unit_scores_", tag, ".tsv")),
              sep = "\t", quote = FALSE, row.names = FALSE)
  write.table(cbind(sc, t(E)), file.path(out, paste0("a7a_unit_logcpm_", tag, ".tsv")),
              sep = "\t", quote = FALSE, row.names = FALSE)

  # donor-level paired deltas
  dons <- sort(unique(sc$donor))
  drows <- list()
  grows_h <- matrix(NA_real_, nrow = sum(rec$recovered), ncol = 0)
  grows_m <- matrix(NA_real_, nrow = sum(rec$recovered), ncol = 0)
  rownames(grows_h) <- rec$symbol
  rownames(grows_m) <- rec$symbol
  for (d in dons) {
    sl <- sc[sc$donor == d, ]
    getv <- function(lin, col) {
      x <- sl[sl$lineage == lin, col]
      if (length(x)) as.numeric(x[1]) else NA_real_
    }
    recd <- data.frame(
      donor = d,
      n_HSPC = getv("HSC/MPP/HSPC", "n_cells"),
      n_CLP = getv("lymphoid progenitor/CLP", "n_cells"),
      n_T = getv("T lineage", "n_cells"),
      n_Myeloid = getv("Myeloid", "n_cells"),
      HSPC_Z1 = getv("HSC/MPP/HSPC", "ZEB1_side50"),
      CLP_Z1 = getv("lymphoid progenitor/CLP", "ZEB1_side50"),
      T_Z1 = getv("T lineage", "ZEB1_side50"),
      My_Z1 = getv("Myeloid", "ZEB1_side50"),
      HSPC_Z2 = getv("HSC/MPP/HSPC", "ZEB2_side50"),
      CLP_Z2 = getv("lymphoid progenitor/CLP", "ZEB2_side50"),
      T_Z2 = getv("T lineage", "ZEB2_side50"),
      My_Z2 = getv("Myeloid", "ZEB2_side50"),
      stringsAsFactors = FALSE
    )
    recd$d_HSPC_minus_T_Z1 <- recd$HSPC_Z1 - recd$T_Z1
    recd$d_CLP_minus_T_Z1 <- recd$CLP_Z1 - recd$T_Z1
    recd$d_My_minus_T_Z1 <- recd$My_Z1 - recd$T_Z1
    recd$d_HSPC_minus_T_Z2 <- recd$HSPC_Z2 - recd$T_Z2
    recd$d_CLP_minus_T_Z2 <- recd$CLP_Z2 - recd$T_Z2
    recd$d_My_minus_T_Z2 <- recd$My_Z2 - recd$T_Z2
    drows[[d]] <- recd
    # gene deltas on logCPM, donors with both
    cols <- colnames(E)
    hcol <- paste(d, "HSC/MPP/HSPC", sep = "||")
    tcol <- paste(d, "T lineage", sep = "||")
    mcol <- paste(d, "Myeloid", sep = "||")
    if (hcol %in% cols && tcol %in% cols) {
      grows_h <- cbind(grows_h, E[, hcol] - E[, tcol])
    }
    if (mcol %in% cols && tcol %in% cols) {
      grows_m <- cbind(grows_m, E[, mcol] - E[, tcol])
    }
  }
  don <- do.call(rbind, drows)
  write.table(don, file.path(out, paste0("a7a_donor_", tag, ".tsv")),
              sep = "\t", quote = FALSE, row.names = FALSE)

  ge <- rec
  ge$delta_HSPC_minus_T <- if (ncol(grows_h)) rowMeans(grows_h) else NA_real_
  ge$delta_Myeloid_minus_T <- if (ncol(grows_m)) rowMeans(grows_m) else NA_real_
  ge$neg_beta <- -ge$beta_R
  ge$exp_HSPC <- sign(ge$delta_HSPC_minus_T) == sign(ge$neg_beta)
  ge$exp_My <- sign(ge$delta_Myeloid_minus_T) == sign(ge$neg_beta)
  ge$exp_HSPC[ge$delta_HSPC_minus_T == 0] <- FALSE
  ge$exp_My[ge$delta_Myeloid_minus_T == 0] <- FALSE
  write.table(ge, file.path(out, paste0("a7a_genes_", tag, ".tsv")),
              sep = "\t", quote = FALSE, row.names = FALSE)

  w <- function(x) {
    x <- x[is.finite(x)]
    if (length(x) < 4) return(c(n = length(x), median = median(x), p = NA_real_))
    p <- tryCatch(wilcox.test(x, alternative = "two.sided", exact = TRUE)$p.value, error = function(e) NA_real_)
    c(n = length(x), median = median(x), p = p)
  }
  rho <- function(a, b) {
    ok <- is.finite(a) & is.finite(b)
    if (sum(ok) < 5) return(c(rho = NA_real_, p = NA_real_))
    s <- suppressWarnings(cor.test(a[ok], b[ok], method = "spearman", exact = FALSE))
    c(rho = unname(s$estimate), p = s$p.value)
  }
  frac <- function(x) paste0(sum(x, na.rm = TRUE), "/", sum(is.finite(x)))
  ntrue <- function(x) paste0(sum(x > 0, na.rm = TRUE), "/", sum(is.finite(x)))

  rh <- rho(ge$delta_HSPC_minus_T, ge$neg_beta)
  rm_ <- rho(ge$delta_Myeloid_minus_T, ge$neg_beta)
  rh1 <- rho(ge$delta_HSPC_minus_T[z1], ge$neg_beta[z1])
  rh2 <- rho(ge$delta_HSPC_minus_T[z2], ge$neg_beta[z2])
  rm1 <- rho(ge$delta_Myeloid_minus_T[z1], ge$neg_beta[z1])
  rm2 <- rho(ge$delta_Myeloid_minus_T[z2], ge$neg_beta[z2])

  wh1 <- w(don$d_HSPC_minus_T_Z1)
  wh2 <- w(don$d_HSPC_minus_T_Z2)
  wm1 <- w(don$d_My_minus_T_Z1)
  wm2 <- w(don$d_My_minus_T_Z2)
  wc1 <- w(don$d_CLP_minus_T_Z1)
  wc2 <- w(don$d_CLP_minus_T_Z2)

  med <- function(x) median(x, na.rm = TRUE)
  sm <- data.frame(
    item = c(
      "min_n", "recovered_Z1", "recovered_Z2",
      "n_donor_HSPC_and_T", "n_donor_CLP_and_T", "n_donor_My_and_T",
      "HSPC_Z1_median", "CLP_Z1_median", "T_Z1_median", "My_Z1_median",
      "HSPC_Z2_median", "CLP_Z2_median", "T_Z2_median", "My_Z2_median",
      "donors_HSPC_gt_T_Z1", "donors_HSPC_gt_T_Z2",
      "donors_My_gt_T_Z1", "donors_My_gt_T_Z2",
      "donors_CLP_gt_T_Z1", "donors_CLP_gt_T_Z2",
      "median_d_HSPC_T_Z1", "wilcox_HSPC_T_Z1_p",
      "median_d_HSPC_T_Z2", "wilcox_HSPC_T_Z2_p",
      "median_d_My_T_Z1", "wilcox_My_T_Z1_p",
      "median_d_My_T_Z2", "wilcox_My_T_Z2_p",
      "median_d_CLP_T_Z1", "wilcox_CLP_T_Z1_p",
      "median_d_CLP_T_Z2", "wilcox_CLP_T_Z2_p",
      "gene_n_HSPC", "gene_exp_HSPC", "rho_HSPC_vs_negbeta", "rho_HSPC_p",
      "gene_exp_HSPC_Z1", "gene_exp_HSPC_Z2",
      "gene_n_My", "gene_exp_My", "rho_My_vs_negbeta", "rho_My_p",
      "gene_exp_My_Z1", "gene_exp_My_Z2"
    ),
    value = c(
      min_n,
      paste0(sum(z1), "/50"), paste0(sum(z2), "/50"),
      sum(is.finite(don$d_HSPC_minus_T_Z1)),
      sum(is.finite(don$d_CLP_minus_T_Z1)),
      sum(is.finite(don$d_My_minus_T_Z1)),
      med(don$HSPC_Z1), med(don$CLP_Z1), med(don$T_Z1), med(don$My_Z1),
      med(don$HSPC_Z2), med(don$CLP_Z2), med(don$T_Z2), med(don$My_Z2),
      ntrue(don$d_HSPC_minus_T_Z1), ntrue(don$d_HSPC_minus_T_Z2),
      ntrue(don$d_My_minus_T_Z1), ntrue(don$d_My_minus_T_Z2),
      ntrue(don$d_CLP_minus_T_Z1), ntrue(don$d_CLP_minus_T_Z2),
      unname(wh1["median"]), unname(wh1["p"]),
      unname(wh2["median"]), unname(wh2["p"]),
      unname(wm1["median"]), unname(wm1["p"]),
      unname(wm2["median"]), unname(wm2["p"]),
      unname(wc1["median"]), unname(wc1["p"]),
      unname(wc2["median"]), unname(wc2["p"]),
      ncol(grows_h), mean(ge$exp_HSPC, na.rm = TRUE), unname(rh["rho"]), unname(rh["p"]),
      mean(ge$exp_HSPC[z1], na.rm = TRUE), mean(ge$exp_HSPC[z2], na.rm = TRUE),
      ncol(grows_m), mean(ge$exp_My, na.rm = TRUE), unname(rm_["rho"]), unname(rm_["p"]),
      mean(ge$exp_My[z1], na.rm = TRUE), mean(ge$exp_My[z2], na.rm = TRUE)
    ),
    stringsAsFactors = FALSE
  )
  write.table(sm, file.path(out, paste0("a7a_summary_", tag, ".tsv")),
              sep = "\t", quote = FALSE, row.names = FALSE)
  print(sm)
  invisible(list(don = don, ge = ge, sc = sc))
}

cat("primary min20\n"); flush.console()
run_scores(20L, "min20")
cat("sensitivity min10\n"); flush.console()
run_scores(10L, "min10")
cat("done\n")

# Developmental residual of ZEB balance.
# A gene enters the frozen score only if its early-versus-DN3/DP direction
# agrees in at least 3 of the 4 normal thymus datasets. ZEB1, ZEB2, LMO2,
# and T-ALL driver genes are not candidates.
# developmental_score = mean(z of DP-high genes) - mean(z of early-high genes).
# The expected-balance curve is fit on normal stage units only, then applied
# to the Pölönen 1309. A second model asks whether subtype still explains
# balance after the developmental score.

suppressPackageStartupMessages({
  library(data.table)
  library(edgeR)
  library(ggplot2)
  library(splines)
})

out_dir <- "D:/_bioinformation/ZEB1/total/data/validation/zeb_developmental_residual"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
fig_dir <- file.path(out_dir, "figures")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

early_cand <- c("CD34", "LYL1", "MEF2C", "SPI1", "KIT")
dp_cand <- c("PTCRA", "RAG1", "RAG2", "CD1A", "CD1B", "LEF1", "TRAC")
candidates <- c(early_cand, dp_cand)

zsd <- function(x) {
  x <- as.numeric(x)
  s <- stats::sd(x, na.rm = TRUE)
  if (!is.finite(s) || s == 0) return(rep(NA_real_, length(x)))
  (x - mean(x, na.rm = TRUE)) / s
}

num <- function(x) suppressWarnings(as.numeric(x))

# ---- stage matrices: genes in columns, one row per stage unit ----
read_195 <- function() {
  a <- fread("D:/_bioinformation/ZEB1/modules/module2/tables/M2_r2_GSE195812_stage_gene_means.tsv")
  order_195 <- c("DN1", "DN2", "DN3", "ISP", "DP_CD3neg", "DP_CD3pos", "CD4SP", "CD8SP")
  a <- a[facs_stage %in% order_195]
  a$stage <- a$facs_stage
  a$order <- match(a$stage, order_195)
  a$window <- ifelse(a$stage %in% c("DN1", "DN2"), "early",
                     ifelse(a$stage %in% c("DN3", "ISP", "DP_CD3neg"), "mid", "late"))
  genes <- sub("^g_", "", grep("^g_", names(a), value = TRUE))
  for (g in genes) a[[g]] <- num(a[[paste0("g_", g)]])
  a$dataset <- "GSE195812"
  a$unit <- a$stage
  a[, c("dataset", "unit", "stage", "order", "window", genes), with = FALSE]
}

read_206 <- function() {
  # donor x fraction, so direction is not a single pooled library
  c <- fread("D:/_bioinformation/ZEB1/modules/module2/tables/M2_r2_GSE206710_sample_means.tsv")
  order_206 <- c("DN_early", "DP_immature", "DP_mature")
  c <- c[stage %in% order_206]
  c$order <- match(c$stage, order_206)
  c$window <- ifelse(c$stage == "DN_early", "early",
                     ifelse(c$stage == "DP_immature", "mid", "late"))
  genes <- sub("^g_", "", grep("^g_", names(c), value = TRUE))
  for (g in genes) c[[g]] <- num(c[[paste0("g_", g)]])
  c$dataset <- "GSE206710"
  # Several libraries from donor D1 occupy the same stage. Use the unique
  # sample-stage identifier; donor remains available as a nesting variable.
  c$unit <- paste(c$sample, c$stage, sep = "|")
  c[, c("dataset", "unit", "stage", "order", "window", "donor", genes), with = FALSE]
}

read_park <- function() {
  p <- fread("D:/_bioinformation/ZEB1/total/data/SourceData_ParkHTA_donor_celltype.tsv")
  keep <- c(
    "double negative thymocyte",
    "double-positive, alpha-beta thymocyte",
    "CD4-positive, alpha-beta T cell",
    "CD8-positive, alpha-beta T cell"
  )
  lab <- c(
    "double negative thymocyte" = "DN",
    "double-positive, alpha-beta thymocyte" = "DP",
    "CD4-positive, alpha-beta T cell" = "CD4 SP",
    "CD8-positive, alpha-beta T cell" = "CD8 SP"
  )
  p <- p[celltype %in% keep & n_cells >= 20]
  p$stage <- unname(lab[p$celltype])
  p$order <- match(p$stage, c("DN", "DP", "CD4 SP", "CD8 SP"))
  p$window <- ifelse(p$stage == "DN", "early", ifelse(p$stage == "DP", "mid", "late"))
  p$dataset <- "Park/HTA"
  p$unit <- paste(p$donor, p$stage, sep = "|")
  p
}

read_142 <- function() {
  b <- fread("D:/_bioinformation/ZEB1/total/data/SourceData_Fig2A_thymus_trajectory.tsv")
  stages <- c("cd34", "isp", "ec", "lc", "sp4", "sp8")
  labels <- c("CD34+", "ISP", "early cortical", "late cortical", "CD4 SP", "CD8 SP")
  windows <- c("early", "other", "mid", "mid", "late", "late")
  wide <- as.data.table(t(as.matrix(b[, ..stages])))
  setnames(wide, b$gene)
  for (g in names(wide)) wide[[g]] <- num(wide[[g]])
  wide$stage <- labels
  wide$order <- seq_len(nrow(wide))
  wide$window <- windows
  wide$dataset <- "GSE142522"
  wide$unit <- wide$stage
  wide
}

g195 <- read_195()
g206 <- read_206()
gpark <- read_park()
g142 <- read_142()

direction_one <- function(df, gene) {
  if (!gene %in% names(df)) {
    return(data.table(gene = gene, dataset = df$dataset[1], n_early = 0L, n_mid = 0L,
                       early_mean = NA_real_, mid_mean = NA_real_, mid_minus_early = NA_real_))
  }
  e <- df[window == "early", num(get(gene))]
  m <- df[window == "mid", num(get(gene))]
  e <- e[is.finite(e)]
  m <- m[is.finite(m)]
  data.table(
    gene = gene,
    dataset = df$dataset[1],
    n_early = length(e),
    n_mid = length(m),
    early_mean = if (length(e)) mean(e) else NA_real_,
    mid_mean = if (length(m)) mean(m) else NA_real_,
    mid_minus_early = if (length(e) && length(m)) mean(m) - mean(e) else NA_real_
  )
}

dirs <- rbindlist(lapply(candidates, function(g) {
  rbindlist(list(direction_one(g195, g), direction_one(g206, g),
                 direction_one(gpark, g), direction_one(g142, g)))
}))
dirs$role <- ifelse(dirs$gene %in% early_cand, "early_high", "DP_high")
dirs$agrees <- ifelse(
  !is.finite(dirs$mid_minus_early), NA,
  ifelse(dirs$role == "early_high", dirs$mid_minus_early < 0, dirs$mid_minus_early > 0)
)
summ <- dirs[, .(
  n_measured = sum(is.finite(mid_minus_early)),
  n_agree = sum(agrees %in% TRUE),
  datasets_agree = paste(dataset[agrees %in% TRUE], collapse = ";"),
  datasets_disagree = paste(dataset[agrees %in% FALSE], collapse = ";"),
  datasets_missing = paste(dataset[!is.finite(mid_minus_early)], collapse = ";")
), by = .(gene, role)]
summ$frozen <- summ$n_agree >= 3 & summ$n_measured >= 3 & summ$n_agree == summ$n_measured
# a gene measured in 4 with 3 agreeing also passes; disagreement in a measured set does not
summ$frozen <- summ$n_agree >= 3
fwrite(dirs, file.path(out_dir, "gene_direction_by_dataset.tsv"), sep = "\t")
fwrite(summ, file.path(out_dir, "gene_freeze_decision.tsv"), sep = "\t")

frozen <- summ[frozen == TRUE]
if (nrow(frozen) < 2 || !any(frozen$role == "early_high") || !any(frozen$role == "DP_high")) {
  stop("Frozen score needs at least one early-high and one DP-high gene with support in 3 datasets")
}
early_genes <- frozen[role == "early_high", gene]
dp_genes <- frozen[role == "DP_high", gene]
message("FROZEN early: ", paste(early_genes, collapse = ", "))
message("FROZEN DP: ", paste(dp_genes, collapse = ", "))

score_units <- function(df, dataset_label) {
  need <- c(early_genes, dp_genes, "ZEB1", "ZEB2")
  miss <- setdiff(need, names(df))
  if (length(miss)) {
    message(dataset_label, " lacks ", paste(miss, collapse = ","), "; skipped for the curve")
    return(NULL)
  }
  for (g in need) df[[g]] <- num(df[[g]])
  if (any(!is.finite(unlist(df[, ..need])))) {
    # drop units with a missing frozen gene or ZEB
    ok <- stats::complete.cases(df[, ..need])
    df <- df[ok]
  }
  if (nrow(df) < 4) return(NULL)
  for (g in need) df[[paste0("z_", g)]] <- zsd(df[[g]])
  df$dev <- rowMeans(df[, paste0("z_", dp_genes), with = FALSE]) -
    rowMeans(df[, paste0("z_", early_genes), with = FALSE])
  df$balance <- df$z_ZEB1 - df$z_ZEB2
  df$dataset <- dataset_label
  df[, .(dataset, unit, stage, order, window, dev, balance)]
}

normal_pts <- rbindlist(list(
  score_units(g195, "GSE195812"),
  score_units(g206, "GSE206710"),
  score_units(gpark, "Park/HTA"),
  score_units(g142, "GSE142522")
))
fwrite(normal_pts, file.path(out_dir, "normal_stage_units_scores.tsv"), sep = "\t")

# Natural spline on the pooled normal units. df=3 matches a rise-then-fall shape
# without fitting one parameter per point.
normal_pts <- normal_pts[is.finite(dev) & is.finite(balance)]
stopifnot(nrow(normal_pts) >= 8)
fit_normal <- lm(balance ~ ns(dev, df = 3), data = normal_pts)
normal_pts$expected_balance <- predict(fit_normal, normal_pts)
normal_pts$residual <- normal_pts$balance - normal_pts$expected_balance
r2_normal <- summary(fit_normal)$r.squared

# ---- Pölönen: same TMM logCPM definition as the frozen balance ----
message("TMM logCPM for frozen genes")
e <- new.env(parent = emptyenv())
load("D:/_bioinformation/ZEB1/data/polonen_syn54032669/Data_1309Samples.RData", envir = e)
annot <- e$annot
ids <- rownames(annot)
gexp <- e$gexp
sub_f <- e$M3.subtype.factor
counts <- fread("D:/_bioinformation/ZEB1/data/TALL_X01_counts.tsv")
gene_col <- names(counts)[1]
setnames(counts, gene_col, "ensembl")
stopifnot(all(ids %in% names(counts)))
mat <- as.matrix(counts[, ..ids])
storage.mode(mat) <- "double"
rownames(mat) <- counts$ensembl
sym <- rownames(gexp)
stopifnot(nrow(mat) == length(sym))
want <- c(early_genes, dp_genes, "ZEB1", "ZEB2")
hit <- sapply(want, function(g) which(sym == g))
stopifnot(all(lengths(hit) == 1))
y <- DGEList(counts = mat)
y <- calcNormFactors(y, method = "TMM")
lcpm <- cpm(y, log = TRUE, prior.count = 1)
# sapply puts one gene per column and one sample per row
expr <- as.data.table(sapply(want, function(g) lcpm[hit[[g]], ]))
expr$sample_id <- ids
for (g in want) expr[[paste0("z_", g)]] <- zsd(expr[[g]])
expr$dev <- rowMeans(expr[, paste0("z_", dp_genes), with = FALSE]) -
  rowMeans(expr[, paste0("z_", early_genes), with = FALSE])
expr$balance <- expr$z_ZEB1 - expr$z_ZEB2

pt <- fread("D:/_bioinformation/ZEB1/total/data/validation/polonen_round1/patient_level_frozen_balance.tsv")
# confirm this run reproduces the frozen balance
mrg <- merge(expr[, .(sample_id, balance)], pt[, .(sample_id, balance_frozen = balance, subtype)], by = "sample_id")
stopifnot(cor(mrg$balance, mrg$balance_frozen) > 0.999)

expr$expected_balance <- predict(fit_normal, newdata = data.frame(dev = expr$dev))
expr$residual_normal <- expr$balance - expr$expected_balance
expr$dev_outside_normal_range <- expr$dev < min(normal_pts$dev) | expr$dev > max(normal_pts$dev)

clin <- data.table(
  sample_id = ids,
  subtype = as.character(sub_f$Reviewed.subtype)
)
expr <- merge(expr, clin, by = "sample_id", sort = FALSE)
expr <- expr[!is.na(subtype) & subtype != ""]

# within-cohort spline, used only for the variance split
fit_dev <- lm(balance ~ ns(dev, df = 3), data = expr)
fit_sub <- lm(balance ~ subtype, data = expr)
fit_both <- lm(balance ~ ns(dev, df = 3) + subtype, data = expr)
expr$expected_within <- predict(fit_dev)
expr$residual_within <- expr$balance - expr$expected_within

r2_tab <- data.table(
  model = c(
    "normal spline: balance ~ ns(dev,3), stage units",
    "Polonen: balance ~ ns(dev,3)",
    "Polonen: balance ~ subtype",
    "Polonen: balance ~ ns(dev,3) + subtype"
  ),
  n = c(nrow(normal_pts), nrow(expr), nrow(expr), nrow(expr)),
  r_squared = c(
    summary(fit_normal)$r.squared,
    summary(fit_dev)$r.squared,
    summary(fit_sub)$r.squared,
    summary(fit_both)$r.squared
  )
)
r2_tab$incremental_r2_subtype_after_dev <- c(NA, NA, NA, r2_tab$r_squared[4] - r2_tab$r_squared[2])
anova_add <- anova(fit_dev, fit_both)
r2_tab$anova_p_subtype_after_dev <- c(NA, NA, NA, anova_add$`Pr(>F)`[2])
fwrite(r2_tab, file.path(out_dir, "variance_partition.tsv"), sep = "\t")

# subtype summary on the normal-anchored residual
sub_sum <- expr[, .(
  n = .N,
  dev_median = median(dev),
  balance_median = median(balance),
  residual_normal_median = median(residual_normal),
  residual_within_median = median(residual_within),
  n_outside_normal_dev_range = sum(dev_outside_normal_range)
), by = subtype]
# one-sample Wilcoxon of normal-anchored residual against 0
wil <- expr[, {
  p <- if (.N >= 5 && length(unique(residual_normal)) > 1) {
    stats::wilcox.test(residual_normal, mu = 0)$p.value
  } else NA_real_
  .(wilcox_p_residual_normal = p)
}, by = subtype]
sub_sum <- merge(sub_sum, wil, by = "subtype")
sub_sum$wilcox_fdr <- p.adjust(sub_sum$wilcox_p_residual_normal, method = "BH")
sub_sum$call <- ifelse(
  is.na(sub_sum$wilcox_fdr), "too few",
  ifelse(sub_sum$wilcox_fdr >= 0.05, "residual compatible with 0",
         ifelse(sub_sum$residual_normal_median < 0, "ZEB2-side rewiring", "ZEB1-side rewiring"))
)
setorder(sub_sum, residual_normal_median)
fwrite(sub_sum, file.path(out_dir, "subtype_residual_summary.tsv"), sep = "\t")

keep_cols <- c("sample_id", "subtype", want, paste0("z_", want),
               "dev", "balance", "expected_balance", "residual_normal",
               "expected_within", "residual_within", "dev_outside_normal_range")
fwrite(expr[, ..keep_cols], file.path(out_dir, "polonen_patient_residual.tsv"), sep = "\t")

meta <- data.table(
  item = c("early_genes", "dp_genes", "normal_units", "normal_r2",
           "polonen_n", "n_dev_outside_normal_range", "definition"),
  value = c(
    paste(early_genes, collapse = ","),
    paste(dp_genes, collapse = ","),
    as.character(nrow(normal_pts)),
    sprintf("%.4f", r2_normal),
    as.character(nrow(expr)),
    as.character(sum(expr$dev_outside_normal_range)),
    "dev = mean(z DP genes) - mean(z early genes); z within the dataset; expected balance from lm(balance ~ ns(dev, df=3)) fit on normal units only"
  )
)
fwrite(meta, file.path(out_dir, "frozen_score_definition.tsv"), sep = "\t")

# figures
curve_line <- normal_pts[order(dev), .(dev, expected_balance)]
p_curve <- ggplot(normal_pts, aes(dev, balance, colour = dataset)) +
  geom_point(size = 1.6, alpha = 0.85) +
  geom_line(data = curve_line, aes(dev, expected_balance), colour = "grey20", linewidth = 0.7, inherit.aes = FALSE) +
  labs(title = "Normal thymus: ZEB balance given developmental score",
       subtitle = "Curve is ns(dev, df=3) fit on these points only. ZEB genes are not in the score.",
       x = "Developmental score (within-dataset z)", y = "Balance z(ZEB1)-z(ZEB2)") +
  theme_bw(base_size = 11)
ggsave(file.path(fig_dir, "normal_balance_vs_dev.pdf"), p_curve, width = 8, height = 5.5)
ggsave(file.path(fig_dir, "normal_balance_vs_dev.png"), p_curve, width = 8, height = 5.5, dpi = 150)

ord <- sub_sum$subtype
plot_s <- expr[, .(subtype, residual_normal, balance, dev)]
plot_s$subtype <- factor(plot_s$subtype, levels = ord)
p_res <- ggplot(plot_s, aes(subtype, residual_normal)) +
  geom_hline(yintercept = 0, linewidth = 0.3, colour = "grey40") +
  geom_boxplot(outlier.size = 0.4, fill = "#e7eef5") +
  coord_flip() +
  labs(title = "Residual ZEB balance after the normal developmental curve",
       subtitle = "Negative means more ZEB2-dominant than the normal score predicts",
       x = NULL, y = "Actual balance minus normal-curve prediction") +
  theme_bw(base_size = 11)
ggsave(file.path(fig_dir, "subtype_residual_normal.pdf"), p_res, width = 8, height = 7)
ggsave(file.path(fig_dir, "subtype_residual_normal.png"), p_res, width = 8, height = 7, dpi = 150)

message("DONE")
print(summ)
print(r2_tab)
print(sub_sum)

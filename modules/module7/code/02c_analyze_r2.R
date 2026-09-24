# Module 7 round-2: pharmacotype PCA / clustering / ZEB1-core / MRD.
# mixOmics / NMF / pcaMethods / missMDA not installed.
# Backbone drugs = n>=50 complete-case. ZEB1 stays continuous.
suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(cluster)
})
has_cht <- requireNamespace("ComplexHeatmap", quietly = TRUE)
has_circlize <- requireNamespace("circlize", quietly = TRUE)

root <- "/data-b/liangfuhua/projects/ZEB1_TALL_analysis/module7"
m4p <- "/data-b/liangfuhua/projects/ZEB1_TALL_analysis/module4/processed"
m4t <- "/data-b/liangfuhua/projects/ZEB1_TALL_analysis/module4/tables"
m3p <- "/data-b/liangfuhua/projects/ZEB1_TALL_analysis/module3/processed"
fig <- file.path(root, "figures"); tab <- file.path(root, "tables"); proc <- file.path(root, "processed")
dir.create(fig, FALSE, TRUE); dir.create(tab, FALSE, TRUE)

theme_sci <- function(base_size = 11) {
  theme_classic(base_size = base_size) +
    theme(axis.text = element_text(color = "black"), axis.title = element_text(color = "black"),
          axis.line = element_line(linewidth = 0.45, color = "black"),
          plot.title = element_text(face = "bold", size = base_size + 1),
          plot.subtitle = element_text(size = base_size - 1, color = "grey25"),
          strip.background = element_blank())
}
save_both <- function(plot, file, width, height) {
  tryCatch({
    pdf(file, width = width, height = height, useDingbats = FALSE); print(plot); dev.off()
    png(sub("\\.pdf$", ".png", file), width = width * 160, height = height * 160, res = 160)
    print(plot); dev.off(); message("saved ", basename(file))
  }, error = function(e) message("FAIL ", basename(file), ": ", conditionMessage(e)))
}
run <- function(label, fn) tryCatch({ fn(); message("OK ", label) },
                                    error = function(e) message("FAIL ", label, ": ", conditionMessage(e)))
fmt_p <- function(p) ifelse(is.na(p), "NA", ifelse(p < 0.001, "P < 0.001", sprintf("P = %.3f", p)))
short_drug <- function(x) gsub(" \\(.*\\)$", "", gsub("_normalized$", "", x))

read_mat <- function(path) {
  d <- as.data.frame(fread(path), stringsAsFactors = FALSE)
  rownames(d) <- d[[1]]; d[[1]] <- NULL
  as.matrix(d)
}

ph <- fread(file.path(proc, "pharmacotype_tall_drugs.tsv"))
raw_drugs <- grep("\\((IU/ml|nM|µM|uM|μM)\\)$", names(ph), value = TRUE)
if ("Age at diagnosis (years)" %in% names(ph) && !"age_years" %in% names(ph))
  ph[, age_years := `Age at diagnosis (years)`]
if (!"etp" %in% names(ph) && "Immunophenotype" %in% names(ph))
  ph[, etp := fifelse(grepl("ETP", Immunophenotype, ignore.case = TRUE), "ETP", "notETP")]
if (!"subtype" %in% names(ph) && file.exists(file.path(m3p, "pharmacotype_tall_subtype.tsv"))) {
  lab <- fread(file.path(m3p, "pharmacotype_tall_subtype.tsv"))
  add <- intersect(c("sample", "subtype", "etp"), names(lab))
  ph <- merge(ph, lab[, ..add], by = "sample", all.x = TRUE, suffixes = c("", "_m3"))
}

drug_class <- data.table(
  drug_short = c("Vincristine", "Daunorubicin", "Prednisolone", "Dexamethasone",
                 "Asparaginase", "Mercaptopurine", "Thioguanine", "Cytarabine",
                 "Nelarabine", "Bortezomib", "Panobinostat", "Vorinostat",
                 "Ruxolitinib", "Dasatinib", "Ibrutinib", "Trametinib",
                 "Venetoclax", "CHZ868"),
  class = c("vinca", "anthracycline", "steroid", "steroid",
            "asparaginase", "antimetabolite", "antimetabolite", "antimetabolite",
            "nucleoside", "proteasome", "HDAC", "HDAC",
            "JAK", "TKI", "BTK", "MEK", "BCL2", "CK2")
)

# missingness
miss <- rbindlist(lapply(raw_drugs, function(d) {
  x <- as.numeric(ph[[d]])
  data.table(drug = d, drug_short = short_drug(d), n = sum(is.finite(x)),
             n_missing = sum(!is.finite(x)), frac_missing = mean(!is.finite(x)))
}))
fwrite(miss, file.path(tab, "M7_r2_missingness.tsv"), sep = "\t")

run("7r2 missingness", function() {
  d <- copy(miss)
  setorder(d, -frac_missing)
  d[, drug_short := factor(drug_short, levels = drug_short)]
  p <- ggplot(d, aes(drug_short, frac_missing, fill = n >= 50)) +
    geom_col(width = 0.75, color = "white") +
    coord_flip() +
    scale_fill_manual(values = c("TRUE" = "#4DBBD5", "FALSE" = "#E64B35"),
                      labels = c("TRUE" = "n>=50 backbone", "FALSE" = "n<50"), name = NULL) +
    labs(title = "Pharmacotype LC50 missingness in T-ALL",
         x = NULL, y = "Fraction missing") + theme_sci()
  save_both(p, file.path(fig, "M7_r2_7.1_missingness.pdf"), 7.2, 6.2)
})

backbone <- miss[n >= 50, drug]
message("backbone drugs: ", paste(short_drug(backbone), collapse = ", "))

# median-impute for clustering only; tests stay complete-case
X <- as.matrix(ph[, ..backbone])
storage.mode(X) <- "double"
Xz <- apply(X, 2, function(v) {
  m <- median(v, na.rm = TRUE)
  v[!is.finite(v)] <- m
  as.numeric(scale(v))
})
rownames(Xz) <- ph$sample
pca <- prcomp(Xz, center = FALSE, scale. = FALSE)
varp <- (pca$sdev^2) / sum(pca$sdev^2)
sc <- as.data.table(pca$x[, 1:2])
sc$sample <- ph$sample
sc$ZEB1_log <- ph$ZEB1_log
if ("subtype" %in% names(ph)) sc$subtype <- ph$subtype
if ("etp" %in% names(ph)) sc$etp <- ph$etp
fwrite(sc, file.path(tab, "M7_r2_drug_pca.tsv"), sep = "\t")

run("7r2 drug PCA", function() {
  d <- copy(sc)
  fill <- if ("subtype" %in% names(d)) "subtype" else if ("etp" %in% names(d)) "etp" else NULL
  p <- ggplot(d, aes(PC1, PC2)) +
    geom_point(aes_string(fill = fill, size = "ZEB1_log"), shape = 21, color = "grey20", alpha = 0.9) +
    scale_size(range = c(1.6, 5)) +
    labs(title = "Pharmacotype PCA of backbone drugs (median-imputed z)",
         subtitle = sprintf("PC1 %.1f%%  PC2 %.1f%%; n = %d; mixOmics/NMF not installed",
                            100 * varp[1], 100 * varp[2], nrow(d)),
         x = sprintf("PC1 (%.1f%%)", 100 * varp[1]),
         y = sprintf("PC2 (%.1f%%)", 100 * varp[2])) + theme_sci()
  save_both(p, file.path(fig, "M7_r2_7.2_drug_PCA.pdf"), 7.2, 5.8)
})

# k-means on PC1-3
kfit <- kmeans(pca$x[, 1:min(3, ncol(pca$x))], centers = 3, nstart = 25)
ph$drug_cluster <- paste0("C", kfit$cluster)
sc$drug_cluster <- ph$drug_cluster
fwrite(ph[, .(sample, ZEB1_log, drug_cluster)], file.path(tab, "M7_r2_drug_clusters.tsv"), sep = "\t")

run("7r2 cluster vs ZEB1", function() {
  d <- copy(ph)
  kw <- kruskal.test(ZEB1_log ~ drug_cluster, data = d)
  p <- ggplot(d, aes(drug_cluster, ZEB1_log, fill = drug_cluster)) +
    geom_boxplot(width = 0.55, outlier.shape = NA, alpha = 0.85) +
    geom_jitter(width = 0.12, size = 1.2, color = "grey20", alpha = 0.6) +
    scale_fill_manual(values = c(C1 = "#E64B35", C2 = "#4DBBD5", C3 = "#00A087"), guide = "none") +
    labs(title = "ZEB1 across drug-response k-means clusters",
         subtitle = sprintf("k = 3 on backbone PCA; Kruskal %s", fmt_p(kw$p.value)),
         x = NULL, y = "log2(FPKM+1) ZEB1") + theme_sci()
  save_both(p, file.path(fig, "M7_r2_7.7_ZEB1_by_drugcluster.pdf"), 6.2, 5.4)
})

# ZEB1 cores on pharmacotype expression
cores_pos <- character()
cores_neg <- character()
adjf <- file.path(m4t, "M4_r2_subtype_adjusted_ZEB1.tsv")
if (file.exists(adjf)) {
  st <- fread(adjf)
  cores_pos <- st[fdr_adj < 0.05 & t_adj > 0][order(-t_adj)][1:80, gene]
  cores_neg <- st[fdr_adj < 0.05 & t_adj < 0][order(t_adj)][1:80, gene]
}
pexpr_f <- file.path(m4p, "pharmacotype_log2fpkm_filtered.tsv.gz")
if (file.exists(pexpr_f) && length(cores_pos)) {
  pexpr <- read_mat(pexpr_f)
  score_one <- function(genes) {
    gg <- intersect(genes, rownames(pexpr))
    if (length(gg) < 8) return(setNames(rep(NA_real_, ncol(pexpr)), colnames(pexpr)))
    z <- t(scale(t(pexpr[gg, , drop = FALSE])))
    colMeans(z, na.rm = TRUE)
  }
  scp <- data.table(sample = colnames(pexpr),
                    core_pos = as.numeric(score_one(cores_pos)),
                    core_neg = as.numeric(score_one(cores_neg)))
  ph <- merge(ph, scp, by = "sample", all.x = TRUE)
  fwrite(ph[, .(sample, ZEB1_log, core_pos, core_neg, subtype, etp, age_years,
                `Day 15 MRD (%)`, `Day 42 or 46 MRD (%)`)],
         file.path(tab, "M7_r2_core_scores.tsv"), sep = "\t")
}

run("7r2 core vs ZEB1", function() {
  if (!"core_pos" %in% names(ph)) stop("no core")
  d <- ph[is.finite(ZEB1_log) & is.finite(core_pos)]
  ct <- cor.test(d$ZEB1_log, d$core_pos, method = "spearman", exact = FALSE)
  p <- ggplot(d, aes(ZEB1_log, core_pos, fill = etp)) +
    geom_point(shape = 21, size = 2.2, color = "grey20") +
    geom_smooth(method = "lm", se = TRUE, linewidth = 0.5, color = "grey20", inherit.aes = FALSE,
                aes(ZEB1_log, core_pos)) +
    scale_fill_manual(values = c(ETP = "#E64B35", notETP = "#4DBBD5"), na.value = "grey70") +
    labs(title = "Pharmacotype: M4 ZEB1-positive core vs ZEB1 RNA",
         subtitle = sprintf("n = %d; rho = %.2f, %s", nrow(d), unname(ct$estimate), fmt_p(ct$p.value)),
         x = "ZEB1", y = "Positive core", fill = NULL) + theme_sci()
  save_both(p, file.path(fig, "M7_r2_7.4_core_vs_ZEB1.pdf"), 6.6, 5.4)
})

# LC50 ~ core + subtype/ETP + age  (complete-case per drug)
reg_rows <- list()
for (d in raw_drugs) {
  dat <- copy(ph)
  dat[, lc50 := as.numeric(dat[[d]])]
  dat <- dat[is.finite(lc50) & is.finite(ZEB1_log)]
  if (nrow(dat) < 12) next
  fit1 <- tryCatch(lm(lc50 ~ ZEB1_log, dat), error = function(e) NULL)
  if (!is.null(fit1)) {
    s <- summary(fit1)$coefficients
    reg_rows[[length(reg_rows) + 1]] <- data.table(drug = d, model = "ZEB1", n = nrow(dat),
                                                   term = "ZEB1_log", beta = s["ZEB1_log", 1], p = s["ZEB1_log", 4])
  }
  if ("core_pos" %in% names(dat) && sum(is.finite(dat$core_pos)) >= 12) {
    d2 <- dat[is.finite(core_pos)]
    form <- "lc50 ~ core_pos"
    if ("etp" %in% names(d2) && uniqueN(d2$etp) > 1) form <- paste(form, "+ etp")
    if ("age_years" %in% names(d2) && sum(is.finite(d2$age_years)) >= 12) form <- paste(form, "+ age_years")
    fitc <- tryCatch(lm(as.formula(form), d2), error = function(e) NULL)
    if (!is.null(fitc) && "core_pos" %in% rownames(summary(fitc)$coefficients)) {
      s <- summary(fitc)$coefficients
      reg_rows[[length(reg_rows) + 1]] <- data.table(drug = d, model = form, n = nrow(d2),
                                                     term = "core_pos", beta = s["core_pos", 1], p = s["core_pos", 4])
    }
  }
}
reg <- rbindlist(reg_rows, fill = TRUE)
if (nrow(reg)) {
  reg[, fdr := p.adjust(p, "BH"), by = term]
  reg[, drug_short := short_drug(drug)]
  fwrite(reg, file.path(tab, "M7_r2_core_drug_lm.tsv"), sep = "\t")
}

run("7r2 core-drug heatmap of -logp", function() {
  if (!nrow(reg)) stop("no reg")
  d <- reg[term == "core_pos"]
  if (!nrow(d)) stop("no core models")
  setorder(d, beta)
  d[, drug_short := factor(drug_short, levels = unique(drug_short))]
  p <- ggplot(d, aes(beta, drug_short, fill = p < 0.05)) +
    geom_vline(xintercept = 0, linetype = 2) +
    geom_point(shape = 21, size = 3.2, color = "grey20") +
    scale_fill_manual(values = c("TRUE" = "#E64B35", "FALSE" = "#8491B4"), guide = "none") +
    labs(title = "LC50 ~ ZEB1-core (+ETP/age when available)",
         subtitle = "Complete-case; filled = nominal P < 0.05. MixOmics sPLS not installed.",
         x = "beta (core_pos)", y = NULL) + theme_sci()
  save_both(p, file.path(fig, "M7_r2_7.8_core_drug_beta.pdf"), 7.4, 6.4)
})

# drug-class: mean rho of ZEB1 vs LC50
scr <- fread(file.path(tab, "M7_7.7_7.8_drug_spearman.tsv"))
scr <- scr[gene == "ZEB1_log"]
scr[, drug_short := short_drug(drug)]
scr <- merge(scr, drug_class, by = "drug_short", all.x = TRUE)
fwrite(scr, file.path(tab, "M7_r2_ZEB1_rho_with_class.tsv"), sep = "\t")

run("7r2 drug class", function() {
  d <- scr[!is.na(class)]
  p <- ggplot(d, aes(class, rho, fill = class)) +
    geom_hline(yintercept = 0, linetype = 2, color = "grey50") +
    geom_boxplot(width = 0.55, outlier.shape = NA, alpha = 0.8) +
    geom_jitter(width = 0.12, size = 2, color = "grey20") +
    labs(title = "ZEB1–LC50 Spearman by drug class",
         subtitle = "Each point is one drug; class n is small",
         x = NULL, y = "rho (ZEB1 vs LC50)") + theme_sci() +
    theme(axis.text.x = element_text(angle = 35, hjust = 1), legend.position = "none")
  save_both(p, file.path(fig, "M7_r2_7.10_drug_class.pdf"), 8.0, 5.4)
})

run("7r2 MRD ~ core", function() {
  if (!"core_pos" %in% names(ph) || !"Day 15 MRD (%)" %in% names(ph)) stop("no MRD/core")
  d <- copy(ph)
  d[, mrd15 := as.numeric(`Day 15 MRD (%)`)]
  d <- d[is.finite(mrd15) & is.finite(core_pos)]
  fit <- lm(mrd15 ~ core_pos + age_years, d)
  s <- summary(fit)$coefficients
  ct <- cor.test(d$core_pos, d$mrd15, method = "spearman", exact = FALSE)
  fwrite(data.table(n = nrow(d), spearman_rho = unname(ct$estimate), spearman_p = ct$p.value,
                    lm_beta = s["core_pos", 1], lm_p = s["core_pos", 4]),
         file.path(tab, "M7_r2_MRD15_core.tsv"), sep = "\t")
  p <- ggplot(d, aes(core_pos, mrd15)) +
    geom_point(size = 2, alpha = 0.7, color = "#3C5488") +
    geom_smooth(method = "lm", se = TRUE, linewidth = 0.5, color = "grey20") +
    labs(title = "Day 15 MRD vs ZEB1-positive core",
         subtitle = sprintf("n = %d; Spearman rho = %.2f, %s", nrow(d), unname(ct$estimate), fmt_p(ct$p.value)),
         x = "ZEB1-positive core", y = "Day 15 MRD (%)") + theme_sci()
  save_both(p, file.path(fig, "M7_r2_7.15_MRD_core.pdf"), 6.4, 5.4)
})

writeLines(c(
  "mixOmics / NMF / pcaMethods / missMDA not installed.",
  "PCA used prcomp on median-imputed z of backbone drugs (n>=50).",
  "Clusters: k-means k=3 on PC1-3.",
  "sPLS skipped. PRISM still absent (see Module 8)."
), file.path(tab, "M7_r2_methods_note.txt"))
message("M7 R2 DONE")

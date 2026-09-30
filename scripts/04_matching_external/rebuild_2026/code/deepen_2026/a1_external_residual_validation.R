# Frozen within-cohort residual applied to independent pediatric T-ALL RNA.
# Do not invent a new score. Each cohort is standardized and spline-fit alone.
# Labels from clinTALL S9 are classifier-assigned and cover only the
# fusion-annotated AIEOP subset. Full 17-subtype rank concordance waits for
# RNA-only clinTALL predictions on all 120 / 108 samples.
suppressPackageStartupMessages({
  library(data.table)
  library(edgeR)
  library(splines)
  library(readr)
})

argv <- grep("^--file=", commandArgs(FALSE), value = TRUE)
script <- normalizePath(sub("^--file=", "", argv), winslash = "/")
code <- dirname(script)
rebuild <- normalizePath(file.path(code, "../.."), winslash = "/")
root <- normalizePath(file.path(rebuild, "../.."), winslash = "/")
out <- file.path(rebuild, "data/deepen_2026/a1_external_residual")
dir.create(out, recursive = TRUE, showWarnings = FALSE)

genes <- c(ZEB1 = "ENSG00000148516", ZEB2 = "ENSG00000169554",
           CD1A = "ENSG00000158477", CD34 = "ENSG00000174059",
           LYL1 = "ENSG00000104903")
zsd <- function(x) {
  x <- as.numeric(x)
  s <- stats::sd(x, na.rm = TRUE)
  if (!is.finite(s) || s == 0) return(rep(NA_real_, length(x)))
  (x - mean(x, na.rm = TRUE)) / s
}
score <- function(dt) {
  dt[, `:=`(z_ZEB1 = zsd(ZEB1), z_ZEB2 = zsd(ZEB2), z_CD1A = zsd(CD1A),
            z_CD34 = zsd(CD34), z_LYL1 = zsd(LYL1))]
  dt[, `:=`(
    balance = z_ZEB1 - z_ZEB2,
    dev = z_CD1A - (z_CD34 + z_LYL1) / 2
  )]
  fit <- lm(balance ~ ns(dev, df = 3), data = dt)
  dt[, expected := as.numeric(predict(fit, newdata = dt))]
  dt[, residual := balance - expected]
  list(data = dt, fit = fit, r2 = summary(fit)$r.squared)
}
boot_median <- function(x, B = 3000L, seed = 20260926L) {
  x <- x[is.finite(x)]
  if (!length(x)) return(c(NA_real_, NA_real_, NA_real_))
  if (length(x) == 1L) return(c(x[1], NA_real_, NA_real_))
  set.seed(seed)
  m <- replicate(B, median(sample(x, length(x), replace = TRUE)))
  c(median(x), quantile(m, 0.025, names = FALSE), quantile(m, 0.975, names = FALSE))
}

# ---- discovery reference ----
disc <- fread(file.path(root, "total/data/validation/zeb_developmental_residual/polonen_patient_residual.tsv"),
              encoding = "UTF-8")
stopifnot(nrow(disc) == 1309L, uniqueN(disc$sample_id) == 1309L)
disc_med <- disc[, .(n = .N, residual_median = median(residual_within)), by = subtype]
setorder(disc_med, residual_median)
fwrite(disc_med, file.path(out, "polonen_discovery_subtype_medians.tsv"), sep = "\t")

# ---- AIEOP counts: TMM logCPM, then the frozen within-cohort residual ----
cnt_path <- file.path(rebuild, "data/deepen_2026/clintall/counts_matrix.csv")
cnt <- fread(cnt_path, sep = ";")
stopifnot(names(cnt)[1] == "Geneid")
sid <- names(cnt)[-1]
stopifnot(length(sid) == 120L)
mat <- as.matrix(cnt[, ..sid])
storage.mode(mat) <- "integer"
rownames(mat) <- cnt$Geneid
miss <- setdiff(unname(genes), rownames(mat))
if (length(miss)) stop("AIEOP missing genes: ", paste(miss, collapse = ","))
y <- DGEList(counts = mat)
y <- calcNormFactors(y, method = "TMM")
lcpm <- cpm(y, log = TRUE, prior.count = 1)
aie <- data.table(sample_id = sid)
for (g in names(genes)) aie[[g]] <- as.numeric(lcpm[genes[[g]], sid])
aie_scored <- score(aie)
aie <- aie_scored$data
aie[, cohort := "AIEOP-BFM_ALL2017"]
aie[, expression_scale := "TMM_log2CPM_prior1"]

lab <- fread(file.path(rebuild, "data/deepen_2026/clintall/aieop_s9_labels.tsv"))
aie <- merge(aie, lab, by = "sample_id", all.x = TRUE, sort = FALSE)
aie[, predicted_subtype := fifelse(is.na(predicted_subtype) | predicted_subtype == "",
                                   NA_character_, predicted_subtype)]
fwrite(aie, file.path(out, "aieop_patient_residual.tsv"), sep = "\t")

labeled <- aie[!is.na(predicted_subtype)]
aie_med <- labeled[, {
  ci <- boot_median(residual)
  .(n = .N, residual_median = ci[1], ci_low = ci[2], ci_high = ci[3],
    balance_median = median(balance), dev_median = median(dev))
}, by = predicted_subtype]
setnames(aie_med, "predicted_subtype", "subtype")
setorder(aie_med, residual_median)
fwrite(aie_med, file.path(out, "aieop_s9_subtype_residual_medians.tsv"), sep = "\t")

# Pole tests that are actually available in S9.
pole <- list()
add_pole <- function(name, x, expected_sign) {
  ci <- boot_median(x)
  pole[[length(pole) + 1L]] <<- data.table(
    contrast = name, n = sum(is.finite(x)), estimate = ci[1],
    ci_low = ci[2], ci_high = ci[3], expected_sign = expected_sign,
    interval_excludes_zero = is.finite(ci[2]) && (ci[2] * ci[3] > 0),
    direction_matches = is.finite(ci[1]) && sign(ci[1]) == expected_sign
  )
}
add_pole("AIEOP_S9_TLX3_median_residual", labeled[predicted_subtype == "TLX3", residual], 1)
add_pole("AIEOP_S9_ETP_like_median_residual", labeled[predicted_subtype == "ETP-like", residual], 0)
etp <- labeled[predicted_subtype == "ETP-like"]
if (nrow(etp)) {
  add_pole("AIEOP_S9_ETP_like_abs_residual_median", abs(etp$residual), -1)
}

# ---- NOPHO VST: same z/spline algebra on the deposited scale ----
vst_path <- file.path(root, "total/data/validation/GSE272023_MatrixRNAseqVSTcountsFilted.txt.gz")
meta <- fread(file.path(root, "total/data/validation/gse272023_sample_map.tsv"))
tum <- meta[tissue_type == "Tumor"]
stopifnot(nrow(tum) == 108L, !anyDuplicated(tum$sample))
vst <- fread(vst_path)
setnames(vst, 1L, "ID")
vst[, ID := as.character(ID)]
# GEO matrix also carries CD3/CD34/Lymphnode annotation columns; keep tumors only.
keep <- intersect(tum$sample, names(vst))
stopifnot(length(keep) == 108L)
hit <- vst[ID %in% unname(genes), c("ID", keep), with = FALSE]
stopifnot(nrow(hit) == 5L, !anyNA(hit))
ord <- match(unname(genes), hit$ID)
stopifnot(!anyNA(ord))
expr <- as.matrix(hit[ord, ..keep])
storage.mode(expr) <- "double"
rownames(expr) <- names(genes)
nop <- data.table(sample_id = keep, gsm = tum$gsm[match(keep, tum$sample)],
                  cimp = tum$cimp[match(keep, tum$sample)])
for (g in names(genes)) nop[[g]] <- as.numeric(expr[g, keep])
nop_scored <- score(nop)
nop <- nop_scored$data
nop[, `:=`(cohort = "NOPHO_ALL2008_GSE272023",
           expression_scale = "deposited_VST",
           predicted_subtype = NA_character_)]
fwrite(nop, file.path(out, "nopho_patient_residual.tsv"), sep = "\t")

# ---- concordance is undefined until full classifier labels exist ----
ov <- merge(disc_med, aie_med[, .(subtype, n_external = n, residual_median_external = residual_median)],
            by = "subtype")
if (nrow(ov) >= 3L) {
  ov$spearman_note <- "fusion-annotated S9 subset only; not a 17-subtype replication"
  fwrite(ov, file.path(out, "aieop_s9_vs_polonen_overlap.tsv"), sep = "\t")
}

qc <- rbindlist(pole, fill = TRUE)
qc_head <- data.table(
  item = c("aieop_n", "aieop_s9_labeled_n", "aieop_s9_has_BCL11B",
           "aieop_dev_r2", "nopho_n", "nopho_dev_r2",
           "nopho_predicted_subtype", "full_17subtype_rank_ready"),
  value = c(nrow(aie), nrow(labeled),
            as.character(any(labeled$predicted_subtype == "BCL11B", na.rm = TRUE)),
            sprintf("%.4f", aie_scored$r2), nrow(nop),
            sprintf("%.4f", nop_scored$r2), "not_published", "FALSE")
)
fwrite(qc, file.path(out, "available_pole_tests.tsv"), sep = "\t")
fwrite(qc_head, file.path(out, "a1_status.tsv"), sep = "\t")
writeLines(c(
  "A1 applies the frozen residual algebra inside each external cohort.",
  "B = z(ZEB1)-z(ZEB2); D = z(CD1A)-[z(CD34)+z(LYL1)]/2; R = B - ns(D, df=3).",
  "Z-scores and the spline are fit inside the external cohort. Discovery coefficients are not reused.",
  "AIEOP uses featureCounts integers and TMM log2CPM (prior.count=1), matching the Pölönen residual script.",
  "NOPHO uses deposited VST values. That scale is not TMM logCPM; only within-cohort ranks are comparable.",
  "clinTALL S9 publishes predictions for the fusion-annotated AIEOP subset, not all 120 patients, and contains no BCL11B call.",
  "Therefore this run cannot claim independent 17-subtype residual ranking. It prepares the frozen residual layer.",
  "Next required input: RNA-only clinTALL predictions for AIEOP n=120 and, if integer counts become available, NOPHO n=108.",
  "Write labels as classifier-assigned external subtype validation, not genomic confirmation."
), file.path(out, "METHODS.txt"))
writeLines(capture.output(sessionInfo()), file.path(out, "R_sessionInfo.txt"))
cat("A1 residual layer complete\n")
print(qc_head)
print(aie_med)

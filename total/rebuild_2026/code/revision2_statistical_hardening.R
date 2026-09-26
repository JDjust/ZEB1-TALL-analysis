# Revision 2: patient-stratified bootstrap of the frozen Pölönen analysis.
# No genes, samples, outcomes or subtypes are added. Run with Rscript [script] [B].
suppressPackageStartupMessages({library(splines); library(sandwich); library(readr)})
args <- commandArgs(trailingOnly = TRUE)
B <- if (length(args)) as.integer(args[[1]]) else 3000L
stopifnot(is.finite(B), B >= 100L)
argv <- grep("^--file=", commandArgs(FALSE), value = TRUE)
script <- normalizePath(sub("^--file=", "", argv), winslash = "/")
base <- normalizePath(file.path(dirname(script), ".."), winslash = "/")
input <- file.path(base, "../data/validation/zeb_developmental_residual/polonen_patient_residual.tsv")
out <- file.path(base, "data/source_data_rebuilt/revision2_statistics")
dir.create(out, recursive = TRUE, showWarnings = FALSE)
d <- as.data.frame(readr::read_tsv(input, locale = readr::locale(encoding = "UTF-8"),
                                   show_col_types = FALSE, progress = FALSE,
                                   name_repair = "minimal"))
stopifnot(nrow(d) == 1309L, length(unique(d$sample_id)) == 1309L,
          length(unique(d$subtype)) == 17L, all(is.finite(d$dev)),
          all(is.finite(d$balance)))
lev <- sort(unique(d$subtype))
d$subtype <- factor(d$subtype, levels = lev)
contrasts(d$subtype) <- contr.sum(length(lev))

fit_dev <- lm(balance ~ ns(dev, df = 3), data = d)
fit_both <- lm(balance ~ ns(dev, df = 3) + subtype, data = d)
stopifnot(max(abs(residuals(fit_dev) - d$residual_within)) < 1e-7,
          abs(summary(fit_dev)$r.squared - 0.255563) < 1e-5,
          abs(summary(fit_both)$r.squared - 0.394922) < 1e-5)

# Sum contrasts yield effects relative to the equal-weighted average of the
# 17 subtype intercepts, conditional on the common developmental spline.
effect_matrix <- function(model) {
  cn <- names(coef(model))
  C <- matrix(0, nrow = length(lev), ncol = length(cn),
              dimnames = list(lev, cn))
  sc <- grep("^subtype", cn)
  stopifnot(length(sc) == length(lev) - 1L)
  C[seq_len(length(lev) - 1L), sc] <- diag(length(lev) - 1L)
  C[length(lev), sc] <- -1
  C
}
C <- effect_matrix(fit_both)
beta <- as.vector(C %*% coef(fit_both))
V <- sandwich::vcovHC(fit_both, type = "HC3")
se <- sqrt(diag(C %*% V %*% t(C)))
p <- 2 * pt(abs(beta / se), df = df.residual(fit_both), lower.tail = FALSE)
names(beta) <- names(se) <- names(p) <- lev
point_median <- tapply(residuals(fit_dev), d$subtype, median)[lev]
n <- as.integer(table(d$subtype)[lev])

set.seed(20260926L)
idx <- split(seq_len(nrow(d)), d$subtype)
boot_med <- matrix(NA_real_, B, length(lev), dimnames = list(NULL, lev))
boot_eff <- matrix(NA_real_, B, length(lev), dimnames = list(NULL, lev))
for (b in seq_len(B)) {
  draw <- unlist(lapply(idx, function(ii) sample(ii, length(ii), replace = TRUE)),
                 use.names = FALSE)
  x <- d[draw, , drop = FALSE]
  contrasts(x$subtype) <- contr.sum(length(lev))
  fd <- lm(balance ~ ns(dev, df = 3), data = x)
  fm <- lm(balance ~ ns(dev, df = 3) + subtype, data = x)
  boot_med[b, ] <- tapply(residuals(fd), x$subtype, median)[lev]
  boot_eff[b, ] <- as.vector(effect_matrix(fm) %*% coef(fm))
}
stopifnot(!anyNA(boot_med), !anyNA(boot_eff))
ci <- function(m) t(apply(m, 2, quantile, probs = c(.025, .975), names = FALSE))
med_ci <- ci(boot_med)
eff_ci <- ci(boot_eff)
ans <- data.frame(subtype = lev, n = n,
  residual_within_median = as.numeric(point_median),
  residual_within_ci_low = med_ci[, 1], residual_within_ci_high = med_ci[, 2],
  adjusted_effect = beta,
  adjusted_effect_ci_low = eff_ci[, 1],
  adjusted_effect_ci_high = eff_ci[, 2],
  adjusted_effect_hc3_se = se,
  adjusted_effect_hc3_p = p,
  adjusted_effect_bh_fdr = p.adjust(p, method = "BH"),
  stringsAsFactors = FALSE)
ans <- ans[order(ans$residual_within_median), ]
readr::write_tsv(ans, file.path(out, "subtype_effects_bootstrap.tsv"), na = "NA")
rep <- data.frame(replicate = rep(seq_len(B), each = length(lev)),
  subtype = rep(lev, times = B),
  residual_within_median = as.vector(t(boot_med)),
  adjusted_effect = as.vector(t(boot_eff)))
readr::write_tsv(rep, file.path(out, "patient_bootstrap_replicates.tsv"), na = "NA")
global <- anova(fit_dev, fit_both)
qc <- data.frame(n_patients = nrow(d), n_subtypes = length(lev),
  bootstrap_replicates = B, seed = 20260926L,
  dev_r2 = summary(fit_dev)$r.squared,
  combined_r2 = summary(fit_both)$r.squared,
  incremental_r2 = summary(fit_both)$r.squared - summary(fit_dev)$r.squared,
  global_subtype_p = global$`Pr(>F)`[2],
  max_frozen_residual_difference = max(abs(residuals(fit_dev) - d$residual_within)))
write.table(qc, file.path(out, "model_qc.tsv"), sep = "\t", quote = FALSE,
            row.names = FALSE)
writeLines(c(
  "Primary inferential population: 1,309 unique diagnostic Pölönen patients; one RNA sample per patient.",
  "Balance and developmental coordinate are the frozen, separately cohort-standardized quantities.",
  "Developmental model: balance ~ ns(dev, df=3), fitted within the Pölönen cohort.",
  "Primary descriptive residual: observed balance minus that cohort-fitted expected balance.",
  "Bootstrap: within-subtype patient resampling with replacement; fixed subtype sample sizes; refit both models each replicate; percentile 95% CI; fixed seed 20260926.",
  "Adjusted effect: coefficient from balance ~ ns(dev, df=3) + subtype with sum-to-zero contrasts, relative to the unweighted average of 17 subtype intercepts.",
  "Adjusted-effect P values: two-sided HC3 robust t tests; BH correction across 17 subtype effects. Confidence intervals are patient-bootstrap percentiles.",
  "The global nested-model F test is distinct from 17 individual effects. The median residual and adjusted coefficient are different estimands and are labelled separately.",
  "Normal-thymus spline projection is retained as a cross-cohort sensitivity and must not be interpreted as an absolute biological offset."),
  file.path(out, "METHODS.txt"))
writeLines(capture.output(sessionInfo()), file.path(out, "R_sessionInfo.txt"))
cat(sprintf("Revision 2 completed: B=%d; n=%d; global P=%.4g\n", B, nrow(d),
            global$`Pr(>F)`[2]))

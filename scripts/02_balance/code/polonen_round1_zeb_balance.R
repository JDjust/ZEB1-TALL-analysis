# Pölönen syn54032669 round 1.
# Frozen metric: balance = z(logCPM ZEB1) - z(logCPM ZEB2), z within the 1309-patient set.
# No reweighting. Primary binary endpoint is the author Induction Failure label.

suppressPackageStartupMessages({
  library(data.table)
  library(edgeR)
  library(logistf)
  library(survival)
})

seed <- 20260925
set.seed(seed)
out_dir <- "D:/_bioinformation/ZEB1/data/analysis/total/data/validation/polonen_round1"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
data_dir <- "D:/_bioinformation/ZEB1/data/analysis/data/polonen_syn54032669"

e <- new.env(parent = emptyenv())
load(file.path(data_dir, "Data_1309Samples.RData"), envir = e)
annot <- e$annot
gexp <- e$gexp
sub_f <- e$M3.subtype.factor
ids <- rownames(annot)
stopifnot(length(unique(annot$subject_sjuid)) == nrow(annot))

message("reading counts")
ann_g <- fread("D:/_bioinformation/ZEB1/data/analysis/data/polonen_syn54032669/TALL_X01_gene_annotations_unfiltered.tsv")
counts <- fread("D:/_bioinformation/ZEB1/data/analysis/data/TALL_X01_counts.tsv")
gene_col <- names(counts)[1]
setnames(counts, gene_col, "ensembl")
stopifnot(all(ids %in% names(counts)))
mat <- as.matrix(counts[, ..ids])
storage.mode(mat) <- "double"
rownames(mat) <- counts$ensembl
# symbol from gexp row order if lengths match, else annotation
sym <- NULL
if (nrow(mat) == nrow(gexp)) {
  # confirm first symbols via annotation if possible
  sym <- rownames(gexp)
}
want <- c(ZEB1 = "ZEB1", ZEB2 = "ZEB2", LMO2 = "LMO2", ZBTB16 = "ZBTB16")
hit <- sapply(want, function(g) which(sym == g))
stopifnot(all(lengths(hit) == 1))

message("TMM logCPM")
y <- DGEList(counts = mat)
y <- calcNormFactors(y, method = "TMM")
lcpm <- cpm(y, log = TRUE, prior.count = 1)
logcpm <- sapply(want, function(g) lcpm[hit[[g]], ])
logcpm <- as.data.frame(logcpm)
logcpm$sample_id <- ids

zsd <- function(x) (x - mean(x)) / sd(x)
logcpm$z_ZEB1 <- zsd(logcpm$ZEB1)
logcpm$z_ZEB2 <- zsd(logcpm$ZEB2)
logcpm$balance <- logcpm$z_ZEB1 - logcpm$z_ZEB2
logcpm$z_LMO2 <- zsd(logcpm$LMO2)
logcpm$z_ZBTB16 <- zsd(logcpm$ZBTB16)

# author matrix on the same patients, z within cohort, sensitivity only
gx <- data.frame(sample_id = ids, stringsAsFactors = FALSE)
for (g in c("ZEB1", "ZEB2", "LMO2", "ZBTB16")) gx[[g]] <- as.numeric(gexp[g, ids])
gx$balance_gexp <- zsd(gx$ZEB1) - zsd(gx$ZEB2)

d <- data.frame(
  sample_id = ids,
  subject_sjuid = annot$subject_sjuid,
  specimen = annot$Tumor.Specimen.Type,
  age = annot$Age.at.Diagnosis.in.Years,
  sex = annot$Sex,
  wbc = annot$WBC.at.Diagnosis,
  risk = annot$Risk.Group,
  mrd = annot$Day.29.MRD,
  morph = annot$Day.29.morphologic.Response,
  etp = as.character(e$M2.ETP_status.factor),
  subtype = as.character(sub_f$Reviewed.subtype),
  protocol = annot$Clinical.Trial.Protocol.Name,
  event_type = annot$Event.Type,
  efs_time = annot$EFS,
  efs_status = annot$EFS.status,
  os_time = annot$OS,
  os_status = annot$OS.status,
  stringsAsFactors = FALSE
)
d <- merge(d, logcpm, by = "sample_id", sort = FALSE)
d <- merge(d, gx[, c("sample_id", "balance_gexp")], by = "sample_id", sort = FALSE)
d$balance_gexp_z <- d$balance_gexp  # already z(ZEB1)-z(ZEB2) on gexp

# WBC: cohort values are typically tens to hundreds; one 325000 is a unit break
d$wbc_use <- d$wbc
d$wbc_use[d$wbc_use > 2000] <- NA
d$log10_wbc <- log10(d$wbc_use)

d$if_author <- as.integer(d$risk == "Induction Failure")
d$m3 <- as.integer(d$morph == "M3")
d$poor_morph <- ifelse(d$morph %in% c("M2", "M3"), 1L, ifelse(d$morph == "M1", 0L, NA_integer_))
d$m3_known <- ifelse(d$morph == "M3", 1L, ifelse(d$morph %in% c("M1", "M2"), 0L, NA_integer_))

mrd_max <- max(d$mrd, na.rm = TRUE)
# Stored numbers match percent (0.01, 0.1, ...), and the paper's outcome cut is 0.1%.
# Residual-disease cut in Fig. 5 is 0.01%. Both are author thresholds.
d$mrd_pos_0.1 <- ifelse(is.na(d$mrd), NA_integer_, as.integer(d$mrd >= 0.1))
d$mrd_pos_0.01 <- ifelse(is.na(d$mrd), NA_integer_, as.integer(d$mrd >= 0.01))

# primary analysis set: author IF label is defined (not the outcome itself missing)
# Inevaluable / Not Assigned / PH+ stay in if their label is not Induction Failure,
# because the author still classified them. Unknown morph is reported separately.
d$primary_ok <- !is.na(d$if_author) & d$risk != ""

fwrite(d, file.path(out_dir, "patient_level_frozen_balance.tsv"), sep = "\t")

audit <- data.frame(
  item = c(
    "expression samples in TALL_X01_counts",
    "Data_1309Samples patients",
    "unique subjects in 1309",
    "subjects with >1 RNA row in 1309",
    "1309 samples present in counts",
    "ZEB1 ZEB2 LMO2 ZBTB16 in matrix",
    "author  subtype assigned",
    "Day 29 MRD non-missing",
    "author Induction Failure",
    "Day 29 morphologic M3",
    "Day 29 morphologic M2",
    "EFS complete (time and status)",
    "OS complete (time and status)",
    "MRD stored maximum",
    "TMM norm factor range"
  ),
  n = c(
    ncol(counts) - 1,
    nrow(d),
    length(unique(d$subject_sjuid)),
    0,
    sum(ids %in% names(counts)),
    1309,
    sum(!is.na(d$subtype) & d$subtype != ""),
    sum(!is.na(d$mrd)),
    sum(d$if_author == 1),
    sum(d$morph == "M3"),
    sum(d$morph == "M2"),
    sum(!is.na(d$efs_time) & !is.na(d$efs_status)),
    sum(!is.na(d$os_time) & !is.na(d$os_status)),
    mrd_max,
    NA_real_
  ),
  stringsAsFactors = FALSE
)
audit$n[nrow(audit)] <- paste(round(range(y$samples$norm.factors), 3), collapse = " to ")
# fix type
audit$n <- as.character(audit$n)
fwrite(audit, file.path(out_dir, "sample_audit.tsv"), sep = "\t")

cor_note <- cor(d$balance, d$balance_gexp, method = "spearman")
message("spearman balance TMM vs gexp ", round(cor_note, 4))
message("IF vs M3 agreement ", sum(d$if_author == 1 & d$morph == "M3"), " / IF ", sum(d$if_author == 1), " / M3 ", sum(d$morph == "M3"))
message("MRD max ", mrd_max, " MRD>=0.1 ", sum(d$mrd_pos_0.1 == 1, na.rm = TRUE), " MRD>=0.01 ", sum(d$mrd_pos_0.01 == 1, na.rm = TRUE))

or_row <- function(fit, term, model, n, events, method) {
  co <- coef(fit)
  se <- sqrt(diag(vcov(fit)))
  b <- unname(co[term])
  s <- unname(se[term])
  data.frame(
    model = model, method = method, term = term, n = n, events = events,
    log_or = b, or = exp(b), ci_low = exp(b - 1.96 * s), ci_high = exp(b + 1.96 * s),
    p = 2 * pnorm(-abs(b / s)),
    stringsAsFactors = FALSE
  )
}

boot_or <- function(df, formula, term, B = 2000) {
  rng <- sample.int(.Machine$integer.max, 1)
  set.seed(seed)
  n <- nrow(df)
  est <- rep(NA_real_, B)
  for (i in seq_len(B)) {
    ii <- sample.int(n, n, replace = TRUE)
    sub <- df[ii, , drop = FALSE]
    if (length(unique(sub[[all.vars(formula)[1]]])) < 2) next
    fit <- try(glm(formula, data = sub, family = binomial()), silent = TRUE)
    if (inherits(fit, "try-error")) next
    est[i] <- unname(coef(fit)[term])
  }
  est <- est[is.finite(est)]
  quantile(exp(est), c(0.025, 0.975), na.rm = TRUE)
}

prep <- function(df) {
  df <- df[is.finite(df$balance) & !is.na(df$if_author), ]
  df$subtype <- factor(df$subtype)
  df$sex <- factor(df$sex)
  df
}

base <- prep(d)
# drop subtype levels that are unused
base$subtype <- droplevels(base$subtype)

message("model 0")
m0 <- glm(if_author ~ balance, data = base, family = binomial())
r0 <- or_row(m0, "balance", "M0 failure ~ balance", nrow(base), sum(base$if_author), "glm")
b0 <- boot_or(base, if_author ~ balance, "balance", 2000)
r0$boot_low <- unname(b0[1])
r0$boot_high <- unname(b0[2])

message("model 1 Firth + subtype")
m1 <- logistf(if_author ~ balance + subtype, data = base, pl = FALSE)
r1 <- or_row(m1, "balance", "M1 failure ~ balance + subtype", nrow(base), sum(base$if_author), "firth")

message("model 2")
b2 <- base[is.finite(base$age) & !is.na(base$sex), ]
m2 <- logistf(if_author ~ balance + subtype + age + sex, data = b2, pl = FALSE)
r2 <- or_row(m2, "balance", "M2 failure ~ balance + subtype + age + sex", nrow(b2), sum(b2$if_author), "firth")

message("model 3")
b3 <- b2[is.finite(b2$log10_wbc), ]
m3 <- logistf(if_author ~ balance + subtype + age + sex + log10_wbc, data = b3, pl = FALSE)
r3 <- or_row(m3, "balance", "M3 failure ~ balance + subtype + age + sex + log10 WBC", nrow(b3), sum(b3$if_author), "firth")

# same-n attenuation: refit M0 on M3 complete cases
m0_on3 <- glm(if_author ~ balance, data = b3, family = binomial())
r0b <- or_row(m0_on3, "balance", "M0 on M3 complete cases", nrow(b3), sum(b3$if_author), "glm")

rows <- rbindlist(list(r0, r1, r2, r3, r0b), fill = TRUE)
# single-gene Firth unadjusted, frozen, not for model selection
for (term in c("z_ZEB1", "z_ZEB2", "z_LMO2")) {
  f <- glm(reformulate(term, "if_author"), data = base, family = binomial())
  rows <- rbindlist(list(rows, or_row(f, term, paste0("M0 failure ~ ", term), nrow(base), sum(base$if_author), "glm")), fill = TRUE)
}

# MRD endpoints
mrd_df <- d[!is.na(d$mrd) & is.finite(d$balance), ]
sp <- cor.test(mrd_df$balance, mrd_df$mrd, method = "spearman", exact = FALSE)
mrd_rows <- list()
for (ep in c("mrd_pos_0.1", "mrd_pos_0.01")) {
  md <- mrd_df
  md$y <- md[[ep]]
  md$subtype <- droplevels(factor(md$subtype))
  f0 <- glm(y ~ balance, data = md, family = binomial())
  f1 <- logistf(y ~ balance + subtype, data = md, pl = FALSE)
  a <- or_row(f0, "balance", paste0(ep, " ~ balance"), nrow(md), sum(md$y), "glm")
  bci <- boot_or(md, y ~ balance, "balance", 1000)
  a$boot_low <- unname(bci[1]); a$boot_high <- unname(bci[2])
  c <- or_row(f1, "balance", paste0(ep, " ~ balance + subtype"), nrow(md), sum(md$y), "firth")
  mrd_rows[[ep]] <- rbindlist(list(a, c), fill = TRUE)
}
# poor morph sensitivity
pm <- d[!is.na(d$poor_morph) & is.finite(d$balance), ]
pm$y <- pm$poor_morph
pm$subtype <- droplevels(factor(pm$subtype))
fp <- glm(y ~ balance, data = pm, family = binomial())
fp1 <- logistf(y ~ balance + subtype, data = pm, pl = FALSE)
poor <- rbindlist(list(
  or_row(fp, "balance", "M2orM3 ~ balance", nrow(pm), sum(pm$y), "glm"),
  or_row(fp1, "balance", "M2orM3 ~ balance + subtype", nrow(pm), sum(pm$y), "firth")
), fill = TRUE)

# EFS / OS Cox, per 1 SD balance. Author outcome models used Firth Cox; survival::coxph is the available implementation.
cox_row <- function(time, status, df, model) {
  df <- df[is.finite(df[[time]]) & !is.na(df[[status]]) & is.finite(df$balance), ]
  df$subtype <- droplevels(factor(df$subtype))
  df$sex <- factor(df$sex)
  f0 <- coxph(as.formula(paste0("Surv(", time, ",", status, ") ~ balance")), data = df)
  f1 <- coxph(as.formula(paste0("Surv(", time, ",", status, ") ~ balance + subtype")), data = df)
  grab <- function(fit, label) {
    b <- coef(fit)["balance"]
    s <- sqrt(vcov(fit)["balance", "balance"])
    data.frame(model = label, n = nrow(df), events = sum(df[[status]] == 1),
               hr = exp(b), ci_low = exp(b - 1.96 * s), ci_high = exp(b + 1.96 * s),
               p = 2 * pnorm(-abs(b / s)), stringsAsFactors = FALSE)
  }
  rbind(grab(f0, paste0(model, " ~ balance")), grab(f1, paste0(model, " ~ balance + subtype")))
}
surv_tab <- rbind(cox_row("efs_time", "efs_status", d, "EFS"), cox_row("os_time", "os_status", d, "OS"))

# ETP decomposition
etp_ok <- d[d$etp %in% c("ETP", "Near-ETP", "Non-ETP") & is.finite(d$balance), ]
welch <- function(a, b) {
  t <- t.test(a, b)
  data.frame(n_a = length(a), n_b = length(b), mean_a = mean(a), mean_b = mean(b),
             diff = mean(a) - mean(b), ci_low = t$conf.int[1], ci_high = t$conf.int[2], p = t$p.value)
}
etp_contrast <- rbind(
  cbind(contrast = "ETP vs Non-ETP", welch(etp_ok$balance[etp_ok$etp == "ETP"], etp_ok$balance[etp_ok$etp == "Non-ETP"])),
  cbind(contrast = "Near-ETP vs Non-ETP", welch(etp_ok$balance[etp_ok$etp == "Near-ETP"], etp_ok$balance[etp_ok$etp == "Non-ETP"])),
  cbind(contrast = "ETP-like subtype vs other", welch(d$balance[d$subtype == "ETP-like"], d$balance[d$subtype != "ETP-like"]))
)
# within ETP / non-ETP outcome
within_etp <- list()
for (lab in c("ETP", "Non-ETP", "Near-ETP")) {
  keep <- d$etp == lab & is.finite(d$balance) & !is.na(d$if_author)
  keep[is.na(keep)] <- FALSE
  w <- d[keep, ]
  evn <- sum(w$if_author)
  if (evn < 3 || evn > nrow(w) - 3) {
    within_etp[[lab]] <- data.frame(group = lab, n = nrow(w), events = sum(w$if_author),
                                    or = NA, ci_low = NA, ci_high = NA, p = NA)
  } else {
    f <- glm(if_author ~ balance, data = w, family = binomial())
    rr <- or_row(f, "balance", lab, nrow(w), sum(w$if_author), "glm")
    within_etp[[lab]] <- data.frame(group = lab, n = rr$n, events = rr$events, or = rr$or,
                                    ci_low = rr$ci_low, ci_high = rr$ci_high, p = rr$p)
  }
}
within_etp <- rbindlist(within_etp)

# does ETP label add variance in balance after subtype?
dd <- d[d$etp %in% c("ETP", "Near-ETP", "Non-ETP"), ]
dd$subtype <- factor(dd$subtype)
dd$etp <- factor(dd$etp)
a0 <- anova(lm(balance ~ subtype, dd))
a1 <- anova(lm(balance ~ subtype + etp, dd))
etp_var <- data.frame(
  r2_subtype = summary(lm(balance ~ subtype, dd))$r.squared,
  r2_subtype_etp = summary(lm(balance ~ subtype + etp, dd))$r.squared,
  p_etp_after_subtype = a1$`Pr(>F)`[2]
)

# subtype landscape
hedges_g <- function(a, b) {
  n1 <- length(a); n2 <- length(b)
  if (n1 < 2 || n2 < 2) return(NA_real_)
  sp <- sqrt(((n1 - 1) * var(a) + (n2 - 1) * var(b)) / (n1 + n2 - 2))
  if (!is.finite(sp) || sp == 0) return(NA_real_)
  j <- 1 - 3 / (4 * (n1 + n2) - 9)
  j * (mean(a) - mean(b)) / sp
}
metrics <- c("ZEB1", "ZEB2", "LMO2", "balance")
subs <- sort(unique(d$subtype))
land <- list()
for (metric in metrics) {
  for (s in subs) {
    a <- d[[metric]][d$subtype == s]
    b <- d[[metric]][d$subtype != s]
    wt <- t.test(a, b)
    land[[paste(metric, s)]] <- data.frame(
      metric = metric, subtype = s, n = length(a),
      median = median(a), mean = mean(a),
      hedges_g_vs_rest = hedges_g(a, b),
      diff = mean(a) - mean(b),
      ci_low = wt$conf.int[1], ci_high = wt$conf.int[2],
      p = wt$p.value, stringsAsFactors = FALSE
    )
  }
}
land <- rbindlist(land)
land[, p_bh := p.adjust(p, "BH"), by = metric]

# old four-group map, secondary only
map4 <- function(s) {
  if (s %in% c("TLX1", "TLX3")) return("TLX")
  if (s %in% c("TAL1 αβ-like", "TAL1 DP-like")) return("TAL")
  if (s %in% c("HOXA9 TCR", "KMT2A", "MLLT10")) return("HOXA")
  if (s == "ETP-like") return("ETP-like")
  "other"
}
d$group4 <- vapply(d$subtype, map4, "")
g4 <- list()
for (metric in metrics) {
  for (s in c("TLX", "TAL", "HOXA", "ETP-like")) {
    a <- d[[metric]][d$group4 == s]
    b <- d[[metric]][d$group4 != s]
    wt <- t.test(a, b)
    g4[[paste(metric, s)]] <- data.frame(
      metric = metric, group = s, n = length(a), median = median(a),
      hedges_g_vs_rest = hedges_g(a, b),
      ci_low = wt$conf.int[1], ci_high = wt$conf.int[2], p = wt$p.value
    )
  }
}
g4 <- rbindlist(g4)
g4[, p_bh := p.adjust(p, "BH"), by = metric]

# events by subtype
ev <- as.data.table(d)[, .(
  n = .N,
  induction_failure = sum(if_author),
  mrd_ge_0.1 = sum(mrd_pos_0.1 == 1, na.rm = TRUE),
  mrd_known = sum(!is.na(mrd)),
  mrd_ge_0.01 = sum(mrd_pos_0.01 == 1, na.rm = TRUE),
  median_balance = median(balance),
  etp_n = sum(etp == "ETP", na.rm = TRUE)
), by = subtype]

fwrite(rows, file.path(out_dir, "model_induction_failure.tsv"), sep = "\t")
fwrite(rbindlist(mrd_rows), file.path(out_dir, "model_mrd.tsv"), sep = "\t")
fwrite(poor, file.path(out_dir, "model_poor_morphology.tsv"), sep = "\t")
fwrite(surv_tab, file.path(out_dir, "model_survival.tsv"), sep = "\t")
fwrite(as.data.frame(etp_contrast), file.path(out_dir, "etp_balance_contrast.tsv"), sep = "\t")
fwrite(within_etp, file.path(out_dir, "etp_within_failure.tsv"), sep = "\t")
fwrite(etp_var, file.path(out_dir, "etp_variance_after_subtype.tsv"), sep = "\t")
fwrite(land, file.path(out_dir, "subtype_landscape.tsv"), sep = "\t")
fwrite(g4, file.path(out_dir, "group4_correspondence.tsv"), sep = "\t")
fwrite(ev, file.path(out_dir, "events_by_subtype.tsv"), sep = "\t")

ref <- data.frame(
  gene = c("ZEB1", "ZEB2"),
  mean_logcpm = c(mean(d$ZEB1), mean(d$ZEB2)),
  sd_logcpm = c(sd(d$ZEB1), sd(d$ZEB2)),
  spearman_balance_tmm_vs_author_gexp = cor_note,
  mrd_max_stored = mrd_max,
  mrd_unit = "stored values match percent; author outcome cut MRD>=0.1; author residual-disease cut MRD>=0.01",
  endpoint = "author Risk.Group Induction Failure, identical patients to Day 29 M3 if counts match",
  if_and_m3 = sum(d$if_author == 1 & d$morph == "M3"),
  wbc_excluded_gt_2000 = sum(d$wbc > 2000, na.rm = TRUE),
  spearman_balance_mrd_rho = unname(sp$estimate),
  spearman_balance_mrd_p = sp$p.value
)
fwrite(ref, file.path(out_dir, "frozen_definition_check.tsv"), sep = "\t")
message("done")
print(rows[, c("model", "or", "ci_low", "ci_high", "p", "n", "events")])

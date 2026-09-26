# Scope-limited clinical sensitivity. These models do not define a biomarker.
suppressPackageStartupMessages({library(readr); library(logistf); library(survival)})
argv <- grep("^--file=", commandArgs(FALSE), value = TRUE)
script <- normalizePath(sub("^--file=", "", argv), winslash = "/")
base <- normalizePath(file.path(dirname(script), ".."), winslash = "/")
input <- file.path(base, "../data/validation/polonen_round1/patient_level_frozen_balance.tsv")
out <- file.path(base, "data/source_data_rebuilt/revision2_statistics")
d <- as.data.frame(read_tsv(input, locale = locale(encoding = "UTF-8"),
                            show_col_types = FALSE, progress = FALSE,
                            name_repair = "minimal"))
stopifnot(nrow(d) == 1309L, length(unique(d$sample_id)) == 1309L)
d$subtype <- factor(d$subtype)
d$sex <- factor(d$sex)

models <- list(
  list(endpoint="Induction failure", outcome="if_author", kind="binary"),
  list(endpoint="M2/M3 morphology", outcome="poor_morph", kind="binary"),
  list(endpoint="MRD >=0.1%", outcome="mrd_pos_0.1", kind="binary"),
  list(endpoint="MRD >=0.01%", outcome="mrd_pos_0.01", kind="binary"),
  list(endpoint="Event-free survival", outcome=c("efs_time","efs_status"), kind="survival"),
  list(endpoint="Overall survival", outcome=c("os_time","os_status"), kind="survival"))
results <- list()
for (spec in models) {
  needed <- c(spec$outcome, "balance", "subtype", "age", "sex", "log10_wbc")
  x <- d[complete.cases(d[, needed]), needed, drop = FALSE]
  if (spec$kind == "binary") {
    formula <- as.formula(paste0("`", spec$outcome, "` ~ balance + subtype + age + sex + log10_wbc"))
    fit <- tryCatch(logistf(formula, data=x, pl=FALSE), error=function(e) e)
    if (inherits(fit,"error")) {
      results[[length(results)+1L]] <- data.frame(endpoint=spec$endpoint,
        method="Firth logistic", n=nrow(x), events=sum(x[[spec$outcome]]==1),
        estimate=NA_real_, ci_low=NA_real_, ci_high=NA_real_, p=NA_real_,
        ph_balance_p=NA_real_, ph_global_p=NA_real_, status=conditionMessage(fit))
    } else {
      results[[length(results)+1L]] <- data.frame(endpoint=spec$endpoint,
        method="Firth logistic", n=nrow(x), events=sum(x[[spec$outcome]]==1),
        estimate=exp(fit$coefficients["balance"]),
        ci_low=exp(fit$ci.lower["balance"]),
        ci_high=exp(fit$ci.upper["balance"]), p=fit$prob["balance"],
        ph_balance_p=NA_real_, ph_global_p=NA_real_, status="converged")
    }
  } else {
    formula <- as.formula(paste0("Surv(",spec$outcome[1],",",spec$outcome[2],
      ") ~ balance + subtype + age + sex + log10_wbc"))
    warning_messages <- character()
    fit <- withCallingHandlers(
      tryCatch(coxph(formula,data=x,ties="efron",x=TRUE), error=function(e) e),
      warning=function(w) {
        warning_messages <<- c(warning_messages, conditionMessage(w))
        invokeRestart("muffleWarning")
      })
    if (inherits(fit,"error")) {
      results[[length(results)+1L]] <- data.frame(endpoint=spec$endpoint,
        method="Cox PH", n=nrow(x), events=sum(x[[spec$outcome[2]]]==1),
        estimate=NA_real_,ci_low=NA_real_,ci_high=NA_real_,p=NA_real_,
        ph_balance_p=NA_real_,ph_global_p=NA_real_,status=conditionMessage(fit))
    } else {
      s <- summary(fit)
      z <- cox.zph(fit, transform="km")$table
      results[[length(results)+1L]] <- data.frame(endpoint=spec$endpoint,
        method="Cox PH", n=nrow(x), events=sum(x[[spec$outcome[2]]]==1),
        estimate=exp(coef(fit)["balance"]),
        ci_low=s$conf.int["balance","lower .95"],
        ci_high=s$conf.int["balance","upper .95"],
        p=s$coefficients["balance","Pr(>|z|)"],
        ph_balance_p=z["balance","p"], ph_global_p=z["GLOBAL","p"],
        status=if (length(warning_messages)) paste(unique(warning_messages),collapse=" | ") else "converged")
    }
  }
}
ans <- do.call(rbind, results)
ans$endpoint_bh_fdr <- p.adjust(ans$p, method="BH")
readr::write_tsv(ans, file.path(out,"clinical_age_sex_wbc_sensitivity.tsv"), na="NA")
writeLines(c(
  "Exploratory clinical sensitivity using the frozen 1309-patient table.",
  "Each model includes balance, all 17 molecular subtypes, age, sex and log10(WBC).",
  "Binary endpoints use Firth logistic regression; survival endpoints use Efron-ties Cox regression.",
  "Complete cases are used per endpoint and counts are reported. No cutoff or new score is fitted.",
  "Cox proportional-hazards tests use cox.zph with KM time transform; p values test balance and the global model, respectively.",
  "Cox rare-subtype coefficient warnings are retained in the status column; finite balance estimates do not remove this caution.",
  "BH FDR is calculated across these six sensitivity endpoints.",
  "These sensitivity estimates are descriptive and do not establish independent prognostic utility."),
  file.path(out,"CLINICAL_METHODS.txt"))
writeLines(capture.output(sessionInfo()), file.path(out,"clinical_R_sessionInfo.txt"))
print(ans)

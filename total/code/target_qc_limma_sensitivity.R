#!/usr/bin/env Rscript
# Diagnostic sensitivity: does broad ZEB1 association survive STAR QC proxies?
# These QC fields may carry biology and are not prespecified confounders.
suppressPackageStartupMessages({
  library(data.table)
  library(limma)
})
root <- "/data-b/liangfuhua/projects/ZEB1_TALL_analysis/module4"
proc <- file.path(root, "processed")
out <- file.path(root, "tables", "target_star_raw_qc")
dir.create(out, recursive=TRUE, showWarnings=FALSE)
read_mat <- function(path) {
  d <- as.data.frame(fread(path), stringsAsFactors=FALSE)
  rownames(d) <- d[[1]]; d[[1]] <- NULL
  as.matrix(d)
}
expr <- read_mat(file.path(proc, "target_log2tpm_filtered.tsv.gz"))
lab <- fread(file.path(proc, "target_samples.tsv"))
qc <- fread(file.path(out, "target_star_raw_sample_qc.tsv"))
stopifnot(!anyDuplicated(lab$usi), !anyDuplicated(qc$usi))
lab <- merge(lab, qc[, .(usi, assigned_fraction_all_reported,
                          N_noFeature, mitochondrial_count_fraction)],
             by="usi", all.x=TRUE, sort=FALSE)
lab <- lab[usi %in% colnames(expr) & subtype != "Unknown" & !is.na(subtype)]
lab <- lab[match(colnames(expr), lab$usi)]
lab <- lab[!is.na(subtype) & subtype != "Unknown" & is.finite(age_years) &
           is.finite(assigned_fraction_all_reported) &
           is.finite(N_noFeature) & is.finite(mitochondrial_count_fraction)]
stopifnot(nrow(lab) == 242L, !anyDuplicated(lab$usi))
expr <- expr[, lab$usi, drop=FALSE]
stopifnot(identical(colnames(expr), lab$usi))
lab[, ZEB1_z := as.numeric(scale(log2(ZEB1+1)))]
lab[, assigned_z := as.numeric(scale(assigned_fraction_all_reported))]
lab[, nofeature_z := as.numeric(scale(log1p(N_noFeature)))]
lab[, mito_z := as.numeric(scale(mitochondrial_count_fraction))]
lab[, subtype := factor(subtype)]
lab[, sample_type := factor(sample_type)]

forms <- list(
  original = ~ ZEB1_z + subtype + age_years,
  source = ~ ZEB1_z + subtype + age_years + sample_type,
  assigned = ~ ZEB1_z + subtype + age_years + sample_type + assigned_z,
  nofeature = ~ ZEB1_z + subtype + age_years + sample_type + nofeature_z,
  mito = ~ ZEB1_z + subtype + age_years + sample_type + mito_z
)
summary <- list(); key <- list()
focus <- c("GATA3", "TCF7", "BCL11B", "IL7R", "LMO2", "LYL1", "ZEB2", "CD34")
for (name in names(forms)) {
  design <- model.matrix(forms[[name]], data=lab)
  if (qr(design)$rank != ncol(design)) stop("rank-deficient design: ", name)
  fit <- eBayes(lmFit(expr, design), trend=TRUE, robust=TRUE)
  tt <- as.data.table(topTable(fit, coef="ZEB1_z", number=Inf, sort.by="none"),
                      keep.rownames="gene")
  tt[, model := name]
  setcolorder(tt, c("model", "gene"))
  fwrite(tt, file.path(out, paste0("zeb1_genomewide_", name, ".tsv")), sep="\t")
  summary[[name]] <- data.table(model=name, n=nrow(lab), genes=nrow(tt),
                                fdr05=sum(tt$adj.P.Val < .05, na.rm=TRUE),
                                positive_coef=sum(tt$logFC > 0, na.rm=TRUE),
                                formula=paste(deparse(forms[[name]]), collapse=" "))
  key[[name]] <- tt[gene %in% focus, .(model, gene, logFC, t, P.Value, adj.P.Val)]
  message("done ", name, ": FDR<.05 ", summary[[name]]$fdr05)
}
fwrite(rbindlist(summary), file.path(out, "qc_sensitivity_summary.tsv"), sep="\t")
fwrite(rbindlist(key), file.path(out, "qc_sensitivity_keygenes.tsv"), sep="\t")
writeLines(c(capture.output(sessionInfo()), "", paste("samples:", nrow(lab))),
           file.path(out, "qc_sensitivity_sessionInfo.txt"))

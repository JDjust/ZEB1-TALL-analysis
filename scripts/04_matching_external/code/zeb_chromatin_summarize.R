suppressPackageStartupMessages(library(data.table))
out <- "D:/_bioinformation/ZEB1/data/analysis/total/data/validation/zeb_chromatin_reciprocal"
new <- fread(file.path(out, "hichip_zeb1_zeb2_contacts.tsv"))
old <- fread(file.path(out, "M6_6.13_6.16_ZEB1_contacts.tsv"))
atac <- fread(file.path(out, "scatac_zeb1_zeb2_peaks.tsv"))

chk <- merge(new[, .(sample, ZEB1, ZEB1_per_M)], old[, .(sample, ZEB1_old = ZEB1_locus, ZEB1_per_M_old = ZEB1_per_M)], by = "sample")
chk[, abs_diff := ZEB1_old - ZEB1]
chk[, pct := 100 * abs_diff / ZEB1_old]
fwrite(chk, file.path(out, "zeb1_count_crosscheck.tsv"), sep = "\t")

use <- new[group %in% c("non-ETP T-ALL", "ETP", "CD34", "thymus")]
summ <- use[, .(
  n = .N,
  ZEB1_per_M_median = median(ZEB1_per_M),
  ZEB1_per_M_min = min(ZEB1_per_M),
  ZEB1_per_M_max = max(ZEB1_per_M),
  ZEB2_per_M_median = median(ZEB2_per_M),
  ZEB1_per_Mkb_median = median(ZEB1_per_M_per_kb),
  ZEB2_per_Mkb_median = median(ZEB2_per_M_per_kb),
  ZEB2_over_ZEB1_median = median(ZEB2_over_ZEB1)
), by = group]
fwrite(summ, file.path(out, "hichip_group_summary.tsv"), sep = "\t")

mw <- function(a, b) {
  if (length(a) < 2 || length(b) < 2) return(list(p = NA_real_, n_a = length(a), n_b = length(b)))
  t <- wilcox.test(a, b, exact = TRUE)
  list(p = unname(t$p.value), n_a = length(a), n_b = length(b))
}
etp <- use[group == "ETP"]
tall <- use[group == "non-ETP T-ALL"]
tests <- rbindlist(list(
  data.table(metric = "ZEB1_per_M", contrast = "ETP vs non-ETP", etp_median = median(etp$ZEB1_per_M), tall_median = median(tall$ZEB1_per_M), p = mw(etp$ZEB1_per_M, tall$ZEB1_per_M)$p),
  data.table(metric = "ZEB2_per_M", contrast = "ETP vs non-ETP", etp_median = median(etp$ZEB2_per_M), tall_median = median(tall$ZEB2_per_M), p = mw(etp$ZEB2_per_M, tall$ZEB2_per_M)$p),
  data.table(metric = "ZEB1_per_M_per_kb", contrast = "ETP vs non-ETP", etp_median = median(etp$ZEB1_per_M_per_kb), tall_median = median(tall$ZEB1_per_M_per_kb), p = mw(etp$ZEB1_per_M_per_kb, tall$ZEB1_per_M_per_kb)$p),
  data.table(metric = "ZEB2_per_M_per_kb", contrast = "ETP vs non-ETP", etp_median = median(etp$ZEB2_per_M_per_kb), tall_median = median(tall$ZEB2_per_M_per_kb), p = mw(etp$ZEB2_per_M_per_kb, tall$ZEB2_per_M_per_kb)$p),
  data.table(metric = "ZEB2_over_ZEB1", contrast = "ETP vs non-ETP", etp_median = median(etp$ZEB2_over_ZEB1), tall_median = median(tall$ZEB2_over_ZEB1), p = mw(etp$ZEB2_over_ZEB1, tall$ZEB2_over_ZEB1)$p)
))
fwrite(tests, file.path(out, "hichip_etp_vs_nontall.tsv"), sep = "\t")

atac_sum <- atac[, .(
  n = .N,
  ZEB1_peaks_median = median(ZEB1_peaks),
  ZEB2_peaks_median = median(ZEB2_peaks),
  ZEB1_per_10k_median = median(ZEB1_per_10k_peaks),
  ZEB2_per_10k_median = median(ZEB2_per_10k_peaks)
), by = group]
fwrite(atac_sum, file.path(out, "scatac_group_summary.tsv"), sep = "\t")
a_etp <- atac[group == "ETP"]
a_tall <- atac[group == "non-ETP T-ALL"]
atac_tests <- rbindlist(list(
  data.table(metric = "ZEB1_per_10k", p = mw(a_etp$ZEB1_per_10k_peaks, a_tall$ZEB1_per_10k_peaks)$p, etp_n = nrow(a_etp), tall_n = nrow(a_tall)),
  data.table(metric = "ZEB2_per_10k", p = mw(a_etp$ZEB2_per_10k_peaks, a_tall$ZEB2_per_10k_peaks)$p, etp_n = nrow(a_etp), tall_n = nrow(a_tall))
))
fwrite(atac_tests, file.path(out, "scatac_etp_vs_nontall.tsv"), sep = "\t")

cat("ZEB1 crosscheck vs prior table:\n"); print(chk)
cat("\nHiChIP group:\n"); print(summ)
cat("\nHiChIP tests:\n"); print(tests)
cat("\nscATAC group:\n"); print(atac_sum)
cat("\nscATAC tests:\n"); print(atac_tests)

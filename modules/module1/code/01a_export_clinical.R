# Robust TARGET clinical merge + pharmacotype lineage counts
suppressPackageStartupMessages({
  library(readxl)
  library(data.table)
})
clin_dir <- "/data-b/liangfuhua/projects/TALL_dataset_download/data/TARGET-ALL-P2/Clinical_Supplement/na"
out <- "/data-b/liangfuhua/projects/ZEB1_TALL_analysis/module1/processed"
dir.create(out, recursive=TRUE, showWarnings=FALSE)

clean_names <- function(x) {
  x <- gsub("[\r\n]+", " ", x)
  x <- gsub("[[:space:]]+", " ", x)
  trimws(x)
}

files <- list.files(clin_dir, pattern="ClinicalData_.*\\.xlsx$", full.names=TRUE)
lst <- lapply(files, function(f) {
  d <- as.data.frame(read_excel(f, sheet=1))
  names(d) <- clean_names(names(d))
  d$source_file <- basename(f)
  as.data.table(d)
})
d <- rbindlist(lst, fill=TRUE)
setnames(d, names(d), clean_names(names(d)))
cat("combined", nrow(d), "x", ncol(d), "\n")
cat("cols:", paste(names(d), collapse=" | "), "\n")
cat("Protocol:\n"); print(table(as.character(d$Protocol), useNA="ifany"))
if ("Cell of Origin" %in% names(d)) {
  cat("Cell of Origin:\n"); print(table(as.character(d[["Cell of Origin"]]), useNA="ifany"))
}
cat("unique USI", uniqueN(d[["TARGET USI"]]), "\n")
# prefer validation then discovery if duplicated USI
pref <- c(
  "TARGET_ALL_ClinicalData_Phase_II_Validation_20230727.xlsx",
  "TARGET_ALL_ClinicalData_Phase_II_Discovery_20230727.xlsx",
  "TARGET_ALL_ClinicalData_Phase_I_20230727.xlsx",
  "TARGET_ALL_ClinicalData_Dicentric_20230727.xlsx",
  "TARGET_ALL_ClinicalData_Xenografts_20230727.xlsx"
)
src_col <- if ("source_file" %in% names(d)) "source_file" else grep("source", names(d), value=TRUE, ignore.case=TRUE)[1]
d[, pref_rank := match(get(src_col), pref)]
d[is.na(pref_rank), pref_rank := 99]
setorder(d, pref_rank)
d1 <- d[!duplicated(`TARGET USI`)]
cat("dedup USI", nrow(d1), "\n")
fwrite(d, file.path(out, "target_clinical_combined_allrows.tsv"), sep="\t")
fwrite(d1, file.path(out, "target_clinical_combined.tsv"), sep="\t")

sup <- as.data.frame(read_excel(file.path(clin_dir, "TARGET_ALL_ClinicalData_Validation_Supplement_20230727.xlsx")))
names(sup) <- clean_names(names(sup))
fwrite(as.data.table(sup), file.path(out, "target_t_all_supplement.tsv"), sep="\t")

ph <- "/data-b/liangfuhua/projects/TALL_dataset_download/data/StJude_ALL_Pharmacotype/41591_2022_2112_MOESM3_ESM.xlsx"
s1 <- as.data.table(read_excel(ph, sheet=1, skip=1))
s2 <- as.data.table(read_excel(ph, sheet=2, skip=1))
cat("Immunophenotype:\n"); print(table(as.character(s1$Immunophenotype), useNA="ifany"))
cat("Molecular subtype top:\n"); print(head(sort(table(as.character(s1[["Molecular subtype"]]), useNA="ifany"), decreasing=TRUE), 20))
fwrite(s1, file.path(out, "pharmacotype_supp_table1.tsv"), sep="\t")
fwrite(s2, file.path(out, "pharmacotype_supp_table2.tsv"), sep="\t")
cat("s2 sample id examples:", paste(head(s2[["Sample ID"]], 8), collapse=", "), "\n")
cat("s2 patient id examples:", paste(head(s1[["Patient ID"]], 8), collapse=", "), "\n")

# Read-only extract of curated GNAB flags vs ZEB1. Does not modify RData.
fm_path <- "D:/_bioinformation/ZEB1/.project-artifacts/work/20261009T052031Z-04811c83/BCL11B_QUALIFICATION/TALL_X01_FeatureMatrix_genomics_subtype.Rdata"
clin_path <- "D:/_bioinformation/ZEB1/02_analyses/03_q1_bioinfo_paper/cox_analysis_table.tsv"
out_dir <- "D:/_bioinformation/ZEB1/02_analyses/03_q1_bioinfo_paper"
genes <- c(
  "NOTCH1","FBXW7","PTEN","PHF6","WT1","FLT3","NRAS","KRAS","NF1","PTPN11",
  "IL7R","JAK1","JAK3","STAT5B","PIK3R1","BCL11B","EZH2","SUZ12","RUNX1",
  "GATA3","DNM2","LEF1","CDKN2A","TAL1","LMO2","TLX1","TLX3","MYB","USP7",
  "KDM6A","CTCF","RPL5","DDX3X","BCOR","ETV6","CNOT3","MED12"
)
e <- new.env(parent = emptyenv())
load(fm_path, envir = e)
fm <- e$fm
rows <- rownames(fm)
pick <- paste0("B:GNAB:", genes, ":Significant_genes")
have <- pick[pick %in% rows]
miss <- genes[!paste0("B:GNAB:", genes, ":Significant_genes") %in% rows]
cat("found", length(have), "missing", paste(miss, collapse=", "), "\n")
mat <- t(as.matrix(fm[have, , drop = FALSE]))
colnames(mat) <- sub("^B:GNAB:", "", sub(":Significant_genes$", "", colnames(mat)))
mat <- as.data.frame(mat, check.names = FALSE)
mat$sample_id <- rownames(mat)
clin <- read.delim(clin_path, check.names = FALSE, stringsAsFactors = FALSE)
df <- merge(clin[, c("sample_id","ZEB1","etp","subtype","zeb1_tertile","os_status","efs_status")], mat, by = "sample_id")
write.table(df, file.path(out_dir, "mutation_join.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
cat("n join", nrow(df), "\n")

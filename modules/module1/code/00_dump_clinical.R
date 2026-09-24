# dump TARGET clinical column names and lineage-like fields
clin_dir <- "/data-b/liangfuhua/projects/TALL_dataset_download/data/TARGET-ALL-P2/Clinical_Supplement/na"
out <- "/data-b/liangfuhua/projects/ZEB1_TALL_analysis/module1/inspect"
dir.create(out, recursive=TRUE, showWarnings=FALSE)
files <- list.files(clin_dir, pattern="\\.xlsx$", full.names=TRUE)
for (f in files) {
  cat("\n====", basename(f), "====\n")
  d <- tryCatch(readxl::read_excel(f, sheet=1), error=function(e) { cat("ERR", e$message, "\n"); NULL })
  if (is.null(d)) next
  cat("dim", paste(dim(d), collapse=" x "), "\n")
  cat("cols:", paste(colnames(d), collapse=" | "), "\n")
  hit <- grep("lineage|immuno|phenotype|diagnos|subtype|T-ALL|B-ALL|ALL Type|WHO|disease|protocol", colnames(d), ignore.case=TRUE, value=TRUE)
  cat("hit:", paste(hit, collapse=" | "), "\n")
  for (h in hit) {
    tb <- sort(table(as.character(d[[h]])), decreasing=TRUE)
    cat("  ", h, "\n")
    print(utils::head(tb, 15))
  }
  # also dump first 2 rows of key ID columns
  idhit <- grep("TARGET|USI|case|submitter|barcode|patient", colnames(d), ignore.case=TRUE, value=TRUE)
  cat("ids:", paste(idhit, collapse=" | "), "\n")
}

# pharmacotype excel
ph <- "/data-b/liangfuhua/projects/TALL_dataset_download/data/StJude_ALL_Pharmacotype/41591_2022_2112_MOESM3_ESM.xlsx"
cat("\n==== PHARM XLSX ====\n")
shts <- readxl::excel_sheets(ph)
cat("sheets", paste(shts, collapse=" | "), "\n")
for (s in shts) {
  d <- readxl::read_excel(ph, sheet=s, n_max=5)
  cat("SHEET", s, "cols", paste(colnames(d), collapse=" | "), "\n")
}

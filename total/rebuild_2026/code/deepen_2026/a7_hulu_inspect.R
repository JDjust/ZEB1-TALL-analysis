# A7-A step 1 on hulu: annotation / assay only. No side50. No scores.
suppressPackageStartupMessages(library(Seurat))

dir <- "/data-b/liangfuhua/zeb1_a7"
p <- file.path(dir, "GSE253355_Normal_Bone_Marrow_Atlas_Seurat_SB_v2.rds.gz")
out <- file.path(dir, "inspect")
dir.create(out, recursive = TRUE, showWarnings = FALSE)

cat("file_bytes", file.info(p)$size, "\n"); flush.console()
cat("readRDS start (double gzip)\n"); flush.console()
obj <- readRDS(gzcon(gzfile(p, "rb")))
cat("class", paste(class(obj), collapse = "/"), "ncells", ncol(obj), "nfeatures", nrow(obj), "\n")
flush.console()

sink(file.path(out, "a7_object_structure.txt"))
cat("class:\n"); print(class(obj))
cat("\nAssays:\n"); print(Assays(obj))
cat("\nDefaultAssay:\n"); print(DefaultAssay(obj))
cat("\nReductions:\n"); print(Reductions(obj))
cat("\nncells nfeatures:\n")
cat(ncol(obj), nrow(obj), "\n")
for (a in Assays(obj)) {
  cat("\n==== assay", a, "====\n")
  asy <- tryCatch(obj[[a]], error = function(e) { cat("assay_error", conditionMessage(e), "\n"); NULL })
  if (is.null(asy)) next
  cat("class", paste(class(asy), collapse = "/"), "\n")
  tryCatch({
    if (inherits(asy, "Assay5") && exists("Layers", mode = "function")) {
      cat("layers:\n"); print(Layers(asy))
    } else {
      for (s in intersect(c("counts", "data", "scale.data"), slotNames(asy))) {
        m <- tryCatch(slot(asy, s), error = function(e) NULL)
        if (!is.null(m)) cat("slot", s, "dim", paste(dim(m), collapse = "x"), "\n")
      }
    }
  }, error = function(e) cat("assay_detail_error", conditionMessage(e), "\n"))
}
cat("\nmeta columns:\n")
print(colnames(obj[[]]))
sink()

md <- obj[[]]
write.table(
  data.frame(
    column = colnames(md),
    class = vapply(md, function(x) paste(class(x), collapse = "/"), ""),
    n_unique = vapply(md, function(x) length(unique(x)), 1L)
  ),
  file.path(out, "a7_meta_columns.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

for (k in colnames(md)) {
  x <- md[[k]]
  nu <- length(unique(x))
  ok <- nu <= 200 && (is.character(x) || is.factor(x) || is.logical(x) || (is.numeric(x) && nu <= 50))
  if (!ok) next
  tab <- sort(table(as.character(x), useNA = "ifany"), decreasing = TRUE)
  write.table(
    data.frame(level = names(tab), n = as.integer(tab)),
    file.path(out, paste0("a7_meta_", gsub("[^A-Za-z0-9]+", "_", k), ".tsv")),
    sep = "\t", quote = FALSE, row.names = FALSE
  )
}

saveRDS(md, file.path(out, "a7_meta_only.rds"))
cat("wrote", out, "\n")

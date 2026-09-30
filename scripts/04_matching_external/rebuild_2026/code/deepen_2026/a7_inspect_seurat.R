# A7-A step 1: annotation / assay only. No side50. No scores.
suppressPackageStartupMessages({
  library(Seurat)
})

p <- "D:/_bioinformation/ZEB1/data/analysis/total/rebuild_2026/data/deepen_2026/a7_lineage/GSE253355_Normal_Bone_Marrow_Atlas_Seurat_SB_v2.rds.gz"
out <- "D:/_bioinformation/ZEB1/data/analysis/total/rebuild_2026/data/deepen_2026/a7_lineage"

cat("readRDS start (double gzip)\n"); flush.console()
obj <- readRDS(gzcon(gzfile(p, "rb")))
cat("class", paste(class(obj), collapse = "/"), "\n"); flush.console()

sink(file.path(out, "a7_object_structure.txt"))
cat("class:\n"); print(class(obj))
cat("\nAssays:\n"); print(Assays(obj))
cat("\nDefaultAssay:\n"); print(DefaultAssay(obj))
cat("\nReductions:\n"); print(Reductions(obj))
cat("\nncells nfeatures:\n")
cat(ncol(obj), nrow(obj), "\n")
for (a in Assays(obj)) {
  cat("\n==== assay", a, "====\n")
  asy <- obj[[a]]
  cat("class", paste(class(asy), collapse = "/"), "\n")
  if (inherits(asy, "Assay5")) {
    cat("layers:\n"); print(Layers(asy))
    for (ly in Layers(asy)) {
      cat("layer", ly, "dim", paste(dim(asy[layer = ly]), collapse = "x"), "\n")
    }
  } else {
    slots <- intersect(c("counts", "data", "scale.data"), slotNames(asy))
    for (s in slots) {
      m <- slot(asy, s)
      cat("slot", s, "dim", paste(dim(m), collapse = "x"), "nnz_guess_ok\n")
    }
  }
}
cat("\nmeta columns:\n")
print(colnames(obj[[]]))
sink()

md <- obj[[]]
write.table(
  data.frame(column = colnames(md), class = vapply(md, function(x) paste(class(x), collapse = "/"), ""),
             n_unique = vapply(md, function(x) length(unique(x)), 1L)),
  file.path(out, "a7_meta_columns.tsv"), sep = "\t", quote = FALSE, row.names = FALSE
)

# write n-per-level for likely annotation / donor columns
keys <- colnames(md)
want <- keys[grepl("donor|sample|patient|orig|ident|cluster|anno|cell|type|lineage|library|nCount|nFeature|percent",
                   keys, ignore.case = TRUE)]
if (!length(want)) want <- keys
writeLines(want, file.path(out, "a7_meta_priority_cols.txt"))

for (k in keys) {
  x <- md[[k]]
  nu <- length(unique(x))
  if (nu <= 200 && (is.character(x) || is.factor(x) || is.logical(x) || (is.numeric(x) && nu <= 50))) {
    tab <- sort(table(as.character(x), useNA = "ifany"), decreasing = TRUE)
    write.table(
      data.frame(level = names(tab), n = as.integer(tab)),
      file.path(out, paste0("a7_meta_", gsub("[^A-Za-z0-9]+", "_", k), ".tsv")),
      sep = "\t", quote = FALSE, row.names = FALSE
    )
  }
}

cat("wrote annotation tables\n")
cat("ncells", ncol(obj), "nfeatures", nrow(obj), "\n")

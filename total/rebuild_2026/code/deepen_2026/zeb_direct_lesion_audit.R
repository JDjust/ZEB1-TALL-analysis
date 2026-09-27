# Ask whether subtype-associated ZEB states are explained by recurrent
# direct ZEB1/ZEB2 genomic lesions in the Pölönen feature matrix.
suppressPackageStartupMessages({library(data.table)})
argv <- grep("^--file=", commandArgs(FALSE), value = TRUE)
script <- normalizePath(sub("^--file=", "", argv), winslash = "/")
rebuild <- normalizePath(file.path(dirname(script), "../.."), winslash = "/")
root <- normalizePath(file.path(rebuild, "../.."), winslash = "/")
out <- file.path(rebuild, "data/deepen_2026/zeb_lesion_audit")
dir.create(out, recursive = TRUE, showWarnings = FALSE)

pt <- fread(file.path(root, "total/data/validation/zeb_developmental_residual/polonen_patient_residual.tsv"),
            encoding = "UTF-8")
e <- new.env(parent = emptyenv())
load(file.path(root, "data/polonen_syn54032669/Data_1309Samples.RData"), envir = e)
objs <- ls(e)
fwrite(data.table(object = objs, class = vapply(objs, function(n) paste(class(e[[n]]), collapse = "/"), "")),
       file.path(out, "rdata_objects.tsv"), sep = "\t")

search_matrix <- function(x, name) {
  if (!is.matrix(x) && !is.data.frame(x)) return(NULL)
  rn <- rownames(x)
  if (is.null(rn)) return(NULL)
  hit <- grepl("ZEB1|ZEB2|ZEB-1|ZEB-2", rn, ignore.case = TRUE)
  if (!any(hit)) return(NULL)
  m <- as.matrix(x[hit, , drop = FALSE])
  data.table(source = name, feature = rownames(m), n_nonzero = rowSums(m != 0, na.rm = TRUE),
             n_na = rowSums(is.na(m)), n_samples = ncol(m))
}

rows <- rbindlist(lapply(objs, function(n) search_matrix(e[[n]], n)), fill = TRUE)
if (!nrow(rows)) {
  rows <- data.table(source = NA, feature = NA, n_nonzero = NA, n_na = NA, n_samples = NA,
                     note = "No rownames containing ZEB1/ZEB2 in loaded objects.")
}
fwrite(rows, file.path(out, "zeb_named_features.tsv"), sep = "\t")

# Also dump any character/factor columns that mention ZEB lesions.
text_hits <- list()
for (n in objs) {
  x <- e[[n]]
  if (is.data.frame(x)) {
    for (cn in names(x)) {
      if (!is.character(x[[cn]]) && !is.factor(x[[cn]])) next
      v <- as.character(x[[cn]])
      hit <- grepl("ZEB1|ZEB2", v, ignore.case = TRUE)
      if (any(hit, na.rm = TRUE)) {
        tab <- as.data.table(table(v[hit], useNA = "ifany"))
        setnames(tab, c("value", "n"))
        tab[, `:=`(object = n, column = cn)]
        text_hits[[length(text_hits) + 1L]] <- tab
      }
    }
  }
}
if (length(text_hits)) {
  fwrite(rbindlist(text_hits, fill = TRUE), file.path(out, "zeb_text_annotations.tsv"), sep = "\t")
}

# M5 lesion-gene matrix: count any ZEB-named rows and the most recurrent genes.
if (!is.null(e$M5.AllLesions.genes)) {
  m5 <- as.matrix(e$M5.AllLesions.genes)
  m5n <- data.table(feature = rownames(m5), n_patients = rowSums(m5 != 0, na.rm = TRUE))
  setorder(m5n, -n_patients)
  fwrite(m5n[n_patients > 0][1:100], file.path(out, "m5_top100_lesion_genes.tsv"), sep = "\t")
  fwrite(m5n[grepl("ZEB", feature, ignore.case = TRUE)],
         file.path(out, "m5_zeb_lesion_genes.tsv"), sep = "\t")
}
if (!is.null(e$M1.classifying.driver)) {
  d1 <- as.matrix(e$M1.classifying.driver)
  d1n <- data.table(feature = rownames(d1), n_patients = rowSums(d1 != 0, na.rm = TRUE))
  fwrite(d1n[grepl("ZEB", feature, ignore.case = TRUE)],
         file.path(out, "m1_zeb_classifying_driver.tsv"), sep = "\t")
}

writeLines(c(
  "This audit only asks whether ZEB1/ZEB2 appear as recurrent named genomic features.",
  "Absence of a named row is not proof that no cryptic structural variant exists.",
  "It is sufficient to say the deposited feature matrix does not annotate recurrent direct ZEB lesions.",
  paste("Patients in residual table:", nrow(pt))
), file.path(out, "METHODS.txt"))
cat("Lesion audit objects:", paste(objs, collapse = ", "), "\n")
print(rows)

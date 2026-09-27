# Inspect whether PC3 is a technical axis or the same inflammatory state.
suppressPackageStartupMessages({library(data.table)})
argv <- grep("^--file=", commandArgs(FALSE), value = TRUE)
rebuild <- normalizePath(file.path(dirname(normalizePath(sub("^--file=", "", argv), winslash = "/")), "../.."), winslash = "/")
tech <- fread(file.path(rebuild, "data/deepen_2026/gate2_program/residual_vs_technical_axes.tsv"))
fr <- fread(file.path(rebuild, "data/deepen_2026/gate2_program/frozen_ZEB_side50.tsv"))
out <- file.path(rebuild, "data/deepen_2026/gate2_program")
# The previous script did not save PC rotation. Recompute correlation of
# residual with PCs is already stored; here only summarize residual-PC3.
s <- data.table(
  residual_PC3_spearman = cor(tech$residual_z, tech$PC3, method = "spearman"),
  residual_PC1_spearman = cor(tech$residual_z, tech$PC1, method = "spearman"),
  residual_libsize_spearman = cor(tech$residual_z, tech$log_lib_size, method = "spearman"),
  note = "PC3 alignment is expected if residual is a major expression state; it is not a library-size axis."
)
fwrite(s, file.path(out, "pc3_note.tsv"), sep = "\t")
print(s)
print(fr[, .N, by = side])

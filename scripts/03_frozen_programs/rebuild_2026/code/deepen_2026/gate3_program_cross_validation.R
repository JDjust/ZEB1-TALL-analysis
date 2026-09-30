# Gate 3A/3B: frozen side50 programs on matched BCL11B-ETP and GSE146901.
# Lim 3C is skipped unless a genome-wide patient/state matrix is present.
suppressPackageStartupMessages({
  library(data.table); library(limma)
})
argv <- grep("^--file=", commandArgs(FALSE), value = TRUE)
rebuild <- normalizePath(file.path(dirname(normalizePath(sub("^--file=", "", argv), winslash = "/")), "../.."), winslash = "/")
root <- normalizePath(file.path(rebuild, "../.."), winslash = "/")
prog_dir <- file.path(rebuild, "data/deepen_2026/gate2_program")
out <- file.path(rebuild, "data/deepen_2026/gate3_validation")
dir.create(out, recursive = TRUE, showWarnings = FALSE)

fr <- fread(file.path(prog_dir, "frozen_ZEB_side50.tsv"))
stopifnot(nrow(fr[side == "ZEB1-side50"]) == 50L, nrow(fr[side == "ZEB2-side50"]) == 50L)
sets <- split(fr$symbol, fr$side)

# 3A: rank enrichment on the already pair-blocked matched contrast.
deg <- fread(file.path(rebuild, "data/source_data_rebuilt/revision4_targeted/matched_BCL11B_ETP_gene_results.tsv"))
stat <- deg$t
names(stat) <- deg$symbol
stat <- stat[is.finite(stat) & !duplicated(names(stat))]
idx <- lapply(sets, function(s) which(names(stat) %in% unique(s)))
cam <- cameraPR(stat, index = idx)
cam$program <- rownames(cam)
cam <- as.data.table(cam)[, .(program, NGenes, Direction, PValue, FDR)]
# Direction is relative to BCL11B-minus-ETP. ZEB2-side should be Up in BCL11B.
cam[, expected := fifelse(program == "ZEB2-side50", "Up", "Down")]
cam[, direction_matches := Direction == expected]
fwrite(cam, file.path(out, "matched_cameraPR_frozen_programs.tsv"), sep = "\t")

# 3B is scored in gate3_gse146901_program_scores.py because the GEO
# workbook is a gzipped xls that readxl cannot open.

lim_note <- data.table(
  gate = "3C",
  status = "blocked",
  reason = "Lim genome-wide ZBTB16-state matrix is not on this machine; stored tables have only ZEB1/ZEB2/LMO2/ZBTB16."
)
fwrite(lim_note, file.path(out, "lim_status.tsv"), sep = "\t")

writeLines(c(
  "3A uses the frozen 18-pair matched t-statistics and CAMERA PR. No threshold was retuned.",
  "ZEB2-side50 is expected Up in BCL11B minus ETP-like; ZEB1-side50 is expected Down.",
  "3B scores each program as the mean within-leukemia z of available FPKM genes.",
  "ETP is expected to have a higher ZEB2-side50 score and a lower ZEB1-side50 score.",
  "3C is not run without a genome-wide Lim state matrix.",
  "OXPHOS is not tested."
), file.path(out, "METHODS.txt"))
cat("Gate 3A complete\n")
print(cam)

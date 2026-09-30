# A6-3: cameraPR on frozen matched BCL11B vs ETP t-statistics. No new DE.
suppressPackageStartupMessages(library(limma))
out <- "D:/_bioinformation/ZEB1/data/analysis/total/rebuild_2026/data/deepen_2026/a6_yayon"
deg <- read.delim(
  "D:/_bioinformation/ZEB1/data/analysis/total/rebuild_2026/data/source_data_rebuilt/revision4_targeted/matched_BCL11B_ETP_gene_results.tsv",
  check.names = FALSE
)
frz <- read.delim(
  "D:/_bioinformation/ZEB1/data/analysis/total/rebuild_2026/data/deepen_2026/gate2_program/frozen_ZEB_side50.tsv",
  check.names = FALSE
)
conc <- read.delim(file.path(out, "a6_3_zeb2_early_concordant.tsv"), check.names = FALSE)
disc <- read.delim(file.path(out, "a6_3_zeb2_discordant_or_neutral.tsv"), check.names = FALSE)
tstat <- deg$t
names(tstat) <- deg$symbol
idx <- list(
  `ZEB2-side50` = intersect(frz$symbol[frz$side == "ZEB2-side50"], names(tstat)),
  `ZEB2-side50_early_concordant` = intersect(conc$symbol, names(tstat)),
  `ZEB2-side50_discordant_or_neutral` = intersect(disc$symbol, names(tstat)),
  `ZEB1-side50` = intersect(frz$symbol[frz$side == "ZEB1-side50"], names(tstat))
)
cam <- cameraPR(tstat, idx)
cam$program <- rownames(cam)
cam$n_in_index <- vapply(idx[cam$program], length, 1L)
write.table(cam, file.path(out, "a6_3_cameraPR.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
print(cam)

# 1.17 substitute: adult (GSE142522) vs pediatric (TARGET) ZEB1, different platforms.
suppressPackageStartupMessages({ library(data.table); library(ggplot2) })
m1 <- "/data-b/liangfuhua/projects/ZEB1_TALL_analysis/module1"
m2 <- "/data-b/liangfuhua/projects/ZEB1_TALL_analysis/module2"
m3 <- "/data-b/liangfuhua/projects/ZEB1_TALL_analysis/module3"
fig <- file.path(m1, "figures")
theme_sci <- function() {
  theme_classic(base_size = 11) +
    theme(axis.text = element_text(color = "black"), plot.title = element_text(face = "bold"))
}
save_both <- function(plot, file, width, height) {
  pdf(file, width = width, height = height, useDingbats = FALSE); print(plot); dev.off()
  png(sub("\\.pdf$", ".png", file), width = width * 160, height = height * 160, res = 160)
  print(plot); dev.off()
}

ped <- fread(file.path(m3, "processed/target_tall_subtype.tsv"))
ped <- ped[, .(cohort = "TARGET pediatric RNA-seq", ZEB1 = as.numeric(ZEB1))]
# GSE142522 adult from module2
cand <- c(file.path(m2, "processed/adult_tall.tsv"),
          file.path(m2, "processed/sample_scores.tsv"),
          file.path(m2, "processed/gse142522_samples.tsv"))
ad <- NULL
for (f in cand) if (file.exists(f) && is.null(ad)) {
  tmp <- fread(f)
  if ("ZEB1_log" %in% names(tmp)) tmp[, ZEB1 := ZEB1_log]
  if ("ZEB1" %in% names(tmp)) {
    if ("group" %in% names(tmp)) tmp <- tmp[grepl("T-ALL|TALL", as.character(group), ignore.case = TRUE)]
    if ("class" %in% names(tmp)) tmp <- tmp[grepl("T-ALL|TALL", as.character(class), ignore.case = TRUE)]
    if ("is_tall" %in% names(tmp)) {
      v <- as.character(tmp$is_tall)
      tmp <- tmp[v %in% c("TRUE", "True", "true", "1")]
    }
    if (nrow(tmp)) ad <- tmp[, .(cohort = "GSE142522 adult SOLiD", ZEB1 = as.numeric(ZEB1))]
  }
}
if (is.null(ad) || !nrow(ad)) {
  p <- ggplot() + annotate("text", 0, 0, label = "1.17  Same-matrix child vs adult is not available.\nTARGET = pediatric STAR TPM; GSE142522 = adult SOLiD log2 counts.\nDistributions plotted only if adult table is found.") + theme_void()
  save_both(p, file.path(fig, "M1_1.17_child_vs_adult_note.pdf"), 8.0, 3.4)
} else {
  d <- rbind(ped, ad, fill = TRUE)
  d[, z := as.numeric(scale(ZEB1)), by = cohort]
  p <- ggplot(d, aes(cohort, z, fill = cohort)) +
    geom_boxplot(width = 0.55, outlier.shape = NA) +
    geom_point(position = position_jitter(width = 0.08, seed = 1), size = 1.2, alpha = 0.55) +
    scale_fill_manual(values = c("TARGET pediatric RNA-seq" = "#E64B35",
                                 "GSE142522 adult SOLiD" = "#4DBBD5"), guide = "none") +
    labs(title = "1.17  Child vs adult ZEB1 (within-cohort z)",
         subtitle = "Not a same-matrix contrast; platforms differ. Descriptive only.",
         x = NULL, y = "Within-cohort z(ZEB1)") + theme_sci() +
    theme(axis.text.x = element_text(angle = 20, hjust = 1))
  save_both(p, file.path(fig, "M1_1.17_child_vs_adult_z.pdf"), 7.2, 5.6)
}
cat("M1 1.17 DONE\n")

# Module 10 figures: PDX L-IC vs DP (GSE260697)
suppressPackageStartupMessages({
  library(data.table); library(ggplot2)
})
root <- "/data-b/liangfuhua/projects/ZEB1_TALL_analysis/module10"
fig <- file.path(root, "figures"); tab <- file.path(root, "tables")
dir.create(fig, FALSE, TRUE)
theme_sci <- function(base_size = 11) {
  theme_classic(base_size = base_size) +
    theme(axis.text = element_text(color = "black"),
          plot.title = element_text(face = "bold"))
}
save_both <- function(plot, file, width, height) {
  tryCatch({
    pdf(file, width = width, height = height, useDingbats = FALSE); print(plot); dev.off()
    png(sub("\\.pdf$", ".png", file), width = width * 160, height = height * 160, res = 160)
    print(plot); dev.off()
  }, error = function(e) message("FAIL ", e$message))
}
long <- fread(file.path(tab, "M10_10.1_keygene_values.tsv"))
tst <- fread(file.path(tab, "M10_10.2_LIC_vs_DP.tsv"))
pal <- c(LIC = "#E64B35", DP = "#4DBBD5", other = "#8491B4")
d <- long[kind == "fpkm" & group %in% c("LIC", "DP") & gene %in% c("ZEB1", "LMO2", "CD1A", "TAL1", "MKI67", "sig_ETP", "sig_cycle")]
if (nrow(d)) {
  p <- ggplot(d, aes(group, value, fill = group)) +
    geom_boxplot(width = 0.55, outlier.shape = NA, alpha = 0.85) +
    geom_point(position = position_jitter(width = 0.12, seed = 1), size = 1.8, alpha = 0.75) +
    facet_grid(patient ~ gene, scales = "free_y") +
    scale_fill_manual(values = pal, guide = "none") +
    labs(title = "GSE260697 PDX: L-IC (CD7+CD1a-) vs DP",
         subtitle = "REL-13 / REL-14 FPKM; n is sorted replicates within patient",
         x = NULL, y = "FPKM or signature z") + theme_sci()
  save_both(p, file.path(fig, "M10_10.1_LIC_vs_DP_keys.pdf"), 11.5, 6.5)
}
if (nrow(tst)) {
  z <- tst[gene == "ZEB1"]
  fwrite(tst, file.path(tab, "M10_10.2_LIC_vs_DP.tsv"), sep = "\t")
  p <- ggplot(tst[gene %in% c("ZEB1", "LMO2", "CD1A", "MKI67", "sig_ETP", "sig_cycle")],
              aes(reorder(gene, delta_LIC_minus_DP), delta_LIC_minus_DP, fill = patient)) +
    geom_hline(yintercept = 0, linetype = 2) +
    geom_col(position = position_dodge(width = 0.7), width = 0.65) +
    coord_flip() +
    labs(title = "L-IC minus DP (FPKM mean)", x = NULL, y = "delta") + theme_sci()
  save_both(p, file.path(fig, "M10_10.2_delta_bar.pdf"), 7.2, 5.2)
}
message("M10 figures done")

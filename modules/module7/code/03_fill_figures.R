# Fill missing Module 7 numbered figures from existing tables.
suppressPackageStartupMessages({ library(data.table); library(ggplot2) })
root <- "/data-b/liangfuhua/projects/ZEB1_TALL_analysis/module7"
fig <- file.path(root, "figures"); tab <- file.path(root, "tables")
theme_sci <- function(base_size = 11) {
  theme_classic(base_size = base_size) +
    theme(axis.text = element_text(color = "black"), axis.title = element_text(color = "black"),
          axis.line = element_line(linewidth = 0.45, color = "black"),
          plot.title = element_text(face = "bold", size = base_size + 1),
          plot.subtitle = element_text(size = base_size - 1, color = "grey25"))
}
save_both <- function(plot, file, width, height) {
  pdf(file, width = width, height = height, useDingbats = FALSE); print(plot); dev.off()
  png(sub("\\.pdf$", ".png", file), width = width * 160, height = height * 160, res = 160)
  print(plot); dev.off()
}

reg <- fread(file.path(tab, "M7_7.5_7.13_regressions.tsv"))
sp <- fread(file.path(tab, "M7_7.7_7.8_drug_spearman.tsv"))
q14 <- fread(file.path(tab, "M7_7.6_Q1Q4.tsv"))
short <- function(x) gsub(" \\(.*\\)$", "", x)

# 7.3-7.4 MTX absent
p <- ggplot() +
  annotate("text", x = 0, y = 0, label = "7.3 / 7.4 / 7.19  MTX is not in the 18-drug Pharmacotype panel.\nNo LC50, Spearman, or MTX-MRD mediation was computed.\nDo not rewrite the paper around MTX.") +
  theme_void() + xlim(-1, 1) + ylim(-1, 1)
save_both(p, file.path(fig, "M7_7.3_7.4_7.19_no_MTX.pdf"), 8.2, 3.6)

# 7.5 robust regression
d <- reg[model == "rlm_ZEB1"]
d[, drug := factor(short(drug), levels = short(drug)[order(estimate)])]
p <- ggplot(d, aes(estimate, drug, fill = p < 0.05)) +
  geom_vline(xintercept = 0, linetype = 2) +
  geom_col(width = 0.7, color = "white") +
  scale_fill_manual(values = c("TRUE" = "#E64B35", "FALSE" = "#8491B4")) +
  labs(title = "7.5  Robust regression ZEB1 vs LC50",
       subtitle = "All ZEB1-drug FDR > 0.05 in Spearman screen",
       x = "rlm coefficient", y = NULL, fill = "P<0.05") + theme_sci()
save_both(p, file.path(fig, "M7_7.5_robust_regression.pdf"), 8.0, 6.6)

# 7.6 Q1 vs Q4 dumbbell
q14[, drug := factor(short(drug), levels = short(drug)[order(median_q4 - median_q1)])]
long <- melt(q14, id.vars = c("drug", "wilcox_p", "fdr"),
             measure.vars = c("median_q1", "median_q4"), variable.name = "q", value.name = "lc50")
p <- ggplot(long, aes(lc50, drug, color = q)) +
  geom_line(aes(group = drug), color = "grey70") +
  geom_point(size = 2.8) +
  scale_x_log10() +
  scale_color_manual(values = c(median_q1 = "#3C5488", median_q4 = "#E64B35"),
                     labels = c("Q1 ZEB1-low", "Q4 ZEB1-high")) +
  labs(title = "7.6  Median LC50 in ZEB1 Q1 vs Q4", x = "Median LC50", y = NULL, color = NULL) + theme_sci()
save_both(p, file.path(fig, "M7_7.6_Q1Q4_dumbbell.pdf"), 8.2, 6.6)

# 7.8 FDR bar
z <- sp[gene == "ZEB1_log"]
setorder(z, rho)
z[, drug_short := factor(drug_short, levels = drug_short)]
p <- ggplot(z, aes(rho, drug_short, fill = fdr < 0.05)) +
  geom_col(width = 0.7, color = "white") +
  geom_text(aes(label = sprintf("FDR=%.2f", fdr)), hjust = ifelse(z$rho >= 0, -0.05, 1.05), size = 2.6) +
  scale_fill_manual(values = c("TRUE" = "#E64B35", "FALSE" = "#8491B4")) +
  labs(title = "7.8  ZEB1 vs LC50 Spearman (BH-FDR)", x = "rho", y = NULL, fill = "FDR<0.05") + theme_sci()
save_both(p, file.path(fig, "M7_7.8_FDR_bars.pdf"), 8.0, 6.6)

# 7.11-7.13 adjusted / interaction
for (mod in c("ZEB1+ETP", "ZEB1+age", "ZEB1xETP")) {
  dd <- reg[model == mod]
  dd[, drug := factor(short(drug), levels = short(drug)[order(estimate)])]
  p <- ggplot(dd, aes(estimate, drug, fill = p < 0.05)) +
    geom_vline(xintercept = 0, linetype = 2) +
    geom_col(width = 0.7, color = "white") +
    scale_fill_manual(values = c("TRUE" = "#E64B35", "FALSE" = "#8491B4"), guide = "none") +
    labs(title = paste("7.11-7.13 ", mod), x = "coefficient", y = NULL) + theme_sci()
  fn <- switch(mod,
               "ZEB1+ETP" = "M7_7.11_ETP_adjusted.pdf",
               "ZEB1+age" = "M7_7.12_age_adjusted.pdf",
               "M7_7.13_ZEB1xETP.pdf")
  save_both(p, file.path(fig, fn), 8.0, 6.4)
}

# 7.14-7.16 gene comparison heatmap of rho
mm <- dcast(sp, drug_short ~ gene, value.var = "rho")
long <- melt(mm, id.vars = "drug_short", variable.name = "gene", value.name = "rho")
p <- ggplot(long, aes(gene, drug_short, fill = rho)) +
  geom_tile(color = "white") +
  scale_fill_gradient2(low = "#3C5488", mid = "white", high = "#E64B35", midpoint = 0) +
  labs(title = "7.14-7.17  LC50 Spearman: ZEB1 / ZEB2 / LMO2 / ratio",
       x = NULL, y = NULL, fill = "rho") + theme_sci() + theme(axis.line = element_blank())
save_both(p, file.path(fig, "M7_7.14_7.17_gene_drug_rho.pdf"), 7.4, 6.8)

# 7.20 pointer
p <- ggplot() +
  annotate("text", x = 0, y = 0, label = "7.20  Independent cell-line drug validation is Module 8\n(DepMap CRISPR + Pharmacotype target mapping; no PRISM file).") +
  theme_void() + xlim(-1, 1) + ylim(-1, 1)
save_both(p, file.path(fig, "M7_7.20_see_module8.pdf"), 8.0, 3.2)
cat("M7 FILL DONE\n")

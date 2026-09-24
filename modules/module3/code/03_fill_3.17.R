# 3.17-3.18 ZEB1 mutation/fusion are 0/265 — still draw the count.
suppressPackageStartupMessages({ library(data.table); library(ggplot2) })
root <- "/data-b/liangfuhua/projects/ZEB1_TALL_analysis/module3"
fig <- file.path(root, "figures"); tab <- file.path(root, "tables")
d <- fread(file.path(tab, "M3_3.17_3.18_ZEB1_mutation_fusion.tsv"))
long <- melt(d, id.vars = "gene", measure.vars = c("n_mut", "n_fusion"),
             variable.name = "event", value.name = "count")
long[, event := fifelse(event == "n_mut", "coding mutation", "fusion")]
long[is.na(count), count := 0]
p <- ggplot(long, aes(gene, count, fill = event)) +
  geom_col(position = position_dodge(0.7), width = 0.65, color = "white") +
  geom_text(aes(label = paste0(count, "/", unique(d$n))), position = position_dodge(0.7), vjust = -0.3, size = 3) +
  scale_fill_manual(values = c("coding mutation" = "#E64B35", fusion = "#3C5488")) +
  labs(title = "3.17-3.18  ZEB1 locus alteration in TARGET T-ALL",
       subtitle = "0 coding mutations and 0 fusions in 265 patients",
       x = NULL, y = "Events", fill = NULL) +
  theme_classic(base_size = 11)
pdf(file.path(fig, "M3_3.17_3.18_mutation_fusion.pdf"), 6.8, 5.0, useDingbats = FALSE)
print(p); dev.off()
png(file.path(fig, "M3_3.17_3.18_mutation_fusion.png"), width = 6.8 * 160, height = 5.0 * 160, res = 160)
print(p); dev.off()
p2 <- ggplot() +
  annotate("text", 0, 0,
           label = "3.20  CIMP x MRD\nGSE272023 GEO metadata has no MRD table.\nNot filled. CIMP vs ZEB1 is 3.19 / Module 6.18-6.20.",
           size = 4.2, lineheight = 1.15) +
  theme_void() + xlim(-1, 1) + ylim(-1, 1)
pdf(file.path(fig, "M3_3.20_CIMP_MRD_note.pdf"), 8.0, 3.4, useDingbats = FALSE)
print(p2); dev.off()
png(file.path(fig, "M3_3.20_CIMP_MRD_note.png"), width = 8.0 * 160, height = 3.4 * 160, res = 160)
print(p2); dev.off()
cat("M3 3.17 DONE\n")

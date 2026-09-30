suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})
out_dir <- "D:/_bioinformation/ZEB1/data/analysis/total/data/validation/polonen_round1"
fig_dir <- file.path(out_dir, "figures")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
d <- fread(file.path(out_dir, "patient_level_frozen_balance.tsv"))
d$subtype <- gsub("γ", "g", d$subtype, fixed = TRUE)
d$subtype <- gsub("α", "a", d$subtype, fixed = TRUE)
d$subtype <- gsub("δ", "d", d$subtype, fixed = TRUE)
ord <- d[, .(med = median(balance)), by = subtype][order(med)]
d$subtype <- factor(d$subtype, levels = ord$subtype)

long <- melt(d, id.vars = "subtype", measure.vars = c("ZEB1", "ZEB2", "LMO2", "balance"),
             variable.name = "metric", value.name = "value")
long$metric <- factor(long$metric, levels = c("ZEB1", "ZEB2", "LMO2", "balance"))

p <- ggplot(long, aes(subtype, value)) +
  geom_hline(data = data.frame(metric = factor("balance", levels = levels(long$metric)), y = 0),
             aes(yintercept = y), linewidth = 0.3, colour = "grey40") +
  geom_jitter(width = 0.18, height = 0, size = 0.35, alpha = 0.35, colour = "#3d4f66") +
  stat_summary(fun = median, geom = "crossbar", width = 0.55, linewidth = 0.35, colour = "#9a3412") +
  facet_wrap(~metric, scales = "free_x", ncol = 2) +
  coord_flip() +
  labs(
    title = "Pölönen 1,309 pediatric T-ALL: ZEB-family state by author subtype",
    subtitle = "Points are patients. Bars are medians. Balance = z(TMM logCPM ZEB1) − z(TMM logCPM ZEB2), z-scored inside this cohort.",
    x = NULL, y = "log2 CPM  or  balance (z ZEB1 − z ZEB2)"
  ) +
  theme_bw(base_size = 11) +
  theme(strip.background = element_rect(fill = "white"),
        panel.grid.major.y = element_blank())
ggsave(file.path(fig_dir, "subtype_zeb_landscape.pdf"), p, width = 11, height = 9)
ggsave(file.path(fig_dir, "subtype_zeb_landscape.png"), p, width = 11, height = 9, dpi = 160)
message("wrote landscape")

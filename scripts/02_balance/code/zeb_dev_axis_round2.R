# ZEB1-ZEB2 developmental axis.
# Balance = z(ZEB1) - z(ZEB2), z computed inside each normal dataset.
# Pölönen: author immunophenotype state versus 17-subtype variance.

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

out_dir <- "D:/_bioinformation/ZEB1/data/analysis/total/data/validation/zeb_developmental_axis"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
fig_dir <- file.path(out_dir, "figures")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

zsd <- function(x) {
  s <- stats::sd(x)
  if (!is.finite(s) || s == 0) return(rep(NA_real_, length(x)))
  (x - mean(x)) / s
}

add_balance <- function(df, zeb1, zeb2) {
  df$z_ZEB1 <- zsd(df[[zeb1]])
  df$z_ZEB2 <- zsd(df[[zeb2]])
  df$balance <- df$z_ZEB1 - df$z_ZEB2
  df$raw_diff <- df[[zeb1]] - df[[zeb2]]
  df
}

# ---- GSE142522 sorted bulk, one library per stage ----
b <- fread("D:/_bioinformation/ZEB1/data/analysis/total/data/SourceData_Fig2A_thymus_trajectory.tsv")
stages_142 <- c("cd34", "isp", "ec", "lc", "sp4", "sp8")
labels_142 <- c("CD34+", "ISP", "early cortical", "late cortical", "CD4 SP", "CD8 SP")
row <- b[gene %in% c("ZEB1", "ZEB2", "LMO2", "LYL1")]
wide <- as.data.table(t(as.matrix(row[, ..stages_142])))
setnames(wide, row$gene)
wide$stage <- labels_142
wide$order <- seq_len(nrow(wide))
wide$dataset <- "GSE142522"
wide <- add_balance(wide, "ZEB1", "ZEB2")
wide$unit <- "one sorted bulk library"
g142 <- wide

# ---- GSE195812 pooled FACS libraries ----
a <- fread("D:/_bioinformation/ZEB1/data/analysis/modules/module2/tables/M2_r2_GSE195812_sample_means.tsv")
order_195 <- c("DN1", "DN2", "DN3", "ISP", "DP_CD3neg", "DP_CD3pos", "CD4SP", "CD8SP")
a <- a[stage %in% order_195]
a$order <- match(a$stage, order_195)
setorder(a, order)
a$dataset <- "GSE195812"
# z across the eight stage libraries, not across genes
tmp <- data.table(ZEB1 = a$g_ZEB1, ZEB2 = a$g_ZEB2, LMO2 = a$g_LMO2, LYL1 = a$g_LYL1,
                  MEF2C = a$g_MEF2C, stage = a$stage, order = a$order, dataset = a$dataset)
tmp <- add_balance(tmp, "ZEB1", "ZEB2")
tmp$unit <- "one pooled FACS library"
g195 <- tmp

# ---- GSE206710 donor x fraction; axis starts at DN, not at a CD34 progenitor ----
c <- fread("D:/_bioinformation/ZEB1/data/analysis/modules/module2/tables/M2_r2_GSE206710_sample_means.tsv")
order_206 <- c("DN_early", "DP_immature", "DP_mature")
c <- c[stage %in% order_206]
# within-donor z, then average
bits <- list()
for (don in unique(c$donor)) {
  sub <- c[donor == don]
  if (nrow(sub) < 2) next
  sub$z_ZEB1 <- zsd(sub$g_ZEB1)
  sub$z_ZEB2 <- zsd(sub$g_ZEB2)
  sub$balance <- sub$z_ZEB1 - sub$z_ZEB2
  sub$raw_diff <- sub$g_ZEB1 - sub$g_ZEB2
  bits[[don]] <- sub
}
c2 <- rbindlist(bits)
c2$order <- match(c2$stage, order_206)
g206_donor <- c2[, .(
  ZEB1 = mean(g_ZEB1), ZEB2 = mean(g_ZEB2), LMO2 = mean(g_LMO2), LYL1 = mean(g_LYL1),
  MEF2C = mean(g_MEF2C), balance = mean(balance), raw_diff = mean(raw_diff),
  n_donors = uniqueN(donor)
), by = .(stage, order)]
g206_donor$dataset <- "GSE206710"
g206_donor$unit <- "within-donor z, then mean across donors"

# ---- Park/HTA donor x annotated cell type ----
p <- fread("D:/_bioinformation/ZEB1/data/analysis/total/data/SourceData_ParkHTA_donor_celltype.tsv")
keep <- c(
  "double negative thymocyte",
  "double-positive, alpha-beta thymocyte",
  "CD4-positive, alpha-beta T cell",
  "CD8-positive, alpha-beta T cell"
)
lab <- c(
  "double negative thymocyte" = "DN",
  "double-positive, alpha-beta thymocyte" = "DP",
  "CD4-positive, alpha-beta T cell" = "CD4 SP",
  "CD8-positive, alpha-beta T cell" = "CD8 SP"
)
p <- p[celltype %in% keep & n_cells >= 20]
p$stage <- unname(lab[p$celltype])
pbits <- list()
for (don in unique(p$donor)) {
  sub <- p[donor == don]
  if (uniqueN(sub$stage) < 3) next
  sub$z_ZEB1 <- zsd(sub$ZEB1)
  sub$z_ZEB2 <- zsd(sub$ZEB2)
  sub$balance <- sub$z_ZEB1 - sub$z_ZEB2
  sub$raw_diff <- sub$ZEB1 - sub$ZEB2
  pbits[[don]] <- sub
}
pd <- rbindlist(pbits)
pd$order <- match(pd$stage, c("DN", "DP", "CD4 SP", "CD8 SP"))
gpark <- pd[, .(
  ZEB1 = mean(ZEB1), ZEB2 = mean(ZEB2), LMO2 = mean(LMO2),
  balance = mean(balance), raw_diff = mean(raw_diff),
  n_donors = uniqueN(donor), n_cells = sum(n_cells)
), by = .(stage, order)]
gpark$dataset <- "Park/HTA"
gpark$unit <- "within-donor z among donors with at least 3 of DN/DP/SP"

traj <- rbindlist(list(
  g142[, .(dataset, stage, order, ZEB1, ZEB2, LMO2, LYL1, balance, raw_diff, unit)],
  g195[, .(dataset, stage, order, ZEB1, ZEB2, LMO2, LYL1, MEF2C, balance, raw_diff, unit)],
  g206_donor[, .(dataset, stage, order, ZEB1, ZEB2, LMO2, LYL1, MEF2C, balance, raw_diff, unit, n_donors)],
  gpark[, .(dataset, stage, order, ZEB1, ZEB2, LMO2, balance, raw_diff, unit, n_donors, n_cells)]
), fill = TRUE)
fwrite(traj, file.path(out_dir, "normal_thymus_balance_trajectory.tsv"), sep = "\t")

# pattern calls: early vs cortical/DP balance, and whether ZEB2 falls into that window
call_one <- function(df, early, mid, late) {
  e <- mean(df$balance[df$stage %in% early])
  m <- mean(df$balance[df$stage %in% mid])
  l <- mean(df$balance[df$stage %in% late])
  z2e <- mean(df$ZEB2[df$stage %in% early])
  z2m <- mean(df$ZEB2[df$stage %in% mid])
  z1e <- mean(df$ZEB1[df$stage %in% early])
  z1m <- mean(df$ZEB1[df$stage %in% mid])
  z1l <- mean(df$ZEB1[df$stage %in% late])
  z2l <- mean(df$ZEB2[df$stage %in% late])
  data.table(
    dataset = df$dataset[1],
    balance_early = e, balance_mid = m, balance_late = l,
    mid_minus_early = m - e,
    ZEB2_falls_into_mid = z2m < z2e,
    ZEB1_higher_or_held_at_mid = z1m >= z1e * 0.9,
    both_lower_at_late_than_peak = (z1l < max(df$ZEB1) * 0.9) && (z2l <= z2m * 1.25 || z2l < z2e)
  )
}
calls <- rbindlist(list(
  call_one(g142, "CD34+", c("early cortical", "late cortical"), c("CD4 SP", "CD8 SP")),
  call_one(g195, c("DN1", "DN2"), c("DN3", "ISP", "DP_CD3neg"), c("CD4SP", "CD8SP")),
  call_one(g206_donor, "DN_early", "DP_immature", "DP_mature"),
  call_one(gpark, "DN", "DP", c("CD4 SP", "CD8 SP"))
))
fwrite(calls, file.path(out_dir, "trajectory_pattern_calls.tsv"), sep = "\t")

# plot
plot_df <- traj[, .(dataset, stage, order, balance, ZEB1, ZEB2)]
plot_df$stage <- factor(plot_df$stage, levels = unique(plot_df$stage[order(plot_df$dataset, plot_df$order)]))
longp <- melt(plot_df, id.vars = c("dataset", "stage", "order"),
              measure.vars = c("ZEB1", "ZEB2", "balance"),
              variable.name = "metric", value.name = "value")
# free y, but balance and expression are different scales: facet metric x dataset is messy.
# two-row: expression trajectories are dataset-native; balance is comparable in shape only.
p1 <- ggplot(longp[metric != "balance"], aes(order, value, colour = metric, group = metric)) +
  geom_line(linewidth = 0.6) +
  geom_point(size = 1.8) +
  facet_wrap(~dataset, scales = "free_y", ncol = 2) +
  scale_x_continuous(breaks = function(x) seq(min(x), max(x), by = 1)) +
  scale_colour_manual(values = c(ZEB1 = "#9a3412", ZEB2 = "#1d4e89")) +
  labs(title = "Normal thymus: ZEB1 and ZEB2 along developmental order",
       subtitle = "Each panel uses that dataset's own expression scale. X is stage order, not a shared time unit.",
       x = "Developmental order within the dataset", y = "Expression (dataset-native)") +
  theme_bw(base_size = 11) +
  theme(strip.background = element_rect(fill = "white"))
p2 <- ggplot(longp[metric == "balance"], aes(order, value)) +
  geom_hline(yintercept = 0, linewidth = 0.3, colour = "grey40") +
  geom_line(linewidth = 0.6, colour = "#3d4f66") +
  geom_point(size = 1.8, colour = "#3d4f66") +
  facet_wrap(~dataset, scales = "free_x", ncol = 2) +
  labs(title = "Within-dataset balance: z(ZEB1) - z(ZEB2)",
       subtitle = "z is inside the dataset. GSE206710 and Park/HTA are within-donor z, then averaged.",
       x = "Developmental order within the dataset", y = "Balance") +
  theme_bw(base_size = 11) +
  theme(strip.background = element_rect(fill = "white"))
ggsave(file.path(fig_dir, "normal_zeb_trajectories.pdf"), p1, width = 9, height = 7)
ggsave(file.path(fig_dir, "normal_zeb_trajectories.png"), p1, width = 9, height = 7, dpi = 150)
ggsave(file.path(fig_dir, "normal_zeb_balance.pdf"), p2, width = 9, height = 7)
ggsave(file.path(fig_dir, "normal_zeb_balance.png"), p2, width = 9, height = 7, dpi = 150)

# stage labels for the plot axes
labs <- traj[, .(dataset, order, stage)]
fwrite(labs, file.path(out_dir, "stage_order_labels.tsv"), sep = "\t")

# ---- Pölönen author immunophenotype ----
e <- new.env(parent = emptyenv())
load("D:/_bioinformation/ZEB1/data/analysis/data/polonen_syn54032669/Data_1309Samples.RData", envir = e)
ip <- data.table(sample_id = names(e$M7.IP.factor), ip = as.character(e$M7.IP.factor))
# factor may be a data.frame
if (is.data.frame(e$M7.IP.factor)) {
  ip <- data.table(sample_id = rownames(e$M7.IP.factor), ip = as.character(e$M7.IP.factor[[1]]))
}
pt <- fread("D:/_bioinformation/ZEB1/data/analysis/total/data/validation/polonen_round1/patient_level_frozen_balance.tsv")
pt <- merge(pt, ip, by = "sample_id", all.x = TRUE)
pt$subtype <- factor(pt$subtype)
pt$ip <- factor(pt$ip)
r_sub <- summary(lm(balance ~ subtype, pt))$r.squared
r_ip <- summary(lm(balance ~ ip, pt))$r.squared
r_both <- summary(lm(balance ~ subtype + ip, pt))$r.squared
a_ip <- anova(lm(balance ~ subtype, pt), lm(balance ~ subtype + ip, pt))
a_sub <- anova(lm(balance ~ ip, pt), lm(balance ~ ip + subtype, pt))
var_tab <- data.table(
  model = c("subtype", "immunophenotype", "subtype + immunophenotype"),
  r_squared = c(r_sub, r_ip, r_both),
  p_added_term = c(NA, a_sub$`Pr(>F)`[2], a_ip$`Pr(>F)`[2])
)
# wait I swapped. Fix explicitly.
m_sub <- lm(balance ~ subtype, pt)
m_ip <- lm(balance ~ ip, pt)
m_both_s <- lm(balance ~ subtype + ip, pt)
m_both_i <- lm(balance ~ ip + subtype, pt)
var_tab <- data.table(
  model = c("17 subtype", "author immunophenotype (4)", "subtype then immunophenotype", "immunophenotype then subtype"),
  r_squared = c(summary(m_sub)$r.squared, summary(m_ip)$r.squared, summary(m_both_s)$r.squared, summary(m_both_i)$r.squared),
  extra_p = c(NA, NA, anova(m_sub, m_both_s)$`Pr(>F)`[2], anova(m_ip, m_both_i)$`Pr(>F)`[2])
)
ip_sum <- pt[, .(
  n = .N,
  median_balance = median(balance),
  mean_balance = mean(balance),
  mean_ZEB1 = mean(ZEB1),
  mean_ZEB2 = mean(ZEB2),
  induction_failure = sum(if_author)
), by = ip]
fwrite(var_tab, file.path(out_dir, "polonen_balance_variance.tsv"), sep = "\t")
fwrite(ip_sum, file.path(out_dir, "polonen_balance_by_immunophenotype.tsv"), sep = "\t")

# ---- BCL11B subtype: which enhancer-hijack partners are on in those 18 ----
fm_env <- new.env(parent = emptyenv())
load("D:/_bioinformation/ZEB1/data/analysis/data/polonen_syn54032669/TALL_X01_FeatureMatrix_genomics_subtype.Rdata", envir = fm_env)
fm <- fm_env$fm
bcl <- pt[subtype == "BCL11B", sample_id]
rows <- grep("BCL11B", rownames(fm), value = TRUE)
rows <- rows[grepl("EnhHj|Translocation|Fusion|ZEB2", rows)]
subm <- fm[rows, bcl, drop = FALSE]
# keep partners present in at least one BCL11B-subtype sample
present <- rows[rowSums(subm != 0, na.rm = TRUE) > 0]
partner <- data.table(feature = present, n_of_18 = as.integer(rowSums(fm[present, bcl, drop = FALSE] != 0)))
setorder(partner, -n_of_18)
fwrite(partner, file.path(out_dir, "bcl11b_subtype_lesion_partners.tsv"), sep = "\t")

message("trajectory calls:")
print(calls)
message("variance:")
print(var_tab)
message("immunophenotype:")
print(ip_sum)
message("BCL11B partners with any hit:")
print(partner)

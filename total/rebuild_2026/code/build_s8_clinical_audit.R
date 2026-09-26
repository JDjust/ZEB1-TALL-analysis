# Supplementary Figure S8: clinical attenuation and Lim response sensitivity.
suppressPackageStartupMessages({
  library(ggplot2); library(patchwork); library(readr); library(dplyr)
})
argv <- grep('^--file=', commandArgs(FALSE), value = TRUE)
script <- normalizePath(sub('^--file=', '', argv), winslash = '/')
root <- normalizePath(file.path(dirname(script), '../../..'), winslash = '/')
out <- file.path(root, 'total/rebuild_2026/edition_20260925/figures/supplementary')
dir.create(out, recursive = TRUE, showWarnings = FALSE)
rd <- function(path) as.data.frame(read_tsv(path, show_col_types = FALSE,
                                            progress = FALSE))
table_dir <- file.path(root, 'total/rebuild_2026/data/source_data_rebuilt')
gate <- file.path(root, 'total/data/validation/gateA_lim2025')
base <- rd(file.path(table_dir, 'supplementary_tables/Table_S1_clinical_models.tsv'))
full <- rd(file.path(table_dir, 'revision2_statistics/clinical_age_sex_wbc_sensitivity.tsv'))
lim <- rd(file.path(gate, 'q1_baseline.tsv'))
lim_adj <- rd(file.path(gate, 'q1_sensitivity_ols.tsv'))
day28 <- rd(file.path(gate, 'q2_longitudinal.tsv'))
stopifnot(nrow(base) == 12L, nrow(full) == 6L)
full$adjustment <- 'Subtype + age + sex + WBC'
full$measure <- ifelse(full$method == 'Cox PH', 'HR', 'OR')
clinical <- rbind(base[, c('endpoint', 'measure', 'adjustment', 'estimate',
                           'ci_low', 'ci_high', 'p', 'endpoint_bh_fdr')],
                  full[, c('endpoint', 'measure', 'adjustment', 'estimate',
                           'ci_low', 'ci_high', 'p', 'endpoint_bh_fdr')])
clinical$adjustment <- factor(clinical$adjustment,
  levels = c('Unadjusted', 'Subtype-adjusted', 'Subtype + age + sex + WBC'))
ink <- '#203744'; grey <- '#9eb0b8'; purple <- '#776298'
teal <- '#168b82'; orange <- '#c78a32'
model_cols <- c('Unadjusted' = orange, 'Subtype-adjusted' = teal,
                'Subtype + age + sex + WBC' = purple)
theme_plate <- theme_classic(base_size = 8, base_family = 'Arial') +
  theme(plot.title = element_text(size = 9, face = 'bold', colour = ink),
        plot.subtitle = element_text(size = 7.1, colour = ink),
        axis.title = element_text(size = 7.8, colour = ink),
        axis.text = element_text(size = 7.1, colour = ink),
        legend.text = element_text(size = 7.0, colour = ink),
        panel.grid.major.y = element_line(colour = '#e5edf0', linewidth = .28),
        plot.margin = margin(4, 5, 4, 5))

one_endpoint <- function(endpoint, letter, short) {
  d <- clinical[clinical$endpoint == endpoint, ]
  ggplot(d, aes(estimate, adjustment, colour = adjustment)) +
    geom_vline(xintercept = 1, colour = grey, linewidth = .4) +
    geom_segment(aes(x = ci_low, xend = ci_high, yend = adjustment),
                 linewidth = .7) + geom_point(size = 2.0) +
    scale_colour_manual(values = model_cols, guide = 'none') +
    scale_x_log10(limits = c(.4, 1.6), breaks = c(.5, 1, 1.5)) +
    labs(title = paste0(letter, '  ', short),
         subtitle = paste0(d$measure[1], ' per unit balance; 95% CI'),
         x = d$measure[1], y = NULL) + theme_plate
}
pA <- one_endpoint('Induction failure', 'A', 'Induction failure')
pB <- one_endpoint('MRD >=0.1%', 'B', 'MRD 0.1%')
pC <- one_endpoint('MRD >=0.01%', 'C', 'MRD 0.01%')

surv <- clinical[clinical$endpoint %in%
  c('Event-free survival', 'Overall survival'), ]
surv$endpoint <- factor(surv$endpoint,
  levels = c('Event-free survival', 'Overall survival'),
  labels = c('EFS', 'OS'))
pD <- ggplot(surv, aes(estimate, adjustment, colour = adjustment)) +
  geom_vline(xintercept = 1, colour = grey, linewidth = .4) +
  geom_segment(aes(x = ci_low, xend = ci_high, yend = adjustment),
               linewidth = .65) + geom_point(size = 1.7) +
  facet_wrap(~endpoint, ncol = 2) +
  scale_colour_manual(values = model_cols, guide = 'none') +
  scale_x_log10(limits = c(.55, 1.45), breaks = c(.6, 1, 1.4)) +
  labs(title = 'D  EFS and OS attenuation',
       subtitle = 'Hazard ratio per unit balance', x = 'HR', y = NULL) +
  theme_plate + theme(strip.text = element_text(size = 7))

full_plot <- full
full_plot$endpoint <- factor(full_plot$endpoint,
  levels = rev(c('Induction failure', 'M2/M3 morphology', 'MRD >=0.1%',
                 'MRD >=0.01%', 'Event-free survival', 'Overall survival')))
pE <- ggplot(full_plot, aes(estimate, endpoint)) +
  geom_vline(xintercept = 1, colour = grey, linewidth = .4) +
  geom_segment(aes(x = ci_low, xend = ci_high, yend = endpoint),
               colour = grey, linewidth = .7) +
  geom_point(aes(fill = endpoint_bh_fdr < .05), shape = 21,
             size = 2.0, colour = ink) +
  scale_fill_manual(values = c('TRUE' = purple, 'FALSE' = 'white'),
                    guide = 'none') +
  scale_x_log10(limits = c(.55, 1.25), breaks = c(.6, .8, 1, 1.2)) +
  labs(title = 'E  All six fully adjusted endpoints',
       subtitle = 'Filled: six-endpoint BH FDR < 0.05',
       x = 'OR or HR per unit balance', y = NULL) + theme_plate

lim_balance <- subset(lim, metric == 'balance')
lim_balance$cohort <- factor(lim_balance$cohort,
  levels = c('combined', 'discovery', 'extension'))
pF <- ggplot(lim_balance, aes(diff_a_minus_b, cohort)) +
  geom_vline(xintercept = 0, colour = grey, linewidth = .4) +
  geom_segment(aes(x = boot_ci95_lo, xend = boot_ci95_hi,
                   yend = cohort), colour = grey, linewidth = .7) +
  geom_point(size = 2, colour = teal) +
  labs(title = 'F  Lim baseline cohorts',
       subtitle = 'IF - responsive; bootstrap 95% CI',
       x = 'Mean balance difference', y = NULL) + theme_plate

lim_adj <- subset(lim_adj, y == 'balance')
lim_adj$model <- factor(lim_adj$model,
  levels = rev(c('unadjusted', 'plus_subtype_l1',
                 'plus_age_bin', 'plus_etp')),
  labels = rev(c('Unadjusted', '+ subtype', '+ age', '+ ETP')))
pG <- ggplot(lim_adj, aes(coef_refractory, model)) +
  geom_vline(xintercept = 0, colour = grey, linewidth = .4) +
  geom_segment(aes(x = ci95_lo, xend = ci95_hi, yend = model),
               colour = grey, linewidth = .7) +
  geom_point(size = 2, colour = purple) +
  labs(title = 'G  Lim adjustment models',
       subtitle = 'One-covariate sensitivity; 95% CI',
       x = 'IF coefficient', y = NULL) + theme_plate

day28 <- subset(day28, group == 'all_paired')
day28$metric <- factor(day28$metric,
  levels = rev(c('ZEB1_logcpm', 'ZEB2_logcpm', 'LMO2_logcpm', 'balance')),
  labels = rev(c('ZEB1', 'ZEB2', 'LMO2', 'Balance')))
pH <- ggplot(day28, aes(mean, metric)) +
  geom_vline(xintercept = 0, colour = grey, linewidth = .4) +
  geom_segment(aes(x = boot_ci95_lo, xend = boot_ci95_hi,
                   yend = metric), colour = grey, linewidth = .7) +
  geom_point(size = 2, colour = purple) +
  labs(title = 'H  Day28 null boundary',
       subtitle = '12 paired patients; bootstrap 95% CI',
       x = 'Mean Day28 - Day0', y = NULL) + theme_plate

fig <- (pA | pB) / (pC | pD) / (pE | pF) / (pG | pH) +
  plot_layout(heights = c(1, 1, 1.2, 1))
ggsave(file.path(out, 'S8.pdf'), fig, width = 7.09, height = 8.35,
       units = 'in', device = cairo_pdf, bg = 'white', limitsize = FALSE)
cat('S8 clinical audit: eight panels\n')

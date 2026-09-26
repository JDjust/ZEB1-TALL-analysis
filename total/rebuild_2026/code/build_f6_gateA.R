# Figure 6: patient-level Lim Gate A results at the final print dimensions.
# All effect estimates and intervals are read from the frozen Gate A outputs.
suppressPackageStartupMessages({
  library(ggplot2); library(patchwork); library(readr); library(dplyr); library(tidyr)
})
argv <- grep('^--file=', commandArgs(FALSE), value = TRUE)
script <- normalizePath(sub('^--file=', '', argv), winslash = '/')
root <- normalizePath(file.path(dirname(script), '../../..'), winslash = '/')
source_dir <- file.path(root, 'total/data/validation/gateA_lim2025')
out <- file.path(root, 'total/rebuild_2026/edition_20260925/figures/main')
dir.create(out, recursive = TRUE, showWarnings = FALSE)
rd <- function(name) as.data.frame(read_tsv(file.path(source_dir, name),
  show_col_types = FALSE, progress = FALSE))
day0 <- rd('eligible_patient_timepoint.tsv')
day0 <- subset(day0, timepoint == 'Day0')
pairs <- rd('q3_state_day0_pairs.tsv')
baseline <- rd('q1_baseline.tsv')
state <- rd('q3_state.tsv')
adjust <- rd('q1_sensitivity_ols.tsv')
stopifnot(nrow(day0) == 54L, nrow(pairs) == 41L,
          sum(day0$induction == 'IF') == 19L,
          sum(day0$induction == 'responsive') == 35L)

ink <- '#203744'; grid <- '#e5edf0'; teal <- '#168b82'
orange <- '#c78a32'; purple <- '#776298'; blue <- '#246b96'
red <- '#bc3446'; grey <- '#9eb0b8'
theme_plate <- theme_classic(base_size = 8, base_family = 'Arial') +
  theme(plot.title = element_text(size = 9, face = 'bold', colour = ink,
                                  margin = margin(b = 3)),
        plot.subtitle = element_text(size = 7.2, colour = ink),
        axis.title = element_text(size = 7.8, colour = ink),
        axis.text = element_text(size = 7.2, colour = ink),
        legend.text = element_text(size = 7, colour = ink),
        legend.title = element_text(size = 7, colour = ink),
        panel.grid.major.y = element_line(colour = grid, linewidth = .28),
        plot.margin = margin(4, 5, 4, 5))
cols <- c(IF = orange, responsive = teal)
metric_labels <- c(ZEB1_logcpm = 'ZEB1', ZEB2_logcpm = 'ZEB2',
                   LMO2_logcpm = 'LMO2', balance = 'Balance')

pairs$induction <- factor(pairs$induction_pos,
                          levels = c('IF', 'responsive'))
pairs$rank <- rank(pairs$d_balance, ties.method = 'first')
pairs$fraction_pos <- pairs$n_cells_pos /
  (pairs$n_cells_pos + pairs$n_cells_neg)

pA <- ggplot(pairs, aes(balance_neg, balance_pos, colour = induction)) +
  geom_abline(slope = 1, intercept = 0, colour = grey,
              linetype = 'dashed', linewidth = .35) +
  geom_point(size = 1.35, alpha = .82) +
  scale_colour_manual(values = cols, guide = 'none') +
  coord_equal() +
  labs(title = 'A  Paired ZBTB16 states',
       subtitle = '41 patients; paired Day0 states',
       x = 'ZBTB16 neg balance', y = 'ZBTB16 pos balance') + theme_plate

pB <- ggplot(pairs, aes(rank, d_balance)) +
  geom_hline(yintercept = 0, colour = grey, linewidth = .4) +
  geom_segment(aes(xend = rank, y = 0, yend = d_balance,
                   colour = d_balance < 0), linewidth = .55) +
  geom_point(aes(colour = d_balance < 0), size = 1.25) +
  scale_colour_manual(values = c('TRUE' = teal, 'FALSE' = grey), guide = 'none') +
  labs(title = 'B  Ranked state difference',
       subtitle = 'Positive minus negative state',
       x = 'Patient rank', y = 'Balance difference') + theme_plate

state_all <- subset(state, group == 'all_day0_pos_minus_neg')
state_all$label <- factor(metric_labels[state_all$metric],
  levels = rev(c('ZEB1', 'ZEB2', 'LMO2', 'Balance')))
pC <- ggplot(state_all, aes(mean, label)) +
  geom_vline(xintercept = 0, colour = grey, linewidth = .4) +
  geom_segment(aes(x = boot_ci95_lo, xend = boot_ci95_hi, yend = label),
               colour = grey, linewidth = .65) +
  geom_point(aes(colour = label), size = 2.0) +
  scale_colour_manual(values = c(ZEB1 = red, ZEB2 = blue,
                                 LMO2 = orange, Balance = purple), guide = 'none') +
  labs(title = 'C  State component changes',
       subtitle = 'Mean paired difference; bootstrap 95% CI',
       x = 'Positive minus negative state', y = NULL) + theme_plate

pD <- ggplot(pairs, aes(fraction_pos, d_balance, colour = induction)) +
  geom_hline(yintercept = 0, colour = grey, linewidth = .4) +
  geom_point(size = 1.4, alpha = .85) +
  scale_colour_manual(values = cols, name = NULL) +
  labs(title = 'D  State fraction and effect',
       subtitle = 'Each dot is one of 41 paired patients',
       x = 'ZBTB16 positive fraction', y = 'Paired balance difference') +
  theme_plate + theme(legend.position = 'bottom')

day0$induction <- factor(day0$induction,
  levels = c('IF', 'responsive'), labels = c('Induction failure', 'Responsive'))
pE <- ggplot(day0, aes(induction, balance, fill = induction,
                       colour = induction)) +
  geom_boxplot(width = .48, outlier.shape = NA, alpha = .3,
               linewidth = .4) +
  geom_jitter(width = .13, height = 0, size = 1.2, alpha = .86) +
  scale_fill_manual(values = c('Induction failure' = orange,
                               Responsive = teal), guide = 'none') +
  scale_colour_manual(values = c('Induction failure' = orange,
                                 Responsive = teal), guide = 'none') +
  labs(title = 'E  Baseline treatment response',
       subtitle = '19 failure; 35 responsive; Welch P = 0.0050',
       x = NULL, y = 'Day0 ZEB balance') + theme_plate +
  theme(axis.text.x = element_text(size = 6.8))

combined <- subset(baseline, cohort == 'combined')
combined$label <- factor(metric_labels[combined$metric],
  levels = rev(c('ZEB1', 'ZEB2', 'LMO2', 'Balance')))
pF <- ggplot(combined, aes(diff_a_minus_b, label)) +
  geom_vline(xintercept = 0, colour = grey, linewidth = .4) +
  geom_segment(aes(x = boot_ci95_lo, xend = boot_ci95_hi, yend = label),
               colour = grey, linewidth = .7) +
  geom_point(size = 2.0, colour = purple) +
  geom_text(aes(x = boot_ci95_hi + .12,
                label = sprintf('g=%+.2f', hedges_g_a_minus_b)),
            hjust = 0, size = 2.35, colour = ink, family = 'Arial') +
  scale_x_continuous(limits = c(-2.5, 5.4)) +
  labs(title = 'F  Four baseline effects',
       subtitle = 'Mean difference (CI); Hedges g labels',
       x = 'Mean difference', y = NULL) + theme_plate

cohort <- subset(baseline, metric == 'balance')
cohort$cohort <- factor(cohort$cohort,
  levels = c('combined', 'discovery', 'extension'))
pG <- ggplot(cohort, aes(diff_a_minus_b, cohort)) +
  geom_vline(xintercept = 0, colour = grey, linewidth = .4) +
  geom_segment(aes(x = boot_ci95_lo, xend = boot_ci95_hi, yend = cohort),
               colour = grey, linewidth = .7) +
  geom_point(aes(colour = cohort), size = 2.1) +
  scale_colour_manual(values = c(combined = purple, discovery = teal,
                                 extension = orange), guide = 'none') +
  labs(title = 'G  Discovery and extension',
       subtitle = 'Balance difference; patient bootstrap 95% CI',
       x = 'IF - responsive', y = NULL) + theme_plate

adj <- subset(adjust, y == 'balance')
adj$model <- factor(adj$model,
  levels = rev(c('unadjusted', 'plus_subtype_l1',
                 'plus_age_bin', 'plus_etp')),
  labels = rev(c('Unadjusted', '+ subtype', '+ age', '+ ETP')))
pH <- ggplot(adj, aes(coef_refractory, model)) +
  geom_vline(xintercept = 0, colour = grey, linewidth = .4) +
  geom_segment(aes(x = ci95_lo, xend = ci95_hi, yend = model),
               colour = grey, linewidth = .7) +
  geom_point(size = 2.0, colour = purple) +
  labs(title = 'H  Adjustment sensitivity',
       subtitle = 'Separate one-covariate sensitivity models',
       x = 'IF coefficient; 95% CI', y = NULL) + theme_plate

comp <- bind_rows(
  day0 %>% count(induction, category = etp_status) %>%
    mutate(type = 'ETP status'),
  day0 %>% count(induction, category = subtype_l1) %>%
    mutate(type = 'Broad subtype'))
comp$type <- factor(comp$type, levels = c('ETP status', 'Broad subtype'))
comp_cols <- c(ETP = teal, nonETP = '#a8b5bb', HOXA = purple,
               LMO = red, TAL = orange, TLX = blue, other = '#b4a381')
pI <- ggplot(comp, aes(induction, n, fill = category)) +
  geom_col(position = 'fill', width = .65, colour = 'white', linewidth = .25) +
  facet_wrap(~type, nrow = 1) +
  scale_fill_manual(values = comp_cols, name = NULL) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  labs(title = 'I  Diagnostic composition',
       subtitle = 'ETP and broad subtype groups',
       x = NULL, y = 'Patients') + theme_plate +
  theme(axis.text.x = element_text(angle = 30, hjust = 1, size = 6.4),
        strip.text = element_text(size = 7),
        legend.position = 'bottom', legend.key.size = grid::unit(2.2, 'mm')) +
  guides(fill = guide_legend(nrow = 2))

full <- (pA | pB | pC) / (pD | pE | pF) / (pG | pH | pI) +
  plot_layout(heights = c(1, 1, 1.05))
base <- file.path(out, 'F6')
ggsave(paste0(base, '.pdf'), full, width = 7.09, height = 8.35,
       units = 'in', device = cairo_pdf, bg = 'white', limitsize = FALSE)
write.table(data.frame(panel = LETTERS[1:9],
  source = c('q3_state_day0_pairs.tsv', 'q3_state_day0_pairs.tsv',
             'q3_state.tsv', 'q3_state_day0_pairs.tsv',
             'eligible_patient_timepoint.tsv + q1_baseline.tsv',
             'q1_baseline.tsv', 'q1_baseline.tsv',
             'q1_sensitivity_ols.tsv', 'eligible_patient_timepoint.tsv')),
  paste0(base, '_panel_manifest.tsv'), sep = '\t', quote = FALSE,
  row.names = FALSE)
cat('F6 Gate A: 54 baseline patients, 41 paired states, 9 panels\n')

# Supplementary Figure S3: Lim eligibility, technical decisions and Day28 boundary.
suppressPackageStartupMessages({
  library(ggplot2); library(patchwork); library(readr); library(dplyr)
})
argv <- grep('^--file=', commandArgs(FALSE), value = TRUE)
script <- normalizePath(sub('^--file=', '', argv), winslash = '/')
root <- normalizePath(file.path(dirname(script), '../../..'), winslash = '/')
src <- file.path(root, 'total/data/validation/gateA_lim2025')
out <- file.path(root, 'total/rebuild_2026/edition_20260925/figures/supplementary')
dir.create(out, recursive = TRUE, showWarnings = FALSE)
rd <- function(n) as.data.frame(read_tsv(file.path(src, n),
  show_col_types = FALSE, progress = FALSE))
inv <- rd('inventory.tsv')
all <- rd('pseudobulk_all_timepoints.tsv')
eligible <- rd('eligible_patient_timepoint.tsv')
states <- rd('q3_state_day0_pairs.tsv')
long <- rd('q2_pairs.tsv')
long_stats <- rd('q2_longitudinal.tsv')
ref <- rd('balance_reference.tsv')
stopifnot(nrow(all) == 81L, nrow(eligible) == 66L,
          nrow(states) == 41L, nrow(long) == 12L)
all$balance_frozen <-
  (all$ZEB1_logcpm - ref$logcpm_mean[ref$gene == 'ZEB1']) /
    ref$logcpm_sd[ref$gene == 'ZEB1'] -
  (all$ZEB2_logcpm - ref$logcpm_mean[ref$gene == 'ZEB2']) /
    ref$logcpm_sd[ref$gene == 'ZEB2']
match_idx <- match(eligible$biological_sample_id, all$biological_sample_id)
stopifnot(!anyNA(match_idx),
          max(abs(eligible$balance - all$balance_frozen[match_idx])) < 1e-8)

ink <- '#203744'; teal <- '#168b82'; orange <- '#c78a32'
blue <- '#246b96'; purple <- '#776298'; grey <- '#9eb0b8'
theme_plate <- theme_classic(base_size = 8, base_family = 'Arial') +
  theme(plot.title = element_text(size = 9, face = 'bold', colour = ink),
        plot.subtitle = element_text(size = 7.3, colour = ink),
        axis.title = element_text(size = 8, colour = ink),
        axis.text = element_text(size = 7.3, colour = ink),
        legend.text = element_text(size = 7.1, colour = ink),
        panel.grid.major.y = element_line(colour = '#e5edf0', linewidth = .28),
        plot.margin = margin(4, 6, 4, 6))
cols <- c(IF = orange, responsive = teal)

units <- data.frame(
  label = factor(c('Metadata patients', 'Malignant biological samples',
                   'Eligible Day0 samples', 'Eligible state pairs',
                   'Day0/Day28 patient pairs'),
                 levels = rev(c('Metadata patients', 'Malignant biological samples',
                                'Eligible Day0 samples', 'Eligible state pairs',
                                'Day0/Day28 patient pairs'))),
  n = c(58, 81, 54, 41, 12),
  unit = c('patients', 'samples', 'samples', 'patients', 'patients'))
pA <- ggplot(units, aes(n, label)) +
  geom_segment(aes(x = 0, xend = n, yend = label),
               colour = '#bfd0d5', linewidth = .8) +
  geom_point(size = 2.2, colour = teal) +
  geom_text(aes(label = paste(n, unit)), hjust = -0.1,
            size = 2.5, family = 'Arial', colour = ink) +
  scale_x_continuous(limits = c(0, 104)) +
  labs(title = 'A  Analysis units and eligibility',
       subtitle = 'Counts have different biological denominators',
       x = 'Count', y = NULL) + theme_plate

day0 <- subset(eligible, timepoint == 'Day0')
pB <- ggplot(day0, aes(cohort, n_cells, colour = induction)) +
  geom_point(position = position_jitter(width = .13, height = 0),
             size = 1.5, alpha = .85) +
  geom_hline(yintercept = 50, colour = grey, linetype = 'dashed') +
  scale_y_log10() + scale_colour_manual(values = cols, name = NULL) +
  labs(title = 'B  Day0 malignant blast counts',
       subtitle = '54 eligible patients; dashed: 50-blast threshold',
       x = NULL, y = 'Malignant blasts (log scale)') + theme_plate +
  theme(legend.position = 'bottom')

pC <- ggplot(states, aes(n_cells_neg, n_cells_pos,
                         colour = induction_pos)) +
  geom_vline(xintercept = 20, colour = grey, linetype = 'dashed') +
  geom_hline(yintercept = 20, colour = grey, linetype = 'dashed') +
  geom_point(size = 1.45, alpha = .85) +
  scale_x_log10() + scale_y_log10() +
  scale_colour_manual(values = cols, guide = 'none') +
  labs(title = 'C  Paired-state cell counts',
       subtitle = 'Both states require at least 20 cells',
       x = 'ZBTB16 negative cells', y = 'ZBTB16 positive cells') + theme_plate

decisions <- data.frame(
  item = factor(c('Matrix', 'L086', 'P058', 'Scale'),
                levels = rev(c('Matrix', 'L086', 'P058', 'Scale'))),
  decision = c('Integer UMI in raw/X', 'R1 + R2 summed before logCPM',
               'Overlapping P058 h5ad omitted', 'log2(CPM + 1)'))
pD <- ggplot(decisions, aes(y = item)) +
  geom_text(aes(x = 0, label = decision), hjust = 0,
            size = 2.5, colour = ink, family = 'Arial') +
  scale_x_continuous(limits = c(0, 3.0)) +
  labs(title = 'D  Technical handling',
       subtitle = 'Frozen decisions in Gate A inventory', x = NULL, y = NULL) +
  theme_plate + theme(axis.line = element_blank(), axis.ticks = element_blank(),
                      axis.text.x = element_blank(), panel.grid = element_blank())

long$induction <- factor(long$induction, levels = c('IF', 'responsive'))
long_line <- rbind(
  data.frame(patient = long$patient_id, induction = long$induction,
             time = 'Day0', balance = long$day0_balance),
  data.frame(patient = long$patient_id, induction = long$induction,
             time = 'Day28', balance = long$day28_balance))
long_line$time <- factor(long_line$time, levels = c('Day0', 'Day28'))
pE <- ggplot(long_line, aes(time, balance, group = patient,
                            colour = induction)) +
  geom_line(alpha = .72, linewidth = .45) +
  geom_point(size = 1.45) +
  scale_colour_manual(values = cols, name = NULL) +
  labs(title = 'E  Day0 to Day28 trajectories',
       subtitle = '12 patients; 9 IF and 3 responsive',
       x = NULL, y = 'Frozen ZEB balance') + theme_plate +
  theme(legend.position = 'bottom')

loo <- data.frame(omitted = long$patient_id,
                  mean_delta = vapply(seq_len(nrow(long)), function(i)
                    mean(long$d_balance[-i]), numeric(1)))
loo$omitted <- factor(loo$omitted,
  levels = loo$omitted[order(loo$mean_delta)])
pF <- ggplot(loo, aes(mean_delta, omitted)) +
  geom_vline(xintercept = 0, colour = grey, linewidth = .4) +
  geom_point(size = 1.8, colour = purple) +
  labs(title = 'F  Day28 leave-one-out',
       subtitle = 'Mean paired change after omitting one patient',
       x = 'Mean Day28 - Day0', y = 'Omitted patient') + theme_plate

thresholds <- do.call(rbind, lapply(c(30, 50, 100), function(t) {
  d <- subset(all, timepoint == 'Day0' & n_cells >= t)
  data.frame(threshold = t, n = nrow(d),
             diff = mean(d$balance_frozen[d$induction == 'IF']) -
               mean(d$balance_frozen[d$induction == 'responsive']))
}))
pG <- ggplot(thresholds, aes(factor(threshold), diff)) +
  geom_hline(yintercept = 0, colour = grey, linewidth = .4) +
  geom_col(width = .55, fill = teal) +
  geom_text(aes(label = paste0('n=', n)), vjust = -0.5,
            size = 2.5, colour = ink, family = 'Arial') +
  labs(title = 'G  Blast threshold sensitivity',
       subtitle = 'Frozen balance reference; 30/50/100 blasts',
       x = 'Minimum malignant blasts', y = 'IF - responsive mean') + theme_plate

eligible$timepoint <- factor(eligible$timepoint,
                             levels = c('Day0', 'Day28'))
pH <- ggplot(eligible, aes(lib_size, balance, colour = timepoint)) +
  geom_point(size = 1.45, alpha = .82) +
  scale_x_log10() +
  scale_colour_manual(values = c(Day0 = teal, Day28 = orange), name = NULL) +
  labs(title = 'H  Library-size boundary',
       subtitle = '66 eligible biological samples; descriptive QC',
       x = 'Library UMI count (log scale)', y = 'Frozen ZEB balance') +
  theme_plate + theme(legend.position = 'bottom')

day28_stats <- subset(long_stats, group == 'all_paired')
metric_label <- c(ZEB1_logcpm = 'ZEB1', ZEB2_logcpm = 'ZEB2',
                  LMO2_logcpm = 'LMO2', balance = 'Balance')
day28_stats$label <- factor(metric_label[day28_stats$metric],
  levels = rev(c('ZEB1', 'ZEB2', 'LMO2', 'Balance')))
pI <- ggplot(day28_stats, aes(mean, label)) +
  geom_vline(xintercept = 0, colour = grey, linewidth = .4) +
  geom_segment(aes(x = boot_ci95_lo, xend = boot_ci95_hi, yend = label),
               linewidth = .75, colour = grey) +
  geom_point(size = 2, colour = purple) +
  labs(title = 'I  Day28 paired-effect boundary',
       subtitle = '12 patients; mean Day28 - Day0; bootstrap 95% CI',
       x = 'Paired difference', y = NULL) + theme_plate

page1 <- (pA | pB) / (pC | pD) + plot_layout(heights = c(1, 1))
page2 <- (pE | pF) / (pG | pH) / pI +
  plot_layout(heights = c(1, 1, .75))
ggsave(file.path(out, 'S3_p1.pdf'), page1, width = 7.09, height = 7.85,
       units = 'in', device = cairo_pdf, bg = 'white', limitsize = FALSE)
ggsave(file.path(out, 'S3_p2.pdf'), page2, width = 7.09, height = 8.2,
       units = 'in', device = cairo_pdf, bg = 'white', limitsize = FALSE)
write.table(thresholds, file.path(out, 'S3_threshold_sensitivity.tsv'),
            sep = '\t', quote = FALSE, row.names = FALSE)
write.table(loo, file.path(out, 'S3_day28_leave_one_out.tsv'),
            sep = '\t', quote = FALSE, row.names = FALSE)
cat('S3 audit: A-I across two native-size pages\n')

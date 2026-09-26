# Figure 6: frozen HiChIP support and cross-assay boundaries.
suppressPackageStartupMessages({library(ggplot2); library(patchwork); library(ragg)})
argv <- grep("^--file=", commandArgs(FALSE), value = TRUE)
script <- normalizePath(sub("^--file=", "", argv), winslash = "/")
root <- normalizePath(file.path(dirname(script), "../../.."), winslash = "/")
valid <- file.path(root, "total/data/validation")
out <- file.path(root, "total/rebuild_2026/edition_20260925/figures/main")
sdir <- file.path(root, "total/rebuild_2026/data/source_data_rebuilt/F6")
dir.create(out, recursive = TRUE, showWarnings = FALSE)
dir.create(sdir, recursive = TRUE, showWarnings = FALSE)
rd <- function(x) as.data.frame(readr::read_tsv(file.path(valid, x),
  locale = readr::locale(encoding = "UTF-8"), show_col_types = FALSE,
  progress = FALSE, name_repair = "minimal"))
h <- rd("zeb_chromatin_reciprocal/hichip_zeb1_zeb2_contacts.tsv")
hstat <- rd("zeb_chromatin_reciprocal/hichip_etp_vs_nontall.tsv")
atac <- rd("zeb_chromatin_reciprocal/scatac_zeb1_zeb2_peaks.tsv")
loop <- rd("gse146901/zeb_merged_loop_rates.tsv")
stopifnot(nrow(h) == 10L, nrow(atac) == 7L, nrow(loop) == 3L)

ink <- "#152b39"; zeb1 <- "#bc3446"; zeb2 <- "#246b96"
teal <- "#168b82"; gold <- "#c28830"; purple <- "#776298"
font <- "Arial"
theme_pub <- theme_classic(base_family = font, base_size = 11.5) +
  theme(plot.title = element_text(colour = ink, face = "bold", size = 12.5,
                                  margin = margin(b = 5)),
        plot.subtitle = element_text(colour = ink, size = 9.2),
        plot.caption = element_text(colour = ink, size = 9),
        axis.title = element_text(colour = ink, size = 10),
        axis.text = element_text(colour = ink, size = 9.2),
        legend.title = element_text(colour = ink, size = 9),
        legend.text = element_text(colour = ink, size = 9),
        plot.margin = margin(5, 5, 5, 5),
        panel.grid.major.y = element_line(colour = "#e6edef", linewidth = .3))

h$group <- factor(h$group, levels = c("non-ETP T-ALL", "ETP", "CD34", "thymus"))
gcols <- c("non-ETP T-ALL" = gold, "ETP" = teal,
           "CD34" = "#7196aa", "thymus" = purple)
set.seed(2026)
pA <- ggplot(h, aes(group, n_contacts / 1e6)) +
  geom_point(aes(colour = group), position = position_jitter(width = .10),
             size = 2.5, alpha = .9) +
  scale_colour_manual(values = gcols, guide = "none") +
  labs(title = "A  HiChIP sample context",
       subtitle = "4 non-ETP, 3 ETP, 2 CD34, 1 thymus",
       x = NULL, y = "Filtered contacts (millions)") +
  theme_pub + theme(axis.text.x = element_text(angle = 23, hjust = 1))

pB <- ggplot(h, aes(group, ZEB1_per_M)) +
  geom_point(aes(colour = group), position = position_jitter(width = .10),
             size = 2.6, alpha = .9) +
  scale_colour_manual(values = gcols, guide = "none") +
  labs(title = "B  ZEB1 locus contacts",
       subtitle = "Sample-level normalized contacts",
       x = NULL, y = "ZEB1 contacts / million") +
  theme_pub + theme(axis.text.x = element_text(angle = 23, hjust = 1))

pC <- ggplot(h, aes(group, ZEB2_over_ZEB1)) +
  geom_hline(yintercept = 1, colour = "#b4c2c8", linewidth = .4) +
  geom_point(aes(colour = group), position = position_jitter(width = .10),
             size = 2.6, alpha = .9) +
  scale_colour_manual(values = gcols, guide = "none") +
  labs(title = "C  Relative ZEB locus organization",
       subtitle = "Higher ratios in ETP, CD34 and thymus",
       x = NULL, y = "ZEB2 / ZEB1 contacts") +
  theme_pub + theme(axis.text.x = element_text(angle = 23, hjust = 1))

leuk <- h[h$group %in% c("non-ETP T-ALL", "ETP"), ]
leuk <- leuk[order(leuk$ZEB2_over_ZEB1), ]
leuk$sample <- factor(leuk$sample, levels = leuk$sample)
pD <- ggplot(leuk, aes(sample, ZEB2_over_ZEB1, colour = group)) +
  geom_hline(yintercept = 1, colour = "#a9bac1", linewidth = .4) +
  geom_point(size = 3.0) +
  geom_segment(aes(xend = sample, y = 1, yend = ZEB2_over_ZEB1),
               colour = "#bacbd0", linewidth = .65) +
  scale_colour_manual(values = gcols, guide = "none") +
  labs(title = "D  Complete rank separation in this series",
       subtitle = "3 ETP vs 4 non-ETP; exact two-sided MW P = 0.057",
       x = "HiChIP sample, ordered by ratio", y = "ZEB2 / ZEB1") +
  theme_pub + theme(axis.text.x = element_text(angle = 35, hjust = 1))

h_long <- rbind(data.frame(sample = h$sample, group = h$group,
                           gene = "ZEB1", value = h$ZEB1_per_M),
                data.frame(sample = h$sample, group = h$group,
                           gene = "ZEB2", value = h$ZEB2_per_M))
h_long$gene <- factor(h_long$gene, levels = c("ZEB1", "ZEB2"))
pE <- ggplot(h_long, aes(gene, value, group = sample)) +
  geom_line(aes(colour = group), linewidth = .65, alpha = .55) +
  geom_point(aes(fill = gene), shape = 21, colour = "white", stroke = .2,
             size = 2.5) +
  scale_colour_manual(values = gcols, guide = "none") +
  scale_fill_manual(values = c(ZEB1 = zeb1, ZEB2 = zeb2), guide = "none") +
  labs(title = "E  Paired ZEB contact profiles",
       subtitle = "Both loci measured in every HiChIP sample",
       x = NULL, y = "Contacts / million") +
  theme_pub

loop_long <- rbind(
 data.frame(group = loop$group, feature = "ZEB1 body", rate = loop$ZEB1_body_per_1000),
 data.frame(group = loop$group, feature = "ZEB2 body", rate = loop$ZEB2_body_per_1000),
 data.frame(group = loop$group, feature = "ZEB1 +/-1 Mb", rate = loop$ZEB1_locus_1Mb_per_1000),
 data.frame(group = loop$group, feature = "ZEB2 +/-1 Mb", rate = loop$ZEB2_locus_1Mb_per_1000))
loop_long$z <- ave(loop_long$rate, loop_long$feature,
  FUN = function(x) if (sd(x) == 0) rep(0, length(x)) else as.numeric(scale(x)))
loop_long$feature <- factor(loop_long$feature,
  levels = c("ZEB1 body", "ZEB2 body", "ZEB1 +/-1 Mb", "ZEB2 +/-1 Mb"))
loop_long$group <- factor(loop_long$group, levels = c("ETP", "non-ETP", "normal T"))
pF <- ggplot(loop_long, aes(feature, group, fill = z)) +
  geom_tile(colour = "white", linewidth = .7) +
  geom_text(aes(label = sprintf("%.2f", rate)), colour = ink,
            size = 3.0, fontface = "bold", family = font) +
  scale_fill_gradient2(low = "#7094ad", mid = "#f6f7f5", high = "#c77d73",
                       midpoint = 0, limits = c(-1.6, 1.6), guide = "none") +
  labs(title = "F  Pooled loops: no polarity",
       subtitle = "GSE146901; calls per 1,000 loops",
       x = NULL, y = NULL) +
  theme_pub + theme(panel.grid = element_blank(), axis.line = element_blank(),
                    axis.ticks = element_blank(),
                    axis.text.x = element_text(angle = 30, hjust = 1))

atac$group <- factor(atac$group, levels = c("non-ETP T-ALL", "ETP"))
pG <- ggplot(atac, aes(ZEB1_per_10k_peaks, ZEB2_per_10k_peaks,
                        colour = group)) +
  geom_abline(slope = 1, intercept = 0, colour = "#a8bac1", linetype = "dashed") +
  geom_point(size = 2.8) +
  scale_colour_manual(values = gcols, name = NULL) +
  labs(title = "G  scATAC is not an independent replication",
       subtitle = "Seven peak sets / distinct assay and normalization",
       x = "ZEB1 peaks / 10,000", y = "ZEB2 peaks / 10,000") +
  theme_pub + theme(legend.position = "bottom", legend.key.width = grid::unit(7, "pt"))

prom <- rbind(data.frame(sample = h$sample, group = h$group,
                         gene = "ZEB1", value = h$ZEB1_prom / h$n_contacts * 1e6),
              data.frame(sample = h$sample, group = h$group,
                         gene = "ZEB2", value = h$ZEB2_prom / h$n_contacts * 1e6))
prom$gene <- factor(prom$gene, levels = c("ZEB1", "ZEB2"))
pH <- ggplot(prom, aes(gene, value, group = sample)) +
  geom_line(aes(colour = group), linewidth = .65, alpha = .5) +
  geom_point(aes(fill = gene), shape = 21, colour = "white", stroke = .2,
             size = 2.6) +
  scale_colour_manual(values = gcols, guide = "none") +
  scale_fill_manual(values = c(ZEB1 = zeb1, ZEB2 = zeb2), guide = "none") +
  labs(title = "H  Promoter-window contacts",
       subtitle = "Same ten samples; depth-normalized promoter counts",
       x = NULL, y = "Promoter contacts / million") + theme_pub

top <- wrap_plots(pA, pB, pC, ncol = 3)
middle <- wrap_plots(pD, pE, pF, ncol = 3, widths = c(1, .8, 1.2))
bottom <- wrap_plots(pG, pH, ncol = 2, widths = c(.9, 1.3))
full <- wrap_plots(top, middle, bottom, ncol = 1, heights = c(.95, 1, .95))
wt <- function(x, name) write.table(x, file.path(sdir, name), sep = "\t",
  quote = FALSE, row.names = FALSE, na = "NA")
wt(h, "F6A-E_HiChIP_samples.tsv")
wt(hstat, "F6B-D_HiChIP_frozen_tests.tsv")
wt(loop_long, "F6F_GSE146901_pooled_loops.tsv")
wt(atac, "F6G_scATAC_peak_sets.tsv")
wt(prom, "F6H_promoter_contacts.tsv")
base <- file.path(out, "F6")
ggsave(paste0(base, ".pdf"), full, width = 14.4, height = 11.4,
       units = "in", device = cairo_pdf, bg = "white", limitsize = FALSE)
ggsave(paste0(base, ".svg"), full, width = 14.4, height = 11.4,
       units = "in", device = svg, bg = "white", limitsize = FALSE)
ggsave(paste0(base, ".png"), full, width = 14.4, height = 11.4,
       units = "in", device = agg_png, dpi = 300, bg = "white", limitsize = FALSE)
ggsave(paste0(base, ".tiff"), full, width = 14.4, height = 11.4,
       units = "in", device = agg_tiff, dpi = 300, bg = "white", limitsize = FALSE)
writeLines(capture.output(sessionInfo()), paste0(base, "_R_sessionInfo.txt"))
manifest <- data.frame(
  figure = "F6", panel = LETTERS[1:8],
  source = c(rep("zeb_chromatin_reciprocal/hichip_zeb1_zeb2_contacts.tsv", 5),
    "gse146901/zeb_merged_loop_rates.tsv",
    "zeb_chromatin_reciprocal/scatac_zeb1_zeb2_peaks.tsv",
    "zeb_chromatin_reciprocal/hichip_zeb1_zeb2_contacts.tsv"),
  displayed_data = c(rep("F6A-E_HiChIP_samples.tsv", 5),
    "F6F_GSE146901_pooled_loops.tsv", "F6G_scATAC_peak_sets.tsv",
    "F6H_promoter_contacts.tsv"),
  display_transform = c("Sample depth and group annotation, frozen values",
    "Sample-level ZEB1 contacts per million, frozen values",
    "Sample-level ZEB2/ZEB1 contacts, frozen values",
    "Seven leukemia samples ordered by frozen contact ratio; exact P from frozen table",
    "Paired frozen ZEB1 and ZEB2 contacts per sample",
    "Frozen pooled-loop rates; within-feature display z and printed rates",
    "Seven frozen scATAC peak sets, sample-level rates",
    "Depth-normalized promoter contacts from frozen sample-level counts"),
  script = "total/rebuild_2026/code/build_f6_ggplot.R")
write.table(manifest, paste0(base, "_panel_manifest.tsv"), sep = "\t",
            quote = FALSE, row.names = FALSE)
cat(base, "\n")

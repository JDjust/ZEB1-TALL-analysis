# Module 6 figures M6_6.1-6.25. One numbered file per checklist ID.
# Recipes: code2 346 bars, 238 bar+points, 88 forest, 193 scatter.
# Biological n for HiChIP is sample, not contact.

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

root <- "/data-b/liangfuhua/projects/ZEB1_TALL_analysis/module6"
m3proc <- "/data-b/liangfuhua/projects/ZEB1_TALL_analysis/module3/processed"
m3tab <- "/data-b/liangfuhua/projects/ZEB1_TALL_analysis/module3/tables"
fig <- file.path(root, "figures"); tab <- file.path(root, "tables")
dir.create(fig, FALSE, TRUE)

theme_sci <- function(base_size = 11) {
  theme_classic(base_size = base_size) +
    theme(axis.text = element_text(color = "black"), axis.title = element_text(color = "black"),
          axis.line = element_line(linewidth = 0.45, color = "black"),
          plot.title = element_text(face = "bold", size = base_size + 1),
          plot.subtitle = element_text(size = base_size - 1, color = "grey25"),
          strip.background = element_blank())
}
save_both <- function(plot, file, width, height) {
  tryCatch({
    pdf(file, width = width, height = height, useDingbats = FALSE); print(plot); dev.off()
    png(sub("\\.pdf$", ".png", file), width = width * 160, height = height * 160, res = 160)
    print(plot); dev.off()
    message("saved ", basename(file))
  }, error = function(e) message("FAIL ", basename(file), ": ", conditionMessage(e)))
}
run <- function(label, fn) tryCatch({ fn(); message("OK ", label) },
                                    error = function(e) message("FAIL ", label, ": ", conditionMessage(e)))
rd <- function(f) {
  p <- if (file.exists(f)) f else file.path(tab, f)
  if (!file.exists(p)) return(data.table())
  fread(p)
}
note_fig <- function(id, text, fn, width = 8.2, height = 3.5) {
  p <- ggplot() +
    annotate("text", x = 0, y = 0, label = paste0(id, "\n", text), size = 4.0, lineheight = 1.15) +
    theme_void() + xlim(-1, 1) + ylim(-1, 1)
  save_both(p, file.path(fig, fn), width, height)
}

pk <- rd(file.path(tab, "M6_6.1_6.12_ZEB1_peak_overlap.tsv"))
ct <- rd(file.path(tab, "M6_6.13_6.16_ZEB1_contacts.tsv"))
cimp <- rd(file.path(m3proc, "gse272023_cimp.tsv"))
pal_g <- c(TALL = "#E64B35", ETP = "#4DBBD5", CD34 = "#00A087", THY = "#3C5488")
pal_m <- c(ATAC = "#E64B35", H3K4me1 = "#4DBBD5", H3K4me3 = "#00A087")

if (nrow(pk)) {
  pk[, distal := pmax(ZEB1_genebody - ZEB1_promoter, 0)]
  pk_long <- melt(pk, id.vars = c("file", "mark", "n_peaks"),
                  measure.vars = c("ZEB1_promoter", "ZEB1_genebody", "distal"),
                  variable.name = "window", value.name = "n")
}
if (nrow(ct)) {
  ct[, ZEB1_per_M := as.numeric(ZEB1_per_M)]
  ct[, group := factor(group, levels = c("TALL", "ETP", "CD34", "THY"))]
}
if (nrow(cimp)) {
  d_cimp <- copy(cimp)
  if (!"cimp_g" %in% names(d_cimp)) {
    cand <- names(d_cimp)[grepl("cimp", names(d_cimp), ignore.case = TRUE)]
    if (length(cand)) d_cimp[, cimp_g := d_cimp[[cand[1]]]]
  }
  d_cimp[, cimp_g := fifelse(grepl("high", cimp_g, ignore.case = TRUE), "CIMP-high",
                             fifelse(grepl("low", cimp_g, ignore.case = TRUE), "CIMP-low", NA_character_))]
  d_cimp <- d_cimp[cimp_g %in% c("CIMP-high", "CIMP-low") & is.finite(as.numeric(ZEB1))]
  d_cimp[, ZEB1 := as.numeric(ZEB1)]
} else {
  d_cimp <- data.table()
}

peak_mark <- function(mark_name, id, title, fn) {
  if (!nrow(pk)) stop("no peaks")
  d <- pk[mark == mark_name]
  if (!nrow(d)) stop(paste("no", mark_name))
  long <- melt(d, id.vars = "mark",
               measure.vars = c("ZEB1_promoter", "ZEB1_genebody"),
               variable.name = "window", value.name = "n")
  long[, window := fifelse(window == "ZEB1_promoter", "promoter +/- 2 kb", "gene body")]
  p <- ggplot(long, aes(window, n, fill = window)) +
    geom_col(width = 0.55, color = "white") +
    geom_text(aes(label = n), vjust = -0.3, size = 4) +
    scale_fill_manual(values = c("promoter +/- 2 kb" = "#E64B35", "gene body" = "#3C5488"), guide = "none") +
    labs(title = title,
         subtitle = "ALL-SIL (GSE110637); one cell line, not a patient cohort",
         x = NULL, y = "Overlapping peaks") + theme_sci()
  save_both(p, file.path(fig, fn), 6.4, 5.2)
}

run("6.1 TSS/promoter ATAC", function() {
  if (!nrow(pk)) stop("no peaks")
  d <- pk[mark == "ATAC"]
  p <- ggplot(d, aes(x = "ATAC", y = ZEB1_promoter)) +
    geom_col(width = 0.45, fill = "#E64B35", color = "white") +
    geom_text(aes(label = ZEB1_promoter), vjust = -0.3, size = 4.5) +
    labs(title = "6.1  ZEB1 TSS / promoter accessibility",
         subtitle = sprintf("ALL-SIL ATAC: %d promoter peaks / %d total", d$ZEB1_promoter[1], d$n_peaks[1]),
         x = NULL, y = "Peaks overlapping ZEB1 promoter +/- 2 kb") + theme_sci()
  save_both(p, file.path(fig, "M6_6.1_ZEB1_TSS_ATAC.pdf"), 5.6, 5.2)
})

run("6.2 promoter peaks", function() {
  if (!nrow(pk)) stop("no peaks")
  p <- ggplot(pk, aes(mark, ZEB1_promoter, fill = mark)) +
    geom_col(width = 0.6, color = "white") +
    geom_text(aes(label = ZEB1_promoter), vjust = -0.3, size = 4) +
    scale_fill_manual(values = pal_m, guide = "none") +
    labs(title = "6.2  ZEB1 promoter peaks by mark",
         subtitle = "ALL-SIL GSE110637", x = NULL, y = "Promoter-overlapping peaks") + theme_sci()
  save_both(p, file.path(fig, "M6_6.2_promoter_peaks.pdf"), 6.2, 5.2)
})

run("6.3 distal peaks", function() {
  if (!nrow(pk)) stop("no peaks")
  p <- ggplot(pk, aes(mark, distal, fill = mark)) +
    geom_col(width = 0.6, color = "white") +
    geom_text(aes(label = distal), vjust = -0.3, size = 4) +
    scale_fill_manual(values = pal_m, guide = "none") +
    labs(title = "6.3  Distal / gene-body peaks at ZEB1",
         subtitle = "Gene-body minus promoter as distal proxy (no enhancer annotation file)",
         x = NULL, y = "Peaks") + theme_sci()
  save_both(p, file.path(fig, "M6_6.3_distal_peaks.pdf"), 6.2, 5.2)
})

run("6.4 ETP vs non-ETP chromatin", function() {
  if (nrow(ct) < 4) stop(sprintf("contacts n=%d", nrow(ct)))
  d <- ct[group %in% c("TALL", "ETP")]
  wt <- tryCatch(wilcox.test(ZEB1_per_M ~ group, data = d)$p.value, error = function(e) NA_real_)
  p <- ggplot(d, aes(group, ZEB1_per_M, fill = group)) +
    geom_boxplot(width = 0.55, outlier.shape = NA, alpha = 0.9) +
    geom_point(position = position_jitter(width = 0.08, seed = 1), size = 2.6) +
    scale_fill_manual(values = pal_g, guide = "none") +
    labs(title = "6.4  ETP vs non-ETP ZEB1 HiChIP contact rate",
         subtitle = sprintf("GSE243915 H3K27ac HiChIP; Wilcoxon P=%.4g (n_TALL=%d, n_ETP=%d)",
                            wt, d[group == "TALL", .N], d[group == "ETP", .N]),
         x = NULL, y = "ZEB1-locus contacts / million") + theme_sci()
  save_both(p, file.path(fig, "M6_6.4_ETP_vs_TALL_HiChIP.pdf"), 6.4, 5.6)
})

run("6.5 high vs low ZEB1 ATAC substitute", function() {
  note_fig("6.5",
           "No patient ATAC paired with ZEB1 RNA in the download inventory.\nSubstitute: GSE243915 H3K27ac HiChIP contact rate by group (6.4 / 6.16).\nGSE110637 ATAC is a single ALL-SIL bed, so high vs low ZEB1 is undefined.",
           "M6_6.5_no_patient_ATAC.pdf")
})

run("6.6 peak-to-gene", function() {
  note_fig("6.6",
           "Peak-to-gene correlation needs a multi-sample ATAC count matrix plus matched RNA.\nNot available. Closest substitute is HiChIP ZEB1 vs LMO2/GATA3 locus counts (6.15).",
           "M6_6.6_no_peak_to_gene.pdf")
})

run("6.7 chromVAR", function() {
  note_fig("6.7",
           "chromVAR needs an ATAC count matrix across samples.\nGSE110637 ATAC is one ALL-SIL peak bed. Motif activity not computed.",
           "M6_6.7_no_chromVAR.pdf")
})

run("6.8 motifs", function() {
  note_fig("6.8",
           "LMO2 / TAL1 / GATA / RUNX motif activity was not computed (no chromVAR matrix).\nHiChIP LMO2 and GATA3 locus overlapping pairs are shown in 6.15.",
           "M6_6.8_no_motif_activity.pdf")
})

run("6.9 H3K4me3", function() peak_mark("H3K4me3", "6.9", "6.9  H3K4me3 at ZEB1 (ALL-SIL)", "M6_6.9_H3K4me3_ZEB1.pdf"))
run("6.10 H3K4me1", function() peak_mark("H3K4me1", "6.10", "6.10  H3K4me1 at ZEB1 (ALL-SIL)", "M6_6.10_H3K4me1_ZEB1.pdf"))

run("6.11 H3K27ac ChIP substitute", function() {
  note_fig("6.11",
           "No H3K27ac ChIP peak beds in GSE110637 (ATAC + H3K4me1 + H3K4me3 only).\nSubstitute: GSE243915 is H3K27ac HiChIP; contact rate is 6.13-6.16.",
           "M6_6.11_no_H3K27ac_ChIP.pdf")
  if (nrow(ct) >= 4) {
    p <- ggplot(ct, aes(ZEB1_per_M, reorder(paste(group, sample), ZEB1_per_M), fill = group)) +
      geom_col(width = 0.7, color = "white") +
      scale_fill_manual(values = pal_g) +
      labs(title = "6.11  H3K27ac HiChIP at ZEB1 (substitute for ChIP)",
           x = "Contacts / million", y = NULL, fill = NULL) + theme_sci()
    save_both(p, file.path(fig, "M6_6.11_H3K27ac_HiChIP_substitute.pdf"), 8.2, 5.8)
  }
})

run("6.12 enhancer score", function() {
  if (!nrow(pk)) stop("no peaks")
  d <- pk[, .(mark, promoter = ZEB1_promoter, genebody = ZEB1_genebody,
              enhancer_proxy = fifelse(mark == "H3K4me1", ZEB1_genebody, ZEB1_promoter))]
  p <- ggplot(d, aes(mark, enhancer_proxy, fill = mark)) +
    geom_col(width = 0.6, color = "white") +
    geom_text(aes(label = enhancer_proxy), vjust = -0.3, size = 4) +
    scale_fill_manual(values = pal_m, guide = "none") +
    labs(title = "6.12  Active-enhancer proxy at ZEB1",
         subtitle = "H3K4me1 gene-body peaks as enhancer proxy; no H3K27ac ChIP",
         x = NULL, y = "Peaks") + theme_sci()
  save_both(p, file.path(fig, "M6_6.12_enhancer_proxy.pdf"), 6.4, 5.2)
})

run("6.13 promoter-enhancer loops", function() {
  if (nrow(ct) < 4) stop(sprintf("contacts n=%d", nrow(ct)))
  p <- ggplot(ct, aes(ZEB1_per_M, reorder(paste(group, sample), ZEB1_per_M), fill = group)) +
    geom_col(width = 0.7, color = "white") +
    scale_fill_manual(values = pal_g) +
    labs(title = "6.13  HiChIP contacts overlapping ZEB1 (hg38)",
         subtitle = "GSE243915 filtered pairs; rate per million. THY-2124 gzip stream failed.",
         x = "ZEB1-locus contacts / million", y = NULL, fill = NULL) + theme_sci()
  save_both(p, file.path(fig, "M6_6.13_ZEB1_contact_rate.pdf"), 8.2, 5.8)
})

run("6.14 H3K27ac HiChIP", function() {
  if (nrow(ct) < 4) stop("no contacts")
  agg <- ct[, .(n = .N, mean = mean(ZEB1_per_M, na.rm = TRUE),
                se = sd(ZEB1_per_M, na.rm = TRUE) / sqrt(.N)), by = group]
  p <- ggplot(agg, aes(mean, group, fill = group)) +
    geom_col(width = 0.6, color = "white") +
    geom_errorbar(aes(xmin = mean - se, xmax = mean + se), width = 0.15) +
    geom_point(data = ct, aes(ZEB1_per_M, group), inherit.aes = FALSE,
               position = position_jitter(height = 0.08, seed = 1), size = 2.4) +
    scale_fill_manual(values = pal_g, guide = "none") +
    labs(title = "6.14  H3K27ac HiChIP ZEB1 contact rate by group",
         x = "Contacts / million", y = NULL) + theme_sci()
  save_both(p, file.path(fig, "M6_6.14_H3K27ac_HiChIP_group.pdf"), 7.4, 5.0)
})

run("6.15 loop strength / loci", function() {
  if (nrow(ct) < 4) stop("no contacts")
  long <- melt(ct, id.vars = c("group", "sample"),
               measure.vars = c("ZEB1_locus", "LMO2_locus", "GATA3_locus", "ZEB1_LMO2_pairs"),
               variable.name = "locus", value.name = "n")
  p <- ggplot(long, aes(group, n, fill = group)) +
    geom_boxplot(width = 0.55, outlier.shape = NA, alpha = 0.85) +
    geom_point(position = position_jitter(width = 0.08, seed = 1), size = 2.2) +
    facet_wrap(~locus, scales = "free_y") +
    scale_fill_manual(values = pal_g, guide = "none") +
    labs(title = "6.15  Locus-overlapping HiChIP pairs",
         subtitle = "ZEB1-LMO2 trans pairs are 0 in every finished sample",
         x = NULL, y = "Count") + theme_sci()
  save_both(p, file.path(fig, "M6_6.15_locus_contact_counts.pdf"), 9.0, 6.2)
})

run("6.16 ETP/non-ETP comparison", function() {
  if (nrow(ct) < 4) stop("no contacts")
  d <- ct[group %in% c("TALL", "ETP", "CD34", "THY")]
  p <- ggplot(d, aes(group, ZEB1_per_M, fill = group)) +
    geom_boxplot(width = 0.55, outlier.shape = NA) +
    geom_point(position = position_jitter(width = 0.08, seed = 1), size = 2.5) +
    scale_fill_manual(values = pal_g, guide = "none") +
    labs(title = "6.16  ZEB1 HiChIP by group including CD34 / thymus",
         subtitle = "T-ALL ~2-fold above ETP / CD34 / THY; n=10 (THY-2124 failed)",
         x = NULL, y = "Contacts / million") + theme_sci()
  save_both(p, file.path(fig, "M6_6.16_group_comparison.pdf"), 6.8, 5.6)
})

run("6.17 promoter CpG", function() {
  b <- rd(file.path(tab, "M6_6.17_ZEB1_sample_region_mean.tsv"))
  if (!nrow(b)) {
    note_fig("6.17",
             "GSE69954 450k promoter table missing.",
             "M6_6.17_no_450k_beta.pdf")
    return(invisible(NULL))
  }
  d <- b[region == "promoter" & cimp %in% c("CIMP+", "CIMP-")]
  d[, cimp := factor(cimp, levels = c("CIMP-", "CIMP+"))]
  wt <- tryCatch(wilcox.test(mean_beta ~ cimp, data = d)$p.value, error = function(e) NA_real_)
  p <- ggplot(d, aes(cimp, mean_beta, fill = cimp)) +
    geom_boxplot(width = 0.55, outlier.shape = NA, alpha = 0.9) +
    geom_point(position = position_jitter(width = 0.08, seed = 1), size = 1.6, alpha = 0.75) +
    scale_fill_manual(values = c("CIMP-" = "#4DBBD5", "CIMP+" = "#E64B35"), guide = "none") +
    labs(title = "6.17  ZEB1 promoter methylation (GSE69954)",
         subtitle = sprintf("NOPHO T-ALL n=65; Wilcoxon P=%.3g. Promoter remains hypomethylated.", wt),
         x = NULL, y = "Mean promoter beta") + theme_sci()
  save_both(p, file.path(fig, "M6_6.17_ZEB1_promoter_beta.pdf"), 5.6, 5.4)
})

run("6.18 methylation vs expression", function() {
  if (!nrow(d_cimp)) stop("no CIMP")
  p <- ggplot(d_cimp, aes(cimp_g, ZEB1, fill = cimp_g)) +
    geom_boxplot(width = 0.55, outlier.shape = NA, alpha = 0.9) +
    geom_point(position = position_jitter(width = 0.08, seed = 1), size = 1.4, alpha = 0.7) +
    scale_fill_manual(values = c("CIMP-low" = "#4DBBD5", "CIMP-high" = "#E64B35"), guide = "none") +
    labs(title = "6.18  ZEB1 vs CIMP (methylation surrogate)",
         subtitle = "GSE272023 RNA+CIMP; 450k promoter betas unavailable",
         x = NULL, y = "ZEB1") + theme_sci()
  save_both(p, file.path(fig, "M6_6.18_CIMP_vs_ZEB1.pdf"), 6.4, 5.6)
})

run("6.19 CIMP-high/low counts", function() {
  if (!nrow(d_cimp)) stop("no CIMP")
  sm <- d_cimp[, .N, by = cimp_g]
  p <- ggplot(sm, aes(cimp_g, N, fill = cimp_g)) +
    geom_col(width = 0.55, color = "white") +
    geom_text(aes(label = N), vjust = -0.3, size = 4) +
    scale_fill_manual(values = c("CIMP-low" = "#4DBBD5", "CIMP-high" = "#E64B35"), guide = "none") +
    labs(title = "6.19  CIMP-high / CIMP-low tumors (GSE272023)",
         x = NULL, y = "Tumors") + theme_sci()
  save_both(p, file.path(fig, "M6_6.19_CIMP_counts.pdf"), 6.0, 5.2)
})

run("6.20 CIMP x ZEB1", function() {
  if (!nrow(d_cimp)) stop("no CIMP")
  wt <- tryCatch(wilcox.test(ZEB1 ~ cimp_g, data = d_cimp)$p.value, error = function(e) NA_real_)
  a <- d_cimp[cimp_g == "CIMP-high", ZEB1]
  b <- d_cimp[cimp_g == "CIMP-low", ZEB1]
  p <- ggplot(d_cimp, aes(cimp_g, ZEB1, fill = cimp_g)) +
    geom_violin(width = 0.7, alpha = 0.35, color = NA) +
    geom_boxplot(width = 0.35, outlier.shape = NA, alpha = 0.9) +
    geom_point(position = position_jitter(width = 0.08, seed = 1), size = 1.2, alpha = 0.65) +
    scale_fill_manual(values = c("CIMP-low" = "#4DBBD5", "CIMP-high" = "#E64B35"), guide = "none") +
    labs(title = "6.20  CIMP x ZEB1 (Wilcoxon)",
         subtitle = sprintf("n_high=%d n_low=%d P=%.3g. No MRD table for CIMP x MRD (3.20).",
                            length(a), length(b), wt),
         x = NULL, y = "ZEB1") + theme_sci()
  save_both(p, file.path(fig, "M6_6.20_CIMP_ZEB1_wilcox.pdf"), 6.4, 5.6)
})

run("6.21 RNA-ATAC", function() {
  note_fig("6.21",
           "RNA-ATAC integration is limited to ALL-SIL peaks vs T-ALL bulk RNA from other modules.\nNo matched patient ATAC+RNA matrix. Peak overlap is 6.1-6.3.",
           "M6_6.21_RNA_ATAC_note.pdf")
})

run("6.22 RNA-H3K27ac", function() {
  note_fig("6.22",
           "GSE243915 HiChIP samples are not the TARGET/Pharmacotype RNA patients.\nGroup-level substitute: T-ALL HiChIP ZEB1 contact rate is higher than ETP (6.4),\nwhile bulk ETP RNA is ZEB1-low (Module 3). Not a within-patient correlation.",
           "M6_6.22_RNA_H3K27ac_note.pdf")
})

run("6.23 RNA-methylation", function() {
  if (!nrow(d_cimp)) {
    note_fig("6.23", "CIMP table missing.", "M6_6.23_RNA_methylation_note.pdf")
  } else {
    p <- ggplot(d_cimp, aes(cimp_g, ZEB1, fill = cimp_g)) +
      geom_boxplot(width = 0.5, outlier.shape = NA) +
      geom_point(position = position_jitter(width = 0.08, seed = 1), size = 1.3, alpha = 0.7) +
      scale_fill_manual(values = c("CIMP-low" = "#4DBBD5", "CIMP-high" = "#E64B35"), guide = "none") +
      labs(title = "6.23  RNA vs methylation surrogate (CIMP)",
           x = NULL, y = "ZEB1") + theme_sci()
    save_both(p, file.path(fig, "M6_6.23_RNA_methylation.pdf"), 6.2, 5.4)
  }
})

run("6.24 accessibility-regulon", function() {
  note_fig("6.24",
           "Accessibility-regulon needs chromVAR or a matched ATAC+RNA regulon matrix.\nNot available. TF-activity substitutes are in Module 4 (C3 TFT).",
           "M6_6.24_no_regulon_ATAC.pdf")
})

run("6.25 LMO2-chromatin-ZEB1 mediation", function() {
  if (!nrow(ct)) {
    note_fig("6.25", "No HiChIP table.", "M6_6.25_no_mediation.pdf")
  } else {
    sm <- ct[, .(n_pairs = sum(ZEB1_LMO2_pairs), n_samples = .N), by = group]
    p <- ggplot(sm, aes(group, n_pairs, fill = group)) +
      geom_col(width = 0.6, color = "white") +
      geom_text(aes(label = n_pairs), vjust = -0.3, size = 4) +
      scale_fill_manual(values = pal_g, guide = "none") +
      labs(title = "6.25  ZEB1-LMO2 trans HiChIP pairs",
           subtitle = "Zero pairs in all samples; LMO2-chromatin-ZEB1 mediation was not fit.",
           x = NULL, y = "ZEB1-LMO2 contacts") + theme_sci()
    save_both(p, file.path(fig, "M6_6.25_ZEB1_LMO2_pairs.pdf"), 6.6, 5.4)
  }
})

run("overview peak overlap", function() {
  if (!nrow(pk)) stop("no peaks")
  long <- melt(pk, id.vars = c("file", "mark", "n_peaks"),
               measure.vars = c("ZEB1_promoter", "ZEB1_genebody"),
               variable.name = "window", value.name = "n")
  long[, window := fifelse(window == "ZEB1_promoter", "promoter +/- 2 kb", "gene body")]
  p <- ggplot(long, aes(n, mark, fill = window)) +
    geom_col(width = 0.65, position = position_dodge(0.75), color = "white") +
    geom_text(aes(label = n), position = position_dodge(0.75), hjust = -0.15, size = 3.2) +
    scale_fill_manual(values = c("promoter +/- 2 kb" = "#E64B35", "gene body" = "#3C5488")) +
    labs(title = "6.1-6.12 overview  ZEB1 locus peaks in ALL-SIL",
         x = "Overlapping peaks", y = NULL, fill = NULL) + theme_sci()
  save_both(p, file.path(fig, "M6_6.1_6.12_ZEB1_peak_overlap.pdf"), 8.0, 5.0)
})

fwrite(data.table(
  id = c("6.5", "6.6", "6.7", "6.8", "6.11", "6.21", "6.22", "6.24", "6.25"),
  status = c(
    "No patient ATAC paired with ZEB1 RNA; HiChIP group contrast used.",
    "No multi-sample ATAC count matrix.",
    "chromVAR not run; single ALL-SIL ATAC bed.",
    "Motif activity not computed.",
    "No H3K27ac ChIP beds; H3K27ac HiChIP used.",
    "No matched RNA-ATAC patients.",
    "HiChIP and bulk RNA are different cohorts; group-level only.",
    "No accessibility-regulon matrix.",
    "ZEB1-LMO2 trans pairs = 0; mediation not fit."
  )
), file.path(tab, "M6_missing_or_substitute.tsv"), sep = "\t")

cat("MODULE6 ANALYZE DONE\n")
cat("FIG", length(list.files(fig, pattern = "\\.(pdf|png)$")), "\n")
cat("CONTACTS_N", nrow(ct), "\n")

# Module 6 round-2 figures: ChIP occupancy, H3K27ac, scATAC peaks, mouse ATAC.
# Occupancy is not repression.
suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})
has_cht <- requireNamespace("ComplexHeatmap", quietly = TRUE)
has_circlize <- requireNamespace("circlize", quietly = TRUE)

root <- "/data-b/liangfuhua/projects/ZEB1_TALL_analysis/module6"
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
    print(plot); dev.off(); message("saved ", basename(file))
  }, error = function(e) message("FAIL ", basename(file), ": ", conditionMessage(e)))
}
run <- function(label, fn) tryCatch({ fn(); message("OK ", label) },
                                    error = function(e) message("FAIL ", label, ": ", conditionMessage(e)))
rd <- function(f) {
  p <- file.path(tab, f)
  if (!file.exists(p)) return(data.table())
  fread(p)
}

occ <- rd("M6_r2_GSE154675_occupancy.tsv")
sig <- rd("M6_r2_GSE154675_promoter_signal.tsv")
h3 <- rd("M6_r2_GSE70734_H3K27ac_ZEB1.tsv")
sc <- rd("M6_r2_GSE280728_ZEB1_peaks.tsv")
mm <- rd("M6_r2_GSE287757_Zeb1_peaks.tsv")
pal_ab <- c(LMO2 = "#E64B35", TAL1 = "#4DBBD5", LDB1 = "#00A087", GATA2 = "#3C5488")

run("6r2 occupancy heatmap-like bars", function() {
  if (!nrow(occ)) stop("no occupancy")
  d <- copy(occ)
  d[, prom := fifelse(is.finite(hg38_prom), hg38_prom, hg19_prom)]
  d <- d[is.finite(prom)]
  d[, antibody := factor(antibody, levels = c("LMO2", "TAL1", "LDB1", "GATA2"))]
  p <- ggplot(d, aes(antibody, line, fill = log1p(pmax(prom, 0)))) +
    geom_tile(color = "white", linewidth = 0.4) +
    geom_text(aes(label = sprintf("%.2f", prom)), size = 3) +
    scale_fill_gradient(low = "#F5F5F5", high = "#E64B35") +
    labs(title = "GSE154675 mean bigWig at the ZEB1 promoter",
         subtitle = "4 T-ALL lines × 4 antibodies. Occupancy ≠ repression. No input.",
         x = NULL, y = NULL, fill = "log1p(mean)") + theme_sci()
  save_both(p, file.path(fig, "M6_r2_6.12_occupancy_heatmap.pdf"), 7.2, 5.2)
})

run("6r2 occupancy grouped bars", function() {
  if (!nrow(occ)) stop("no occupancy")
  d <- copy(occ)
  d[, prom := fifelse(is.finite(hg38_prom), hg38_prom, hg19_prom)]
  p <- ggplot(d, aes(line, prom, fill = antibody)) +
    geom_col(position = position_dodge(0.8), width = 0.75, color = "white") +
    scale_fill_manual(values = pal_ab) +
    labs(title = "ZEB1 promoter occupancy by antibody",
         x = NULL, y = "Mean bigWig signal", fill = NULL) + theme_sci()
  save_both(p, file.path(fig, "M6_r2_6.7_occupancy_bars.pdf"), 7.6, 5.2)
})

run("6r2 H3K27ac ALL-SIL", function() {
  if (!nrow(h3)) stop("no H3K27ac")
  d <- melt(h3, id.vars = "file",
            measure.vars = intersect(c("ZEB1_prom_hg19", "ZEB1_gene_hg19", "ZEB1_prom_hg38", "ZEB1_gene_hg38"), names(h3)),
            variable.name = "window", value.name = "n")
  p <- ggplot(d, aes(window, n, fill = window)) +
    geom_col(width = 0.7, color = "white") +
    geom_text(aes(label = n), vjust = -0.3, size = 3.5) +
    scale_fill_manual(values = c("#E64B35", "#4DBBD5", "#00A087", "#3C5488"), guide = "none") +
    labs(title = "ALL-SIL H3K27ac peaks overlapping ZEB1 (GSE70734)",
         subtitle = "Single cell line; n = 1. Not a high vs low comparison.",
         x = NULL, y = "Overlapping peaks") + theme_sci() +
    theme(axis.text.x = element_text(angle = 25, hjust = 1))
  save_both(p, file.path(fig, "M6_r2_6.18_H3K27ac_ZEB1.pdf"), 7.0, 5.2)
})

run("6r2 scATAC ZEB1 peaks", function() {
  if (!nrow(sc)) stop("no scATAC")
  d <- copy(sc)
  d[, n_zeb := pmax(ZEB1_gene_hg38, ZEB1_gene_hg19, na.rm = TRUE)]
  d[, frac := n_zeb / pmax(n_peaks, 1)]
  p <- ggplot(d, aes(reorder(sample, n_zeb), n_zeb)) +
    geom_col(fill = "#3C5488", width = 0.7, color = "white") +
    coord_flip() +
    labs(title = "GSE280728 scATAC: peaks overlapping ZEB1",
         subtitle = "ArchR not installed; n = 7 peak-barcode libraries",
         x = NULL, y = "n ZEB1-overlapping peaks") + theme_sci()
  save_both(p, file.path(fig, "M6_r2_6.1_scATAC_ZEB1_peaks.pdf"), 7.2, 5.4)
})

run("6r2 mouse ATAC Zeb1", function() {
  if (!nrow(mm)) stop("no mouse ATAC")
  d <- copy(mm)
  d[, n_zeb := pmax(Zeb1_mm39, Zeb1_mm10, na.rm = TRUE)]
  d[, stage := factor(stage, levels = intersect(c("ETP", "ETP_trans", "DN2a", "DN2b", "DN3", "DN4", "unk"), unique(stage)))]
  p <- ggplot(d, aes(stage, n_zeb, fill = stage)) +
    geom_boxplot(width = 0.55, outlier.shape = NA, alpha = 0.85) +
    geom_jitter(width = 0.12, size = 2, color = "grey20") +
    scale_fill_manual(values = c(ETP = "#3C5488", ETP_trans = "#4DBBD5", DN2a = "#00A087",
                                 DN2b = "#91D1C2", DN3 = "#F39B7F", DN4 = "#E64B35", unk = "grey70"),
                      guide = "none") +
    labs(title = "GSE287757 in-vivo thymocyte ATAC peaks at Zeb1",
         subtitle = "narrowPeak overlap; biological n = replicate",
         x = NULL, y = "n peaks in Zeb1 window") + theme_sci()
  save_both(p, file.path(fig, "M6_r2_6.6_mouse_ATAC_Zeb1.pdf"), 7.4, 5.4)
})

message("M6 R2 FIGURES DONE")

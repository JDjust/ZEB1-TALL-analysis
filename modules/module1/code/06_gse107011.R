# Module 1 r2: GSE107011 29 immune cell types — T-lineage heritage vs T-ALL feature.
suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})
has_org <- requireNamespace("org.Hs.eg.db", quietly = TRUE)
has_ann <- requireNamespace("AnnotationDbi", quietly = TRUE)
has_cht <- requireNamespace("ComplexHeatmap", quietly = TRUE)
has_circlize <- requireNamespace("circlize", quietly = TRUE)
has_vp <- requireNamespace("variancePartition", quietly = TRUE)
has_limma <- requireNamespace("limma", quietly = TRUE)

root <- "/data-b/liangfuhua/projects/ZEB1_TALL_analysis/module1"
fig <- file.path(root, "figures"); tab <- file.path(root, "tables"); proc <- file.path(root, "processed")
dir.create(fig, FALSE, TRUE); dir.create(tab, FALSE, TRUE)
tpm_f <- "/data-b/liangfuhua/projects/TALL_dataset_download/data/GSE107011/suppl/GSE107011_Processed_data_TPM.txt.gz"
soft_f <- "/data-b/liangfuhua/projects/TALL_dataset_download/data/GSE107011/soft/GSE107011_family.soft.gz"

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

KEYS <- c("ZEB1", "ZEB2", "LMO2", "LYL1", "TCF7", "GATA3", "BCL11B", "MEF2C",
          "IL7R", "CD34", "CD3D", "CD3E", "CD4", "CD8A", "CD19", "MS4A1",
          "NKG7", "CD14", "SPI1", "PAX5", "EBF1", "TBX21", "RORC", "FOXP3",
          "LEF1", "NOTCH1", "MYC", "TAL1", "CD7", "CD79A")
FORCE_ENS <- c(
  ZEB1 = "ENSG00000148516", ZEB2 = "ENSG00000169554", LMO2 = "ENSG00000135363",
  LYL1 = "ENSG00000104903", TCF7 = "ENSG00000081059", GATA3 = "ENSG00000107485",
  BCL11B = "ENSG00000127152", MEF2C = "ENSG00000081189", IL7R = "ENSG00000168685",
  CD34 = "ENSG00000174059", CD3D = "ENSG00000167286", CD3E = "ENSG00000198851",
  CD4 = "ENSG00000010610", CD8A = "ENSG00000153563", CD19 = "ENSG00000177455",
  MS4A1 = "ENSG00000156738", NKG7 = "ENSG00000105374", CD14 = "ENSG00000170458",
  SPI1 = "ENSG00000066336", PAX5 = "ENSG00000196092", EBF1 = "ENSG00000164330",
  TBX21 = "ENSG00000073861", RORC = "ENSG00000143365", FOXP3 = "ENSG00000049768",
  LEF1 = "ENSG00000138795", NOTCH1 = "ENSG00000148400", MYC = "ENSG00000136997",
  TAL1 = "ENSG00000162367", CD7 = "ENSG00000173762", CD79A = "ENSG00000105369",
  CD2 = "ENSG00000116824", TRAC = "ENSG00000277734", LCK = "ENSG00000182866"
)
T_SIG <- c("CD3D", "CD3E", "CD7", "LCK", "TRAC", "BCL11B", "TCF7", "GATA3", "CD2")
B_SIG <- c("CD19", "MS4A1", "CD79A", "PAX5", "EBF1", "CD74")

message("read TPM")
tpm <- fread(cmd = paste("gzip -dc", shQuote(tpm_f)))
ens <- tpm[[1]]
ens0 <- sub("\\.\\d+$", "", ens)
mat <- as.matrix(tpm[, -1, with = FALSE])
storage.mode(mat) <- "double"
samples <- colnames(tpm)[-1]
message("matrix ", nrow(mat), " x ", ncol(mat))

ens2gene <- setNames(names(FORCE_ENS), unname(FORCE_ENS))
sym <- ens0
hit_force <- ens0 %in% names(ens2gene)
sym[hit_force] <- unname(ens2gene[ens0[hit_force]])
message("forced symbols ", sum(hit_force))
if (has_org && has_ann) {
  map <- AnnotationDbi::select(org.Hs.eg.db::org.Hs.eg.db, keys = unique(ens0),
                               columns = "SYMBOL", keytype = "ENSEMBL")
  map <- as.data.table(map)
  map <- map[!is.na(SYMBOL) & SYMBOL != ""]
  map <- map[!duplicated(ENSEMBL)]
  setkey(map, ENSEMBL)
  hit <- map[ens0, SYMBOL]
  fill <- !is.na(hit) & !hit_force
  sym[fill] <- hit[fill]
  message("org.Hs mapped extra ", sum(fill), " / ", length(ens0))
} else {
  message("org.Hs.eg.db missing; using forced key-gene map")
}

# collapse duplicate symbols by max mean
logm <- log2(mat + 1)
agg <- data.table(symbol = sym, idx = seq_along(sym), mu = rowMeans(logm))
agg <- agg[order(-mu)]
keep <- agg[!duplicated(symbol)]
logm <- logm[keep$idx, , drop = FALSE]
rownames(logm) <- keep$symbol

meta <- data.table(sample = samples)
meta[, donor := sub("_.*$", "", sample)]
meta[, cell_raw := sub("^[^_]+_", "", sample)]
lineage_of <- function(x) {
  s <- tolower(x)
  if (grepl("progenitor", s)) return("Progenitor")
  if (grepl("^b_|plasmablast|memory_b|naive_b|exhausted_b|switched|nsm", s) || grepl("^b-", s) || grepl("b_naive|b_nsm|b_ex|b_sm", s) || grepl("^bnaive|^bnsm|^bex|^bsm", s)) return("B")
  if (grepl("nk|natural_killer|natural.killer", s)) return("NK")
  if (grepl("mono|neutro|baso|dc|dendritic|myeloid", s)) return("Myeloid")
  if (grepl("cd4|cd8|th1|th2|th17|treg|tfh|mait|vd2|gd|t_naive|naive_cd|cm|em|te", s) ||
      grepl("^th|^treg|^tfh", s)) return("T")
  "Other"
}
# GSE107011 sample suffixes from Monaco 2019
meta[, lineage := vapply(cell_raw, lineage_of, "")]
# refine B types that use B_ prefix
meta[grepl("^B_", cell_raw) | grepl("^Bnaive|^B_NSM|^B_Ex|^B_SM|^Plasmablasts|^B_n", cell_raw, ignore.case = TRUE), lineage := "B"]
meta[grepl("Plasmablast", cell_raw, ignore.case = TRUE), lineage := "B"]
meta[grepl("Progenitor", cell_raw, ignore.case = TRUE), lineage := "Progenitor"]
meta[grepl("NK|Natural", cell_raw, ignore.case = TRUE), lineage := "NK"]
meta[grepl("mono|DC|neutro|baso|pDC|mDC", cell_raw, ignore.case = TRUE), lineage := "Myeloid"]
fwrite(meta, file.path(tab, "M1_r2_GSE107011_sample_meta.tsv"), sep = "\t")
message(paste(capture.output(print(meta[, .N, lineage][order(-N)])), collapse = "\n"))

score_mean <- function(genes) {
  g <- intersect(genes, rownames(logm))
  if (length(g) < 2) return(rep(NA_real_, ncol(logm)))
  z <- t(scale(t(logm[g, , drop = FALSE])))
  colMeans(z, na.rm = TRUE)
}
meta[, T_score := score_mean(T_SIG)]
meta[, B_score := score_mean(B_SIG)]
meta[, TB_axis := T_score - B_score]
for (g in KEYS) {
  if (g %in% rownames(logm)) meta[[g]] <- as.numeric(logm[g, meta$sample])
  else meta[[g]] <- NA_real_
}
fwrite(meta, file.path(tab, "M1_r2_GSE107011_scores.tsv"), sep = "\t")
message("ZEB1 finite ", sum(is.finite(meta$ZEB1)), " / ", nrow(meta),
        " lineage table:")
print(meta[, .(n = .N, zeb1 = mean(ZEB1, na.rm = TRUE)), lineage])
if (!sum(is.finite(meta$ZEB1))) stop("ZEB1 still all NA after mapping")

pal_lin <- c(T = "#E64B35", B = "#4DBBD5", NK = "#00A087", Myeloid = "#3C5488",
             Progenitor = "#F39B7F", Other = "#8491B4")

p <- ggplot(meta, aes(reorder(cell_raw, ZEB1, FUN = median), ZEB1, fill = lineage)) +
  geom_boxplot(width = 0.7, outlier.size = 0.6) +
  scale_fill_manual(values = pal_lin) +
  coord_flip() +
  labs(title = "ZEB1 across 29 healthy immune cell types",
       subtitle = "GSE107011 TPM, log2(TPM+1); 4 donors",
       x = NULL, y = "ZEB1") + theme_sci()
save_both(p, file.path(fig, "M1_r2_1.1_ZEB1_celltypes.pdf"), 8.2, 9.5)

p <- ggplot(meta[lineage %in% c("T", "B", "NK", "Myeloid", "Progenitor")],
            aes(lineage, ZEB1, fill = lineage)) +
  geom_boxplot(width = 0.55, outlier.shape = NA, alpha = 0.9) +
  geom_point(position = position_jitter(width = 0.12, seed = 1), size = 1.5, alpha = 0.7) +
  scale_fill_manual(values = pal_lin, guide = "none") +
  labs(title = "ZEB1 is a T-lineage gene in healthy blood",
       x = NULL, y = "ZEB1 log2(TPM+1)") + theme_sci()
save_both(p, file.path(fig, "M1_r2_1.2_ZEB1_lineage.pdf"), 6.4, 5.4)

long <- melt(meta, id.vars = c("sample", "donor", "cell_raw", "lineage"),
             measure.vars = intersect(c("ZEB1", "ZEB2", "LMO2", "TCF7", "GATA3", "BCL11B"), names(meta)),
             variable.name = "gene", value.name = "expr")
p <- ggplot(long[lineage %in% c("T", "B", "NK", "Myeloid")],
            aes(lineage, expr, fill = lineage)) +
  geom_boxplot(width = 0.65, outlier.size = 0.4) +
  facet_wrap(~gene, scales = "free_y", nrow = 2) +
  scale_fill_manual(values = pal_lin, guide = "none") +
  labs(title = "Lineage TFs in GSE107011", x = NULL, y = "log2(TPM+1)") + theme_sci()
save_both(p, file.path(fig, "M1_r2_1.3_TF_lineage_boxes.pdf"), 9.2, 6.2)

p <- ggplot(meta, aes(TB_axis, ZEB1, color = lineage)) +
  geom_point(size = 2.2, alpha = 0.85) +
  geom_smooth(method = "lm", se = TRUE, color = "grey30", linewidth = 0.6, inherit.aes = FALSE,
              aes(TB_axis, ZEB1)) +
  scale_color_manual(values = pal_lin) +
  labs(title = "ZEB1 tracks the T-minus-B lineage axis",
       x = "T-lineage score − B-lineage score", y = "ZEB1") + theme_sci()
save_both(p, file.path(fig, "M1_r2_1.4_ZEB1_vs_TBaxis.pdf"), 6.8, 5.6)

p <- ggplot(meta, aes(ZEB1, ZEB2, color = lineage)) +
  geom_hline(yintercept = median(meta$ZEB2, na.rm = TRUE), linetype = 2, color = "grey70") +
  geom_vline(xintercept = median(meta$ZEB1, na.rm = TRUE), linetype = 2, color = "grey70") +
  geom_point(size = 2.4, alpha = 0.9) +
  scale_color_manual(values = pal_lin) +
  labs(title = "ZEB1/ZEB2 phase in healthy immune cells",
       subtitle = "T cells occupy ZEB1-high / ZEB2-variable space") + theme_sci()
save_both(p, file.path(fig, "M1_r2_1.5_ZEB1_ZEB2_phase.pdf"), 6.8, 5.6)

# heatmap of key genes x cell type means
if (has_cht && has_circlize) {
  guse <- intersect(KEYS, rownames(logm))
  ct <- unique(meta$cell_raw)
  hm <- sapply(ct, function(cl) rowMeans(logm[guse, meta$sample[meta$cell_raw == cl], drop = FALSE]))
  hm <- t(scale(t(hm)))
  lin <- meta[, .(lineage = lineage[1]), cell_raw]
  setkey(lin, cell_raw)
  col_lin <- pal_lin[lin[colnames(hm), lineage]]
  ha <- ComplexHeatmap::HeatmapAnnotation(lineage = lin[colnames(hm), lineage],
                                          col = list(lineage = pal_lin))
  pdf(file.path(fig, "M1_r2_1.6_celltype_heatmap.pdf"), width = 12, height = 6.2)
  print(ComplexHeatmap::Heatmap(hm, name = "z", top_annotation = ha,
                                cluster_rows = TRUE, cluster_columns = TRUE,
                                column_names_gp = grid::gpar(fontsize = 7),
                                row_names_gp = grid::gpar(fontsize = 9),
                                column_title = "GSE107011 cell-type means"))
  dev.off()
}

if (has_vp && "ZEB1" %in% names(meta)) {
  form <- ~ (1 | lineage) + (1 | donor) + (1 | cell_raw)
  vp <- tryCatch(variancePartition::fitExtractVarPartModel(t(as.matrix(meta[, .(ZEB1)])), form, meta),
                 error = function(e) { message("vp fail ", e$message); NULL })
  if (!is.null(vp)) {
    fwrite(as.data.table(vp, keep.rownames = "gene"), file.path(tab, "M1_r2_1.7_variancePartition.tsv"), sep = "\t")
    pdf(file.path(fig, "M1_r2_1.7_variancePartition.pdf"), width = 6.2, height = 4.2)
    print(variancePartition::plotPercentBars(vp))
    dev.off()
  }
}

# limma lineage (T vs each)
if (has_limma && "ZEB1" %in% rownames(logm)) {
  d <- meta[lineage %in% c("T", "B", "NK", "Myeloid")]
  y <- logm["ZEB1", d$sample]
  fit <- lm(y ~ lineage, data = d)
  sm <- summary(fit)
  coefs <- as.data.table(sm$coefficients, keep.rownames = "term")
  fwrite(coefs, file.path(tab, "M1_r2_1.8_ZEB1_lm_lineage.tsv"), sep = "\t")
}

# Spearman ZEB1 vs T_score
ok <- is.finite(meta$ZEB1) & is.finite(meta$T_score)
ct <- suppressWarnings(cor.test(meta$ZEB1[ok], meta$T_score[ok], method = "spearman"))
fwrite(data.table(n = sum(ok), rho = unname(ct$estimate), p = ct$p.value),
       file.path(tab, "M1_r2_1.4_ZEB1_Tscore_spearman.tsv"), sep = "\t")

# cell-type median ranking
rk <- meta[, .(n = .N, median_ZEB1 = median(ZEB1, na.rm = TRUE),
               mean_ZEB1 = mean(ZEB1, na.rm = TRUE), lineage = lineage[1]), cell_raw]
rk <- rk[order(-median_ZEB1)]
fwrite(rk, file.path(tab, "M1_r2_1.1_celltype_ZEB1_rank.tsv"), sep = "\t")

p <- ggplot(rk, aes(reorder(cell_raw, median_ZEB1), median_ZEB1, fill = lineage)) +
  geom_col(width = 0.8) + coord_flip() +
  scale_fill_manual(values = pal_lin) +
  labs(title = "Median ZEB1 by cell type", x = NULL, y = "median ZEB1") + theme_sci()
save_both(p, file.path(fig, "M1_r2_1.9_median_rank.pdf"), 8.0, 9.2)

fwrite(data.table(
  note = "GSE107011 is healthy immune RNA-seq (Monaco et al.). T-ALL is not on this platform. Compare lineage effect here with Module 1 T-ALL vs B-ALL Hedges g; do not merge TPM with TARGET."
), file.path(tab, "M1_r2_limitation.txt"), sep = "\t")

message("M1 r2 DONE fig=", length(list.files(fig, pattern = "^M1_r2_")))

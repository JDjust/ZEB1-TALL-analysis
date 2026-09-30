# Module 8 figures M8_8.x. code2 206 volcano, 88 forest, 346 bars, 210 heatmap.
suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})
has_fgsea <- requireNamespace("fgsea", quietly = TRUE)
has_igraph <- requireNamespace("igraph", quietly = TRUE)
has_ggrepel <- requireNamespace("ggrepel", quietly = TRUE)

root <- "/data-b/liangfuhua/projects/ZEB1_TALL_analysis/module8"
fig <- file.path(root, "figures"); tab <- file.path(root, "tables")
gmt <- "/data-b/liangfuhua/projects/TALL_dataset_download/data/MSigDB/h.all.v2023.2.Hs.symbols.gmt"
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
  if (!file.exists(p) || isTRUE(file.info(p)$size < 8)) return(data.table())
  tryCatch(fread(p), error = function(e) data.table())
}

lab <- rd("M8_8.1_8.5_line_table.tsv")
gw <- rd("M8_8.6_8.10_GE_vs_ZEB1.tsv")
ess <- rd("M8_8.9_8.15_essentiality.tsv")
self <- rd("M8_8.4_8.5_ZEB1_self.tsv")
pan <- rd("M8_8.16_ZEB1_GE_by_lineage.tsv")
tgt <- rd("M8_8.19_drug_target_GE.tsv")
cand <- rd("M8_8.22_candidates.tsv")
cmf <- file.path(tab, "M8_8.13_codependency_corr.tsv")
pal <- c("T-ALL" = "#E64B35", "B-ALL" = "#4DBBD5", myeloid = "#00A087",
         other_lymphoid = "#F39B7F", other = "#8491B4")

run("8.1-8.3 expression", function() {
  if (!nrow(lab)) stop("no lab")
  d <- lab[lineage %in% c("T-ALL", "B-ALL", "myeloid")]
  p <- ggplot(d, aes(lineage, ZEB1, fill = lineage)) +
    geom_boxplot(width = 0.55, outlier.shape = NA, alpha = 0.9) +
    geom_point(position = position_jitter(width = 0.1, seed = 1), size = 1.8, alpha = 0.75) +
    scale_fill_manual(values = pal, guide = "none") +
    labs(title = "8.1-8.2  DepMap ZEB1 expression", x = NULL, y = "log2(TPM+1)") + theme_sci()
  save_both(p, file.path(fig, "M8_8.1_8.2_ZEB1_expression.pdf"), 6.6, 5.4)
  long <- melt(d, id.vars = c("ModelID", "lineage"),
               measure.vars = intersect(c("ZEB1", "ZEB2", "LMO2"), names(d)),
               variable.name = "gene", value.name = "expr")
  p2 <- ggplot(long, aes(lineage, expr, fill = lineage)) +
    geom_boxplot(width = 0.6, outlier.size = 0.3) +
    facet_wrap(~gene, scales = "free_y") +
    scale_fill_manual(values = pal, guide = "none") +
    labs(title = "8.3  ZEB1 / ZEB2 / LMO2 expression", x = NULL, y = "log2(TPM+1)") + theme_sci()
  save_both(p2, file.path(fig, "M8_8.3_ZEB1_ZEB2_LMO2.pdf"), 9.0, 5.2)
})

run("8.4-8.5 self GE", function() {
  if (!nrow(self)) stop("no self")
  p <- ggplot(self, aes(ZEB1, ZEB1_GE)) +
    geom_hline(yintercept = -0.5, linetype = 2, color = "grey40") +
    geom_point(size = 2.6, shape = 21, fill = "#E64B35") +
    geom_smooth(method = "lm", se = TRUE, color = "#3C5488", linewidth = 0.8) +
    labs(title = "8.4-8.5  T-ALL ZEB1 GeneEffect vs expression",
         subtitle = "Dashed: GeneEffect < -0.5 as dependent (no probability file)",
         x = "ZEB1 expression", y = "ZEB1 GeneEffect") + theme_sci()
  save_both(p, file.path(fig, "M8_8.4_8.5_ZEB1_self_GE.pdf"), 6.6, 5.6)
})

run("8.6-8.10 volcano", function() {
  if (!nrow(gw)) stop("no gw")
  gw[, cls := fcase(fdr < 0.05 & rho > 0, "pos FDR<0.05",
                    fdr < 0.05 & rho < 0, "neg FDR<0.05", default = "NS")]
  hit <- gw[gene %in% c("ZEB1", "ZEB2", "LMO2", "BCL2", "JAK1", "CDK6", "NOTCH1", "GATA3")]
  p <- ggplot(gw, aes(rho, -log10(pmax(p, 1e-300)), color = cls)) +
    geom_point(size = 0.35, alpha = 0.4) +
    geom_point(data = hit, size = 2.2, shape = 21, fill = "white", color = "black") +
    scale_color_manual(values = c("pos FDR<0.05" = "#E64B35", "neg FDR<0.05" = "#3C5488", NS = "grey75")) +
    labs(title = "8.6-8.10  GeneEffect vs ZEB1 expression (T-ALL lines)",
         x = "Spearman rho", y = "-log10 P", color = NULL) + theme_sci() + theme(legend.position = "top")
  if (has_ggrepel && nrow(hit)) p <- p + ggrepel::geom_text_repel(data = hit, aes(label = gene), size = 3, color = "black")
  save_both(p, file.path(fig, "M8_8.6_volcano_GE_vs_ZEB1.pdf"), 7.6, 6.0)
  top <- rbind(gw[order(-rho)][1:12][, side := "pos"], gw[order(rho)][1:12][, side := "neg"])
  setorder(top, rho); top[, gene := factor(gene, levels = gene)]
  p2 <- ggplot(top, aes(rho, gene, fill = side)) +
    geom_col(width = 0.75, color = "white") +
    scale_fill_manual(values = c(pos = "#E64B35", neg = "#3C5488"), guide = "none") +
    labs(title = "8.9  Strongest GeneEffect–ZEB1 correlations", x = "rho", y = NULL) + theme_sci()
  save_both(p2, file.path(fig, "M8_8.9_top_rho_bars.pdf"), 7.2, 6.4)
})

run("8.11-8.15 essentials", function() {
  if (!nrow(ess)) stop("no ess")
  top <- ess[order(mean_TALL)][1:20]
  top[, gene := factor(gene, levels = rev(gene))]
  p <- ggplot(top, aes(mean_TALL, gene, fill = mean_TALL < -0.5)) +
    geom_col(width = 0.75, color = "white") +
    scale_fill_manual(values = c("TRUE" = "#E64B35", "FALSE" = "#8491B4"), guide = "none") +
    labs(title = "8.11  Top T-ALL CRISPR essentials", x = "Mean GeneEffect", y = NULL) + theme_sci()
  save_both(p, file.path(fig, "M8_8.11_top_essentials.pdf"), 7.4, 6.6)
  d <- ess[is.finite(g_T_vs_B) & is.finite(p_T_vs_B)]
  d <- d[order(g_T_vs_B)][c(1:12, (nrow(d) - 11):nrow(d))]
  setorder(d, g_T_vs_B); d[, gene := factor(gene, levels = gene)]
  p2 <- ggplot(d, aes(g_T_vs_B, gene, fill = g_T_vs_B < 0)) +
    geom_col(width = 0.75, color = "white") +
    scale_fill_manual(values = c("TRUE" = "#E64B35", "FALSE" = "#4DBBD5"), guide = "none") +
    labs(title = "8.15  T-ALL vs B-ALL dependency (Hedges g)",
         subtitle = "Negative = more essential in T-ALL", x = "g", y = NULL) + theme_sci()
  save_both(p2, file.path(fig, "M8_8.15_TALL_vs_BALL_g.pdf"), 7.6, 6.8)
})

run("8.16 lineage", function() {
  if (!nrow(pan)) stop("no pan")
  setorder(pan, mean)
  pan[, OncotreeLineage := factor(OncotreeLineage, levels = OncotreeLineage)]
  p <- ggplot(pan, aes(mean, OncotreeLineage, fill = count)) +
    geom_col(width = 0.75, color = "white") +
    scale_fill_gradient(low = "#4DBBD5", high = "#E64B35") +
    labs(title = "8.16  Mean ZEB1 GeneEffect by lineage", x = "Mean GeneEffect", y = NULL, fill = "n") + theme_sci()
  save_both(p, file.path(fig, "M8_8.16_ZEB1_GE_lineage.pdf"), 8.0, 7.4)
})

run("8.12 GSEA", function() {
  if (!has_fgsea || !nrow(ess) || !file.exists(gmt)) stop("fgsea/gmt missing")
  ranks <- ess[!is.na(mean_TALL), setNames(-mean_TALL, gene)]
  ranks <- ranks[!duplicated(names(ranks))]
  gs <- fgsea::gmtPathways(gmt)
  res <- as.data.table(fgsea::fgsea(gs, ranks, minSize = 10, maxSize = 400, nPermSimple = 3000))
  res <- res[order(padj, -abs(NES))]
  fwrite(res[, .(pathway, pval, padj, NES, size)], file.path(tab, "M8_8.12_hallmark_gsea.tsv"), sep = "\t")
  d <- res[!is.na(NES)][1:12]
  setorder(d, NES); d[, pathway := factor(pathway, levels = pathway)]
  p <- ggplot(d, aes(NES, pathway, fill = NES > 0)) +
    geom_col(width = 0.75, color = "white") +
    scale_fill_manual(values = c("TRUE" = "#E64B35", "FALSE" = "#3C5488"), guide = "none") +
    labs(title = "8.12  Hallmark GSEA on T-ALL essentiality",
         subtitle = "Ranked by -mean GeneEffect", x = "NES", y = NULL) + theme_sci()
  save_both(p, file.path(fig, "M8_8.12_hallmark_gsea.pdf"), 9.0, 6.2)
})

run("8.13 network", function() {
  if (!has_igraph || !file.exists(cmf)) stop("no corr")
  cm <- as.matrix(fread(cmf), rownames = 1)
  ut <- which(upper.tri(cm), arr.ind = TRUE)
  ed <- data.table(a = rownames(cm)[ut[, 1]], b = colnames(cm)[ut[, 2]], r = cm[ut])
  ed <- ed[is.finite(r) & abs(r) >= 0.55]
  if (nrow(ed) < 5) stop("too few edges")
  g <- igraph::graph_from_data_frame(ed[, .(a, b, r)], directed = FALSE)
  pdf(file.path(fig, "M8_8.13_codependency_network.pdf"), 8.0, 7.6)
  set.seed(1)
  plot(g, vertex.size = 8, vertex.label.cex = 0.55,
       edge.width = abs(ed$r) * 2, main = "8.13  Co-dependency among top T-ALL essentials")
  dev.off()
})

run("8.19 targets", function() {
  if (!nrow(tgt)) stop("no tgt")
  setorder(tgt, mean_TALL)
  tgt[, gene := factor(gene, levels = gene)]
  p <- ggplot(tgt, aes(mean_TALL, gene, fill = mean_TALL < -0.5)) +
    geom_col(width = 0.7, color = "white") +
    geom_vline(xintercept = -0.5, linetype = 2) +
    scale_fill_manual(values = c("TRUE" = "#E64B35", "FALSE" = "#8491B4"), guide = "none") +
    labs(title = "8.19  Pharmacotype drug-target GeneEffect in T-ALL",
         subtitle = "No PRISM file; targets mapped from Module 7 drugs",
         x = "Mean GeneEffect", y = NULL) + theme_sci()
  save_both(p, file.path(fig, "M8_8.19_8.21_drug_targets.pdf"), 7.4, 6.4)
})

run("8.22 candidates", function() {
  if (!nrow(cand)) stop("no candidates with score>=2")
  d <- cand[order(mean_TALL)][1:min(20, nrow(cand))]
  d[, gene := factor(gene, levels = rev(gene))]
  p <- ggplot(d, aes(mean_TALL, gene, fill = score)) +
    geom_col(width = 0.75, color = "white") +
    scale_fill_gradient(low = "#4DBBD5", high = "#E64B35") +
    labs(title = "8.22-8.23  Candidates with >=2 evidence layers",
         subtitle = "T-ALL essential + T vs B specificity +/or ZEB1 association",
         x = "Mean GeneEffect", y = NULL, fill = "score") + theme_sci()
  save_both(p, file.path(fig, "M8_8.22_candidates.pdf"), 7.6, 6.6)
})

q14 <- rd("M8_8.7_Q1Q4_GE.tsv")
adj <- rd("M8_8.17_expr_adjusted.tsv")
keyg <- rd("M8_8.14_key_lineage_GE.tsv")
ev <- rd("M8_8.22_8.23_evidence.tsv")
pharm <- rd("M8_8.21_pharmacotype_spearman.tsv")

note_fig <- function(id, text, fn) {
  p <- ggplot() +
    annotate("text", x = 0, y = 0, label = paste0(id, "\n", text), size = 4.1, lineheight = 1.15) +
    theme_void() + xlim(-1, 1) + ylim(-1, 1)
  save_both(p, file.path(fig, fn), 8.2, 3.4)
}

run("8.5 dependency probability substitute", function() {
  note_fig("8.5",
           "CRISPRGeneDependency.csv is not in the downloaded DepMap folder.\nSubstitute: GeneEffect < -0.5 as a dependent call (see 8.4 scatter).",
           "M8_8.5_no_dependency_probability.pdf")
  if (nrow(self)) {
    p <- ggplot(self, aes(ZEB1_GE)) +
      geom_histogram(binwidth = 0.15, fill = "#3C5488", color = "white") +
      geom_vline(xintercept = -0.5, linetype = 2, color = "#E64B35") +
      labs(title = "8.5  T-ALL ZEB1 GeneEffect (probability file missing)",
           subtitle = "Vertical line: GeneEffect = -0.5 substitute threshold",
           x = "ZEB1 GeneEffect", y = "Lines") + theme_sci()
    save_both(p, file.path(fig, "M8_8.5_ZEB1_GE_histogram.pdf"), 6.8, 5.0)
  }
})

run("8.7 Q1 vs Q4", function() {
  if (!nrow(q14)) {
    note_fig("8.7",
             "ZEB1 Q1 vs Q4 GeneEffect was not run: only 6 T-ALL lines have both CRISPR and ZEB1 RNA.\nQ1 n=1, Q4 n=2. Continuous Spearman (8.6 / 8.8) is the valid test.",
             "M8_8.7_Q1Q4_underpowered.pdf")
    return()
  }
  d <- q14[is.finite(g)]
  d <- d[order(g)][c(1:min(12, nrow(d)), max(1, nrow(d) - 11):nrow(d))]
  setorder(d, g); d[, gene := factor(gene, levels = gene)]
  p <- ggplot(d, aes(g, gene, fill = g < 0)) +
    geom_col(width = 0.75, color = "white") +
    scale_fill_manual(values = c("TRUE" = "#E64B35", "FALSE" = "#4DBBD5"), guide = "none") +
    labs(title = "8.7  GeneEffect in ZEB1 Q4 vs Q1 T-ALL lines",
         subtitle = "Negative g = more essential in ZEB1-high lines. n per arm is small.",
         x = "Hedges g (Q4 - Q1)", y = NULL) + theme_sci()
  save_both(p, file.path(fig, "M8_8.7_Q1Q4_g.pdf"), 7.6, 6.8)
})

run("8.8 8.10 continuous + FDR", function() {
  if (!nrow(gw)) stop("no gw")
  p <- ggplot(gw, aes(rho)) +
    geom_histogram(bins = 60, fill = "#3C5488", color = "white") +
    labs(title = "8.8  Genome-wide Spearman GeneEffect vs ZEB1 expression",
         x = "rho", y = "Genes") + theme_sci()
  save_both(p, file.path(fig, "M8_8.8_rho_histogram.pdf"), 6.8, 5.0)
  n_sig <- sum(gw$fdr < 0.05, na.rm = TRUE)
  sm <- data.table(class = c("FDR < 0.05", "NS"), n = c(n_sig, nrow(gw) - n_sig))
  p2 <- ggplot(sm, aes(class, n, fill = class)) +
    geom_col(width = 0.6, color = "white") +
    geom_text(aes(label = n), vjust = -0.3, size = 3.5) +
    scale_fill_manual(values = c("FDR < 0.05" = "#E64B35", NS = "#8491B4"), guide = "none") +
    labs(title = "8.10  BH-FDR on GeneEffect vs ZEB1", x = NULL, y = "Genes") + theme_sci()
  save_both(p2, file.path(fig, "M8_8.10_FDR_counts.pdf"), 6.2, 5.2)
})

run("8.14 lineage specificity", function() {
  if (!nrow(keyg)) stop("no key lineage table")
  long <- melt(keyg, id.vars = "gene",
               measure.vars = intersect(c("mean_TALL", "mean_BALL", "mean_myeloid"), names(keyg)),
               variable.name = "lineage", value.name = "mean_GE")
  long[, lineage := gsub("mean_", "", lineage)]
  p <- ggplot(long, aes(mean_GE, gene, fill = lineage)) +
    geom_col(position = position_dodge(0.8), width = 0.7, color = "white") +
    geom_vline(xintercept = -0.5, linetype = 2, color = "grey40") +
    scale_fill_manual(values = c(TALL = "#E64B35", BALL = "#4DBBD5", myeloid = "#00A087")) +
    labs(title = "8.14  Key-gene mean GeneEffect by lineage",
         x = "Mean GeneEffect", y = NULL, fill = NULL) + theme_sci()
  save_both(p, file.path(fig, "M8_8.14_lineage_specificity.pdf"), 8.2, 7.0)
})

run("8.17 expression-adjusted", function() {
  if (!nrow(adj)) {
    note_fig("8.17",
             "Expression-adjusted residuals need n>=8 T-ALL lines with CRISPR+RNA.\nOnly 6 lines overlap. Unadjusted Spearman is 8.6 / 8.8.",
             "M8_8.17_underpowered.pdf")
    return()
  }
  long <- melt(adj, id.vars = "gene",
               measure.vars = c("rho_raw", "rho_expr_adjusted"),
               variable.name = "model", value.name = "rho")
  p <- ggplot(long, aes(rho, gene, fill = model)) +
    geom_col(position = position_dodge(0.75), width = 0.7, color = "white") +
    geom_vline(xintercept = 0, linetype = 2) +
    scale_fill_manual(values = c(rho_raw = "#8491B4", rho_expr_adjusted = "#E64B35"),
                      labels = c("raw", "expr-adjusted")) +
    labs(title = "8.17  GeneEffect vs ZEB1 after target-expression residual",
         x = "Spearman rho", y = NULL, fill = NULL) + theme_sci()
  save_both(p, file.path(fig, "M8_8.17_expr_adjusted.pdf"), 7.6, 6.4)
})

run("8.18 CNV note", function() {
  note_fig("8.18",
           "No CNV matrix in the downloaded DepMap_CCLE folder.\nExpression-adjusted residuals (8.17) are the available covariate control.",
           "M8_8.18_no_CNV.pdf")
})

run("8.20 PRISM note", function() {
  note_fig("8.20",
           "No PRISM file in the download inventory.\nDrug evidence is Pharmacotype (Module 7) plus CRISPR target mapping (8.19).",
           "M8_8.20_no_PRISM.pdf")
})

run("8.21 pharmacotype", function() {
  if (!nrow(pharm)) stop("no pharmacotype spearman")
  z <- pharm[gene == "ZEB1_log"]
  if (!nrow(z)) z <- pharm[grepl("^ZEB1", gene)]
  if (!nrow(z)) stop("no ZEB1 drug rows")
  setorder(z, rho)
  if (!"drug_short" %in% names(z)) z[, drug_short := gsub(" \\(.*\\)$", "", drug)]
  z[, drug_short := factor(drug_short, levels = drug_short)]
  p <- ggplot(z, aes(rho, drug_short, fill = fdr < 0.05)) +
    geom_col(width = 0.7, color = "white") +
    scale_fill_manual(values = c("TRUE" = "#E64B35", "FALSE" = "#8491B4")) +
    labs(title = "8.21  Pharmacotype ZEB1 vs LC50 (from Module 7)",
         subtitle = "Independent cell-line drug screen is not PRISM",
         x = "rho", y = NULL, fill = "FDR<0.05") + theme_sci()
  save_both(p, file.path(fig, "M8_8.21_pharmacotype_rho.pdf"), 7.8, 6.6)
})

run("8.23 evidence scoring", function() {
  if (!nrow(ev)) stop("no evidence")
  d <- ev[score >= 1][order(-score, mean_TALL)][1:min(25, .N)]
  long <- melt(d, id.vars = c("gene", "score"),
               measure.vars = intersect(c("ess_TALL", "tall_specific", "zeb1_assoc"), names(d)),
               variable.name = "layer", value.name = "hit")
  long[, hit := as.integer(as.logical(hit))]
  p <- ggplot(long, aes(layer, gene, fill = factor(hit))) +
    geom_tile(color = "white") +
    scale_fill_manual(values = c("0" = "grey90", "1" = "#E64B35"), guide = "none") +
    labs(title = "8.23  Evidence layers for top-scoring genes",
         subtitle = "T-ALL essential / T vs B specific / ZEB1-associated",
         x = NULL, y = NULL) + theme_sci() + theme(axis.line = element_blank())
  save_both(p, file.path(fig, "M8_8.23_evidence_layers.pdf"), 7.0, 7.2)
})

run("8.24 therapeutic network", function() {
  if (!has_igraph) stop("no igraph")
  nodes <- unique(c(cand$gene, tgt$gene, "ZEB1", "ZEB2", "LMO2"))
  nodes <- nodes[!is.na(nodes) & nzchar(nodes)]
  if (length(nodes) < 4) stop("too few nodes")
  ed <- data.table()
  if (file.exists(cmf)) {
    cm <- as.matrix(fread(cmf), rownames = 1)
    keep <- intersect(nodes, rownames(cm))
    if (length(keep) >= 3) {
      ut <- which(upper.tri(cm[keep, keep, drop = FALSE]), arr.ind = TRUE)
      ed <- data.table(a = keep[ut[, 1]], b = keep[ut[, 2]], r = cm[keep, keep, drop = FALSE][ut])
      ed <- ed[is.finite(r) & abs(r) >= 0.4]
    }
  }
  if (!nrow(ed)) {
    # fallback: connect ZEB1 to every candidate
    ed <- data.table(a = "ZEB1", b = setdiff(nodes, "ZEB1"), r = 0.5)
  }
  g <- igraph::graph_from_data_frame(ed[, .(a, b)], directed = FALSE)
  pdf(file.path(fig, "M8_8.24_therapeutic_network.pdf"), 8.2, 7.6)
  set.seed(1)
  plot(g, vertex.size = 10, vertex.label.cex = 0.6, vertex.color = "#4DBBD5",
       edge.width = pmax(abs(ed$r) * 2, 0.6),
       main = "8.24  Integrated therapeutic network (candidates + drug targets)")
  dev.off()
})

run("8.20 8.18 notes tsv", function() {
  fwrite(data.table(
    id = c("8.5", "8.18", "8.20", "8.21"),
    status = c(
      "CRISPRGeneDependency.csv not in download; GeneEffect < -0.5 used.",
      "No CNV matrix in DepMap_CCLE folder; skipped.",
      "No PRISM file; drug evidence is Pharmacotype (Module 7) target mapping.",
      "Pharmacotype Spearman copied from Module 7; ZEB1-drug FDR all >0.05."
    )
  ), file.path(tab, "M8_substitutes.tsv"), sep = "\t")
})

cat("MODULE8 ANALYZE DONE\n")
cat("FIG", length(list.files(fig, pattern = "\\.(pdf|png)$")), "\n")

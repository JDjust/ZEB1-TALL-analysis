# Module 10 GSE287751 scRNA figures. n is HTO group, not cell.
suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})
root <- "/data-b/liangfuhua/projects/ZEB1_TALL_analysis/module10"
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

means <- rd("M10_r2_HTO_group_means.tsv")
zl <- rd("M10_r2_HTO_Zeb1_Lmo2.tsv")
cells <- rd("M10_r2_scrna_plot_cells.tsv.gz")

run("10r2 group Lmo2 Zeb1", function() {
  d <- if (nrow(zl)) copy(zl) else copy(means)
  if (!nrow(d)) stop("no group means")
  if ("Zeb1_mean" %in% names(d)) {
    long <- melt(d, id.vars = intersect(c("group", "n"), names(d)),
                 measure.vars = intersect(c("Zeb1_mean", "Lmo2_mean"), names(d)),
                 variable.name = "gene", value.name = "mean")
    long[, gene := gsub("_mean", "", gene)]
  } else {
    gcol <- setdiff(names(d), c("n_cells", "geno"))[1]
    long <- melt(d, id.vars = intersect(c(gcol, "n_cells", "geno"), names(d)),
                 measure.vars = intersect(c("g_Zeb1", "g_Lmo2"), names(d)),
                 variable.name = "gene", value.name = "mean")
    setnames(long, gcol, "group")
    long[, gene := gsub("^g_", "", gene)]
  }
  p <- ggplot(long, aes(group, mean, fill = gene)) +
    geom_col(position = position_dodge(0.8), width = 0.7, color = "white") +
    scale_fill_manual(values = c(Zeb1 = "#E64B35", Lmo2 = "#4DBBD5",
                                 Zeb1_mean = "#E64B35", Lmo2_mean = "#4DBBD5")) +
    labs(title = "GSE287751 hashed scRNA: Lmo2 and Zeb1 by HTO group",
         subtitle = "Biological n is the HTO group, not the cell",
         x = NULL, y = "Mean log-normalized expression", fill = NULL) +
    theme_sci() + theme(axis.text.x = element_text(angle = 35, hjust = 1))
  save_both(p, file.path(fig, "M10_r2_HTO_Lmo2_Zeb1.pdf"), 8.2, 5.4)
})

run("10r2 UMAP", function() {
  if (!nrow(cells) || !"UMAP1" %in% names(cells)) stop("no cells")
  gcol <- setdiff(names(cells), c("UMAP1", "UMAP2", grep("^g_", names(cells), value = TRUE)))[1]
  p <- ggplot(cells, aes(UMAP1, UMAP2, color = .data[[gcol]])) +
    geom_point(size = 0.25, alpha = 0.8) +
    labs(title = "GSE287751 Lmo2 KO hashed scRNA", color = gcol) +
    theme_sci() + theme(axis.ticks = element_blank(), axis.text = element_blank())
  save_both(p, file.path(fig, "M10_r2_scrna_umap_group.pdf"), 7.2, 5.6)
})

run("10r2 UMAP Zeb1", function() {
  if (!nrow(cells) || !"g_Zeb1" %in% names(cells)) stop("no Zeb1")
  p <- ggplot(cells, aes(UMAP1, UMAP2, color = g_Zeb1)) +
    geom_point(size = 0.25, alpha = 0.85) +
    scale_color_gradientn(colors = c("#3C5488", "#EEEEEE", "#E64B35")) +
    labs(title = "Zeb1 in GSE287751 hashed scRNA", color = "Zeb1") +
    theme_sci() + theme(axis.ticks = element_blank(), axis.text = element_blank())
  save_both(p, file.path(fig, "M10_r2_scrna_umap_Zeb1.pdf"), 7.0, 5.6)
})

message("M10 R2 FIGURES DONE")

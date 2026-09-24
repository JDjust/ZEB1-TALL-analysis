pkgs <- c("ggplot2","ggpubr","pheatmap","pROC","limma","metafor","reshape2","plyr",
          "GEOquery","hgu133plus2.db","AnnotationDbi","cowplot","ggsci","data.table")
for (p in pkgs) {
  cat(p, as.character(requireNamespace(p, quietly=TRUE)), "\n")
}

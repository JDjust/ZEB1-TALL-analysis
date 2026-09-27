pkgs <- c("rlang","caret","glmnet","LiblineaR","kknn","randomForest",
          "ranger","caTools","elasticnet","singscore","ggplot2","ggrepel","umap")
for (p in pkgs) {
  cat(p, requireNamespace(p, quietly = TRUE), "\n")
}
cat("R", R.version.string, "\n")

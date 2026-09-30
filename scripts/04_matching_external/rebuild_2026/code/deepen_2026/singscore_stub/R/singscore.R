# Official singscore rank-mean (Foroutan et al. 2018).
# ALLCatchR2 uses only column 1 of simpleScore (TotalScore).

rankGenes <- function(expreMatrix, tiesMethod = "min", stableGenes = NULL) {
  m <- as.matrix(expreMatrix)
  if (is.null(rownames(m))) stop("rankGenes: expression matrix needs rownames")
  rankData <- apply(m, 2, rank, ties.method = tiesMethod)
  rownames(rankData) <- rownames(m)
  colnames(rankData) <- colnames(m)
  rankData
}

.score_one <- function(normRanks, genes, centerScore) {
  genes <- intersect(as.character(genes), rownames(normRanks))
  if (!length(genes)) {
    warning("None of the genes in the geneset were found")
    return(rep(if (isTRUE(centerScore)) -0.5 else 0, ncol(normRanks)))
  }
  s <- colMeans(normRanks[genes, , drop = FALSE], na.rm = TRUE)
  if (isTRUE(centerScore)) s <- s - 0.5
  s
}

.disp_one <- function(normRanks, genes, dispersionFun) {
  genes <- intersect(as.character(genes), rownames(normRanks))
  if (length(genes) < 2) return(rep(NA_real_, ncol(normRanks)))
  apply(normRanks[genes, , drop = FALSE], 2, dispersionFun, na.rm = TRUE)
}

simpleScore <- function(rankData, upSet, downSet = NULL, subSamples = NULL,
                        centerScore = TRUE, dispersionFun = stats::mad,
                        knownDirection = TRUE) {
  rankData <- as.matrix(rankData)
  if (!is.null(subSamples)) rankData <- rankData[, subSamples, drop = FALSE]
  n <- nrow(rankData)
  norm <- rankData / (n + 1)
  up <- .score_one(norm, upSet, centerScore)
  upD <- .disp_one(norm, upSet, dispersionFun)
  if (is.null(downSet) || !length(downSet)) {
    TotalScore <- up
    down <- rep(NA_real_, length(up))
    downD <- rep(NA_real_, length(up))
    TotalDispersion <- upD
  } else {
    down <- .score_one(norm, downSet, centerScore)
    downD <- .disp_one(norm, downSet, dispersionFun)
    TotalScore <- up - down
    TotalDispersion <- upD + downD
  }
  data.frame(
    TotalScore = TotalScore,
    TotalDispersion = TotalDispersion,
    UpScore = up,
    UpDispersion = upD,
    DownScore = down,
    DownDispersion = downD,
    row.names = colnames(rankData),
    check.names = FALSE
  )
}

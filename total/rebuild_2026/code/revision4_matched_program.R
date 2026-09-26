# Exploratory, pair-adjusted BCL11B versus developmental-score-matched ETP-like
# expression programs. This does not imply a ZEB-driven causal mechanism.
suppressPackageStartupMessages({library(data.table); library(edgeR); library(limma); library(msigdbr)})
argv <- grep('^--file=',commandArgs(FALSE),value=TRUE)
script <- normalizePath(sub('^--file=','',argv),winslash='/')
root <- normalizePath(file.path(dirname(script),'../../..'),winslash='/')
out <- file.path(root,'total/rebuild_2026/data/source_data_rebuilt/revision4_targeted')
dir.create(out,recursive=TRUE,showWarnings=FALSE)
pairs <- fread(file.path(root,'total/rebuild_2026/data/source_data_rebuilt/revision3_robustness/bcl11b_etp_matched_pairs.tsv'))
stopifnot(nrow(pairs)==18L,uniqueN(c(pairs$BCL11B_sample,pairs$ETP_sample))==36L)
ids <- as.vector(t(as.matrix(pairs[,.(BCL11B_sample,ETP_sample)])))
cn <- names(fread(file.path(root,'data/TALL_X01_counts.tsv'),nrows=0))
stopifnot(all(ids %in% cn))
z <- fread(file.path(root,'data/TALL_X01_counts.tsv'),select=c(cn[1],ids))
gene_id <- z[[1]]
mat <- as.matrix(z[,..ids]); storage.mode(mat)<-'integer'
rownames(mat) <- gene_id
e <- new.env(parent=emptyenv())
load(file.path(root,'data/polonen_syn54032669/Data_1309Samples.RData'),envir=e)
symbols <- rownames(e$gexp)
stopifnot(length(symbols)==nrow(mat),all(nzchar(symbols)))
rm(e,z); gc()
pair <- factor(rep(seq_len(nrow(pairs)),each=2))
group <- factor(rep(c('BCL11B','ETP'),nrow(pairs)),levels=c('ETP','BCL11B'))
design <- model.matrix(~pair+group)
y <- DGEList(counts=mat)
keep <- filterByExpr(y,design=design)
y <- y[keep,,keep.lib.sizes=FALSE]
sym <- symbols[keep]
y <- calcNormFactors(y,method='TMM')
v <- voom(y,design,plot=FALSE)
fit <- eBayes(lmFit(v,design),robust=TRUE)
t <- topTable(fit,coef='groupBCL11B',number=Inf,sort.by='none')
t[,c('gene_id','symbol')] <- list(rownames(t),sym)
res <- as.data.table(t)[,.(gene_id,symbol,logFC,AveExpr,t,P.Value,adj.P.Val)]
setorder(res,P.Value)
fwrite(res,file.path(out,'matched_BCL11B_ETP_gene_results.tsv'),sep='\t')
# Hallmark gene sets, package-pinned and traceable; camera accounts for
# intergene correlation and tests each program against the remaining genes.
gs <- msigdbr(species='Homo sapiens',collection='H')
sets <- split(gs$gene_symbol,gs$gs_name)
indices <- lapply(sets,function(s) which(sym %in% unique(s)))
indices <- indices[lengths(indices)>=10L & lengths(indices)<=500L]
cam <- camera(v,index=indices,design=design,contrast='groupBCL11B')
cam$pathway <- rownames(cam)
cam <- as.data.table(cam)[,.(pathway,NGenes,Direction,PValue,FDR)]
setorder(cam,FDR,PValue)
fwrite(cam,file.path(out,'matched_BCL11B_ETP_hallmark_camera.tsv'),sep='\t')
writeLines(c('Analysis: 18 fixed developmental-score-matched BCL11B/ETP-like pairs.',
  'Counts: source patient count matrix; TMM/voom; pair-blocked limma model.',
  'Contrasts: BCL11B minus ETP-like; exploratory phenotype comparison.',
  'Gene set: MSigDB Hallmark via msigdbr, camera competitive test with intergene correlation.',
  'Neither matching nor enrichment establishes causality or direct ZEB regulation.',
  paste('edgeR',as.character(packageVersion('edgeR'))),
  paste('limma',as.character(packageVersion('limma'))),
  paste('msigdbr',as.character(packageVersion('msigdbr')))),
  file.path(out,'matched_program_methods.txt'))
writeLines(capture.output(sessionInfo()),file.path(out,'matched_program_R_sessionInfo.txt'))
cat('Matched program complete:',nrow(res),'genes;',nrow(cam),'Hallmark sets\n')

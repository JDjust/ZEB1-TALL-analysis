# Post-review sensitivity analyses. All results are exploratory and use the
# already frozen 1,309-patient cohort. Run from any directory with Rscript.
suppressPackageStartupMessages({
  library(data.table); library(edgeR); library(splines)
})
argv <- grep("^--file=", commandArgs(FALSE), value=TRUE)
script <- normalizePath(sub("^--file=", "", argv), winslash="/")
root <- normalizePath(file.path(dirname(script), "../../.."), winslash="/")
out <- file.path(root, "total/rebuild_2026/data/source_data_rebuilt/revision3_robustness")
dir.create(out, recursive=TRUE, showWarnings=FALSE)
read <- function(p) fread(file.path(root,p), encoding="UTF-8")
write <- function(x,n) fwrite(x,file.path(out,n),sep="\t")
d <- read("total/data/validation/zeb_developmental_residual/polonen_patient_residual.tsv")
frozen <- read("total/data/validation/polonen_round1/patient_level_frozen_balance.tsv")
stopifnot(nrow(d)==1309L, uniqueN(d$sample_id)==1309L, uniqueN(d$subtype)==17L)
d <- merge(d, frozen[,.(sample_id, etp)], by="sample_id", sort=FALSE)
stopifnot(nrow(d)==1309L, max(abs(d$dev-(d$z_CD1A-(d$z_CD34+d$z_LYL1)/2)))<1e-8)
subtypes <- sort(unique(d$subtype))

# Group-held-out predictions keep the tested subtype out of its own baseline.
cf <- copy(d[,.(sample_id,subtype,dev,balance)])
cf[, `:=`(expected_crossfit=NA_real_, training_min=NA_real_, training_max=NA_real_)]
for(s in subtypes) {
  train <- d[subtype!=s]; test <- d[subtype==s]
  fit <- lm(balance~ns(dev,df=3),data=train)
  cf[subtype==s, `:=`(expected_crossfit=as.numeric(predict(fit,newdata=test)),
                      training_min=min(train$dev), training_max=max(train$dev))]
}
cf[, residual_crossfit:=balance-expected_crossfit]
cf[, out_of_training_range:=dev<training_min | dev>training_max]
cf_summary <- cf[,.(n=.N,median_crossfit=median(residual_crossfit),
  mean_crossfit=mean(residual_crossfit),n_out_of_training_range=sum(out_of_training_range)),by=subtype]
# Match directly by subtype to avoid dependence on input row order.
primary_medians <- d[,.(subtype,median_primary=median(residual_within)),by=subtype]
cf_summary[,median_primary:=primary_medians$median_primary[match(subtype,primary_medians$subtype)]]
setorder(cf_summary,median_crossfit)
write(cf,"crossfit_patient_residuals.tsv")
write(cf_summary,"crossfit_subtype_summary.tsv")

# Alternative marker definitions are expression proxies, not independent
# measurements of maturation. The expanded set is based on early/cortical
# directions in two normal references and excludes ZEB genes, LMO2 and SPI1.
d[, `:=`(minus_CD34=z_CD1A-z_LYL1,
          minus_LYL1=z_CD1A-z_CD34,
          minus_CD1A=-(z_CD34+z_LYL1)/2)]
extra <- c("MEF2C","PTCRA","RAG1","RAG2","CD1B","LEF1","TRAC")
e <- new.env(parent=emptyenv())
load(file.path(root,"data/polonen_syn54032669/Data_1309Samples.RData"),envir=e)
ids <- rownames(e$annot); sym <- rownames(e$gexp)
ip_raw <- e$M7.IP.factor
if(is.data.frame(ip_raw)) {
  ip <- data.table(sample_id=rownames(ip_raw),ip=as.character(ip_raw[[1]]))
} else {
  ip <- data.table(sample_id=names(ip_raw),ip=as.character(ip_raw))
}
stopifnot(nrow(ip)==1309L, uniqueN(ip$sample_id)==1309L)
counts <- read("data/TALL_X01_counts.tsv")
stopifnot(all(ids %in% names(counts)),nrow(counts)==length(sym))
mat <- as.matrix(counts[,..ids]); storage.mode(mat)<-"double"
y <- DGEList(counts=mat); y <- calcNormFactors(y,method="TMM")
gene_all <- c("CD34","LYL1","CD1A","ZEB1","ZEB2",extra)
hits <- lapply(gene_all,function(g) which(sym==g)); names(hits)<-gene_all
stopifnot(all(lengths(hits)==1L))
small <- cpm(y, log=TRUE, prior.count=1, normalized.lib.sizes=TRUE)[unlist(hits),,drop=FALSE]
rownames(small)<-gene_all
check <- merge(data.table(sample_id=ids,CD34_recalc=as.numeric(small["CD34",])),
               d[,.(sample_id,CD34)],by="sample_id")
stopifnot(max(abs(check$CD34_recalc-check$CD34))<1e-6)
for(g in extra) {
  x <- as.numeric(small[g,match(d$sample_id,ids)])
  d[[paste0("z_",g)]] <- as.numeric(scale(x))
}
d[, expanded:= (z_CD1A+z_PTCRA+z_RAG1+z_RAG2+z_CD1B+z_LEF1+z_TRAC)/7-
                 (z_CD34+z_LYL1+z_MEF2C)/3]
rm(mat,y,small,counts,e); gc()
scores <- c(primary="dev",minus_CD34="minus_CD34",minus_LYL1="minus_LYL1",
            minus_CD1A="minus_CD1A",expanded="expanded")
alt_pat <- data.table()
alt_sum <- data.table()
for(nm in names(scores)) {
  xx <- copy(d[,.(sample_id,subtype,balance)])
  xx[,coordinate:=d[[scores[[nm]]]]]
  fit <- lm(balance~ns(coordinate,df=3),data=xx)
  xx[,residual:=as.numeric(resid(fit))]
  xx[,definition:=nm]
  alt_pat <- rbind(alt_pat,xx)
  sm <- xx[,.(n=.N,median_residual=median(residual)),by=subtype]
  sm[,`:=`(definition=nm,model_r2=summary(fit)$r.squared,
           rank_negative=rank(median_residual,ties.method="min"))]
  alt_sum <- rbind(alt_sum,sm)
}
write(alt_pat,"alternative_coordinates_patient_residuals.tsv")
write(alt_sum,"alternative_coordinates_subtype_summary.tsv")

# Fair model comparison: raw/adjusted and repeated subtype-stratified 5-fold
# predictive R2. Each patient enters one held-out fold per repeat.
v <- merge(d[,.(sample_id,balance,dev,subtype,CD34,LYL1,CD1A,ZEB1,ZEB2)],ip,by="sample_id")
v[,ip:=factor(ip)]
v[,subtype:=factor(subtype,levels=subtypes)]
forms <- list(immunophenotype=balance~ip,development=balance~ns(dev,df=3),
              subtype=balance~subtype,combined=balance~ns(dev,df=3)+subtype)
in_sample <- rbindlist(lapply(names(forms),function(nm) {
  m<-lm(forms[[nm]],data=v); s<-summary(m)
  data.table(model=nm,n=nrow(v),df_model=length(coef(m))-1L,
             raw_r2=s$r.squared,adjusted_r2=s$adj.r.squared)
}))
set.seed(20260926L)
cv <- data.table()
for(rep in 1:5) {
  fold <- integer(nrow(v))
  for(s in subtypes) {
    ii<-which(v$subtype==s)
    fold[ii]<-sample(rep(1:5,length.out=length(ii)))
  }
  # A single author-labelled 'Other' patient cannot define a held-out IP
  # factor level. Retain that patient in each training set for every model.
  fold[as.character(v$ip)=="Other"]<-0L
  for(k in 1:5) {
    tr<-copy(v[fold!=k]); te<-copy(v[fold==k])
    # Re-estimate all expression standardization parameters using training
    # patients only. The held-out patients never set score or outcome scales.
    for(g in c("CD34","LYL1","CD1A","ZEB1","ZEB2")) {
      mu<-mean(tr[[g]]); sig<-sd(tr[[g]])
      tr[[paste0("z_",g)]]<-(tr[[g]]-mu)/sig
      te[[paste0("z_",g)]]<-(te[[g]]-mu)/sig
    }
    tr[,`:=`(dev=z_CD1A-(z_CD34+z_LYL1)/2,balance=z_ZEB1-z_ZEB2)]
    te[,`:=`(dev=z_CD1A-(z_CD34+z_LYL1)/2,balance=z_ZEB1-z_ZEB2)]
    for(nm in names(forms)) {
      m<-lm(forms[[nm]],data=tr)
      p<-as.numeric(predict(m,newdata=te))
      cv<-rbind(cv,data.table(replicate_id=rep,fold=k,model=nm,
                              sample_id=te$sample_id,observed=te$balance,predicted=p))
    }
  }
}
cvsum <- cv[,.(cv_r2=1-sum((observed-predicted)^2)/
                   sum((observed-mean(observed))^2)),by=.(replicate_id,model)]
cvavg <- cvsum[,.(cv_r2_mean=mean(cv_r2),cv_r2_min=min(cv_r2),
                  cv_r2_max=max(cv_r2)),by=model]
in_sample<-merge(in_sample,cvavg,by="model")
full_dev<-lm(forms$development,data=v); full_both<-lm(forms$combined,data=v)
incremental <- summary(full_both)$r.squared-summary(full_dev)$r.squared
partial <- incremental/(1-summary(full_dev)$r.squared)
in_sample[,`:=`(incremental_r2_after_development=ifelse(model=="combined",incremental,NA_real_),
                partial_r2_subtype_given_development=ifelse(model=="combined",partial,NA_real_))]
write(in_sample,"model_comparison.tsv");write(cvsum,"cv_repeats.tsv")

# Exploratory, fixed set of pairwise median-residual contrasts. Bootstrap at
# the patient level within each of the two subtypes; no post-hoc rank P-values.
set.seed(20260926L)
comparators<-c("ETP-like","SPI1",subtypes[grepl("^LMO2",subtypes) & grepl("like$",subtypes)])
contrast<-rbindlist(lapply(comparators,function(s) {
  a<-d[subtype=="BCL11B",residual_within]; b<-d[subtype==s,residual_within]
  boot<-replicate(3000,median(sample(a,length(a),TRUE))-median(sample(b,length(b),TRUE)))
  data.table(comparison=paste("BCL11B vs",s),n_BCL11B=length(a),n_comparator=length(b),
             median_difference=median(a)-median(b),
             ci_low=unname(quantile(boot,.025)),ci_high=unname(quantile(boot,.975)))
}))
write(contrast,"bcl11b_pairwise_bootstrap.tsv")

# One-to-one nearest-neighbour contrast with ETP-like inside common coordinate
# support. Greedy no-replacement matching is deterministic, so label as an
# exploratory matched sensitivity rather than a causal matched analysis.
a<-d[subtype=="BCL11B"][order(dev)]
b<-d[subtype=="ETP-like"]
used<-rep(FALSE,nrow(b)); pairs<-data.table()
for(i in seq_len(nrow(a))) {
  dist<-abs(b$dev-a$dev[i]); dist[used]<-Inf
  j<-which.min(dist); used[j]<-TRUE
  pairs<-rbind(pairs,data.table(BCL11B_sample=a$sample_id[i],ETP_sample=b$sample_id[j],
    BCL11B_dev=a$dev[i],ETP_dev=b$dev[j],distance=dist[j],
    BCL11B_balance=a$balance[i],ETP_balance=b$balance[j],
    BCL11B_residual=a$residual_within[i],ETP_residual=b$residual_within[j]))
}
pairs[,`:=`(balance_difference=BCL11B_balance-ETP_balance,
            residual_difference=BCL11B_residual-ETP_residual)]
write(pairs,"bcl11b_etp_matched_pairs.tsv")
set.seed(20260927L)
boots<-replicate(3000,mean(sample(pairs$residual_difference,nrow(pairs),TRUE)))
match_sum<-data.table(n_pairs=nrow(pairs),median_distance=median(pairs$distance),
  max_distance=max(pairs$distance),mean_residual_difference=mean(pairs$residual_difference),
  ci_low=unname(quantile(boots,.025)),ci_high=unname(quantile(boots,.975)),
  mean_balance_difference=mean(pairs$balance_difference))
write(match_sum,"bcl11b_etp_matched_summary.tsv")
writeLines(capture.output(sessionInfo()),file.path(out,"R_sessionInfo.txt"))
cat("Revision 3 robustness complete; output:",out,"\n")

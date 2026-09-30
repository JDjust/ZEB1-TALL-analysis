# Prespecified, scope-limited sensitivity analyses for journal variants.
suppressPackageStartupMessages({library(data.table); library(splines); library(survival); library(glmnet)})
argv <- grep('^--file=', commandArgs(FALSE), value=TRUE)
script <- normalizePath(sub('^--file=', '', argv), winslash='/')
root <- normalizePath(file.path(dirname(script), '../../..'), winslash='/')
out <- file.path(root, 'total/rebuild_2026/data/source_data_rebuilt/revision4_targeted')
dir.create(out, recursive=TRUE, showWarnings=FALSE)
write <- function(x, name) fwrite(x, file.path(out,name), sep='\t')

# Keep the patient set and the three-gene score fixed while changing only df.
d <- fread(file.path(root,'total/data/validation/zeb_developmental_residual/polonen_patient_residual.tsv'))
pairs <- fread(file.path(root,'total/rebuild_2026/data/source_data_rebuilt/revision3_robustness/bcl11b_etp_matched_pairs.tsv'))
stopifnot(nrow(d)==1309L, uniqueN(d$sample_id)==1309L, nrow(pairs)==18L)
df_rows <- rbindlist(lapply(2:4, function(k) {
  f <- lm(balance ~ ns(dev,df=k),data=d)
  g <- lm(balance ~ ns(dev,df=k)+subtype,data=d)
  r <- copy(d[,.(sample_id,subtype,dev,balance)])
  r[,residual:=as.numeric(residuals(f))]
  a <- r[subtype %in% c('BCL11B','TLX3','ETP-like'),
         .(n=.N,median_residual=median(residual)),by=subtype]
  ix1 <- match(pairs$BCL11B_sample,r$sample_id)
  ix2 <- match(pairs$ETP_sample,r$sample_id)
  a[,`:=`(spline_df=k,
    matched_BCL11B_minus_ETP_mean=mean(r$residual[ix1]-r$residual[ix2]),
    subtype_incremental_r2=summary(g)$r.squared-summary(f)$r.squared,
    developmental_r2=summary(f)$r.squared)]
  a
}))
setcolorder(df_rows,c('spline_df','subtype','n','median_residual','matched_BCL11B_minus_ETP_mean','developmental_r2','subtype_incremental_r2'))
write(df_rows,'spline_df_2_3_4.tsv')

# The measured normal-stage units are descriptive; GSE206710 units are nested
# within donors and its D1 label includes repeated sample-stage aggregates.
n <- fread(file.path(root,'total/rebuild_2026/data/source_data_rebuilt/S3/S3C_normal_stage_units.tsv'))
stopifnot(nrow(n)==20L)
n[, donor:=ifelse(dataset=='GSE206710',sub('\\|.*','',unit),NA_character_)]
normal_summary <- rbindlist(list(
  n[dataset=='GSE195812',.(dataset='GSE195812',n_measured=.N,n_donors=NA_integer_,
    rho_stage_score=cor(order,dev,method='spearman'),
    median_early=median(dev[window=='early']),median_cortical=median(dev[grepl('^DP',stage)]))],
  n[dataset=='GSE206710',.(dataset='GSE206710',n_measured=.N,n_donors=uniqueN(donor),
    rho_stage_score=cor(order,dev,method='spearman'),
    median_early=median(dev[window=='early']),median_cortical=median(dev[grepl('^DP',stage)]))]
),fill=TRUE)
write(normal_summary,'normal_proxy_stage_summary.tsv')
donor_stage <- n[dataset=='GSE206710',.(dev=mean(dev),n_aggregates=.N),by=.(donor,stage,order)]
donor_summary <- donor_stage[,.(n_stages=.N,rho_stage_score=cor(order,dev,method='spearman'),
  early_to_immature_DP=dev[stage=='DP_immature']-dev[stage=='DN_early'],
  early_to_mature_DP=dev[stage=='DP_mature']-dev[stage=='DN_early']),by=donor]
write(donor_stage,'normal_proxy_donor_stage.tsv')
write(donor_summary,'normal_proxy_donor_summary.tsv')

# Ridge Cox: penalize subtype dummies only; keep balance, age, sex and WBC
# unpenalized. CV chooses a common penalty per endpoint before fixed-lambda
# patient bootstrap. This is sensitivity analysis, not a prediction model.
clin <- fread(file.path(root,'total/data/validation/polonen_round1/patient_level_frozen_balance.tsv'))
clin[,subtype:=factor(subtype)]
clin[,sex:=factor(sex)]
cox_rows <- list()
set.seed(20260926L)
for (ep in c('efs','os')) {
  cols <- c(paste0(ep,c('_time','_status')),'balance','subtype','age','sex','log10_wbc')
  x <- clin[complete.cases(clin[,..cols])]
  mm <- model.matrix(~ balance + subtype + age + sex + log10_wbc,data=x)[,-1,drop=FALSE]
  penalty <- ifelse(grepl('^subtype',colnames(mm)),1,0)
  stopifnot(sum(penalty)==nlevels(x$subtype)-1, 'balance' %in% colnames(mm))
  surv <- Surv(x[[paste0(ep,'_time')]],x[[paste0(ep,'_status')]])
  cv <- cv.glmnet(mm,surv,family='cox',alpha=0,penalty.factor=penalty,
                  nfolds=5,type.measure='deviance',standardize=TRUE,cox.ties='breslow')
  lambda <- cv$lambda.min
  fit <- glmnet(mm,surv,family='cox',alpha=0,penalty.factor=penalty,
                lambda=lambda,standardize=TRUE,cox.ties='breslow')
  beta <- as.numeric(coef(fit)['balance',1])
  path_lambdas <- sort(unique(c(cv$lambda.min,cv$lambda.1se,0.01,0.1,1,10)))
  path <- rbindlist(lapply(path_lambdas,function(v) {
    q <- glmnet(mm,surv,family='cox',alpha=0,penalty.factor=penalty,
                lambda=v,standardize=TRUE,cox.ties='breslow')
    cf <- as.numeric(coef(q)[,1]); names(cf) <- rownames(coef(q))
    data.table(endpoint=ep,lambda=v,balance_HR=exp(cf['balance']),
      subtype_coefficient_L2=sqrt(sum(cf[grepl('^subtype',names(cf))]^2)),
      max_abs_subtype_coefficient=max(abs(cf[grepl('^subtype',names(cf))])))
  }))
  write(path,paste0('ridge_cox_penalty_path_',ep,'.tsv'))
  B <- 500L
  boot <- rep(NA_real_,B)
  for (b in seq_len(B)) {
    ix <- sample.int(nrow(x),nrow(x),replace=TRUE)
    xb <- mm[ix,,drop=FALSE]
    yb <- surv[ix]
    z <- try(glmnet(xb,yb,family='cox',alpha=0,penalty.factor=penalty,
                    lambda=lambda,standardize=TRUE,cox.ties='breslow'),silent=TRUE)
    if (!inherits(z,'try-error')) {
      zz <- as.numeric(coef(z)[,1]); names(zz) <- rownames(coef(z))
      # glmnet can return an all-zero placeholder after a convergence warning.
      if (all(is.finite(zz)) && any(abs(zz)>1e-10)) boot[b] <- zz['balance']
    }
  }
  ok <- is.finite(boot)
  stopifnot(sum(ok)>=450L)
  ci <- quantile(exp(boot[ok]),c(.025,.975),names=FALSE)
  cox_rows[[ep]] <- data.table(endpoint=ep,n=nrow(x),events=sum(x[[paste0(ep,'_status')]]),
    method='ridge Cox (Breslow ties); subtype dummies penalized only',lambda_min=lambda,
    lambda_1se=cv$lambda.1se,
    balance_HR=exp(beta),bootstrap_CI_low=ci[1],bootstrap_CI_high=ci[2],
    bootstrap_success=sum(ok),bootstrap_replicates=B,
    bootstrap_fraction_HR_below_1=mean(boot[ok]<0))
}
write(rbindlist(cox_rows),'ridge_cox_sensitivity.tsv')
writeLines(capture.output(warnings()),file.path(out,'R_warnings.txt'))
writeLines(capture.output(sessionInfo()),file.path(out,'R_sessionInfo.txt'))
cat('Revision 4 targeted sensitivity complete:',out,'\n')

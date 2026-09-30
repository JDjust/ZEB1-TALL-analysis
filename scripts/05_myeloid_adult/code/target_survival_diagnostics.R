# Independent survival-package diagnostics on the exported patient-level data.
args <- commandArgs(trailingOnly=TRUE)
root <- if(length(args)) args[1] else "D:/_bioinformation/ZEB1/data/analysis/total"
suppressPackageStartupMessages(library(survival))
out <- file.path(root, "data/validation/target_survival")
dir.create(out, recursive=TRUE, showWarnings=FALSE)
d <- read.delim(file.path(root,"data/F7_target_surv_merged.tsv"), check.names=FALSE)
stopifnot(!anyDuplicated(d[["_k"]]))
# Match the original population-SD scaling; do not refit cutoffs.
popz <- function(x) (x-mean(x,na.rm=TRUE))/sqrt(mean((x-mean(x,na.rm=TRUE))^2,na.rm=TRUE))
d$ZEB1_z <- popz(d$ZEB1)
d$age_z <- popz(d$age_days)
fits <- ph <- follow <- list()
for(ep in c("OS","EFS")) {
  d$time <- d[[paste0(tolower(ep),"_time")]]
  d$event <- d[[paste0(tolower(ep),"_event")]]
  a <- d[complete.cases(d[,c("time","event","ZEB1_z","age_z")]) & !is.na(d$time) & d$time>0,]
  stopifnot(all(a$event %in% c(0,1)))
  revkm <- survfit(Surv(time,1-event)~1,data=a)
  st <- summary(revkm)$table
  follow[[ep]] <- data.frame(endpoint=ep,n=nrow(a),events=sum(a$event),
    reverse_KM_median_days=unname(st["median"]),lo=unname(st["0.95LCL"]),hi=unname(st["0.95UCL"]),
    min_days=min(a$time),max_days=max(a$time))
  for(model in c("ZEB1_unadj","ZEB1_age")) {
    f <- if(model=="ZEB1_unadj") Surv(time,event)~ZEB1_z else Surv(time,event)~ZEB1_z+age_z
    for(ties in c("breslow","efron")) {
      fit <- coxph(f,data=a,ties=ties,x=TRUE,y=TRUE)
      s <- summary(fit)
      fits[[paste(ep,model,ties)]] <- data.frame(endpoint=ep,model=model,ties=ties,term=rownames(s$coefficients),
        n=fit$n,n_event=fit$nevent,hr=s$conf.int[,1],ci_lo=s$conf.int[,3],ci_hi=s$conf.int[,4],p=s$coefficients[,5])
      z <- cox.zph(fit,transform="km")$table
      ph[[paste(ep,model,ties)]] <- data.frame(endpoint=ep,model=model,ties=ties,term=rownames(z),chisq=z[,1],df=z[,2],p=z[,3])
    }
  }
}
write.table(do.call(rbind,fits),file.path(out,"cox_refit.tsv"),sep="\t",row.names=FALSE,quote=FALSE)
write.table(do.call(rbind,ph),file.path(out,"proportional_hazards.tsv"),sep="\t",row.names=FALSE,quote=FALSE)
write.table(do.call(rbind,follow),file.path(out,"followup.tsv"),sep="\t",row.names=FALSE,quote=FALSE)
capture.output(sessionInfo(),file=file.path(out,"R_sessionInfo.txt"))
print(do.call(rbind,follow))
print(do.call(rbind,ph))

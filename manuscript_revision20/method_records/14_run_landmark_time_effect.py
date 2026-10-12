from pathlib import Path
import pandas as pd,numpy as np,json,hashlib,datetime,subprocess
HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[1];OUT=HERE/'landmark_time_effect';OUT.mkdir(exist_ok=True)
contract_path=HERE/'LANDMARK_TIME_EFFECT_CONTRACT.json'
qualification_path=HERE/'POLONEN_TIME_UNIT_QUALIFICATION.json'
assert json.loads(qualification_path.read_text(encoding='utf-8'))['qualified_days']
clinical=ROOT/'02_analyses/03_q1_bioinfo_paper/cox_analysis_table.tsv'
dev_path=ROOT/'tasks/ZEB1_TARGET_MIN_PROBE/inputs/PROBE_PATIENT_SCORES.tsv'
df=pd.read_csv(clinical,sep='\t').merge(pd.read_csv(dev_path,sep='\t')[['sample_id','dev']],on='sample_id',validate='one_to_one')
assert df.subject_sjuid.is_unique and df.subject_sjuid.notna().all()
df['dev_centered']=df.dev-df.dev.mean();df['dev_squared']=df.dev_centered**2
assert set(df.sex)=={'Male','Female'};assert df.protocol.nunique()==2
df['sex_male']=(df.sex=='Male').astype(int);df['protocol_binary']=(df.protocol==sorted(df.protocol.unique())[1]).astype(int)
ledgers={};r_jobs=[]
for ep,duration,event in [('OS','os_time','os_status'),('EFS','efs_time','efs_status')]:
 needed=[duration,event,'ZEB1','dev_centered','dev_squared','age_years','sex_male','log_wbc','protocol_binary','etp3','subtype_g','log_mrd']
 complete=df.dropna(subset=needed)
 eligible=complete.loc[complete[duration]>29].copy();eligible['duration']=eligible[duration]-29;eligible['event']=eligible[event].astype(int)
 early=complete.loc[complete[duration]<=29].copy()
 early.to_csv(OUT/f'{ep}_landmark_early_exclusions.tsv',sep='\t',index=False)
 df.loc[~df.index.isin(complete.index)].to_csv(OUT/f'{ep}_missing_exclusions.tsv',sep='\t',index=False)
 eligible[['sample_id','subject_sjuid',duration,event,'duration','event']].to_csv(OUT/f'{ep}_included_subjects.tsv',sep='\t',index=False)
 ledgers[ep]={'n_complete_before_landmark':len(complete),'n_landmark':len(eligible),'n_events':int(eligible.event.sum()),'n_early_excluded':len(early),'n_early_events':int(early[event].sum()),'n_missing_excluded':len(df)-len(complete),'at_risk':{str(t):int((eligible.duration>=t).sum()) for t in [0,90,365.25,1095.75,1826.25]},'n_future_events':{str(t):int(((eligible.duration>=t)&(eligible.event==1)).sum()) for t in [0,90,365.25,1095.75,1826.25]}}
 cols=['duration','event','ZEB1','age_years','sex_male','log_wbc','protocol_binary','etp3','subtype_g','log_mrd','dev_centered','dev_squared']
 csv=eligible[cols].to_csv(sep='\t',index=False);(OUT/f'{ep}_MODEL_INPUT.tsv').write_text(csv,encoding='utf-8')
 r_jobs.append('data <- read.delim(textConnection('+json.dumps(csv,ensure_ascii=False)+'),check.names=FALSE)\nresults[["'+ep+'"]] <- analyze(data,"'+ep+'")')
r_code=r'''suppressPackageStartupMessages(library(survival))
suppressPackageStartupMessages(library(jsonlite))
base <- "age_years + sex_male + log_wbc + protocol_binary + etp3 + log_mrd + dev_centered + dev_squared + strata(subtype_g)"
ttfun <- function(x,t,...) x*log1p(t/365.25)
fit <- function(data,formula){
 warns <- character()
 value <- withCallingHandlers(coxph(as.formula(paste("Surv(duration,event) ~",formula)),data=data,ties="efron",tt=ttfun,x=TRUE),warning=function(w){warns <<- c(warns,conditionMessage(w));invokeRestart("muffleWarning")})
 list(fit=value,warnings=warns)
}
summarize <- function(obj){
 z <- summary(obj$fit)$coefficients
 list(coefficients=data.frame(term=rownames(z),z,check.names=FALSE,row.names=NULL),warnings=obj$warnings,log_likelihood=as.numeric(logLik(obj$fit)),n_model_rows=obj$fit$n,events=obj$fit$nevent)
}
lrt <- function(full,reduced){
 statistic <- 2*(as.numeric(logLik(full$fit))-as.numeric(logLik(reduced$fit)))
 df <- attr(logLik(full$fit),"df")-attr(logLik(reduced$fit),"df")
 list(statistic=statistic,df=df,p=pchisq(statistic,df,lower.tail=FALSE))
}
analyze <- function(data,endpoint){
 static_reduced <- fit(data,base)
 static_full <- fit(data,paste(base,"+ ZEB1"))
 extended_base <- paste(base,"+ tt(log_mrd) + tt(sex_male)")
 extended_reduced <- fit(data,extended_base)
 extended_constant <- fit(data,paste(extended_base,"+ ZEB1"))
 extended_full <- fit(data,paste(extended_base,"+ ZEB1 + tt(ZEB1)"))
 ph <- cox.zph(static_full$fit,transform="rank")
 coeff <- coef(extended_full$fit);covariance <- vcov(extended_full$fit)
 hr <- lapply(c(0,90,365.25,1095.75,1826.25),function(t){
  v <- setNames(rep(0,length(coeff)),names(coeff));v['ZEB1'] <- 1;v['tt(ZEB1)'] <- log1p(t/365.25)
  loghr <- sum(v*coeff);se <- sqrt(as.numeric(t(v)%*%covariance%*%v))
  list(elapsed_days=t,HR=exp(loghr),lower95=exp(loghr-1.96*se),upper95=exp(loghr+1.96*se))
 })
 list(n_original_subjects=nrow(data),events=sum(data$event),static_reduced=summarize(static_reduced),static_full=summarize(static_full),static_LRT=lrt(static_full,static_reduced),static_PH=data.frame(term=rownames(ph$table),ph$table,row.names=NULL),extended_reduced=summarize(extended_reduced),extended_constant=summarize(extended_constant),extended_full=summarize(extended_full),ZEB1_block_LRT=lrt(extended_full,extended_reduced),ZEB1_time_increment_LRT=lrt(extended_full,extended_constant),HR_by_elapsed_time=hr)
}
results <- list()
'''
r_code+='\n'.join(r_jobs)
r_code+='\ncat("RESULT_BEGIN\\n");cat(toJSON(list(environment=list(R=R.version.string,survival=as.character(packageVersion("survival"))),results=results),auto_unbox=TRUE,digits=16,na="null"));cat("\\nRESULT_END\\n")\n'
script=OUT/'landmark_time_effect.R';script.write_text(r_code,encoding='utf-8')
(OUT/'LANDMARK_PATIENT_LEDGER.json').write_text(json.dumps(ledgers,indent=2),encoding='utf-8')
print('PATIENT_LEDGER',json.dumps(ledgers),flush=True)
process=subprocess.run([r'D:\Git\Git\usr\bin\ssh.exe','-o','BatchMode=yes','-o','ConnectTimeout=15','hulu','/opt/bioinfo/envs/r-bio/bin/Rscript','-'],input=r_code,capture_output=True,text=True,encoding='utf-8',timeout=240)
(OUT/'REMOTE_MODEL_LOG.txt').write_text(process.stdout+'\nSTDERR\n'+process.stderr,encoding='utf-8')
if process.returncode:raise RuntimeError(process.stderr[-4000:])
payload=process.stdout.split('RESULT_BEGIN\n',1)[1].split('\nRESULT_END',1)[0];result=json.loads(payload)
result['source_sha256']={str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in [clinical,dev_path,contract_path,qualification_path,script]}
result['patient_ledger']=ledgers;result['created_utc']=datetime.datetime.now(datetime.timezone.utc).isoformat()
from statsmodels.stats.multitest import multipletests
q=multipletests([result['results'][ep]['ZEB1_block_LRT']['p'] for ep in ['OS','EFS']],method='fdr_bh')[1]
for ep,adj in zip(['OS','EFS'],q):
 result['results'][ep]['ZEB1_block_LRT']['BH_q_two_endpoints']=float(adj)
 info=result['results'][ep]
 print(ep,'STATIC_LRT',info['static_LRT'],'ZEB1_BLOCK',info['ZEB1_block_LRT'],'TIME_INCREMENT',info['ZEB1_time_increment_LRT'],'HR_TIME',info['HR_by_elapsed_time'],'WARNINGS',info['extended_full']['warnings'],flush=True)
(OUT/'RESULTS.json').write_text(json.dumps(result,indent=2,ensure_ascii=False),encoding='utf-8')

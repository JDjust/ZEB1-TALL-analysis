from pathlib import Path
import pandas as pd,numpy as np,json,hashlib,subprocess,datetime
from statsmodels.stats.multitest import multipletests
HERE=Path(__file__).resolve().parent;OUT=HERE/'multigene_maturity_proxy';qualification=OUT/'QUALIFICATION.json';contract=HERE/'MULTIGENE_MATURITY_PROXY_CONTRACT.json'
assert json.loads(qualification.read_text())['status']=='QUALIFIED_FOR_EXPLORATORY_CLINICAL_SENSITIVITY'
scorepath=OUT/'ALL1309_MULTIGENE_SCORES.tsv';scores=pd.read_csv(scorepath,sep='\t');assert scores.sample_id.is_unique and len(scores)==1309
jobs=[];inputs=[scorepath,qualification,contract];ledgers={}
for ep in ['OS','EFS']:
    oldinput=HERE/'landmark_time_effect'/f'{ep}_MODEL_INPUT.tsv';oldledger=HERE/'landmark_time_effect'/f'{ep}_included_subjects.tsv';inputs.extend([oldinput,oldledger])
    data=pd.read_csv(oldinput,sep='\t');ledger=pd.read_csv(oldledger,sep='\t');assert len(data)==len(ledger) and ledger.subject_sjuid.is_unique
    assert np.allclose(data.duration,ledger.duration) and np.array_equal(data.event,ledger.event)
    indexed=pd.concat([ledger[['sample_id','subject_sjuid']],data],axis=1).merge(scores,on='sample_id',how='left',validate='one_to_one',sort=False)
    assert np.array_equal(indexed.subject_sjuid,ledger.subject_sjuid) and np.isfinite(indexed[['cortical_vs_early','SP_vs_cortical']]).all().all()
    indexed['cortical_squared']=indexed.cortical_vs_early**2;indexed['SP_squared']=indexed.SP_vs_cortical**2
    indexed.to_csv(OUT/f'{ep}_PATIENT_MODEL_INPUT.tsv',sep='\t',index=False)
    ledger.to_csv(OUT/f'{ep}_included_subjects.tsv',sep='\t',index=False)
    ledgers[ep]={'n_patients':len(indexed),'events':int(indexed.event.sum()),'same_order_as_existing_landmark':True,'missing_proxy_patients':0,'ZEB1_score_correlations':indexed[['ZEB1','cortical_vs_early','SP_vs_cortical','dev_centered']].corr(method='spearman').to_dict()}
    jobs.append('data<-read.delim(textConnection('+json.dumps(indexed.to_csv(sep='\t',index=False),ensure_ascii=False)+'),check.names=FALSE)\nresults[["'+ep+'"]]<-analyze(data)')
code=r'''suppressPackageStartupMessages(library(survival));suppressPackageStartupMessages(library(jsonlite))
base<-"age_years + sex_male + log_wbc + protocol_binary + etp3 + log_mrd + cortical_vs_early + cortical_squared + SP_vs_cortical + SP_squared + strata(subtype_g)"
ttfun<-function(x,t,...) x*log1p(t/365.25)
fit<-function(data,formula){
 warns<-character();value<-withCallingHandlers(coxph(as.formula(paste("Surv(duration,event) ~",formula)),data=data,ties="efron",tt=ttfun,x=TRUE),warning=function(w){warns<<-c(warns,conditionMessage(w));invokeRestart("muffleWarning")});list(fit=value,warnings=warns)
}
summarize<-function(obj){
 z<-summary(obj$fit)$coefficients;ci<-summary(obj$fit)$conf.int
 list(coefficients=data.frame(term=rownames(z),z,lower95=ci[,"lower .95"],upper95=ci[,"upper .95"],check.names=FALSE,row.names=NULL),warnings=obj$warnings,log_likelihood=as.numeric(logLik(obj$fit)),df=attr(logLik(obj$fit),"df"),n_model_rows=obj$fit$n,events=obj$fit$nevent)
}
lrt<-function(full,reduced){statistic<-2*(as.numeric(logLik(full$fit))-as.numeric(logLik(reduced$fit)));df<-attr(logLik(full$fit),"df")-attr(logLik(reduced$fit),"df");list(statistic=statistic,df=df,p=pchisq(statistic,df,lower.tail=FALSE))}
analyze<-function(data){
 sr<-fit(data,base);sf<-fit(data,paste(base,"+ ZEB1"));eb<-paste(base,"+ tt(log_mrd) + tt(sex_male)")
 er<-fit(data,eb);ec<-fit(data,paste(eb,"+ ZEB1"));ef<-fit(data,paste(eb,"+ ZEB1 + tt(ZEB1)"))
 ph<-cox.zph(sf$fit,transform="rank");coef<-coef(ef$fit);cov<-vcov(ef$fit)
 hr<-lapply(c(0,90,365.25,1095.75,1826.25),function(t){v<-setNames(rep(0,length(coef)),names(coef));v['ZEB1']<-1;v['tt(ZEB1)']<-log1p(t/365.25);lh<-sum(v*coef);se<-sqrt(as.numeric(t(v)%*%cov%*%v));list(elapsed_days=t,HR=exp(lh),lower95=exp(lh-1.96*se),upper95=exp(lh+1.96*se))})
 design<-model.matrix(~ZEB1+cortical_vs_early+cortical_squared+SP_vs_cortical+SP_squared,data=data)
 list(n_patients=nrow(data),events=sum(data$event),static_reduced=summarize(sr),static_full=summarize(sf),static_LRT=lrt(sf,sr),static_PH=data.frame(term=rownames(ph$table),ph$table,row.names=NULL),extended_reduced=summarize(er),extended_constant=summarize(ec),extended_full=summarize(ef),ZEB1_block_LRT=lrt(ef,er),ZEB1_time_increment_LRT=lrt(ef,ec),HR_by_elapsed_time=hr,maturity_design_condition_number=kappa(design,exact=TRUE))
}
results<-list()
'''+ '\n'.join(jobs)
code+='\ncat("RESULT_BEGIN\\n");cat(toJSON(list(environment=list(R=R.version.string,survival=as.character(packageVersion("survival"))),results=results),auto_unbox=TRUE,digits=16,na="null"));cat("\\nRESULT_END\\n")\n'
script=OUT/'MULTIGENE_CLINICAL_SENSITIVITY.R';script.write_text(code,encoding='utf-8');print('MULTIGENE_COX_START',json.dumps({ep:{k:v for k,v in info.items() if k!='ZEB1_score_correlations'} for ep,info in ledgers.items()}),flush=True)
p=subprocess.run([r'D:\Git\Git\usr\bin\ssh.exe','-o','BatchMode=yes','-o','ConnectTimeout=15','hulu','/opt/bioinfo/envs/r-bio/bin/Rscript','-'],input=code,capture_output=True,text=True,encoding='utf-8',timeout=240)
(OUT/'CLINICAL_REMOTE_LOG.txt').write_text(p.stdout+'\nSTDERR\n'+p.stderr,encoding='utf-8')
if p.returncode:raise RuntimeError(p.stderr[-3500:])
result=json.loads(p.stdout.split('RESULT_BEGIN\n',1)[1].split('\nRESULT_END',1)[0]);q=multipletests([result['results'][ep]['ZEB1_block_LRT']['p'] for ep in ['OS','EFS']],method='fdr_bh')[1]
for ep,qq in zip(['OS','EFS'],q):
    obj=result['results'][ep];obj['ZEB1_block_LRT']['BH_q_two_endpoints']=float(qq)
    for name in ['static_reduced','static_full','extended_reduced','extended_constant','extended_full']:pd.DataFrame(obj[name]['coefficients']).to_csv(OUT/f'{ep}_{name}_coefficients.tsv',sep='\t',index=False)
    pd.DataFrame(obj['static_PH']).to_csv(OUT/f'{ep}_STATIC_PH.tsv',sep='\t',index=False)
    print('MULTIGENE_COX_RESULT',ep,'static_ZEB1',[x for x in obj['static_full']['coefficients'] if x['term']=='ZEB1'],'primary_block',obj['ZEB1_block_LRT'],'timeincrement',obj['ZEB1_time_increment_LRT'],'warnings',{x:obj[x]['warnings'] for x in ['static_full','extended_reduced','extended_full']},'condition',obj['maturity_design_condition_number'],flush=True)
result.update(patient_ledger=ledgers,source_sha256={str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in inputs+[Path(__file__),script]},created_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),interpretation='Exploratory same-cohort measured-proxy sensitivity, normal-reference multi-gene similarities qualified by internal donor holdout. No complete developmental trajectory, external survival replication, causal mechanism or validated risk model claim.')
(OUT/'CLINICAL_RESULTS.json').write_text(json.dumps(result,indent=2),encoding='utf-8');print('MULTIGENE_CLINICAL_SENSITIVITY_COMPLETE',flush=True)

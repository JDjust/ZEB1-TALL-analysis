from pathlib import Path
import json,subprocess,datetime,hashlib
import pandas as pd
from statsmodels.stats.multitest import multipletests
root=Path(__file__).resolve().parent
dest=root/'author_BMP17_conditioned_clinical_sensitivity'
contract=json.loads((dest/'PREVALUE_CONTRACT.json').read_text(encoding='utf-8'))
qualification=json.loads((dest/'INPUT_QUALIFICATION.json').read_text(encoding='utf-8'))
assert qualification['status']=='QUALIFIED_SAME_PATIENT_AUTHOR_BMP17_CONDITIONING'
code=r'''suppressPackageStartupMessages(library(survival));suppressPackageStartupMessages(library(jsonlite))
base<-"age_years + sex_male + log_wbc + protocol_binary + etp3 + log_mrd + cortical_vs_early + cortical_squared + SP_vs_cortical + SP_squared + strata(subtype_g) + blast_fraction_centered + blast_fraction_squared + tt(log_mrd) + tt(sex_male)"
ttfun<-function(x,t,...)x*log1p(t/365.25)
fitmodel<-function(d,f){warnings<-character();m<-withCallingHandlers(coxph(as.formula(paste("Surv(duration,event)~",f)),data=d,ties="efron",tt=ttfun,x=TRUE),warning=function(w){warnings<<-c(warnings,conditionMessage(w));invokeRestart("muffleWarning")});list(fit=m,warnings=warnings)}
summarize<-function(obj){s<-summary(obj$fit);a<-s$coefficients;ci<-s$conf.int;list(coefficients=data.frame(term=rownames(a),a,lower95=ci[,"lower .95"],upper95=ci[,"upper .95"],check.names=FALSE,row.names=NULL),warnings=obj$warnings,log_likelihood=as.numeric(logLik(obj$fit)),df=attr(logLik(obj$fit),"df"),n_model_rows=obj$fit$n,events=obj$fit$nevent,coef_names=names(coef(obj$fit)),covariance=unname(vcov(obj$fit)),all_coefficients_finite=all(is.finite(coef(obj$fit))))}
lrt<-function(f,r){statistic<-2*(as.numeric(logLik(f$fit))-as.numeric(logLik(r$fit)));df<-attr(logLik(f$fit),"df")-attr(logLik(r$fit),"df");list(statistic=statistic,df=df,p=pchisq(statistic,df,lower.tail=FALSE))}
hrtime<-function(obj){b<-coef(obj$fit);cov<-vcov(obj$fit);lapply(c(0,90,365.25,1095.75,1826.25),function(days){v<-setNames(rep(0,length(b)),names(b));v['ZEB1']<-1;v['tt(ZEB1)']<-log1p(days/365.25);estimate<-sum(v*b);se<-sqrt(as.numeric(t(v)%*%cov%*%v));list(elapsed_landmark_days=days,logHR=estimate,SE=se,HR=exp(estimate),lower95=exp(estimate-1.96*se),upper95=exp(estimate+1.96*se))})}
analyze<-function(d){pre.r<-fitmodel(d,base);pre.f<-fitmodel(d,paste(base,"+ ZEB1 + tt(ZEB1)"));postbase<-paste(base,"+ BMP17_z + BMP17_squared + tt(BMP17_z)");post.r<-fitmodel(d,postbase);post.c<-fitmodel(d,paste(postbase,"+ ZEB1"));post.f<-fitmodel(d,paste(postbase,"+ ZEB1 + tt(ZEB1)"));staticbase<-"age_years + sex_male + log_wbc + protocol_binary + etp3 + log_mrd + cortical_vs_early + cortical_squared + SP_vs_cortical + SP_squared + strata(subtype_g) + blast_fraction_centered + blast_fraction_squared + BMP17_z + BMP17_squared + ZEB1";static<-fitmodel(d,staticbase);ph<-cox.zph(static$fit,transform="rank");design<-model.matrix(~ZEB1+cortical_vs_early+cortical_squared+SP_vs_cortical+SP_squared+blast_fraction_centered+blast_fraction_squared+BMP17_z+BMP17_squared,data=d);list(n_patients=nrow(d),unique_patients=length(unique(d$subject_sjuid)),events=sum(d$event),samecase_pre_reduced=summarize(pre.r),samecase_pre_full=summarize(pre.f),post_reduced=summarize(post.r),post_constant=summarize(post.c),post_full=summarize(post.f),static_full=summarize(static),static_PH=data.frame(term=rownames(ph$table),ph$table,row.names=NULL),samecase_pre_ZEB1_block=lrt(pre.f,pre.r),primary_post_ZEB1_block=lrt(post.f,post.r),post_ZEB1_time_increment=lrt(post.f,post.c),samecase_pre_HR_by_time=hrtime(pre.f),post_HR_by_time=hrtime(post.f),continuous_design_condition_number=kappa(design,exact=TRUE),descriptive_continuous_correlations=unname(cor(d[,c("ZEB1","BMP17_z","cortical_vs_early","SP_vs_cortical")])),correlation_variable_order=c("ZEB1","BMP17_z","cortical_vs_early","SP_vs_cortical"))}
results<-list()
'''
for ep in ['OS','EFS']:
    p=dest/(ep+'_FROZEN_CONDITIONED_INPUT.tsv')
    code+='d<-read.delim(textConnection('+json.dumps(p.read_text(encoding='utf-8'),ensure_ascii=False)+'),check.names=FALSE)\nresults[["'+ep+'"]]<-analyze(d)\n'
code+='cat("BMP17_RESULT_BEGIN\\n");cat(toJSON(list(environment=list(R=R.version.string,survival=as.character(packageVersion("survival"))),results=results),auto_unbox=TRUE,digits=16,na="null"));cat("\\nBMP17_RESULT_END\\n")\n'
rp=dest/'AUTHOR_BMP17_CONDITIONED_MODELS.R';rp.write_text(code,encoding='utf-8')
attempts=[];result=None
for attempt in range(2):
    r=subprocess.run([r'D:\Git\Git\usr\bin\ssh.exe','-o','BatchMode=yes','-o','ConnectTimeout=15','hulu','/opt/bioinfo/envs/r-bio/bin/Rscript','-'],input=code,encoding='utf-8',capture_output=True,timeout=240)
    (dest/f'MODEL_REMOTE_ATTEMPT{attempt+1}.log').write_text(r.stdout+'\nSTDERR\n'+r.stderr,encoding='utf-8');attempts.append({'attempt':attempt+1,'returncode':r.returncode,'stderr_tail':r.stderr[-1500:]})
    if r.returncode==0:
        result=json.loads(r.stdout.split('BMP17_RESULT_BEGIN\n',1)[1].split('\nBMP17_RESULT_END',1)[0]);break
    print('REMOTE_ATTEMPT_FAILED',attempts[-1],flush=True)
    if 'Connection' not in r.stderr and 'ssh' not in r.stderr:break
(dest/'MODEL_EXECUTION_RECEIPT.json').write_text(json.dumps({'created_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'attempts':attempts},ensure_ascii=False,indent=2),encoding='utf-8')
if result is None:raise RuntimeError('No successful model run; inspect preserved source/logs.')
q=multipletests([result['results'][ep]['primary_post_ZEB1_block']['p'] for ep in ['OS','EFS']],method='fdr_bh')[1]
for ep,qq in zip(['OS','EFS'],q):
    x=result['results'][ep];assert [x['n_patients'],x['events']]==contract['retain_population'][ep] and x['unique_patients']==x['n_patients'];x['primary_post_ZEB1_block']['BH2_q']=float(qq)
    for name in ['samecase_pre_reduced','samecase_pre_full','post_reduced','post_constant','post_full','static_full']:pd.DataFrame(x[name]['coefficients']).to_csv(dest/(ep+'_'+name+'_coefficients.tsv'),sep='\t',index=False)
    pd.DataFrame(x['static_PH']).to_csv(dest/(ep+'_STATIC_PH.tsv'),sep='\t',index=False)
    print('AUTHOR_BMP17_CONDITIONED_RESULT',ep,'patients/events',[x['n_patients'],x['events']],'ZEB1_block',x['primary_post_ZEB1_block'],'time_increment',x['post_ZEB1_time_increment'],'1year_HR',x['post_HR_by_time'][2],'warnings',x['post_full']['warnings'],flush=True)
result.update(created_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),contract=contract,input_qualification=qualification,scope='Same cohort post-result exploratory sensitivity; not independentOS, causalproof or predictionvalidation.',source_sha256={str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in [dest/'PREVALUE_CONTRACT.json',dest/'INPUT_QUALIFICATION.json',rp,dest/'OS_FROZEN_CONDITIONED_INPUT.tsv',dest/'EFS_FROZEN_CONDITIONED_INPUT.tsv',Path(__file__)]})
(dest/'RESULTS.json').write_text(json.dumps(result,ensure_ascii=False,indent=2),encoding='utf-8')

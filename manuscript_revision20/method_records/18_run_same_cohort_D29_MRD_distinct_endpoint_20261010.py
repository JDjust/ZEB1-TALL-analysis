from pathlib import Path
from datetime import datetime,timezone
import json,hashlib,subprocess,sys
import pandas as pd,numpy as np,patsy,statsmodels.api as sm,statsmodels
from statsmodels.stats.multitest import multipletests
R=Path(__file__).resolve().parent;O=R/'same_cohort_D29_MRD_distinct_endpoint'
cp=O/'PREVALUE_CONTRACT.json';ip=O/'ALL_FOUR_FROZEN_SAME_CASE_MRD_INPUT.tsv'
qual=json.loads((O/'INPUT_QUALIFICATION.json').read_text(encoding='utf-8'))
assert hashlib.sha256(cp.read_bytes()).hexdigest()==qual['contract_sha256']
assert hashlib.sha256(ip.read_bytes()).hexdigest()==qual['input_sha256']
d=pd.read_csv(ip,sep='\t');assert len(d)==qual['n_complete_same_cases'] and d.subject_sjuid.is_unique
base='ZEB1 + age_years + sex_male + log_wbc + protocol_binary + C(etp3) + C(subtype_g) + cortical_vs_early + cortical_squared + SP_vs_cortical + SP_squared + blast_fraction_centered + blast_fraction_squared'
extra=' + BMP17_z + BMP17_squared + ancestry_native_american_25pct + ancestry_eastasian_25pct + ancestry_african_25pct + ancestry_southasian_25pct'
results=[];allcoeff=[];allloo=[];jobs=[]
for context,formula in [('fixed_maturity_clinical_blast',base),('plus_BMP17_ancestry',base+extra)]:
    X=patsy.dmatrix(formula,d,return_type='dataframe');a=X.to_numpy();n,k=a.shape
    assert np.isfinite(a).all() and np.linalg.matrix_rank(a)==k
    X.to_csv(O/(context+'_DESIGN_MATRIX.tsv'),sep='\t',index=False)
    y=d.MRD_log1p_native_percent.to_numpy()
    fit=sm.OLS(y,X).fit();robust=fit.get_robustcov_results(cov_type='HC3',use_t=True)
    j=X.columns.get_loc('ZEB1');ci=robust.conf_int(alpha=.05)[j];beta=float(robust.params[j]);p=float(robust.pvalues[j])
    inv=np.linalg.inv(a.T@a);res=y-a@fit.params;hat=np.einsum('ij,jk,ik->i',a,inv,a);loo=float(fit.params.iloc[j])-(a@inv[:,j])*res/(1-hat)
    for pid,v in zip(d.subject_sjuid,loo):allloo.append({'context':context,'excluded_subject_sjuid':pid,'analytic_OLS_ZEB1_beta':float(v)})
    results.append({'endpoint':'continuous_ln1p_MRD_percent','context':context,'n_patients':n,'MRD_ge_0p01_percent_n':qual['MRD_ge_0p01_percent_n'],'target_beta':beta,'target_SE':float(robust.bse[j]),'CI95_low':float(ci[0]),'CI95_high':float(ci[1]),'P':p,'CI_method':'HC3 residual-df t Wald','residual_df':float(robust.df_resid),'n_coefficients':k,'design_rank':int(np.linalg.matrix_rank(a)),'design_condition_number':float(np.linalg.cond(a)),'max_leverage':float(hat.max()),'all_patient_leave_one_out_beta_min':float(loo.min()),'all_patient_leave_one_out_beta_max':float(loo.max()),'all_patient_leave_one_out_same_sign':bool(np.all(np.sign(loo)==np.sign(beta))),'statsmodels_version':statsmodels.__version__})
    for l,b,se,pv,c in zip(X.columns,robust.params,robust.bse,robust.pvalues,robust.conf_int()):allcoeff.append({'endpoint':'continuous_ln1p_MRD_percent','context':context,'term':l,'coef':float(b),'SE':float(se),'P':float(pv),'CI95_low':float(c[0]),'CI95_high':float(c[1])})
    cov=pd.DataFrame(robust.cov_params(),index=X.columns,columns=X.columns);cov.to_csv(O/(context+'_OLS_HC3_COVARIANCE.tsv'),sep='\t')
    jobs.append({'context':context,'design':a.tolist(),'column_labels':X.columns.tolist(),'y':d.MRD_positive_ge_0p01_percent.astype(int).tolist(),'target_index_one_based':int(j+1)})
    print('OLS_COMPLETED',context,'N',n,'K',k,flush=True)
# Firth mean-bias-reduced fits on the identical numeric designs, using the maintained brglm2 fitter.
payload=json.dumps(jobs,ensure_ascii=False,separators=(',',':'))
rcode='''.libPaths(c("/home/liangfuhua/zeb1_maturation_clinical_strengthen_20261009/runtime_R_library",.libPaths()))
suppressPackageStartupMessages(library(jsonlite))
suppressPackageStartupMessages(library(brglm2))
jobs <- fromJSON(INPUT_JSON, simplifyVector=FALSE)
out <- list()
for(job in jobs){
 X <- do.call(rbind,lapply(job$design,unlist))
 y <- as.integer(unlist(job$y))
 labels <- unlist(job$column_labels)
 dat <- as.data.frame(X[,-1,drop=FALSE]);names(dat)<-paste0("x",seq_len(ncol(dat)));dat$y<-y
 warns <- character()
 fit <- withCallingHandlers(glm(y~.,data=dat,family=binomial(link="logit"),method=brglm2::brglmFit,type="AS_mean",control=glm.control(maxit=1000,epsilon=1e-10)),warning=function(w){warns<<-c(warns,conditionMessage(w));invokeRestart("muffleWarning")})
 b<-coef(fit);V<-vcov(fit);se<-sqrt(diag(V));z<-b/se;pv<-2*pnorm(abs(z),lower.tail=FALSE);ci<-cbind(b-qnorm(.975)*se,b+qnorm(.975)*se)
 names(b)<-labels;dimnames(V)<-list(labels,labels)
 target<-as.integer(job$target_index_one_based)
 out[[job$context]]<-list(context=job$context,n_patients=length(y),events=sum(y),n_coefficients=length(b),rank=fit$rank,converged=fit$converged,iterations=fit$iter,warnings=warns,target_beta=unname(b[target]),target_SE=unname(se[target]),CI95_low=unname(ci[target,1]),CI95_high=unname(ci[target,2]),P=unname(pv[target]),coefficients=data.frame(term=labels,coef=unname(b),SE=unname(se),P=unname(pv),CI95_low=unname(ci[,1]),CI95_high=unname(ci[,2])),covariance=unname(V),fitted_probability=unname(fitted(fit)))
}
cat(toJSON(list(environment=list(R=R.version.string,brglm2=as.character(packageVersion("brglm2"))),results=out),auto_unbox=TRUE,digits=17,null="null"))
'''.replace('INPUT_JSON',json.dumps(payload,ensure_ascii=False))
rp=O/'FROZEN_TWO_FIRTH_MRD_DESIGNS_AND_FITS.R';rp.write_text(rcode,encoding='utf-8')
cmd=['D:/Git/Git/usr/bin/ssh.exe','-o','BatchMode=yes','-o','ConnectTimeout=15','hulu','/opt/bioinfo/envs/r-bio/bin/Rscript -']
process=subprocess.run(cmd,input=rcode,text=True,encoding='utf-8',capture_output=True,timeout=180)
if (O/'FIRTH_REMOTE_STDERR.txt').exists() and 'no package called' in (O/'FIRTH_REMOTE_STDERR.txt').read_text(encoding='utf-8'):
    (O/'FAILED_INITIAL_MISSING_LIBRARY_STDERR.txt').write_text((O/'FIRTH_REMOTE_STDERR.txt').read_text(encoding='utf-8'),encoding='utf-8')
(O/'FIRTH_REMOTE_STDOUT.json').write_text(process.stdout,encoding='utf-8');(O/'FIRTH_REMOTE_STDERR.txt').write_text(process.stderr,encoding='utf-8')
if process.returncode!=0:
    print('REMOTE_FITTER_FAILED',process.returncode,process.stderr[-2500:]);sys.exit(1)
firth=json.loads(process.stdout)
for context,x in firth['results'].items():
    assert x['converged'] and x['rank']==x['n_coefficients'] and not x['warnings']
    assert x['n_patients']==len(d) and x['events']==qual['MRD_ge_0p01_percent_n']
    results.append({'endpoint':'binary_MRD_ge_0p01_percent','context':context,'n_patients':x['n_patients'],'MRD_ge_0p01_percent_n':x['events'],'target_beta':x['target_beta'],'target_SE':x['target_SE'],'CI95_low':x['CI95_low'],'CI95_high':x['CI95_high'],'OR':float(np.exp(x['target_beta'])),'OR_CI95_low':float(np.exp(x['CI95_low'])),'OR_CI95_high':float(np.exp(x['CI95_high'])),'P':x['P'],'CI_method':'Mean bias-reduced logit Wald normal','n_coefficients':x['n_coefficients'],'design_rank':x['rank'],'converged':x['converged'],'iterations':x['iterations'],'warnings':x['warnings']})
    for row in x['coefficients']:allcoeff.append({'endpoint':'binary_MRD_ge_0p01_percent','context':context,**row})
    pd.DataFrame(x['covariance'],index=[v['term'] for v in x['coefficients']],columns=[v['term'] for v in x['coefficients']]).to_csv(O/(context+'_FIRTH_INFORMATION_COVARIANCE.tsv'),sep='\t')
qs=multipletests([r['P'] for r in results],method='fdr_bh')[1]
for row,q in zip(results,qs):row['BH4_q']=float(q)
out={'created_utc':datetime.now(timezone.utc).isoformat(),'contract_sha256':qual['contract_sha256'],'input_sha256':qual['input_sha256'],'R_environment':firth['environment'],'results':results,'evidence_role':'New exploratory same-cohort response endpoint; not independent OS or causal/temporal uncoupling.','new_OLS_models':2,'new_bias_reduced_logit_models':2,'accepted_Cox_models_rerun':0,'formal_manuscript_or_figure_sync':False}
(O/'RESULTS.json').write_text(json.dumps(out,indent=2),encoding='utf-8');pd.DataFrame(results).to_csv(O/'ALL_FOUR_MRD_TARGET_RESULTS.tsv',sep='\t',index=False);pd.DataFrame(allcoeff).to_csv(O/'ALL_FOUR_MODEL_COEFFICIENTS.tsv',sep='\t',index=False);pd.DataFrame(allloo).to_csv(O/'ALL_PATIENT_ANALYTIC_OLS_LEAVE_ONE_OUT.tsv',sep='\t',index=False)
print('NEW_SAME_COHORT_MRD_RESULTS',json.dumps(results),flush=True)

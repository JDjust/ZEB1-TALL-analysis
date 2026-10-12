from pathlib import Path
import json,hashlib,datetime
import pandas as pd,numpy as np
from scipy.stats import rankdata
H=Path(__file__).resolve().parent;S=H/'independent_AIIMS35_maturation';O=H/'independent_AIIMS35_native_state_conditional_maturation';O.mkdir(exist_ok=True)
axes=['SP_vs_cortical_raw','cortical_vs_early_raw'];paths=[S/'ALL35_FIXED_PATIENT_SCORES.tsv',S/'ALL35_DIAGNOSTIC_EXPRESSION_UNIT_LEDGER.tsv']
contract={'created_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'question':'Does the observed ZEB1-programme coupling remain after conditioning on native Immature/Cortical/Mature labels?','analysis_role':'Exploratory discrimination after marginal results were already known; not preregistered confirmation','independent_unit':'original35 diagnostic discovery patients','sources':[{'file':str(p),'sha256':hashlib.sha256(p.read_bytes()).hexdigest()} for p in paths],'input_rule':'exact one-to-one expression_unit_id join; retain all35; no clinical D1 mapping','covariate':'original GEO three-category disease state; no merging or assumed ordered clock','estimand':'correlation of global midrank residuals after projection onto native state indicator matrix','test_family':axes,'permutations':20000,'permutation_rule':'permute ZEB1 ranks only within original state; same permutations across axes; two-sided absolute statistic plus1 correction','bootstrap':10000,'bootstrap_rule':'resample paired patients within each of the three original categories; recompute ranks/projection per bootstrap','seed':202610101,'multiplicity':'BH2 across the two prespecified conditional tests','sensitivity':'all35 leave-one-patient-out correlations; no subgroup selection or retesting','scope_limit':'conditional cross-sectional expression association; neither independent OS nor native ZEB1 perturbation nor exact developmental time'}
(O/'EXPLORATORY_CONDITIONAL_STATE_CONTRACT.json').write_text(json.dumps(contract,ensure_ascii=False,indent=2),encoding='utf-8')
a=pd.read_csv(paths[0],sep='\t');b=pd.read_csv(paths[1],sep='\t');assert a['expression_unit_id'].is_unique and b['expression_unit_id'].is_unique
d=a.merge(b,on='expression_unit_id',how='outer',validate='one_to_one',indicator=True);assert len(d)==35 and (d['_merge']=='both').all();assert d[['ZEB1_log2_FPKM_plus1']+axes].notna().all().all()
labels=d['source_disease_state'].to_numpy();levels=['Immature_TALL','Cortical_TALL','Mature_TALL'];groups=[np.flatnonzero(labels==x) for x in levels];assert list(map(len,groups))==[13,17,5]
y=d['ZEB1_log2_FPKM_plus1'].to_numpy();xs=d[axes].to_numpy();rng=np.random.default_rng(contract['seed'])
def residuals(v,g):
 r=rankdata(v,method='average');return r-np.array([r[g==z].mean() for z in g])
def corr(v,w):return float(v@w/np.sqrt(v@v*(w@w)))
yr=residuals(y,labels);xr=np.column_stack([residuals(xs[:,i],labels) for i in range(2)]);observed=np.array([corr(yr,xr[:,i]) for i in range(2)])
null=np.empty((contract['permutations'],2))
for i in range(contract['permutations']):
 perm=np.arange(35)
 for g in groups:perm[g]=rng.permutation(g)
 yp=yr[perm];null[i]=[corr(yp,xr[:,k]) for k in range(2)]
pv=(1+np.sum(np.abs(null)>=np.abs(observed)-1e-14,axis=0))/(len(null)+1);order=np.argsort(pv);adjusted=np.empty(2);adjusted[order]=np.minimum(1,np.minimum.accumulate((pv[order]*2/np.arange(1,3))[::-1])[::-1])
boot=np.empty((contract['bootstrap'],2))
for i in range(contract['bootstrap']):
 idx=np.concatenate([rng.choice(g,size=len(g),replace=True) for g in groups]);g=labels[idx];yy=residuals(y[idx],g)
 for k in range(2):boot[i,k]=corr(yy,residuals(xs[idx,k],g))
loo=[]
for i in range(35):
 keep=np.arange(35)!=i;g=labels[keep];yy=residuals(y[keep],g)
 for k,ax in enumerate(axes):loo.append({'left_out':d['expression_unit_id'].iloc[i],'axis':ax,'rho':corr(yy,residuals(xs[keep,k],g))})
rows=[]
for k,ax in enumerate(axes):
 ci=np.quantile(boot[:,k],[.025,.975]);v=[r['rho'] for r in loo if r['axis']==ax]
 rows.append({'axis':ax,'n':35,'native_state_counts':dict(zip(levels,map(len,groups))),'conditional_rank_rho':float(observed[k]),'two_sided_within_state_permutation_p':float(pv[k]),'BH2_q':float(adjusted[k]),'stratified_patient_bootstrap_CI95_low':float(ci[0]),'stratified_patient_bootstrap_CI95_high':float(ci[1]),'leave_one_out_rho_min':min(v),'leave_one_out_rho_max':max(v),'supported_BH2':bool(adjusted[k]<.05)})
result={'contract':contract,'results':rows,'input_members_unchanged':True,'native_state_time_clock_assumed':False,'new_clinical_models':0,'new_native_perturbation_models':0}
(O/'RESULTS.json').write_text(json.dumps(result,ensure_ascii=False,indent=2),encoding='utf-8');pd.DataFrame(rows).to_csv(O/'ALL2_CONDITIONAL_ASSOCIATION_RESULTS.tsv',sep='\t',index=False);pd.DataFrame(loo).to_csv(O/'ALL70_LEAVE_ONE_PATIENT_OUT_READBACK.tsv',sep='\t',index=False);d.drop(columns='_merge').to_csv(O/'EXACT35_JOINED_INPUT.tsv',sep='\t',index=False);np.savez_compressed(O/'CONDITIONAL_NULL_AND_STRATIFIED_BOOTSTRAP.npz',null=null,bootstrap=boot)
print('AIIMS35_CONDITIONAL_NATIVE_STATE_RESULTS',json.dumps(rows,ensure_ascii=False),flush=True)

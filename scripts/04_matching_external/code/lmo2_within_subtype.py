"""Bounded R04 supplement from existing expression exports; no upstream refit."""
from pathlib import Path
import hashlib
import json
import warnings
import numpy as np
import pandas as pd
from scipy.stats import rankdata
from statsmodels.stats.multitest import multipletests
import statsmodels.formula.api as smf

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / 'data/validation/lmo2_within_subtype'
OUT.mkdir(exist_ok=True)
module = ROOT.parent / 'modules/module3/tables/independent_keygenes.tsv'
target = ROOT / 'data/validation/lmo2_subgroups_existing/TALL_clinical_annotation_full.csv'
SEED, BOOT, PERM, MIN_N = 20260924, 2000, 9999, 10
g = pd.read_csv(module, sep='\t')
g = g[(g.cohort != 'TARGET-ALL-P2') & g.disease.fillna('').str.contains('T-ALL',case=False)].copy()
g = g.rename(columns={'gsm':'sample','level2':'subtype'})
t = pd.read_csv(target).rename(columns={'Sample_ID':'sample','MolecularSubtype':'subtype'})
t['cohort'] = 'TARGET-ALL-P2'
d = pd.concat([g, t],ignore_index=True)[['cohort','sample','subtype','LMO2','ZEB1']]
assert not d.duplicated(['cohort','sample']).any()
d.to_csv(OUT/'sample_inputs.tsv',sep='\t',index=False)

def cor_rows(a,b):
    a=rankdata(a,axis=-1); b=rankdata(b,axis=-1)
    a=a-a.mean(axis=-1,keepdims=True); b=b-b.mean(axis=-1,keepdims=True)
    with np.errstate(invalid='ignore',divide='ignore'):
        return (a*b).sum(axis=-1)/np.sqrt((a*a).sum(axis=-1)*(b*b).sum(axis=-1))

rows=[]; models=[]
for i, ((cohort,subtype),part) in enumerate(d.groupby(['cohort','subtype'],sort=True)):
    part=part.dropna(subset=['LMO2','ZEB1']); n=len(part)
    x=part.LMO2.to_numpy(); y=part.ZEB1.to_numpy()
    variable = n>=3 and len(np.unique(x))>1 and len(np.unique(y))>1
    # Unknown molecular labels remain descriptive, not an interaction category.
    eligible=variable and n>=MIN_N and 'unknown' not in subtype.lower()
    r={'cohort':cohort,'subtype':subtype,'n':n,'rho':float(cor_rows(x,y)) if variable else np.nan,
       'inference':eligible,'ci_low':np.nan,'ci_high':np.nan,'p_perm':np.nan,'valid_boot':0}
    if eligible:
        rng=np.random.default_rng(SEED+i)
        ix=rng.integers(n,size=(BOOT,n)); bs=cor_rows(x[ix],y[ix]); bs=bs[np.isfinite(bs)]
        r['ci_low'],r['ci_high']=np.quantile(bs,[.025,.975]); r['valid_boot']=len(bs)
        xp=np.broadcast_to(x,(PERM,n)); yp=np.vstack([rng.permutation(y) for _ in range(PERM)])
        null=cor_rows(xp,yp)
        r['p_perm']=(1+np.sum(np.abs(null)>=abs(r['rho'])-1e-12))/(PERM+1)
        p=part.copy()
        for gene,label in [('LMO2','xrank'),('ZEB1','yrank')]:
            v=rankdata(p[gene]); p[label]=(v-v.mean())/v.std(ddof=0)
        models.append(p)
    rows.append(r)
res=pd.DataFrame(rows); mask=res.inference
res.loc[mask,'q_bh']=multipletests(res.loc[mask,'p_perm'],method='fdr_bh')[1]
res.to_csv(OUT/'within_subtype.tsv',sep='\t',index=False)
tests=[]
for cohort,part in pd.concat(models).groupby('cohort'):
    k=part.subtype.nunique()
    rec={'cohort':cohort,'n':len(part),'subtypes':k,'p_hc3':np.nan,'status':'fewer_than_two_eligible_subtypes'}
    if k>=2:
        fit=smf.ols('yrank ~ xrank * C(subtype)',part).fit(cov_type='HC3')
        terms=[j for j,name in enumerate(fit.params.index) if ':' in name]
        restriction=np.eye(len(fit.params))[terms]
        w=fit.wald_test(restriction,scalar=True,use_f=False)
        rec.update(p_hc3=float(w.pvalue),wald_chi2=float(w.statistic),df=len(terms),status='exploratory_asymptotic')
    tests.append(rec)
tests=pd.DataFrame(tests); mask=tests.p_hc3.notna()
tests.loc[mask,'q_bh']=multipletests(tests.loc[mask,'p_hc3'],method='fdr_bh')[1]
tests.to_csv(OUT/'cohort_interactions.tsv',sep='\t',index=False)
(OUT/'manifest.json').write_text(json.dumps({'seed':SEED,'bootstrap':BOOT,'permutations':PERM,
    'min_n':MIN_N,'rules':'Original subtype labels; no cross-cohort subtype pooling. Unknown and n<10 descriptive only. Percentile bootstrap; two-sided permutation plus-one P; BH across all eligible strata. HC3 omnibus rank-slope interactions with separate BH across evaluable cohorts; asymptotic exploratory tests.',
    'inputs':[{'path':str(p),'sha256':hashlib.sha256(p.read_bytes()).hexdigest()} for p in [module,target]],
    'n_strata':len(res),'n_inferential':int(res.inference.sum())},indent=2),encoding='utf-8')
print(res.loc[res.inference,['cohort','subtype','n','rho','ci_low','ci_high','p_perm','q_bh']].to_string(index=False))
print(tests.to_string(index=False))

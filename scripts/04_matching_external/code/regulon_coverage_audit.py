"""Re-run the archived CollecTRI network and export measurement coverage and siRNA effects."""
from pathlib import Path
import argparse
import hashlib
import json
import sys
import importlib.metadata
import numpy as np
import pandas as pd
from scipy import stats

p=argparse.ArgumentParser()
p.add_argument('--project',required=True)
p.add_argument('--network',required=True)
p.add_argument('--output',required=True)
a=p.parse_args()
root=Path(a.project); out=Path(a.output); out.mkdir(parents=True,exist_ok=True)
sys.path.insert(0,str(root/'total_fig467F/pydeps'))
import decoupler as dc
paths={'network':Path(a.network),'counts':root/'module9/processed/gse110635_symbol_counts.tsv.gz',
       'samples':root/'module9/processed/gse110635_symbol_coldata.tsv'}
raw=pd.read_csv(paths['network'],sep='\t')
inh=raw.is_inhibition.astype(str).isin(['1','True','true'])
net=pd.DataFrame({'source':raw.source_genesymbol.astype(str),'target':raw.target_genesymbol.astype(str),'weight':np.where(inh,-1.,1.)})
net=net[(net.source!='nan')&(net.target!='nan')].groupby(['source','target'],as_index=False).weight.mean()
net.to_csv(out/'frozen_network.tsv.gz',sep='\t',index=False)
counts=pd.read_csv(paths['counts'],sep='\t',index_col=0)
cd=pd.read_csv(paths['samples'],sep='\t').set_index('gsm',drop=False)
assert counts.index.is_unique and counts.columns.is_unique and cd.index.is_unique
mat=np.log1p(counts.div(counts.sum(axis=0),axis=1)*1e6).T
assert set(mat.index)==set(cd.index)
cd=cd.loc[mat.index]
cd.to_csv(out/'samples.tsv',sep='\t',index=False)
net['in_matrix']=net.target.isin(counts.index)
net['detected_any']=net.target.isin(counts.index[counts.gt(0).any(axis=1)])
net['detected_all']=net.target.isin(counts.index[counts.gt(0).all(axis=1)])
net['nonzero_weight']=net.weight.ne(0)
coverage=net.groupby('source').agg(network_targets=('target','size'),matrix_targets=('in_matrix','sum'),
    detected_any=('detected_any','sum'),detected_all=('detected_all','sum'),nonzero_weight=('nonzero_weight','sum')).reset_index()
acts,_=dc.mt.ulm(mat,net[['source','target','weight']],tmin=5)
acts.to_csv(out/'sample_ulm_activities.tsv',sep='\t')
coverage['activity_returned']=coverage.source.isin(acts.columns)
coverage.to_csv(out/'coverage_all.tsv',sep='\t',index=False)
net[net.source.isin(['TLX1','ZEB1'])].to_csv(out/'tlx1_zeb1_target_coverage.tsv',sep='\t',index=False)
rows=[]
ctrl=acts.loc[cd.group.eq('control')]
assert len(ctrl)==3
for label,mask in [('pooled',cd.group.isin(['siTLX1_1','siTLX1_2'])),('siTLX1_1',cd.group.eq('siTLX1_1')),('siTLX1_2',cd.group.eq('siTLX1_2'))]:
    kd=acts.loc[mask]
    for tf in acts.columns:
        test=stats.ttest_ind(kd[tf],ctrl[tf],equal_var=False)
        rows.append(dict(contrast=label,tf=tf,mean_ctrl=ctrl[tf].mean(),mean_kd=kd[tf].mean(),
                         delta_kd_minus_ctrl=kd[tf].mean()-ctrl[tf].mean(),p=test.pvalue,n_ctrl=len(ctrl),n_kd=len(kd)))
res=pd.DataFrame(rows)
for label,ix in res.groupby('contrast').groups.items():
    res.loc[ix,'fdr']=stats.false_discovery_control(res.loc[ix,'p'].to_numpy())
res.to_csv(out/'activity_contrasts.tsv',sep='\t',index=False)
old=pd.read_csv(root/'total_fig467F/tables/F6_decoupler_ulm_all.tsv',sep='\t')
comp=res[res.contrast.eq('pooled')].merge(old,on='tf',validate='one_to_one',suffixes=('_rerun','_original'))
comp['absolute_delta_difference']=abs(comp.delta_kd_minus_ctrl_rerun-comp.delta_kd_minus_ctrl_original)
comp.to_csv(out/'original_comparison.tsv',sep='\t',index=False)
record={'inputs':{k:{'path':str(v),'sha256':hashlib.sha256(v.read_bytes()).hexdigest()} for k,v in paths.items()},
        'versions':{k:importlib.metadata.version(k) for k in ['decoupler','numpy','pandas','scipy']},
        'network_edges':len(net),'network_TFs':net.source.nunique(),'activity_TFs':acts.shape[1],
        'compared_TFs':len(comp),'maximum_absolute_delta_difference':float(comp.absolute_delta_difference.max()),
        'note':'Archived local OmniPath file reused; no newly downloaded network. siRNAs share controls and are not independent cohorts.'}
(out/'run.json').write_text(json.dumps(record,indent=2),encoding='utf-8')
print(json.dumps(record,indent=2),flush=True)
print(coverage[coverage.source.isin(['TLX1','ZEB1'])].to_string(index=False),flush=True)

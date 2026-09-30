"""Preserve the original MRD models and add protocol as a timing sensitivity."""
from pathlib import Path
import pandas as pd
import numpy as np
from sklearn.metrics import roc_auc_score, brier_score_loss
from clinical_mrd_persistence import loo_predict

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT/'data/validation/clinical_endpoint_audit'
d = pd.read_csv(OUT/'patient_protocol_mrd_manifest.tsv', sep='\t')
d = d[d.persistence_model_eligible].copy()
assert len(d)==66 and set(d.Protocol)=={'T15','T16'}
d['early_log10'] = np.log10(d['Day 15 MRD (%)'])
d['protocol_T15'] = d.Protocol.eq('T15').astype(int)
y = d['Day 42 or 46 MRD (%)'].ge(.01).astype(int).to_numpy()
assert y.sum()==18
models = {
    'mid_MRD': ['early_log10'],
    'mid_MRD_ZEB1': ['early_log10','ZEB1_log'],
    'mid_MRD_protocol': ['early_log10','protocol_T15'],
    'mid_MRD_protocol_ZEB1': ['early_log10','protocol_T15','ZEB1_log']}
rows=[]
for name, columns in models.items():
    p=loo_predict(d[columns].to_numpy(float),y)
    d['loo_'+name]=p
    rows.append(dict(model=name,n_patients=len(d),events=int(y.sum()),
                     auc=roc_auc_score(y,p),brier=brier_score_loss(y,p),
                     predictors=','.join(columns),method='LOO training-fold scaling; logistic L2 C=1',
                     scope='internal exploratory protocol sensitivity; no selected model'))
pd.DataFrame(rows).to_csv(OUT/'protocol_model_sensitivity.tsv',sep='\t',index=False)
d.to_csv(OUT/'protocol_model_predictions.tsv',sep='\t',index=False)
print(pd.DataFrame(rows).to_string(index=False))

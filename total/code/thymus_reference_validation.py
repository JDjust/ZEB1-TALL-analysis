"""Donor-held-out centroid validation and exploratory patient reference affinity."""
from pathlib import Path
import json
import numpy as np
import pandas as pd
from scipy.stats import entropy
from sklearn.metrics import confusion_matrix, balanced_accuracy_score, accuracy_score

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'data/validation/thymus_reference'
STAGES=['DN','DP','CD4SP','CD8SP']


def predict(train,test,genes):
    mean=train[genes].mean();sd=train[genes].std(ddof=1).replace(0,1)
    z=(train[genes]-mean)/sd
    centroid=z.assign(stage=train.stage.to_numpy()).groupby('stage')[genes].mean().loc[STAGES].to_numpy()
    query=((test[genes]-mean)/sd).to_numpy()
    cn=np.linalg.norm(centroid,axis=1,keepdims=True)
    qn=np.linalg.norm(query,axis=1,keepdims=True)
    assert (cn>0).all() and (qn>0).all()
    sim=(query/qn)@(centroid/cn).T
    dist=1-sim
    affinity=np.exp((sim-sim.max(axis=1,keepdims=True))/.2)
    affinity/=affinity.sum(axis=1,keepdims=True)
    result=test[['donor','stage']].copy() if 'donor' in test else test[['patient','n_cells']].copy()
    result['predicted_stage']=[STAGES[i] for i in sim.argmax(axis=1)]
    result['nearest_cosine_distance']=dist.min(axis=1)
    result['similarity_margin']=np.sort(sim,axis=1)[:,-1]-np.sort(sim,axis=1)[:,-2]
    result['affinity_entropy']=entropy(affinity,axis=1)/np.log(len(STAGES))
    for i,stage in enumerate(STAGES):
        result['similarity_'+stage]=sim[:,i]
        result['affinity_'+stage]=affinity[:,i]
    return result


def main():
    settings=json.loads((OUT/'aggregate_run.json').read_text())
    genes=settings['markers']
    assert not set(genes)&set(settings['excluded_genes'])
    original=pd.read_csv(OUT/'normal_donor_stage_markers.tsv',sep='\t')
    ref=original[original.n_cells>=20].copy()
    assert not ref[['donor','stage']].duplicated().any()
    predictions=[]
    for donor in sorted(ref.donor.unique()):
        train=ref[ref.donor.ne(donor)];test=ref[ref.donor.eq(donor)]
        assert set(train.stage)==set(STAGES)
        predictions.append(predict(train,test,genes))
    cv=pd.concat(predictions,ignore_index=True)
    cutoff=float(cv.nearest_cosine_distance.quantile(.95))
    cv.to_csv(OUT/'heldout_donor_predictions.tsv',sep='\t',index=False)
    matrix=confusion_matrix(cv.stage,cv.predicted_stage,labels=STAGES)
    pd.DataFrame(matrix,index=STAGES,columns=STAGES).to_csv(OUT/'heldout_confusion.tsv',sep='\t')
    donor_metrics=cv.assign(correct=cv.stage.eq(cv.predicted_stage)).groupby('donor').agg(
        n_stages=('correct','size'),accuracy=('correct','mean')).reset_index()
    donor_metrics.to_csv(OUT/'heldout_donor_accuracy.tsv',sep='\t',index=False)
    query=pd.read_csv(OUT/'query_patient_markers.tsv',sep='\t')
    result=predict(ref,query,genes)
    result['outside_heldout_95pct_distance']=result.nearest_cosine_distance>cutoff
    result.to_csv(OUT/'query_patient_reference_affinity.tsv',sep='\t',index=False)
    summary=dict(n_donors=int(ref.donor.nunique()),n_donor_stage_profiles=len(ref),
        min_cells=20,n_excluded_profiles=int((original.n_cells<20).sum()),
        n_features=len(genes),excluded_genes=settings['excluded_genes'],
        heldout_profile_accuracy=accuracy_score(cv.stage,cv.predicted_stage),
        heldout_balanced_accuracy=balanced_accuracy_score(cv.stage,cv.predicted_stage),
        donor_mean_accuracy=float(donor_metrics.accuracy.mean()),
        heldout_distance95=cutoff,n_query_patients=len(result),
        n_query_outside_distance95=int(result.outside_heldout_95pct_distance.sum()),
        affinity_temperature=.2,scope='donor-stage mean profiles, not per-cell classification',
        caveat='Softmax affinity is not calibrated probability; distance flag is exploratory, not validated OOD detection. No cross-cohort ZEB1 residual is computed.')
    (OUT/'validation_summary.json').write_text(json.dumps(summary,indent=2))
    print(json.dumps(summary,indent=2));print(result.to_string(index=False))


if __name__=='__main__':main()

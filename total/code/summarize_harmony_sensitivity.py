"""Score the two label sensitivities using frozen original patient scaling."""
from pathlib import Path
import json
import numpy as np
import pandas as pd
from scipy.stats import spearmanr
from three_gene_patient_validation import perm_p, boot_ci

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'data/validation/harmony_label_sensitivity'


def main():
    run=json.loads((OUT/'run.json').read_text())
    assert run['convergence_return'] and run['iterations'] < 50
    d=pd.read_csv(OUT/'patient_label_pseudobulk.tsv',sep='\t').assign(run_id='1736')
    initial=pd.read_csv(OUT/'initial_run_1735/patient_label_pseudobulk.tsv',sep='\t').assign(run_id='1735')
    d=pd.concat([initial,d],ignore_index=True)
    base=pd.read_csv(ROOT/'data/validation/three_gene_patient/patient_three_gene_score.tsv',sep='\t')
    base=base[base.cohort.eq('GSE227122')].set_index('patient')
    genes=['GATA3','TCF7','BCL11B']
    means,sd=base[genes].mean(),base[genes].std(ddof=1)
    d['three_gene_score']=((d[genes]-means)/sd).mean(axis=1)
    d['jaccard']=d.n_intersection/d.n_union
    rows=[]
    for (run_id,resolution),g in d.groupby(['run_id','resolution']):
        g=g[g.n_new>=30].copy()
        rng=np.random.default_rng(20260924+int(resolution*100))
        x,y=g.ZEB1.to_numpy(),g.three_gene_score.to_numpy()
        lo,hi=boot_ci(x,y,2000,rng)
        rows.append(dict(run_id=run_id,resolution=resolution,n_patients=len(g),
            rho=spearmanr(x,y).statistic,permutation_p=perm_p(x,y,10000,rng),
            bootstrap_ci95_lo=lo,bootstrap_ci95_hi=hi,
            score_scaling='original GSE227122 patient means/SDs frozen',
            min_patient_jaccard=g.jaccard.min(),max_patient_jaccard=g.jaccard.max(),
            interpretation='same heuristic under altered clustering; not independent malignant ground truth'))
    d.to_csv(OUT/'patient_scores_and_overlap.tsv',sep='\t',index=False)
    pd.DataFrame(rows).to_csv(OUT/'score_association_sensitivity.tsv',sep='\t',index=False)
    print(pd.DataFrame(rows).to_string(index=False))


if __name__=='__main__':main()

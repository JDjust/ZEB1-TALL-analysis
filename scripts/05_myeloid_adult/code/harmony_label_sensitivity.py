"""Recreate original Omicverse PCA/Harmony path with an extended iteration cap.

Keep original object immutable. This checks heuristic-label stability, not
malignancy ground truth. Parameters are never chosen using ZEB1 association.
"""
from pathlib import Path
import argparse
import hashlib
import inspect
import json
import importlib.metadata
import numpy as np
import pandas as pd
import anndata as ad
import scanpy as sc
import omicverse as ov
from scipy.stats import entropy, spearmanr
from sklearn.metrics import adjusted_rand_score
from omicverse.external.harmony import run_harmony


def label(obs, clusters):
    tmp = obs.copy()
    tmp['cluster'] = np.asarray(clusters).astype(str)
    details = []
    for name, g in tmp.groupby('cluster', observed=True):
        e = entropy(g.patient.value_counts(normalize=True), base=2) if len(g)>=5 else np.nan
        t, my, b = g.score_T.median(), g.score_Myeloid.median(), g.score_B.median()
        selected = t > my and t > b and ((np.isfinite(e) and e < 1.6) or t > .4)
        details.append(dict(cluster=name, n_cells=len(g), patient_entropy=e,
                            median_T=t, median_myeloid=my, median_B=b, selected=selected))
    selected = [v['cluster'] for v in details if v['selected']]
    return (tmp.cluster.isin(selected) & tmp.lineage.isin(['T','Unknown'])).to_numpy(), details


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--input',type=Path,required=True)
    ap.add_argument('--output',type=Path,required=True);args=ap.parse_args()
    args.output.mkdir(parents=True,exist_ok=False)
    a=ad.read_h5ad(args.input)
    old=a.obs.malignant.eq('malignant').to_numpy()
    reproduced, original_details=label(a.obs,a.obs.leiden)
    assert np.array_equal(reproduced,old), 'Original heuristic could not be reproduced'
    print('Original malignant heuristic reproduced exactly',flush=True)
    # batch_correction default use_rep is scaled|original|X_pca; the original
    # wrapper calculated this on a temporary copy. Recreate that path explicitly.
    if 'scaled|original|X_pca' not in a.obsm:
        ov.pp.scale(a)
        ov.pp.pca(a,layer='scaled',n_pcs=30)
    representation=a.obsm['scaled|original|X_pca']
    print('Recreated default representation',representation.shape,flush=True)
    h=run_harmony(representation,a.obs,'sample',max_iter_harmony=50,random_state=0,device='cpu')
    a.obsm['X_harmony_extended']=h.result()
    objectives=[float(v) for v in h.objective_harmony]
    pd.DataFrame({'iteration':range(len(objectives)),'objective':objectives}).to_csv(
        args.output/'harmony_objectives.tsv',sep='\t',index=False)
    sc.pp.neighbors(a,use_rep='X_harmony_extended',n_neighbors=15,n_pcs=30,random_state=0)
    result_rows=[];patient_rows=[];cell=a.obs[['sample','patient','timepoint','malignant','leiden']].copy()
    gene_idx=[a.var_names.get_loc(g) for g in ['ZEB1','GATA3','TCF7','BCL11B']]
    for resolution in [.6,1.0]:
        key=f'extended_{resolution}'
        sc.tl.leiden(a,resolution=resolution,key_added=key,flavor='igraph',n_iterations=2,
                     directed=False,random_state=0)
        new,detail=label(a.obs,a.obs[key]);cell[key]=a.obs[key].astype(str)
        cell[key+'_malignant']=new
        pd.DataFrame(detail).to_csv(args.output/f'clusters_{resolution}.tsv',sep='\t',index=False)
        dx=a.obs.timepoint.eq('Dx').to_numpy()
        for patient in sorted(a.obs.loc[dx,'patient'].unique()):
            pm=dx & a.obs.patient.eq(patient).to_numpy()
            selected=pm & new
            raw=np.asarray(a.layers['counts'][selected,:].sum(axis=0)).ravel()
            lib=raw.sum()
            row=dict(resolution=resolution,patient=patient,n_old=int((pm&old).sum()),
                     n_new=int(selected.sum()),n_intersection=int((pm&old&new).sum()),
                     n_union=int((pm&(old|new)).sum()),library_umi=float(lib))
            for gene,idx in zip(['ZEB1','GATA3','TCF7','BCL11B'],gene_idx):
                row[gene]=float(np.log2(1+raw[idx]*1e6/lib)) if lib>0 else np.nan
            patient_rows.append(row)
        result_rows.append(dict(resolution=resolution,ari_vs_original=adjusted_rand_score(a.obs.leiden,a.obs[key]),
            n_original=int(old.sum()),n_extended=int(new.sum()),n_label_changed=int((old!=new).sum()),
            malignant_jaccard=float((old&new).sum()/(old|new).sum())))
    cell.to_csv(args.output/'cell_label_comparison.tsv.gz',sep='\t',index=True)
    pd.DataFrame(patient_rows).to_csv(args.output/'patient_label_pseudobulk.tsv',sep='\t',index=False)
    pd.DataFrame(result_rows).to_csv(args.output/'label_stability_summary.tsv',sep='\t',index=False)
    run=dict(input=str(args.input),n_cells=a.n_obs,iterations=len(objectives)-1,
        convergence_return=bool(h._check_convergence(1) if hasattr(h, '_check_convergence')
                                else h.check_convergence(1)),objective_history=objectives,
        representation='recreated original wrapper default scaled|original|X_pca',
        representation_sha256=hashlib.sha256(np.ascontiguousarray(representation).tobytes()).hexdigest(),
        harmony_source_sha256=hashlib.sha256(Path(inspect.getfile(run_harmony)).read_bytes()).hexdigest(),
        packages={p:importlib.metadata.version(p) for p in ['omicverse','scanpy','anndata','numpy','torch']},
        scope='extended-cap sensitivity with current installed code; no claim of ground-truth malignant annotation')
    (args.output/'run.json').write_text(json.dumps(run,indent=2))
    print(json.dumps(run,indent=2),flush=True);print(pd.DataFrame(result_rows).to_string(index=False),flush=True)


if __name__=='__main__':main()

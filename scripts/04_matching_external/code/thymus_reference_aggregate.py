"""Export donor-stage normal and patient malignant marker means from raw counts."""
import argparse
import json
from pathlib import Path
import anndata as ad
import numpy as np
import pandas as pd

MARKERS=['CD34','KIT','HES1','RAG1','RAG2','DNTT','PTCRA','CD1A','CD1B',
         'CD4','CD8A','CD8B','CCR7','CD3D','CD3E','CD3G','TRAC','LCK']
EXCLUDED=['ZEB1','ZEB2','LMO2','GATA3','TCF7','BCL11B','IL7R']
STAGES={'double negative thymocyte':'DN','double-positive, alpha-beta thymocyte':'DP',
        'CD4-positive, alpha-beta T cell':'CD4SP','CD8-positive, alpha-beta T cell':'CD8SP'}


def gene_indices(var):
    names=var.feature_name.astype(str) if 'feature_name' in var else pd.Series(var.index,index=var.index)
    out=[]
    for gene in MARKERS:
        hit=np.flatnonzero(names.to_numpy()==gene)
        if len(hit)!=1:raise ValueError(f'{gene}: {len(hit)} matches')
        out.append(int(hit[0]))
    return out


def normalized_markers(raw, idx, n):
    out=np.empty((n,len(idx)),dtype=np.float32)
    for start in range(0,n,4096):
        x=raw[start:start+4096,:].tocsr()
        total=np.asarray(x.sum(axis=1)).ravel()
        if (total<=0).any():raise ValueError('empty raw-count cell')
        vals=x[:,idx].toarray()
        if not np.allclose(vals,np.round(vals)):raise ValueError('noninteger marker counts')
        out[start:start+len(total)]=np.log1p(vals/total[:,None]*10000)
    return out


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--reference',type=Path,required=True)
    ap.add_argument('--query',type=Path,required=True);ap.add_argument('--output',type=Path,required=True)
    args=ap.parse_args();args.output.mkdir(parents=True,exist_ok=False)
    a=ad.read_h5ad(args.reference,backed='r')
    assert a.raw is not None
    x=normalized_markers(a.raw.X,gene_indices(a.raw.var),a.n_obs)
    data=pd.DataFrame(x,columns=MARKERS,index=a.obs_names)
    data['donor']=a.obs.donor_id.astype(str)
    data['stage']=a.obs.cell_type.astype(str).map(STAGES)
    data['age_label']=a.obs.development_stage.astype(str)
    data=data[data.stage.notna()]
    groups=data.groupby(['donor','stage'],observed=True)
    agg=groups[MARKERS].mean().join(groups.size().rename('n_cells')).reset_index()
    agg=agg.merge(groups.age_label.first().reset_index(),on=['donor','stage'],validate='one_to_one')
    agg.to_csv(args.output/'normal_donor_stage_markers.tsv',sep='\t',index=False)
    citation=str(a.uns.get('citation',''));a.file.close()
    q=ad.read_h5ad(args.query)
    keep=(q.obs.timepoint.eq('Dx') & q.obs.malignant.eq('malignant')).to_numpy()
    counts=q.layers['counts'][keep,:]
    vals=normalized_markers(counts,gene_indices(q.var),counts.shape[0])
    frame=pd.DataFrame(vals,columns=MARKERS)
    frame['patient']=q.obs.loc[keep,'patient'].astype(str).to_numpy()
    groups=frame.groupby('patient')
    query=groups[MARKERS].mean().join(groups.size().rename('n_cells')).reset_index()
    query.to_csv(args.output/'query_patient_markers.tsv',sep='\t',index=False)
    (args.output/'aggregate_run.json').write_text(json.dumps(dict(reference=str(args.reference),
        query=str(args.query),reference_citation=citation,markers=MARKERS,excluded_genes=EXCLUDED,
        normalization='per-cell log1p(10000 * raw UMI / retained raw library total), then cell mean',
        query_label='original project heuristic diagnosis malignant',
        query_caveat='saved count layer excludes a small fraction of genes; normalization does not remove cohort batch effects'),indent=2))
    print('Normal donor-stage groups',len(agg),'donors',agg.donor.nunique(),'query patients',len(query),flush=True)


if __name__=='__main__':main()

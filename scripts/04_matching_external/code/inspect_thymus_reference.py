from pathlib import Path
import anndata as ad
import numpy as np
import json

p=Path('/data-b/liangfuhua/projects/TALL_dataset_download/data/external_fig467F/Park_HTA/Human_T_cells.h5ad')
a=ad.read_h5ad(p,backed='r')
print('shape',a.shape,'obs',list(a.obs),'var',list(a.var),'layers',list(a.layers),flush=True)
print('uns',str(a.uns)[:6500],flush=True)
print('raw',a.raw.shape if a.raw is not None else None,flush=True)
print(a.var.head().to_string(),flush=True)
for name in ['donor_id','cell_type','tissue','development_stage']:
    if name in a.obs:print(name,a.obs[name].value_counts().to_string(),flush=True)
x=a.X[:3,:]
print('X first rows sum/max',np.asarray(x.sum(axis=1)).ravel(),x.max(),flush=True)
if a.raw is not None:
    x=a.raw.X[:3,:]
    print('raw first rows sum/max',np.asarray(x.sum(axis=1)).ravel(),x.max(),flush=True)

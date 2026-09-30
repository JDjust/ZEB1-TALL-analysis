"""Repair S4D's post-freeze t-statistic pairing from existing saved outputs."""
from pathlib import Path
from hashlib import sha256
import json

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

ROOT=Path(__file__).resolve().parents[1]
DATA=ROOT/'data/deepen_2026/gate2_program'
STAGE=DATA
STAGE.mkdir(parents=True,exist_ok=True)
def sha(path):return sha256(path.read_bytes()).hexdigest()

primary_path=DATA/'all_genes_residual_effects.tsv'
legacy_path=DATA/'residual_t_primary_vs_pc15.tsv'
frozen_path=DATA/'frozen_ZEB_side50.tsv'
before={p.name:sha(p) for p in (primary_path,legacy_path,frozen_path)}
primary=pd.read_csv(primary_path,sep='\t',dtype={'gene_id':str,'symbol':str})
legacy=pd.read_csv(legacy_path,sep='\t',dtype={'symbol':str})
exclude={'ZEB1','ZEB2','CD34','LYL1','CD1A'}
primary=primary.loc[~primary.symbol.str.upper().isin(exclude),['gene_id','symbol','t']].copy()
assert len(primary)==len(legacy)==24356
assert primary.gene_id.notna().all() and primary.gene_id.is_unique
assert primary.symbol.notna().all() and primary.symbol.is_unique
assert legacy.symbol.notna().all() and legacy.symbol.is_unique
assert legacy.t_pc.notna().all() and primary.t.notna().all()
assert set(primary.symbol)==set(legacy.symbol)

# The legacy sensitivity output lacks gene_id. Use its unique symbol only as a
# checked fallback to recover gene_id, then retain gene_id as the canonical key.
joined=legacy.merge(primary,on='symbol',how='left',validate='one_to_one',indicator=True)
assert len(joined)==24356 and joined._merge.eq('both').all()
assert joined.gene_id.is_unique and joined.symbol.is_unique
joined=joined.rename(columns={'t_primary':'t_primary_old_positional','t':'t_primary'})
joined=joined[['gene_id','symbol','t_primary','t_pc','t_primary_old_positional']]
assert np.isfinite(joined[['t_primary','t_pc','t_primary_old_positional']].to_numpy()).all()

old_match=int(np.isclose(joined.t_primary_old_positional,joined.t_primary,rtol=1e-9,atol=1e-9).sum())
old_rho=float(spearmanr(joined.t_primary_old_positional,joined.t_pc).statistic)
new_rho=float(spearmanr(joined.t_primary,joined.t_pc).statistic)
assert old_match==4
assert abs(old_rho-0.009936501542298415)<1e-10
assert abs(new_rho-0.4827668014127643)<1e-10

fixed=STAGE/'residual_t_primary_vs_pc15_KEYED.tsv'
joined.to_csv(fixed,sep='\t',index=False,float_format='%.15g')
after={p.name:sha(p) for p in (primary_path,legacy_path,frozen_path)}
assert before==after
evidence={'input_sha256':before,'output_sha256':sha(fixed),
          'key_policy':'primary gene_id retained; legacy unique symbol validated as fallback',
          'rows':len(joined),'gene_id_unique':True,'symbol_unique':True,'unmatched':0,
          'old_primary_t_matching_rows':old_match,'old_spearman_rho':old_rho,
          'fixed_spearman_rho':new_rho,'t_pc_provenance':'tt2$t from PC1-5 limma model, paired to symbol before original faulty merge; source code inspected, no independent re-fit',
          'primary_t_provenance':'all_genes_residual_effects.tsv t column, gene_id + symbol, produced before freeze',
          'frozen_gene_membership_changed':False,'primary_model_refit':False}
(STAGE/'S4D_KEY_REPAIR_EVIDENCE.json').write_text(json.dumps(evidence,ensure_ascii=False,indent=2),encoding='utf-8')
print(json.dumps({k:evidence[k] for k in ['rows','old_primary_t_matching_rows','old_spearman_rho','fixed_spearman_rho','unmatched']},ensure_ascii=False))

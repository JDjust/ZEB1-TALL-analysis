"""Thin adapter of accepted probe_NOPHO108_fixed_maturation_replication.py.
Rank scoring is the accepted project definition; inference uses SciPy/statsmodels.
"""
from pathlib import Path
import json, hashlib, sys, platform
import pandas as pd
import numpy as np
import scipy
from scipy import stats
import statsmodels
import statsmodels.api as sm
from statsmodels.stats.multitest import multipletests
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
P=Path('D:/_bioinformation/ZEB1');O=Path(__file__).resolve().parent
H=P/'tasks/ZEB1_MATURATION_CLINICAL_STRENGTHEN_20261009';L=P/'tasks/ZEB1_PUBLIC_MAINLINE_NO_SKILL_20261009'
members_path=H/'multigene_maturity_proxy/FIXED200_MODULE_MEMBERS.tsv'
normal_path=P/'tasks/ZEB1_TARGET_DISCOVERY/inputs/NORMAL_DONOR_STATE_COUNTS.tsv.gz'
members=pd.read_csv(members_path,sep='\t');normal=set(pd.read_csv(normal_path,sep='\t',usecols=[0]).iloc[:,0].astype(str))
assert len(members)==200 and members.symbol.nunique()==166 and not members.symbol.isin(['ZEB1','ZEB2']).any()
assert members.groupby(['axis','direction']).size().eq(50).all()
inputs=[members_path,normal_path,O/'EXECUTION_CONTRACT.json',H/'probe_NOPHO108_fixed_maturation_replication.py']
results=[];banks=[];diagnostics=[];coverage=[];cohort_patients={};unit_ledgers=[]
rng=np.random.default_rng(20261010)
for cohort in ['discovery','validation']:
    cp=P/f'tasks/ZEB1_TARGET_DISCOVERY/inputs/Lim_{cohort}_pseudobulk_counts.tsv.gz'
    mp=L/f'LIM_{cohort.upper()}_D0_METADATA.tsv';inputs.extend([cp,mp])
    full=pd.read_csv(cp,sep='\t',index_col=0);md=pd.read_csv(mp,sep='\t').set_index('pseudobulk_id')
    assert md.patient.is_unique and md.timepoint.eq('D0').all() and md.tissue.eq('BM').all() and md.mapped_state.eq('ALL').all()
    unit_ledger=md.reset_index().copy();unit_ledger['cohort_analysis']=cohort;unit_ledger['eligibility']=np.where(unit_ledger.n_cells.ge(100),'ELIGIBLE','UNESTIMATED_LT100');unit_ledgers.append(unit_ledger)
    md=md.loc[md.n_cells.ge(100)].copy()
    assert md.n_cells.ge(100).all() and len(md)>=15 and full.index.is_unique and full.columns.is_unique
    cohort_patients[cohort]=md.patient.tolist()
    d=full.loc[:,md.index];assert np.isfinite(d.to_numpy()).all() and (d.to_numpy()>=0).all() and np.allclose(d,np.round(d))
    common=sorted(set(d.index)&normal);assert len(common)>=20000 and 'ZEB1' in common
    ranks=d.loc[common].rank(axis=0,method='average',pct=True)
    # Library totals are accepted source metadata; log discrepancies rather than silently replacing.
    lib=md.lib_size;assert (lib>0).all()
    scores=md.copy();scores['cohort_analysis']=cohort;scores['ZEB1_log2_CPM_plus1']=np.log2(1+d.loc['ZEB1']/lib*1e6)
    scores['sum_retained_counts']=d.sum(axis=0);scores['library_gap']=lib-d.sum(axis=0)
    scores['log2_mean_umi_plus1']=np.log2(1+scores.mean_umi)
    scores['zero_fraction_shared']=d.loc[common].eq(0).mean(axis=0)
    for axis,g in members.groupby('axis'):
        up=g.loc[g.direction.eq('up') & g.symbol.isin(common),'symbol'];down=g.loc[g.direction.eq('down') & g.symbol.isin(common),'symbol']
        assert len(up)>=45 and len(down)>=45 and len(up)+len(down)>=90
        coverage.append({'cohort':cohort,'axis':axis,'n_up':len(up),'n_down':len(down),'n_ranking_symbols':len(common)})
        y=(ranks.loc[up].mean()-ranks.loc[down].mean()).reindex(md.index).to_numpy();scores[axis]=y
        x=scores.ZEB1_log2_CPM_plus1.to_numpy();rho,p=stats.spearmanr(x,y)
        boot=[]
        for _ in range(10000):
            ix=rng.integers(0,len(x),len(x));rr=stats.spearmanr(x[ix],y[ix]).statistic
            if np.isfinite(rr):boot.append(rr)
        ci=np.quantile(boot,[.025,.975]);loo=[stats.spearmanr(np.delete(x,i),np.delete(y,i)).statistic for i in range(len(x))]
        results.append({'cohort':cohort,'axis':axis,'model':'primary_spearman','n_patients':len(x),'effect':float(rho),'ci_low':float(ci[0]),'ci_high':float(ci[1]),'p':float(p),'loo_min':float(min(loo)),'loo_max':float(max(loo)),'n_valid_bootstrap':len(boot),'status':'ESTIMATED'})
        def zr(v):
            v=stats.rankdata(v,method='average');return (v-v.mean())/v.std(ddof=1)
        X=sm.add_constant(pd.DataFrame({'ZEB1_rank_z':zr(x),'cycle_rank_z':zr(scores.cycle_fraction),'depth_rank_z':zr(scores.log2_mean_umi_plus1)}))
        assert np.isfinite(X.to_numpy()).all() and np.linalg.matrix_rank(X)==X.shape[1] and len(x)-X.shape[1]>=10
        fit=sm.OLS(zr(y),X).fit(cov_type='HC3',use_t=True);fc=fit.conf_int().loc['ZEB1_rank_z']
        adjusted_loo=[float(sm.OLS(np.delete(zr(y),i),X.drop(i)).fit().params.ZEB1_rank_z) for i in range(len(x))]
        results.append({'cohort':cohort,'axis':axis,'model':'cycle_depth_rank_HC3','n_patients':len(x),'effect':float(fit.params.ZEB1_rank_z),'ci_low':float(fc.iloc[0]),'ci_high':float(fc.iloc[1]),'p':float(fit.pvalues.ZEB1_rank_z),'loo_min':min(adjusted_loo),'loo_max':max(adjusted_loo),'status':'ESTIMATED'})
        for term in X.columns:
            diagnostics.append({'cohort':cohort,'axis':axis,'term':term,'beta':float(fit.params[term]),'se_HC3':float(fit.bse[term]),'residual_df':float(fit.df_resid),'condition_number':float(np.linalg.cond(X)),'max_leverage':float(fit.get_influence().hat_matrix_diag.max()),'residual_variance':float(np.var(fit.resid,ddof=1))})
    banks.append(scores.reset_index())
assert not set(cohort_patients['discovery']) & set(cohort_patients['validation'])
res=pd.DataFrame(results);res['BH8_joint_exploratory']=multipletests(res.p,method='fdr_bh')[1]
primary=res.model.eq('primary_spearman');res.loc[primary,'BH4_primary']=multipletests(res.loc[primary,'p'],method='fdr_bh')[1]
res.to_csv(O/'LIM_FIXED_DEVELOPMENT_RESULTS.tsv',sep='\t',index=False)
bank=pd.concat(banks,ignore_index=True);bank.to_csv(O/'LIM_ELIGIBLE49_PATIENT_SCORES.tsv',sep='\t',index=False)
pd.concat(unit_ledgers,ignore_index=True).to_csv(O/'LIM_ALL50_INPUT_UNIT_LEDGER.tsv',sep='\t',index=False)
pd.DataFrame(diagnostics).to_csv(O/'LIM_MODEL_DIAGNOSTICS.tsv',sep='\t',index=False);pd.DataFrame(coverage).to_csv(O/'LIM_COVERAGE.tsv',sep='\t',index=False)
fig,axs=plt.subplots(2,2,figsize=(8,6.5),layout='constrained')
axes=['cortical_vs_early','SP_vs_cortical']
for i,cohort in enumerate(['discovery','validation']):
    b=bank[bank.cohort_analysis.eq(cohort)]
    for j,axis in enumerate(axes):
        ax=axs[i,j];rr=res[(res.cohort==cohort)&(res.axis==axis)&primary].iloc[0]
        ax.scatter(b[axis],b.ZEB1_log2_CPM_plus1,s=28,alpha=.85,color=['#c27836','#417b78'][j],edgecolor='white',linewidth=.35)
        ax.set(xlabel=axis.replace('_',' '),ylabel='ZEB1 log2(CPM+1)',title=f'{cohort.capitalize()} n={len(b)}; rho={rr.effect:.2f}; joint BH8 q={rr.BH8_joint_exploratory:.3g}')
        ax.spines[['top','right']].set_visible(False)
fig.savefig(O/'LIM_PATIENT_DEVELOPMENT_NATIVE.png',dpi=180);fig.savefig(O/'LIM_PATIENT_DEVELOPMENT_NATIVE.pdf');plt.close(fig)
out={'execution_status':'COMPLETE','scientific_status':'INTERPRET_RESULTS_WITH_BOUNDARIES','n_patients':cohort_patients,'models':json.loads(res.to_json(orient='records')),'same_study_two_disjoint_patient_sets':True,'library_gap_range':[float(bank.library_gap.min()),float(bank.library_gap.max())],'zero_fraction_range':[float(bank.zero_fraction_shared.min()),float(bank.zero_fraction_shared.max())],'source_sha256':{str(p.relative_to(P)):hashlib.sha256(p.read_bytes()).hexdigest() for p in inputs},'environment':{'executable':sys.executable,'python':platform.python_version(),'pandas':pd.__version__,'numpy':np.__version__,'scipy':scipy.__version__,'statsmodels':statsmodels.__version__,'matplotlib':matplotlib.__version__},'implementation':'Adapter of existing project fixed-rank score protocol; scipy.stats.spearmanr and statsmodels.OLS HC3; statsmodels multipletests; no custom mature algorithm replacement','boundaries':['Author malignant assignment reused, not new independent genotype identity certification','Response-enriched diagnosis BM patient subsets, not population-representative T-ALL samples','Sparsity and different normal/bulk/single-cell capture; no causal regulation inference','No survival or new clinical biomarker model fitted','Primary and sensitivity are different estimands; HC3 standardized rank slope is not partial Spearman rho']}
(O/'LIM_EXECUTION_RECEIPT.json').write_text(json.dumps(out,indent=2),encoding='utf-8');print(res.to_string(index=False));print('LIBRARY_GAP',out['library_gap_range'],'ZERO_FRACTION',out['zero_fraction_range'])

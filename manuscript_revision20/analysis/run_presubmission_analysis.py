#!/usr/bin/env python3
"""Authorized manuscript reanalysis. Frozen inputs and original result files are read only.

Maintained scipy/statsmodels implement rank, resampling, OLS, HC3 and Wald inference.
The adapter preserves accepted score definitions and exposes rank projections/formulas.
"""
from __future__ import annotations
import argparse, datetime, hashlib, json, math, platform, shutil, warnings
from pathlib import Path
import numpy as np
import pandas as pd
import scipy
from scipy import stats
import statsmodels
import statsmodels.api as sm
from statsmodels.stats.multitest import multipletests

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
AXES = ['cortical_vs_early', 'SP_vs_cortical']
SEED = 20261012

def sha(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as f:
        for block in iter(lambda:f.read(1024*1024),b''): h.update(block)
    return h.hexdigest()

def write_json(path,obj):
    Path(path).write_text(json.dumps(obj,indent=2,ensure_ascii=False,allow_nan=False)+'\n',encoding='utf-8')

def clean(obj):
    if isinstance(obj,dict): return {k:clean(v) for k,v in obj.items()}
    if isinstance(obj,list): return [clean(v) for v in obj]
    if isinstance(obj,(np.floating,float)): return float(obj) if np.isfinite(obj) else None
    if isinstance(obj,np.integer): return int(obj)
    if isinstance(obj,np.bool_): return bool(obj)
    return obj

def save(df,name):
    df.to_csv(HERE/'results'/name,sep='\t',index=False,float_format='%.12g')

def zrank(a):
    a=stats.rankdata(np.asarray(a),method='average')
    return (a-a.mean())/a.std(ddof=0)

def residual(a,group):
    return np.asarray(a)-pd.Series(a).groupby(np.asarray(group),sort=False).transform('mean').to_numpy()

def cor(a,b):
    a=np.asarray(a)-np.mean(a); b=np.asarray(b)-np.mean(b)
    den=np.sqrt(a@a*(b@b))
    return float(a@b/den) if den>0 else np.nan

def rankstat(x,y,axis=-1):
    x=stats.rankdata(x,axis=axis); y=stats.rankdata(y,axis=axis)
    x=x-x.mean(axis=axis,keepdims=True); y=y-y.mean(axis=axis,keepdims=True)
    den=np.sqrt(np.sum(x*x,axis=axis)*np.sum(y*y,axis=axis))
    with np.errstate(invalid='ignore',divide='ignore'): return np.sum(x*y,axis=axis)/den

def bootstrap_rho(x,y,seed):
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        r=stats.bootstrap((np.asarray(x),np.asarray(y)),rankstat,paired=True,vectorized=True,
                          n_resamples=10000,batch=100,method='percentile',rng=np.random.default_rng(seed))
    d=r.bootstrap_distribution; finite=d[np.isfinite(d)]
    lo,hi=np.quantile(finite,[.025,.975]) if len(finite) else [np.nan,np.nan]
    return float(lo),float(hi),len(finite),d

def local_permutation(x,y,seed):
    xr=zrank(x); yr=zrank(y)
    def st(yp,axis=-1): return np.sum(yp*xr,axis=axis)/len(xr)
    n=len(x); nr=math.factorial(n) if n<=8 else 19999
    r=stats.permutation_test((yr,),st,permutation_type='pairings',vectorized=True,
                             n_resamples=nr,batch=256,alternative='two-sided',rng=np.random.default_rng(seed))
    return float(r.pvalue),'exact' if n<=8 else 'Monte_Carlo_19999'

def fit_rank(x,y,group=None):
    X=pd.DataFrame({'intercept':np.ones(len(x)), 'score_rank_sd':x})
    if group is not None:
        d=pd.get_dummies(pd.Series(group),prefix='subtype',drop_first=True,dtype=float)
        X=pd.concat([X,d],axis=1)
    f=sm.OLS(np.asarray(y),X).fit(cov_type='HC3',use_t=True)
    j=list(X.columns).index('score_rank_sd'); ci=f.conf_int().iloc[j]
    return {'effect':float(f.params.iloc[j]),'ci_low':float(ci.iloc[0]),'ci_high':float(ci.iloc[1]),
            'p':float(f.pvalues.iloc[j]),'se_HC3':float(f.bse.iloc[j]),'df_resid':float(f.df_resid),
            'n_patients':int(f.nobs),'r_squared':float(f.rsquared),
            'max_leverage':float(f.get_influence().hat_matrix_diag.max()),
            'condition_number':float(f.condition_number)},f,X

def freeze_inputs(project):
    old=project/'03_publication_outputs/17_figure6_revision_20261011'
    sf=json.loads((old/'SOURCE_FREEZE.json').read_text(encoding='utf-8'))
    records=[]
    names=['bulk','maturity','fixed_members','normal_heldout','normal_state_scores','lim_scores',
           'lim_results','lim_units','lim_modules','lim_coverage','aiims_native','aiims_conditional',
           'nopho','aiims','india','nopho_results','aiims_results','india_results']
    for r in sf:
        if r['alias'] not in names: continue
        p=old/'source_data'/r['frozen']; actual=sha(p); assert actual==r['sha256'],r['alias']
        dest=HERE/'source_inputs'/p.name;shutil.copy2(p,dest)
        records.append({'source':str(p),'snapshot':str(dest.relative_to(HERE)),'sha256':actual,'source_freeze_match':True})
    m=project/'tasks/ZEB1_MATURATION_CLINICAL_STRENGTHEN_20261009/multigene_maturity_proxy'
    q=json.loads((m/'QUALIFICATION.json').read_text(encoding='utf-8'))
    src=project/'tasks/ZEB1_TARGET_DISCOVERY/inputs'
    for name in ['NORMAL_DONOR_STATE_COUNTS.tsv.gz','NORMAL_DONOR_STATE_METADATA.tsv']:
        p=src/name; actual=sha(p)
        expected=next(v for k,v in q['source_sha256'].items() if Path(k).name==name)
        assert actual==expected,name
        shutil.copy2(p,HERE/'source_inputs'/name)
        records.append({'source':str(p),'snapshot':'source_inputs/'+name,'sha256':actual,'accepted_qualification_match':True})
    for name in ['GENE_MAPPING_LEDGER.tsv','cortical_vs_early_NORMAL_SELECTION_LEDGER.tsv','SP_vs_cortical_NORMAL_SELECTION_LEDGER.tsv','QUALIFICATION.json']:
        p=m/name;shutil.copy2(p,HERE/'source_inputs'/name)
        records.append({'source':str(p),'snapshot':'source_inputs/'+name,'sha256':sha(p)})
    for name in ['09_MULTIGENE_MATURITY_PROXY_CONTRACT.json','24_NORMAL_AGGREGATION_CERTIFICATE.json','25_NORMAL_EXPORT_IDENTITY_QA.json']:
        p=old/'method_records'/name;shutil.copy2(p,HERE/'source_inputs'/name)
        records.append({'source':str(p),'snapshot':'source_inputs/'+name,'sha256':sha(p)})
    write_json(HERE/'SOURCE_INPUT_READBACK.json',records)
    return records

def verify_snapshots():
    records=json.loads((HERE/'SOURCE_INPUT_READBACK.json').read_text(encoding='utf-8'))
    for r in records: assert sha(HERE/r['snapshot'])==r['sha256'],r['snapshot']
    rank_receipt=HERE/'BULK_MEMBER_RANK_EXTRACTION_RECEIPT.json'
    if rank_receipt.exists():
        r=json.loads(rank_receipt.read_text(encoding='utf-8'))
        assert sha(HERE/r['snapshot'])==r['snapshot_sha256'],'member rank snapshot'
    return {'status':'PASS','n_snapshots':len(records),'member_rank_snapshot_verified':rank_receipt.exists()}

def bulk_stage():
    s=HERE/'source_inputs'; b=pd.read_csv(s/'bulk.tsv',sep='\t'); m=pd.read_csv(s/'maturity.tsv',sep='\t')
    assert len(b)==1309 and b.sample_id.is_unique and b.subject_sjuid.is_unique
    d=b[['sample_id','subject_sjuid','subtype','ZEB1']].merge(m,on='sample_id',validate='one_to_one',how='outer')
    assert len(d)==1309 and d[['ZEB1']+AXES+['subtype']].notna().all().all() and d.subtype.nunique()==17
    d['ZEB1_rank_sd']=zrank(d.ZEB1);d['ZEB1_subtype_residual']=residual(d.ZEB1_rank_sd,d.subtype)
    rows=[]; within=[]; interactions=[];decomp=[];fwl=[]
    for k,ax in enumerate(AXES):
        d[ax+'_rank_sd']=zrank(d[ax]);d[ax+'_subtype_residual']=residual(d[ax+'_rank_sd'],d.subtype)
        x=d[ax+'_rank_sd'].to_numpy(); y=d.ZEB1_rank_sd.to_numpy(); g=d.subtype.to_numpy()
        for label,groups in [('marginal',None),('subtype_adjusted',g)]:
            out,f,X=fit_rank(x,y,groups);out.update(axis=ax,model=label,estimand='SD_global_ZEB1_midrank_per_SD_global_score_midrank')
            rows.append(out)
            pd.DataFrame({'term':X.columns,'coefficient':f.params,'se_HC3':f.bse,'p':f.pvalues,
                          'ci_low':f.conf_int().iloc[:,0],'ci_high':f.conf_int().iloc[:,1]}).to_csv(HERE/'results'/f'{ax}_{label}_model_terms.tsv',sep='\t',index=False)
        xr=residual(x,g);yr=residual(y,g);slope=float(xr@yr/(xr@xr))
        assert abs(slope-rows[-1]['effect'])<1e-12
        fwl.append({'axis':ax,'adjusted_slope_FWL':slope,'adjusted_slope_statsmodels':rows[-1]['effect'],
                    'absolute_difference':abs(slope-rows[-1]['effect']),'partial_rank_rho':cor(xr,yr)})
        total=float(x@y/len(x));ww=float(xr@yr/len(x));between=total-ww
        decomp.append({'axis':ax,'total_rank_covariance':total,'within_subtype_covariance':ww,
                       'between_subtype_covariance':between,'subtype_R2_score_rank':float(1-np.var(xr)/np.var(x)),
                       'subtype_R2_ZEB1_rank':float(1-np.var(yr)/np.var(y)),'partial_rank_rho':cor(xr,yr)})
        dummy=pd.get_dummies(d.subtype,drop_first=True,dtype=float)
        X=pd.concat([pd.DataFrame({'intercept':np.ones(len(x)),'score':x}),dummy,
                     dummy.mul(x,axis=0).rename(columns=lambda c:'interaction_'+c)],axis=1)
        f=sm.OLS(y,X).fit(cov_type='HC3',use_t=True)
        ids=[i for i,c in enumerate(X) if str(c).startswith('interaction_')]
        R=np.eye(len(X.columns))[ids];w=f.wald_test(R,use_f=True,scalar=True)
        interactions.append({'axis':ax,'n_patients':len(d),'n_subtypes':17,'F_HC3':float(w.statistic),
                             'df_num':int(w.df_num),'df_denom':float(w.df_denom),'p':float(w.pvalue),
                             'max_leverage':float(f.get_influence().hat_matrix_diag.max()),
                             'estimand':'heterogeneity_of_global_standardized_rank_slopes_across_17_categories'})
        for j,(sub,z) in enumerate(d.groupby('subtype',sort=True)):
            xx=z[ax].to_numpy();yy=z.ZEB1.to_numpy();rho=float(stats.spearmanr(xx,yy).statistic)
            lo,hi,nv,_=bootstrap_rho(xx,yy,SEED+k*100+j)
            pp,method=local_permutation(xx,yy,SEED+k*100+j+1000)
            glob,_,_=fit_rank(z[ax+'_rank_sd'].to_numpy(),z.ZEB1_rank_sd.to_numpy())
            within.append({'axis':ax,'subtype':sub,'n_patients':len(z),'spearman_rho':rho,'ci_low':lo,'ci_high':hi,
                           'p':pp,'p_method':method,'n_valid_bootstrap':nv,'global_rank_slope':glob['effect'],
                           'global_rank_slope_CI_low':glob['ci_low'],'global_rank_slope_CI_high':glob['ci_high'],
                           'status':'ESTIMATED_SMALL_N_DESCRIPTIVE' if len(z)<10 else 'ESTIMATED',
                           'interval_method':'10000_patient_percentile_bootstrap_recomputed_local_midranks'})
    leave_groups=[]
    for ax in AXES:
        for sub in sorted(d.subtype.unique()):
            keep=d.subtype.ne(sub);z=d[keep]
            out,_,_=fit_rank(z[ax+'_rank_sd'].to_numpy(),z.ZEB1_rank_sd.to_numpy(),z.subtype.to_numpy())
            out.update(axis=ax,omitted_subtype=sub,n_omitted_patients=int((~keep).sum()),
                       role='ALL17_LEAVE_ONE_SOURCE_GROUP_OUT_SENSITIVITY_GLOBAL_RANKS_HELD_FIXED_NO_PRIMARY_REDEFINITION')
            leave_groups.append(out)
    save(pd.DataFrame(leave_groups),'BULK34_LEAVE_ONE_SUBTYPE_OUT_SENSITIVITY.tsv')
    rr=pd.DataFrame(rows);rr['BH4_all_display_tests']=multipletests(rr.p,method='fdr_bh')[1]
    mask=rr.model.eq('subtype_adjusted');rr.loc[mask,'BH2_primary_adjusted']=multipletests(rr.loc[mask,'p'],method='fdr_bh')[1]
    ww=pd.DataFrame(within);ww['BH34_within_subtype']=multipletests(ww.p,method='fdr_bh')[1]
    ii=pd.DataFrame(interactions);ii['BH2_interactions']=multipletests(ii.p,method='fdr_bh')[1]
    save(rr,'BULK_MARGINAL_ADJUSTED_EFFECTS.tsv');save(ww,'ALL34_WITHIN_SUBTYPE_ASSOCIATIONS.tsv')
    save(ii,'SUBTYPE_INTERACTION_HETEROGENEITY.tsv');save(pd.DataFrame(decomp),'RANK_COVARIANCE_DECOMPOSITION.tsv')
    save(pd.DataFrame(fwl),'FWL_NUMERIC_READBACK.tsv')
    d.to_csv(HERE/'source_data'/'ALL1309_RANKS_AND_SUBTYPE_RESIDUALS.tsv',sep='\t',index=False,float_format='%.12g')
    write_json(HERE/'results'/'BULK_RESULT_GATE.json',{'status':'PASS','n_patients':1309,'n_unique_subjects':1309,
        'n_subtypes':17,'missing_values_in_primary_inputs':0,'n_primary_adjusted_tests':2,'n_within_subtype_tests':34,
        'source_labels_unchanged':True,'FWL_numeric_readback_pass':True,'original_clinical_models_refit':False,
        'HC3_use_t':True,'rank_scaling_ddof':0,'small_subtype_intervals':'Exploratory, discrete and unstable for n<10; retained without filtering.'})
    print('BULK_COMPLETE',rr[['axis','model','effect','ci_low','ci_high','p']].to_json(orient='records'),flush=True)

def normal_stage():
    s=HERE/'source_inputs';md=pd.read_csv(s/'NORMAL_DONOR_STATE_METADATA.tsv',sep='\t')
    full=pd.read_csv(s/'NORMAL_DONOR_STATE_COUNTS.tsv.gz',sep='\t',index_col=0)
    assert full.index.is_unique and set(full.columns)==set(md.pseudobulk_id)
    assert np.allclose(full.sum().to_numpy(),md.set_index('pseudobulk_id').loc[full.columns,'lib_size'].to_numpy())
    universe=pd.read_csv(s/'cortical_vs_early_NORMAL_SELECTION_LEDGER.tsv',sep='\t').symbol.tolist()
    assert len(universe)==25288 and len(set(universe))==25288
    counts=full.loc[universe];ranks=counts.rank(axis=0,method='average')/len(universe)
    members=pd.read_csv(s/'fixed_members.tsv',sep='\t');held=pd.read_csv(s/'normal_heldout.tsv',sep='\t')
    original=pd.read_csv(s/'normal_state_scores.tsv',sep='\t').set_index('pseudobulk_id')
    normal=md.copy().set_index('pseudobulk_id')
    normal['ZEB1_UMI_sum']=full.loc['ZEB1'];normal['ZEB1_log2_CPM_plus1']=np.log2(full.loc['ZEB1']/normal.lib_size*1e6+1)
    pairrows=[];gene_rows=[];membership=[];member_sets=[];qa=[]
    st=md.loc[md.n_cells>=30].set_index(['donor','mapped_state']).pseudobulk_id.to_dict()
    for ax in AXES:
        up=members.loc[(members.axis==ax)&(members.direction=='up'),'symbol'].tolist()
        dn=members.loc[(members.axis==ax)&(members.direction=='down'),'symbol'].tolist()
        score=ranks.loc[up].mean()-ranks.loc[dn].mean();normal[ax]=score
        delta=float(np.max(np.abs(score.reindex(original.index)-original[ax])));assert delta<1e-12
        qa.append({'axis':ax,'max_score_difference_from_frozen':delta})
        pairs=[]
        for donor in sorted(md.donor.unique()):
            required=['DP','Early_DN'] if ax==AXES[0] else ['DP','CD4SP','CD8SP']
            eligible=all((donor,x) in st for x in required)
            if not eligible:
                pairrows.append({'axis':ax,'donor':donor,'status':'UNESTIMATED_INELIGIBLE_STATE_LT30_CELLS','fixed_member_difference':np.nan});continue
            plus=[st[(donor,'DP')]] if ax==AXES[0] else [st[(donor,'CD4SP')],st[(donor,'CD8SP')]]
            minus=[st[(donor,'Early_DN')]] if ax==AXES[0] else [st[(donor,'DP')]]
            diff=float(score[plus].mean()-score[minus].mean());pairs.append((donor,plus,minus))
            z=normal.ZEB1_log2_CPM_plus1
            pairrows.append({'axis':ax,'donor':donor,'status':'ELIGIBLE','fixed_member_difference':diff,
                             'ZEB1_log2_CPM_plus1_paired_difference':float(z[plus].mean()-z[minus].mean()),
                             'n_cells_plus':';'.join(str(int(normal.loc[x,'n_cells'])) for x in plus),
                             'n_cells_minus':';'.join(str(int(normal.loc[x,'n_cells'])) for x in minus)})
        for direction,genes in [('up',up),('down',dn)]:
            fixed=set(genes)
            for h in held.loc[held.axis==ax].itertuples():
                fold=set(str(getattr(h,direction+'_symbols')).split(';'));inter=len(fixed&fold);union=len(fixed|fold)
                membership.append({'axis':ax,'direction':direction,'heldout_donor':h.heldout_donor,
                                   'fixed_n':len(fixed),'fold_n':len(fold),'overlap_n':inter,'fraction_fixed_retained':inter/len(fixed),
                                   'jaccard':inter/union,'heldout_difference_existing':h.heldout_difference})
                for gene in sorted(fixed|fold): member_sets.append({'axis':ax,'direction':direction,'heldout_donor':h.heldout_donor,'symbol':gene,'in_fixed':gene in fixed,'in_fold':gene in fold})
            for gene in genes:
                plus=ranks.loc[[g for g in up if g!=gene]].mean();minus=ranks.loc[[g for g in dn if g!=gene]].mean();leave=plus-minus
                diffs=[float(leave[hi].mean()-leave[lo].mean()) for _,hi,lo in pairs]
                gene_rows.append({'axis':ax,'omitted_direction':direction,'omitted_gene':gene,
                                  'n_up_used':49 if direction=='up' else 50,'n_down_used':49 if direction=='down' else 50,
                                  'n_eligible_donors':len(pairs),'minimum_paired_difference':min(diffs),
                                  'median_paired_difference':float(np.median(diffs)),'all_donor_directions_preserved':all(x>0 for x in diffs),
                                  'spearman_with_full_score_across30_profiles':float(stats.spearmanr(leave,score).statistic)})
    normal.reset_index().to_csv(HERE/'source_data'/'NORMAL30_DONOR_STATE_SCORES_AND_ZEB1.tsv',sep='\t',index=False)
    save(pd.DataFrame(pairrows),'NORMAL10_DONOR_CONTRAST_LEDGER.tsv');save(pd.DataFrame(gene_rows),'NORMAL200_FIXED_MEMBER_LEAVE_ONE_OUT.tsv')
    save(pd.DataFrame(membership),'NORMAL18_HELDOUT_MEMBERSHIP_STABILITY.tsv');save(pd.DataFrame(member_sets),'NORMAL_HELDOUT_MEMBER_PRESENCE_LEDGER.tsv')
    save(pd.DataFrame(qa),'NORMAL_SCORE_NUMERIC_READBACK.tsv')
    corrrows=[]
    for label,q in [('all30_profiles_descriptive',normal),('eligible_ge30_profiles_descriptive',normal[normal.n_cells>=30])]:
        corrrows.append({'resource':label,'n_profiles':len(q),'n_donors':q.donor.nunique(),
                         'spearman_rho_between_axes':float(stats.spearmanr(q[AXES[0]],q[AXES[1]]).statistic),
                         'inference':'DESCRIPTIVE_REPEATED_DONOR_STATES_NO_CELL_OR_PROFILE_LEVEL_P'})
    save(pd.DataFrame(corrrows),'NORMAL_INTERSCORE_CORRELATION.tsv')
    overlap=[]
    for da in ['up','down']:
        for db in ['up','down']:
            a=set(members.loc[(members.axis==AXES[0])&(members.direction==da),'symbol'])
            b=set(members.loc[(members.axis==AXES[1])&(members.direction==db),'symbol'])
            overlap.append({'early_direction':da,'later_direction':db,'n_shared_genes':len(a&b),'shared_genes':';'.join(sorted(a&b))})
    save(pd.DataFrame(overlap),'FIXED_MODULE_CROSS_AXIS_GENE_OVERLAP.tsv')
    nsum=[]
    pr=pd.DataFrame(pairrows)
    for ax in AXES:
        q=pr[(pr.axis==ax)&(pr.status=='ELIGIBLE')]
        for field in ['fixed_member_difference','ZEB1_log2_CPM_plus1_paired_difference']:
            v=q[field].to_numpy();n=len(v);mean=float(v.mean());se=float(stats.sem(v));half=float(stats.t.ppf(.975,n-1)*se)
            nsum.append({'axis':ax,'quantity':field,'n_donors':n,'mean_difference':mean,'ci_low':mean-half,'ci_high':mean+half,
                'median_difference':float(np.median(v)),'min_difference':float(v.min()),'max_difference':float(v.max()),
                'n_positive_donor_differences':int((v>0).sum()),'interval_method':'descriptive_donor_mean_t_interval_not_independent_signature_validation'})
    save(pd.DataFrame(nsum),'NORMAL_DONOR_CONTRAST_SUMMARY.tsv')
    print('NORMAL_COMPLETE',pd.DataFrame(pairrows).groupby('axis').fixed_member_difference.agg(['count','min','median','max']).to_json(),flush=True)

def cross_stage(project):
    s=HERE/'source_inputs';lim=pd.read_csv(s/'lim_scores.tsv',sep='\t'); lr=pd.read_csv(s/'lim_results.tsv',sep='\t')
    primary=lr[lr.model.eq('primary_spearman')];rows=[];long=[];compositions=[];sensitivity=[];correlations=[]
    b=pd.read_csv(s/'bulk.tsv',sep='\t').merge(pd.read_csv(s/'maturity.tsv',sep='\t'),on='sample_id',validate='one_to_one')
    frames={'main1309':(b,'sample_id','ZEB1',AXES),
            'Lim_discovery21':(lim[lim.cohort_analysis.eq('discovery')],'patient','ZEB1_log2_CPM_plus1',AXES),
            'Lim_validation28':(lim[lim.cohort_analysis.eq('validation')],'patient','ZEB1_log2_CPM_plus1',AXES)}
    for label,(df,pid,ycol,xcols) in list(frames.items()):
        assert df[pid].is_unique
        for ax in AXES:
            a=primary[(primary.cohort==('discovery' if 'discovery' in label else 'validation'))&(primary.axis==ax)]
            if label.startswith('Lim'):
                r=a.iloc[0];rho=float(stats.spearmanr(df[ycol],df[ax]).statistic);assert abs(rho-r.effect)<1e-9
                rows.append({'cohort':label,'axis':ax,'n_patients':len(df),'effect':rho,'ci_low':r.ci_low,'ci_high':r.ci_high,'p':r.p,
                             'original_q':r.BH4_primary,'original_family':'BH4_Lim_primary','uncertainty':'saved_10000_patient_percentile_bootstrap',
                             'assay':'malignant_patient_pseudobulk_UMI','score_members_used':100,'status':'REUSED_SOURCE_CHECKED'})
    for label,alias,ycol,assay in [('NOPHO108','nopho','ZEB1_log2_CPM_plus1','bulk_raw_read_counts'),
                                  ('AIIMS35','aiims','ZEB1_log2_FPKM_plus1','bulk_FPKM'),
                                  ('India25','india','ZEB1_log2_counts_plus1','bulk_author_normalized_counts')]:
        df=pd.read_csv(s/(alias+'.tsv'),sep='\t');j=json.loads((s/(alias+'_results.json')).read_text(encoding='utf-8'))
        pid=df.columns[0];assert df[pid].is_unique
        frames[label]=(df,pid,ycol,[a+'_raw' for a in AXES])
        for r in j['results']:
            ax=r['axis'];rho=float(stats.spearmanr(df[ycol],df[ax+'_raw']).statistic);assert abs(rho-r['spearman_rho'])<1e-9
            rows.append({'cohort':label,'axis':ax,'n_patients':len(df),'effect':rho,'ci_low':r['bootstrap95CI_low'],
                         'ci_high':r['bootstrap95CI_high'],'p':r['p_two_sided'],'original_q':r['BH2_q'],
                         'original_family':'BH2_'+label,'uncertainty':'saved_10000_patient_percentile_bootstrap',
                         'assay':assay,'score_members_used':r['members_up_used']+r['members_down_used'],'status':'REUSED_SOURCE_CHECKED'})
    out=pd.DataFrame(rows);out['BH10_exploratory_harmonization']=multipletests(out.p,method='fdr_bh')[1]
    out['direction']=np.where(out.effect>0,'positive','negative');out['CI_excludes_zero']=(out.ci_low*out.ci_high)>0
    save(out,'CROSS10_COHORT_EFFECTS_CI.tsv')
    for label,(df,pid,ycol,xcols) in frames.items():
        correlations.append({'cohort':label,'n_patients':len(df),'spearman_rho_between_axes':float(stats.spearmanr(df[xcols[0]],df[xcols[1]]).statistic)})
        for ax,xc in zip(AXES,xcols):
            for r in df.itertuples(index=False):
                long.append({'cohort':label,'patient_id':str(getattr(r,pid)),'axis':ax,'ZEB1_observed':float(getattr(r,ycol)),
                             'score_observed':float(getattr(r,xc))})
    save(pd.DataFrame(correlations),'PATIENT_RESOURCE_INTERSCORE_CORRELATION.tsv')
    pd.DataFrame(long).to_csv(HERE/'source_data'/'PATIENT_RESOURCE_SCORE_VALUES_LONG.tsv',sep='\t',index=False)
    for cohort,df in lim.groupby('cohort_analysis'):
        for field in ['subtype','response','tissue','timepoint']:
            for category,n in df[field].value_counts(dropna=False).items():
                compositions.append({'cohort':'Lim_'+cohort,'annotation':field,'category':str(category),'n_patients':int(n),'fraction':n/len(df)})
        for ax in AXES:
            x=zrank(df[ax]);y=zrank(df.ZEB1_log2_CPM_plus1)
            for field in ['subtype','response']:
                groups=df[field].fillna('MISSING_SOURCE_LABEL').to_numpy();cnt=pd.Series(groups).value_counts();inform=np.array([cnt[g]>1 for g in groups])
                # A singleton fixed-effect group has zero score information and unit leverage.
                # Preserve its original ranks and report its exclusion from the HC3 fit only.
                f,_,_=fit_rank(x[inform],y[inform],groups[inform]);f.update(cohort='Lim_'+cohort,axis=ax,adjustment=field,
                    n_source_patients=len(df),n_singleton_patients_without_slope_information=int((~inform).sum()),n_source_categories=len(cnt),
                    partial_rank_rho=cor(residual(x,groups),residual(y,groups)),
                    status='EXPLORATORY_SOURCE_CATEGORY_CONDITIONING_NO_LABEL_MERGING')
                sensitivity.append(f)
    ai=pd.read_csv(s/'aiims_native.tsv',sep='\t')
    for category,n in ai.source_disease_state.value_counts().items():
        compositions.append({'cohort':'AIIMS35','annotation':'source_disease_state','category':category,'n_patients':int(n),'fraction':n/len(ai)})
    for ax in AXES:
        f,_,_=fit_rank(zrank(ai[ax+'_raw']),zrank(ai.ZEB1_log2_FPKM_plus1),ai.source_disease_state.to_numpy())
        f.update(cohort='AIIMS35',axis=ax,adjustment='source_disease_state',n_source_patients=35,
            n_singleton_patients_without_slope_information=0,n_source_categories=3,
            partial_rank_rho=cor(residual(zrank(ai[ax+'_raw']),ai.source_disease_state),residual(zrank(ai.ZEB1_log2_FPKM_plus1),ai.source_disease_state)),
            status='EXPLORATORY_SOURCE_CATEGORY_CONDITIONING_NO_LABEL_MERGING');sensitivity.append(f)
    npth=s/'GSE272023_ALL_SAMPLE_HEADERS.json'
    original_header=project/'tasks/ZEB1_MATURATION_CLINICAL_STRENGTHEN_20261009/public_source_qualification/recent_clinical2025_2026/GSE272023_ALL_SAMPLE_HEADERS.json'
    if not npth.exists():shutil.copy2(original_header,npth)
    header=json.loads(npth.read_text(encoding='utf-8'));linked=[]
    for sample in header['samples']:
        fields=sample['fields']
        pid=next((x.split(' = ',1)[1] for x in fields if x.startswith('!Sample_description = TALL')),None)
        cimp=next((x.split('cimp subgroup: ',1)[1] for x in fields if 'cimp subgroup: ' in x),None)
        if pid is not None:linked.append({'patient_id':pid,'source_CIMP':cimp,'gsm':sample['gsm']})
    ann=pd.DataFrame(linked);assert len(ann)==108 and ann.patient_id.is_unique and ann.source_CIMP.notna().all()
    npdf=pd.read_csv(s/'nopho.tsv',sep='\t').merge(ann,on='patient_id',how='outer',validate='one_to_one')
    assert len(npdf)==108 and npdf.notna().all().all()
    npdf.to_csv(HERE/'source_data'/'NOPHO108_EXACT_SOURCE_CIMP_JOIN.tsv',sep='\t',index=False)
    for category,n in npdf.source_CIMP.value_counts().items():
        compositions.append({'cohort':'NOPHO108','annotation':'source_CIMP','category':category,'n_patients':int(n),'fraction':n/len(npdf)})
    for ax in AXES:
        f,_,_=fit_rank(zrank(npdf[ax+'_raw']),zrank(npdf.ZEB1_log2_CPM_plus1),npdf.source_CIMP.to_numpy())
        f.update(cohort='NOPHO108',axis=ax,adjustment='source_CIMP',n_source_patients=108,
            n_singleton_patients_without_slope_information=0,n_source_categories=2,
            partial_rank_rho=cor(residual(zrank(npdf[ax+'_raw']),npdf.source_CIMP),residual(zrank(npdf.ZEB1_log2_CPM_plus1),npdf.source_CIMP)),
            status='EXPLORATORY_SOURCE_PHENOTYPE_CONDITIONING_NOT_MOLECULAR_SUBTYPE_ADJUSTMENT');sensitivity.append(f)
    save(pd.DataFrame(compositions),'SOURCE_LINKED_COHORT_COMPOSITIONS.tsv')
    sn=pd.DataFrame(sensitivity);sn['BH12_composition_sensitivity']=multipletests(sn.p,method='fdr_bh')[1]
    save(sn,'CROSS_SOURCE_COMPOSITION_CONDITIONED_EFFECTS.tsv')
    # A same-study split comparison can be formalized without pooling different platforms.
    dv=[];gen=np.random.default_rng(SEED+500)
    for ax in AXES:
        a=lim[lim.cohort_analysis.eq('discovery')];b=lim[lim.cohort_analysis.eq('validation')]
        x=np.concatenate([a[ax],b[ax]]);y=np.concatenate([a.ZEB1_log2_CPM_plus1,b.ZEB1_log2_CPM_plus1]);n=len(a)
        observed=float(stats.spearmanr(b[ax],b.ZEB1_log2_CPM_plus1).statistic-stats.spearmanr(a[ax],a.ZEB1_log2_CPM_plus1).statistic)
        # Paired-patient bootstrap separately within each original cohort; no patient pooling.
        _,_,_,ba=bootstrap_rho(a[ax],a.ZEB1_log2_CPM_plus1,SEED+700+AXES.index(ax))
        _,_,_,bb=bootstrap_rho(b[ax],b.ZEB1_log2_CPM_plus1,SEED+800+AXES.index(ax))
        boot=bb-ba;ci=np.quantile(boot[np.isfinite(boot)],[.025,.975])
        # Fisher transform difference is an asymptotic sensitivity for independent Spearman coefficients.
        ra=float(stats.spearmanr(a[ax],a.ZEB1_log2_CPM_plus1).statistic);rb=float(stats.spearmanr(b[ax],b.ZEB1_log2_CPM_plus1).statistic)
        z=(np.arctanh(rb)-np.arctanh(ra))/np.sqrt(1/(len(a)-3)+1/(len(b)-3))
        dv.append({'axis':ax,'discovery_n':len(a),'validation_n':len(b),'validation_minus_discovery_rho':observed,
                   'ci_low':float(ci[0]),'ci_high':float(ci[1]),'p_Fisher_approximation':float(2*stats.norm.sf(abs(z))),
                   'method_boundary':'separate_cohort_patient_bootstrap_difference; Fisher_z_for_Spearman_is_approximate; composition_and_gene_universe_differ'})
    dv=pd.DataFrame(dv);dv['BH2_difference_approximate']=multipletests(dv.p_Fisher_approximation,method='fdr_bh')[1]
    save(dv,'LIM_DISCOVERY_VALIDATION_EFFECT_DIFFERENCE.tsv')
    table=np.array([[13,8],[5,23]]);fr=stats.fisher_exact(table)
    write_json(HERE/'results'/'LIM_RESPONSE_COMPOSITION_DIFFERENCE.json',{'discovery_IF':13,'discovery_responsive':8,
        'validation_IF':5,'validation_responsive':23,'Fisher_two_sided_p':float(fr.pvalue),'odds_ratio':float(fr.statistic),
        'role':'Exploratory evidence that source cohorts have different response composition; does not explain the association causally.'})
    print('CROSS_COMPLETE',out[['cohort','axis','effect','ci_low','ci_high']].to_json(orient='records'),flush=True)

def bulk_member_stage(project,bulk_tpm=None):
    s=HERE/'source_inputs';rankpath=s/'BULK_FIXED_AND_HELDOUT_MEMBER_PERCENTILE_RANKS.tsv.gz'
    common=pd.read_csv(s/'cortical_vs_early_NORMAL_SELECTION_LEDGER.tsv',sep='\t').symbol.tolist()
    mem=pd.read_csv(s/'fixed_members.tsv',sep='\t');h=pd.read_csv(s/'normal_heldout.tsv',sep='\t')
    targets=set(mem.symbol)
    for col in ['up_symbols','down_symbols']:
        for x in h[col]: targets.update(x.split(';'))
    d=pd.read_csv(s/'bulk.tsv',sep='\t').merge(pd.read_csv(s/'maturity.tsv',sep='\t'),on='sample_id',validate='one_to_one')
    ids=d.sample_id.tolist();path=Path(bulk_tpm) if bulk_tpm else project/'00_data/00_raw/TALL_X01_tpm.tsv'
    if not rankpath.exists():
        gm=pd.read_csv(s/'GENE_MAPPING_LEDGER.tsv',sep='\t');gm=gm[gm.eligible & gm.symbol.isin(common)]
        first=pd.read_csv(path,sep='\t',nrows=0).columns[0];blocks=[]
        for block in pd.read_csv(path,sep='\t',index_col=0,usecols=[first]+ids,chunksize=400):
            block.index=block.index.astype(str).str.split('.').str[0]
            blocks.append(block.loc[block.index.isin(set(gm.stable))])
        bulk=pd.concat(blocks)
        if bulk.index.has_duplicates:bulk=bulk.groupby(level=0,sort=False).max()
        sy=gm.set_index('stable').symbol.reindex(bulk.index);assert sy.notna().all()
        bulk.index=sy.to_numpy();bulk=bulk.loc[common,ids];assert bulk.index.is_unique
        selected=[]
        for start in range(0,len(ids),32):
            r=bulk.iloc[:,start:start+32].rank(axis=0,method='average')/len(common)
            selected.append(r.loc[sorted(targets)])
        rr=pd.concat(selected,axis=1);rr.to_csv(rankpath,sep='\t',compression='gzip')
        receipt={'input':str(path),'sha256':sha(path),'n_original_ranking_symbols':len(common),'n_patient_columns':len(ids),
                 'n_fixed_or_existing_fold_genes':len(rr),'exact_normalization':'average_within_sample_TPM_ranks/25288',
                 'snapshot':str(rankpath.relative_to(HERE)),'snapshot_sha256':sha(rankpath),
                 'reference_source_script':'qualify_multigene_maturity_proxy.py','duplicates':'maximum_TPM_per_stable_ID_as_original'}
        accepted=json.loads((s/'QUALIFICATION.json').read_text(encoding='utf-8'))
        expected=next(v for k,v in accepted['source_sha256'].items() if Path(k).name==path.name)
        assert receipt['sha256']==expected
        write_json(HERE/'BULK_MEMBER_RANK_EXTRACTION_RECEIPT.json',receipt)
    rr=pd.read_csv(rankpath,sep='\t',index_col=0).loc[:,ids]
    yr=zrank(d.ZEB1);g=d.subtype.to_numpy();yrres=residual(yr,g);rows=[];fold=[];readback=[]
    for ax in AXES:
        up=mem.loc[(mem.axis==ax)&(mem.direction=='up'),'symbol'].tolist();dn=mem.loc[(mem.axis==ax)&(mem.direction=='down'),'symbol'].tolist()
        full=rr.loc[up].mean()-rr.loc[dn].mean();diff=float(np.max(np.abs(full.to_numpy()-d[ax+'_raw'].to_numpy())))
        assert diff<1e-12;readback.append({'axis':ax,'n_patients':len(d),'max_raw_score_difference_from_frozen':diff})
        for direction,genes in [('up',up),('down',dn)]:
            for gene in genes:
                score=rr.loc[[z for z in up if z!=gene]].mean()-rr.loc[[z for z in dn if z!=gene]].mean()
                xr=zrank(score);res=residual(xr,g)
                rows.append({'axis':ax,'omitted_direction':direction,'omitted_gene':gene,'n_patients':len(d),
                    'rho_ZEB1_marginal':float(stats.spearmanr(score,d.ZEB1).statistic),
                    'adjusted_global_rank_slope':float(res@yrres/(res@res)),
                    'partial_rank_rho':cor(res,yrres),'rho_with_full_fixed_score':float(stats.spearmanr(score,full).statistic),
                    'max_absolute_raw_score_change':float(np.max(np.abs(score-full))),
                    'inference':'DESCRIPTIVE_SINGLE_FIXED_MEMBER_OMISSION_NO_RESELECTION_OR_NEW_P_FAMILY'})
        for r in h.loc[h.axis==ax].itertuples():
            score=rr.loc[r.up_symbols.split(';')].mean()-rr.loc[r.down_symbols.split(';')].mean();xr=zrank(score)
            f,_,_=fit_rank(xr,yr,g);f.update(axis=ax,heldout_normal_donor=r.heldout_donor,
                normal_heldout_difference=r.heldout_difference,rho_with_full_fixed_score=float(stats.spearmanr(score,full).statistic),
                partial_rank_rho=cor(residual(xr,g),yrres),role='EXISTING_INTERNAL_NORMAL_HOLDOUT_MEMBERS_APPLIED_WITHOUT_NEW_SELECTION')
            fold.append(f)
    save(pd.DataFrame(rows),'BULK200_FIXED_MEMBER_LEAVE_ONE_OUT.tsv');save(pd.DataFrame(fold),'BULK9_EXISTING_DONOR_FOLD_SCORE_SENSITIVITY.tsv')
    save(pd.DataFrame(readback),'BULK_SCORE_RECONSTRUCTION_READBACK.tsv')
    print('BULK_MEMBER_COMPLETE',pd.DataFrame(rows).groupby('axis')[['adjusted_global_rank_slope','rho_with_full_fixed_score']].agg(['min','max']).to_json(),flush=True)

def main():
    p=argparse.ArgumentParser();p.add_argument('--project',type=Path,default=ROOT)
    p.add_argument('--stage',choices=['freeze','bulk','normal','cross','bulk-members','all'],default='all');p.add_argument('--bulk-tpm',type=Path)
    a=p.parse_args();start=datetime.datetime.now(datetime.timezone.utc).isoformat()
    for x in ['source_inputs','source_data','results']: (HERE/x).mkdir(exist_ok=True)
    if a.stage=='freeze' or (a.stage=='all' and not (HERE/'SOURCE_INPUT_READBACK.json').exists()):freeze_inputs(a.project)
    snapshot_gate=verify_snapshots()
    write_json(HERE/'INPUT_SNAPSHOT_GATE.json',snapshot_gate)
    if a.stage in ['bulk','all']:bulk_stage()
    if a.stage in ['normal','all']:normal_stage()
    if a.stage in ['cross','all']:cross_stage(a.project)
    if a.stage in ['bulk-members','all']:bulk_member_stage(a.project,a.bulk_tpm)
    write_json(HERE/f'EXECUTION_RECEIPT_{a.stage}.json',{'status':'COMPLETE','stage':a.stage,'started_UTC':start,
        'finished_UTC':datetime.datetime.now(datetime.timezone.utc).isoformat(),'seed':SEED,'script_sha256':sha(__file__),
        'environment':{'python':platform.python_version(),'numpy':np.__version__,'pandas':pd.__version__,
                       'scipy':scipy.__version__,'statsmodels':statsmodels.__version__},
        'design_gate':'DESIGN_GATE.json','original_science_modified':False,'clinical_models_refit':False,
        'source_implementation':['https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.rankdata.html',
             'https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.bootstrap.html',
             'https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.permutation_test.html',
             'https://www.statsmodels.org/stable/generated/statsmodels.regression.linear_model.OLS.html',
             'https://www.statsmodels.org/stable/generated/statsmodels.regression.linear_model.OLSResults.HC3_se.html']})

if __name__=='__main__':main()

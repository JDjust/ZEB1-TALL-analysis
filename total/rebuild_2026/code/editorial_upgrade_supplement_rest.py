"""Data-rich one-page S1/S2/S4/S5/S7/S8 plates from frozen inputs."""
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from matplotlib.patches import Rectangle
from matplotlib.colors import TwoSlopeNorm

import editorial_rebuild_figures as b
import editorial_rebuild_supplement as s
from editorial_upgrade_1 import _coverage
from editorial_upgrade_347 import CMAP

D,SD,DATA,OUT=b.D,b.SD,b.DATA,s.OUT
Z1,Z2,ETP,GREY,LIGHT,DARK=b.Z1,b.Z2,b.ETP,b.GREY,b.LIGHT,b.DARK
read,panel,grid=b.read,b.panel,b.clean_grid

def save(f,n):
    for t in f.findobj(match=plt.Text):
        if t.get_text().strip() and t.get_fontsize()<6.5:t.set_fontsize(6.5)
    f.savefig(OUT/f'FigureS{n}.pdf');f.savefig(OUT/f'FigureS{n}.png',dpi=300)
    plt.close(f)

def detection_matrix():
    path=D/'a4_yayon/merged_thymus_visium_pediatric.h5ad'
    import anndata as ad
    x=ad.read_h5ad(path,backed='r')
    genes=['ZEB1','ZEB2','CD1A','CD34','LYL1']
    vi=x.var.feature_name.astype(str)
    ids=[x.var_names[vi.eq(g)][0] for g in genes]
    q=x[:,ids].X
    if hasattr(q,'toarray'):q=q.toarray()
    q=np.asarray(q)
    libs=x.obs.library_id.astype(str).to_numpy()
    tcols=[c for c in x.obs.columns if str(c).startswith('T_')]
    tfrac=x.obs[tcols].apply(pd.to_numeric,errors='coerce').sum(axis=1)
    total=pd.to_numeric(x.obs['tot_cell_abundance'],errors='coerce')
    t_enriched=(tfrac/total.replace(0,np.nan)>=.5).to_numpy()
    assert q.shape==(25593,5) and len(np.unique(libs))==16
    rows=[]
    for lib in sorted(np.unique(libs)):
        mask=(libs==lib)&t_enriched
        assert mask.sum()>=10
        row={'library_id':lib,'n_spots':int(mask.sum())}
        for j,gene in enumerate(genes):row[gene]=float((q[mask,j]>0).mean())
        rows.append(row)
    out=pd.DataFrame(rows)
    assert len(out)==16
    # The published 10%-of-spots eligibility summary must be reproduced.
    assert int((out.ZEB2>=.10).sum())==2
    x.file.close()
    folder=DATA/'editorial_visual_sources';folder.mkdir(parents=True,exist_ok=True)
    out.to_csv(folder/'S1D_Visium_library_detection.tsv',sep='\t',index=False)
    return out,genes

def sc1():
    det,genes=detection_matrix()
    states=read(D/'a4_yayon/a4a_stage_balance.tsv').merge(
        read(D/'a4_yayon/a4ab_stage_balance_and_cma.tsv')[['stage','cma_weighted_median']],
        on='stage',validate='one_to_one')
    assert len(states)==33
    f=b.mmfig(180,215);g=GridSpec(2,3,figure=f,height_ratios=[.72,1.62])
    a=f.add_subplot(g[0,:]);_coverage(a)
    bb=f.add_subplot(g[1,:2]);right=g[1,2].subgridspec(2,1,hspace=.65)
    c=f.add_subplot(right[0,0]);dd=f.add_subplot(right[1,0])
    states=states.sort_values(['pole','balance_median']).reset_index(drop=True)
    cols=['balance_median','dev_median','cma_weighted_median','n_cells','n_donors']
    z=states[cols].copy();z['n_cells']=np.log10(z.n_cells+1)
    z=(z-z.mean())/z.std(ddof=0)
    bb.imshow(z.to_numpy(),aspect='auto',cmap=CMAP,vmin=-2,vmax=2,interpolation='nearest')
    bb.set_yticks(range(33),states.stage.str.replace('T_','',regex=False),fontsize=6.5)
    bb.set_xticks(range(5),['B','D','CMA','log n cells','n donors'],
                  rotation=25,ha='right',fontsize=6.5)
    panel(bb,'B','All 33 author states','Column-standardized display; grey CMA cells unavailable')
    marker=read(SD/'S3/S3B_marker_directions.tsv')
    marker=marker[marker.gene.isin(['CD34','LYL1','CD1A'])]
    mat=marker.pivot(index='gene',columns='dataset',values='mid_minus_early')
    assert mat.shape[0]==3
    c.imshow(mat.to_numpy(float),aspect='auto',cmap=CMAP,vmin=-.9,vmax=.9)
    c.set_yticks(range(len(mat)),mat.index,fontsize=6.5)
    c.set_xticks(range(len(mat.columns)),[q.replace('GSE','G') for q in mat.columns],
                 rotation=30,ha='right',fontsize=6.5)
    panel(c,'C','Proxy-marker directions','Middle − early; grey = NA')
    dd.imshow(det[genes].to_numpy(),aspect='auto',cmap='Blues',vmin=0,vmax=1,
              interpolation='nearest')
    dd.set_yticks([0,3,7,11,15],[f'Lib {i+1}' for i in [0,3,7,11,15]],fontsize=6.5)
    dd.set_xticks(range(5),genes,rotation=35,ha='right',fontsize=6.5)
    panel(dd,'D','T-enriched Visium spots','16 libraries; ZEB2 ≥10% in 2/16')
    f.subplots_adjust(left=.19,right=.965,top=.94,bottom=.07,wspace=.50,hspace=.66)
    save(f,1)

def sc2():
    f=b.mmfig(180,205);g=GridSpec(3,2,figure=f)
    a,bb,c,d,e,ff=[f.add_subplot(g[i,j]) for i in range(3) for j in range(2)]
    x=read(SD/'revision3_robustness/crossfit_subtype_summary.tsv').sort_values('median_primary')
    for i,r in enumerate(x.itertuples()):
        a.plot([r.median_primary,r.median_crossfit],[i,i],color=LIGHT,lw=.8)
        a.scatter([r.median_primary,r.median_crossfit],[i,i],color=[GREY,Z2],s=9)
    a.set_yticks(range(17),[b.short_subtype(v) for v in x.subtype],fontsize=6.5)
    a.invert_yaxis();a.axvline(0,color=LIGHT,lw=.5);a.set_xlabel('Median R')
    panel(a,'A','Leave-subtype-out comparison','Shared baseline versus subtype-held-out spline');grid(a,'x')
    alt=read(SD/'revision3_robustness/alternative_coordinates_subtype_summary.tsv')
    definitions=['primary','minus_CD34','minus_LYL1','minus_CD1A','expanded']
    names=x.subtype.tolist();mat=alt.pivot(index='subtype',columns='definition',values='median_residual').loc[names,definitions]
    ranks=mat.rank(axis=0,ascending=True)
    bb.imshow(ranks.to_numpy(),aspect='auto',cmap='Blues',vmin=1,vmax=17)
    bb.set_yticks(range(17),[b.short_subtype(v) for v in names],fontsize=6.5)
    bb.set_xticks(range(5),['Primary','−CD34','−LYL1','−CD1A','Expanded'],
                   rotation=25,ha='right',fontsize=6.5)
    for i,name in enumerate(names):
        if name in ['BCL11B','TLX3']:
            bb.add_patch(Rectangle((-.48,i-.48),4.96,.96,fill=False,
                                   edgecolor=b.color_subtype(name),lw=.9))
    panel(bb,'B','Seventeen-subtype rank heatmap','Rank 1 = most ZEB2-skewed; middle ranks may move')
    models=read(SD/'revision3_robustness/model_comparison.tsv')
    models=models[models.model.isin(['development','subtype','combined'])]
    for k,(col,color,lab) in enumerate([('raw_r2',GREY,'Raw'),('adjusted_r2',ETP,'Adjusted'),('cv_r2_mean',DARK,'CV')]):
        c.scatter(models[col]*100,np.arange(3)+k*.16,s=17,color=color,label=lab)
    c.set_yticks(np.arange(3)+.16,models.model);c.set_xlabel('R² (%)')
    c.legend(frameon=False,fontsize=6.5,ncol=3)
    panel(c,'C','Model-complexity audit','Patient-level repeated five-fold CV');grid(c,'x')
    df=read(SD/'revision4_targeted/spline_df_2_3_4.tsv')
    for name,col in [('BCL11B',Z2),('ETP-like',ETP),('TLX3',Z1)]:
        q=df[df.subtype.eq(name)].sort_values('spline_df')
        d.plot(q.spline_df,q.median_residual,marker='o',ms=3,color=col,label=name)
    d.set_xticks([2,3,4]);d.set_xlabel('Natural spline df');d.set_ylabel('Median R')
    d.legend(frameon=False,fontsize=6.5)
    panel(d,'D','Prespecified spline flexibility','No setting selected for favorable result');grid(d)
    coef=read(D/'gate2_sensitivity/frozen_genes_batch_adjusted.tsv')
    e.scatter(coef.beta_primary,coef.beta_batch,s=7,color=Z2,alpha=.7)
    e.plot([-1.2,1.2],[-1.2,1.2],color=GREY,lw=.6,ls='--')
    e.set_xlabel('Primary βR');e.set_ylabel('Batch-adjusted βR')
    panel(e,'E','Batch effect audit','100/100 frozen coefficients preserve sign');grid(e)
    pt=read(SD/'F3/F3B-C_patient_level.tsv')
    normal=read(SD/'S3/S3C_normal_stage_units.tsv')
    ff.hist(pt.dev,bins=28,color='#D8DEE1',density=True,label='T-ALL 1309')
    ff.hist(pt.loc[pt.subtype.eq('BCL11B'),'dev'],bins=10,color=Z2,
            density=True,alpha=.55,label='BCL11B 18')
    lo,hi=normal.dev.min(),normal.dev.max()
    ff.axvspan(lo,hi,color=ETP,alpha=.15)
    ff.scatter(normal.dev,np.zeros(len(normal.dev))+.015,marker='|',s=27,color=ETP)
    ff.set_xlabel('Developmental proxy D');ff.set_ylabel('Density')
    ff.legend(frameon=False,fontsize=6.5)
    panel(ff,'F','Normal-reference domain','20 units are not 20 donors; 14/18 BCL11B outside')
    grid(ff,'y')
    f.subplots_adjust(left=.19,right=.97,top=.945,bottom=.07,wspace=.50,hspace=.75)
    save(f,2)

def sc4():
    orig=read(D/'gate2_program/frozen_ZEB_side50.tsv')
    hold=read(D/'gate2_sensitivity/heldout_ZEB_side50.tsv')
    assert len(orig)==len(hold)==100
    f=b.mmfig(180,175);g=GridSpec(2,2,figure=f)
    a,bb,c,d=[f.add_subplot(g[i,j]) for i in range(2) for j in range(2)]
    allg=read(D/'gate2_program/all_genes_residual_effects.tsv')
    a.scatter(allg.beta_R,allg.partial_R2,s=2,color=LIGHT,alpha=.55,rasterized=True)
    for side,col in [('ZEB1-side50',Z1),('ZEB2-side50',Z2)]:
        q=orig[orig.side.eq(side)]
        a.scatter(q.beta_R,q.partial_R2,s=9,color=col)
    a.set_xlabel('Residual βR');a.set_ylabel('Partial R²')
    panel(a,'A','Frozen 100-gene selection','24,356 eligible; 3,560 threshold-pass genes');grid(a)
    membership=[]
    for side in ['ZEB2-side50','ZEB1-side50']:
        one=set(orig.loc[orig.side.eq(side),'symbol']);two=set(hold.loc[hold.side.eq(side),'symbol'])
        assert len(one)==len(two)==50
        membership+=sorted(one|two)
    membership=list(dict.fromkeys(membership));assert len(membership)>=120
    columns=[('Original Z2',orig,'ZEB2-side50'),('Held-out Z2',hold,'ZEB2-side50'),
             ('Original Z1',orig,'ZEB1-side50'),('Held-out Z1',hold,'ZEB1-side50')]
    m=np.column_stack([[int(gene in set(t.loc[t.side.eq(side),'symbol'])) for gene in membership]
                       for _,t,side in columns])
    bb.imshow(m,aspect='auto',cmap='Greys',vmin=0,vmax=1,interpolation='nearest')
    bb.set_xticks(range(4),[q[0] for q in columns],rotation=26,ha='right',fontsize=6.5)
    bb.set_yticks([]);bb.axhline(len(set(orig.loc[orig.side.eq('ZEB2-side50'),'symbol'])|
                                 set(hold.loc[hold.side.eq('ZEB2-side50'),'symbol']))-.5,
                                 color=Z2,lw=1)
    panel(bb,'B','Frozen membership overlap','Binary matrix of original and held-out gene sets')
    cam=read(D/'gate2_sensitivity/heldout_matched_cameraPR.tsv')
    q=cam[cam.program.isin(['ZEB1-side50','ZEB2-side50'])]
    for i,r in enumerate(q.itertuples()):
        col=Z1 if 'ZEB1' in r.program else Z2
        c.plot([0,-np.log10(r.PValue)],[i,i],color=LIGHT,lw=1)
        c.scatter(-np.log10(r.PValue),i,s=25,color=col)
    c.set_yticks(range(len(q)),q.program);c.set_xlabel('−log10 CAMERA P')
    panel(c,'C','Held-out CAMERA in matched pairs','Within-cohort phenotype bridge');grid(c,'x')
    bt=read(D/'gate2_program/residual_t_primary_vs_pc15.tsv').dropna(subset=['t_primary','t_pc'])
    d.hexbin(bt.t_primary,bt.t_pc,gridsize=42,mincnt=1,cmap='Blues',linewidths=0,rasterized=True)
    d.set_xlim(*np.quantile(bt.t_primary,[.01,.99]));d.set_ylim(*np.quantile(bt.t_pc,[.01,.99]))
    d.set_xlabel('Primary gene t');d.set_ylabel('PC1–5-adjusted gene t')
    panel(d,'D','Technical-axis overadjustment','Diagnostic, not robustness confirmation');grid(d)
    f.subplots_adjust(left=.14,right=.975,top=.94,bottom=.09,wspace=.46,hspace=.75)
    save(f,4)

def sc5():
    pairs=read(SD/'revision3_robustness/bcl11b_etp_matched_pairs.tsv')
    cases=read(SD/'F4/F4C_GSE162280_cases.tsv').sort_values('log2_ZEB1_over_ZEB2')
    assert len(pairs)==18 and len(cases)==12
    f=b.mmfig(180,185);g=GridSpec(2,2,figure=f)
    a,bb,c,d=[f.add_subplot(g[i,j]) for i in range(2) for j in range(2)]
    q=pairs.sort_values('BCL11B_dev').reset_index(drop=True)
    right={sample:i for i,sample in enumerate(q.sort_values('ETP_dev').ETP_sample)}
    for i,r in q.iterrows():
        j=right[r.ETP_sample]
        a.plot([0,1],[i,j],color=Z2,alpha=.25+.6*(1-r.distance/.10),lw=.85)
        a.scatter([0,1],[i,j],s=6,color=[Z2,ETP])
    a.set_xticks([0,1],['BCL11B','Matched ETP']);a.set_yticks([])
    a.set_ylim(17.5,-.5)
    panel(a,'A','Bipartite matching','18 measured pairs; line opacity reflects |ΔD|')
    for i,r in q.iterrows():
        bb.plot([0,1],[r.BCL11B_residual,r.ETP_residual],color=LIGHT,lw=.8)
        bb.scatter([0,1],[r.BCL11B_residual,r.ETP_residual],s=9,color=[Z2,ETP])
    bb.set_xticks([0,1],['BCL11B','ETP-like']);bb.set_ylabel('Patient R')
    panel(bb,'B','Matched residual geometry','Same 18 fixed pairs');grid(bb)
    h=read(SD/'revision4_targeted/matched_BCL11B_ETP_hallmark_camera.tsv').sort_values('FDR').head(7)
    h=h.sort_values('FDR',ascending=False)
    yy=np.arange(len(h));c.scatter(-np.log10(h.FDR),yy,s=h.NGenes*.18,
                                  color=[Z2 if v=='Up' else Z1 for v in h.Direction])
    c.set_yticks(yy,[v.replace('HALLMARK_','').replace('_',' ')[:17] for v in h.pathway],fontsize=6.5)
    c.set_xlabel('−log10 BH FDR');panel(c,'C','Hallmark CAMERA context','Dot area = pathway genes; not mechanism');grid(c,'x')
    col1=cases.ZEB1_cpm.to_numpy();col2=cases.ZEB2_cpm.to_numpy()
    vals=np.column_stack([np.log2(col1+1),np.log2(col2+1),cases.log2_ZEB1_over_ZEB2])
    norm=np.clip(vals/5,-1,1)
    d.imshow(norm,aspect='auto',cmap=CMAP,vmin=-1,vmax=1,interpolation='nearest')
    d.set_yticks(range(12),[f'{r.case} · {str(r.phenotype_class)[:4]}' for r in cases.itertuples()],fontsize=6.5)
    d.set_xticks(range(3),['log2 ZEB1 CPM','log2 ZEB2 CPM','log2 ratio'],
                 rotation=25,ha='right',fontsize=6.5)
    for i,r in enumerate(cases.itertuples()):
        d.text(2.55,i,'fusion' if r.zeb2_partner=='yes' else 'enh.',
               ha='left',va='center',fontsize=6.5,color=DARK)
    d.set_xlim(-.5,2.9)
    panel(d,'D','Complete 12-case lesion matrix','Column-specific units; all ratios < 0')
    f.subplots_adjust(left=.18,right=.97,top=.94,bottom=.075,wspace=.64,hspace=.73)
    save(f,5)

def sc7():
    units=read(D/'a7_lineage/scores/a7a_unit_logcpm_min20.tsv')
    genes=read(D/'a7_lineage/scores/a7a_genes_min20.tsv')
    assert len(genes)==93 and genes.symbol.is_unique
    scores=read(D/'a7_lineage/scores/a7a_unit_scores_min20.tsv')
    f=b.mmfig(180,185);g=GridSpec(2,2,figure=f)
    a,bb,c,d=[f.add_subplot(g[i,j]) for i in range(2) for j in range(2)]
    order=['HSC/MPP/HSPC','lymphoid progenitor/CLP','T lineage','Myeloid']
    donors=sorted(scores.donor.unique())
    mat=scores.pivot_table(index='donor',columns='lineage',values='ZEB2_side50').reindex(index=donors,columns=order)
    a.imshow(mat.to_numpy(),aspect='auto',cmap=CMAP,vmin=-1.5,vmax=1.5)
    a.set_yticks(range(len(mat)),donors,fontsize=6.5)
    a.set_xticks(range(4),['HSPC','CLP','T','Myeloid'],rotation=20,ha='right',fontsize=6.5)
    panel(a,'A','Donor × lineage ZEB2 scores','Eligibility ≥20 cells; grey = unavailable')
    med=units.groupby('lineage')[genes.symbol.tolist()].median().reindex(order)
    z=med.T.loc[genes.sort_values(['side','rank']).symbol]
    z=(z-z.mean(axis=1).to_numpy()[:,None])/z.std(axis=1,ddof=0).replace(0,np.nan).to_numpy()[:,None]
    bb.imshow(z.to_numpy(),aspect='auto',cmap=CMAP,vmin=-2,vmax=2,interpolation='nearest')
    bb.set_yticks([22,68],['ZEB1 side','ZEB2 side'],fontsize=6.5)
    bb.set_xticks(range(4),['HSPC','CLP','T','Myeloid'],rotation=20,ha='right',fontsize=6.5)
    panel(bb,'B','Ninety-three frozen genes','Row-standardized lineage expression')
    paired=read(D/'a7_lineage/scores/a7a_donor_min20.tsv')
    # Use the existing Figure 6 paired source if shape differs between releases.
    if {'contrast','difference'}.issubset(paired.columns):
        for i,(name,q) in enumerate(paired.groupby('contrast')):
            c.scatter(q.difference,np.full(len(q),i),s=13,color=Z2)
        c.set_yticks(range(paired.contrast.nunique()),paired.contrast.unique(),fontsize=6.5)
    else:
        for i,lineage in enumerate(['Myeloid','HSC/MPP/HSPC','lymphoid progenitor/CLP']):
            v=(mat[lineage]-mat['T lineage']).dropna()
            c.scatter(v,np.full(len(v),i)+b.jitter(len(v),.04,i),s=11,color=Z2 if i==0 else ETP)
        c.set_yticks(range(3),['Myeloid−T','HSPC−T','CLP−T'])
    c.axvline(0,color=LIGHT,lw=.6);c.set_xlabel('Within-donor ZEB2-side difference')
    panel(c,'C','Paired lineage differences','Donor is the independent unit');grid(c,'x')
    for i,(name,file) in enumerate([('≥20',D/'a7_lineage/scores/a7a_unit_scores_min20.tsv'),
                                    ('≥10',D/'a7_lineage/scores/a7a_unit_scores_min10.tsv')]):
        q=read(file).groupby('lineage').ZEB2_side50.median().reindex(order)
        d.plot(range(4),q,marker='o',ms=3,color=Z2 if i==0 else GREY,label=name)
    d.set_xticks(range(4),['HSPC','CLP','T','Myeloid'],rotation=20,ha='right')
    d.set_ylabel('Median ZEB2-side score');d.legend(frameon=False,fontsize=6.5)
    panel(d,'D','Cell-threshold sensitivity','Same frozen genes under ≥20 / ≥10 cells');grid(d)
    f.subplots_adjust(left=.16,right=.97,top=.94,bottom=.08,wspace=.50,hspace=.73)
    save(f,7)

def sc8():
    x=read(D/'a7b_malignant/a7b_patient_scores.tsv').sort_values('B')
    assert len(x)==15
    paired=read(SD/'F5/F5D-E_Lim_Day0_state_pairs.tsv')
    day28=read(SD/'S7/S7D_Day0_Day28_pairs.tsv')
    assert len(paired)==41 and len(day28)==12
    f=b.mmfig(180,190);g=GridSpec(2,3,figure=f)
    a=f.add_subplot(g[0,0]);bg=g[0,1].subgridspec(1,2,wspace=.35)
    b1=f.add_subplot(bg[0,0]);b2=f.add_subplot(bg[0,1]);c=f.add_subplot(g[0,2])
    d=f.add_subplot(g[1,:2]);e=f.add_subplot(g[1,2])
    mat=x[['B','ZEB1_side50','ZEB2_side50']]
    a.imshow(mat.to_numpy(),aspect='auto',cmap=CMAP,vmin=-2,vmax=2)
    a.set_yticks(range(15),x.patient,fontsize=6.5)
    a.set_xticks(range(3),['Balance B','Z1 score','Z2 score'],rotation=30,ha='right',fontsize=6.5)
    panel(a,'A','Malignant patient matrix','15 author-labelled malignant pseudobulks')
    for ax,col,color in [(b1,'ZEB1_side50',Z1),(b2,'ZEB2_side50',Z2)]:
        ax.scatter(x.B,x[col],s=14,color=color)
        ax.axhline(0,color=LIGHT,lw=.5);ax.set_xlabel('Malignant B');grid(ax)
    b1.set_ylabel('Z1 score');b2.set_ylabel('Z2 score')
    panel(b1,'B','Z1','ρ +.73')
    panel(b2,'','Z2','ρ −.21')
    c.barh(range(15),x.n_cells,color=ETP)
    c.set_yticks(range(15),x.patient,fontsize=6.5);c.set_xscale('log')
    c.set_xticks([100,1000],['100','1000'])
    c.set_xlabel('Malignant cells');panel(c,'C','Patient cell counts','15 biological units');grid(c,'x')
    p=paired.sort_values('d_balance').reset_index(drop=True)
    for i,r in p.iterrows():
        d.plot([r.balance_neg,r.balance_pos],[i,i],color=ETP,alpha=.55,lw=.65)
        d.scatter([r.balance_neg,r.balance_pos],[i,i],s=5,color=[GREY,Z2],rasterized=True)
    d.set_yticks([]);d.set_xlabel('Day-0 malignant-state balance')
    panel(d,'D','ZBTB16 paired states','41 within-patient positive versus negative state pairs');grid(d,'x')
    for _,r in day28.iterrows():
        e.plot([0,1],[r.day0_balance,r.day28_balance],color=GREY,alpha=.65,lw=.75)
        e.scatter([0,1],[r.day0_balance,r.day28_balance],s=9,color=[ETP,Z2])
    e.set_xticks([0,1],['Day 0','Day 28']);e.set_ylabel('Patient balance')
    panel(e,'E','Day-28 boundary','12 pairs; no stable selection');grid(e)
    f.subplots_adjust(left=.17,right=.975,top=.94,bottom=.08,wspace=.55,hspace=.72)
    save(f,8)

if __name__=='__main__':
    for n,fn in [(1,sc1),(2,sc2),(4,sc4),(5,sc5),(7,sc7),(8,sc8)]:
        print('S',n,flush=True);fn()

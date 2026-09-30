"""One-page supplementary evidence plates S1–S10 from locked tables."""
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from matplotlib.patches import Rectangle
from scipy.stats import spearmanr
import editorial_rebuild_figures as m

SD=m.SD;D=m.D;read=m.read;panel=m.panel;grid=m.clean_grid
Z1=m.Z1;Z2=m.Z2;ETP=m.ETP;GREY=m.GREY;LIGHT=m.LIGHT;DARK=m.DARK;SOFT=m.SOFT
OUT=m.OUT.parent/'supplementary'

def save(fig,n):
    fig.subplots_adjust(left=.14,right=.97,top=.94,bottom=.07,wspace=.42,hspace=.63)
    for label in fig.findobj(match=matplotlib.text.Text):
        if label.get_text().strip() and label.get_fontsize()<6.5:label.set_fontsize(6.5)
    fig.savefig(OUT/f'FigureS{n}.pdf')
    fig.savefig(OUT/f'FigureS{n}.png',dpi=300)
    plt.close(fig)

def fig(n,h=180,rows=2,cols=2):
    f=m.mmfig(180,h);g=GridSpec(rows,cols,figure=f)
    return f,g

def sc1():
    f=m.mmfig(180,195);g=GridSpec(2,3,figure=f,height_ratios=[.85,1.65]);a=f.add_subplot(g[0,:]);b=f.add_subplot(g[1,:2]);c=f.add_subplot(g[1,2])
    x=read(SD/'Supplement_S1/source_S2_S2A_complete_stage_balance.tsv')
    for i,(name,t) in enumerate(x.groupby('dataset')):
        t=t.sort_values('order');a.plot(t['order'],t['balance'],lw=1.2,marker='o',ms=2.6,label=name,color=['#4A7189','#C0713F','#749486','#A39963'][i])
    panel(a,'A','Complete four-atlas stage profiles','Source-specific author stages; no pooled donor inference');a.set_xlabel('Within-source stage order');a.set_ylabel('Displayed balance');a.legend(frameon=False,fontsize=6.5,ncol=4,loc='upper right');grid(a)
    stage=read(D/'a4_yayon/a4a_stage_balance.tsv').sort_values('balance_median');yy=np.arange(len(stage))
    cols=[Z2 if q=='early_DN' else Z1 if q=='cortical_DP' else GREY for q in stage.pole]
    b.scatter(stage.balance_median,yy,c=cols,s=13);b.set_yticks(yy,[str(q).replace('T_','') for q in stage.stage],fontsize=6.5);b.axvline(0,color=LIGHT,lw=.5)
    panel(b,'B','All 33 eligible Yayon author states','Sorted by observed balance, not developmental time');b.set_xlabel('Stage median ZEB balance');grid(b,'x')
    panel(c,'C','Proxy eligibility','3/3 marker agreement; Visium 2/16')
    marker=read(SD/'S3/S3A_marker_freeze_decision.tsv');marker=marker[marker['gene'].isin(['CD34','LYL1','CD1A'])]
    det=read(D/'a4_yayon/a4c_detection_eligibility.tsv')
    ze=det[det.gene=='ZEB2'].iloc[0]
    labels=[]; numerators=[]; denominators=[]
    for r in marker.itertuples():
        labels.append(r.gene);numerators.append(int(r.n_agree));denominators.append(int(r.n_measured))
    labels.append('Visium ZEB2');numerators.append(int(ze.n_sections_ge10pct));denominators.append(int(ze.n_sections_tested))
    for i,(lab,num,den) in enumerate(zip(labels,numerators,denominators)):
        frac=num/den
        c.hlines(i,0,frac,color=Z2 if i<3 else GREY,lw=2.2)
        c.scatter(frac,i,color=Z2 if i<3 else GREY,s=23,zorder=3)
        c.text(min(frac+.03,.83),i-.18,f'{num}/{den}',fontsize=6.6,color=DARK)
    c.set_xlim(0,1.12);c.set_ylim(-.5,3.6);c.set_yticks(range(4),labels,fontsize=6.5);c.invert_yaxis()
    c.set_xlabel('Eligible fraction');grid(c,'x')
    save(f,1)

def sc2():
    f=m.mmfig(180,195);g=GridSpec(3,2,figure=f);axes=[f.add_subplot(g[i,j]) for i in range(3) for j in range(2)]
    a,b,c,d,e,h=axes
    x=read(SD/'revision3_robustness/crossfit_subtype_summary.tsv').sort_values('median_primary')
    y=np.arange(len(x));a.plot(x.median_primary,y,'o',color=GREY,ms=4,label='Primary');a.plot(x.median_crossfit,y,'o',color=Z2,ms=4,label='Leave-subtype-out')
    for i,r in enumerate(x.itertuples()):a.plot([r.median_primary,r.median_crossfit],[i,i],color=LIGHT,lw=.7)
    a.set_yticks(y,[m.short_subtype(q) for q in x.subtype],fontsize=6.5);a.invert_yaxis();a.axvline(0,color=LIGHT,lw=.5)
    panel(a,'A','Leave-one-subtype-out residual','17 classes; target subtype excluded from spline fit');a.set_xlabel('Median R');a.legend(frameon=False,fontsize=6.5);grid(a,'x')
    x=read(SD/'revision3_robustness/alternative_coordinates_subtype_summary.tsv');x=x[x.subtype.isin(['BCL11B','ETP-like','TLX3'])]
    defs=list(x.definition.drop_duplicates());for_colors={'BCL11B':Z2,'ETP-like':ETP,'TLX3':Z1}
    for j,s in enumerate(['BCL11B','ETP-like','TLX3']):
        xx=x[x.subtype==s];b.plot(np.arange(len(xx)),xx.median_residual,marker='o',ms=3,lw=1,color=for_colors[s],label=s)
    b.set_xticks(np.arange(len(defs)),[q.replace('leave_','−')[:13] for q in defs],rotation=30,ha='right',fontsize=6.5)
    panel(b,'B','Alternative proxy definitions','BCL11B and TLX3 remain opposite poles');b.set_ylabel('Median R');b.legend(frameon=False,fontsize=6.5);grid(b)
    x=read(SD/'revision3_robustness/model_comparison.tsv')
    x=x[x.model.isin(['development','subtype','combined'])]
    xx=np.arange(len(x));w=.22
    for off,col,color,lab in [(-w,'raw_r2',GREY,'Raw'),(0,'adjusted_r2',ETP,'Adjusted'),(w,'cv_r2_mean',DARK,'CV')]:c.bar(xx+off,100*x[col],width=w,color=color,label=lab)
    c.set_xticks(xx,x.model,fontsize=6.5);c.set_ylabel('R² (%)');c.legend(frameon=False,fontsize=6.5,ncol=3)
    panel(c,'C','Model-complexity check','Raw, adjusted and repeated cross-validated R²');grid(c)
    x=read(SD/'revision4_targeted/spline_df_2_3_4.tsv')
    for s,col in [('BCL11B',Z2),('ETP-like',ETP),('TLX3',Z1)]:
        xx=x[x.subtype==s].sort_values('spline_df');d.plot(xx.spline_df,xx.median_residual,marker='o',color=col,lw=1.2,label=s)
    d.axhline(0,color=LIGHT,lw=.5);d.set_xticks([2,3,4]);d.set_xlabel('Natural spline df');d.set_ylabel('Median R');d.legend(frameon=False,fontsize=6.5)
    panel(d,'D','Prespecified spline flexibility','df = 2, 3, 4; no model selection');grid(d)
    x=read(D/'gate2_sensitivity/frozen_genes_batch_adjusted.tsv').dropna(subset=['beta_primary','beta_batch'])
    if len(x):e.scatter(x.beta_primary,x.beta_batch,s=10,color=Z2,alpha=.6);e.plot([-1,1],[-1,1],color=LIGHT,lw=.6)
    else:
        x=read(D/'gate2_sensitivity/batch_sensitivity_overall.tsv').iloc[0]
        e.axis('off');e.text(.07,.62,f'{int(x.n_same_sign)}/100 same sign',transform=e.transAxes,fontsize=11,color=Z2,weight='bold');e.text(.07,.39,f'ρ = {x.spearman_beta:.4f}',transform=e.transAxes,fontsize=9,color=DARK)
    panel(e,'E','Batch-adjusted frozen effects','Direction preserved for all 100 frozen genes');grid(e)
    p=read(SD/'F3/F3B-C_patient_level.tsv');bc=p[p.subtype=='BCL11B'];outside=int((bc.dev_outside_normal_range.astype(str).str.lower()=='true').sum())
    h.bar([0,1],[outside,len(bc)-outside],color=[Z2,LIGHT])
    h.set_xticks([0,1],['Outside normal D','Within domain'],fontsize=6.5);h.set_ylabel('BCL11B patients')
    panel(h,'F','Reference-domain boundary','14/18 BCL11B outside measured normal D');grid(h)
    save(f,2)

def sc3():
    f,g=fig(3,175);a=f.add_subplot(g[0,0]);b=f.add_subplot(g[0,1]);c=f.add_subplot(g[1,0]);d=f.add_subplot(g[1,1])
    p=read(D/'aieop120_official/aieop120_residual_with_official_labels.tsv')
    count=p.predicted_subtype.value_counts().sort_values();a.barh(np.arange(len(count)),count.values,color=GREY)
    a.set_yticks(np.arange(len(count)),[m.short_subtype(x) for x in count.index],fontsize=6.5);a.set_xlabel('Officially assigned patients')
    panel(a,'A','Official AIEOP120 labels','120 patients; 16 classes');grid(a,'x')
    q=read(D/'aieop120_official/aieop120_vs_s9.tsv');b.scatter(np.arange(len(q)),q.max_probability,s=19,color=ETP)
    b.set_xlabel('Fusion-anchored case');b.set_ylabel('Official max probability');b.set_ylim(0,1.05)
    panel(b,'B','Fusion-anchor classifier audit','33/33 exact label matches; no added probability cutoff');grid(b)
    x=read(D/'aieop120_official/aieop120_subtype_residual_medians.tsv');x=x[x.n_aieop>=5].sort_values('residual_median_aieop');yy=np.arange(len(x))
    c.hlines(yy,x.residual_median_lo,x.residual_median_hi,color=GREY,lw=1.1);c.scatter(x.residual_median_aieop,yy,c=[m.color_subtype(s) for s in x.predicted_subtype],s=18)
    c.set_yticks(yy,[m.short_subtype(s) for s in x.predicted_subtype],fontsize=6.5);c.invert_yaxis();c.axvline(0,color=LIGHT,lw=.5);c.set_xlabel('External median R and bootstrap 95% CI')
    panel(c,'C','Six supported subtype summaries','n ≥ 5 only; no 17-class rank claim');grid(c,'x')
    for s,col in [('BCL11B',Z2),('ETP-like',ETP),('TLX3',Z1)]:
        t=p[p.predicted_subtype==s];d.scatter(t.dev,t.residual,s=20,color=col,label=f'{s} n={len(t)}')
    d.axhline(0,color=LIGHT,lw=.5);d.set_xlabel('External developmental proxy D');d.set_ylabel('External within-cohort R');d.legend(frameon=False,fontsize=6.5)
    panel(d,'D','Patient-level external geometry','Selected poles; all 120 patients in Table S4');grid(d)
    save(f,3)

def sc4():
    f,g=fig(4,170);a=f.add_subplot(g[0,0]);b=f.add_subplot(g[0,1]);c=f.add_subplot(g[1,0]);d=f.add_subplot(g[1,1])
    x=read(D/'gate2_program/all_genes_residual_effects.tsv');a.scatter(x.beta_R,-np.log10(x['adj.P.Val'].clip(lower=1e-300)),s=2,color=LIGHT,alpha=.5,rasterized=True)
    fr=read(D/'gate2_program/frozen_ZEB_side50.tsv')
    for s,col in [('ZEB1-side50',Z1),('ZEB2-side50',Z2)]:
        t=fr[fr.side==s];a.scatter(t.beta_R,-np.log10(t['adj.P.Val'].clip(lower=1e-300)),s=8,color=col,label=s)
    a.axvline(0,color=LIGHT,lw=.5);a.set_xlabel('Residual effect βR');a.set_ylabel('−log10 BH FDR');a.legend(frameon=False,fontsize=6.5)
    panel(a,'A','Frozen selection audit','3,560 effect-threshold genes → 50 per sign');grid(a)
    x=read(D/'gate2_sensitivity/heldout_vs_original_overlap.tsv');b.barh(x.side,x.n_overlap,color=[Z1 if 'ZEB1' in s else Z2 for s in x.side],height=.5);b.set_xlim(0,55);b.set_xlabel('Genes retained in held-out rederivation')
    for i,r in enumerate(x.itertuples()):b.text(r.n_overlap+.5,i,f'{int(r.n_overlap)}/50',va='center',fontsize=6.5)
    panel(b,'B','Held-out gene-set overlap','Original versus leave BCL11B/ETP-out');grid(b,'x')
    x=read(D/'gate2_sensitivity/heldout_matched_cameraPR.tsv');x=x[x.program.isin(['ZEB1-side50','ZEB2-side50'])]
    c.barh(x.program,-np.log10(x.PValue),color=[Z1 if 'ZEB1' in s else Z2 for s in x.program],height=.5);c.set_xlabel('−log10 CAMERA P')
    panel(c,'C','Held-out pair-set enrichment','Same-cohort phenotype bridge');grid(c,'x')
    pc=read(D/'gate2_program/residual_t_primary_vs_pc15.tsv').dropna(subset=['t_primary','t_pc'])
    x=read(D/'gate2_program/residual_pc_sensitivity.tsv').iloc[0]
    bounds=[np.quantile(pc[q],[.01,.99]) for q in ['t_primary','t_pc']]
    d.hexbin(pc.t_primary,pc.t_pc,gridsize=45,mincnt=1,cmap='Blues',linewidths=0,rasterized=True)
    d.axhline(0,color=LIGHT,lw=.6);d.axvline(0,color=LIGHT,lw=.6)
    d.set_xlim(*bounds[0]);d.set_ylim(*bounds[1]);d.set_xlabel('Primary gene t');d.set_ylabel('PC1–5 adjusted gene t')
    panel(d,'D','Technical-axis boundary',f'24,356 genes · rank ρ={x.spearman_t_primary_vs_pc:.2f} · central 98% axes');grid(d)
    save(f,4)

def sc5():
    f,g=fig(5,185);a=f.add_subplot(g[0,0]);b=f.add_subplot(g[0,1]);c=f.add_subplot(g[1,0]);d=f.add_subplot(g[1,1])
    p=read(SD/'revision3_robustness/bcl11b_etp_matched_pairs.tsv');a.hist(p.distance,bins=8,color=ETP,edgecolor='white');a.set_xlabel('|Δ developmental proxy|');a.set_ylabel('Matched pairs')
    panel(a,'A','Match quality','18 one-to-one pairs; median distance .016');grid(a)
    b.scatter(p.BCL11B_residual,p.ETP_residual,s=16,color=Z2);ll=[min(p.BCL11B_residual.min(),p.ETP_residual.min()),max(p.BCL11B_residual.max(),p.ETP_residual.max())];b.plot(ll,ll,color=GREY,lw=.7,ls='--')
    b.set_xlabel('BCL11B patient R');b.set_ylabel('Matched ETP-like R');panel(b,'B','Pair-level residual audit','Difference estimated with patient-pair bootstrap');grid(b)
    h=read(SD/'revision4_targeted/matched_BCL11B_ETP_hallmark_camera.tsv').sort_values('FDR').head(7).sort_values('FDR',ascending=False)
    c.barh(np.arange(len(h)),-np.log10(h.FDR),color=[ETP if q=='Up' else GREY for q in h.Direction]);c.set_yticks(np.arange(len(h)),[q.replace('HALLMARK_','').replace('_',' ').title()[:17] for q in h.pathway],fontsize=6.5)
    c.set_xlabel('−log10 Hallmark FDR');panel(c,'C','Matched transcriptional context','Directional enrichment; no pathway mechanism inferred');grid(c,'x')
    x=read(SD/'F4/F4C_GSE162280_cases.tsv').sort_values('log2_ZEB1_over_ZEB2');yy=np.arange(len(x))
    d.hlines(yy,0,x.log2_ZEB1_over_ZEB2,color=LIGHT,lw=1.3)
    d.scatter(x.log2_ZEB1_over_ZEB2,yy,color=Z2,s=19)
    d.set_yticks(yy,x['case'],fontsize=6.5);d.invert_yaxis();d.axvline(0,color=GREY,lw=.5)
    d.set_xlabel('log2(ZEB1 CPM / ZEB2 CPM)');panel(d,'D','Complete 12-case lesion audit','Partner, phenotype and ratio; no ordinary T-ALL control');grid(d,'x')
    save(f,5)

def sc6():
    f,g=fig(6,180);a=f.add_subplot(g[0,0]);b=f.add_subplot(g[0,1]);c=f.add_subplot(g[1,0]);d=f.add_subplot(g[1,1])
    cohorts=[('GSE146901',read(D/'gate3_validation/gse146901_program_scores.tsv'),'ZEB1_side50','ZEB2_side50','group'),('GSE234608',read(D/'a5_gse234608/a5_patient_frozen_scores.tsv'),'ZEB1_side50','ZEB2_side50','group'),('GSE243914',read(D/'a7b_malignant/a7b05_sample_scores.tsv'),'ZEB1_side50','ZEB2_side50','group')]
    for ax,letter,(name,x,z1,z2,grp) in zip([a,b,c],'ABC',cohorts):
        for i,(group,t) in enumerate(x.groupby(grp)):
            ax.scatter(t[z1],t[z2],s=17,color=[Z2,ETP,GREY][i],alpha=.8,label=f'{group} n={len(t)}')
        ax.axhline(0,color=LIGHT,lw=.5);ax.axvline(0,color=LIGHT,lw=.5);ax.set_xlabel('Frozen ZEB1-side score');ax.set_ylabel('Frozen ZEB2-side score');ax.legend(frameon=False,fontsize=6.5)
        panel(ax,letter,f'{name} patient audit','Frozen scores without gene reselection');grid(ax)
    summ=read(D/'a65_external/a65_summary.tsv');summ=summ[summ['set']=='all_recovered']
    others=[('GSE243914',.918,.77),('GSE248287',.777,.62),('GSE280250',.987,.86)]
    vals=[(r.cohort,r.expected_frac,r.rho_delta_vs_neg_beta) for r in summ.itertuples()]+others
    yy=np.arange(len(vals));d.barh(yy,[v[1]*100 for v in vals],color=[GREY,GREY,Z2,ETP,DARK]);d.set_yticks(yy,[v[0] for v in vals],fontsize=6.5);d.invert_yaxis();d.set_xlim(0,112)
    for i,v in enumerate(vals):d.text(v[1]*100+1,i,f'ρ {v[2]:.2f}',fontsize=6.5,va='center')
    d.set_xlabel('Expected-direction genes (%)');panel(d,'D','Gene-level preservation','Five external assay contexts; recovered n in Table S7');grid(d,'x')
    save(f,6)

def sc7():
    f,g=fig(7,180);a=f.add_subplot(g[0,0]);b=f.add_subplot(g[0,1]);c=f.add_subplot(g[1,0]);d=f.add_subplot(g[1,1])
    x=read(D/'a7_lineage/scores/a7a_n_cells.tsv');
    if {'donor','official','n_cells'}.issubset(x.columns):
        x=x.rename(columns={'official':'broad'})
        x=x[x.broad.isin(['HSC/MPP/HSPC','lymphoid progenitor/CLP','T lineage','Myeloid'])]
        z=x.pivot_table(index='donor',columns='broad',values='n_cells',aggfunc='sum').fillna(0)
        im=a.imshow(np.log10(z+1),aspect='auto',cmap='Blues',vmin=0);a.set_xticks(range(len(z.columns)),[q.replace('lymphoid progenitor/','').replace('HSC/MPP/','') for q in z.columns],rotation=25,ha='right',fontsize=6.5);a.set_yticks(range(len(z.index)),z.index,fontsize=6.5)
        cb=f.colorbar(im,ax=a,fraction=.035,pad=.02);cb.set_label('log10(cells+1)',fontsize=6.5);cb.ax.tick_params(labelsize=6.5)
    panel(a,'A','Marrow donor × lineage eligibility','Cell counts; ≥20 primary, ≥10 sensitivity')
    order=['HSC/MPP/HSPC','lymphoid progenitor/CLP','T lineage','Myeloid']
    for i,(name,path) in enumerate([('≥20',D/'a7_lineage/scores/a7a_unit_scores_min20.tsv'),('≥10',D/'a7_lineage/scores/a7a_unit_scores_min10.tsv')]):
        t=read(path).groupby('lineage')['ZEB2_side50'].median().reindex(order)
        b.plot(np.arange(4)+i*.09,t.values,marker='o',ms=4,lw=1,color=Z2 if i==0 else GREY,label=name)
    b.set_xticks(np.arange(4),['HSPC','CLP','T','Myeloid'],rotation=20,ha='right');b.set_ylabel('Median ZEB2-side score');b.legend(frameon=False,fontsize=6.5)
    panel(b,'B','Cell-threshold sensitivity','Myeloid and HSPC directions retained');grid(b)
    x=read(D/'a7_lineage/a7_0_patient_scores.tsv');
    for i,(spec,t) in enumerate(x.groupby('specimen')):
        vals=t.ZEB2_side50.dropna().to_numpy();c.scatter(np.full(len(vals),i)+m.jitter(len(vals),.09,i+7),vals,s=7,color=Z2 if 'blood' in spec.lower() else GREY,alpha=.45,rasterized=True)
        c.plot([i-.22,i+.22],[np.median(vals)]*2,color=DARK,lw=1.7)
    c.set_xticks(range(len(x.specimen.unique())),sorted(x.specimen.unique()));c.set_ylabel('Frozen ZEB2-side score');panel(c,'C','BM/PB source comparison','Specimen P=.13; patient scores shown');grid(c)
    x=read(D/'a7c_decompose/a7c_zeb2_split.tsv');x=x.dropna(subset=['delta_Myeloid_minus_T','beta_R'])
    for group,t in x.groupby('subset'):
        d.scatter(t.delta_Myeloid_minus_T,t.beta_R,s=18,color=Z2 if 'myeloid_concordant'==group else Z1,label=f'{group} n={len(t)}')
    d.axhline(0,color=LIGHT,lw=.5);d.axvline(0,color=LIGHT,lw=.5);d.set_xlabel('Normal myeloid−T Δ');d.set_ylabel('Leukemia βR');d.legend(frameon=False,fontsize=6.5)
    panel(d,'D','45 + 3 gene decomposition','Three non-myeloid genes are a slice, not a new program');grid(d)
    save(f,7)

def sc8():
    f,g=fig(8,180);a=f.add_subplot(g[0,0]);b=f.add_subplot(g[0,1]);c=f.add_subplot(g[1,0]);d=f.add_subplot(g[1,1])
    x=read(D/'a7b_malignant/a7b_patient_scores.tsv').sort_values('n_cells');a.barh(np.arange(len(x)),x.n_cells,color=ETP);a.set_yticks(np.arange(len(x)),x.patient,fontsize=6.5);a.set_xscale('log');a.set_xticks([100,1000,10000],['100','1,000','10,000']);a.set_xlabel('Author malignant cells per patient')
    panel(a,'A','Malignant-cell pseudobulk eligibility','15 diagnosis patients · author Celltypes_all label');grid(a,'x')
    b.scatter(x.B,x.ZEB1_side50,s=17,color=Z1,label='ZEB1-side');b.scatter(x.B,x.ZEB2_side50,s=17,color=Z2,label='ZEB2-side');b.axvline(0,color=LIGHT,lw=.5);b.set_xlabel('Malignant pseudobulk B');b.set_ylabel('Frozen score');b.legend(frameon=False,fontsize=6.5)
    panel(b,'B','Malignant score boundary','ZEB1-side ρ+.73; ZEB2-side ρ−.21, CI crosses 0');grid(b)
    x=read(SD/'S7/S7A_eligible_units.tsv');x=x[x.stage.isin(['Eligible Day0 samples','Day0 paired-state samples','Eligible Day28 samples'])]
    order=['Eligible Day0 samples','Day0 paired-state samples','Eligible Day28 samples'];x=x.set_index('stage').loc[order].reset_index()
    c.barh(np.arange(len(x)),x.n,color=[GREY,ETP,Z2]);c.set_yticks(np.arange(len(x)),['Day0 54','ZBTB16 pairs 41','Day28 pairs 12'],fontsize=6.5);c.set_xlabel('Eligible patients / samples')
    panel(c,'C','Lim denominator audit','Day0, paired states and Day28 are distinct denominators');grid(c,'x')
    x=read(SD/'S7/S7D_Day0_Day28_pairs.tsv');yy=np.arange(len(x))
    d.plot([0,0],[0,0],alpha=0)
    for i,r in enumerate(x.itertuples()):d.plot([0,1],[r.day0_balance,r.day28_balance],color=GREY,lw=.7);d.scatter([0,1],[r.day0_balance,r.day28_balance],color=[ETP,Z2],s=12)
    d.set_xticks([0,1],['Day 0','Day 28']);d.set_ylabel('Patient pseudobulk balance')
    panel(d,'D','Day28 negative boundary','12 paired patients; no reproducible state selection');grid(d)
    save(f,8)

def sc9():
    f,g=fig(9,175);a=f.add_subplot(g[0,0]);b=f.add_subplot(g[0,1]);c=f.add_subplot(g[1,0]);d=f.add_subplot(g[1,1])
    p=read(D/'a8_adult/a8_patient_scores.tsv');a.hist(p.age,bins=np.arange(10,81,10),color=ETP,edgecolor='white');a.axvline(18,color=Z2,lw=.9,ls='--');a.set_xlabel('Age, years');a.set_ylabel('Diagnosis samples')
    panel(a,'A','Author adult-series age spectrum','79 cases, seven younger than 18');grid(a)
    g=read(D/'a8_adult/a8_genes.tsv');g=g.dropna(subset=['beta_adult','beta_R']);b.scatter(g.beta_R,g.beta_adult,c=[Z1 if 'ZEB1' in s else Z2 for s in g.side],s=16);b.axhline(0,color=LIGHT,lw=.5);b.axvline(0,color=LIGHT,lw=.5)
    b.set_xlabel('Discovery βR');b.set_ylabel('Adult βR');panel(b,'B','Adult constituent genes','86 recovered, 78 tested, 98.7% expected');grid(b)
    x=read(D/'a8_adult/a8_2_merged.tsv');col='T-ALL main-cluster high-confidence';vc=x[col].value_counts(dropna=True).head(10).sort_values()
    c.barh(np.arange(len(vc)),vc.values,color=GREY);c.set_yticks(np.arange(len(vc)),[str(s).split(' (')[0] for s in vc.index],fontsize=6.5);c.set_xlabel('High-confidence patients')
    panel(c,'C','Official ALLCatchR2 calls','Subtype context only, not the A8-1 endpoint');grid(c,'x')
    ff=[('BCL11B',7,-2.34),('ETP-like',47,-.06),('TAL1 DP',7,.32),('TLX3 imm.',2,1.54)]
    for i,(name,n,v) in enumerate(ff):d.scatter(v,i,color=m.color_subtype(name),s=22);d.text(v+.08,i,f'n={n}',va='center',fontsize=6.5)
    d.set_yticks(range(4),[q[0] for q in ff]);d.invert_yaxis();d.axvline(0,color=LIGHT,lw=.5);d.set_xlabel('Adult median R')
    panel(d,'D','Prespecified adult poles','Only three uniquely mapped n≥5 classes enter rank');grid(d,'x')
    save(f,9)

def sc10():
    f,g=fig(10,185);a=f.add_subplot(g[0,0]);b=f.add_subplot(g[0,1]);c=f.add_subplot(g[1,0]);d=f.add_subplot(g[1,1])
    x=read(D/'a3_gse165209/gse173432/a3_cameraPR_frozen_side50.tsv')
    y=np.arange(len(x));direction=np.where(x.Direction.astype(str)=='Down',-1,1)
    a.barh(y,direction*(-np.log10(x['PValue'].clip(lower=1e-300))),color=[Z1 if 'ZEB1' in q else Z2 for q in x.program]);a.set_yticks(y,x.program,fontsize=6.5);a.axvline(0,color=LIGHT,lw=.5);a.set_xlabel('Signed −log10 CAMERA P (left = Down)')
    panel(a,'A','Acute BCL11B overexpression','Frozen ZEB2-side direction not induced in CD34+ cells');grid(a,'x')
    x=read(SD/'S10/S10A-C_HiChIP_sample_contacts.tsv');x=x.drop_duplicates('sample')
    for group, color, label in [(True,Z2,'ETP'),(False,GREY,'non-ETP')]:
        mask=x.group.astype(str).str.contains('ETP',case=False) & ~x.group.astype(str).str.contains('non',case=False)
        t=x[mask if group else ~mask]
        b.scatter(t.ZEB2_over_ZEB1,t.n_contacts/1e6,color=color,s=25,label=label)
    b.set_xlabel('ZEB2/ZEB1 contact ratio');b.set_ylabel('Contacts (million)');panel(b,'B','Exploratory HiChIP polarity','3 ETP versus 4 non-ETP; exact two-sided P=.057');b.legend(frameon=False,fontsize=6.5,loc='upper right');grid(b)
    x=read(SD/'S10/S10D_scATAC_peak_sets.tsv');c.bar(np.arange(len(x)),x.ZEB2_per_10k_peaks,color=[Z2 if 'ETP' in str(q) and 'non' not in str(q) else GREY for q in x.group]);c.set_xticks(np.arange(len(x)),x['sample'],rotation=40,ha='right',fontsize=6.5);c.set_ylabel('ZEB2 peaks per 10k')
    panel(c,'C','Cross-assay negative boundary','Seven scATAC peak sets do not reproduce polarity');grid(c)
    x=read(SD/'S6/S6A_all_endpoint_models.tsv')
    endpoint_order=['Induction failure','MRD >=0.1%']
    model_order=[('Unadjusted',Z2,-.08),('Subtype-adjusted',GREY,.08)]
    for i,ep in enumerate(endpoint_order):
        t=x[x.endpoint==ep].set_index('adjustment')
        for label,color,offset in model_order:
            r=t.loc[label]
            d.scatter(r.estimate,i+offset,color=color,s=20,label=label if i==0 else None)
            d.plot([r.ci_low,r.ci_high],[i+offset]*2,color=color,lw=1)
    d.set_yticks(range(len(endpoint_order)),endpoint_order,fontsize=6.5)
    d.legend(frameon=False,fontsize=6.5,loc='center right')
    d.axvline(1,color=LIGHT,lw=.6);d.set_xlabel('OR per unit balance (95% CI)')
    panel(d,'D','Clinical attenuation','Subtype adjustment; sparse Cox warning');grid(d,'x')
    save(f,10)

def main():
    for n,fn in enumerate([sc1,sc2,sc3,sc4,sc5,sc6,sc7,sc8,sc9,sc10],1):
        print('S',n,flush=True);fn()
    print('Exported',OUT)

if __name__=='__main__':main()

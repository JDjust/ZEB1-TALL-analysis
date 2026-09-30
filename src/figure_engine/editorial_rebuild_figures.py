"""Rebuild the seven evidence-led plates from frozen source tables.

Only visualization joins and descriptive summaries are performed here. No
gene selection, model fitting beyond drawing the already saved spline curve,
or hypothesis testing is done in this script.
"""
from pathlib import Path
import os
import math
import numpy as np
import pandas as pd
from scipy.stats import spearmanr
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from matplotlib.lines import Line2D
from matplotlib.patches import FancyArrowPatch, Rectangle, Circle

ROOT = Path(os.environ.get('ZEB1_DATA_ROOT',str(Path(__file__).resolve().parents[2]/'data/processed'))) / 'total/rebuild_2026'
DATA = ROOT / 'data'
SD = DATA / 'source_data_rebuilt'
D = DATA / 'deepen_2026'
OUT = Path(__file__).resolve().parents[2] / 'results/figures'

Z2 = '#23658A'      # negative residual / ZEB2-side
Z1 = '#C46D35'      # positive residual / ZEB1-side
ETP = '#638D8A'
BCL = '#275B78'
TLX = '#B45E27'
GREY = '#7C858B'
LIGHT = '#D6DBDE'
DARK = '#172B39'
SOFT = '#F4F6F7'

plt.rcParams.update({
    'font.family': 'Arial', 'font.size': 7.4, 'axes.titlesize': 8.2,
    'axes.labelsize': 7.5, 'xtick.labelsize': 6.8, 'ytick.labelsize': 6.8,
    'axes.linewidth': .55, 'pdf.fonttype': 42, 'ps.fonttype': 42,
    'savefig.facecolor': 'white', 'figure.facecolor': 'white',
    'axes.spines.top': False, 'axes.spines.right': False,
})

def read(path):
    return pd.read_csv(path, sep='\t', low_memory=False)

def mmfig(width, height):
    return plt.figure(figsize=(width / 25.4, height / 25.4), constrained_layout=False)

def panel(ax, letter, title, subtitle=None):
    ax.text(-.02, 1.10, letter, transform=ax.transAxes, fontsize=10.5,
            weight='bold', color=DARK, va='bottom', clip_on=False)
    ax.text(.07, 1.10, title, transform=ax.transAxes, fontsize=8.2,
            weight='bold', color=DARK, va='bottom', clip_on=False)
    if subtitle:
        ax.text(.07, 1.02, subtitle, transform=ax.transAxes, fontsize=6.7,
                color=GREY, va='bottom', clip_on=False)
    ax.tick_params(length=2.2, width=.45, pad=2)
    ax.spines['left'].set_color('#66747d')
    ax.spines['bottom'].set_color('#66747d')

def save(fig, num):
    fig.subplots_adjust(left=.14 if num==2 else .085, right=.975, top=.94, bottom=.065,
                        wspace=.42 if num==2 else .32, hspace=.43 if num==5 else .62)
    for label in fig.findobj(match=matplotlib.text.Text):
        if label.get_text().strip() and label.get_fontsize()<6.5:
            label.set_fontsize(6.5)
    fig.savefig(OUT / f'Figure{num}.pdf', bbox_inches=None)
    fig.savefig(OUT / f'Figure{num}.png', dpi=300, bbox_inches=None)
    plt.close(fig)

def color_subtype(s):
    q=str(s).lower()
    if 'bcl11b' in q: return BCL
    if 'etp' in q: return ETP
    if 'tlx3' in q: return TLX
    return '#AAB3B8'

def short_subtype(s):
    s=str(s)
    for old,new in [('<U+03B3><U+03B4>','gd'),('<U+03B1><U+03B2>','ab'),('<U+03B3>','g'),('<U+03B4>','d'),('<U+03B1>','a'),('<U+03B2>','b')]:
        s=s.replace(old,new)
    s=s.replace('LMO2 gd-like','LMO2 gd').replace('TAL1 ab-like','TAL1 ab').replace('TAL1 DP-like','TAL1 DP')
    return s

def jitter(n, scale=.11, seed=2):
    return np.random.default_rng(seed).normal(0, scale, n)

def clean_grid(ax, axis='y'):
    ax.grid(axis=axis, color='#E5E9EB', linewidth=.38, zorder=0)
    ax.set_axisbelow(True)

def fig1():
    fig=mmfig(180,190); g=GridSpec(3, 6, figure=fig, height_ratios=[.9,1.2,1.1])
    a=fig.add_subplot(g[0,:2]); b=fig.add_subplot(g[0,2:]);
    c=fig.add_subplot(g[1,:4]); d=fig.add_subplot(g[1,4:]);
    e=fig.add_subplot(g[2,:3]); f=fig.add_subplot(g[2,3:])
    for ax in (a,): ax.axis('off')
    panel(a,'A','Thymic compartments')
    steps=[('ENTRY','ETP / early DN',Z2),('CORTEX','DP selection',Z1),('MEDULLA','SP return',GREY)]
    for i,(lab,sub,col) in enumerate(steps):
        y=.81-i*.31
        a.add_patch(Rectangle((.04,y-.10),.85,.19,facecolor=SOFT,edgecolor='none',transform=a.transAxes))
        a.add_patch(Rectangle((.04,y-.10),.025,.19,facecolor=col,edgecolor='none',transform=a.transAxes))
        a.text(.10,y+.015,lab,transform=a.transAxes,weight='bold',fontsize=8,color=DARK)
        a.text(.10,y-.065,sub,transform=a.transAxes,fontsize=6.8,color=GREY)
        if i<2: a.add_patch(FancyArrowPatch((.47,y-.11),(.47,y-.20),transform=a.transAxes,arrowstyle='-|>',mutation_scale=7,color=GREY,lw=.8))
    tr=read(SD/'F1/F1D_stage_tracks.tsv')
    bt=tr[tr['track'].astype(str).str.contains('balance',case=False,na=False)]
    if bt.empty: bt=tr[tr['track'].isin(['ZEB1','ZEB2'])].copy()
    atlas_colors={'GSE142522':'#465D6B','GSE195812':'#9A6A4C','GSE206710':'#6F8F83','Park/HTA':'#A99055'}
    for k,x in bt.groupby(['dataset','track']):
        x=x.sort_values('order')
        col=atlas_colors.get(k[0],GREY)
        b.plot(x['order'],x['value'],marker='o',ms=2.8,lw=1.25,label=k[0],color=col,alpha=.9)
    panel(b,'B','Four normal-thymus references','Aligned stage trajectories; source-specific scales')
    b.set_xlabel('Author-defined developmental order'); b.set_ylabel('Displayed expression / balance')
    b.legend(loc='lower right',fontsize=6.2,frameon=False,ncol=2)
    clean_grid(b)
    st=read(D/'a4_yayon/a4a_stage_balance.tsv')
    stage_order=['T_ETP','T_DN(early)','T_DP(Q)-early','T_DP(Q)-CD99','T_DP(Q)','T_DP(P)','T_CD4','T_CD8']
    st=st.set_index('stage').loc[stage_order].reset_index()
    st['rank']=np.arange(len(st)); st['pole']=st['pole'].astype(str)
    col=[Z2 if s in stage_order[:2] else Z1 if s in stage_order[2:6] else GREY for s in st['stage']]
    c.axvspan(-.5,1.5,color='#E7EFF2',zorder=0)
    c.axvspan(1.5,5.5,color='#F7EDE5',zorder=0)
    c.axvspan(5.5,7.5,color='#F2F3F4',zorder=0)
    c.scatter(st['rank'],st['balance_median'],c=col,s=25,zorder=2,edgecolors='white',linewidths=.3)
    panel(c,'C','Yayon author-state refinement','Eight prespecified poles among 33 eligible states')
    c.set_ylabel('Median ZEB balance'); c.set_xlabel('Author stage order')
    c.set_xticks(st['rank'],st['stage'].str.replace('T_','',regex=False),rotation=28,ha='right',fontsize=6.3)
    c.annotate('ETP: one donor',xy=(0,st.loc[0,'balance_median']),xytext=(.5,-2.55),fontsize=6.1,color=Z2,arrowprops=dict(arrowstyle='-',lw=.5,color=Z2))
    c.axhline(0,color=LIGHT,lw=.6);clean_grid(c)
    dd=read(D/'a4_yayon/a4a_donor_level_balance.tsv')
    for i,r in dd.iterrows():
        d.plot([0,1],[r['early'],r['cortical_DP']],color=Z1 if r['DP_minus_early']>0 else GREY,lw=1.1,alpha=.9)
        d.scatter([0,1],[r['early'],r['cortical_DP']],color=[Z2,Z1],s=17,zorder=3)
    panel(d,'D','Donor-paired early → DP','4/5 concordant; P = 0.063')
    d.set_xticks([0,1],['Early','Cortical DP']);d.set_xlim(-.15,1.15);d.set_ylabel('ZEB balance');clean_grid(d)
    cm=read(D/'a4_yayon/a4b_tstate_cma_by_donor.tsv')
    selected=['T_ETP','T_DP(Q)','T_CD4','T_CD8']
    cm=cm[cm.state.isin(selected)]
    for i,s in enumerate(selected):
        x=cm.loc[cm.state==s,'cma_weighted_median'].dropna().to_numpy()
        e.scatter(x,np.repeat(i,len(x))+jitter(len(x),.07,i+1),color=Z2 if i==0 else Z1 if i==1 else GREY,s=20,alpha=.85)
        if len(x): e.plot([np.median(x),np.median(x)],[i-.2,i+.2],color=DARK,lw=2)
    panel(e,'E','Spatial CMA localization','Six donor maps; author cell2location states')
    e.set_yticks(range(len(selected)),selected);e.set_xlabel('Corticomedullary axis');e.invert_yaxis();clean_grid(e,'x')
    topo=read(D/'a4_yayon/a4ab_stage_balance_and_cma.tsv')
    pc={'early_DN':Z2,'cortical_DP':Z1,'SP':GREY,'early':Z2,'DP':Z1}
    for pole,x in topo.groupby('pole'):
        f.scatter(x['cma_weighted_median'],x['balance_median'],s=np.clip(x['n_cells']/250,15,70),color=pc.get(pole,GREY),label=pole,alpha=.78,edgecolors='white',linewidths=.35)
    panel(f,'F','Development × spatial topology','Stage medians; circle area reflects sampled cells')
    f.set_xlabel('Corticomedullary axis');f.set_ylabel('Stage median ZEB balance');f.axhline(0,color=LIGHT,lw=.5);f.legend(frameon=False,fontsize=6,loc='best');clean_grid(f)
    save(fig,1)

def fig2():
    fig=mmfig(180,205);g=GridSpec(3,6,figure=fig,height_ratios=[.82,1.28,1.1])
    a=fig.add_subplot(g[0,:2]);b=fig.add_subplot(g[0,2:]);c=fig.add_subplot(g[1,:2]);d=fig.add_subplot(g[1,2:]);e=fig.add_subplot(g[2,:3]);f=fig.add_subplot(g[2,3:])
    st=read(SD/'F3/F3D-F_subtype_statistics.tsv').sort_values('residual_within_median')
    names=st.subtype.astype(str).tolist(); display=[short_subtype(s) for s in names]; yy=np.arange(len(st))
    a.barh(yy,st['n'],color=[color_subtype(s) for s in names],height=.63)
    a.set_yticks(yy,['']*len(names));a.invert_yaxis();a.set_xlabel('Patients; order as D/E');
    panel(a,'A','Discovery cohort','1,309 diagnostic patients · 17 subtypes')
    for y,n in zip(yy,st['n']):a.text(n+.5,y,str(int(n)),va='center',fontsize=5.2,color=DARK)
    a.set_xlim(0,max(st['n'])*1.16);clean_grid(a,'x')
    pt=read(SD/'F3/F3B-C_patient_level.tsv')
    b.scatter(pt['dev'],pt['balance'],s=3,c='#BEC6CA',alpha=.55,rasterized=True)
    for s in ['BCL11B','ETP-like','TLX3']:
        x=pt[pt.subtype==s]
        b.scatter(x['dev'],x['balance'],s=11,color=color_subtype(s),label=f'{s} (n={len(x)})',alpha=.8,rasterized=True)
    curve=read(SD/'F3/F3B_within_cohort_curve.tsv')
    xc='dev' if 'dev' in curve else curve.columns[0]
    yc='expected_within' if 'expected_within' in curve else curve.columns[1]
    b.plot(curve[xc],curve[yc],color=DARK,lw=1.5)
    panel(b,'B','Balance versus developmental proxy','Within-cohort natural spline, df = 3')
    b.set_xlabel('Developmental-expression proxy D');b.set_ylabel('ZEB balance B');b.legend(frameon=False,fontsize=6,loc='lower right');clean_grid(b)
    var=read(SD/'F2/F2F_model_variance.tsv')
    c.barh(np.arange(len(var)),var['r2']*100,color=[GREY,Z2,ETP,DARK],height=.55)
    c.set_yticks(np.arange(len(var)),['Immunophen.','Development','Subtype','D + subtype'],fontsize=6.4);c.invert_yaxis()
    for i,v in enumerate(var['r2']*100):
        c.text(v-1.0,i,f'{v:.1f}%',va='center',ha='right',fontsize=6.5,
               color=DARK if i==0 else 'white',weight='bold')
    c.set_xlim(0,48);c.set_xlabel('Variance explained, raw R² (%)')
    panel(c,'C','Complementary variance','Subtype adds 13.9 percentage points over D')
    clean_grid(c,'x')
    rng=np.random.default_rng(71)
    for i,s in enumerate(names):
        x=pt.loc[pt.subtype==s,'residual_within'].dropna().to_numpy()
        d.scatter(x,np.full(len(x),i)+rng.normal(0,.09,len(x)),s=2.5,color=color_subtype(s),alpha=.50,rasterized=True)
        d.plot([np.median(x),np.median(x)],[i-.25,i+.25],color=DARK,lw=1.3)
    d.set_yticks(yy,display,fontsize=6.5);d.invert_yaxis();d.axvline(0,color=GREY,lw=.6)
    d.set_xlabel('Patient within-cohort residual R');panel(d,'D','All 17 subtype residual distributions','Points are patients; dark ticks are medians');clean_grid(d,'x')
    e.hlines(yy,st['residual_within_ci_low'],st['residual_within_ci_high'],color=[color_subtype(s) for s in names],lw=1.5)
    e.scatter(st['residual_within_median'],yy,c=[color_subtype(s) for s in names],s=16,zorder=3)
    e.set_yticks(yy,display,fontsize=6.4);e.invert_yaxis();e.axvline(0,color=GREY,lw=.6);e.set_xlabel('Median R (bootstrap 95% CI)')
    panel(e,'E','Subtype poles','BCL11B −1.76 · ETP-like −0.23 · TLX3 +0.64');clean_grid(e,'x')
    ex=read(D/'aieop120_official/aieop120_subtype_residual_medians.tsv')
    ex=ex[ex['n_aieop']>=5].sort_values('residual_median_aieop')
    fy=np.arange(len(ex))
    f.hlines(fy,ex['residual_median_lo'],ex['residual_median_hi'],color=[color_subtype(s) for s in ex['predicted_subtype']],lw=1.3)
    f.scatter(ex['residual_median_aieop'],fy,c=[color_subtype(s) for s in ex['predicted_subtype']],s=20,zorder=3)
    ext_labels=[]
    for subtype in ex.predicted_subtype:
        label=short_subtype(subtype)
        if label.startswith('TAL1') and 'like' in label: label='TAL1 ab'
        if label=='ETP-like': label='ETP'
        ext_labels.append(label)
    f.set_yticks(fy,ext_labels,fontsize=6.2);f.invert_yaxis();f.axvline(0,color=GREY,lw=.6)
    f.set_xlabel('AIEOP within-cohort median R');panel(f,'F','Moderate external architecture support','120 cases · official labels · six classes with n ≥ 5');clean_grid(f,'x')
    save(fig,2)

def fig3():
    fig=mmfig(180,210);g=GridSpec(3,6,figure=fig,height_ratios=[.8,1.0,1.27])
    a=fig.add_subplot(g[0,:2]);b=fig.add_subplot(g[0,2:]);c=fig.add_subplot(g[1,:3]);d=fig.add_subplot(g[1,3:]);e=fig.add_subplot(g[2,:4]);f=fig.add_subplot(g[2,4:])
    a.axis('off');panel(a,'A','Frozen program derivation')
    lines=[r'$E_g \sim ns(D,3) + subtype + R$',r'$FDR < 0.01\;\;\; |\beta_R| \geq 0.25$',r'Top 50 by partial R2 per sign',r'100 genes frozen before validation']
    for i,s in enumerate(lines):
        a.text(.04,.81-i*.19,s,transform=a.transAxes,fontsize=9.4 if i<2 else 7.7,color=DARK if i<3 else Z1,weight='bold' if i==3 else 'normal')
    eff=read(D/'gate2_program/all_genes_residual_effects.tsv');fr=read(D/'gate2_program/frozen_ZEB_side50.tsv')
    b.scatter(eff['beta_R'],eff['partial_R2'],s=2,color='#CCD2D6',alpha=.5,rasterized=True)
    for side,col in [('ZEB1-side50',Z1),('ZEB2-side50',Z2)]:
        x=fr[fr.side==side];b.scatter(x['beta_R'],x['partial_R2'],s=10,color=col,alpha=.88,rasterized=True,label=side)
    b.axvline(0,color=LIGHT,lw=.6);b.set_xlabel('Frozen leukemia residual effect βR');b.set_ylabel('Partial R²');b.legend(frameon=False,fontsize=6)
    panel(b,'B','Genome-wide residual effects','100 frozen genes emphasized');clean_grid(b)
    ps=read(D/'a7_lineage/a7_0_patient_scores.tsv');rr=read(SD/'F3/F3B-C_patient_level.tsv')[['sample_id','residual_within','subtype']]
    both=ps.merge(rr,on='sample_id',validate='one_to_one');assert len(both)==1309
    for side,col in [('ZEB1_side50',Z1),('ZEB2_side50',Z2)]:
        c.scatter(both.residual_within,both[side],s=3,alpha=.30,color=col,rasterized=True,label=side)
        x=both[['residual_within',side]].sort_values('residual_within')
        roll=x[side].rolling(125,min_periods=50,center=True).median()
        c.plot(x.residual_within,roll,color=col,lw=1.5)
    c.axvline(0,color=LIGHT,lw=.6);c.set_xlabel('Patient R');c.set_ylabel('Frozen side score');c.legend(frameon=False,fontsize=6)
    panel(c,'C','Continuous patient architecture','1,309 discovery patients; frozen scores');clean_grid(c)
    no=read(D/'a6_yayon/a6_1_donor_program.tsv')
    for side,col in [('ZEB1',Z1),('ZEB2',Z2)]:
        cols=[f'early_{side}_side50',f'DP_{side}_side50',f'SP_{side}_side50']
        for _,row in no.iterrows(): d.plot([0,1,2],[row[q] for q in cols],color=col,alpha=.18,lw=.7)
        d.plot([0,1,2],[no[q].median() for q in cols],color=col,lw=2,marker='o',ms=3,label=f'{side}-side')
    d.set_xticks([0,1,2],['Early','DP','SP']);d.legend(frameon=False,fontsize=6)
    panel(d,'D','Normal thymus programs','Five donor-level trajectories');clean_grid(d)
    q=read(D/'a6_yayon/a6_2_gene_effects.tsv')
    for cls,marker,alpha in [('concordant','o',.92),('discordant','x',.85),('neutral','.',.6),('weak','.',.6)]:
        x=q[q['class'].astype(str).str.lower()==cls]
        if len(x):
            e.scatter(x.delta_normal,x.beta_R,s=25 if marker!='.' else 11,
                      c=[Z1 if str(side).startswith('ZEB1') else Z2 for side in x.side],
                      marker=marker,alpha=alpha,linewidths=.8 if marker=='x' else .2)
    e.axhline(0,color=LIGHT,lw=.6);e.axvline(0,color=LIGHT,lw=.6)
    e.set_xlabel('Normal Δ early → DP');e.set_ylabel('Leukemia βR')
    panel(e,'E','Normal development ≠ leukemia program','13 concordant · 32 discordant · 45 weak · ρ = −0.50');clean_grid(e)
    cam=read(D/'a6_yayon/a6_3_cameraPR.tsv')
    subset=cam[cam['program'].astype(str).str.contains('concordant|neutral',case=False,regex=True)].copy()
    subset['name']=subset.program.map(lambda x:'Discordant / neutral' if 'discordant' in x else 'Normal-early concordant')
    subset['power']=-np.log10(subset.PValue.clip(lower=1e-20))
    f.barh(np.arange(len(subset)),subset.power,color=[Z2 if 'Discordant' in x else GREY for x in subset.name],height=.55)
    f.set_yticks(np.arange(len(subset)),['']*len(subset));f.invert_yaxis()
    for i,(name,n) in enumerate(zip(subset.name,subset.NGenes)):
        f.text(.5,i,f'{name} (n={int(n)})',fontsize=6.4,va='center',color='white' if 'Discordant' in name else DARK)
    f.set_xlabel('−log10 CAMERA P');panel(f,'F','Matched BCL11B–ETP partition','38-gene block: P = 3.9e−16');clean_grid(f,'x')
    save(fig,3)

def fig4():
    fig=mmfig(180,195);g=GridSpec(3,6,figure=fig,height_ratios=[.85,.95,1.25])
    a=fig.add_subplot(g[0,:3]);b=fig.add_subplot(g[0,3:]);c=fig.add_subplot(g[1,:2]);d=fig.add_subplot(g[1,3:]);e=fig.add_subplot(g[2,:4]);f=fig.add_subplot(g[2,4:])
    p=read(SD/'revision3_robustness/bcl11b_etp_matched_pairs.tsv')
    for i,r in p.iterrows():
        a.plot([0,1],[r.BCL11B_dev,r.ETP_dev],color=LIGHT,lw=.7)
        a.scatter([0,1],[r.BCL11B_dev,r.ETP_dev],color=[BCL,ETP],s=12)
        b.plot([0,1],[r.BCL11B_residual,r.ETP_residual],color=LIGHT,lw=.8)
        b.scatter([0,1],[r.BCL11B_residual,r.ETP_residual],color=[BCL,ETP],s=13)
    for ax in (a,b): ax.set_xticks([0,1],['BCL11B','Matched ETP-like']);ax.set_xlim(-.15,1.15);clean_grid(ax)
    panel(a,'A','Developmentally matched cases','18 pairs; median |ΔD| = 0.016');a.set_ylabel('Developmental proxy D')
    panel(b,'B','Paired residual contrast','Mean difference −1.92 (95% CI −2.98 to −0.67)');b.set_ylabel('Within-cohort R')
    h=read(D/'gate2_sensitivity/heldout_matched_cameraPR.tsv')
    h=h[h.program.isin(['ZEB1-side50','ZEB2-side50'])]
    c.barh(np.arange(len(h)),-np.log10(h.PValue),color=[Z1 if 'ZEB1' in x else Z2 for x in h.program],height=.55)
    c.set_yticks(np.arange(len(h)),h.program,fontsize=6);c.invert_yaxis();c.set_xlabel('−log10 CAMERA P')
    panel(c,'C','Held-out within-cohort bridge','Leave BCL11B and ETP-like out of program derivation');clean_grid(c,'x')
    gene=read(SD/'revision4_targeted/matched_BCL11B_ETP_gene_results.tsv')
    gene=gene.dropna(subset=['symbol','logFC']).sort_values('adj.P.Val').head(10).sort_values('logFC')
    d.barh(np.arange(len(gene)),gene.logFC,color=[Z2 if x>0 else Z1 for x in gene.logFC],height=.58)
    d.set_yticks(np.arange(len(gene)),gene.symbol,fontsize=6);d.axvline(0,color=GREY,lw=.55);d.set_xlabel('BCL11B − ETP log2 effect')
    panel(d,'D','Matched transcriptome','Top adjusted-P genes; Hallmark context in S5');clean_grid(d,'x')
    cases=read(SD/'F4/F4C_GSE162280_cases.tsv').sort_values('log2_ZEB1_over_ZEB2')
    yy=np.arange(len(cases))
    lesion_colors={'ZEB2 fusion':'#466D88','ARID1B enhancer':'#658A7F','CDK6 enhancer':'#AE8A55'}
    for i,r in enumerate(cases.itertuples()):
        group='ZEB2 fusion' if 'ZEB2' in str(r.lesion) else 'ARID1B enhancer' if 'ARID1B' in str(r.lesion) else 'CDK6 enhancer'
        e.plot([r.ZEB1_cpm,r.ZEB2_cpm],[i,i],color=lesion_colors[group],lw=1.0,alpha=.55)
        e.scatter([r.ZEB1_cpm],[i],color=Z1,s=18,zorder=3)
        e.scatter([r.ZEB2_cpm],[i],color=Z2,s=18,zorder=3)
    e.set_yticks(yy,cases['case'],fontsize=6.6)
    e.invert_yaxis();e.set_xlabel('Expression, CPM');e.set_xscale('symlog',linthresh=20)
    e.set_xticks([0,10,100,300]);e.set_xticklabels(['0','10','100','300'])
    panel(e,'E','Lesion-defined convergence','12/12 ZEB2-dominant; 7/7 non-ZEB2 partners');clean_grid(e,'x')
    adult=read(D/'a8_adult/a8_patient_scores.tsv');m=read(D/'a8_adult/a8_2_merged.tsv')
    # Official high-confidence BCL11B calls are identified in the locked A8-2 merged table.
    call='T-ALL main-cluster high-confidence'
    sel=m[call].astype(str).str.contains('BCL11B',case=False,na=False)
    xx=m.loc[sel,'R'].dropna().to_numpy()
    if len(xx)!=7:
        # The lock is authoritative; show the locked summary if the encoded classifier field changes.
        f.text(.5,.57,'7 high-confidence BCL11B cases',ha='center',transform=f.transAxes,fontsize=8,color=BCL)
    else:
        f.scatter(xx,np.arange(len(xx)),color=BCL,s=22);f.set_yticks([]);f.axvline(0,color=LIGHT,lw=.6);f.set_xlabel('Adult within-cohort R')
    panel(f,'F','Adult BCL11B context','n=7 · median R −2.34 · 7/7 negative');clean_grid(f,'x')
    save(fig,4)

def fig5():
    fig=mmfig(180,190);g=GridSpec(3,6,figure=fig,height_ratios=[.32,1.20,1.30])
    a=fig.add_subplot(g[0,:]);b=fig.add_subplot(g[1,:2]);c=fig.add_subplot(g[1,2:4]);d=fig.add_subplot(g[1,4:]);e=fig.add_subplot(g[2,:4]);f=fig.add_subplot(g[2,4:])
    a.axis('off');panel(a,'A','Independent evidence levels')
    levels=[('Bulk T-ALL','GSE146901'),('ETP / MPAL','GSE234608'),('Blast-enriched','GSE243914'),('Author malignant','GSE248287'),('Adult bulk','GSE280250')]
    for i,(lab,sub) in enumerate(levels):
        x=.01+i*.198
        a.add_patch(Rectangle((x,.30),.185,.34,transform=a.transAxes,facecolor=SOFT,edgecolor='none'))
        a.text(x+.09,.51,lab,ha='center',transform=a.transAxes,fontsize=6.8,weight='bold',color=DARK)
        a.text(x+.09,.35,sub,ha='center',transform=a.transAxes,fontsize=5.8,color=GREY)
    x1=read(D/'gate3_validation/gse146901_program_scores.tsv')
    x1['sample']=x1['sample'].astype(str).str.zfill(3)
    base146=read(SD/'F5/F5A-C_GSE146901_samples.tsv')[['sample','balance']]
    base146['sample']=base146['sample'].astype(str).str.zfill(3)
    x1=x1.merge(base146,on='sample',validate='one_to_one')
    assert len(x1)==18 and x1['balance'].notna().all()
    x2=read(D/'a5_gse234608/a5_patient_frozen_scores.tsv').rename(columns={'balance':'B'})
    x3=read(D/'a7b_malignant/a7b05_sample_scores.tsv')
    sources=[(b,x1,'GSE146901','balance'),(c,x2,'GSE234608','B'),(d,x3,'GSE243914','B')]
    metric_cols=['ZEB1_side50','ZEB2_side50']
    for letter,(ax,x,title,balance_col) in zip('BCD',sources):
        groups=x['group'].astype(str); uniq=groups.unique().tolist()
        for mi,(metric,col) in enumerate([(metric_cols[0],Z1),(metric_cols[1],Z2),(balance_col,DARK)]):
            for gi,u in enumerate(uniq):
                vals=pd.to_numeric(x.loc[groups==u,metric],errors='coerce').dropna().to_numpy()
                offset=(gi-(len(uniq)-1)/2)*.23
                xpos=np.full(len(vals),mi+offset)+jitter(len(vals),.028,mi*5+gi+1)
                ax.scatter(xpos,vals,s=12,c=col,marker='o' if gi==0 else 's',alpha=.52,linewidths=0)
                if len(vals):ax.hlines(np.median(vals),mi+offset-.09,mi+offset+.09,color=col,lw=1.8)
        ax.axhline(0,color=LIGHT,lw=.55)
        ax.set_xticks(range(3),['ZEB1-side','ZEB2-side','Balance'],fontsize=6.1,rotation=25,ha='right')
        if letter=='B':ax.set_ylabel('Within-cohort standardized value')
        labels='  ·  '.join(f'{"○" if i==0 else "□"} {u.replace("conventional_TALL","T-ALL").replace("ETP_or_MPAL","ETP/MPAL")} n={int((groups==u).sum())}' for i,u in enumerate(uniq))
        panel(ax,letter,title,labels);clean_grid(ax)
    names=['GSE146901','GSE234608','GSE243914','GSE248287','GSE280250']
    pct=[87,80,91.8,77.7,98.7];rho=[.60,.64,.77,.62,.86]
    y=np.arange(5)
    e.barh(y,pct,color=[GREY,GREY,Z2,ETP,DARK],height=.58)
    for i,(p,r) in enumerate(zip(pct,rho)):e.text(p+1,i,f'{p:g}%   ρ={r:.2f}',va='center',fontsize=6.2,color=DARK)
    e.set_xlim(0,120);e.set_yticks(y,names);e.invert_yaxis();e.set_xlabel('Frozen genes in expected direction (%)')
    panel(e,'E','Constituent-gene preservation','Recovered/tested n reported in legend and Table S7');clean_grid(e,'x')
    mal=read(D/'a7b_malignant/a7b_patient_scores.tsv')
    f.scatter(mal['B'],mal['ZEB1_side50'],color=Z1,s=18,label='ZEB1-side ρ=+.73')
    f.scatter(mal['B'],mal['ZEB2_side50'],color=Z2,s=18,label='ZEB2-side ρ=−.21')
    f.axvline(0,color=LIGHT,lw=.5);f.set_xlabel('Malignant pseudobulk B');f.set_ylabel('Frozen score')
    panel(f,'F','Author-defined malignant cells','15 patients; Z1 ρ=+.73; Z2 ρ=−.21 (CI spans 0)');clean_grid(f)
    save(fig,5)

def fig6():
    fig=mmfig(180,185);g=GridSpec(3,6,figure=fig,height_ratios=[1.05,1.1,.62])
    a=fig.add_subplot(g[0,:3]);b=fig.add_subplot(g[0,3:]);c=fig.add_subplot(g[1,:3]);d=fig.add_subplot(g[1,3:]);e=fig.add_subplot(g[2,:])
    score=read(D/'a7_lineage/scores/a7a_unit_scores_min20.tsv')
    groups=['HSC/MPP/HSPC','lymphoid progenitor/CLP','T lineage','Myeloid'];lab=['HSPC','CLP','T','Myeloid']
    for i,(gname,l) in enumerate(zip(groups,lab)):
        sub=score[score.lineage==gname]
        for j,col in enumerate(['ZEB1_side50','ZEB2_side50']):
            v=sub[col].dropna().to_numpy();a.scatter(np.full(len(v),i+j*.30)+jitter(len(v),.045,i+j),v,color=Z1 if j==0 else Z2,s=14,alpha=.75)
            if len(v):a.hlines(np.median(v),i+j*.30-.10,i+j*.30+.10,color=DARK,lw=1.5)
    a.set_xticks(np.arange(4)+.15,lab);a.set_ylabel('Donor×lineage frozen score')
    panel(a,'A','Normal marrow lineage affinity','≥20 cells per donor×lineage unit');clean_grid(a)
    don=read(D/'a7_lineage/scores/a7a_donor_min20.tsv')
    contrasts=[('Myeloid−T','d_My_minus_T_Z2',Z2),('HSPC−T','d_HSPC_minus_T_Z2',ETP),('CLP−T','d_CLP_minus_T_Z2',GREY)]
    for i,(name,col,color) in enumerate(contrasts):
        x=don[col].dropna().to_numpy();b.scatter(x,np.full(len(x),i)+jitter(len(x),.06,i),s=16,color=color,alpha=.8)
        if len(x):b.plot([np.median(x),np.median(x)],[i-.2,i+.2],color=DARK,lw=1.6)
    b.axvline(0,color=LIGHT,lw=.6);b.set_yticks(range(3),['']*3);b.invert_yaxis();b.set_xlim(-.55,1.95);b.set_xlabel('Paired donor ZEB2-side difference')
    for i,label in enumerate(['Myeloid 12/12','HSPC 11/11','CLP 0/8']):
        b.text(-.52,i-.20,label,fontsize=6.5,color=DARK)
    panel(b,'B','Donor-paired lineage contrasts','Myeloid 12/12 · HSPC 11/11 · CLP 0/8');clean_grid(b,'x')
    genes=read(D/'a7_lineage/scores/a7a_genes_min20.tsv').dropna(subset=['delta_Myeloid_minus_T','neg_beta'])
    c.scatter(genes.delta_Myeloid_minus_T,genes.neg_beta,s=16,color=[Z2 if str(s).startswith('ZEB2') else Z1 for s in genes.side],alpha=.7,edgecolors='white',linewidths=.2)
    c.axhline(0,color=LIGHT,lw=.5);c.axvline(0,color=LIGHT,lw=.5)
    c.set_xlabel('Normal Δ myeloid − T');c.set_ylabel('−βR in leukemia')
    panel(c,'C','Gene-level lineage affinity','75% expected direction · ρ = +0.60');clean_grid(c)
    adj=read(D/'a7_lineage/a7_0_frozen_genes_specimen_adjusted.tsv')
    d.scatter(adj.beta_primary,adj.beta_specimen,s=13,color=[Z2 if str(s).startswith('ZEB2') else Z1 for s in adj.side],alpha=.72)
    lim=max(abs(adj.beta_primary).max(),abs(adj.beta_specimen).max())*1.08
    d.plot([-lim,lim],[-lim,lim],color=GREY,lw=.7,ls='--')
    d.set_xlim(-lim,lim);d.set_ylim(-lim,lim);d.set_aspect('equal',adjustable='box')
    d.set_xlabel('Primary βR');d.set_ylabel('Specimen-adjusted βR')
    panel(d,'D','BM/PB specimen audit','100/100 same sign · ρ = .9995');clean_grid(d)
    e.axis('off');panel(e,'E','Alternative-explanation audit')
    steps=[('Normal marrow','myeloid > T in 12/12',Z2),('BM/PB source','100/100 β signs stable',GREY),('Blast-enriched','91.8% genes expected',ETP),('Author malignant','15 patient pseudobulks',DARK)]
    for i,(h,s,col) in enumerate(steps):
        x=.01+i*.246
        e.add_patch(Rectangle((x,.18),.229,.56,transform=e.transAxes,facecolor=SOFT,edgecolor='none'))
        e.add_patch(Rectangle((x,.18),.015,.56,transform=e.transAxes,facecolor=col,edgecolor='none'))
        e.text(x+.025,.53,h,transform=e.transAxes,fontsize=6.8,weight='bold',color=DARK)
        e.text(x+.025,.33,s,transform=e.transAxes,fontsize=6.6,color=GREY)
        if i<3:e.annotate('',xy=(x+.243,.46),xytext=(x+.232,.46),xycoords='axes fraction',arrowprops=dict(arrowstyle='->',lw=.8,color=GREY))
    save(fig,6)

def fig7():
    fig=mmfig(180,215);g=GridSpec(4,6,figure=fig,height_ratios=[.9,1.0,1.12,.88])
    a=fig.add_subplot(g[0,:]);b=fig.add_subplot(g[1,:3]);c=fig.add_subplot(g[1,3:]);d=fig.add_subplot(g[2,:4]);e=fig.add_subplot(g[2,4:]);f=fig.add_subplot(g[3,:3]);h=fig.add_subplot(g[3,3:])
    pt=read(D/'a8_adult/a8_patient_scores.tsv')
    a.scatter(pt.D,pt.B,s=15,color='#BAC4C9',alpha=.8)
    x=pt.sort_values('D');a.plot(x.D,x.B-x.R,color=DARK,lw=1.5)
    a.set_xlabel('Adult developmental proxy D');a.set_ylabel('Adult balance B')
    panel(a,'A','Adult cohort-internal developmental relationship','GSE280250 · 79 diagnosis samples; spline df = 3');clean_grid(a)
    for ax,col,color,letter,r,ci in [(b,'ZEB1_side50',Z1,'B','+.60','.43 to .73'),(c,'ZEB2_side50',Z2,'C','−.66','−.76 to −.52')]:
        ax.scatter(pt.R,pt[col],s=20,color=color,alpha=.75)
        z=np.polyfit(pt.R,pt[col],1);xx=np.linspace(pt.R.min(),pt.R.max(),100);ax.plot(xx,np.polyval(z,xx),color=DARK,lw=1)
        ax.axvline(0,color=LIGHT,lw=.55);ax.set_xlabel('Adult within-cohort R');ax.set_ylabel(col.replace('_side50','-side score'))
        panel(ax,letter,f'{col.replace("_side50", "-side")} versus R',f'Spearman ρ = {r}; bootstrap 95% CI {ci}');clean_grid(ax)
    ge=read(D/'a8_adult/a8_genes.tsv').dropna(subset=['beta_adult','beta_R'])
    d.scatter(ge.beta_R,ge.beta_adult,s=19,color=[Z1 if str(s).startswith('ZEB1') else Z2 for s in ge.side],alpha=.77,edgecolors='white',linewidths=.25)
    d.axhline(0,color=LIGHT,lw=.55);d.axvline(0,color=LIGHT,lw=.55)
    d.set_xlabel('Pediatric discovery βR');d.set_ylabel('Adult cohort βR')
    panel(d,'D','Constituent-gene effects replicate','86/100 recovered · 78 model-tested · 98.7% expected · ρ=.86');clean_grid(d)
    panel(e,'E','Age sensitivity','Descriptive Spearman correlations')
    for i,(label,subset) in enumerate([('All  n=79',pt),('Age ≥18  n=72',pt[pt.age>=18]),('Age ≥21  n=70',pt[pt.age>=21])]):
        for metric,color,offset in [('ZEB1_side50',Z1,-.12),('ZEB2_side50',Z2,.12)]:
            rho=spearmanr(subset.R,subset[metric]).statistic
            e.scatter(rho,i+offset,s=21,color=color,zorder=3)
        e.plot([spearmanr(subset.R,subset.ZEB2_side50).statistic,spearmanr(subset.R,subset.ZEB1_side50).statistic],[i,i],color=LIGHT,lw=.7)
    e.axvline(0,color=GREY,lw=.55);e.set_xlim(-.85,.85);e.set_yticks(range(3),['All 79','≥18 72','≥21 70'],fontsize=6.3);e.invert_yaxis();e.set_xlabel('ρ with adult R');clean_grid(e,'x')
    ff=[('BCL11B',7,-2.34),('ETP-like',47,-.06),('TAL1 DP-like',7,.32),('TLX3 immature',2,1.54)]
    for i,(name,n,val) in enumerate(ff):
        f.scatter(val,i,color=color_subtype(name),s=22);f.text(val+.07,i,f'n={n}',va='center',fontsize=6)
    f.set_yticks(range(4),['BCL11B','ETP-like','TAL1 DP','TLX3 imm.'],fontsize=6.5);f.invert_yaxis();f.axvline(0,color=LIGHT,lw=.55);f.set_xlabel('Adult median R')
    panel(f,'F','ALLCatchR2 subtype context','Only three unique mapped n≥5 classes enter ordering');clean_grid(f,'x')
    h.axis('off');panel(h,'G','Evidence chain and boundaries')
    nodes=[(.16,'4 + 6','thymus refs +\nspatial donors',Z1),
           (.50,'1,309','diagnostic patients\n17 subtypes',DARK),
           (.84,'5','external program\ncohorts',Z2)]
    for i,(x,number,description,col) in enumerate(nodes):
        h.add_patch(Circle((x,.65),.095,transform=h.transAxes,facecolor='white',edgecolor=col,lw=1.5))
        h.text(x,.65,number,transform=h.transAxes,ha='center',va='center',fontsize=8.6,weight='bold',color=col)
        h.text(x,.42,description,transform=h.transAxes,ha='center',va='top',fontsize=6.5,color=DARK)
        if i<2:
            h.add_patch(FancyArrowPatch((x+.11,.65),(nodes[i+1][0]-.11,.65),
                          transform=h.transAxes,arrowstyle='-|>',mutation_scale=8,lw=.85,color=GREY))
    h.text(.5,.10,'No cell-of-origin, acute-OE mechanism or clinical-utility claim',
           transform=h.transAxes,ha='center',va='center',fontsize=6.5,color=GREY)
    save(fig,7)

def main():
    for n,fn in enumerate([fig1,fig2,fig3,fig4,fig5,fig6,fig7],1):
        print('Figure',n,flush=True);fn()
    print('Exported to',OUT)

if __name__=='__main__': main()

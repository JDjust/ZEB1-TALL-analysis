"""Data-led second-pass plates for Figures 3, 4 and 7.

All quantitative marks originate in frozen source tables. This module changes
figure grammar and descriptive aggregation only; it does not refit an inference
model or choose genes from validation cohorts.
"""
from pathlib import Path
import math
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from matplotlib.patches import Rectangle, PathPatch, Circle
from matplotlib.path import Path as MplPath
from matplotlib.colors import LinearSegmentedColormap, Normalize
from scipy.stats import spearmanr

import editorial_rebuild_figures as base

SD, D, DATA, OUT = base.SD, base.D, base.DATA, base.OUT
read, panel, grid = base.read, base.panel, base.clean_grid
Z1, Z2, ETP, GREY, LIGHT, DARK = base.Z1, base.Z2, base.ETP, base.GREY, base.LIGHT, base.DARK
CMAP = LinearSegmentedColormap.from_list('zeb_diverging', [Z2, '#F7F8F8', Z1], N=256)
CMAP.set_bad('#E9EDEF')


def save(fig, number):
    fig.subplots_adjust(left=.088, right=.975, top=.945, bottom=.06,
                        wspace=.60 if number==5 else .34,
                        hspace=.42 if number==5 else .62)
    for t in fig.findobj(match=matplotlib.text.Text):
        if t.get_text().strip() and t.get_fontsize() < 6.5:
            t.set_fontsize(6.5)
    fig.savefig(OUT / f'Figure{number}.pdf', bbox_inches=None)
    fig.savefig(OUT / f'Figure{number}.png', dpi=300, bbox_inches=None)
    plt.close(fig)


def _ranked_program_selection(ax):
    """Rank curve is the panel; the measured selection flow is a small inset."""
    s = read(D / 'gate2_program/freeze_summary.tsv').iloc[0]
    assert (s.n_eligible, s.n_pass, s.n_pass_ZEB1_side,
            s.n_pass_ZEB2_side, s.n_frozen_ZEB1, s.n_frozen_ZEB2) == (24356, 3560, 499, 3061, 50, 50)
    passed=read(D/'gate2_program/effect_threshold_pass.tsv')
    assert len(passed)==3560
    panel(ax,'A','Frozen-program selection',
          'Partial-R² ranking; inset shows 24,356 → 3,560 → 499 / 3,061')
    flow=ax.inset_axes([.71,.43,.27,.43],facecolor='white')
    flow.set_xlim(0,1);flow.set_ylim(0,1);flow.axis('off')
    eligible=24356; selected=3560; rejected=eligible-selected
    pass_h=.73*selected/eligible; reject_h=.73-pass_h
    # Source and gate bars share one global count scale; rejected genes remain visible.
    flow.add_patch(Rectangle((.03,.12),.09,.73,facecolor=DARK,edgecolor='none'))
    flow.add_patch(Rectangle((.39,.12),.09,pass_h,facecolor=GREY,edgecolor='none'))
    flow.add_patch(Rectangle((.39,.25),.09,reject_h,facecolor='#DEE4E7',edgecolor='none'))
    def ribbon(y0,h0,y1,h1,color,alpha):
        x0,x1=.12,.39
        verts=[(x0,y0),(x0+.12,y0),(x1-.12,y1),(x1,y1),
               (x1,y1+h1),(x1-.12,y1+h1),(x0+.12,y0+h0),(x0,y0+h0),(x0,y0)]
        codes=[MplPath.MOVETO,MplPath.CURVE4,MplPath.CURVE4,MplPath.CURVE4,
               MplPath.LINETO,MplPath.CURVE4,MplPath.CURVE4,MplPath.CURVE4,MplPath.CLOSEPOLY]
        flow.add_patch(PathPatch(MplPath(verts,codes),facecolor=color,alpha=alpha,edgecolor='none'))
    ribbon(.12,pass_h,.12,pass_h,GREY,.48)
    ribbon(.12+pass_h,reject_h,.25,reject_h,'#DEE4E7',.75)
    # At the side split, the vertical total is normalized to the passing set.
    z1_h=.62*499/3560;z2_h=.62*3061/3560
    flow.add_patch(Rectangle((.75,.14),.09,z2_h,facecolor=Z2,edgecolor='none'))
    flow.add_patch(Rectangle((.75,.14+z2_h),.09,z1_h,facecolor=Z1,edgecolor='none'))
    flow.annotate('',xy=(.73,.43),xytext=(.50,.18),
                  arrowprops=dict(arrowstyle='->',lw=.8,color=GREY))
    flow.text(.10,.89,'24k',ha='center',fontsize=6.5,weight='bold',color=DARK)
    flow.text(.45,.02,'3.6k',ha='center',fontsize=6.5,weight='bold',color=DARK)
    flow.text(.83,.89,'499',ha='center',fontsize=6.5,color=Z1,weight='bold')
    flow.text(.83,.03,'3061',ha='center',fontsize=6.5,color=Z2,weight='bold')
    for sign,color,label,n in [(1,Z1,'ZEB1 side',499),(-1,Z2,'ZEB2 side',3061)]:
        q=passed[passed.beta_R*sign>0].sort_values('partial_R2',ascending=False)
        assert len(q)==n
        ranks=np.arange(1,n+1)
        ax.plot(ranks,q.partial_R2,color=color,alpha=.30,lw=1.0)
        ax.plot(ranks[:50],q.partial_R2.iloc[:50],color=color,lw=1.8,
                label=f'{label}: 50 / {n}')
        ax.scatter([50],[q.partial_R2.iloc[49]],s=13,color=color,zorder=3)
    ax.axvline(50,color=GREY,ls=':',lw=.6)
    ax.set_xscale('log');ax.set_xlim(.85,3400)
    ax.set_xticks([1,10,50,500,3000],['1','10','50','500','3k'])
    ax.set_xlabel('Partial R² rank')
    ax.set_ylabel('Partial R²')
    grid(ax,'y')


def _donor_trajectories(ax):
    normal=read(D/'a6_yayon/a6_1_donor_program.tsv')
    assert len(normal)==5
    for side,col in [('ZEB1',Z1),('ZEB2',Z2)]:
        keys=[f'early_{side}_side50',f'DP_{side}_side50',f'SP_{side}_side50']
        for _,r in normal.iterrows():
            ax.plot(range(3),[r[k] for k in keys],color=col,alpha=.25,lw=.7)
        ax.plot(range(3),[normal[k].median() for k in keys],color=col,lw=1.8,
                marker='o',ms=3,label=f'{side} side')
    ax.set_xticks(range(3),['Early','DP','SP'])
    ax.set_xlim(-.10,2.1)
    ax.set_ylabel('Score',fontsize=6.5)
    ax.text(.02,.98,'D  Normal donors (n=5)',transform=ax.transAxes,
            va='top',fontsize=6.5,weight='bold',color=DARK)
    grid(ax,'y')


def _decile_heatmap(ax):
    x = read(DATA / 'editorial_visual_sources/F3C_frozen100_residual_deciles.tsv')
    ann = read(DATA / 'editorial_visual_sources/F3C_residual_decile_annotations.tsv')
    assert len(x) == 100 and len(ann) == 10 and ann.n_patients.sum() == 1309
    x['order_side'] = np.where(x.side.eq('ZEB2-side50'), 0, 1)
    x = x.sort_values(['order_side','rank']).reset_index(drop=True)
    cols = [f'decile_{i}' for i in range(1,11)]
    z = x[cols].to_numpy(float)
    im = ax.imshow(z, aspect='auto', interpolation='nearest', origin='upper',
                   extent=(0,10,100,0), cmap=CMAP, vmin=-1.25, vmax=1.25,
                   rasterized=True)
    ax.set_xlim(-1.72,10.15); ax.set_ylim(101,-18.2)
    ax.axhline(50,color='white',lw=1.2)
    ax.set_xticks(np.arange(.5,10.5,1), [f'D{i}' for i in range(1,11)])
    ax.set_yticks([25,75],['ZEB2 · 50','ZEB1 · 50'],fontsize=6.6)
    ax.tick_params(axis='y',length=0)
    # Annotation columns describe distinct quantities and have their own keys.
    for i,r in ann.sort_values('decile').reset_index(drop=True).iterrows():
        vals = [(r.pct_BCL11B,Z2),(r.pct_ETP_like,ETP),(r.pct_TLX3,Z1)]
        left=i; total=0.
        for pct,col in vals:
            frac=min(max(float(pct)/100,0),1-total)
            if frac>0:
                ax.add_patch(Rectangle((left,-16),frac,3.0,facecolor=col,edgecolor='none'))
                left+=frac; total+=frac
        if total<1:
            ax.add_patch(Rectangle((left,-16),1-total,3.0,facecolor='#E5E9EB',edgecolor='none'))
        ax.add_patch(Rectangle((i,-11),1,3.0,facecolor=CMAP((np.clip(r.median_R,-2,2)+2)/4),edgecolor='white',lw=.3))
        ax.add_patch(Rectangle((i,-6),1,3.0,facecolor=CMAP((np.clip(r.median_D,-2,2)+2)/4),edgecolor='white',lw=.3))
    for y,label in [(-14.5,'Subtype'),(-9.5,'R'),(-4.5,'D')]:
        ax.text(-1.66,y,label,va='center',fontsize=6.5,color=DARK)
    ax.text(10.05,26,'ZEB2',rotation=90,ha='center',va='center',fontsize=6.5,color=Z2)
    ax.text(10.05,76,'ZEB1',rotation=90,ha='center',va='center',fontsize=6.5,color=Z1)
    ax.set_xlabel('Patient residual decile (low → high)')
    panel(ax,'C','Frozen 100-gene expression across residual deciles',
          'Discovery patients · median per-gene expression z score · 130–131 per decile')
    cb = ax.inset_axes([.72,-.17,.24,.035])
    plt.colorbar(im,cax=cb,orientation='horizontal',ticks=[-1,0,1])
    cb.set_xlabel('Median expression z',fontsize=6.5,labelpad=1)
    cb.tick_params(labelsize=6.5,length=1,pad=1)


def fig3():
    fig=base.mmfig(180,215)
    g=GridSpec(3,7,figure=fig,height_ratios=[.80,1.22,1.30])
    a=fig.add_subplot(g[0,:4]); b=fig.add_subplot(g[0,4:])
    c=fig.add_subplot(g[1,:]); e=fig.add_subplot(g[2,:4]); f=fig.add_subplot(g[2,4:])
    _ranked_program_selection(a)
    eff=read(D/'gate2_program/all_genes_residual_effects.tsv')
    fr=read(D/'gate2_program/frozen_ZEB_side50.tsv')
    b.scatter(eff.beta_R,eff.partial_R2,s=2,color='#C9D1D5',alpha=.55,rasterized=True)
    for side,col in [('ZEB1-side50',Z1),('ZEB2-side50',Z2)]:
        q=fr[fr.side==side]
        b.scatter(q.beta_R,q.partial_R2,s=11,color=col,alpha=.85,label=side,rasterized=True)
    b.axvline(0,color=LIGHT,lw=.6)
    b.set_xlabel('Leukemia residual coefficient βR');b.set_ylabel('Partial R²')
    b.legend(frameon=False,fontsize=6.5,loc='upper right')
    panel(b,'B','Genome-wide residual effects','100 frozen genes highlighted among 24,356 eligible genes');grid(b)
    _decile_heatmap(c)

    q=read(D/'a6_yayon/a6_2_gene_effects.tsv')
    for cls,marker,alpha in [('concordant','o',.9),('discordant','x',.9),('neutral','.',.55)]:
        t=q[q['class']==cls]
        e.scatter(t.delta_normal,t.beta_R,s=24 if marker!='.' else 12,
                  color=[Z1 if u.startswith('ZEB1') else Z2 for u in t.side],
                  marker=marker,alpha=alpha,linewidths=.7 if marker=='x' else .2,
                  rasterized=True)
    e.axhline(0,color=LIGHT,lw=.6);e.axvline(0,color=LIGHT,lw=.6)
    e.set_xlabel('Normal thymus early → DP effect');e.set_ylabel('Leukemia βR')
    panel(e,'E','Normal change versus leukemia association',
          '90 genes · 13 concordant · 32 discordant · 45 weak · ρ = −0.50');grid(e)
    d=e.inset_axes([.68,.56,.30,.40],facecolor='white')
    _donor_trajectories(d)
    z=q[q.side.eq('ZEB2-side50')].copy()
    matched=read(SD/'revision4_targeted/matched_BCL11B_ETP_gene_results.tsv')
    z=z.merge(matched[['symbol','logFC']],on='symbol',validate='one_to_one')
    assert len(z)==48 and (z['class'].eq('concordant').sum(),z['class'].isin(['neutral','discordant']).sum())==(10,38)
    z=z.sort_values('logFC').reset_index(drop=True)
    y=np.arange(len(z))
    for cls,col,label in [('concordant',GREY,'Normal-early concordant (10)'),
                          ('neutral',Z2,'Normal discordant/neutral (38)'),
                          ('discordant',Z2,None)]:
        mask=z['class'].eq(cls)
        f.hlines(y[mask],0,z.loc[mask,'logFC'],color=col,lw=.65,alpha=.65)
        f.scatter(z.loc[mask,'logFC'],y[mask],s=9,color=col,label=label,zorder=2,rasterized=True)
    f.axvline(0,color=LIGHT,lw=.6);f.set_ylim(-1,len(z));f.set_yticks([])
    f.set_xlabel('BCL11B − matched ETP log2 effect')
    # Group colors are defined in the panel subtitle and legend; preserve all marks.
    panel(f,'F','Matched gene effects','38 blue / 10 grey genes');grid(f,'x')
    save(fig,3)


def _overlap_and_heldout(ax):
    overlap=read(D/'gate2_sensitivity/heldout_vs_original_overlap.tsv')
    camera=read(D/'gate2_sensitivity/heldout_matched_cameraPR.tsv')
    assert set(overlap.n_overlap)=={39,29}
    assert overlap.n_orig.eq(50).all() and overlap.n_heldout.eq(50).all()
    ax.set_xlim(0,2.0);ax.set_ylim(0,1.03)
    ax.set_aspect('equal',adjustable='box',anchor='N');ax.axis('off')
    panel(ax,'C','Frozen-set overlap',
          'Original / held-out sets · 50 genes each')
    for center,side,color,dist in [(.48,'ZEB1-side50',Z1,.20),
                                   (1.52,'ZEB2-side50',Z2,.29)]:
        q=overlap.loc[overlap.side.eq(side)].iloc[0]
        shared=int(q.n_overlap);only=50-shared
        lx=center-dist/2;rx=center+dist/2;cy=.55;radius=.27
        # Two overlapping vector circles are mandatory. The printed counts are
        # exact; circle overlap area is illustrative rather than area-proportional.
        ax.add_patch(Circle((lx,cy),radius,facecolor=color,edgecolor=color,
                            alpha=.37,lw=.9))
        ax.add_patch(Circle((rx,cy),radius,facecolor=color,edgecolor=color,
                            alpha=.20,lw=1.0,linestyle='--'))
        for xpos,value in [(lx-radius*.67,only),(center,shared),(rx+radius*.67,only)]:
            ax.text(xpos,cy,str(value),ha='center',va='center',fontsize=7.0,
                    weight='bold',color=DARK)
        ax.text(center,.92,'ZEB1 side' if side.startswith('ZEB1') else 'ZEB2 side',
                ha='center',fontsize=7.0,weight='bold',color=color)
        p=float(camera.loc[camera.program.eq(side),'PValue'].iloc[0])
        ptxt='$0.68$' if side.startswith('ZEB1') else '$2.0\\times10^{-11}$'
        ax.text(center,.21,'CAMERA P',ha='center',fontsize=6.5,color=color)
        ax.text(center,.11,ptxt,ha='center',fontsize=6.5,color=color)
    ax.text(1.0,.01,'solid = original   ·   dashed = held-out',
            ha='center',fontsize=6.5,color=GREY)


def _matched_volcano_hallmark(ax, fig, sub_spec):
    nested=sub_spec.subgridspec(1,2,width_ratios=[2.25,1],wspace=.48)
    v=fig.add_subplot(nested[0,0]); h=fig.add_subplot(nested[0,1])
    genes=read(SD/'revision4_targeted/matched_BCL11B_ETP_gene_results.tsv')
    frozen=read(D/'gate2_program/frozen_ZEB_side50.tsv')[['symbol','side']]
    genes=genes.merge(frozen,on='symbol',how='left',validate='one_to_one')
    yy=-np.log10(genes['adj.P.Val'].clip(lower=1e-300))
    v.scatter(genes.logFC,yy,s=2,color='#CAD1D5',alpha=.45,rasterized=True)
    for side,col in [('ZEB1-side50',Z1),('ZEB2-side50',Z2)]:
        t=genes.side.eq(side)
        v.scatter(genes.loc[t,'logFC'],yy[t],s=9,color=col,alpha=.88,rasterized=True)
    v.axhline(-np.log10(.05),color=GREY,lw=.5,ls='--')
    v.axvline(0,color=LIGHT,lw=.5)
    v.set_xlim(-7.7,7.4);v.set_ylim(0,7.5)
    v.set_xlabel('Pair-blocked BCL11B − ETP log2 effect');v.set_ylabel('−log10 BH FDR')
    panel(v,'D','Matched transcriptome',
          '18 pairs; frozen genes highlighted');grid(v)
    pathways=read(SD/'revision4_targeted/matched_BCL11B_ETP_hallmark_camera.tsv')
    pathways=pathways[pathways.FDR<.05].sort_values('FDR').head(6).copy()
    label_map={
      'HALLMARK_INTERFERON_GAMMA_RESPONSE':'IFNγ',
      'HALLMARK_IL6_JAK_STAT3_SIGNALING':'IL6/JAK',
      'HALLMARK_HEME_METABOLISM':'Heme',
      'HALLMARK_ALLOGRAFT_REJECTION':'Allograft',
      'HALLMARK_INTERFERON_ALPHA_RESPONSE':'IFNα',
      'HALLMARK_OXIDATIVE_PHOSPHORYLATION':'OxPhos',
    }
    py=np.arange(len(pathways))
    px=-np.log10(pathways.FDR)
    h.scatter(px,py,s=pathways.NGenes*.18,color='#626F77',alpha=.85,zorder=2)
    h.hlines(py,0,px,color=LIGHT,lw=.8)
    h.set_yticks(py,[label_map.get(t,t.replace('HALLMARK_','')) for t in pathways.pathway],fontsize=6.5)
    h.invert_yaxis();h.set_xlabel('−log10 FDR')
    h.text(.02,1.03,'Hallmark CAMERA',transform=h.transAxes,fontsize=6.5,weight='bold',color=DARK)
    grid(h,'x')


def _adult_bcl_matrix(ax):
    m=read(D/'a8_adult/a8_2_merged.tsv')
    call='T-ALL main-cluster high-confidence'
    z=m[m[call].astype(str).str.contains('BCL11B',case=False,na=False)].copy()
    assert len(z)==7 and z.R.lt(0).all() and z.age.between(31,75).all()
    assert z['T-ALL immature high-confidence'].isna().all()
    z=z.sort_values('R').reset_index(drop=True)
    ax.set_xlim(0,1.06);ax.set_ylim(0,8);ax.axis('off')
    panel(ax,'F','Adult BCL11B patient matrix',
          '7/7 negative R · ages 31–75 · 0/7 ETP-like')
    for xpos,head in [(.07,'Age'),(.28,'R'),(.65,'B'),(.77,'Z1'),(.89,'Z2'),(1.01,'ETP')]:
        ax.text(xpos,7.55,head,transform=ax.transData,ha='center',fontsize=6.5,
                weight='bold',color=DARK)
    ax.plot([.56,.56],[.5,7.35],color=GREY,lw=.7)
    for i,r in z.iterrows():
        yy=6.9-i*.86
        ax.text(.07,yy,f'{int(r.age)}',ha='center',va='center',fontsize=6.5)
        xx=.56+.047*float(r.R)
        ax.plot([xx,.56],[yy,yy],color=Z2,lw=1)
        ax.scatter([xx],[yy],s=13,color=Z2,zorder=3)
        ax.text(.28,yy,f'{r.R:.2f}',ha='center',va='center',fontsize=6.5,color=DARK)
        for xpos,val in [(.65,r.B),(.77,r.ZEB1_side50),(.89,r.ZEB2_side50)]:
            col=CMAP((np.clip(float(val),-1.5,1.5)+1.5)/3)
            ax.add_patch(Rectangle((xpos-.052,yy-.30),.104,.60,facecolor=col,edgecolor='white',lw=.4))
            ax.text(xpos,yy,f'{val:+.1f}',ha='center',va='center',fontsize=6.5,color=DARK)
        ax.add_patch(Rectangle((.96,yy-.30),.10,.60,facecolor='#E9ECEE',
                               edgecolor='white',lw=.4))
        ax.text(1.01,yy,'No',ha='center',va='center',fontsize=6.5,color=GREY)


def fig4():
    fig=base.mmfig(180,195)
    g=GridSpec(3,6,figure=fig,height_ratios=[.90,1.08,1.32])
    a=fig.add_subplot(g[0,:3]);b=fig.add_subplot(g[0,3:])
    c=fig.add_subplot(g[1,:2]);e=fig.add_subplot(g[2,:4]);f=fig.add_subplot(g[2,4:])
    p=read(SD/'revision3_robustness/bcl11b_etp_matched_pairs.tsv')
    assert len(p)==18
    for _,r in p.iterrows():
        a.plot([0,1],[r.BCL11B_dev,r.ETP_dev],color=LIGHT,lw=.7)
        a.scatter([0,1],[r.BCL11B_dev,r.ETP_dev],color=[Z2,ETP],s=12)
        b.plot([0,1],[r.BCL11B_residual,r.ETP_residual],color=LIGHT,lw=.7)
        b.scatter([0,1],[r.BCL11B_residual,r.ETP_residual],color=[Z2,ETP],s=13)
    for ax in (a,b):
        ax.set_xticks([0,1],['BCL11B','Matched ETP-like']);ax.set_xlim(-.15,1.15);grid(ax)
    a.set_ylabel('Developmental-expression proxy D')
    panel(a,'A','Developmentally matched cases','18 pairs; median |ΔD| = 0.016')
    b.set_ylabel('Within-cohort residual R')
    panel(b,'B','Paired residual contrast','Mean difference −1.92; bootstrap 95% CI −2.98 to −0.67')
    _overlap_and_heldout(c)
    _matched_volcano_hallmark(None,fig,g[1,2:])
    cases=read(SD/'F4/F4C_GSE162280_cases.tsv').sort_values('log2_ZEB1_over_ZEB2')
    assert len(cases)==12 and cases.log2_ZEB1_over_ZEB2.lt(0).all()
    yy=np.arange(12)
    lesion_colors={'ZEB2 fusion':'#466D88','ARID1B enhancer':'#658A7F','CDK6 enhancer':'#AE8A55'}
    for i,r in enumerate(cases.itertuples()):
        group='ZEB2 fusion' if 'ZEB2' in str(r.lesion) else 'ARID1B enhancer' if 'ARID1B' in str(r.lesion) else 'CDK6 enhancer'
        e.plot([r.ZEB1_cpm,r.ZEB2_cpm],[i,i],color=lesion_colors[group],lw=1.0,alpha=.55)
        e.scatter([r.ZEB1_cpm],[i],color=Z1,s=16,zorder=3)
        e.scatter([r.ZEB2_cpm],[i],color=Z2,s=16,zorder=3)
    e.set_yticks(yy,cases['case'],fontsize=6.5);e.invert_yaxis()
    e.set_xscale('symlog',linthresh=20);e.set_xticks([0,10,100,300],['0','10','100','300'])
    e.set_xlabel('Measured expression (CPM)')
    panel(e,'E','Lesion-defined convergence','12/12 ZEB2-dominant; 7/7 without ZEB2 partner');grid(e,'x')
    _adult_bcl_matrix(f)
    save(fig,4)


def _adult_main_clusters(ax, label_ax):
    patients=read(D/'a8_adult/a8_2_merged.tsv')
    main='T-ALL main-cluster high-confidence'
    immature='T-ALL immature high-confidence'
    def group(v):
        if pd.isna(v): return 'No HC main'
        s=str(v)
        if ';' in s: return 'Multiple HC'
        code=s.split(' ',1)[0]
        if code in ('C9','C14','C15'): return 'Other HC (n=1 classes)'
        if code=='C8': return 'C8 BCL11B'
        if code=='C2': return 'C2 TAL1 DP'
        if code=='C17': return 'C17 TLX3 imm.'
        return code
    patients['display_group']=patients[main].map(group)
    assert patients.display_group.notna().all() and len(patients)==79
    stats=patients.groupby('display_group').agg(n=('R','size'),median=('R','median'),
                                                 etp=(immature,lambda x:x.notna().mean())).reset_index()
    stats=stats.sort_values('median').reset_index(drop=True)
    assert stats.n.sum()==79 and len(stats)==11
    for i,r in stats.iterrows():
        vals=patients.loc[patients.display_group==r.display_group,'R'].to_numpy()
        yy=np.full(len(vals),i)+base.jitter(len(vals),.065,i+11)
        color=Z2 if 'BCL11B' in r.display_group else Z1 if 'TLX3' in r.display_group else '#9EAAAF'
        ax.scatter(vals,yy,s=9,color=color,alpha=.75,rasterized=True,zorder=2)
        ax.plot([r['median']-.07,r['median']+.07],[i,i],color=DARK,lw=1.8,zorder=3)
        ax.scatter([5.05],[i],s=18+45*r.etp,color=ETP,alpha=.20+.75*r.etp,edgecolors='none')
        ax.text(5.32,i,f'{int(round(100*r.etp))}%',fontsize=6.5,va='center',ha='left',color=DARK)
        label_ax.text(.98,i,f'{r.display_group} ({int(r.n)})',fontsize=6.5,
                      va='center',ha='right',color=DARK)
    ax.set_yticks([])
    ax.invert_yaxis();ax.set_xlim(-3.3,6.2);ax.axvline(0,color=LIGHT,lw=.5)
    label_ax.set_xlim(0,1);label_ax.set_ylim(ax.get_ylim());label_ax.axis('off')
    ax.set_xlabel('Adult within-cohort residual R')
    panel(ax,'E','Adult molecular clusters',
          '79 patients; right bubbles show ETP-like fraction')
    grid(ax,'x')


def _integrated_gene_matrix(ax):
    fr=read(D/'gate2_program/frozen_ZEB_side50.tsv')[['symbol','side','rank','beta_R']]
    assert len(fr)==100 and fr.symbol.is_unique
    pieces=[('Discovery βR',fr.set_index('symbol').beta_R),
            ('Thymus E−DP',read(D/'a6_yayon/a6_2_gene_effects.tsv').set_index('symbol').delta_normal),
            ('Marrow My−T',read(D/'a7_lineage/scores/a7a_genes_min20.tsv').set_index('symbol').delta_Myeloid_minus_T),
            ('GSE146901 non−ETP',-read(D/'a65_external/a65_GSE146901_genes.tsv').set_index('symbol').delta_external),
            ('GSE234608 T−ETP',-read(D/'a65_external/a65_GSE234608_genes.tsv').set_index('symbol').delta_external),
            ('GSE243914 non−ETP',-read(D/'a7b_malignant/a7b05_genes.tsv').set_index('symbol').delta),
            ('Malignant ρ(E,B)',read(D/'a7b_malignant/a7b_genes.tsv').set_index('symbol').rho_E_vs_B),
            ('Adult βR',read(D/'a8_adult/a8_genes.tsv').set_index('symbol').beta_adult)]
    out=fr.copy()
    for name,series in pieces[1:]:
        assert series.index.is_unique,name
        out[name]=out.symbol.map(series)
    out=out.sort_values(['side','rank'],ascending=[False,True]).reset_index(drop=True)
    # The alphabetical descending order puts ZEB2-side first, ZEB1-side second.
    assert out.side.iloc[:50].eq('ZEB2-side50').all() and out.side.iloc[50:].eq('ZEB1-side50').all()
    raw=np.column_stack([out['beta_R'].to_numpy(float)]+[out[name].to_numpy(float) for name,_ in pieces[1:]])
    assert raw.shape==(100,8)
    scales=np.nanquantile(np.abs(raw),.9,axis=0)
    assert np.isfinite(scales).all() and (scales>0).all()
    display=np.clip(raw/scales,-1,1)
    source=out[['symbol','side','rank']].copy()
    for j,(name,_) in enumerate(pieces):
        source[name]=raw[:,j]
        source[f'scaled_{j+1}']=display[:,j]
    source.to_csv(DATA/'editorial_visual_sources/F7G_frozen100_context_effects.tsv',sep='\t',index=False)
    im=ax.imshow(display,aspect='auto',interpolation='nearest',origin='upper',
                 extent=(0,8,100,0),cmap=CMAP,vmin=-1,vmax=1,rasterized=True)
    ax.set_xlim(-.95,8);ax.set_ylim(101,-6.4)
    ax.axhline(50,color='white',lw=1.2)
    labels=['Pediatric\nβR','Thymus\nE→DP','Marrow\nMy−T','GSE146901\nnon−ETP',
            'GSE234608\nT−ETP','GSE243914\nnon−ETP','Malignant\nρ(E,B)','Adult\nβR']
    ax.set_xticks(np.arange(.5,8.5),labels,fontsize=6.5)
    ax.set_yticks([25,75],['ZEB2 · 50','ZEB1 · 50'],fontsize=6.5)
    ax.tick_params(axis='y',length=0)
    for i,(name,_) in enumerate(pieces):
        typ='NORMAL' if i in (1,2) else 'LEUKEMIA'
        ax.add_patch(Rectangle((i,-5.9),1,2.2,facecolor='#DDE4E7' if typ=='NORMAL' else '#D9E5EA',edgecolor='white',lw=.4))
    # Fixed biological annotations: A6 normal class and A7C myeloid affinity.
    a6=read(D/'a6_yayon/a6_2_gene_effects.tsv').set_index('symbol')['class']
    a7=read(D/'a7c_decompose/a7c_zeb2_split.tsv').set_index('symbol')['subset']
    for i,r in out.iterrows():
        ax.add_patch(Rectangle((-.19,i),.13,1,facecolor=Z2 if i<50 else Z1,edgecolor='none'))
        cls=a6.get(r.symbol,None)
        col={'concordant':Z1,'discordant':Z2,'neutral':'#E6EAEC'}.get(cls,'#C7CED1')
        ax.add_patch(Rectangle((-.39,i),.13,1,facecolor=col,edgecolor='none'))
        typ=a7.get(r.symbol,None)
        col='#4E8277' if typ=='myeloid_concordant' else '#AAB4B8' if typ=='non_myeloid_concordant' else 'white'
        ax.add_patch(Rectangle((-.59,i),.13,1,facecolor=col,edgecolor='none'))
    # Compact, explicit scale key: each column is divided by its own 90th |effect|.
    for j,level in enumerate(np.linspace(-1,1,9)):
        ax.add_patch(Rectangle((5.85+j*.23,-3.15),.23,1.25,facecolor=CMAP((level+1)/2),edgecolor='none'))
    ax.text(5.82,-3.6,'−1',fontsize=6.5,ha='left')
    ax.text(6.77,-3.6,'0',fontsize=6.5,ha='center')
    ax.text(7.69,-3.6,'+1',fontsize=6.5,ha='right')
    panel(ax,'F','Frozen gene effects across contexts',
          '100 genes; columnwise 90th-percentile scaling; grey = unavailable')


def fig7():
    fig=base.mmfig(180,220)
    g=GridSpec(4,6,figure=fig,height_ratios=[.72,.92,1.18,1.75])
    a=fig.add_subplot(g[0,:]);b=fig.add_subplot(g[1,:3]);c=fig.add_subplot(g[1,3:])
    d=fig.add_subplot(g[2,:3])
    fg=g[2,3:].subgridspec(1,4,width_ratios=[1.5,1,1,1],wspace=.03)
    fl=fig.add_subplot(fg[0,0]);f=fig.add_subplot(fg[0,1:]);h=fig.add_subplot(g[3,:])
    pt=read(D/'a8_adult/a8_patient_scores.tsv')
    assert len(pt)==79
    a.scatter(pt.D,pt.B,s=14,color='#B8C3C8',alpha=.78,rasterized=True)
    x=pt.sort_values('D');a.plot(x.D,x.B-x.R,color=DARK,lw=1.4)
    a.set_xlabel('Adult developmental-expression proxy D');a.set_ylabel('Adult ZEB balance B')
    panel(a,'A','Adult developmental relationship',
          'GSE280250 · 79 patients · independent df = 3 spline');grid(a)
    for ax,col,color,letter,rho,ci in [(b,'ZEB1_side50',Z1,'B','+.60','+.43 to +.73'),
                                        (c,'ZEB2_side50',Z2,'C','−.66','−.76 to −.52')]:
        ax.scatter(pt.R,pt[col],s=18,color=color,alpha=.75,rasterized=True)
        z=np.polyfit(pt.R,pt[col],1);xx=np.linspace(pt.R.min(),pt.R.max(),100)
        ax.plot(xx,np.polyval(z,xx),color=DARK,lw=.9)
        ax.axvline(0,color=LIGHT,lw=.5)
        ax.set_xlabel('Adult within-cohort R');ax.set_ylabel(col.replace('_side50','-side score'))
        panel(ax,letter,f'{col.replace("_side50","-side")} versus R',
              f'Spearman ρ = {rho}; bootstrap 95% CI {ci}');grid(ax)
    # Small age sensitivity is an inset in the score panel.
    e=c.inset_axes([.70,.56,.28,.40],facecolor='white')
    for i,(name,q) in enumerate([('All 79',pt),('≥18 72',pt[pt.age>=18]),('≥21 70',pt[pt.age>=21])]):
        for metric,col,off in [('ZEB1_side50',Z1,-.11),('ZEB2_side50',Z2,.11)]:
            e.scatter(spearmanr(q.R,q[metric]).statistic,i+off,s=9,color=col)
    e.axvline(0,color=LIGHT,lw=.4);e.set_xlim(-.85,.85)
    e.set_yticks(range(3),['79','72','70'],fontsize=6.5);e.invert_yaxis()
    e.set_xticks([-.5,0,.5]);e.tick_params(labelsize=6.5,length=1,pad=1)
    e.text(.22,.98,'Age',transform=e.transAxes,va='top',fontsize=6.5,
           weight='bold',color=DARK)
    e.spines[['top','right']].set_visible(False)
    ge=read(D/'a8_adult/a8_genes.tsv').dropna(subset=['beta_adult','beta_R'])
    d.scatter(ge.beta_R,ge.beta_adult,s=18,
              color=[Z1 if str(s).startswith('ZEB1') else Z2 for s in ge.side],
              alpha=.78,edgecolors='white',linewidths=.2,rasterized=True)
    d.axhline(0,color=LIGHT,lw=.5);d.axvline(0,color=LIGHT,lw=.5)
    d.set_xlabel('Pediatric discovery βR');d.set_ylabel('Adult cohort βR')
    panel(d,'D','Constituent-gene effects replicate',
          '86/100 recovered · 78 model-tested · 98.7% expected · ρ = .86');grid(d)
    _adult_main_clusters(f,fl)
    _integrated_gene_matrix(h)
    save(fig,7)


def main():
    for n,fn in [(3,fig3),(4,fig4),(7,fig7)]:
        print('Upgrading Figure',n,flush=True)
        fn()
    print('Updated',OUT)


if __name__=='__main__':
    main()

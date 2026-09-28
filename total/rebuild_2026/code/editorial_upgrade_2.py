"""Figure 2: measured subtype landscape and external slopegraph.

All inferential values are read from frozen tables. The only calculations here
are descriptive medians and column-wise color scaling for display.
"""
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from matplotlib.colors import TwoSlopeNorm
from matplotlib.patches import Rectangle
from scipy.stats import gaussian_kde

import editorial_rebuild_figures as b
from editorial_upgrade_347 import OUT, CMAP

D, SD = b.D, b.SD
Z1,Z2,ETP,GREY,LIGHT,DARK = b.Z1,b.Z2,b.ETP,b.GREY,b.LIGHT,b.DARK
read,panel,grid = b.read,b.panel,b.clean_grid

def sources():
    pt=read(SD/'F3/F3B-C_patient_level.tsv')
    st=read(SD/'F3/F3D-F_subtype_statistics.tsv').sort_values('residual_within_median')
    scores=read(D/'a7_lineage/a7_0_patient_scores.tsv')
    ex=read(D/'aieop120_official/aieop120_subtype_residual_medians.tsv')
    assert len(pt)==len(scores)==1309 and len(st)==17 and pt.sample_id.is_unique
    assert scores.sample_id.is_unique and pt.subtype.nunique()==17
    x=pt.merge(scores[['sample_id','ZEB1_side50','ZEB2_side50']],on='sample_id',validate='one_to_one')
    assert len(x)==1309 and x[['ZEB1_side50','ZEB2_side50']].notna().all().all()
    return x,st,ex

def landscape(ax,pt,st):
    names=st.subtype.tolist(); n=len(names)
    med=pt.groupby('subtype').agg(B=('balance','median'),Z1=('ZEB1_side50','median'),
                                  Z2=('ZEB2_side50','median'))
    mat=pd.DataFrame({'D':st.dev_median.to_numpy(),
                      'B':med.loc[names,'B'].to_numpy(),
                      'R':st.residual_within_median.to_numpy(),
                      'Z1':med.loc[names,'Z1'].to_numpy(),
                      'Z2':med.loc[names,'Z2'].to_numpy()})
    scaled=(mat-mat.mean())/mat.std(ddof=0)
    # The color scale encodes ZEB1-side orientation; raw Z2 values are retained
    # in the exported source table and have their display sign reversed only.
    scaled['Z2']=-scaled['Z2']
    assert scaled.shape==(17,5) and np.isfinite(scaled.to_numpy()).all()
    ax.imshow(scaled.to_numpy(),cmap=CMAP,norm=TwoSlopeNorm(vmin=-2.5,vcenter=0,vmax=2.5),
              aspect='auto',extent=(-.5,4.5,16.5,-.5),interpolation='nearest')
    for y,name in enumerate(names):
        if name in ['BCL11B','ETP-like','TLX3']:
            ax.add_patch(Rectangle((-.49,y-.48),4.98,.96,fill=False,
                                   edgecolor=b.color_subtype(name),lw=.95))
    ax.set_xticks(range(5),['D','B','R','Z1','−Z2'],fontsize=6.8)
    ax.tick_params(axis='x',length=0,pad=3)
    ax.set_yticks(range(n),[b.short_subtype(v) for v in names],fontsize=6.5)
    ax.tick_params(axis='y',length=0,pad=2)
    ax.set_xlim(-.5,6.55);ax.set_ylim(16.5,-.5)
    for y,count in enumerate(st.n):
        ax.plot([5.02,5.02+1.0*count/st.n.max()],[y,y],color='#76848A',lw=2.1,
                solid_capstyle='round')
        ax.text(6.45,y,str(int(count)),va='center',ha='right',fontsize=6.5,color=DARK)
    panel(ax,'A','Seventeen-subtype landscape',
          'Column z scores; Z2 sign reversed for common color orientation')
    ax.spines[:].set_visible(False)
    out=b.DATA/'editorial_visual_sources';out.mkdir(parents=True,exist_ok=True)
    raw=mat.copy();raw.insert(0,'subtype',names);raw['n']=st.n.to_numpy()
    raw.to_csv(out/'F2A_subtype_descriptive_medians.tsv',sep='\t',index=False)

def patient_geometry(ax,pt):
    ax.scatter(pt.dev,pt.balance,s=2.2,color='#B9C2C7',alpha=.42,rasterized=True)
    for name in ['BCL11B','ETP-like','TLX3']:
        x=pt[pt.subtype.eq(name)]
        ax.scatter(x.dev,x.balance,s=10,color=b.color_subtype(name),alpha=.80,
                   label=f'{name} (n={len(x)})',rasterized=True)
    q=read(SD/'F3/F3B_within_cohort_curve.tsv')
    ax.plot(q.dev,q.expected_within,color=DARK,lw=1.45)
    panel(ax,'B','Patient balance versus developmental proxy',
          'Cohort-internal natural spline, df = 3; n = 1,309')
    ax.set_xlabel('Developmental-expression proxy D');ax.set_ylabel('ZEB balance B')
    ax.legend(frameon=False,fontsize=6.5,loc='lower right');grid(ax)

def variance(ax):
    x=read(SD/'F2/F2F_model_variance.tsv')
    assert len(x)==4
    y=np.arange(4); val=x.r2.to_numpy()*100
    cols=[GREY,Z2,ETP,DARK]
    for i,v in enumerate(val):
        ax.plot([0,v],[i,i],color=LIGHT,lw=1.2)
        ax.scatter(v,i,s=45,color=cols[i],zorder=2)
        ax.text(v+1.1,i,f'{v:.1f}%',va='center',fontsize=6.7,color=DARK)
    ax.set_yticks(y,['Immunophenotype','Developmental proxy',
                     'Molecular subtype','D + subtype'],fontsize=6.5)
    ax.invert_yaxis();ax.set_xlim(0,48)
    ax.set_xlabel('Raw variance explained, R² (%)')
    panel(ax,'C','Model variance',
          f'Subtype adds +{val[3]-val[1]:.1f} pp beyond D; CV R² in S2')
    grid(ax,'x')

def raincloud(ax,pt,st):
    rng=np.random.default_rng(71)
    for i,name in enumerate(st.subtype):
        vals=pt.loc[pt.subtype.eq(name),'residual_within'].dropna().to_numpy()
        col=b.color_subtype(name)
        if len(vals)>=5 and np.std(vals)>0:
            xs=np.linspace(max(-6.4,vals.min()),min(5.5,vals.max()),90)
            den=gaussian_kde(vals)(xs);den=den/den.max()*.30
            ax.fill_between(xs,i,i-den,color=col,alpha=.28,lw=0)
            ax.plot(xs,i-den,color=col,lw=.55)
        jitter=rng.uniform(.05,.29,len(vals))
        ax.scatter(vals,i+jitter,s=1.7,color=col,alpha=.43,rasterized=True)
        q1,q3=np.quantile(vals,[.25,.75]);med=np.median(vals)
        ax.plot([q1,q3],[i+.01,i+.01],color=DARK,lw=1.3)
        ax.scatter(med,i+.01,s=11,color=DARK,zorder=3)
    ax.set_yticks(range(17),[b.short_subtype(v) for v in st.subtype],fontsize=6.5)
    ax.set_ylim(16.6,-.7);ax.axvline(0,color=GREY,lw=.55)
    ax.set_xlabel('Patient within-cohort residual R')
    panel(ax,'D','Patient residual distributions',
          'Half density, individual patients, median and IQR')
    grid(ax,'x')

def subtype_ci(ax,st):
    y=np.arange(17)
    for i,r in st.reset_index(drop=True).iterrows():
        col=b.color_subtype(r.subtype)
        ax.plot([r.residual_within_ci_low,r.residual_within_ci_high],[i,i],
                color=col,lw=1.5)
        ax.scatter(r.residual_within_median,i,color=col,s=15,zorder=2)
    ax.set_yticks(y,[b.short_subtype(v) for v in st.subtype],fontsize=6.5)
    ax.set_ylim(16.5,-.5);ax.axvline(0,color=GREY,lw=.55)
    ax.set_xlabel('Subtype median R (patient bootstrap 95% CI)')
    panel(ax,'E','Subtype poles and uncertainty',
          'BCL11B −1.76; ETP-like −0.23; TLX3 +0.64')
    grid(ax,'x')

def external_slope(ax,ex):
    q=ex[ex.n_aieop>=5].copy()
    assert len(q)==6
    ax.set_xlim(-.45,1.65);ax.set_ylim(-2.25,1.03)
    ax.axhline(0,color=LIGHT,lw=.5)
    labels=[]
    for _,r in q.iterrows():
        col=b.color_subtype(r.predicted_subtype)
        y0=float(r.residual_median_polonen);y1=float(r.residual_median_aieop)
        ax.plot([0,1],[y0,y1],color=col,lw=1.2,alpha=.9)
        ax.scatter([0,1],[y0,y1],color=col,s=19,zorder=2)
        short=b.short_subtype(r.predicted_subtype).replace('TAL1 immature-like','TAL1 imm.')
        labels.append((y1,short,col))
    labels.sort(key=lambda v:v[0]); placed=[]
    for y,name,col in labels:
        target=max(y,placed[-1][0]+.17) if placed else y
        placed.append((target,name,col,y))
    for target,name,col,actual in placed:
        ax.plot([1.02,1.09],[actual,target],color=col,lw=.45,alpha=.6)
        ax.text(1.11,target,name,va='center',fontsize=6.5,color=col)
    ax.set_xticks([0,1],['Pölönen','AIEOP'],fontsize=6.8)
    ax.set_ylabel('Subtype median R')
    panel(ax,'F','Discovery to AIEOP subtype comparison',
          'Six classes with AIEOP n ≥ 5; rank uncertainty remains wide')
    # The prespecified BCL11B anchor has n=4 and is shown as an inset.
    anchor=ex[ex.predicted_subtype.eq('BCL11B')].iloc[0]
    ax.plot([.12,.12],[.02,.18],transform=ax.transAxes,color=b.BCL,lw=1.1)
    ax.text(.17,.10,f'BCL11B n=4: {anchor.residual_median_aieop:+.2f}',
            transform=ax.transAxes,fontsize=6.5,color=b.BCL,va='center')
    grid(ax,'y')

def fig2():
    pt,st,ex=sources()
    fig=b.mmfig(180,213)
    g=GridSpec(3,6,figure=fig,height_ratios=[1.0,1.18,1.0])
    a=fig.add_subplot(g[0,:3]);bp=fig.add_subplot(g[0,3:])
    c=fig.add_subplot(g[1,:2]);d=fig.add_subplot(g[1,2:])
    e=fig.add_subplot(g[2,:3]);f=fig.add_subplot(g[2,3:])
    landscape(a,pt,st);patient_geometry(bp,pt);variance(c);raincloud(d,pt,st)
    subtype_ci(e,st);external_slope(f,ex)
    fig.subplots_adjust(left=.18,right=.975,top=.945,bottom=.055,wspace=.55,hspace=.74)
    for t in fig.findobj(match=plt.Text):
        if t.get_text().strip() and t.get_fontsize()<6.5:t.set_fontsize(6.5)
    fig.savefig(OUT/'Figure2.pdf');fig.savefig(OUT/'Figure2.png',dpi=300)
    plt.close(fig)

if __name__=='__main__':fig2()

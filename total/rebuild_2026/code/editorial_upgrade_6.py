"""Measured lineage and A7-C decomposition plate for Figure 6."""
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from matplotlib.patches import Rectangle

import editorial_rebuild_figures as b
from editorial_upgrade_347 import save, CMAP

D, DATA = b.D, b.DATA
Z1, Z2, ETP, GREY, LIGHT, DARK = b.Z1, b.Z2, b.ETP, b.GREY, b.LIGHT, b.DARK
read, panel, grid = b.read, b.panel, b.clean_grid


def _donor_heatmaps(fig, spec):
    don=read(D/'a7_lineage/scores/a7a_donor_min20.tsv')
    assert len(don)==12 and don.donor.is_unique
    inner=spec.subgridspec(1,2,wspace=.16)
    left,right=[fig.add_subplot(inner[0,i]) for i in range(2)]
    for ax,side,columns in [
            (left,'ZEB1',['HSPC_Z1','CLP_Z1','T_Z1','My_Z1']),
            (right,'ZEB2',['HSPC_Z2','CLP_Z2','T_Z2','My_Z2'])]:
        vals=don[columns].to_numpy(float)
        assert vals.shape==(12,4)
        ax.imshow(vals,aspect='auto',cmap=CMAP,vmin=-1.3,vmax=1.3,
                  interpolation='nearest',rasterized=True)
        ax.set_xticks(range(4),['HSPC','CLP','T','Myeloid'],fontsize=6.5)
        ax.set_yticks(range(12),don.donor if side=='ZEB1' else [],fontsize=6.5)
        ax.tick_params(length=0)
        for row in range(12):
            for col in range(4):
                if not np.isfinite(vals[row,col]):
                    ax.add_patch(Rectangle((col-.5,row-.5),1,1,facecolor='#E9EDEF',edgecolor='white',lw=.3))
    panel(left,'A','Normal-marrow donor × lineage matrix',
          '12 donors · left ZEB1 side / right ZEB2 side · grey = unavailable')


def _paired_donors(ax):
    don=read(D/'a7_lineage/scores/a7a_donor_min20.tsv')
    contrasts=[('Myeloid−T','d_My_minus_T_Z2',Z2),
               ('HSPC−T','d_HSPC_minus_T_Z2',ETP),
               ('CLP−T','d_CLP_minus_T_Z2',GREY)]
    for i,(label,col,color) in enumerate(contrasts):
        x=don[col].dropna().to_numpy()
        ax.scatter(x,np.full(len(x),i)+b.jitter(len(x),.06,i),s=18,color=color,alpha=.8)
        ax.plot([np.median(x)]*2,[i-.22,i+.22],color=DARK,lw=1.6)
    ax.axvline(0,color=LIGHT,lw=.6)
    ax.set_yticks([])
    ax.invert_yaxis();ax.set_xlim(-1.35,1.95)
    for i,label in enumerate(['Myeloid−T 12/12','HSPC−T 11/11','CLP−T 0/8']):
        ax.text(-1.31,i,label,va='center',fontsize=6.5,color=DARK)
    ax.set_xlabel('Within-donor ZEB2-side score difference')
    panel(ax,'B','Paired lineage contrasts','Positive = higher ZEB2-side score than T cells')
    grid(ax,'x')


def _gene_affinity(ax):
    x=read(D/'a7_lineage/scores/a7a_genes_min20.tsv').dropna(
        subset=['delta_Myeloid_minus_T','neg_beta'])
    assert len(x)==93
    ax.scatter(x.delta_Myeloid_minus_T,x.neg_beta,s=17,
               color=[Z2 if str(s).startswith('ZEB2') else Z1 for s in x.side],
               alpha=.7,edgecolors='white',linewidths=.2,rasterized=True)
    ax.axhline(0,color=LIGHT,lw=.5);ax.axvline(0,color=LIGHT,lw=.5)
    ax.set_xlabel('Normal myeloid − T effect')
    ax.set_ylabel('− pediatric residual coefficient βR')
    panel(ax,'C','Gene-level lineage affinity','93 genes · 75% expected direction · ρ = +0.60')
    grid(ax)


def _specimen_audit(ax):
    x=read(D/'a7_lineage/a7_0_frozen_genes_specimen_adjusted.tsv')
    assert len(x)==100 and np.sign(x.beta_primary).eq(np.sign(x.beta_specimen)).all()
    ax.scatter(x.beta_primary,x.beta_specimen,s=13,
               color=[Z2 if str(s).startswith('ZEB2') else Z1 for s in x.side],
               alpha=.72,rasterized=True)
    lim=max(abs(x.beta_primary).max(),abs(x.beta_specimen).max())*1.08
    ax.plot([-lim,lim],[-lim,lim],color=GREY,lw=.8,ls='--')
    ax.set_xlim(-lim,lim);ax.set_ylim(-lim,lim)
    ax.set_xlabel('Primary βR');ax.set_ylabel('Specimen-adjusted βR')
    panel(ax,'D','BM / PB specimen audit','100/100 same sign · ρ = .9995')
    grid(ax)


def _decomposition_heatmap(ax):
    split=read(D/'a7c_decompose/a7c_zeb2_split.tsv')
    assert len(split)==48 and split.symbol.is_unique
    assert split.subset.value_counts().to_dict()=={'myeloid_concordant':45,
                                                     'non_myeloid_concordant':3}
    names=['GSE146901','GSE234608','GSE243914','GSE248287']
    # Negative display values mean ZEB2-side-compatible in every leukemia column.
    sources=[(D/'a7c_decompose/a7c_GSE146901_genes.tsv','delta_external',-1),
             (D/'a7c_decompose/a7c_GSE234608_genes.tsv','delta_external',-1),
             (D/'a7b_malignant/a7b05_genes.tsv','delta',-1),
             (D/'a7b_malignant/a7b_genes.tsv','rho_E_vs_B',1)]
    split['order']=np.where(split.subset.eq('myeloid_concordant'),0,1)
    split=split.sort_values(['order','beta_R']).reset_index(drop=True)
    values=[]
    for path,column,sign in sources:
        x=read(path).set_index('symbol')[column]
        assert x.index.is_unique
        values.append(split.symbol.map(x).to_numpy(float)*sign)
    raw=np.column_stack(values)
    assert raw.shape==(48,4) and np.isfinite(raw).sum()>=185
    scales=np.nanquantile(np.abs(raw),.9,axis=0)
    disp=np.clip(raw/scales,-1,1)
    src=split[['symbol','subset','beta_R']].copy()
    for j,name in enumerate(names):
        src[name+'_ZEB2_compatible']=raw[:,j]
        src[name+'_scaled']=disp[:,j]
    src.to_csv(DATA/'editorial_visual_sources/F6E_A7C_48gene_contexts.tsv',
               sep='\t',index=False)
    ax.imshow(disp,aspect='auto',cmap=CMAP,vmin=-1,vmax=1,
              interpolation='nearest',rasterized=True)
    ax.axhline(44.5,color='white',lw=1.5)
    ax.set_xticks(range(4),['GSE146901','GSE234608','GSE243914','GSE248287'],fontsize=6.5)
    ax.set_yticks([])
    ax.tick_params(length=0)
    for i in range(48):
        ax.add_patch(Rectangle((-.68,i-.5),.10,1,
                     facecolor='#4E8277' if i<45 else DARK,edgecolor='none'))
    nonmyeloid=split.loc[45:,'symbol'].tolist()
    assert set(nonmyeloid)=={'AOAH','MAP3K5','ADRB2'}
    ax.text(4.12,24,'3 non-myeloid:\n'+' / '.join(nonmyeloid),
            ha='left',va='center',fontsize=6.5,color=DARK)
    ax.set_xlim(-.72,5.15)
    panel(ax,'E','A7-C decomposition across leukemia cohorts',
          '45 myeloid-affiliated + 3 named non-myeloid genes; blue = ZEB2-compatible')


def fig6():
    fig=b.mmfig(180,195)
    gs=GridSpec(3,6,figure=fig,height_ratios=[1.04,1.02,1.06])
    _donor_heatmaps(fig,gs[0,:3])
    bb=fig.add_subplot(gs[0,3:]);c=fig.add_subplot(gs[1,:3])
    dd=fig.add_subplot(gs[1,3:]);e=fig.add_subplot(gs[2,:])
    _paired_donors(bb);_gene_affinity(c);_specimen_audit(dd)
    _decomposition_heatmap(e)
    save(fig,6)


if __name__=='__main__':fig6()

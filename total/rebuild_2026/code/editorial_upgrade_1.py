"""Figure 1 coverage, stage and spatial redesign without inferred pseudotime."""
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from matplotlib.patches import Rectangle

import editorial_rebuild_figures as b
from editorial_upgrade_347 import save

D, SD = b.D, b.SD
Z1,Z2,ETP,GREY,LIGHT,DARK=b.Z1,b.Z2,b.ETP,b.GREY,b.LIGHT,b.DARK
read,panel,grid=b.read,b.panel,b.clean_grid


def _coverage(ax):
    tracks=read(SD/'F1/F1D_stage_tracks.tsv')
    states=read(D/'a4_yayon/a4a_stage_balance.tsv')
    spatial=read(D/'a4_yayon/a4b_tstate_cma_by_donor.tsv')
    assert tracks.dataset.nunique()==4 and len(states)==33
    assert spatial.state.nunique()==24 and spatial.donor_id.nunique()==6
    rows=['GSE142522','GSE195812','GSE206710','Park/HTA',
          'Yayon CITE-seq','Yayon spatial']
    cols=['ETP','DN early','DN(P)','DN(Q)','DP early','DP(Q)',
          'DP(P)','selection','CD4 SP','CD8 SP']
    # These are source-annotation coverage windows, not label-transfer calls.
    broad={
      'GSE142522':['ETP','DN early','DP early','selection','CD4 SP','CD8 SP'],
      'GSE195812':['ETP','DN early','DP early','selection','CD4 SP','CD8 SP'],
      'GSE206710':['DN early','DP early','selection'],
      'Park/HTA':['DN early','DP early','CD4 SP','CD8 SP'],
    }
    precise={
      'Yayon CITE-seq':cols,
      'Yayon spatial':['ETP','DN early','DN(P)','DN(Q)','DP early',
                       'DP(Q)','DP(P)','selection','CD4 SP','CD8 SP']}
    x=np.zeros((6,10),int)
    for i,row in enumerate(rows):
        for name in broad.get(row,[]):x[i,cols.index(name)]=1
        for name in precise.get(row,[]):x[i,cols.index(name)]=2
    assert x.shape==(6,10) and np.sum(x>0)==39
    shift=3.50
    ax.set_xlim(0,14.5);ax.set_ylim(-.7,6.35);ax.axis('off')
    panel(ax,'A','Atlas × developmental-window coverage',
          'Filled = annotated window; pale bulk windows are broad mappings')
    for i,row in enumerate(rows):
        y=5-i
        ax.text(.02,y,row,ha='left',va='center',fontsize=6.5,color=DARK)
        for j in range(10):
            col={0:'#FFFFFF',1:'#BFCBD0',2:'#607C8A'}[int(x[i,j])]
            ax.add_patch(Rectangle((shift+j-.43,y-.39),.86,.78,facecolor=col,
                                   edgecolor='#DAE1E4',lw=.35))
        ax.text(13.75,y,'bulk' if i<4 else 'CITE' if i==4 else 'spatial',
                va='center',fontsize=6.5,color=GREY)
    for j,label in enumerate(cols):
        ax.text(shift+j,-.57,label.replace(' ','\n'),rotation=45,ha='right',
                va='top',fontsize=6.5,color=DARK)
    for start,end,label,col in [(0,3,'ENTRY','#DCEAF0'),(4,7,'CORTEX','#F5E5D8'),
                                (8,9,'MEDULLA','#E9ECEE')]:
        ax.add_patch(Rectangle((shift+start-.44,5.68),end-start+.88,.38,
                               facecolor=col,edgecolor='white'))
        ax.text(shift+(start+end)/2,5.87,label,ha='center',va='center',
                fontsize=6.5,weight='bold',color=DARK)


def _atlas_trajectories(fig,spec):
    tracks=read(SD/'F1/F1D_stage_tracks.tsv')
    x=tracks[tracks.track.eq('Balance')]
    assert x.groupby('dataset').size().to_dict()=={
        'GSE142522':6,'GSE195812':8,'GSE206710':3,'Park/HTA':4}
    inner=spec.subgridspec(3,2,wspace=.24,hspace=.95,
                           height_ratios=[.18,1,1])
    header=fig.add_subplot(inner[0,:]);header.axis('off')
    header.text(0,.45,'B',fontsize=10.5,weight='bold',color=DARK)
    header.text(.09,.45,'Four thymus balance trajectories',fontsize=8.2,
                weight='bold',color=DARK)
    for k,dataset in enumerate(['GSE142522','GSE195812','GSE206710','Park/HTA']):
        ax=fig.add_subplot(inner[1+k//2,k%2]);q=x[x.dataset.eq(dataset)].sort_values('order')
        ax.plot(q.order,q.value,color=DARK,lw=.85,zorder=1)
        ax.scatter(q.order,q.value,c=[Z1 if v>=0 else Z2 for v in q.value],
                   s=17,zorder=2)
        ax.axhline(0,color=LIGHT,lw=.5)
        if k<2:
            ax.set_xticks(q.order,[])
        else:
            ax.set_xticks(q.order,[str(v).replace(' cortical',' cort.').replace('DP_','DP ') for v in q.stage],
                          rotation=28,ha='right',fontsize=6.5)
        ax.set_title(dataset,fontsize=7.0,color=DARK,weight='bold',loc='left')
        grid(ax,'y')


def _all_states_with_donor_inset(ax):
    x=read(D/'a4_yayon/a4a_stage_balance.tsv')
    assert len(x)==33 and x.stage.is_unique
    groups=['early_DN','cortical_DP','SP','other']
    x['group']=pd.Categorical(x.pole,categories=groups,ordered=True)
    x=x.sort_values(['group','balance_median']).reset_index(drop=True)
    assert x.group.value_counts().to_dict()=={
        'other':17,'early_DN':7,'cortical_DP':6,'SP':3}
    colors={'early_DN':Z2,'cortical_DP':Z1,'SP':GREY,'other':'#B6C1C5'}
    for i,r in x.iterrows():
        ax.plot([i,i],[0,r.balance_median],color=colors[str(r.group)],lw=.65,alpha=.65)
        ax.scatter(i,r.balance_median,s=np.clip(10+8*np.log10(r.n_cells),12,42),
                   color=colors[str(r.group)],edgecolors='white',lw=.25,zorder=2)
    ax.axhline(0,color=LIGHT,lw=.5)
    for boundary in [7,13,16]:ax.axvline(boundary-.5,color='#E6EAEC',lw=.6)
    for start,end,label in [(0,6,'EARLY / DN'),(7,12,'CORTICAL DP'),
                            (13,15,'SP'),(16,32,'OTHER STATES')]:
        ax.text((start+end)/2,-3.84,label,ha='center',fontsize=6.5,
                color=DARK,weight='bold')
    ax.set_xlim(-.6,32.7);ax.set_ylim(-4.0,4.3)
    ax.set_xticks([]);ax.set_ylabel('Stage median ZEB balance')
    panel(ax,'C','All 33 eligible author states',
          'Grouped by author pole; within-group position is NOT pseudotime')
    grid(ax,'y')
    labels=['T_ETP','T_DN(early)','T_DN(P)','T_DN(Q)',
            'T_DP(Q)-early','T_DP(Q)','T_DP(P)','T_CD4','T_CD8']
    for name in labels:
        q=x[x.stage.eq(name)]
        if len(q):
            i=int(q.index[0]);r=q.iloc[0]
            offset=(2,12) if name=='T_CD4' else (2,-14) if name=='T_CD8' else (3,4)
            ax.annotate(name.replace('T_',''),(i,r.balance_median),
                        xytext=offset,textcoords='offset points',fontsize=6.5,
                        color=DARK)
    # The donor comparison is a clearly readable inset in the unused upper-right
    # data region, not an independent schematic panel.
    d=ax.inset_axes([.69,.59,.29,.35],facecolor='white')
    dd=read(D/'a4_yayon/a4a_donor_level_balance.tsv')
    assert len(dd)==5
    for _,r in dd.iterrows():
        col=Z1 if r.DP_minus_early>0 else GREY
        d.plot([0,1],[r.early,r.cortical_DP],color=col,lw=.9)
        d.scatter([0,1],[r.early,r.cortical_DP],color=[Z2,Z1],s=10,zorder=2)
    d.set_xticks([0,1],['Early','DP'],fontsize=6.5)
    d.set_xlim(-.12,1.12);d.tick_params(labelsize=6.5,length=1)
    d.text(.02,.98,'D  4/5 donors',transform=d.transAxes,va='top',
           fontsize=6.5,weight='bold',color=DARK)
    d.spines[['top','right']].set_visible(False)


def _spatial_all_states(ax):
    x=read(D/'a4_yayon/a4b_tstate_cma_by_donor.tsv')
    assert x.state.nunique()==24 and x.groupby('state').donor_id.nunique().eq(6).all()
    order=x.groupby('state').cma_weighted_median.median().sort_values().index.tolist()
    for i,state in enumerate(order):
        q=x[x.state.eq(state)]
        col=Z2 if 'ETP' in state or 'DN' in state else Z1 if 'DP' in state else GREY
        ax.scatter(q.cma_weighted_median,np.full(len(q),i)+b.jitter(len(q),.045,i),
                   s=np.clip(8+4*np.log10(q.abundance_sum),10,26),
                   color=col,alpha=.55,rasterized=True)
        ax.plot([q.cma_weighted_median.median()]*2,[i-.20,i+.20],color=DARK,lw=1.5)
    ax.axvline(0,color=GREY,lw=.6)
    short=[s.removeprefix('T_').replace('_',' ').replace('(entry)',' entry')
           .replace('(agonist)',' agonist').replace('Treg-diff','Treg diff')
           for s in order]
    ax.set_yticks(range(24),short,fontsize=6.5);ax.invert_yaxis()
    ax.set_xlabel('Corticomedullary axis (CMA)')
    panel(ax,'E','All 24 spatially mapped T states',
          'Six donor estimates per state; dark tick = median')
    grid(ax,'x')


def _topology(ax):
    x=read(D/'a4_yayon/a4ab_stage_balance_and_cma.tsv').dropna(
        subset=['cma_weighted_median'])
    assert len(x)==19  # 19 of 33 RNA states have an exact spatial-name join.
    colors={'early_DN':Z2,'early':Z2,'cortical_DP':Z1,'DP':Z1,'SP':GREY}
    ax.scatter(x.cma_weighted_median,x.balance_median,
               s=np.clip(x.n_cells/230,12,85),
               color=[colors.get(v,'#BBC3C8') for v in x.pole],
               alpha=.72,edgecolors='white',linewidths=.3,rasterized=True)
    keys=['T_ETP','T_DN(early)','T_DP(Q)-early','T_DP(Q)','T_DP(P)','T_CD4']
    path=x.set_index('stage').reindex(keys).dropna(subset=['cma_weighted_median'])
    assert len(path)==len(keys)
    ax.plot(path.cma_weighted_median,path.balance_median,color=DARK,lw=.75,
            alpha=.55,zorder=1)
    for _,r in x[x.stage.isin(['T_ETP','T_DN(early)','T_DP(Q)','T_DP(P)','T_CD4','T_CD8'])].iterrows():
        label=r.stage.replace('T_','')
        offset=(-20,8) if r.stage=='T_CD4' else (-20,-9) if r.stage=='T_CD8' else (3,3)
        ax.annotate(label,(r.cma_weighted_median,r.balance_median),
                    xytext=offset,textcoords='offset points',fontsize=6.5,color=DARK)
    ax.axhline(0,color=LIGHT,lw=.5);ax.axvline(0,color=LIGHT,lw=.5)
    ax.set_xlabel('Median CMA position');ax.set_ylabel('Stage median ZEB balance')
    panel(ax,'F','Balance × spatial topology',
          'Path joins prespecified poles only; not an inferred trajectory')
    grid(ax)


def fig1():
    fig=b.mmfig(180,220)
    gs=GridSpec(3,6,figure=fig,height_ratios=[.95,1.18,1.70])
    a=fig.add_subplot(gs[0,:3]);c=fig.add_subplot(gs[1,:])
    e=fig.add_subplot(gs[2,:3]);f=fig.add_subplot(gs[2,3:])
    _coverage(a);_atlas_trajectories(fig,gs[0,3:])
    _all_states_with_donor_inset(c);_spatial_all_states(e);_topology(f)
    fig.subplots_adjust(left=.17,right=.97,top=.95,bottom=.06,wspace=.58,hspace=.62)
    # save() fixes the physical size and does not rescale any panel geometry.
    from editorial_upgrade_347 import OUT
    for t in fig.findobj(match=plt.Text):
        if t.get_text().strip() and t.get_fontsize()<6.5:t.set_fontsize(6.5)
    fig.savefig(OUT/'Figure1.pdf');fig.savefig(OUT/'Figure1.png',dpi=300)
    plt.close(fig)


if __name__=='__main__':fig1()

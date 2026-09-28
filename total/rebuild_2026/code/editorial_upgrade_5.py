"""Data-preserving Figure 5 redesign from frozen cross-cohort tables."""
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from matplotlib.patches import Rectangle
from matplotlib.colors import LinearSegmentedColormap

import editorial_rebuild_figures as b
from editorial_upgrade_347 import save, CMAP

D, SD, DATA = b.D, b.SD, b.DATA
Z1, Z2, ETP, GREY, LIGHT, DARK = b.Z1, b.Z2, b.ETP, b.GREY, b.LIGHT, b.DARK
read, panel, grid = b.read, b.panel, b.clean_grid


def _cohorts():
    x1=read(D/'gate3_validation/gse146901_program_scores.tsv')
    x1['sample']=x1['sample'].astype(str).str.zfill(3)
    orig=read(SD/'F5/F5A-C_GSE146901_samples.tsv')[['sample','balance']]
    orig['sample']=orig['sample'].astype(str).str.zfill(3)
    x1=x1.merge(orig,on='sample',validate='one_to_one').rename(columns={'balance':'B'})
    x2=read(D/'a5_gse234608/a5_patient_frozen_scores.tsv').rename(columns={'balance':'B'})
    x3=read(D/'a7b_malignant/a7b05_sample_scores.tsv')
    assert [len(x1),len(x2),len(x3)]==[18,28,24]
    return [(x1,'ETP','non-ETP'),(x2,'ETP_or_MPAL','conventional_TALL'),
            (x3,'ETP','nonETP')]


def _validation_matrix(ax):
    rows=[('GSE146901','diagnostic bulk','patient',18,84,.869,.596,'bulk'),
          ('GSE234608','ETP / MPAL bulk','patient',28,90,.800,.640,'bulk'),
          ('GSE243914','blast-enriched','specimen',24,98,.918,.770,'blast'),
          ('GSE248287','author malignant','patient',15,94,.777,.622,'malignant'),
          ('GSE280250','adult bulk','patient',79,86,.987,.860,'adult')]
    ax.set_xlim(0,10);ax.set_ylim(-.1,5.7);ax.axis('off')
    panel(ax,'A','Cross-cohort validation matrix',
          'One biological unit per row of each source cohort; fractions use tested genes')
    heads=[(.05,'Cohort'),(1.8,'Context'),(4.1,'Unit'),(5.15,'n'),
           (6.05,'Genes'),(7.25,'Direction'),(8.75,'Gene ρ')]
    for x,h in heads:ax.text(x,5.10,h,fontsize=6.5,weight='bold',color=DARK)
    seq=LinearSegmentedColormap.from_list('preserved',['#F1F5F5','#23658A'])
    for i,(cohort,context,unit,n,recovered,expected,rho,status) in enumerate(rows):
        y=4.28-i*.86
        ax.axhline(y-.38,color='#E5E9EB',lw=.45)
        ax.text(.05,y,cohort,va='center',fontsize=7,weight='bold',color=DARK)
        ax.text(1.8,y,context,va='center',fontsize=6.5,color=DARK)
        ax.text(4.1,y,unit,va='center',fontsize=6.5,color=GREY)
        ax.text(5.3,y,str(n),va='center',ha='center',fontsize=7,color=DARK)
        ax.scatter([6.35],[y],s=recovered*.58,facecolor='#AFC1CA',edgecolor=DARK,lw=.35)
        ax.text(6.72,y,str(recovered),va='center',fontsize=6.5,color=DARK)
        ax.add_patch(Rectangle((7.20,y-.28),.86,.56,facecolor=seq((expected-.70)/.30),edgecolor='white'))
        ax.text(7.63,y,f'{100*expected:.1f}%',ha='center',va='center',fontsize=6.5,
                color='white' if expected>.88 else DARK,weight='bold')
        ax.add_patch(Rectangle((8.68,y-.28),.70,.56,
                               facecolor=seq((rho-.50)/.40),edgecolor='white'))
        ax.text(9.03,y,f'{rho:.2f}',ha='center',va='center',fontsize=6.5,color=DARK)
        ax.add_patch(Rectangle((9.58,y-.28),.14,.56,
                  facecolor={'bulk':'#B9C6CB','blast':Z2,'malignant':ETP,'adult':DARK}[status],
                  edgecolor='none'))
    ax.text(9.48,5.10,'Assay',fontsize=6.5,weight='bold',color=DARK)


def _cohort_effect_matrix(ax, cohorts):
    metrics=[('ZEB1_side50','ZEB1 side'),('ZEB2_side50','ZEB2 side'),('B','Balance')]
    names=['GSE146901','GSE234608','GSE243914']
    effects=[]
    for row,(x,case,ref) in enumerate(cohorts):
        assert set(x.group.unique())=={case,ref}
        for col,(metric,label) in enumerate(metrics):
            a=x.loc[x.group.eq(ref),metric].median()
            c=x.loc[x.group.eq(case),metric].median()
            effects.append({'cohort':names[row],'metric':label,'median_reference_minus_ETP':a-c,
                            'n_reference':int(x.group.eq(ref).sum()),
                            'n_ETP':int(x.group.eq(case).sum())})
    tab=pd.DataFrame(effects)
    target=DATA/'editorial_visual_sources/F5B_descriptive_median_effects.tsv'
    target.parent.mkdir(parents=True,exist_ok=True);tab.to_csv(target,sep='\t',index=False)
    ax.set_xlim(-.65,2.55);ax.set_ylim(-.75,2.75)
    for row,name in enumerate(names):
        y=2-row
        for col,(metric,_) in enumerate(metrics):
            v=float(tab.loc[(tab.cohort==name)&(tab.metric==metrics[col][1]),
                            'median_reference_minus_ETP'].iloc[0])
            ax.add_patch(Rectangle((col-.43,y-.38),.86,.76,
                         facecolor=CMAP((np.clip(v,-1.5,1.5)+1.5)/3),edgecolor='white',lw=.5))
            ax.text(col,y,f'{v:+.2f}',ha='center',va='center',fontsize=7,weight='bold',color=DARK)
    ax.set_xticks(range(3),[x[1] for x in metrics],fontsize=6.5)
    ax.set_yticks([2,1,0],names,fontsize=6.5)
    ax.tick_params(length=0);ax.spines[['left','bottom']].set_visible(False)
    panel(ax,'B','Independent bulk / blast contrasts',
          'Cell = reference − ETP median; descriptive, no new test')


def _preservation_forest(ax):
    names=['GSE146901','GSE234608','GSE243914','GSE248287','GSE280250']
    recovered=np.array([84,90,98,94,86])
    frac=np.array([.869,.800,.918,.777,.987])
    rho=np.array([.596,.640,.770,.622,.860])
    ys=np.arange(5)
    cmap=LinearSegmentedColormap.from_list('concordant',['#DCE7E9',Z2])
    ax.scatter(rho,ys,s=recovered*.55,c=frac,cmap=cmap,vmin=.72,vmax=1,
               edgecolors=DARK,linewidths=.35,zorder=3)
    for y,x,p in zip(ys,rho,frac):
        ax.plot([.4,x],[y,y],color=LIGHT,lw=.8,zorder=1)
        ax.text(.98,y,f'{p*100:.1f}%',va='center',fontsize=6.5,color=DARK)
    ax.set_yticks([]);ax.invert_yaxis()
    ax.set_xlim(.18,1.11);ax.set_xticks([.4,.6,.8,1])
    for y,name in zip(ys,names):
        ax.text(.20,y,name,va='center',fontsize=6.5,color=DARK)
    ax.set_xlabel('Gene-level effect correlation ρ')
    panel(ax,'C','Constituent-gene preservation',
          'Dot area = recovered genes; shade = expected-direction fraction')
    grid(ax,'x')


def _malignant_scatter(fig, spec):
    gs=spec.subgridspec(1,2,wspace=.25)
    axes=[fig.add_subplot(gs[0,i]) for i in range(2)]
    mal=read(D/'a7b_malignant/a7b_patient_scores.tsv')
    summ=read(D/'a7b_malignant/a7b_summary.tsv').set_index('item')['value']
    assert len(mal)==15 and int(summ['n_patients'])==15
    for ax,metric,color,label,lo,hi,rho in [
            (axes[0],'ZEB1_side50',Z1,'ZEB1 side',summ['boot_Z1_lo'],summ['boot_Z1_hi'],summ['rho_Z1_vs_B']),
            (axes[1],'ZEB2_side50',Z2,'ZEB2 side',summ['boot_Z2_lo'],summ['boot_Z2_hi'],summ['rho_Z2_vs_B'])]:
        ax.scatter(mal.B,mal[metric],s=18,color=color,alpha=.82,rasterized=True)
        ax.axvline(0,color=LIGHT,lw=.5)
        ax.set_xlabel('Malignant-cell balance B')
        ax.set_ylabel('Frozen score' if ax is axes[0] else '')
        ax.text(.02,.96,f'{label}\nρ={rho:+.2f} [{lo:+.2f}, {hi:+.2f}]',
                transform=ax.transAxes,va='top',fontsize=6.5,color=DARK)
        grid(ax)
    panel(axes[0],'D','Malignant-cell pseudobulk',
          '15 independent patient-level units')


def _malignant_gene_heatmap(ax):
    x=read(D/'a7b_malignant/a7b_genes.tsv')
    assert len(x)==94 and set(x.side)=={'ZEB1-side50','ZEB2-side50'}
    x=x.sort_values(['side','rank'],ascending=[False,True]).reset_index(drop=True)
    raw=x[['beta_R','rho_E_vs_B']].to_numpy(float)
    scales=np.nanquantile(np.abs(raw),.9,axis=0)
    z=np.clip(raw/scales,-1,1)
    ax.imshow(z,aspect='auto',interpolation='nearest',cmap=CMAP,vmin=-1,vmax=1,
              rasterized=True)
    n_z2=int(x.side.eq('ZEB2-side50').sum())
    ax.axhline(n_z2-.5,color='white',lw=1.2)
    ax.set_xticks([0,1],['Discovery βR','Malignant ρ'],fontsize=6.5)
    ax.set_yticks([n_z2/2, n_z2+(94-n_z2)/2],
                  [f'Z2 {n_z2}',f'Z1 {94-n_z2}'],fontsize=6.5)
    ax.tick_params(length=0)
    panel(ax,'E','Gene-level malignant preservation',
          '94 recovered genes; columnwise display scaling')


def fig5():
    fig=b.mmfig(180,205)
    gs=GridSpec(3,6,figure=fig,height_ratios=[.78,.95,1.25])
    a=fig.add_subplot(gs[0,:]);bb=fig.add_subplot(gs[1,:3]);c=fig.add_subplot(gs[1,3:])
    e=fig.add_subplot(gs[2,3:])
    cohorts=_cohorts()
    _validation_matrix(a)
    _cohort_effect_matrix(bb,cohorts)
    _preservation_forest(c)
    _malignant_scatter(fig,gs[2,:3])
    _malignant_gene_heatmap(e)
    save(fig,5)


if __name__=='__main__':fig5()

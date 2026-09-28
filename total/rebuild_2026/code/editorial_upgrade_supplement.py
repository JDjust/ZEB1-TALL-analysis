"""Literal S3/S6/S9 visualization types from frozen source data."""
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from matplotlib.path import Path
from matplotlib.patches import PathPatch, Rectangle
from scipy.stats import spearmanr

import editorial_rebuild_figures as b
import editorial_rebuild_supplement as s

D,SD,DATA,OUT=b.D,b.SD,b.DATA,s.OUT
Z1,Z2,ETP,GREY,LIGHT,DARK=b.Z1,b.Z2,b.ETP,b.GREY,b.LIGHT,b.DARK
read,panel,grid=b.read,b.panel,b.clean_grid

def ribbon(ax,x0,x1,y0,y1,width,color,alpha=.45):
    dx=(x1-x0)*.46
    verts=[(x0,y0),(x0+dx,y0),(x1-dx,y1),(x1,y1),
           (x1,y1+width),(x1-dx,y1+width),(x0+dx,y0+width),(x0,y0+width),(x0,y0)]
    codes=[Path.MOVETO,Path.CURVE4,Path.CURVE4,Path.CURVE4,
           Path.LINETO,Path.CURVE4,Path.CURVE4,Path.CURVE4,Path.CLOSEPOLY]
    ax.add_patch(PathPatch(Path(verts,codes),facecolor=color,edgecolor='white',
                           lw=.14,alpha=alpha))

def alluvial(ax,records,layers,gap=.35,label_side=True):
    """Measured categorical ribbons; widths equal observed case counts."""
    n=len(records); assert n>0
    counts=[]; offsets=[]
    for col in layers:
        order=records[col].value_counts().index.tolist()
        vc=records[col].value_counts().reindex(order)
        starts={};z=0
        for key,count in vc.items(): starts[key]=z;z+=int(count)+gap
        counts.append(vc);offsets.append(starts)
    xmax=len(layers)-1
    for k in range(xmax):
        positions=[{key:0 for key in counts[k].index},
                   {key:0 for key in counts[k+1].index}]
        grouped=records.groupby([layers[k],layers[k+1]],sort=False).size()
        # Stable source/target ordering keeps the measured width exact.
        for a in counts[k].index:
            for c in counts[k+1].index:
                width=int(grouped.get((a,c),0))
                if not width:continue
                y0=offsets[k][a]+positions[0][a]
                y1=offsets[k+1][c]+positions[1][c]
                color=b.color_subtype(c if k==0 else a)
                if color==GREY:color='#8DA5AF'
                ribbon(ax,k+.045,k+.955,y0,y1,width,color)
                positions[0][a]+=width;positions[1][c]+=width
    for k,vc in enumerate(counts):
        for key,count in vc.items():
            y=offsets[k][key]
            ax.add_patch(Rectangle((k-.025,y),.05,int(count),facecolor=DARK,
                                   edgecolor='none',alpha=.78,zorder=3))
            textkey=b.short_subtype(str(key)).replace('Immature T-ALL (ETP-like)','ETP-like')
            if k==0:ax.text(k-.045,y+count/2,textkey,ha='right',va='center',fontsize=6.5)
            elif k==xmax:ax.text(k+.06,y+count/2,textkey,ha='left',va='center',fontsize=6.5)
            elif len(vc)<=3:ax.text(k+.06,y+count/2,textkey,ha='left',va='center',fontsize=6.5)
    total=n+max(len(q) for q in counts)*gap
    ax.set_xlim(-.7,xmax+.74);ax.set_ylim(total+.2,-.8)
    ax.set_xticks(range(len(layers)),[x.replace('_',' ') for x in layers],fontsize=6.7)
    ax.set_yticks([]);ax.spines[:].set_visible(False)
    return [dict(v) for v in counts]

def sc3():
    q=read(D/'aieop120_official/aieop120_vs_s9.tsv')
    assert len(q)==33 and q.exact_match.all() and q.sample_id.is_unique
    q=q.rename(columns={'fusion':'Fusion anchor',
                        'predicted_subtype_official120':'Official subtype'})
    f=b.mmfig(180,195);g=GridSpec(2,3,figure=f,height_ratios=[.75,1.65])
    a=f.add_subplot(g[0,0]);c=f.add_subplot(g[0,1]);d=f.add_subplot(g[0,2]);ab=f.add_subplot(g[1,:])
    p=read(D/'aieop120_official/aieop120_residual_with_official_labels.tsv')
    vc=p.predicted_subtype.value_counts().sort_values()
    a.barh(range(len(vc)),vc.values,color=GREY,height=.65)
    a.set_yticks(range(len(vc)),[b.short_subtype(v) for v in vc.index],fontsize=6.5)
    a.set_xlabel('Patients');panel(a,'A','Official call composition','120 cases; 16 classes')
    grid(a,'x')
    ex=read(D/'aieop120_official/aieop120_subtype_residual_medians.tsv')
    x=ex[ex.n_aieop>=5]
    for _,r in x.iterrows():
        c.plot([0,1],[r.residual_median_polonen,r.residual_median_aieop],
               color=b.color_subtype(r.predicted_subtype),lw=1)
        c.scatter([0,1],[r.residual_median_polonen,r.residual_median_aieop],
                  s=9,color=b.color_subtype(r.predicted_subtype))
    c.set_xticks([0,1],['Discovery','AIEOP']);c.set_ylabel('Median R')
    panel(c,'C','Six-class slopegraph','AIEOP n ≥ 5');grid(c,'y')
    d.hist(q.max_probability,bins=np.linspace(0,1,11),color=ETP,edgecolor='white')
    d.set_xlabel('Official max probability');d.set_ylabel('Fusion anchors')
    panel(d,'D','Assignment probability','33/33 labels matched');grid(d,'y')
    alluvial(ab,q,['Fusion anchor','Official subtype'],gap=.40)
    panel(ab,'B','Fusion-anchored classification alluvial',
          '33 measured cases; ribbon width = case count; source fusion is not a subtype call')
    f.subplots_adjust(left=.19,right=.96,top=.95,bottom=.065,wspace=.55,hspace=.68)
    _save(f,3)

def sc6():
    x=read(DATA/'editorial_visual_sources/F7G_frozen100_context_effects.tsv')
    cols=['Discovery βR','GSE146901 non−ETP','GSE234608 T−ETP',
          'GSE243914 non−ETP','Malignant ρ(E,B)','Adult βR']
    assert len(x)==100 and x.symbol.is_unique and all(c in x for c in cols)
    lab=['Discovery','GSE146901','GSE234608','GSE243914','Malignant','Adult']
    corr=np.eye(6);shared=np.zeros((6,6),int)
    for i in range(6):
        for j in range(6):
            v=x[[cols[i],cols[j]]].dropna();shared[i,j]=len(v)
            if i!=j and len(v)>=5:corr[i,j]=spearmanr(v.iloc[:,0],v.iloc[:,1]).statistic
    assert np.isfinite(corr).all()
    f=b.mmfig(180,170);g=GridSpec(2,3,figure=f,height_ratios=[1.15,.85])
    a=f.add_subplot(g[0,:]);axes=[f.add_subplot(g[1,i]) for i in range(3)]
    for i in range(6):
        for j in range(6):
            if i>j:
                col=plt.get_cmap('RdBu_r')((corr[i,j]+1)/2)
                a.add_patch(Rectangle((j-.48,i-.48),.96,.96,facecolor=col,edgecolor='white',lw=.7))
                a.text(j,i,f'{corr[i,j]:+.2f}',ha='center',va='center',fontsize=6.7,color=DARK)
            elif i<j:
                a.text(j,i,str(shared[i,j]),ha='center',va='center',fontsize=6.5,color=GREY)
            else:
                a.add_patch(Rectangle((j-.48,i-.48),.96,.96,facecolor='#EDF0F2',edgecolor='white'))
                a.text(j,i,'—',ha='center',va='center',fontsize=6.7)
    a.set_xlim(-.5,5.5);a.set_ylim(5.5,-.5)
    a.set_xticks(range(6),lab,rotation=30,ha='right',fontsize=6.5)
    a.set_yticks(range(6),lab,fontsize=6.5);a.tick_params(length=0)
    a.spines[:].set_visible(False)
    panel(a,'A','Frozen-gene effect-vector concordance',
          'Lower triangle: pairwise Spearman ρ; upper triangle: shared measured genes')
    for ax,letter,idx in zip(axes,'BCD',[3,4,5]):
        v=x[[cols[0],cols[idx],'side']].dropna()
        for side,col in [('ZEB1-side50',Z1),('ZEB2-side50',Z2)]:
            t=v[v.side.eq(side)]
            ax.scatter(t.iloc[:,0],t.iloc[:,1],s=8,color=col,alpha=.72)
        ax.axhline(0,color=LIGHT,lw=.5);ax.axvline(0,color=LIGHT,lw=.5)
        ax.set_xlabel('Discovery βR');ax.set_ylabel(lab[idx]+' effect')
        panel(ax,letter,f'Discovery × {lab[idx]}',f'{len(v)} shared genes')
        grid(ax)
    f.subplots_adjust(left=.15,right=.975,top=.94,bottom=.085,wspace=.44,hspace=.68)
    _save(f,6)

def sc9():
    x=read(D/'a8_adult/a8_2_merged.tsv')
    assert len(x)==79 and x.geo_accession.is_unique
    col='T-ALL main-cluster high-confidence'
    calls=x[col].fillna('No HC main').astype(str)
    first=calls.str.extract(r'^(C\d+(?:\.\d+)?)',expand=False).fillna('No HC main')
    x['Main cluster']=np.where(calls.str.contains(';',regex=False),'Multiple HC',first)
    x['Immature call']=np.where(x['T-ALL immature high-confidence'].notna(),'ETP-like','Other')
    x['Residual display']=pd.cut(x.R,[-np.inf,-.5,.5,np.inf],
                                 labels=['R < −0.5','−0.5 to 0.5','R > 0.5'])
    assert x['Residual display'].notna().all()
    f=b.mmfig(180,195);g=GridSpec(2,3,figure=f,height_ratios=[1.72,.8])
    a=f.add_subplot(g[0,:]);bb=f.add_subplot(g[1,0]);cc=f.add_subplot(g[1,1]);dd=f.add_subplot(g[1,2])
    alluvial(a,x,['Main cluster','Immature call','Residual display'],gap=.5)
    panel(a,'A','Adult subtype-to-state alluvial',
          '79 cases; residual bins use ±0.5 for visualization only, not classification')
    vc=x['Main cluster'].value_counts().sort_values()
    bb.barh(range(len(vc)),vc.values,color=GREY)
    bb.set_yticks(range(len(vc)),vc.index,fontsize=6.5);bb.set_xlabel('Patients')
    panel(bb,'B','High-confidence clusters','Multiple or absent calls retained');grid(bb,'x')
    cc.scatter(x.age,x.R,s=13,color=[Z2 if r<-.5 else Z1 if r>.5 else GREY for r in x.R],alpha=.8)
    cc.axvline(18,color=LIGHT,lw=.6);cc.axhline(0,color=LIGHT,lw=.6)
    cc.set_xlabel('Age at diagnosis');cc.set_ylabel('Adult-cohort residual R')
    panel(cc,'C','Age and residual','79 cases; descriptive only');grid(cc)
    med=x.groupby('Immature call').R.median()
    for i,(name,v) in enumerate(med.items()):
        vals=x.loc[x['Immature call'].eq(name),'R']
        dd.scatter(vals,np.full(len(vals),i)+b.jitter(len(vals),.07,i),s=8,color=ETP if name=='ETP-like' else GREY,alpha=.6)
        dd.plot([v,v],[i-.18,i+.18],color=DARK,lw=1.6)
    dd.set_yticks(range(len(med)),med.index);dd.axvline(0,color=LIGHT,lw=.5)
    dd.set_xlabel('Patient R');panel(dd,'D','Immature-call context','Patient values and medians');grid(dd,'x')
    f.subplots_adjust(left=.15,right=.97,top=.865,bottom=.07,wspace=.58,hspace=.43)
    _save(f,9)

def _save(f,n):
    for t in f.findobj(match=plt.Text):
        if t.get_text().strip() and t.get_fontsize()<6.5:t.set_fontsize(6.5)
    f.savefig(OUT/f'FigureS{n}.pdf')
    f.savefig(OUT/f'FigureS{n}.png',dpi=300)
    plt.close(f)

if __name__=='__main__':
    sc3();sc6();sc9()

from plot_common import *
def effect_scatter(a,fr,col,symlog=False):
    q=fr[['Discovery βR',col,'side']].dropna()
    for side in ['ZEB1-side50','ZEB2-side50']:
        z=q[q.side==side];a.scatter(z.iloc[:,0],z.iloc[:,1],color=color(side),s=16,alpha=.8,rasterized=True)
    zeros(a);a.set_xlabel('Discovery βR');a.set_xlim(-1.16,1.16)
    if symlog:
        a.set_yscale('symlog',linthresh=2);a.set_yticks([-100,-10,0,10,100],['−100','−10','0','10','100']);bound=q.iloc[:,1].abs().max()*1.5;a.set_ylim(-bound,bound)
    else:a.margins(y=.25)
    a.text(.98,.97,f'n={len(q)} genes',transform=a.transAxes,ha='right',va='top');grid(a)
def main():
    f=fig(211);fr=read(V/'F7G_frozen100_context_effects.tsv');assert len(fr)==100
    cols=['Discovery βR','GSE146901 non−ETP','GSE234608 T−ETP','GSE243914 non−ETP','Malignant ρ(E,B)','Adult βR'];labels=['Discovery','GSE146901','GSE234608','GSE243914','Malignant','Adult']
    a=ax(f,box(f,34,9,131,34),'A','Pairwise frozen-gene correlations (ρ / shared gene n)')
    # Same pairwise Spearman display as the prior plot; no inferential test.
    corr=fr[cols].corr(method='spearman').to_numpy();masked=np.ma.masked_where(np.triu(np.ones((6,6)),0).astype(bool),corr)
    im=a.imshow(masked,cmap=LinearSegmentedColormap.from_list('seq',['white',Z2]),vmin=0,vmax=1,aspect='auto')
    for i in range(6):
        for j in range(i):
            n=len(fr[[cols[i],cols[j]]].dropna());a.text(j,i,f'{corr[i,j]:.2f} / {n}',ha='center',va='center',color='white' if corr[i,j]>.75 else INK)
    a.set_xticks(range(6),['Disc.','146901','234608','243914','Malignant','Adult']);a.set_yticks(range(6),labels);a.tick_params(length=0);a.spines[:].set_visible(False)

    scale(f,im,box(f,137,55,28,1.2),'Spearman ρ',[0,.5,1])
    settings=[('B','GSE146901',cols[1],'GSE146901: non-ETP − ETP\nFPKM (symlog)',True),('C','GSE234608',cols[2],'GSE234608: T − ETP/MPAL\nlog2(CPM+1)',False),('D','GSE243914: blast-enriched',cols[3],'GSE243914: non-ETP − ETP\nlog2(TPM+1)',False),('E','GSE248287: malignant',cols[4],'GSE248287: ρ(gene,B)',False)]
    for i,(l,title,col,ylabel,symlog) in enumerate(settings):
        a=ax(f,box(f,26+(i%2)*88,73+(i//2)*47,60,28),l,title);effect_scatter(a,fr,col,symlog);a.set_ylabel(ylabel)
    mal=read(D/'a7b_malignant/a7b_patient_scores.tsv');assert len(mal)==15 and mal.n_cells.sum()==29454
    for i,(metric,co,rho,ci) in enumerate([('ZEB1_side50',Z1,'+.73',None),('ZEB2_side50',Z2,'−.21','[−.70, +.36]')]):
        a=ax(f,box(f,26+i*88,170,60,29),['F','G'][i],['Malignant ZEB1-side score','Malignant ZEB2-side score'][i]);a.scatter(mal.B,mal[metric],color=co,s=25)
        a.set(xlabel='Malignant patient balance B',ylabel=('ZEB1' if i==0 else 'ZEB2')+'-side score');a.axvline(0,color=LIGHT,lw=.5);grid(a)
        a.text(.03,.98,f'n=15; ρ={rho}'+('\n95% CI '+ci if ci else ''),transform=a.transAxes,va='top')
        a.margins(y=.50)
    f.legend(handles=[Line2D([],[],marker='o',ls='',color=Z2,label='ZEB2 side'),Line2D([],[],marker='o',ls='',color=Z1,label='ZEB1 side')],loc='center',bbox_to_anchor=anchor(f,63,61),ncol=2,frameon=False)
    save(f,5)
if __name__=='__main__':main()

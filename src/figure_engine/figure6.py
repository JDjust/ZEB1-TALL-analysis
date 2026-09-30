from plot_common import *
def main():
    f=fig(170);u=read(D/'a7_lineage/scores/a7a_unit_scores_min20.tsv');order=['T lineage','lymphoid progenitor/CLP','HSC/MPP/HSPC','Myeloid'];rng=np.random.default_rng(16)
    for j,(metric,co) in enumerate([('ZEB1_side50',Z1),('ZEB2_side50',Z2)]):
        a=ax(f,box(f,24+j*89,10,61,36),'A' if j==0 else '',['Normal marrow: ZEB1 side','Normal marrow: ZEB2 side'][j])
        for i,lin in enumerate(order):
            vals=u.loc[u.lineage==lin,metric];a.scatter(i+rng.uniform(-.12,.12,len(vals)),vals,color=co,s=20,alpha=.8);a.plot([i-.2,i+.2],[vals.median()]*2,color=INK,lw=1.2)
        a.set_xticks(range(4),['T','CLP','HSPC','Myeloid']);a.set_ylabel(('ZEB1' if j==0 else 'ZEB2')+'-side\ndonor pseudobulk score');grid(a)

    b=ax(f,box(f,24,77,58,52),'B','Gene-level myeloid affinity');ge=read(D/'a7_lineage/scores/a7a_genes_min20.tsv').dropna(subset=['delta_Myeloid_minus_T','neg_beta']);assert len(ge)==93
    b.scatter(ge.delta_Myeloid_minus_T,ge.neg_beta,c=ge.side.map(color),s=20,alpha=.8);zeros(b);b.set(xlabel='Normal myeloid − T\n(logCPM effect)',ylabel='− discovery coefficient βR');b.text(.02,.97,'n=93 genes; ρ=.60',va='top',transform=b.transAxes);b.margins(y=.25);grid(b)

    c=ax(f,box(f,121,77,53,58),'C','Frozen 45 + 3 genes across leukemia contexts')
    detail=read(V/'F6E_A7C_48gene_contexts.tsv');assert len(detail)==48 and detail.symbol.is_unique
    cols=['GSE146901','GSE234608','GSE243914','GSE248287'];mat=detail[[k+'_scaled' for k in cols]].to_numpy();edges=np.r_[np.arange(46,dtype=float),48,51,54]
    im=c.pcolormesh(np.arange(5)-.5,edges,mat,cmap=CMAP,vmin=-1,vmax=1,edgecolors='none',rasterized=True);c.set(ylim=(54,0),xlim=(-.5,3.5));c.axhline(45,color='white',lw=2)
    c.set_yticks([22.5,47.5,50.5,53],['45 myeloid-\nconcordant','AOAH','MAP3K5','ADRB2']);c.set_xticks(range(4),['146901\nnon-ETP\n− ETP','234608\nT − ETP/\nMPAL','243914\nnon-ETP\n− ETP','248287\nρ(gene,B)'])
    c.tick_params(length=0);c.spines[:].set_visible(False)
    scale(f,im,box(f,113,162,55,1.2),'',[-1,0,1])
    swatchkey(f,[(NA,'Unavailable')],anchor(f,151,156),1)
    sidekey(f,1-64/170)
    save(f,6)
if __name__=='__main__':main()

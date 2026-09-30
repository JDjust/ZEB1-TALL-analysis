from plot_common import *
def main():
    f=fig(210);a=ax(f,box(f,22,10,64,34),'A','Genome-wide residual effects')
    ge=read(D/'gate2_program/all_genes_residual_effects.tsv');ge=ge[~ge.symbol.str.upper().isin(['ZEB1','ZEB2','CD34','LYL1','CD1A'])];fr=read(D/'gate2_program/frozen_ZEB_side50.tsv');assert len(ge)==24356 and len(fr)==100
    a.scatter(ge.beta_R,ge.partial_R2,s=2,color='#D9D9D9',rasterized=True)
    for side in ['ZEB1-side50','ZEB2-side50']:
        q=fr[fr.side==side];a.scatter(q.beta_R,q.partial_R2,s=9,color=color(side),rasterized=True)
    a.set(xlabel='Leukemia residual coefficient βR',ylabel='Ranking score');grid(a)

    b=ax(f,box(f,30,70,140,45),'C','Frozen expression across patient residual deciles')
    x=read(V/'F3C_frozen100_residual_deciles.tsv');x['side_order']=np.where(x.side=='ZEB2-side50',0,1);x=x.sort_values(['side_order','rank']).reset_index(drop=True)
    mat=x[[f'decile_{i}' for i in range(1,11)]].to_numpy();im=b.imshow(mat,aspect='auto',extent=(0,10,100,0),vmin=-1.25,vmax=1.25,cmap=CMAP,interpolation='nearest',rasterized=True)
    ann=read(V/'F3C_residual_decile_annotations.tsv').sort_values('decile');normal=read(D/'a6_yayon/a6_2_gene_effects.tsv').set_index('symbol')['class'];my=read(D/'a7c_decompose/a7c_zeb2_split.tsv').set_index('symbol')['subset']
    nc={'concordant':Z1,'discordant':Z2,'neutral':GREY};mc={'myeloid_concordant':ETP,'non_myeloid_concordant':INK}
    for i,r in x.iterrows():
        for xx,co in [(-.23,color(r.side)),(-.50,nc.get(normal.get(r.symbol),NA)),(-.77,mc.get(my.get(r.symbol),NA))]:b.add_patch(Rectangle((xx,i),.2,1,color=co,lw=0))
    for i,r in enumerate(ann.itertuples()):
        left=i
        for pct,co in [(r.pct_BCL11B,Z2),(r.pct_ETP_like,ETP),(r.pct_TLX3,Z1),(100-r.pct_BCL11B-r.pct_ETP_like-r.pct_TLX3,NA)]:
            b.add_patch(Rectangle((left,-18),pct/100,3.8,color=co,lw=0));left+=pct/100
        for y,v in [(-12,r.median_R),(-6,r.median_D)]:b.add_patch(Rectangle((i,y),1,3.8,color=CMAP((np.clip(v,-2,2)+2)/4),lw=0))
    for y,label in [(-16,'Subtype'),(-10,'Median R'),(-4,'Median D')]:b.text(-.9,y,label,ha='right',va='center')
    b.set(xlim=(-.85,10),ylim=(100,-20));b.axhline(50,color='white',lw=1.2);b.set_yticks([25,75],['ZEB2\n50 genes','ZEB1\n50 genes']);b.set_xticks(np.arange(.5,10.5),range(1,11));b.set_xlabel('Patient residual decile');b.tick_params(length=0);b.spines[:].set_visible(False)
    scale(f,im,box(f,137,129,33,1.2),'Median expression z',[-1,0,1])
    swatchkey(f,[(Z2,'BCL11B'),(ETP,'ETP-like'),(Z1,'TLX3'),(NA,'Other subtype')],anchor(f,100,62),4)
    f.text(*anchor(f,7,131),'Myeloid | Development | Side',fontsize=8)
    swatchkey(f,[(Z1,'Concordant'),(Z2,'Discordant'),(GREY,'Smaller')],anchor(f,65,139),3)
    swatchkey(f,[(ETP,'Myeloid concordant'),(INK,'Nonconcordant'),(NA,'NA')],anchor(f,75,146),3)
    swatchkey(f,[(Z2,'ZEB2 side'),(Z1,'ZEB1 side')],anchor(f,92,3),2)
    c=ax(f,box(f,118,10,55,34),'B','Normal trajectories (n=5 donors)');n=read(D/'a6_yayon/a6_1_donor_program.tsv')
    for side,co in [('ZEB1',Z1),('ZEB2',Z2)]:
        keys=[f'{stage}_{side}_side50' for stage in ['early','DP','SP']]
        for _,r in n.iterrows():c.plot(range(3),r[keys].astype(float),color=co,alpha=.3,lw=.7)
        c.plot(range(3),n[keys].median().to_numpy(),color=co,lw=1.8,marker='o',ms=3)
    c.set_xticks(range(3),['Early','DP','SP']);c.set_ylabel('Frozen program score');c.text(.98,.97,'n=5 donors',ha='right',va='top',transform=c.transAxes);grid(c)
    d=ax(f,box(f,23,156,58,28),'D','Normal versus leukemia effects');q=read(D/'a6_yayon/a6_2_gene_effects.tsv');assert q['class'].value_counts().to_dict()=={'neutral':45,'discordant':32,'concordant':13}
    for cls,marker in [('concordant','o'),('discordant','x'),('neutral','.')]:
        z=q[q['class']==cls];d.scatter(z.delta_normal,z.beta_R,c=z.side.map(color),marker=marker,s=18,rasterized=True)
    zeros(d);d.set(xlabel='Normal DP − early\n(atlas expression)',ylabel='Leukemia βR');d.text(.98,.97,'n=90 genes; ρ=−.50',ha='right',va='top',transform=d.transAxes)
    f.legend(handles=[Line2D([],[],ls='',marker=m,color=INK,label=t) for m,t in [('o','13 concordant'),('x','32 discordant'),('.','45 smaller')]],loc='center',bbox_to_anchor=anchor(f,90,205),ncol=3,frameon=False)
    z=q[q.side=='ZEB2-side50'].merge(read(SD/'revision4_targeted/matched_BCL11B_ETP_gene_results.tsv')[['symbol','logFC']],on='symbol',validate='one_to_one')
    for j,(mask,title,pv) in enumerate([(z['class']!='concordant','38 discordant or smaller','3.9e-16'),(z['class']=='concordant','10 normal-early-concordant','.019')]):
        e=ax(f,box(f,114,154+j*24,59,15),'E' if j==0 else '',title);v=z[mask].sort_values('logFC');assert len(v)==[38,10][j]
        e.hlines(range(len(v)),0,v.logFC,color=Z2 if j==0 else GREY,lw=.5);e.scatter(v.logFC,range(len(v)),s=6,color=Z2 if j==0 else GREY);e.set(xlim=(-1.05,4.2),ylim=(-1,len(v)*1.35),yticks=[],xlabel='BCL11B − ETP-like log2 effect' if j==1 else '');e.axvline(0,color=LIGHT,lw=.5);e.text(.02,.96,('n=38' if j==0 else 'n=10')+'; P='+pv,transform=e.transAxes,va='top')

    save(f,3)
if __name__=='__main__':main()

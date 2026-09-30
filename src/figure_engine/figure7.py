from plot_common import *
def main():
    f=fig(229);pt=read(D/'a8_adult/a8_patient_scores.tsv');assert len(pt)==79
    a=ax(f,box(f,18,10,44,30),'A','Adult B–D fit');a.scatter(pt.D,pt.B,color=GREY,s=9);p=pt.sort_values('D');a.plot(p.D,p.B-p.R,color=INK,lw=1);a.set(xlabel='Adult proxy D',ylabel='Adult balance B');grid(a)
    for j,(metric,co,title,ci) in enumerate([('ZEB1_side50',Z1,'ZEB1 side','ρ=+.60; CI .43–.73'),('ZEB2_side50',Z2,'ZEB2 side','ρ=−.66; CI −.76–−.52')]):
        a=ax(f,box(f,74+j*55,10,44,30),['B','C'][j],title);a.scatter(pt.R,pt[metric],color=co,s=10);a.set(xlabel='Adult residual R',ylabel='Side score');a.axvline(0,color=LIGHT,lw=.5);grid(a);f.text(*anchor(f,73+j*55,54),ci.replace('; CI','\nCI'),fontsize=8)
    d=ax(f,box(f,23,72,55,32),'D','Adult constituent-gene effects');ge=read(D/'a8_adult/a8_genes.tsv').dropna(subset=['beta_R','beta_adult']);assert len(ge)==78
    d.scatter(ge.beta_R,ge.beta_adult,c=ge.side.map(color),s=14);zeros(d);d.set(xlabel='Discovery βR',ylabel='Adult βR');d.margins(y=.25);d.text(.02,.97,'77/78; ρ=.86',transform=d.transAxes,va='top');grid(d)
    e=ax(f,box(f,134,72,39,75),'E','Adult clusters');ad=read(D/'a8_adult/a8_2_merged.tsv');maincol='T-ALL main-cluster high-confidence';imm='T-ALL immature high-confidence'
    def group(v):
        if pd.isna(v):return 'No HC main'
        if ';' in str(v):return 'Multiple HC'
        c=str(v).split(' ',1)[0]
        return 'Other HC' if c in ['C9','C14','C15'] else {'C8':'C8 BCL11B','C2':'C2 TAL1 DP','C17':'C17 TLX3'}.get(c,c)
    ad['group']=ad[maincol].map(group);stats=ad.groupby('group').agg(n=('R','size'),median=('R','median'),immature=(imm,lambda x:x.notna().mean())).sort_values('median')
    assert len(stats)==11 and stats.n.sum()==79
    for i,(name,r) in enumerate(stats.iterrows()):
        vals=ad.loc[ad.group==name,'R'];e.scatter(vals,i+np.random.default_rng(i).uniform(-.13,.13,len(vals)),s=8,color=color(name));e.plot([r['median']-.09,r['median']+.09],[i,i],color=INK,lw=1.5);e.text(5.1,i,f'{100*r.immature:.0f}',va='center',ha='center')
    e.set_yticks(range(11),[f'{name} ({int(r.n)})' for name,r in stats.iterrows()]);e.set(ylim=(10.6,-.6),xlim=(-3.5,6),xlabel='Adult cohort residual R');e.axvline(0,color=LIGHT,lw=.5);e.text(5.1,-1.1,'Immature\ncalls (%)',ha='center',va='bottom',fontsize=8)
    a=ax(f,box(f,36,127,50,20),'F','Prespecified subtype anchors')
    summ=read(D/'a8_adult/a8_2_summary.tsv');s=summ[summ.item=='primary_rank_point'].set_index('label');dis=read(SD/'F3/F3D-F_subtype_statistics.tsv').set_index('subtype')
    anchors=[('BCL11B','C8 (BCL11B)'),('ETP-like','Immature T-ALL (ETP-like)'),('TAL1 DP','C2 (TAL1 DP-like)')];rows=[]
    for name,key in anchors:
        r=s.loc[key];rows.append((name,int(r.n),r.disc_R_median,r.adult_R_median))
    tlx=ad[ad[maincol].astype(str).str.startswith('C17 ')];rows.append(('TLX3',2,dis.loc['TLX3','residual_within_median'],tlx.R.median()))
    for i,(name,n,v1,v2) in enumerate(rows):
        a.plot([v1,v2],[i,i],color=LIGHT,lw=.8);a.scatter(v1,i,color=GREY,s=17);a.scatter(v2,i,facecolor='white' if n==2 else color(name),edgecolor=color(name),s=23)
    a.set_yticks(range(4),[f'{x[0]} (n={x[1]})' for x in rows]);a.set(ylim=(3.5,-.5),xlabel='Cohort-specific median R');a.axvline(0,color=LIGHT,lw=.5)
    h=ax(f,box(f,30,160,140,29),'G','Frozen-gene context effects');x=read(V/'F7G_frozen100_context_effects.tsv');assert x.side.iloc[:50].eq('ZEB2-side50').all()
    mat=x[[f'scaled_{i}' for i in range(1,9)]].to_numpy();im=h.imshow(mat,aspect='auto',extent=(0,8,100,0),cmap=CMAP,vmin=-1,vmax=1,interpolation='nearest',rasterized=True)
    norm=read(D/'a6_yayon/a6_2_gene_effects.tsv').set_index('symbol')['class'];my=read(D/'a7c_decompose/a7c_zeb2_split.tsv').set_index('symbol')['subset']
    for i,r in x.iterrows():
        for xx,co in [(-.19,color(r.side)),(-.39,{'concordant':Z1,'discordant':Z2,'neutral':GREY}.get(norm.get(r.symbol),NA)),(-.59,{'myeloid_concordant':ETP,'non_myeloid_concordant':INK}.get(my.get(r.symbol),NA))]:h.add_patch(Rectangle((xx,i),.14,1,color=co,lw=0))
    for i in range(8):
        h.add_patch(Rectangle((i,-10),1,7,color=GREY if i in [1,2] else INK,lw=0));h.text(i+.5,-6.5,'N' if i in [1,2] else 'L',color='white',ha='center',va='center',fontsize=8)
    h.set(xlim=(-.65,8),ylim=(100,-11));h.axhline(50,color='white',lw=1);h.set_yticks([25,75],['ZEB2\n50 genes','ZEB1\n50 genes']);h.set_xticks(np.arange(.5,8.5),['Discovery βR','Thymus DP−early','Marrow My−T','146901 non−ETP','234608 T−ETP','243914 non−ETP','Malignant ρ(gene,B)','Adult βR'],rotation=42,ha='right');h.tick_params(length=0);h.spines[:].set_visible(False)
    swatchkey(f,[(GREY,'N: normal'),(INK,'L: leukemia'),(Z2,'ZEB2 side'),(Z1,'ZEB1 side')],anchor(f,95,213),4)
    swatchkey(f,[(Z1,'Development: concordant'),(Z2,'Discordant'),(GREY,'Smaller')],anchor(f,95,218),3)
    swatchkey(f,[(ETP,'Myeloid: concordant'),(INK,'Nonconcordant'),(NA,'Unavailable')],anchor(f,90,223),3)
    scale(f,im,box(f,144,221,26,1),'',[-1,0,1])
    key=f.legend(handles=[Line2D([],[],marker='o',ls='',color=GREY,label='Discovery'),Line2D([],[],marker='o',ls='',color=Z2,label='Adult (subtype color)')],loc='center',bbox_to_anchor=anchor(f,60,157),ncol=2,frameon=False);f.add_artist(key)
    save(f,7)
if __name__=='__main__':main()

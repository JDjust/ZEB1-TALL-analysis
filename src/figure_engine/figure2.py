from plot_common import *
def main():
    pt=read(SD/'F3/F3B-C_patient_level.tsv');st=read(SD/'F3/F3D-F_subtype_statistics.tsv').sort_values('residual_within_median')
    assert len(pt)==1309 and len(st)==17
    f=fig(154)
    a=ax(f,box(f,20,10,56,53),'A','Cohort-internal balance fit')
    a.scatter(pt.dev,pt.balance,s=4,color='#D9D9D9',rasterized=True)
    for name in ['BCL11B','ETP-like','TLX3']:
        q=pt[pt.subtype==name];a.scatter(q.dev,q.balance,s=10,color=color(name),label=name,rasterized=True)
    curve=read(SD/'F3/F3B_within_cohort_curve.tsv');a.plot(curve.dev,curve.expected_within,color=INK,lw=1.1,label='Spline df=3')
    a.set(xlabel='Developmental proxy D',ylabel='ZEB balance B');a.legend(frameon=False,ncol=1,loc='lower right',handlelength=1,fontsize=8);grid(a)
    b=ax(f,box(f,134,10,40,78),'B','Patient residuals and subtype medians')
    rng=np.random.default_rng(20260928)
    for i,r in enumerate(st.itertuples()):
        q=pt.loc[pt.subtype==r.subtype,'residual_within'];c=color(r.subtype)
        b.scatter(q,i+rng.uniform(-.22,.22,len(q)),s=4,color=c,alpha=.4,rasterized=True)
        b.plot([r.residual_within_ci_low,r.residual_within_ci_high],[i,i],lw=1.7,color=INK,zorder=4)
        b.scatter(r.residual_within_median,i,marker='D',s=20,facecolor='white',edgecolor=INK,zorder=5)
    b.set_yticks(range(17),[f'{short(s)} (n={int(n)})' for s,n in zip(st.subtype,st.n)]);b.set_ylim(16.6,-.6);b.axvline(0,color=LIGHT,lw=.6);b.set_xlabel('Discovery cohort residual R');grid(b,'x')
    c=ax(f,box(f,27,95,49,35),'C','Explained balance variation')
    c.barh([1,0],[25.6,39.5],color=[GREY,INK],height=.5);c.set_yticks([1,0],['D only','D + subtype']);c.set_xlim(0,48);c.set_xlabel('Raw in-sample R² (%)')
    for y,v in [(1,25.6),(0,39.5)]:c.text(v+1,y,str(v),va='center')
    c.text(0,.97,'+13.9 percentage points',transform=c.transAxes,va='top');grid(c,'x')
    e=ax(f,box(f,134,110,40,30),'D','External subtype medians')
    ex=read(D/'aieop120_official/aieop120_subtype_residual_medians.tsv');ex=ex[(ex.n_aieop>=5)|(ex.predicted_subtype=='BCL11B')].sort_values('residual_median_aieop')
    for i,r in enumerate(ex.itertuples()):
        co=color(r.predicted_subtype);e.plot([r.residual_median_polonen,r.residual_median_aieop],[i,i],color=co,ls='--' if r.n_aieop==4 else '-',lw=1)
        e.scatter(r.residual_median_polonen,i,color=co,s=18);e.scatter(r.residual_median_aieop,i,facecolor='white',edgecolor=co,s=22)
    e.set_yticks(range(len(ex)),[f'{short(r.predicted_subtype)} {int(r.n_polonen)}/{int(r.n_aieop)}' for r in ex.itertuples()],fontsize=8)
    e.set(xlabel='Cohort-specific median R',ylim=(-.7,6.7));e.axvline(0,color=LIGHT,lw=.5)
    f.legend(handles=[Line2D([],[],marker='o',ls='',color=INK,label='Discovery'),Line2D([],[],marker='o',ls='',markerfacecolor='white',markeredgecolor=INK,label='AIEOP')],loc='center',bbox_to_anchor=anchor(f,132,151),ncol=2,frameon=False)
    save(f,2)
if __name__=='__main__':main()

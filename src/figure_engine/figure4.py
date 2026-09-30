from plot_common import *
def main():
    f=fig(187);p=read(SD/'revision3_robustness/bcl11b_etp_matched_pairs.tsv');assert len(p)==18
    a=ax(f,box(f,23,10,57,31),'A','Fixed developmental matching');b=ax(f,box(f,116,10,57,31),'B','Paired residual contrast')
    for r in p.itertuples():
        for a0,vals in [(a,[r.BCL11B_dev,r.ETP_dev]),(b,[r.BCL11B_residual,r.ETP_residual])]:
            a0.plot([0,1],vals,color=LIGHT,lw=.65);a0.scatter([0,1],vals,color=[Z2,ETP],s=18,zorder=3)
    for z in [a,b]:z.set_xticks([0,1],['BCL11B','ETP-like']);z.set_xlim(-.2,1.2);grid(z)
    a.set_ylabel('Developmental\nproxy D');b.set_ylabel('Cohort-internal residual R')
    a.set_ylim(p[['BCL11B_dev','ETP_dev']].min().min()-.1,-1.6)
    a.text(.02,.98,'18 pairs\nMedian |ΔD|=.016',transform=a.transAxes,va='top')
    b.text(.02,.98,'Mean ΔR = −1.92\n95% CI [−2.98, −.67]',transform=b.transAxes,va='top')
    c=ax(f,box(f,23,66,57,33),'C','Internal held-out programs');c.set_xlim(-.5,1.5);c.set_ylim(0,55)
    overlap=read(D/'gate2_sensitivity/heldout_vs_original_overlap.tsv');camera=read(D/'gate2_sensitivity/heldout_matched_cameraPR.tsv')
    for i,side in enumerate(['ZEB1-side50','ZEB2-side50']):
        r=overlap[overlap.side==side].iloc[0];s=camera[camera.program==side].iloc[0]
        c.bar(i,50,color='white',edgecolor=color(side),lw=1);c.bar(i,r.n_overlap,color=color(side),width=.8)
        c.text(i,r.n_overlap/2,f'{int(r.n_overlap)}/50',ha='center',color='white',weight='bold')
        c.text(i,51,f'FDR {s.FDR:.2g}',ha='center')
    c.set_xticks([0,1],['ZEB1 side','ZEB2 side']);c.set_ylabel('Shared genes')
    d=ax(f,box(f,116,66,57,33),'D','Matched transcriptome')
    g=read(SD/'revision4_targeted/matched_BCL11B_ETP_gene_results.tsv');fr=read(D/'gate2_program/frozen_ZEB_side50.tsv')
    g=g.merge(fr[['symbol','side']],on='symbol',how='left',validate='one_to_one');yy=-np.log10(g['adj.P.Val'].clip(lower=1e-300))
    d.scatter(g.logFC,yy,color='#D9D9D9',s=2,rasterized=True)
    for side in ['ZEB1-side50','ZEB2-side50']:
        m=g.side==side;d.scatter(g.loc[m,'logFC'],yy[m],color=color(side),s=10,rasterized=True)
    d.axhline(-np.log10(.05),color=LIGHT,ls='--',lw=.5);d.axvline(0,color=LIGHT,lw=.5);d.set(xlim=(-7.7,7.4),ylim=(0,7.5),xlabel='BCL11B − ETP-like log₂ effect',ylabel='−log₁₀ BH FDR');grid(d)
    e=ax(f,box(f,70,122,103,54),'E','Lesion-defined cases')
    cases=read(SD/'F4/F4C_GSE162280_cases.tsv').sort_values('log2_ZEB1_over_ZEB2');assert len(cases)==12
    for i,r in enumerate(cases.itertuples()):
        e.plot([r.ZEB1_cpm,r.ZEB2_cpm],[i,i],color=LIGHT,lw=.8);e.scatter(r.ZEB1_cpm,i,s=22,color=Z1);e.scatter(r.ZEB2_cpm,i,s=22,color=Z2)
    e.set_yticks(range(12),[f'{r.case} | {r.phenotype_class}' for r in cases.itertuples()]);e.invert_yaxis();e.set_xscale('symlog',linthresh=20)
    e.set_xticks([0,10,100,300],['0','10','100','300']);e.set_xlabel('Expression (CPM; symlog)');grid(e,'x')
    f.legend(handles=[Line2D([],[],marker='o',ls='',color=Z1,label='ZEB1'),Line2D([],[],marker='o',ls='',color=Z2,label='ZEB2')],frameon=False,ncol=2,loc='center',bbox_to_anchor=anchor(f,131,115))
    key=f.legend(handles=[Patch(facecolor=GREY,label='Overlap'),Patch(facecolor='white',edgecolor=INK,label='Set size')],loc='center',bbox_to_anchor=anchor(f,51,112),ncol=2,frameon=False);f.add_artist(key)
    save(f,4)
if __name__=='__main__':main()

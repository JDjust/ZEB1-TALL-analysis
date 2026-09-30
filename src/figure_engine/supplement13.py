from plot_common import *
def main():
    f=fig(200);a=ax(f,[.13,.70,.32,.23],'A','Adult gene effects (moved duplicate)');g=read(D/'a8_adult/a8_genes.tsv').dropna(subset=['beta_R','beta_adult']);a.scatter(g.beta_R,g.beta_adult,c=g.side.map(color),s=18);zeros(a);a.set(xlabel='Discovery βR',ylabel='Adult βR');a.text(.02,.97,'n=78 model-tested genes\n86 recovered; 77/78 expected',transform=a.transAxes,va='top');a.margins(y=.3);grid(a)
    b=ax(f,[.63,.70,.32,.23],'B','BM/PB source adjustment');s=read(D/'a7_lineage/a7_0_frozen_genes_specimen_adjusted.tsv');b.scatter(s.beta_primary,s.beta_specimen,c=s.side.map(color),s=16);lim=max(s.beta_primary.abs().max(),s.beta_specimen.abs().max())*1.1;b.plot([-lim,lim],[-lim,lim],color=GREY,ls='--',lw=.7);b.set(xlabel='Primary βR',ylabel='BM/PB-adjusted βR');b.text(.02,.97,'100/100 same sign; ρ=.9995\nScore source test P=.13',transform=b.transAxes,va='top');b.margins(y=.3);grid(b)
    su=read(D/'a7c_decompose/a7c_summary.tsv')
    for j,(co,stat,title,xlab) in enumerate([('GSE243914','MWU_ETP_greater','Blast-enriched subset scores','ETP − non-ETP score'),('GSE248287','spearman_vs_B','Malignant subset scores','Spearman ρ versus patient B')]):
        a=ax(f,[.26+j*.47,.40,.21,.16],['C','D'][j],title);z=su[(su.cohort==co)&(su.stat==stat)].set_index('subset');vals=[z.loc[k,'estimate'] for k in ['myeloid_concordant','non_myeloid_concordant']];a.scatter(vals,[0,1],s=38,c=[Z2,INK]);a.set_yticks([0,1],['45 concordant','3 residual']);a.set(ylim=(1.6,-.6),xlabel=xlab);a.axvline(0,color=LIGHT,lw=.6);grid(a,'x')
    e=ax(f,[.30,.12,.55,.15],'E','Age-eligibility sensitivity (existing cohort)');pt=read(D/'a8_adult/a8_patient_scores.tsv')
    for i,(label,q) in enumerate([('All (n=79)',pt),('Age ≥18 y (n=72)',pt[pt.age>=18]),('Age ≥21 y (n=70)',pt[pt.age>=21])]):
        for metric,co,off in [('ZEB1_side50',Z1,-.1),('ZEB2_side50',Z2,.1)]:e.scatter(q.R.corr(q[metric],method='spearman'),i+off,s=28,color=co)
    e.set_yticks(range(3),['All (n=79)','Age ≥18 y (n=72)','Age ≥21 y (n=70)']);e.set(xlim=(-.85,.85),ylim=(2.5,-.5),xlabel='Spearman ρ of side score versus cohort residual R');e.axvline(0,color=LIGHT,lw=.6);grid(e,'x')
    sidekey(f,.01);save(f,'S13')
if __name__=='__main__':main()

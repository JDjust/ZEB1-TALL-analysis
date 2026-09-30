from plot_common import *
def main():
    f=fig(225);a=ax(f,[.27,.57,.24,.36],'A','All mapped spatial states')
    cma=read(D/'a4_yayon/a4b_tstate_cma_by_donor.tsv');order=cma.groupby('state').cma_weighted_median.median().sort_values().index;assert len(order)==24
    for i,name in enumerate(order):
        q=cma[cma.state==name];co=Z2 if 'DN' in name or 'ETP' in name else Z1 if 'DP' in name else GREY;a.scatter(q.cma_weighted_median,np.repeat(i,len(q)),s=10,color=co,alpha=.65);a.plot([q.cma_weighted_median.median()]*2,[i-.25,i+.25],color=INK,lw=1)
    a.set_yticks(range(24),[s.removeprefix('T_') for s in order]);a.set(ylim=(23.7,-.7),xlabel='CMA position');a.axvline(0,color=LIGHT,lw=.5);grid(a,'x')
    b=ax(f,[.76,.57,.19,.36],'B','Subtype medians\nand 95% CI');st=read(SD/'F3/F3D-F_subtype_statistics.tsv').sort_values('residual_within_median')
    for i,r in enumerate(st.itertuples()):b.plot([r.residual_within_ci_low,r.residual_within_ci_high],[i,i],color=color(r.subtype),lw=1);b.scatter(r.residual_within_median,i,s=18,color=color(r.subtype))
    b.set_yticks(range(17),st.subtype.map(short));b.set(ylim=(16.5,-.5),xlabel='Median R (95% CI)');b.axvline(0,color=LIGHT,lw=.5);grid(b,'x')
    c=ax(f,[.30,.13,.55,.32],'C','Complete subtype expression summary')
    mini=read(V/'F2A_subtype_descriptive_medians.tsv').set_index('subtype').loc[st.subtype];raw=mini[['D','B','R','Z1','Z2']].to_numpy();z=(raw-raw.mean(axis=0))/raw.std(axis=0,ddof=0)
    im=c.imshow(z,aspect='auto',cmap=CMAP,vmin=-2,vmax=2);c.set_xticks(range(5),['Median D','Median B','Median R','ZEB1 score','ZEB2 score']);c.set_yticks(range(17),[f'{short(name)} (n={int(n)})' for name,n in zip(mini.index,mini.n)]);c.tick_params(length=0);c.spines[:].set_visible(False)
    scale(f,im,[.45,.065,.24,.008],'Column display z',[-2,0,2]);f.text(.5,.02,'A: 24 predefined states × 6 spatial donors. B: original 3,000-bootstrap summaries. C: all 17 × 5 values.',ha='center')
    save(f,'S11')
if __name__=='__main__':main()

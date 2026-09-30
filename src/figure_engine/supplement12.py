from plot_common import *
def main():
    f=fig(205);a=ax(f,[.13,.70,.81,.22],'A','Frozen-program ranking and eligibility')
    passed=read(D/'gate2_program/effect_threshold_pass.tsv');assert len(passed)==3560
    for sign,co,lab in [(1,Z1,'ZEB1 side: 499 eligible effects'),(-1,Z2,'ZEB2 side: 3,061 eligible effects')]:
        q=passed[passed.beta_R*sign>0].sort_values('partial_R2',ascending=False);ranks=np.arange(1,len(q)+1);a.plot(ranks,q.partial_R2,color=co,alpha=.4,lw=.8);a.plot(ranks[:50],q.partial_R2.iloc[:50],color=co,lw=1.8,label=lab)
    a.axvline(50,color=GREY,ls=':',lw=.7);a.set(xscale='log',xlabel='Within-side ranking position',ylabel='Moderated-t-derived ranking score');a.legend(frameon=False,loc='upper right');grid(a)
    a.set_xticks([1,10,50,500,3000],['1','10','50','500','3000']);a.minorticks_off()
    f.text(.13,.613,'24,361 initially testable → exclude ZEB1/ZEB2/CD1A/CD34/LYL1 → 24,356 eligible.\nFixed BH FDR<.01 and |βR|≥.25: 3,560 pass; freeze top 50 per side.',linespacing=1.5)
    b=ax(f,[.30,.355,.58,.19],'B','Matched Hallmark context (original significant display)')
    hall=read(SD/'revision4_targeted/matched_BCL11B_ETP_hallmark_camera.tsv');q=hall[hall.FDR<.05].sort_values('FDR').head(6)
    labels={'INTERFERON_GAMMA_RESPONSE':'IFN-γ response','IL6_JAK_STAT3_SIGNALING':'IL6/JAK/STAT3','HEME_METABOLISM':'Heme metabolism','ALLOGRAFT_REJECTION':'Allograft rejection','INTERFERON_ALPHA_RESPONSE':'IFN-α response','OXIDATIVE_PHOSPHORYLATION':'Oxidative phosphorylation'}
    b.hlines(range(len(q)),0,-np.log10(q.FDR),color=LIGHT,lw=.8);b.scatter(-np.log10(q.FDR),range(len(q)),s=q.NGenes*.18,color=GREY)
    b.set_yticks(range(len(q)),[labels.get(v.replace('HALLMARK_',''),v) for v in q.pathway]);b.invert_yaxis();b.set_xlabel('−log10 CAMERA FDR');grid(b,'x')
    f.text(.30,.29,'Point area = gene count; all 50 terms remain in Table S5b.\nPair-blocked expression context; no pathway-causation inference.',linespacing=1.5)
    c=ax(f,[.13,.075,.81,.12],'C','Seven adult BCL11B individuals (sorted by existing R)');c.axis('off')
    ad=read(D/'a8_adult/a8_2_merged.tsv');z=ad[ad['T-ALL main-cluster high-confidence'].astype(str).str.contains('BCL11B',na=False)].sort_values('R');assert len(z)==7
    rows=[[str(i+1),str(int(r.age)),f'{r.R:.2f}',f'{r.B:.2f}',f'{r.ZEB1_side50:.2f}',f'{r.ZEB2_side50:.2f}','No'] for i,r in enumerate(z.itertuples())]
    table=c.table(cellText=rows,colLabels=['Order','Age (y)','R','B','ZEB1 score','ZEB2 score','ETP-like'],bbox=[0,0,1,1],cellLoc='center');table.auto_set_font_size(False);table.set_fontsize(8)
    for cell in table.get_celld().values():cell.set_edgecolor('#DDDDDD');cell.set_linewidth(.4)
    f.text(.5,.013,'Original measured individuals retained; 7/7 negative R. Small subtype, descriptive only.',ha='center')
    save(f,'S12')
if __name__=='__main__':main()

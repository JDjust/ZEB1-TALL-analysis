"""S10: measured perturbation, raw narrowPeak genomic tracks, assay boundary, clinical attenuation."""
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from matplotlib.patches import Rectangle

import editorial_rebuild_figures as b
import editorial_rebuild_supplement as s
import final_visual_style_lock as style
style.apply_style()
style.lock_legacy_module(b)

D,SD,OUT=b.D,b.SD,s.OUT
Z1,Z2,ETP,GREY,LIGHT,DARK=style.ZEB1,style.ZEB2,style.ETP,style.GREY,style.REF,style.INK
read,panel,grid=b.read,b.panel,b.clean_grid

ROOT=D/'a3_gse165209'
CHIP=ROOT/'gse165016'
# GRCh38 / Ensembl, copied from the frozen A3 locus-analysis manifest.
LOCI={'ZEB1':('chr10',31318495,31529814,'+',31318495),
      'ZEB2':('chr2',144364364,144521057,'-',144521057)}
TRACKS=[('DND41 BCL11B','GSM5024537_DND41_BCL11B_peaks.narrowPeak.gz',Z2),
        ('SJ005006 rep1','GSM5024540_SJTALL005006_BCL11B_rep1_peaks.narrowPeak.gz',Z2),
        ('SJ005006 rep2','GSM5024541_SJTALL005006_BCL11B_rep2_peaks.narrowPeak.gz',Z2),
        ('SJALL068666','GSM5265324_BCL11B_E15804_peaks.narrowPeak.gz',Z2),
        ('DND41 H3K27ac','GSM5024536_DND41_H3K27ac_peaks.narrowPeak.gz',Z1),
        ('SJ005006 H3K27ac','GSM5024539_SJTALL005006_H3K27ac_peaks.narrowPeak.gz',Z1)]

def perturbation(fig,spec):
    x=read(ROOT/'gse173432/a3_paired_donor_deltas.tsv')
    assert len(x)==3 and x.donor.is_unique
    inner=spec.subgridspec(3,2,wspace=.33,hspace=.75,
                           height_ratios=[.78,1,1])
    header=fig.add_subplot(inner[0,:]);header.axis('off')
    header.text(0,.96,'A',fontsize=10,weight='bold',color=DARK,va='top')
    header.text(.10,.96,'Acute BCL11B overexpression',fontsize=7.8,
                weight='semibold',color=DARK,va='top')
    header.text(.10,.05,'Three paired CD34+ donors',fontsize=6.7,
                color=GREY,va='bottom')
    items=[('log2 ratio','log2_ratio_EV','log2_ratio_OE'),
           ('ZEB1 TPM','ZEB1_EV','ZEB1_OE'),
           ('ZEB2 TPM','ZEB2_EV','ZEB2_OE'),
           ('Developmental D','dev_EV','dev_OE')]
    for i,(name,ev,oe) in enumerate(items):
        ax=fig.add_subplot(inner[1+i//2,i%2])
        for _,r in x.iterrows():
            ax.plot([0,1],[r[ev],r[oe]],color=Z2 if i==2 else Z1,lw=.8,alpha=.8)
            ax.scatter([0,1],[r[ev],r[oe]],s=9,color=[GREY,Z2 if i==2 else Z1],zorder=2)
        ax.set_xticks([0,1],['EV','OE'] if i>=2 else ['', ''],fontsize=6.5)
        ax.set_xlim(-.2,1.2);ax.tick_params(labelsize=6.5)
        ax.set_title(name,fontsize=7.2,loc='left',color=DARK,weight='bold')
        ax.spines[['top','right']].set_visible(False)

def _raw_peaks(path,chrom,start,end):
    assert path.exists(),path
    x=pd.read_csv(path,sep='\t',header=None,usecols=[0,1,2,6],
                  names=['chrom','start','end','signal'],compression='gzip')
    return x[(x.chrom.eq(chrom))&(x.end.ge(start-22000))&
             (x.start.le(end+22000))].copy()

def locus(ax,gene):
    chrom,start,end,strand,tss=LOCI[gene]
    length=(end-start)/1000
    ax.axvspan(-2,2,color=style.NA,alpha=.85,zorder=0)
    ax.axvspan(0,length,color=style.MIDPOINT,alpha=.7,zorder=0)
    all_count=0
    for i,(label,file,color) in enumerate(TRACKS):
        q=_raw_peaks(CHIP/file,chrom,start,end)
        all_count+=len(q)
        if strand=='+':
            left=(q.start-tss)/1000;right=(q.end-tss)/1000
        else:
            left=(tss-q.end)/1000;right=(tss-q.start)/1000
        maxsig=max(1,float(q.signal.quantile(.98)) if len(q) else 1)
        for l,r,v in zip(left,right,q.signal):
            ax.add_patch(Rectangle((l,i+.10),max(.22,r-l),
                                   .7*min(1,float(v)/maxsig),
                                   facecolor=color,edgecolor='none',alpha=.85))
    assert all_count>0
    ax.set_xlim(-20,length+20);ax.set_ylim(6.05,-.10)
    ax.set_yticks(np.arange(6)+.42,[q[0] for q in TRACKS],fontsize=6.5)
    ax.set_xlabel(f'{gene} distance from TSS (kb; transcriptional direction →)')
    ax.axvline(0,color=DARK,lw=.7)
    ax.spines[['top','right','left']].set_visible(False)
    ax.tick_params(axis='y',length=0)
    ax.text(.99,.94,f'{gene} · body {length:.0f} kb',
            transform=ax.transAxes,fontsize=6.7,weight='bold',color=DARK,
            va='top',ha='right')
    return all_count

def assay_matrix(ax):
    x=read(ROOT/'locus/a3_chip_zeb_loci.tsv')
    assert len(x)==12
    # Counts are explicit, assay-specific measured promoter/body peak calls.
    m=[]
    for assay in ['BCL11B_ChIP','H3K27ac']:
        row=[]
        for gene in ['ZEB1','ZEB2']:
            q=x[(x.assay.eq(assay))&(x.gene.eq(gene))]
            row.extend([int(q.promoter_hit.sum()),int(q.body_hit.sum())])
        m.append(row)
    labels=['BCL11B ChIP','H3K27ac']
    for i,row in enumerate(m):
        for j,count in enumerate(row):
            den=4 if i==0 else 2
            ax.add_patch(Rectangle((j-.45,i-.4),.9,.8,
                                   facecolor=style.SEQUENTIAL_CMAP(.35 if count<den else .8),
                                   edgecolor='white',lw=.5))
            ax.text(j,i,f'{count}/{den}',ha='center',va='center',fontsize=6.7,color=DARK)
    # Other assays did not offer the same independently replicated promoter/body
    # readout; NA is explicitly distinguished from a negative measured result.
    for i,label in enumerate(['HiChIP','Pooled loops','scATAC'],2):
        labels.append(label)
        for j in range(4):
            ax.add_patch(Rectangle((j-.45,i-.4),.9,.8,facecolor=style.NA,
                                   edgecolor='white',lw=.5))
            ax.text(j,i,'n/a',ha='center',va='center',fontsize=6.5,color=GREY)
    ax.set_xlim(-.5,3.5);ax.set_ylim(4.6,-.65)
    ax.set_xticks(range(4),['Z1 prom.','Z1 body','Z2 prom.','Z2 body'],
                  rotation=25,ha='right',fontsize=6.5)
    ax.set_yticks(range(5),labels,fontsize=6.5);ax.tick_params(length=0)
    ax.spines[:].set_visible(False)
    panel(ax,'C','Assay-specific locus audit',
          'ChIP peak hits; other assay cells are not comparable')

def clinical(ax):
    x=read(SD/'S6/S6A_all_endpoint_models.tsv')
    names=['Induction failure','M2/M3 morphology','MRD >=0.1%',
           'MRD >=0.01%','Event-free survival','Overall survival']
    assert set(names)==set(x.endpoint)
    for i,name in enumerate(names):
        q=x[x.endpoint.eq(name)].set_index('adjustment')
        for lab,col,off in [('Unadjusted',Z2,-.13),('Subtype-adjusted',GREY,.13)]:
            r=q.loc[lab]; y=i+off
            ax.plot([r.ci_low,r.ci_high],[y,y],color=col,lw=1.05)
            if i>=4:
                ax.scatter(r.estimate,y,s=17,facecolor='white',edgecolor=col,
                           linewidth=.85,zorder=3)
            else:ax.scatter(r.estimate,y,s=17,color=col,zorder=3)
    ax.set_yticks(range(6),names,fontsize=6.5);ax.set_ylim(5.55,-.55)
    ax.axvline(1,color=LIGHT,lw=.65);ax.set_xlim(.5,1.23)
    ax.set_xlabel('OR / HR per unit balance (95% CI)')
    panel(ax,'D','Clinical attenuation',
          'Open circles: EFS/OS Cox sparse-subtype warnings; exploratory only')
    grid(ax,'x')

def sc10():
    f=b.mmfig(180,200);g=GridSpec(3,2,figure=f,height_ratios=[1.22,1.62,.98])
    perturbation(f,g[0,0]);c=f.add_subplot(g[0,1]);assay_matrix(c)
    inner=g[1,:].subgridspec(2,1,hspace=.62)
    z1=f.add_subplot(inner[0,0]);z2=f.add_subplot(inner[1,0])
    n1=locus(z1,'ZEB1');n2=locus(z2,'ZEB2')
    panel(z1,'B','BCL11B and H3K27ac locus tracks',
          f'GRCh38 narrowPeak intervals; {n1+n2} peaks, signal scaled within track')
    d=f.add_subplot(g[2,:]);clinical(d)
    f.subplots_adjust(left=.22,right=.965,top=.945,bottom=.065,wspace=.55,hspace=.48)
    style.polish_text(f)
    f.savefig(OUT/'FigureS10.pdf');f.savefig(OUT/'FigureS10.png',dpi=300)
    plt.close(f)

if __name__=='__main__':sc10()

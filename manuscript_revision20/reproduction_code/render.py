"""Reader-led layout revision; frozen-result rendering without new inference."""
from pathlib import Path
import sys,json,types
import numpy as np,pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.colors import ListedColormap
from scipy.stats import rankdata
B=Path(__file__).resolve().parent.parent
old=types.ModuleType('layout1_helpers');old.__file__=str(Path(__file__).with_name('render_layout1.py'))
exec(compile(Path(old.__file__).read_text(encoding='utf-8').split('FUNCTIONS=')[0],old.__file__,'exec'),old.__dict__)
t=old.t;data=t.data;js=t.js;C=t.C;scatter=t.scatter;effect=t.effect
F=B/'figures';F.mkdir(exist_ok=True);t.OUT=F
W,H=240,160
def analysisdata(name, directory='results'):
    return pd.read_csv(B/'analysis'/directory/name, sep='\t')

plt.rcParams.update({'font.family':'Arial','font.size':8.5,'axes.labelsize':8.5,'xtick.labelsize':8.5,'ytick.labelsize':8.5,'legend.fontsize':8,'axes.linewidth':.65,'lines.linewidth':1.0,'text.color':'black','axes.labelcolor':'black','axes.titlecolor':'black','savefig.facecolor':'white'})
PANELS=[];CURRENT=''
def start(name):
    global CURRENT
    CURRENT=name;return plt.figure(figsize=(190.5/25.4,160/25.4))
def ax(f,l,r,kind,sources,role='COMPACT_SUPPORT',parent=''):
    x,y,w,h=r;a=f.add_axes([x/W,y/H,w/W,h/H]);a.spines[['top','right']].set_visible(False)
    f.text(max(2,x-9)/W,min(H-3,y+h+3)/H,l,fontsize=12,fontweight='bold',va='bottom')
    PANELS.append(dict(figure=CURRENT,panel=l,status='SUPP' if CURRENT.startswith('S') else 'MAIN',plot_type=kind,relative_visual_weight='HERO' if role=='HERO' else 'L' if role=='MAJOR_SUPPORT' else 'M',role=role,parent_main=parent,source_aliases=';'.join(sources),font_family='Arial',min_font_pt=8,semantic_palette='ZEB1-development-v2',rect=json.dumps([x/W,y/H,w/W,h/H])))
    a._panel_receipt=PANELS[-1]
    return a
def tag(a,s):a.text(0,1+2.5/(a.get_position().height*H),s,transform=a.transAxes,fontsize=8,va='bottom',color='black')
def save(f,n):
    p=F/n;p.mkdir(exist_ok=True)
    for text in f.findobj(matplotlib.text.Text):
        if text.get_text():text.set_fontsize(max(8,min(12,text.get_fontsize())));text.set_fontfamily('Arial')
    # Keep actual 8-point labels inside the journal-size page; transform every
    # axes position consistently, retaining all observations and coordinates.
    for _ in range(6):
        f.canvas.draw();renderer=f.canvas.get_renderer()
        boxes=[v.get_window_extent(renderer) for v in f.findobj(matplotlib.text.Text) if v.get_visible() and v.get_text()]
        left=min(b.x0 for b in boxes)/f.bbox.width;right=max(b.x1 for b in boxes)/f.bbox.width
        if left>=.007 and right<=.993:break
        lo=min(left,0);hi=max(right,1);scale=.982/(hi-lo);shift=.009-lo*scale
        for a in f.axes:
            r=a.get_position();a.set_position([shift+scale*r.x0,r.y0,scale*r.width,r.height])
        for v in f.texts:
            x,y=v.get_position();v.set_position((shift+scale*x,y))
        for legend in f.legends:
            r=legend.get_bbox_to_anchor().transformed(f.transFigure.inverted())
            legend.set_bbox_to_anchor((shift+scale*r.x0,r.y0,scale*r.width,r.height),transform=f.transFigure)
    frames=[]
    for a in f.axes:
        r=list(a.get_position().bounds);frames.append(dict(rect=r,has_panel=hasattr(a,'_panel_receipt')))
        if hasattr(a,'_panel_receipt'):
            a._panel_receipt['actual_plot_rect']=json.dumps(r)
    (p/'FINAL_LAYOUT.json').write_text(json.dumps(dict(width_mm=190.5,height_mm=160,axes=frames),indent=2),encoding='utf-8')
    for ext in ['pdf','svg','png']:f.savefig(p/(n+'.'+ext),dpi=300)
    plt.close(f);print(n+' RENDERED',flush=True)
def heat(f,a,m,label,vmin=-2,vmax=2,cmap='RdBu_r',cbar=True,rotation=60):
    cm=plt.get_cmap(cmap).copy();cm.set_bad('#DADADA');im=a.imshow(m.to_numpy(dtype=float),aspect='auto',interpolation='nearest',cmap=cm,vmin=vmin,vmax=vmax,rasterized=True)
    a.set_yticks(range(len(m)),m.index);a.set_xticks(range(len(m.columns)),m.columns,rotation=rotation,ha='right' if rotation else 'center')
    if cbar:f.colorbar(im,ax=a,fraction=.035,pad=.04).set_label(label)
    return im
def sourcebox(a,d,g,v,order,colors,points=True):
    vals=[d.loc[d[g]==k,v].dropna().values for k in order];b=a.boxplot(vals,positions=range(len(order)),widths=.56,patch_artist=True,showfliers=False,medianprops={'color':'black','lw':1},boxprops={'lw':.65},whiskerprops={'lw':.65},capprops={'lw':.65});rng=np.random.default_rng(20261010)
    for i,(p,z) in enumerate(zip(b['boxes'],vals)):
        p.set_facecolor(colors[i]);p.set_alpha(.35)
        if points:a.scatter(i+rng.uniform(-.19,.19,len(z)),z,s=4 if len(z)>100 else 10,color=colors[i],alpha=.4,edgecolors='none')
    a.set_xticks(range(len(order)),order)

def raincloud(a,d,g,v,order,colors):
    """Actual patient values: half violin, compact IQR box and all points.

    Density is omitted for groups below ten patients. No inferential test.
    """
    vals=[d.loc[d[g]==k,v].dropna().to_numpy() for k in order]
    rng=np.random.default_rng(20261010)
    for i,(z,c) in enumerate(zip(vals,colors)):
        if len(z)>=10 and np.ptp(z)>0:
            violin=a.violinplot([z],positions=[i-.04],widths=.72,showextrema=False)
            body=violin['bodies'][0];verts=body.get_paths()[0].vertices
            verts[:,0]=np.minimum(verts[:,0],i-.04)
            body.set_facecolor(c);body.set_edgecolor(c);body.set_alpha(.28);body.set_linewidth(.5)
        box=a.boxplot([z],positions=[i+.035],widths=.13,patch_artist=True,showfliers=False,
            medianprops={'color':'black','lw':.8},boxprops={'lw':.6},
            whiskerprops={'lw':.6},capprops={'lw':.6})
        box['boxes'][0].set_facecolor(c);box['boxes'][0].set_alpha(.65)
        a.scatter(i+rng.uniform(.15,.34,len(z)),z,s=3 if len(z)>100 else 9,
            color=c,alpha=.22 if len(z)>100 else .65,edgecolors='none',rasterized=len(z)>500)
    a.set_xticks(range(len(order)),order);a.set_xlim(-.48,len(order)-.48)
def counted_groups(a,d,g,v,order,labels):
    """Source-observed group n on compact horizontal decoding labels."""
    a.set_xticks(range(len(order)),[f'{label}\nn={d.loc[d[g].eq(k),v].notna().sum():,}' for k,label in zip(order,labels)])
    a.tick_params(axis='x',labelsize=7.5)
def pathpoints(a,genes):
    for i,g in enumerate(genes):
        for j,k in enumerate(['hallmark_raw','hallmark_residual']):
            r=t.normalized_pathway(k).set_index('pathway').loc['HALLMARK_'+g];c=C['early' if j==0 else 'ZEB1'];v=-r.NES if j==0 else r.NES
            a.scatter(v,i+(j-.5)*.2,s=24,marker='o' if j==0 else 's',facecolors=c if r.q<.05 else 'white',edgecolors=c,lw=.8)
    a.set_yticks(range(len(genes)),[t.short_pathway('HALLMARK_'+g) for g in genes]);a.invert_yaxis();a.axvline(0,color='black',lw=.6);a.set_xlabel('NES: + = ZEB1-high')
    lo,hi=a.get_xlim();pad=(hi-lo)*.16;a.set_xlim(lo-pad,hi+pad)
def pathwaykey(f,pos):f.legend([Line2D([0],[0],marker=m,ls='',color=c,mfc=fc,ms=4) for m,c,fc in [('o',C['early'],C['early']),('s',C['ZEB1'],C['ZEB1']),('s',C['ZEB1'],'white')]],['Marginal (display −NES)','Subtype residual','Open: source q ≥ .05'],frameon=False,ncol=3,loc='center',bbox_to_anchor=pos,columnspacing=1.2)
def withinphenotype(d,field):
    r=rankdata(d[field].to_numpy(),method='average');g=d.source_disease_state.to_numpy();return r-np.array([r[g==k].mean() for k in g])
def residualscatter(a,axis):
    d=data('aiims_native');x=withinphenotype(d,axis+'_raw');y=withinphenotype(d,'ZEB1_log2_FPKM_plus1');assert len(d)==35
    for g,c in [('Immature_TALL',C['early']),('Cortical_TALL',C['cortical']),('Mature_TALL',C['SP'])]:
        keep=d.source_disease_state.eq(g);a.scatter(x[keep],y[keep],s=22,color=c,edgecolors='none',alpha=.8,label=g.replace('_TALL',''))
    a.axhline(0,color='#AAAAAA',lw=.5);a.axvline(0,color='#AAAAAA',lw=.5);a.set_xlabel(('SP / cortical' if axis=='SP_vs_cortical' else 'Cortical / early')+' rank residual');a.set_ylabel('ZEB1 rank residual');r=data('aiims_conditional').set_index('axis').loc[axis+'_raw'];tag(a,f'AIIMS n=35 · ρ={r.conditional_rank_rho:.2f}; q={r.BH2_q:.3g}')
def cohort(a,k,axis):
    d,y=t.cohort(k);scatter(a,d,axis+'_raw',y,C['SP' if axis=='SP_vs_cortical' else 'cortical'],xlabel='SP / cortical score' if axis=='SP_vs_cortical' else 'Cortical / early score',ylabel={'nopho':'ZEB1 log2(CPM+1)','aiims':'ZEB1 log2(FPKM+1)','india':'ZEB1 log2(counts+1)'}[k]);r=next(v for v in js(k+'_results')['results'] if v['axis']==axis);tag(a,f'{k.upper()} n={len(d)}\nρ={r["spearman_rho"]:.2f}; q={r["BH2_q"]:.2g}')
def donor(a,axis):
    d=data('normal_state_scores');p=d.pseudobulk_id.str.split('|',expand=True);d['donor']=p[0];d['stage']=p[1];order=['Early_DN','DP','CD4SP','CD8SP']
    for i,(k,s) in enumerate(d.groupby('donor')):a.plot(range(4),s.set_index('stage')[axis].reindex(order),'o-',color=plt.get_cmap('tab10')(i),ms=3,label=k.replace('TSC-',''))
    a.set_xticks(range(4),['Early DN','DP','CD4 SP','CD8 SP'],rotation=45,ha='right');a.set_ylabel('Fixed rank score');tag(a,'Cortical / early' if axis=='cortical_vs_early' else 'SP / cortical')
def pair(a,gene):
    p=data('sc_state').pivot(index='patient',columns='mapped_state',values=gene)[['Early_DN','DP']].dropna()
    for k,v in p.iterrows():a.plot([0,1],v.values,color='#999999',lw=.8);a.scatter([0,1],v.values,c=[C['early'],C['cortical']],s=20)
    a.set_xticks([0,1],['Early DN','DP']);a.set_ylabel(gene+' log2(CPM+1)');a.set_xlim(-.2,1.2);return p

def Fig1():
    f=start('Fig1');d=data('bulk').merge(data('maturity'),on='sample_id',validate='one_to_one');order=d.groupby('subtype').ZEB1.median().sort_values().index.tolist()
    a=ax(f,'a',[16,111,217,37],'all-patient subtype distributions',['bulk'],'HERO');sourcebox(a,d,'subtype','ZEB1',order,[C['ZEB1']]*len(order));a.set_ylabel('ZEB1 log2(TPM+1)');a.set_xticklabels(order,rotation=55,ha='right');a.tick_params(axis='x',labelsize=8);a.set_ylim(0,8.5)
    for i,g in enumerate(order):a.text(i,8.05,f'n={d.subtype.eq(g).sum()}',ha='center',va='center',fontsize=7.5)
    tag(a,'1,309 diagnostic patients · orange: patients · box: median and IQR')
    a=ax(f,'b',[30,17,48,55],'subtype developmental context',['bulk','maturity'],'MAJOR_SUPPORT');m=d.groupby('subtype')[['cortical_vs_early','SP_vs_cortical']].mean().reindex(order);m.columns=['Cortical / early','SP / cortical'];heat(f,a,m,'Mean score (discovery z)',-1.5,1.5,rotation=40)
    a=ax(f,'c',[122,62,37,23],'native ETP patient raincloud',['bulk']);raincloud(a,d,'etp3','ZEB1',['ETP','Near-ETP','Non-ETP'],[C['early'],C['cortical'],C['SP']]);counted_groups(a,d,'etp3','ZEB1',['ETP','Near-ETP','Non-ETP'],['ETP','Near-ETP','Non-ETP']);a.set_xticklabels(['ETP\nn=110','Near\nETP\nn=168','Non\nETP\nn=862']);a.set_ylabel('ZEB1 log2(TPM+1)')
    r=analysisdata('ALL1309_RANKS_AND_SUBTYPE_RESIDUALS.tsv','source_data')
    for letter,axis,rect in [('d','cortical_vs_early',[195,62,38,23]),('e','SP_vs_cortical',[122,14,37,20])]:
        a=ax(f,letter,rect,'within-subtype patient rank residual scatter',['bulk_rank_residuals'],'MAJOR_SUPPORT');scatter(a,r,axis+'_subtype_residual','ZEB1_subtype_residual',C['cortical' if letter=='d' else 'SP'],xlabel=('Cortical / early' if letter=='d' else 'SP / cortical')+'\nrank residual',ylabel='ZEB1 rank residual');tag(a,'Within subtype\nn=1,309')
    a=ax(f,'f',[201,14,31,23],'marginal and subtype-adjusted rank effects',['bulk_rank_effects']);r=analysisdata('BULK_MARGINAL_ADJUSTED_EFFECTS.tsv');labels=[];values=[];low=[];high=[];colors=[]
    for axis,label in [('cortical_vs_early','C/E'),('SP_vs_cortical','SP/C')]:
        for model,word in [('marginal','total'),('subtype_adjusted','within')]:
            v=r[(r.axis==axis)&(r.model==model)].iloc[0];labels.append(label+' '+word);values.append(v.effect);low.append(v.ci_low);high.append(v.ci_high);colors.append(C['cortical' if axis=='cortical_vs_early' else 'SP'])
    effect(a,labels,values,low,high,'Rank slope β\n(95% CI)',colors);a.set_xlim(-.4,.5);a.set_xticks([-.3,0,.3]);tag(a,'Total / within');save(f,'Fig1')

def Fig2():
    f=start('Fig2');d=data('normal_umap');labels=data('normal_labels').set_index('author_cell_type').display_label.to_dict();groups=list(d.author_cell_type.unique());palette=plt.get_cmap('tab20b')(np.arange(len(groups))%20);a=ax(f,'a',[13,90,90,59],'normal atlas embedding',['normal_umap','normal_labels'],'HERO');texts=[];handles=[]
    for i,g in enumerate(groups):
        s=d[d.author_cell_type==g];a.scatter(s.UMAP1,s.UMAP2,s=.12,color=palette[i],linewidths=0,rasterized=True);texts.append(a.text(s.UMAP1.median(),s.UMAP2.median(),str(i+1),fontsize=8,color='black',ha='center',bbox={'facecolor':'white','alpha':.7,'edgecolor':'none','pad':.05}));handles.append(Line2D([0],[0],marker='o',ls='',ms=3,color=palette[i],label=str(i+1)+' '+labels.get(g,g).replace('Treg-','Tr-').replace('Treg_','Tr_').replace('SP_fetal_ambiguous','SP fetal amb.')))
    a.set_xlabel('UMAP 1');a.set_ylabel('UMAP 2');a.set_xticks([]);a.set_yticks([]);tag(a,'Normal thymus atlas · all 391,462 coordinate rows');np.random.seed(20261010);old.adjust_text(texts,ax=a,force_text=(.5,.5),force_static=(.2,.2),force_pull=(.015,.015),expand=(1.25,1.25),iter_lim=400,min_arrow_len=2,arrowprops={'arrowstyle':'-','color':'#777777','lw':.4})
    a=ax(f,'b',[138,90,89,59],'normal ZEB1 FeaturePlot',['normal_umap'],'HERO');t.umap(a,d,True);tag(a,'ZEB1 counts · localization only');f.legend(handles=handles,loc='center',bbox_to_anchor=(.50,.45),ncol=6,frameon=False,columnspacing=.75,handletextpad=.25,labelspacing=.45,fontsize=8.2)
    a=ax(f,'c',[14,18,46,31],'normal FACS-stage expression',['normal195']);d=data('normal195').set_index('facs_stage').reindex(['DN1','DN2','DN3','ISP','DP_CD3neg','DP_CD3pos','CD4SP','CD8SP']);a.plot(range(len(d)),d.g_ZEB1,'o-',color=C['ZEB1'],ms=3);a.set_xticks(range(len(d)),['DN1','DN2','DN3','ISP','DP CD3−','DP CD3+','CD4 SP','CD8 SP'],rotation=55,ha='right');a.set_ylabel('Mean normalized ZEB1');tag(a,'GSE195812 · FACS stages')
    a=ax(f,'d',[80,18,32,31],'normal inferred-stage expression',['normal206']);d=data('normal206').set_index('inferred_stage');a.plot(range(len(d)),d.g_ZEB1,'o-',color=C['ZEB1'],ms=3);a.set_xticks(range(len(d)),['DN early','DP imm.','DP mat.'],rotation=50,ha='right');a.set_ylabel('Mean normalized ZEB1');tag(a,'GSE206710')
    a=ax(f,'e',[140,18,32,31],'normal donor trajectories',['normal_state_scores','normal_heldout']);donor(a,'cortical_vs_early');a=ax(f,'f',[200,18,32,31],'normal donor trajectories',['normal_state_scores','normal_heldout']);donor(a,'SP_vs_cortical');f.legend(*a.get_legend_handles_labels(),frameon=False,ncol=5,loc='lower center',bbox_to_anchor=(.76,.005),columnspacing=.9,handlelength=1);save(f,'Fig2')
def Fig3():
    f=start('Fig3');d=data('lim_scores');r=data('lim_results');parts=[]
    for k in ['discovery','validation']:
        s=d[d.cohort_analysis==k].sort_values('ZEB1_log2_CPM_plus1');m=old.zcols(s[['ZEB1_log2_CPM_plus1','cortical_vs_early','SP_vs_cortical']]);m.index=[k[0].upper()+':'+p for p in s.patient];parts.append(m)
    m=pd.concat(parts).T;m.index=['ZEB1','Cortical /\nearly','SP /\ncortical'];a=ax(f,'a',[24,127,206,21],'malignant patient score landscape',['lim_scores'],'HERO');heat(f,a,m,'Within-set display z',rotation=90);a.axvline(20.5,color='black',lw=1);a.tick_params(axis='x',labelsize=7.5);tag(a,'Lim diagnosis BM · 21 discovery and 28 validation patients')
    a=ax(f,'b',[17,58,48,37],'same-patient malignant state contrast',['sc_state','sc_receipt'],'MAJOR_SUPPORT');p=pair(a,'ZEB1');a.set_xlim(-.2,1.45)
    last=-np.inf
    for k,v in p.sort_values('DP').iterrows():
        ly=max(v.DP,last+.105);last=ly;a.annotate(k,(1,v.DP),xytext=(1.10,ly),fontsize=7.5,arrowprops={'arrowstyle':'-','lw':.4,'color':'#777777'})
    a.set_ylim(p.to_numpy().min()-.08,max(p.to_numpy().max(),last)+.10)
    tag(a,'GSE248287 · five patient pairs\nDP − early DN: exact P=.0625')
    for l,k,axis,rect in [('c','discovery','cortical_vs_early',[112,66,40,29]),('d','validation','cortical_vs_early',[193,66,40,29]),('e','discovery','SP_vs_cortical',[112,17,40,29]),('f','validation','SP_vs_cortical',[193,17,40,29])]:
        a=ax(f,l,rect,'malignant patient scatter',['lim_scores','lim_results']);s=d[d.cohort_analysis==k];v=r[(r.cohort==k)&(r.axis==axis)&(r.model=='primary_spearman')].iloc[0];scatter(a,s,axis,'ZEB1_log2_CPM_plus1',C['cortical' if axis=='cortical_vs_early' else 'SP'],xlabel='Cortical / early' if axis=='cortical_vs_early' else 'SP / cortical',ylabel='ZEB1 log2(CPM+1)');tag(a,f'{k.title()} n={len(s)}\nρ={v.effect:.2f}; BH4 q={v.BH4_primary:.3g}')
    a=ax(f,'g',[32,18,36,27],'cycle-depth adjusted rank slopes',['lim_results']);s=r[r.model=='cycle_depth_rank_HC3'].sort_values(['cohort','axis']);labs=[v.cohort[0].upper()+': '+('C/E' if v.axis=='cortical_vs_early' else 'SP/C') for _,v in s.iterrows()];effect(a,labs,s.effect.tolist(),s.ci_low.tolist(),s.ci_high.tolist(),'HC3 rank slope');tag(a,'D: discovery · V: validation');save(f,'Fig3')
def Fig4():
    f=start('Fig4');d=data('mutations');genes=['CDKN2A','NOTCH1','TLX1','TLX3','DNM2','PHF6','IL7R','BCL11B','FLT3','NRAS','JAK3','MED12'];s=d[d.subtype=='TLX1'].sort_values(['ZEB1','sample_id']);a=ax(f,'a',[22,91,51,56],'focused genetic Oncoprint',['mutations'],'MAJOR_SUPPORT');t.onco(a,f,s,genes);tag(a,'TLX1 · all 71 patients');a.tick_params(axis='y',labelsize=8)
    for l in list(f.legends):l.remove()
    if a.get_legend() is not None:a.get_legend().remove()
    f.legend([plt.Rectangle((0,0),1,1,color=C['early']),plt.Rectangle((0,0),1,1,color='#F2F2F2')],['Altered','Absent'],frameon=False,ncol=2,loc='center',bbox_to_anchor=(.17,.977),fontsize=8,columnspacing=1)
    a=ax(f,'b',[99,91,39,56],'representative genetic subtype prevalence',['mutations'],'MAJOR_SUPPORT');backgrounds=['BCL11B','ETP-like','TAL1 DP-like','TAL1 αβ-like','TLX1','TLX3'];fields=['CDKN2A','NOTCH1','TLX1','TLX3','PHF6','MED12'];m=d.groupby('subtype')[fields].mean().reindex(backgrounds).T;heat(f,a,m,'Altered fraction',0,1,'YlOrBr',cbar=False);a.set_xticklabels(m.columns,rotation=60,ha='right');a.tick_params(axis='x',labelsize=8);a.tick_params(axis='y',labelsize=8);tag(a,'Altered fraction')
    f.legend([plt.Rectangle((0,0),1,1,color=plt.get_cmap('YlOrBr')(0.)),plt.Rectangle((0,0),1,1,color=plt.get_cmap('YlOrBr')(1.))],['0','1'],frameon=False,ncol=2,loc='center',bbox_to_anchor=(.49,.977),fontsize=8,columnspacing=1)
    a=ax(f,'c',[183,18,49,129],'all-subtype developmental coupling',['within_subtype_associations'],'HERO');r=analysisdata('ALL34_WITHIN_SUBTYPE_ASSOCIATIONS.tsv');order=sorted(r.subtype.unique());ys=np.arange(len(order))
    for axis,offset,color in [('cortical_vs_early',-.14,C['cortical']),('SP_vs_cortical',.14,C['SP'])]:
        v=r[r.axis==axis].set_index('subtype').loc[order]
        for y,row in zip(ys+offset,v.itertuples()):
            a.plot([row.ci_low,row.ci_high],[y,y],color=color,lw=.85,alpha=.85);a.plot(row.spearman_rho,y,'o',color=color,ms=3)
    ns=r.drop_duplicates('subtype').set_index('subtype').n_patients.to_dict();a.set_yticks(ys,[k+' ('+str(ns[k])+')' for k in order]);a.tick_params(axis='y',labelsize=7.7);a.invert_yaxis();a.axvline(0,color='black',lw=.6,ls='--');a.set_xlim(-1.05,1.05);a.set_xticks([-1,0,1]);a.set_xlabel('Within-group Spearman ρ\n95% patient-bootstrap CI');tag(a,'All 17 source groups · n=5–296')
    for l,x,g in [('d',20,'NOTCH1'),('e',71,'CDKN2A'),('f',122,'MED12')]:
        a=ax(f,l,[x,15,25,25],'lesion-linked patient raincloud',['mutations','mutation_assoc']);raincloud(a,d,g,'ZEB1',[0,1],[C['early'],C['ZEB1']]);counted_groups(a,d,g,'ZEB1',[0,1],['Absent','Present']);a.set_ylabel('ZEB1 log2(TPM+1)');v=data('mutation_assoc').set_index('gene').loc[g];tag(a,f'{g}\n37-field q={v.mw_q:.2g}')
    f.legend([Line2D([0],[0],marker='o',ls='',color=C['cortical']),Line2D([0],[0],marker='o',ls='',color=C['SP'])],['Cortical / early','SP / cortical'],frameon=False,ncol=2,loc='center',bbox_to_anchor=(.77,.977),fontsize=8,columnspacing=1.2);save(f,'Fig4')

def Fig5():
    f=start('Fig5');r=analysisdata('CROSS10_COHORT_EFFECTS_CI.tsv');a=ax(f,'a',[34,69,46,77],'complete external cohort effect intervals',['cross_cohort_effects'],'MAJOR_SUPPORT');cohorts=['NOPHO108','AIIMS35','India25'];ys=np.arange(3)
    for axis,offset,color in [('cortical_vs_early',-.14,C['cortical']),('SP_vs_cortical',.14,C['SP'])]:
        v=r[(r.axis==axis)&r.cohort.isin(cohorts)].set_index('cohort').loc[cohorts]
        for y,row in zip(ys+offset,v.itertuples()):
            a.plot([row.ci_low,row.ci_high],[y,y],color=color,lw=1.2);a.plot(row.effect,y,'o',color=color,ms=4)
    a.set_yticks(ys,['NOPHO (108)','AIIMS (35)','India (25)']);a.invert_yaxis();a.axvline(0,color='black',lw=.65,ls='--');a.set_xlim(-1,1);a.set_xticks([-1,-.5,0,.5,1]);a.set_xlabel('Spearman ρ · 95% CI');tag(a,'All six developmental associations')
    f.legend([Line2D([0],[0],marker='o',ls='',color=C['cortical']),Line2D([0],[0],marker='o',ls='',color=C['SP'])],['Cortical / early','SP / cortical'],frameon=False,fontsize=8,loc='center',bbox_to_anchor=(.29,.975),ncol=2,columnspacing=.8)
    a=ax(f,'b',[119,113,41,31],'NOPHO patient scatter',['nopho','nopho_results']);cohort(a,'nopho','cortical_vs_early');a=ax(f,'c',[193,113,40,31],'NOPHO patient scatter',['nopho','nopho_results']);cohort(a,'nopho','SP_vs_cortical');a=ax(f,'d',[119,62,41,29],'AIIMS patient scatter',['aiims','aiims_results']);cohort(a,'aiims','SP_vs_cortical')
    a=ax(f,'e',[193,62,40,29],'native AIIMS phenotype raincloud',['aiims_native']);s=data('aiims_native');order=['Immature_TALL','Cortical_TALL','Mature_TALL'];raincloud(a,s,'source_disease_state','ZEB1_log2_FPKM_plus1',order,[C['early'],C['cortical'],C['SP']]);counted_groups(a,s,'source_disease_state','ZEB1_log2_FPKM_plus1',order,['Imm.','Cort.','Mat.']);a.set_ylabel('ZEB1 log2(FPKM+1)');tag(a,'AIIMS · original phenotypes')
    a=ax(f,'f',[26,19,100,24],'India complete patient landscape',['india','india_results'],'MAJOR_SUPPORT');d=data('india').sort_values('ZEB1_log2_counts_plus1');m=old.zcols(d[['ZEB1_log2_counts_plus1','cortical_vs_early_raw','SP_vs_cortical_raw']]).T;m.index=['ZEB1','Cortical / early','SP / cortical'];m.columns=range(1,26);heat(f,a,m,'Within-India z',rotation=0);a.set_xticks([0,4,9,14,19,24],[1,5,10,15,20,25]);a.set_xlabel('All 25 patients, ordered by ZEB1');tag(a,'India · both source estimates uncertain')
    a=ax(f,'g',[180,19,53,25],'AIIMS phenotype-residual patient scatter',['aiims_native','aiims_conditional'],'MAJOR_SUPPORT');residualscatter(a,'SP_vs_cortical');f.legend(*a.get_legend_handles_labels(),frameon=False,ncol=3,loc='lower right',bbox_to_anchor=(.98,.003),handletextpad=.25,columnspacing=.8,markerscale=.75,fontsize=7.5);save(f,'Fig5')
def survival(a,end):
    d=data('km_'+end.lower());p=data('bulk');times=[0,2,4,6,8,10]
    for j,(k,c) in enumerate(zip(['T1_low','T2_mid','T3_high'],[C['early'],'#888888',C['ZEB1']])):
        s=d[d.tertile==k]
        a.step(s.time/365.25,s.surv,where='post',color=c,label=['Low','Middle','High'][j],lw=1.15)
        q=p[p.zeb1_tertile==k]
        for x in times:
            a.text(x,-.28-j*.11,str(int((q[end.lower()+'_time']>=x*365.25).sum())),
                   transform=a.get_xaxis_transform(),ha='center',fontsize=8)
        a.text(-.035,-.28-j*.11,['Low','Mid','High'][j],transform=a.transAxes,ha='right',fontsize=8)
    # A shared, explicitly ticked range keeps every observed displayed curve.
    visible=d[d.time/365.25<=10.2]
    assert visible.surv.min()>=.75
    a.set_xlim(0,10.2);a.set_ylim(.75,1.01);a.set_xticks(times)
    a.set_yticks([.75,.80,.85,.90,.95,1.00])
    a.set_xlabel('Years from diagnosis');a.xaxis.set_label_coords(.5,-.63)
    a.set_ylabel(end+' probability');a.legend(frameon=False,ncol=3,loc='lower left',columnspacing=1.4)
    tag(a,'AALL0434 · ZEB1 tertiles · n=1,309')

def Fig6():
    f=start('Fig6')
    a=ax(f,'a',[20,116,92,33],'observed patient OS',['km_os','km_counts_os','bulk'],'HERO');survival(a,'OS')
    a=ax(f,'b',[141,116,92,33],'observed patient EFS',['km_efs','km_counts_efs','bulk'],'HERO');survival(a,'EFS')

    a=ax(f,'c',[43,59,81,25],'diagnosis hazard effects',['cox_OS_mrd','cox_EFS_mrd'])
    vals=[];lo=[];hi=[];labs=[];cols=[]
    for end in ['OS','EFS']:
        z=data('cox_'+end+'_mrd').set_index('term').loc['ZEB1']
        vals.append(z['exp(coef)']);lo.append(z['exp(coef) lower 95%']);hi.append(z['exp(coef) upper 95%'])
        labs.append(end+' diagnosis');cols.append(C[end])
    effect(a,labs,vals,lo,hi,'HR per log2(TPM+1)',cols,1)
    a.set_ylim(1.6,-.6);a.set_xlim(.52,1.10);a.set_xticks([.6,.8,1.0])
    tag(a,'Diagnosis OS/EFS · MRD adjusted')

    a=ax(f,'d',[160,59,73,25],'observed discovery MRD joint distribution',['bulk','main_mrd'],'MAJOR_SUPPORT')
    d=data('bulk');s=d[d.log_mrd.notna()];assert len(s)==1302
    a.set_position([160/W,59/H,64/W,22/H])
    a.scatter(s.ZEB1,s.log_mrd,s=6,color=C['ZEB1'],alpha=.35,edgecolors='none')
    a.set_xlabel('ZEB1 log2(TPM+1)');a.set_ylabel('Day-29 MRD: ln(1 + %)')
    top=f.add_axes([160/W,82/H,64/W,3/H],sharex=a)
    top.hist(s.ZEB1,bins=30,color=C['ZEB1'],alpha=.55)
    right=f.add_axes([226/W,59/H,7/W,22/H],sharey=a)
    right.hist(s.log_mrd,bins=25,orientation='horizontal',color=C['ZEB1'],alpha=.55)
    for marginal in [top,right]:marginal.set_axis_off()
    f.text(160/W,88/H,'AALL0434 · day 29 · n=1,302',fontsize=8,va='bottom')

    # One lettered panel compares two quantities, with no invented C-value CI.
    a=ax(f,'e',[43,16,81,26],'paired apparent C points and delta-C intervals',['nested'],'MAJOR_SUPPORT')
    a.set_position([43/W,16/H,36/W,26/H])
    delta=f.add_axes([91/W,16/H,33/W,26/H],sharey=a)
    delta.spines[['top','right','left']].set_visible(False)
    s=data('nested').set_index('pair').loc[['OS_clinical_mrd_plus_ZEB1','EFS_clinical_mrd_plus_ZEB1']]
    for row,(end,z) in enumerate(zip(['OS','EFS'],s.itertuples())):
        y=1-row
        a.plot([z.C_reduced,z.C_full],[y,y],color='#999999',lw=1.1,zorder=1)
        a.plot(z.C_reduced,y,'o',mfc='white',mec=C['early'],mew=1.2,ms=5,zorder=3)
        a.plot(z.C_full,y,'o',mfc=C['ZEB1'],mec=C['ZEB1'],ms=5,zorder=3)
        a.text(z.C_reduced,y+.23,f'{z.C_reduced:.4f}',ha='center',va='bottom',fontsize=8)
        a.text(z.C_full,y-.23,f'{z.C_full:.4f}',ha='center',va='top',fontsize=8)
        delta.plot([z.boot_dC_lo,z.boot_dC_hi],[y,y],color=C[end],lw=1.2)
        delta.plot([z.boot_dC_lo,z.boot_dC_lo],[y-.065,y+.065],color=C[end],lw=.8)
        delta.plot([z.boot_dC_hi,z.boot_dC_hi],[y-.065,y+.065],color=C[end],lw=.8)
        delta.plot(z.delta_C,y,'o',color=C[end],ms=4.5)
        delta.text(.024,y+.23,f'{z.delta_C:+.4f}',ha='right',va='bottom',fontsize=8)
    a.set_xlim(.710,.750);a.set_ylim(-.65,1.65)
    a.set_xticks([.71,.73,.75]);a.set_yticks([1,0],['OS','EFS'])
    a.set_xlabel('Harrell C');a.spines['left'].set_visible(False);a.tick_params(axis='y',length=0)
    tag(a,'Apparent C · n=1,302')
    delta.axvline(0,color='black',lw=.65,ls='--')
    delta.set_xlim(-.030,.025);delta.set_xticks([-.02,0,.02],['−0.02','0','0.02'])
    delta.tick_params(axis='y',left=False,labelleft=False);delta.set_xlabel('ΔC')
    tag(delta,'Added ZEB1 · 95% CI')
    f.legend([Line2D([0],[0],marker='o',ls='',mfc='white',mec=C['early'],mew=1.2,ms=4),
              Line2D([0],[0],marker='o',ls='',color=C['ZEB1'],ms=4)],
             ['Clinical + MRD','+ ZEB1'],frameon=False,ncol=2,loc='center',
             bbox_to_anchor=(83.5/W,3.5/H),fontsize=7.5,columnspacing=1.4,handletextpad=.3)

    a=ax(f,'f',[160,16,73,26],'independent end-induction MRD patients',['stjude_end_patients','stjude_mrd'],'MAJOR_SUPPORT')
    d=data('stjude_end_patients');assert len(d)==93
    for protocol,c in [('T15',C['early']),('T16',C['cortical'])]:
        s=d[d.Protocol==protocol]
        a.scatter(s.ZEB1,s.MRD_log10_reported_floor,s=13,color=c,alpha=.75,edgecolors='none',label=protocol)
    a.set_xlabel('ZEB1 log2(FPKM+1)');a.set_ylabel('log10 MRD (%)')
    tag(a,'St Jude · end induction · n=93')
    a.legend(frameon=False,ncol=2,loc='upper left',fontsize=7.5,markerscale=.75)
    save(f,'Fig6')

def SFig1():
    f=start('SFig1');genes=sorted(t.normalized_pathway('hallmark_raw').pathway.str.replace('HALLMARK_','').tolist())
    for l,x,part in [('a',49,genes[:25]),('b',127,genes[25:])]:
        a=ax(f,l,[x,20,18,126],'complete discovery pathways',['hallmark_raw','hallmark_residual'],'HERO','Fig1');pathpoints(a,part);a.set_yticklabels([t.short_pathway('HALLMARK_'+g)[:25] for g in part]);a.tick_params(axis='y',labelsize=8)
    genes=['ZEB1','ZEB2','BCL11B','LYL1','CD34','CD1A','RAG1','RAG2','TCF7','GATA3','CD3D','CD8B']
    for l,k,col,r in [('c','normal195','facs_stage',[164,87,30,58]),('d','normal206','inferred_stage',[215,87,17,58])]:
        d=data(k).set_index(col)
        if k=='normal195':d=d.reindex(['DN1','DN2','DN3','ISP','DP_CD3neg','DP_CD3pos','CD4SP','CD8SP'])
        m=d[['g_'+g for g in genes]].T;m.index=genes;m=m.div(m.max(axis=1),axis=0);a=ax(f,l,r,'normal marker profiles',[k],'MAJOR_SUPPORT','Fig2');heat(f,a,m,'Relative mean',0,1,'YlOrBr',cbar=False);a.set_xticklabels(['DN1','DN2','DN3','ISP','DP CD3−','DP CD3+','CD4 SP','CD8 SP'] if k=='normal195' else ['DN early','DP imm.','DP mat.'],rotation=55,ha='right');a.tick_params(axis='x',labelsize=7.5);tag(a,'GSE195812' if k=='normal195' else 'GSE206710')
    d=data('normal_heldout')
    for l,y,k in [('e',50,'cortical_vs_early'),('f',18,'SP_vs_cortical')]:
        a=ax(f,l,[173,y,59,15],'heldout normal donor contrasts',['normal_heldout'],'COMPACT_SUPPORT','Fig2');s=d[d.axis==k];a.scatter(s.heldout_difference,range(len(s)),color=C['cortical' if l=='e' else 'SP'],s=24);a.set_yticks(range(len(s)),s.heldout_donor.str.replace('TSC-',''));a.set_xlabel('Held-out paired score contrast');a.axvline(0,color='black',lw=.6);tag(a,'Cortical / early' if l=='e' else 'SP / cortical')
    pathwaykey(f,(.35,.018));f.legend([plt.Rectangle((0,0),1,1,color=plt.get_cmap('YlOrBr')(0.)),plt.Rectangle((0,0),1,1,color=plt.get_cmap('YlOrBr')(1.))],['0','1: gene-relative mean'],frameon=False,ncol=2,loc='lower right',bbox_to_anchor=(.99,.005),fontsize=7.5,columnspacing=.7);save(f,'SFig1')
def SFig2():
    f=start('SFig2');d=data('sc_patient');genes=['ZEB1','ZEB2','BCL11B','MKI67','CD34','SPINK2','LYL1','CD1A','CD3D','CD8A','RAG1','TCF7'];m=old.zcols(d.set_index('patient')[genes]).T;a=ax(f,'a',[20,73,126,74],'malignant patient marker landscape',['sc_patient'],'HERO','Fig3');heat(f,a,m,'Across-patient gene z',rotation=90);tag(a,'GSE248287 · all 15 diagnostic patients');a=ax(f,'b',[181,96,51,49],'patient BCL11B context',['sc_patient','sc_assoc'],'MAJOR_SUPPORT','Fig3');scatter(a,d,'BCL11B','ZEB1',C['ZEB1'],xlabel='BCL11B log2(CPM+1)',ylabel='ZEB1 log2(CPM+1)');tag(a,'n=15 · source ρ=−.268; q=.533')
    for l,x,g in [('c',18,'BCL11B'),('d',81,'MKI67'),('e',144,'SPINK2'),('f',207,'CD1A')]:a=ax(f,l,[x,16,26,29],'paired malignant marker state',['sc_state'],'COMPACT_SUPPORT','Fig3');pair(a,g)
    save(f,'SFig2')
def SFig3():
    f=start('SFig3');d=data('mutations').sort_values(['subtype','ZEB1','sample_id']);r=data('mutation_assoc').sort_values('gene');genes=r.gene.tolist();a=ax(f,'a',[25,16,111,132],'complete genetic patient landscape',['mutations'],'HERO','Fig4');a.imshow(d[genes].T,aspect='auto',cmap=ListedColormap(['#F1F1F1',C['early']]),vmin=0,vmax=1,interpolation='nearest',rasterized=True);a.set_yticks(range(len(genes)),genes);a.set_xticks([]);a.set_xlabel('All 1,309 patients: subtype, then ZEB1');tag(a,'All 37 supplied alteration fields');offset=0
    for k,s in d.groupby('subtype',sort=False):a.axvline(offset-.5,color='black',lw=.3);offset+=len(s)
    for l,y,s in [('b',92,r.iloc[:19]),('c',18,r.iloc[19:])]:
        a=ax(f,l,[188,y,44,56],'complete lesion-expression family',['mutation_assoc'],'MAJOR_SUPPORT','Fig4');a.scatter(s.delta,range(len(s)),c=np.where(s.mw_q<.05,C['ZEB1'],'#999999'),s=25);a.set_yticks(range(len(s)),[f'{v.gene} ({int(v.n_mut)})' for _,v in s.iterrows()]);a.invert_yaxis();a.axvline(0,color='black',lw=.6);a.set_xlabel('Altered − absent median ZEB1');a.tick_params(axis='y',labelsize=8)
    f.legend([plt.Rectangle((0,0),1,1,color=C['early']),plt.Rectangle((0,0),1,1,color='#F1F1F1'),Line2D([0],[0],marker='o',ls='',color=C['ZEB1']),Line2D([0],[0],marker='o',ls='',color='#999999')],['Altered','Absent','37-field q < .05','q ≥ .05'],frameon=False,ncol=4,loc='lower center',bbox_to_anchor=(.5,.008));save(f,'SFig3')
def SFig4():
    f=start('SFig4');d=data('aiims_native');order=['Immature_TALL','Cortical_TALL','Mature_TALL']
    for l,x,k in [('a',20,'cortical_vs_early_raw'),('b',85,'SP_vs_cortical_raw')]:
        a=ax(f,l,[x,110,32,36],'native phenotype fixed-score distributions',['aiims_native'],'COMPACT_SUPPORT','Fig5');sourcebox(a,d,'source_disease_state',k,order,[C['early'],C['cortical'],C['SP']]);a.set_xticklabels(['Immature','Cortical','Mature'],rotation=50,ha='right');a.set_ylabel('Cortical / early score' if l=='a' else 'SP / cortical score');tag(a,'AIIMS n=35')
    a=ax(f,'c',[150,110,32,36],'complementary phenotype-residual patients',['aiims_native','aiims_conditional'],'COMPACT_SUPPORT','Fig5');residualscatter(a,'cortical_vs_early');a.set_xlabel('Cortical / early\nrank residual')
    a=ax(f,'d',[215,110,18,36],'NOPHO patient BCL11B context',['nopho_genes','nopho_gene_assoc'],'COMPACT_SUPPORT','Fig5');d=data('nopho_genes');scatter(a,d,'ZEB1','BCL11B',C['ZEB1'],xlabel='ZEB1\n(author scale)',ylabel='BCL11B (author scale)');tag(a,'NOPHO n=108\nρ=.27; q=.0078')
    r=data('lim_results');labs=['D: SP/C','D: C/E','V: SP/C','V: C/E'];a=ax(f,'e',[33,20,40,42],'Lim primary bootstrap uncertainty',['lim_results'],'MAJOR_SUPPORT','Fig3');s=r[r.model=='primary_spearman'].sort_values(['cohort','axis']);effect(a,labs,s.effect.tolist(),s.ci_low.tolist(),s.ci_high.tolist(),'Spearman ρ (patient CI)');a.set_xlim(-1,1);tag(a,'Lim: D=21 · V=28')
    a=ax(f,'f',[111,20,40,42],'Lim adjusted individual-deletion range',['lim_results'],'MAJOR_SUPPORT','Fig3');s=r[r.model=='cycle_depth_rank_HC3'].sort_values(['cohort','axis']);effect(a,labs,s.effect.tolist(),s.loo_min.tolist(),s.loo_max.tolist(),'Rank slope (deletion range)');tag(a,'Cycle/depth adjusted')
    a=ax(f,'g',[185,21,48,41],'complete Lim patient capture context',['lim_scores','lim_units'],'MAJOR_SUPPORT','Fig3');d=data('lim_scores').sort_values(['cohort_analysis','ZEB1_log2_CPM_plus1']);m=d[['zero_fraction_shared','n_cells','cycle_fraction']].copy();m.n_cells=np.log1p(m.n_cells);m=old.zcols(m).T;m.index=['Shared genes:\nzero fraction','log(1 +\nmalignant cells)','Measured\ncycle fraction'];m.columns=range(1,50);heat(f,a,m,'Display z',rotation=0);a.set_xticks([0,48],[1,49]);a.axvline(20.5,color='black',lw=.8);a.set_xlabel('Discovery 21 | validation 28');a.tick_params(axis='y',labelsize=7.5);tag(a,'All 49 eligible Lim patients');f.legend(*f.axes[2].get_legend_handles_labels(),frameon=False,ncol=3,loc='center',bbox_to_anchor=(.61,.52),fontsize=8);save(f,'SFig4')
def SFig5():
    f=start('SFig5')
    for l,end,x in [('a','OS',19),('b','EFS',103)]:
        v=js('landmark')['results'][end];d=pd.DataFrame(v['HR_by_elapsed_time']);a=ax(f,l,[x,111,59,35],'landmark time-specific association',['landmark'],'MAJOR_SUPPORT','Fig6');xx=d.elapsed_days/365.25;a.plot(xx,d.HR,'o-',color=C[end],ms=4);a.fill_between(xx,d.lower95,d.upper95,color=C[end],alpha=.18);a.axhline(1,color='black',lw=.6,ls='--');a.set_xlabel('Years after day-29 landmark');a.set_ylabel(end+' pointwise HR');tag(a,f'{end} · block BH2 q={v["ZEB1_block_LRT"]["BH_q_two_endpoints"]:.3g}')
    a=ax(f,'c',[199,112,33,32],'complete apparent-discrimination intervals',['nested'],'COMPACT_SUPPORT','Fig6');d=data('nested');s=d[d.boot_dC_lo.notna()];key={'OS_clinical_mrd_plus_ZEB1':'OS + ZEB1','OS_mrd_increment_noZEB1':'OS + MRD','EFS_clinical_mrd_plus_ZEB1':'EFS + ZEB1'};labels=[key[p] for p in s.pair];effect(a,labels,s.delta_C.tolist(),s.boot_dC_lo.tolist(),s.boot_dC_hi.tolist(),'Apparent ΔC');tag(a,'200 patient resamples')
    d=data('drug_slots');drugs=d.drug.unique().tolist()
    for l,x,endpoint in [('d',35,'native_normalized'),('e',115,'fractional_midrank')]:
        a=ax(f,l,[x,23,34,63],'full drug-response estimation landscape',['drug_slots'],'HERO','Fig6');m=d[d.endpoint==endpoint].pivot(index='drug',columns='model',values='beta').reindex(drugs);heat(f,a,m,'ZEB1 slope',-.3,.3);a.set_xticklabels(['Base','+ proxy'],rotation=40,ha='right');a.tick_params(axis='y',labelsize=8);tag(a,'Normalized LC50' if l=='d' else 'Fractional midrank')
    a=ax(f,'f',[175,56,58,30],'drug patient availability',['drug_slots'],'COMPACT_SUPPORT','Fig6');s=d.drop_duplicates('drug');a.bar(range(len(s)),s.n,color=C['early'],width=.8);a.set_xticks(range(len(s)),s.drug,rotation=90);a.tick_params(axis='x',labelsize=7.5);a.set_ylabel('Complete patients');a.axhline(30,color='black',lw=.6,ls='--');tag(a,'Independent St Jude · all 18 drugs')
    a=ax(f,'g',[186,15,46,13],'same-discovery-cohort MRD slopes',['main_mrd'],'COMPACT_SUPPORT','Fig6');d=data('main_mrd');s=d[d.endpoint.str.contains('continuous')];effect(a,['Fixed proxy','+ BMP17 / ancestry'],s.target_beta.tolist(),s.CI95_low.tolist(),s.CI95_high.tolist(),'ln(1 + MRD %) slope');tag(a,'AALL0434 · MRD models');f.legend([plt.Rectangle((0,0),1,1,color='#DADADA')],['Gray: unestimated (40/72 drug slots)'],frameon=False,loc='lower center',bbox_to_anchor=(.38,.008));save(f,'SFig5')

def SFig6():
    f=start('SFig6');d=analysisdata('NORMAL30_DONOR_STATE_SCORES_AND_ZEB1.tsv','source_data');donors=sorted(d.donor.unique());palette=plt.get_cmap('tab10')(range(5));states=['Early_DN','DP','CD4SP','CD8SP']
    a=ax(f,'a',[20,93,90,54],'depth-corrected normal donor trajectories',['normal_donor_expression'],'HERO','Fig2')
    for color,donorid in zip(palette,donors):
        z=d[d.donor==donorid].set_index('mapped_state').reindex(states);y=z.ZEB1_log2_CPM_plus1.where(z.eligible_ge30);a.plot(range(4),y,'o-',color=color,ms=3,label=donorid.replace('TSC-',''))
    a.set_xticks(range(4),['Early DN','DP','CD4 SP','CD8 SP']);a.set_ylabel('ZEB1 log2(CPM + 1)');a.legend(frameon=False,ncol=5,loc='upper center',bbox_to_anchor=(.5,-.18),fontsize=8,columnspacing=.9);tag(a,'Five donors · original full-library depth correction')
    a=ax(f,'b',[155,93,76,54],'normal paired-score geometry',['normal_donor_expression'],'MAJOR_SUPPORT','Fig2')
    for state,color in zip(states+['DN','T_entry'],[C['early'],C['cortical'],C['SP'],'#8C66A5','#C25E64','#A58B3F']):
        z=d[(d.mapped_state==state)&d.eligible_ge30];a.scatter(z.cortical_vs_early,z.SP_vs_cortical,s=25,color=color,label=state.replace('_',' '))
    a.set_xlabel('Cortical / early fixed score');a.set_ylabel('SP / cortical fixed score');a.legend(frameon=False,ncol=3,loc='upper center',bbox_to_anchor=(.5,-.15),fontsize=8,columnspacing=.6,handletextpad=.3);tag(a,'29 eligible profiles · five repeated donors')
    a=ax(f,'c',[30,17,37,45],'normal holdout member retention',['normal_member_retention'],'MAJOR_SUPPORT','Fig2');r=analysisdata('NORMAL18_HELDOUT_MEMBERSHIP_STABILITY.tsv');g=r.groupby(['axis','heldout_donor'],sort=False).fraction_fixed_retained.agg(['min','max']).reset_index()
    for j,v in enumerate(g.itertuples()):
        color=C['cortical' if v.axis=='cortical_vs_early' else 'SP'];a.plot([v.min,v.max],[j,j],color=color,lw=2);a.plot((v.min+v.max)/2,j,'o',color=color,ms=3)
    a.set_yticks(range(len(g)),[('C/E ' if v.axis=='cortical_vs_early' else 'SP/C ')+v.heldout_donor.replace('TSC-CT','CT') for v in g.itertuples()]);a.invert_yaxis();a.set_xlim(.85,1.015);a.set_xticks([.9,.95,1]);a.set_xlabel('Fixed members retained');tag(a,'Existing donor holdouts')
    r=analysisdata('BULK200_FIXED_MEMBER_LEAVE_ONE_OUT.tsv');a=ax(f,'d',[109,49,47,17],'all single-member bulk omission effects',['bulk_member_omission'],'COMPACT_SUPPORT','Fig1')
    for axis,color in [('cortical_vs_early',C['cortical']),('SP_vs_cortical',C['SP'])]:
        z=r[r.axis==axis];a.scatter(range(1,101),z.adjusted_global_rank_slope,s=5,color=color,alpha=.6)
    a.set_xlabel('Omitted member (1–100 per axis)');a.set_ylabel('Adjusted β');a.axhline(0,color='black',lw=.5);tag(a,'All 200 member omissions')
    r=analysisdata('NORMAL200_FIXED_MEMBER_LEAVE_ONE_OUT.tsv');a=ax(f,'e',[109,14,47,14],'normal minimum donor contrast under member omission',['normal_member_omission'],'COMPACT_SUPPORT','Fig2')
    for axis,color in [('cortical_vs_early',C['cortical']),('SP_vs_cortical',C['SP'])]:
        z=r[r.axis==axis];a.scatter(range(1,101),z.minimum_paired_difference,s=5,color=color,alpha=.6)
    a.set_xlabel('Omitted member (1–100 per axis)');a.set_ylabel('Min. donor contrast');a.axhline(0,color='black',lw=.5);tag(a,'All donor directions preserved')
    r=analysisdata('BULK9_EXISTING_DONOR_FOLD_SCORE_SENSITIVITY.tsv');a=ax(f,'f',[203,17,30,49],'bulk effects using original donor-fold members',['bulk_normal_fold_effects'],'MAJOR_SUPPORT','Fig1')
    labels=[]
    for j,v in enumerate(r.itertuples()):
        color=C['cortical' if v.axis=='cortical_vs_early' else 'SP'];a.plot([v.ci_low,v.ci_high],[j,j],color=color,lw=.85);a.plot(v.effect,j,'o',color=color,ms=3);labels.append(('C/E ' if v.axis=='cortical_vs_early' else 'SP/C ')+v.heldout_normal_donor.replace('TSC-CT','CT'))
    a.set_yticks(range(len(r)),labels);a.invert_yaxis();a.axvline(0,color='black',lw=.5,ls='--');a.set_xlabel('Adjusted β · 95% CI');tag(a,'Nine donor folds');save(f,'SFig6')

FUNCTIONS={v.__name__:v for v in [Fig1,Fig2,Fig3,Fig4,Fig5,Fig6,SFig1,SFig2,SFig3,SFig4,SFig5,SFig6]}
names=sys.argv[1:] or list(FUNCTIONS)
for n in names:FUNCTIONS[n]()
g=pd.DataFrame(PANELS);p=F/'PANEL_GRAMMAR.tsv'
if p.exists():
    previous=pd.read_csv(p,sep='\t');previous=previous[previous.figure.isin(FUNCTIONS)&~previous.figure.isin(names)];g=pd.concat([previous,g],ignore_index=True)
g.to_csv(p,sep='\t',index=False)

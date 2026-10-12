"""Render the selected observational story from exact frozen results.

No inference, feature selection, cutpoint fitting, clustering or embedding.
Reuses the project's established Matplotlib render adapter. Only documented
display transformations are applied. Run with names to resume a bounded batch.
"""
from pathlib import Path
import importlib.util, sys, json, string
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.colors import ListedColormap

B=Path(__file__).resolve().parent.parent
sys.path.insert(0,str(Path(__file__).resolve().parent/'vendor'))
from adjustText import adjust_text
spec=importlib.util.spec_from_file_location('accepted_render',Path(__file__).with_name('render_template.py'))
t=importlib.util.module_from_spec(spec);spec.loader.exec_module(t)
t.BASE=B;t.I=B/'source_data';t.OUT=B/'figures';t.OUT.mkdir(exist_ok=True)
data=t.data;js=t.js;C=t.C;boxes=t.boxes;scatter=t.scatter;effect=t.effect
for p in (B/'fonts').glob('*.ttf'):t.fontManager.addfont(str(p))
plt.rcParams.update({'font.size':8,'axes.labelsize':8,'xtick.labelsize':7.5,'ytick.labelsize':7.5,'legend.fontsize':7.5})
def start(name,h=140,w=210):
    t.CURRENT=name
    return plt.figure(figsize=(w/25.4,h/25.4))
def ax(f,letter,rect,kind,sources,role='COMPACT_SUPPORT',parent=''):
    weight={'HERO':'HERO','MAJOR_SUPPORT':'L','COMPACT_SUPPORT':'M'}[role]
    a=t.ax(f,letter,rect,kind,sources,weight,parent)
    a.texts[0].set_transform(f.transFigure)
    a.texts[0].set_position((max(.012,rect[0]-.045),min(.966,rect[1]+rect[3]+.03)))
    return a
def tag(a,s):a.text(.02,1.035,s,transform=a.transAxes,fontsize=7.5,va='bottom')
def save(f,n):
    p=t.OUT/n;p.mkdir(exist_ok=True)
    for ext in ['pdf','svg','png']:f.savefig(p/(n+'.'+ext),dpi=300)
    plt.close(f);t.FIGS.append(n);print(n+' RENDERED',flush=True)
def hm(f,a,m,label,vmin=-2,vmax=2,cmap='RdBu_r',cbar=True):
    cm=plt.get_cmap(cmap).copy();cm.set_bad('#DADADA')
    im=a.imshow(m.to_numpy(dtype=float),aspect='auto',interpolation='nearest',cmap=cm,vmin=vmin,vmax=vmax,rasterized=True)
    a.set_yticks(range(len(m)),m.index);a.set_xticks(range(len(m.columns)),m.columns,rotation=60,ha='right')
    if cbar:f.colorbar(im,ax=a,fraction=.03,pad=.03).set_label(label)
    return im
def zcols(m):return (m-m.mean())/m.std(ddof=0)
def cohort_plot(f,l,rect,k,axis):
    d,y=t.cohort(k);color=C['cortical' if axis=='cortical_vs_early' else 'SP']
    a=ax(f,l,rect,'patient scatter',[k,k+'_results'])
    xlabel='Cortical / early score' if axis=='cortical_vs_early' else 'SP / cortical score'
    yl={'nopho':'ZEB1 log2(CPM+1)','aiims':'ZEB1 log2(FPKM+1)','india':'ZEB1 log2(counts+1)'}[k]
    scatter(a,d,axis+'_raw',y,color,xlabel=xlabel,ylabel=yl)
    r=next(x for x in js(k+'_results')['results'] if x['axis']==axis)
    tag(a,f"{k.upper()} n={len(d)}\nρ={r['spearman_rho']:.2f}; q={r['BH2_q']:.2g}")
    return a

def Fig1():
    f=start('Fig1');d=data('bulk').merge(data('maturity'),on='sample_id',validate='one_to_one')
    order=d.groupby('subtype').ZEB1.median().sort_values().index.tolist()
    ordered=d.assign(subtype=pd.Categorical(d.subtype,categories=order,ordered=True)).sort_values(['subtype','ZEB1'])
    a=ax(f,'a',[.08,.77,.87,.16],'patient landscape',['bulk'],'HERO')
    offset=0
    for j,g in enumerate(order):
        s=ordered[ordered.subtype==g];x=np.arange(len(s))+offset
        a.scatter(x,s.ZEB1,s=3,color=plt.get_cmap('tab20')(j%20),rasterized=False)
        a.axvline(offset-.5,color='#BBBBBB',lw=.25);offset+=len(s)
    a.set_xticks([]);a.set_ylabel('ZEB1\nlog2(TPM+1)');a.set_xlabel('All 1,309 patients: subtype blocks, then expression; block key in b')
    a=ax(f,'b',[.16,.12,.24,.49],'subtype distribution',['bulk'],'MAJOR_SUPPORT')
    for j,g in enumerate(order):
        s=d[d.subtype==g];a.scatter(s.ZEB1,np.full(len(s),j),s=2,alpha=.20,color=plt.get_cmap('tab20')(j%20),rasterized=False)
        q=s.ZEB1.quantile([.25,.5,.75]);a.plot([q.iloc[0],q.iloc[2]],[j,j],color='black',lw=1.4);a.plot(q.iloc[1],j,'|',color='black',ms=6)
    a.set_yticks(range(len(order)),[g+f' ({sum(d.subtype==g)})' for g in order]);a.invert_yaxis();a.set_xlabel('ZEB1 log2(TPM+1)');tag(a,'Supplied 17-group annotation · all patients')
    a=ax(f,'c',[.56,.43,.16,.20],'ETP distribution',['bulk']);boxes(a,d,'etp3','ZEB1',['ETP','Near-ETP','Non-ETP'],[C['early'],C['cortical'],C['SP']]);a.set_ylabel('ZEB1 log2(TPM+1)')
    a=ax(f,'d',[.81,.43,.15,.20],'developmental scatter',['bulk','maturity']);scatter(a,d,'cortical_vs_early','ZEB1',C['cortical'],xlabel='Cortical / early (z)',ylabel='ZEB1 log2(TPM+1)')
    a=ax(f,'e',[.56,.12,.16,.20],'developmental scatter',['bulk','maturity']);scatter(a,d,'SP_vs_cortical','ZEB1',C['SP'],xlabel='SP / cortical (z)',ylabel='ZEB1 log2(TPM+1)')
    a=ax(f,'f',[.84,.12,.12,.20],'pathway landscape',['hallmark_raw','hallmark_residual']);pathway_compare(a,['E2F_TARGETS','G2M_CHECKPOINT','MYC_TARGETS_V1','PROTEIN_SECRETION','ANGIOGENESIS'])
    f.legend([Line2D([0],[0],color=C[x],marker=m,ls='') for x,m in [('early','o'),('ZEB1','s')]],['Marginal (display −NES)','Subtype residual'],loc='lower right',bbox_to_anchor=(.98,.002),frameon=False,ncol=2)
    save(f,'Fig1')
def pathway_compare(a,genes):
    for j,g in enumerate(genes):
        for n,k in enumerate(['hallmark_raw','hallmark_residual']):
            r=t.normalized_pathway(k).set_index('pathway').loc['HALLMARK_'+g]
            v=-r.NES if n==0 else r.NES;col=C['early' if n==0 else 'ZEB1']
            a.scatter(v,j+(n-.5)*.18,s=22,marker='o' if n==0 else 's',facecolors=col if float(r.q)<.05 else 'white',edgecolors=col,lw=.7)
    a.set_yticks(range(len(genes)),[g.replace('_TARGETS','').replace('_CHECKPOINT','').replace('PROTEIN_SECRETION','Secretion').replace('ANGIOGENESIS','Angiogenesis') for g in genes]);a.invert_yaxis();a.axvline(0,color='black',lw=.5);a.set_xlabel('NES: + = ZEB1-high')

def Fig2():
    f=start('Fig2');d=data('normal_umap');labels=data('normal_labels').set_index('author_cell_type').display_label.to_dict();groups=list(d.author_cell_type.unique());palette=plt.get_cmap('tab20b')(np.arange(len(groups))%20)
    a=ax(f,'a',[.09,.66,.34,.28],'UMAP embedding',['normal_umap','normal_labels'],'HERO')
    handles=[];category_texts=[]
    for j,g in enumerate(groups):
        s=d[d.author_cell_type==g];a.scatter(s.UMAP1,s.UMAP2,s=.12,color=palette[j],linewidths=0,rasterized=True)
        handles.append(Line2D([0],[0],marker='o',ls='',color=palette[j],ms=3,label=f'{j+1} '+labels.get(g,g).replace('Treg-','Tr-').replace('Treg_','Tr_').replace('SP_fetal_ambiguous','SP fetal amb.')))
        category_texts.append(a.text(s.UMAP1.median(),s.UMAP2.median(),str(j+1),fontsize=7.5,color='black',ha='center',bbox={'facecolor':'white','alpha':.68,'edgecolor':'none','pad':.05}))
    a.set_xlabel('UMAP 1');a.set_ylabel('UMAP 2');a.set_xticks([]);a.set_yticks([]);tag(a,f'Frozen author atlas · {len(d):,} rows')
    np.random.seed(20261010)
    adjust_text(category_texts,ax=a,force_text=(.5,.5),force_static=(.2,.2),force_pull=(.015,.015),expand=(1.25,1.25),iter_lim=400,min_arrow_len=2,arrowprops={'arrowstyle':'-','color':'#666666','lw':.4})
    a=ax(f,'b',[.59,.66,.32,.28],'FeaturePlot',['normal_umap'],'HERO');t.umap(a,d,True)
    f.legend(handles=handles,loc='center',bbox_to_anchor=(.50,.515),ncol=7,frameon=False,columnspacing=.8,handletextpad=.2,labelspacing=.25)
    p=data('normal195').set_index('facs_stage').reindex(['DN1','DN2','DN3','ISP','DP_CD3neg','DP_CD3pos','CD4SP','CD8SP'])
    a=ax(f,'c',[.09,.15,.20,.23],'normal stage profile',['normal195']);a.plot(range(len(p)),p.g_ZEB1,'o-',color=C['ZEB1'],ms=3);a.set_xticks(range(len(p)),['DN1','DN2','DN3','ISP','DP−','DP+','CD4','CD8'],rotation=55,ha='right');a.set_ylabel('Mean normalized ZEB1');tag(a,'GSE195812 · FACS stage means')
    p=data('normal206').set_index('inferred_stage');a=ax(f,'d',[.38,.15,.14,.23],'normal stage profile',['normal206']);a.plot(range(len(p)),p.g_ZEB1,'o-',color=C['ZEB1']);a.set_xticks(range(len(p)),['DN early','DP imm.','DP mat.'],rotation=55,ha='right');a.set_ylabel('Mean normalized ZEB1');tag(a,'GSE206710')
    d=data('normal_state_scores');s=d.pseudobulk_id.str.split('|',expand=True);d['donor']=s[0];d['stage']=s[1]
    for letter,x,k in [('e',.65,'cortical_vs_early'),('f',.85,'SP_vs_cortical')]:
        a=ax(f,letter,[x,.15,.11,.23],'donor trajectories',['normal_state_scores','normal_heldout']);order=['Early_DN','DP','CD4SP','CD8SP']
        for j,(donor,p) in enumerate(d.groupby('donor')):
            v=p.set_index('stage')[k].reindex(order);a.plot(range(4),v,'o-',color=plt.get_cmap('tab10')(j),ms=3,alpha=.7,label=donor.replace('TSC-',''))
        a.set_xticks(range(4),['Early','DP','CD4','CD8'],rotation=55,ha='right');a.set_ylabel('Fixed rank score');tag(a,'Cortical / early' if letter=='e' else 'SP / cortical')
    f.text(.61,.045,'Lines: five source donors; missing early state retained.\nHeld-out qualification: 4 and 5 donors (Supplementary Table S2).',fontsize=7.5)
    f.legend(*a.get_legend_handles_labels(),frameon=False,ncol=5,loc='lower center',bbox_to_anchor=(.32,.002),columnspacing=.8,handlelength=1)
    save(f,'Fig2')

def Fig3():
    f=start('Fig3');d=data('lim_scores');r=data('lim_results');pieces=[]
    for k in ['discovery','validation']:
        s=d[d.cohort_analysis==k].sort_values('ZEB1_log2_CPM_plus1');m=zcols(s[['ZEB1_log2_CPM_plus1','cortical_vs_early','SP_vs_cortical']]);m.index=[k[0].upper()+':'+p for p in s.patient];pieces.append(m)
    m=pd.concat(pieces).T;m.index=['ZEB1','Cortical / early','SP / cortical']
    a=ax(f,'a',[.13,.80,.80,.13],'malignant patient landscape',['lim_scores'],'HERO');hm(f,a,m,'Within-set display z');a.tick_params(axis='x',rotation=90);a.axvline(20.5,color='black',lw=1);tag(a,'Lim diagnosis BM: discovery 21 | validation 28 · author malignant pseudobulks')
    a=ax(f,'b',[.09,.36,.24,.24],'paired patient state',['sc_state','sc_receipt'],'MAJOR_SUPPORT');p=data('sc_state').pivot(index='patient',columns='mapped_state',values='ZEB1')[['Early_DN','DP']].dropna()
    last=-np.inf
    for j,(patient,s) in enumerate(p.sort_values('DP').iterrows()):
        ly=max(s.DP,last+.044);last=ly
        a.plot([0,1],s.values,color='#999999');a.scatter([0,1],s.values,c=[C['early'],C['cortical']],s=20);a.annotate(patient,(1,s.DP),xytext=(1.12,ly),fontsize=7.5,arrowprops={'arrowstyle':'-','lw':.4,'color':'#777777'})
    a.set_xlim(-.1,1.6);a.set_xticks([0,1],['Early DN','DP']);a.set_ylabel('ZEB1 log2(CPM+1)');tag(a,'GSE248287 · 5 same-patient pairs\nExact two-sided P=0.0625')
    for letter,k,axis,rect in [('c','discovery','cortical_vs_early',[.50,.43,.16,.17]),('d','validation','cortical_vs_early',[.80,.43,.16,.17]),('e','discovery','SP_vs_cortical',[.50,.12,.16,.17]),('f','validation','SP_vs_cortical',[.80,.12,.16,.17])]:
        a=ax(f,letter,rect,'malignant patient scatter',['lim_scores','lim_results']);s=d[d.cohort_analysis==k];v=r[(r.cohort==k)&(r.axis==axis)&(r.model=='primary_spearman')].iloc[0];color=C['cortical' if axis=='cortical_vs_early' else 'SP'];scatter(a,s,axis,'ZEB1_log2_CPM_plus1',color,xlabel='Cortical / early' if axis=='cortical_vs_early' else 'SP / cortical',ylabel='ZEB1 log2(CPM+1)');tag(a,f"{k.title()} n={len(s)}\nρ={v.effect:.2f}; BH4 q={v.BH4_primary:.4g}")
    a=ax(f,'g',[.17,.10,.16,.15],'adjusted rank slope',['lim_results']);s=r[r.model=='cycle_depth_rank_HC3'].sort_values(['cohort','axis']);labs=[('D' if v.cohort=='discovery' else 'V')+': '+('C/E' if v.axis=='cortical_vs_early' else 'SP/C') for _,v in s.iterrows()];effect(a,labs,s.effect.tolist(),s.ci_low.tolist(),s.ci_high.tolist(),'Cycle/depth HC3 slope');tag(a,'Different estimand: rank OLS')
    save(f,'Fig3')

def Fig4():
    f=start('Fig4');d=data('mutations');genes=['CDKN2A','NOTCH1','TLX1','TLX3','DNM2','PHF6','IL7R','BCL11B','FLT3','NRAS','JAK3','MED12'];s=d[d.subtype=='TLX1'].sort_values(['ZEB1','sample_id'])
    a=ax(f,'a',[.10,.58,.34,.33],'focused Oncoprint',['mutations'],'HERO');t.onco(a,f,s,genes);tag(a,'TLX1 subtype · complete 71-patient group')
    for legend in list(f.legends):legend.remove()
    a.legend([plt.Rectangle((0,0),1,1,color=C['early']),plt.Rectangle((0,0),1,1,color='#EEEEEE')],['Altered','Not altered'],frameon=False,ncol=2,loc='upper center',bbox_to_anchor=(.5,1.20))
    a=ax(f,'b',[.61,.58,.34,.33],'genotype subtype prevalence',['mutations'],'HERO');m=d.groupby('subtype')[genes].mean().T;hm(f,a,m,'Altered fraction',0,1,'YlOrBr');a.set_xticklabels(m.columns,rotation=65,ha='right')
    for l,x,g in [('c',.10,'NOTCH1'),('d',.35,'CDKN2A'),('e',.61,'TLX1'),('f',.86,'MED12')]:
        a=ax(f,l,[x,.12,.12,.18],'lesion expression distribution',['mutations','mutation_assoc']);boxes(a,d,g,'ZEB1',[0,1],[C['early'],C['ZEB1']]);a.set_xticklabels(['Absent','Present'],rotation=35,ha='right');a.set_ylabel('ZEB1 log2(TPM+1)');v=data('mutation_assoc').set_index('gene').loc[g];tag(a,f'{g}\n37-test BH q={v.mw_q:.2g}')
    save(f,'Fig4')

def Fig5():
    f=start('Fig5');rows=[]
    for k in ['nopho','aiims','india']:
        for v in js(k+'_results')['results']:rows.append(dict(cohort=k.upper(),**v))
    d=pd.DataFrame(rows);m=d.pivot(index='cohort',columns='axis',values='spearman_rho').reindex(index=['NOPHO','AIIMS','INDIA'],columns=['cortical_vs_early','SP_vs_cortical'])
    a=ax(f,'a',[.12,.44,.23,.39],'cross-cohort relationship landscape',['nopho_results','aiims_results','india_results'],'HERO');im=hm(f,a,m,'Spearman ρ',-1,1,cbar=False);a.set_xticklabels(['Cortical / early','SP / cortical'],rotation=40,ha='right');a.set_yticklabels(['NOPHO (108)','AIIMS (35)','India (25)'])
    for i,k in enumerate(['NOPHO','AIIMS','INDIA']):
        for j,axis in enumerate(['cortical_vs_early','SP_vs_cortical']):
            v=d[(d.cohort==k)&(d.axis==axis)].iloc[0];a.text(j,i,f'{v.spearman_rho:+.2f}\n'+('BH2 q<.05' if v.BH2_q<.05 else 'uncertain'),ha='center',va='center',fontsize=8,color='black')
    f.colorbar(im,ax=a,fraction=.05,pad=.04).set_label('Spearman ρ')
    cohort_plot(f,'b',[.51,.65,.13,.20],'nopho','cortical_vs_early');cohort_plot(f,'c',[.80,.65,.13,.20],'nopho','SP_vs_cortical')
    cohort_plot(f,'d',[.51,.33,.13,.17],'aiims','SP_vs_cortical')
    a=ax(f,'e',[.80,.33,.13,.17],'native phenotype distribution',['aiims_native']);s=data('aiims_native');order=['Immature_TALL','Cortical_TALL','Mature_TALL'];order=[x for x in order if x in set(s.source_disease_state)] or list(s.source_disease_state.unique());boxes(a,s,'source_disease_state','ZEB1_log2_FPKM_plus1',order,[C['early'],C['cortical'],C['SP']]);a.set_xticklabels([x.replace('_TALL','') for x in order],rotation=50,ha='right');a.set_ylabel('ZEB1 log2(FPKM+1)');tag(a,'AIIMS native phenotype')
    a=ax(f,'f',[.17,.12,.17,.15],'uncertain cohort effects',['india_results']);s=d[d.cohort=='INDIA'].set_index('axis').loc[['cortical_vs_early','SP_vs_cortical']];effect(a,['Cortical / early','SP / cortical'],s.spearman_rho.tolist(),s.bootstrap95CI_low.tolist(),s.bootstrap95CI_high.tolist(),'India 25: ρ and patient CI');a.set_xlim(-1,1)
    a=ax(f,'g',[.58,.08,.32,.12],'native phenotype conditioned association',['aiims_conditional']);s=data('aiims_conditional');effect(a,['SP / cortical','Cortical / early'],s.conditional_rank_rho.tolist(),s.stratified_patient_bootstrap_CI95_low.tolist(),s.stratified_patient_bootstrap_CI95_high.tolist(),'AIIMS within-phenotype rank ρ');a.set_xlim(-1,1)
    save(f,'Fig5')

def survival_view(f,l,rect,end):
    a=ax(f,l,rect,'patient survival',['km_'+end.lower(),'km_counts_'+end.lower(),'bulk'],'HERO');d=data('km_'+end.lower());p=data('bulk');colors=[C['early'],'#888888',C['ZEB1']];times=[0,2,4,6,8,10]
    for j,(k,col) in enumerate(zip(['T1_low','T2_mid','T3_high'],colors)):
        s=d[d.tertile==k];a.step(s.time/365.25,s.surv,where='post',color=col,label=k.replace('_',' '));q=p[p.zeb1_tertile==k];v=[sum(q[end.lower()+'_time']>=x*365.25) for x in times]
        for x,n in zip(times,v):a.text(x,-.30-j*.115,str(n),transform=a.get_xaxis_transform(),ha='center',color='black',fontsize=7.5)
        a.text(-.06,-.30-j*.115,['Low','Mid','High'][j],transform=a.transAxes,ha='right',fontsize=7.5)
    a.set_xlim(0,10.3);a.set_ylim(.5,1.01);a.set_xlabel('Years from diagnosis',labelpad=2);a.xaxis.set_label_coords(.5,-.72);a.set_ylabel(end+' probability');a.legend(frameon=False,ncol=3,loc='lower left');tag(a,'AALL0434 · fixed expression tertiles')
def Fig6():
    f=start('Fig6',144,216);survival_view(f,'a',[.10,.72,.33,.21],'OS');survival_view(f,'b',[.63,.72,.33,.21],'EFS')
    a=ax(f,'c',[.18,.29,.20,.22],'clinical effects',['cox_OS_mrd','cox_EFS_mrd','landmark']);vals=[];lo=[];hi=[];labs=[];cols=[]
    for end in ['OS','EFS']:
        v=data('cox_'+end+'_mrd').set_index('term').loc['ZEB1'];labs.append(end+' diagnosis model');vals.append(v['exp(coef)']);lo.append(v['exp(coef) lower 95%']);hi.append(v['exp(coef) upper 95%']);cols.append(C[end])
        v=min(js('landmark')['results'][end]['HR_by_elapsed_time'],key=lambda x:abs(x['elapsed_days']-365.25));labs.append(end+' day29 + 1 year');vals.append(v['HR']);lo.append(v['lower95']);hi.append(v['upper95']);cols.append(C[end])
    effect(a,labs,vals,lo,hi,'HR per log2(TPM+1)',cols,1);tag(a,'Distinct populations/time origins')
    a=ax(f,'d',[.66,.31,.27,.18],'apparent discrimination',['nested']);s=data('nested').set_index('pair').loc[['OS_clinical_mrd_plus_ZEB1','EFS_clinical_mrd_plus_ZEB1']];effect(a,['OS (1,302; 140 events)','EFS (1,302; 217 events)'],s.delta_C.tolist(),s.boot_dC_lo.tolist(),s.boot_dC_hi.tolist(),'Apparent ΔC: add ZEB1',[C['OS'],C['EFS']]);tag(a,'200-bootstrap intervals · no external OS test')
    a=ax(f,'e',[.21,.09,.20,.12],'conditional clinical block',['joint_maturity','bmp','ancestry']);s=t.clinical_block_table()
    for j,end in enumerate(['OS','EFS']):q=s[s.endpoint==end];a.scatter(q.P,np.arange(3)+(j-.5)*.16,color=C[end],s=22,label=end)
    a.set_yticks(range(3),['Joint normal scores','+ BMP17','+ BMP17 + ancestry']);a.set_xlabel('Joint ZEB1-block LRT P');a.axvline(.05,color='black',ls='--',lw=.5);a.set_xlim(-.01,max(.30,s.P.max()*1.12));f.legend(*a.get_legend_handles_labels(),frameon=False,ncol=2,loc='lower center',bbox_to_anchor=(.55,.002))
    a=ax(f,'f',[.75,.085,.18,.14],'independent MRD effects',['stjude_mrd']);s=data('stjude_mrd');effect(a,[v.outcome.replace('MRD','').replace('_',' ').strip()+' '+('base' if v.model=='base' else '+proxy') for _,v in s.iterrows()],s.ZEB1_beta_log10_MRD_per_log2_FPKM_plus1.tolist(),s.HC3_t_CI95_low.tolist(),s.HC3_t_CI95_high.tolist(),'St Jude log10 MRD slope');tag(a,'Different response endpoint')
    save(f,'Fig6')

def SFig1():
    f=start('SFig1',140,210);genes=sorted(t.normalized_pathway('hallmark_raw').pathway.str.replace('HALLMARK_','').tolist())
    for j,part in enumerate([genes[:25],genes[25:]]):
        a=ax(f,string.ascii_lowercase[j],[.24+j*.48,.12,.24,.79],'complete Hallmark family',['hallmark_raw','hallmark_residual'],'HERO','Fig1');pathway_compare(a,part)
        a.set_yticklabels([t.short_pathway('HALLMARK_'+g)[:25] for g in part]);a.set_xlabel('NES: + = ZEB1-high')
    f.legend([Line2D([0],[0],marker=m,color=C[k],ls='',mfc=fc) for k,m,fc in [('early','o',C['early']),('ZEB1','s',C['ZEB1']),('ZEB1','s','white')]],['Marginal (−source NES)','Subtype residual','Open = recorded q≥.05'],frameon=False,ncol=3,loc='lower center',bbox_to_anchor=(.5,.006));save(f,'SFig1')
def SFig2():
    f=start('SFig2');genes=['ZEB1','ZEB2','BCL11B','LYL1','CD34','CD1A','RAG1','RAG2','TCF7','GATA3','CD3D','CD8B']
    for letter,k,rect,col in [('a','normal195',[.11,.34,.35,.55],'facs_stage'),('b','normal206',[.66,.34,.26,.55],'inferred_stage')]:
        d=data(k).set_index(col);m=d[['g_'+g for g in genes]].T;m.index=genes;m=m.div(m.max(axis=1),axis=0);a=ax(f,letter,rect,'developmental marker matrix',[k],'HERO','Fig2');hm(f,a,m,'Gene-relative mean',0,1,'YlOrBr')
    d=data('normal_heldout')
    for l,x,k in [('c',.13,'cortical_vs_early'),('d',.64,'SP_vs_cortical')]:
        a=ax(f,l,[x,.10,.28,.10],'normal donor holdout',['normal_heldout'],'COMPACT_SUPPORT','Fig2');s=d[d.axis==k];a.scatter(s.heldout_difference,range(len(s)),color=C['cortical' if l=='c' else 'SP']);a.set_yticks(range(len(s)),s.heldout_donor.str.replace('TSC-',''));a.set_xlabel('Held-out paired score contrast');a.axvline(0,color='black',lw=.5);tag(a,k.replace('_',' '))
    save(f,'SFig2')
def legacy_supp(name):
    oldstart=t.start;t.start=lambda n,h=150:start(n,140,210)
    getattr(t,name)();t.start=oldstart
def SFig3():
    f=start('SFig3');d=data('sc_patient');genes=['ZEB1','ZEB2','BCL11B','MKI67','CD34','SPINK2','LYL1','CD1A','CD3D','CD8A','RAG1','TCF7']
    m=zcols(d.set_index('patient')[genes]).T
    a=ax(f,'a',[.11,.48,.44,.41],'malignant patient gene matrix',['sc_patient'],'HERO','Fig3');hm(f,a,m,'Across-patient gene z');a.set_xticklabels(m.columns,rotation=90);tag(a,'15 diagnosis patients · measured malignant pseudobulks')
    a=ax(f,'b',[.76,.60,.19,.28],'patient marker scatter',['sc_patient','sc_receipt'],'MAJOR_SUPPORT','Fig3');scatter(a,d,'BCL11B','ZEB1',C['ZEB1'],xlabel='BCL11B log2(CPM+1)',ylabel='ZEB1 log2(CPM+1)');tag(a,'Patient-level ρ=−0.268\nBH q=0.533')
    d=data('sc_state')
    for j,g in enumerate(['BCL11B','MKI67','SPINK2','CD1A']):
        a=ax(f,string.ascii_lowercase[j+2],[.10+j*.24,.10,.15,.20],'paired patient marker',['sc_state'],'COMPACT_SUPPORT','Fig3');t.paired(a,d,g)
    save(f,'SFig3')
def SFig4():
    f=start('SFig4');d=data('mutation_assoc').sort_values('gene')
    for j,s in enumerate([d.iloc[:19],d.iloc[19:]]):
        a=ax(f,string.ascii_lowercase[j],[.14+j*.47,.12,.30,.78],'complete lesion family',['mutation_assoc'],'HERO','Fig4');a.scatter(s.delta,range(len(s)),c=np.where(s.mw_q<.05,C['ZEB1'],'#999999'),s=23);a.set_yticks(range(len(s)),[f'{v.gene} ({int(v.n_mut)})' for _,v in s.iterrows()]);a.set_xlabel('Altered − absent median ZEB1');a.axvline(0,color='black',lw=.5);a.invert_yaxis()
    f.text(.17,.025,'Counts are altered patients. Orange: Mann–Whitney BH q<.05 across 37 fields; grey: q≥.05.',fontsize=7.5);save(f,'SFig4')
def SFig7():
    f=start('SFig7');d=data('cd4_hallmark');m=d.pivot(index='pathway',columns='condition',values='NES').sort_index();m.index=[t.short_pathway(x)[:25] for x in m.index]
    for j,s in enumerate([m.iloc[:25],m.iloc[25:]]):
        a=ax(f,string.ascii_lowercase[j],[.23+j*.48,.14,.23,.78],'normal CD4 perturbation landscape',['cd4_hallmark'],'HERO','Fig6');hm(f,a,s,'NES',-3,3);a.set_xticklabels(s.columns,rotation=40,ha='right')
    f.text(.10,.015,'Normal peripheral CD4 context: complete 50×3 pathway table. Inferential families retained in source data.',fontsize=7.5);save(f,'SFig7')
def SFig10():
    f=start('SFig10');d=data('lim_scores');r=data('lim_results');s=r[r.model=='primary_spearman'].sort_values(['cohort','axis'])
    a=ax(f,'a',[.18,.62,.27,.27],'patient bootstrap effects',['lim_results'],'HERO','Fig3');labs=[v.cohort[:1].upper()+': '+('C/E' if v.axis=='cortical_vs_early' else 'SP/C') for _,v in s.iterrows()];effect(a,labs,s.effect.tolist(),s.ci_low.tolist(),s.ci_high.tolist(),'Spearman ρ · patient bootstrap');a.set_xlim(-1,1)
    a=ax(f,'b',[.68,.62,.27,.27],'leave one patient out',['lim_results'],'MAJOR_SUPPORT','Fig3');s=r[r.model=='cycle_depth_rank_HC3'].sort_values(['cohort','axis']);effect(a,labs,s.effect.tolist(),s.loo_min.tolist(),s.loo_max.tolist(),'Adjusted slope: individual deletions')
    for l,x,y in [('c',.10,'zero_fraction_shared'),('d',.43,'n_cells'),('e',.76,'cycle_fraction')]:
        a=ax(f,l,[x,.13,.22,.23],'capture and eligibility',['lim_scores','lim_units'],'COMPACT_SUPPORT','Fig3')
        for k,col in [('discovery',C['early']),('validation',C['cortical'])]:s=d[d.cohort_analysis==k];a.scatter(s.ZEB1_log2_CPM_plus1,s[y],s=15,color=col,label=k,alpha=.7)
        a.set_xlabel('ZEB1 log2(CPM+1)');a.set_ylabel({'n_cells':'Malignant cells / patient','cycle_fraction':'Measured cycle fraction','zero_fraction_shared':'Shared genes: zero fraction'}[y]);a.legend(frameon=False)
    save(f,'SFig10')
def SFig11():
    f=start('SFig11',160,240);d=data('mutations').sort_values(['subtype','ZEB1','sample_id']);genes=data('mutation_assoc').gene.tolist();a=ax(f,'a',[.10,.15,.81,.73],'full genomic patient matrix',['mutations'],'HERO','Fig4');im=a.imshow(d[genes].T,aspect='auto',cmap=ListedColormap(['#F1F1F1',C['early']]),vmin=0,vmax=1,interpolation='nearest',rasterized=True);a.set_yticks(range(len(genes)),genes);a.set_xticks([]);a.set_xlabel('All 1,309 patients: subtype, then ZEB1; all 37 supplied alteration fields')
    offset=0
    for subtype,s in d.groupby('subtype',sort=False):a.axvline(offset-.5,color='black',lw=.35);offset+=len(s)
    f.legend([plt.Rectangle((0,0),1,1,color=C['early']),plt.Rectangle((0,0),1,1,color='#F1F1F1')],['Source field = altered','Source field = absent'],frameon=False,ncol=2,loc='lower center');save(f,'SFig11')
def SFig12():
    f=start('SFig12');d=data('depmap');a=ax(f,'a',[.14,.48,.30,.40],'native CRISPR origins',['depmap'],'HERO','Fig6');scatter(a,d,'ZEB1_log2_TPM_plus1','ZEB1_Chronos_effect',xlabel='Cell-line ZEB1 log2(TPM+1)',ylabel='Native Chronos effect');tag(a,'Qualified T-ALL origins · descriptive')
    a=ax(f,'b',[.66,.48,.27,.40],'native origin effect',['depmap'],'MAJOR_SUPPORT','Fig6');s=d[d.ZEB1_Chronos_effect.notna()];a.barh(np.arange(len(s)),s.ZEB1_Chronos_effect,color=C['ZEB1']);a.set_yticks(np.arange(len(s)),s.CellLineName);a.set_xlabel('Native ZEB1 Chronos');a.axvline(0,color='black',lw=.5);tag(a,'2 measured / 13 listed origins')
    a=ax(f,'c',[.21,.12,.62,.20],'normal CD4 pathway context',['cd4_hallmark'],'COMPACT_SUPPORT','Fig6');h=data('cd4_hallmark');genes=['HALLMARK_E2F_TARGETS','HALLMARK_G2M_CHECKPOINT','HALLMARK_MYC_TARGETS_V1']
    for j,(cond,s) in enumerate(h.groupby('condition')):
        s=s.set_index('pathway').reindex(genes);a.scatter(s.NES,np.arange(3)+(j-1)*.14,label=cond,s=25)
    a.set_yticks(range(3),['E2F','G2M','MYC V1']);a.set_xlabel('NES · normal CD4 ZEB1 perturbation');a.legend(frameon=False,ncol=3,loc='lower right');save(f,'SFig12')

FUNCTIONS={v.__name__:v for v in [Fig1,Fig2,Fig3,Fig4,Fig5,Fig6,SFig1,SFig2,SFig3,SFig4,SFig7,SFig10,SFig11,SFig12]}
for n in ['SFig5','SFig6','SFig8','SFig9']:FUNCTIONS[n]=lambda name=n:legacy_supp(name)
for n in sys.argv[1:] or [f'Fig{i}' for i in range(1,7)]+[f'SFig{i}' for i in range(1,13)]:FUNCTIONS[n]()
g=pd.DataFrame(t.PANELS);p=t.OUT/'PANEL_GRAMMAR.tsv'
if p.exists():old=pd.read_csv(p,sep='\t');g=pd.concat([old[~old.figure.isin(t.FIGS)],g],ignore_index=True)
g.to_csv(p,sep='\t',index=False)

from pathlib import Path
import os,json,string,hashlib,sys
import numpy as np,pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
from matplotlib.font_manager import fontManager
BASE=Path(__file__).resolve().parent;I=BASE/'inputs';OUT=BASE/'renders';OUT.mkdir(exist_ok=True)
for p in (BASE/'fonts').glob('*.ttf'):fontManager.addfont(str(p))
plt.rcParams.update({'font.family':'Arial','font.size':7.5,'axes.labelsize':8,'xtick.labelsize':7.5,'ytick.labelsize':7.5,'legend.fontsize':7.5,'axes.linewidth':.6,'lines.linewidth':1,'pdf.fonttype':42,'svg.fonttype':'none','text.color':'black','axes.labelcolor':'black','xtick.color':'black','ytick.color':'black','savefig.facecolor':'white'})
C={'early':'#5779A3','cortical':'#D99142','SP':'#478C79','OS':'#5A62A2','EFS':'#B56E78','ZEB1':'#C47435','ZEB2':'#597FAB'}
PANELS=[];FIGS=[];USED={};CURRENT=''
def data(k):
 p=I/(k+'.tsv');d=pd.read_csv(p,sep='\t');USED[k]=p;return d
def js(k):
 p=I/(k+'.json');USED[k]=p;return json.loads(p.read_text())
def start(name,h=150):
 global CURRENT
 CURRENT=name;fig=plt.figure(figsize=(180/25.4,h/25.4));return fig
def ax(fig,letter,rect,kind,sources,weight='M',parent=''):
 a=fig.add_axes(rect);a.text(-.13,1.045,letter,transform=a.transAxes,fontsize=11,fontweight='bold',ha='left',va='bottom');a.spines[['top','right']].set_visible(False)
 PANELS.append({'figure':CURRENT,'panel':letter,'status':'SUPP' if CURRENT.startswith('S') else 'MAIN','plot_type':kind,'relative_visual_weight':weight,'role':'HERO' if weight=='HERO' else ('MAJOR_SUPPORT' if weight in ['L','XL'] else 'COMPACT_SUPPORT'),'parent_main':parent,'source_aliases':';'.join(sources),'font_family':'Arial','min_font_pt':7.5,'semantic_palette':'ZEB1-development-v1','rect':json.dumps(rect)})
 return a
def stamp(a,s):a.text(.02,1.02,s,transform=a.transAxes,va='bottom',fontsize=7.5)
def scatter(a,d,x,y,color=C['ZEB1'],xlabel=None,ylabel=None):
 ok=d[[x,y]].dropna();a.scatter(ok[x],ok[y],s=8 if len(ok)>200 else 16,color=color,alpha=.55,edgecolors='none',rasterized=len(ok)>2000);a.set_xlabel(xlabel or x);a.set_ylabel(ylabel or y)
def boxes(a,d,group,value,order=None,colors=None,horizontal=False):
 order=order or list(d[group].dropna().unique());vals=[d.loc[d[group]==g,value].dropna().values for g in order];pos=np.arange(len(order))+1
 b=a.boxplot(vals,positions=pos,vert=not horizontal,patch_artist=True,widths=.6,showfliers=False,medianprops={'color':'black','linewidth':1},boxprops={'linewidth':.6},whiskerprops={'linewidth':.6},capprops={'linewidth':.6})
 for n,p in enumerate(b['boxes']):p.set_facecolor((colors or [C['ZEB1']]*len(order))[n]);p.set_alpha(.55)
 if horizontal:a.set_yticks(pos,order);a.set_xlabel(value)
 else:a.set_xticks(pos,order,rotation=35,ha='right');a.set_ylabel(value)
def matrix(fig,a,df,cmap='RdBu_r',vmin=None,vmax=None,key='Value'):
 cm=plt.get_cmap(cmap).copy();cm.set_bad('#DDDDDD');im=a.imshow(df.to_numpy(dtype=float),aspect='auto',interpolation='nearest',cmap=cm,vmin=vmin,vmax=vmax,rasterized=True);a.set_yticks(np.arange(len(df)),df.index);a.set_xticks(np.arange(len(df.columns)),df.columns,rotation=60,ha='right');cb=fig.colorbar(im,ax=a,fraction=.025,pad=.025);cb.set_label(key);return im
def short_pathway(s):
 return s.replace('HALLMARK_','').replace('EPITHELIAL_MESENCHYMAL_TRANSITION','EMT').replace('REACTIVE_OXYGEN_SPECIES_PATHWAY','ROS pathway').replace('WNT_BETA_CATENIN_SIGNALING','WNT / beta-catenin').replace('OXIDATIVE_PHOSPHORYLATION','Oxidative phosphorylation').replace('PI3K_AKT_MTOR_SIGNALING','PI3K / AKT / mTOR').replace('TNFA_SIGNALING_VIA_NFKB','TNF via NF-kB').replace('UNFOLDED_PROTEIN_RESPONSE','Unfolded protein response').replace('CHOLESTEROL_HOMEOSTASIS','Cholesterol homeostasis').replace('_',' ').capitalize()
def save(fig,name):
 folder=OUT/name;folder.mkdir(exist_ok=True);fig.savefig(folder/(name+'.pdf'));fig.savefig(folder/(name+'.svg'));fig.savefig(folder/(name+'.png'),dpi=300);plt.close(fig);FIGS.append(name);print(name+' RENDERED',flush=True)
def cohort(k):
 d=data(k);yc=[c for c in d if c.startswith('ZEB1_log2')][0];return d,yc
def normalized_pathway(k):
 d=data(k);return d.rename(columns={'Term':'pathway','FDR q-val':'q','NOM p-val':'p','author_method_within_condition_padj':'q','BH_q_full150':'q_full'})
def effect(a,labels,vals,lows,highs,key='Effect',colors=None,null=0):
 y=np.arange(len(labels));
 for j in range(len(y)):a.plot([lows[j],highs[j]],[j,j],color=(colors or [C['ZEB1']]*len(y))[j]);a.scatter(vals[j],j,color=(colors or [C['ZEB1']]*len(y))[j],s=20)
 a.axvline(null,color='black',lw=.6,ls='--');a.set_yticks(y,labels);a.invert_yaxis();a.set_xlabel(key)
def paired(a,d,gene):
 p=d.pivot(index='patient',columns='mapped_state',values=gene).reindex(columns=['Early_DN','DP']).dropna()
 for patient,r in p.iterrows():a.plot([0,1],r.values,color='#777777',alpha=.6);a.scatter([0,1],r.values,c=[C['early'],C['cortical']],s=15)
 a.set_xticks([0,1],['Early DN','DP']);a.set_ylabel(gene+' log2(CPM+1)');return len(p)
def umap(a,d,feature=False):
 if feature:
  z=np.log1p(d.ZEB1_count);idx=np.argsort(z);im=a.scatter(d.UMAP1.iloc[idx],d.UMAP2.iloc[idx],c=z.iloc[idx],s=.12,cmap='YlOrBr',rasterized=True,linewidths=0);plt.colorbar(im,ax=a,fraction=.03,pad=.02,label='ZEB1 log1p(raw UMI)')
 else:
  labels=data('normal_labels').set_index('author_cell_type').display_label.to_dict();groups=list(d.author_cell_type.unique());palette=plt.get_cmap('tab20')(np.arange(len(groups))%20)
  for n,g in enumerate(groups):
   t=d[d.author_cell_type==g];a.scatter(t.UMAP1,t.UMAP2,s=.1,color=palette[n],rasterized=True,linewidths=0)
  # Every frozen category is visibly keyed in a separate complete atlas legend.
  from matplotlib.lines import Line2D
  handles=[Line2D([0],[0],marker='o',color='white',markerfacecolor=palette[n],markersize=3,label=labels.get(g,g)) for n,g in enumerate(groups)]
  a.legend(handles=handles,loc='upper left',bbox_to_anchor=(1.015,1),ncol=1,frameon=False,handletextpad=.2,labelspacing=.12)
 a.set_xlabel('UMAP 1');a.set_ylabel('UMAP 2');a.set_aspect('equal',adjustable='datalim');a.set_xticks([]);a.set_yticks([])
def Fig1():
 f=start('Fig1',160);d=data('bulk').merge(data('maturity'),on='sample_id',validate='one_to_one')
 a=ax(f,'a',[.29,.55,.64,.39],'patient distribution',['bulk'],'HERO');order=d.groupby('subtype').ZEB1.median().sort_values().index.tolist();boxes(a,d,'subtype','ZEB1',order,horizontal=False);a.set_ylabel('ZEB1 log2(TPM+1)');a.tick_params(axis='x',rotation=65);stamp(a,'AALL0434 · 1,309 patients')
 a=ax(f,'b',[.08,.58,.14,.33],'box distribution',['bulk']);boxes(a,d,'etp3','ZEB1',['ETP','Near-ETP','Non-ETP'],[C['early'],C['cortical'],C['SP']]);a.set_ylabel('ZEB1 log2(TPM+1)')
 a=ax(f,'c',[.10,.15,.23,.26],'scatter',['bulk','maturity'],'L');scatter(a,d,'cortical_vs_early','ZEB1',xlabel='Cortical vs early score (z)',ylabel='ZEB1 log2(TPM+1)')
 a=ax(f,'d',[.42,.15,.23,.26],'scatter',['bulk','maturity'],'L');scatter(a,d,'SP_vs_cortical','ZEB1',C['SP'],xlabel='SP vs cortical score (z)',ylabel='ZEB1 log2(TPM+1)')
 a=ax(f,'e',[.76,.15,.20,.26],'gsea pathway',['hallmark_raw','hallmark_residual']);raw=normalized_pathway('hallmark_raw');res=normalized_pathway('hallmark_residual');genes=['HALLMARK_E2F_TARGETS','HALLMARK_G2M_CHECKPOINT','HALLMARK_MYC_TARGETS_V1','HALLMARK_PROTEIN_SECRETION'];
 for j,g in enumerate(genes):
  for n,(t,col) in enumerate([(raw,C['early']),(res,C['ZEB1'])]):
   r=t[t.pathway==g].iloc[0];display_nes=-r.NES if n==0 else r.NES;a.scatter(display_nes,j+(n-.5)*.2,s=18,color=col,marker='o' if float(r.q)<.05 else 'x')
 a.set_yticks(range(4),['E2F','G2M','MYC V1','Secretion']);a.set_xlabel('NES (+ = ZEB1-high)');a.axvline(0,color='black',lw=.5);f.legend([plt.Line2D([0],[0],marker='o',color=c,ls='') for c in [C['early'],C['ZEB1']]]+[plt.Line2D([0],[0],marker=m,color='black',ls='') for m in ['o','x']],['Marginal (−original NES)','Subtype residual','q < 0.05','q ≥ 0.05'],frameon=False,loc='lower right',bbox_to_anchor=(.98,.005),ncol=2);save(f,'Fig1')
def Fig2():
 f=start('Fig2',165);d=data('normal_umap');a=ax(f,'a',[.07,.53,.40,.42],'UMAP embedding',['normal_umap','normal_labels'],'HERO');umap(a,d);stamp(a,'Frozen normal-thymus atlas')
 a=ax(f,'b',[.72,.53,.23,.40],'FeaturePlot embedding',['normal_umap'],'L');umap(a,d,True)
 p=data('normal195').set_index('facs_stage').reindex(['DN1','DN2','DN3','ISP','DP_CD3neg','DP_CD3pos','CD4SP','CD8SP'])
 a=ax(f,'c',[.08,.09,.32,.29],'stage trajectory descriptive',['normal195'],'L');a.plot(range(len(p)),p.g_ZEB1,'o-',color=C['ZEB1'],label='ZEB1');a.plot(range(len(p)),p.g_ZEB2,'o-',color=C['ZEB2'],label='ZEB2');a.set_xticks(range(len(p)),p.index,rotation=45,ha='right');a.set_ylabel('Mean frozen normalized expression');a.legend(frameon=False);stamp(a,'GSE195812 · FACS stages')
 p=data('normal206').set_index('inferred_stage');a=ax(f,'d',[.49,.10,.19,.28],'stage trajectory descriptive',['normal206']);a.plot(range(len(p)),p.g_ZEB1,'o-',color=C['ZEB1']);a.set_xticks(range(len(p)),['DN early','DP immature','DP mature'],rotation=40,ha='right');a.set_ylabel('Mean ZEB1');stamp(a,'GSE206710')
 t=data('normal_heldout');a=ax(f,'e',[.79,.12,.17,.26],'paired donor distribution',['normal_heldout']);boxes(a,t,'axis','heldout_difference',['cortical_vs_early','SP_vs_cortical'],[C['cortical'],C['SP']]);a.set_xticklabels(['Cortical / early','SP / cortical'],rotation=40,ha='right');a.set_ylabel('Held-out donor score contrast');a.axhline(0,color='black',lw=.5);save(f,'Fig2')
def Fig3():
 f=start('Fig3',150);d=data('sc_state');a=ax(f,'a',[.09,.55,.30,.37],'paired patient',['sc_state'],'HERO');n=paired(a,d,'ZEB1');stamp(a,f'GSE248287 · {n} paired patients\nExact Wilcoxon P = 0.0625')
 p=data('sc_patient');a=ax(f,'b',[.50,.57,.44,.35],'patient gene heatmap',['sc_patient'],'XL');g=['ZEB1','BCL11B','MKI67','SPINK2','CD34','CD1A'];m=p.set_index('patient')[g];m=m.loc[m.ZEB1.sort_values().index];m=(m-m.mean())/m.std(ddof=0);matrix(f,a,m.T,key='Within-gene display z');a.set_xlabel('Patients · diagnosis malignant pseudobulk')
 for letter,rect,g,col in [('c',[.08,.10,.23,.29],'BCL11B',C['cortical']),('d',[.42,.10,.23,.29],'SPINK2',C['early']),('e',[.76,.10,.20,.29],'MKI67',C['SP'])]:
  a=ax(f,letter,rect,'patient scatter',['sc_patient','sc_receipt']);scatter(a,p,'ZEB1',g,col,xlabel='ZEB1 log2(CPM+1)',ylabel=g+' log2(CPM+1)');r=js('sc_receipt')['patient_level_spearman_vs_ZEB1'][g];stamp(a,f"n = 15 · rho = {r['rho']:.2f}\nBH q = {r['q']:.3f}")
 save(f,'Fig3')
def onco(a,f,d,genes):
 d=d.sort_values(['ZEB1','sample_id']);m=d[genes].T;im=a.imshow(m,aspect='auto',interpolation='nearest',cmap=ListedColormap(['#F2F2F2',C['early']]),vmin=0,vmax=1,rasterized=True);a.set_yticks(range(len(genes)),genes);a.set_xticks([0,len(d)-1],['Low ZEB1','High ZEB1']);a.set_xlabel(f'{len(d):,} patients ordered by measured ZEB1');a.legend([plt.Rectangle((0,0),1,1,color=C['early']),plt.Rectangle((0,0),1,1,color='#F2F2F2')],['Altered','Not altered'],frameon=False,loc='upper right',bbox_to_anchor=(1.02,1.16),ncol=2)
def Fig4():
 f=start('Fig4',155);d=data('mutations');genes=['CDKN2A','NOTCH1','TLX1','TLX3','DNM2','PHF6','IL7R','BCL11B','FLT3','NRAS','JAK3','MED12'];a=ax(f,'a',[.12,.59,.82,.34],'Oncoprint genomic',['mutations'],'HERO');onco(a,f,d,genes)
 for letter,rect,g in [('b',[.08,.12,.20,.29],'NOTCH1'),('c',[.39,.12,.20,.29],'CDKN2A'),('d',[.70,.12,.20,.29],'MED12')]:
  a=ax(f,letter,rect,'box distribution',['mutations','mutation_assoc'],'L');boxes(a,d,g,'ZEB1',[0,1],[C['SP'],C['ZEB1']]);a.set_xticklabels(['Unaltered','Altered'],rotation=30,ha='right');a.set_ylabel('ZEB1 log2(TPM+1)');r=data('mutation_assoc').set_index('gene').loc[g];stamp(a,f'{g}\nMW BH q = {r.mw_q:.2g}')
 save(f,'Fig4')
def Fig5():
 f=start('Fig5',170)
 for row,(k,label) in enumerate([('nopho','NOPHO108'),('aiims','AIIMS35'),('india','India25')]):
  d,y=cohort(k);r=js(k+'_results')['results']
  for col,(axis,clr) in enumerate([('cortical_vs_early',C['cortical']),('SP_vs_cortical',C['SP'])]):
   yl={'nopho':'ZEB1 log2(CPM+1)','aiims':'ZEB1 log2(FPKM+1)','india':'ZEB1 log2(counts+1)'}[k];a=ax(f,string.ascii_lowercase[row*2+col],[.12+col*.49,.72-row*.30,.34,.20],'scatter',[k,k+'_results'],'HERO' if row==0 else 'L');scatter(a,d,axis+'_raw',y,clr,xlabel=('Cortical vs early' if col==0 else 'SP vs cortical')+' fixed rank score',ylabel=yl);v=next(t for t in r if t['axis']==axis);stamp(a,f"{label} · n = {len(d)}\nrho = {v['spearman_rho']:.2f} · BH2 q = {v['BH2_q']:.2g}")
 save(f,'Fig5')
def clinical_block_table():
 rows=[]
 for end in ['OS','EFS']:
  j=js('joint_maturity');v=j['results'][end]['ZEB1_block_LRT'];q=next(x['BH_q_four'] for x in j['BH4_FAMILY'] if x['model']=='joint_sensitivity' and x['endpoint']==end);rows.append({'endpoint':end,'context':'Joint maturity','P':v['p'],'q':q,'family':'BH4'})
  v=js('bmp')['results'][end]['primary_post_ZEB1_block'];rows.append({'endpoint':end,'context':'+ BMP17','P':v['p'],'q':v['BH2_q'],'family':'BH2'})
  v=js('ancestry')['results'][end]['with_BMP17'];z=next(t for key,t in v.items() if isinstance(t,dict) and 'p' in t and ('block' in key.lower() or 'target' in key.lower()));q=next(t for key,t in z.items() if 'q' in key.lower());rows.append({'endpoint':end,'context':'+ BMP17 + ancestry','P':z['p'],'q':q,'family':'BH4'})
 return pd.DataFrame(rows)
def Fig6():
 f=start('Fig6',180)
 for letter,end,rect in [('a','OS',[.09,.70,.36,.23]),('b','EFS',[.59,.70,.36,.23])]:
  a=ax(f,letter,rect,'survival patient distribution',['km_'+end.lower(),'km_counts_'+end.lower()],'HERO');d=data('km_'+end.lower());
  for k,c in zip(['T1_low','T2_mid','T3_high'],[C['early'],'#999999',C['ZEB1']]):
   r=d[d.tertile==k];a.step(r.time/365.25,r.surv,where='post',color=c,label=k.replace('_',' '))
  a.set_ylim(.5,1.01);a.set_xlabel('Years from diagnosis');a.set_ylabel(end+' probability');a.legend(frameon=False,ncol=3,loc='lower left')
 a=ax(f,'c',[.19,.40,.26,.20],'forest effect',[f'cox_{e}_{m}' for e in ['OS','EFS'] for m in ['uni','etp','subtype','mrd']]);vals=[];lo=[];hi=[];labs=[];cols=[]
 for end in ['OS','EFS']:
  for m,label in [('uni','Unadjusted'),('etp','Clinical + ETP'),('subtype','+ subtype'),('mrd','+ MRD')]:
   r=data(f'cox_{end}_{m}').set_index('term').loc['ZEB1'];labs.append(end+' '+label);vals.append(r['exp(coef)']);lo.append(r['exp(coef) lower 95%']);hi.append(r['exp(coef) upper 95%']);cols.append(C[end])
 effect(a,labs,vals,lo,hi,'HR per log2(TPM+1)',cols,1)
 a=ax(f,'d',[.68,.40,.25,.20],'paired concordance scatter',['nested']);d=data('nested');s=d[d.pair.isin(['OS_clinical_mrd_plus_ZEB1','EFS_clinical_mrd_plus_ZEB1'])];
 for n,(_,r) in enumerate(s.iterrows()):a.plot([0,1],[r.C_reduced,r.C_full],color=C['OS' if n==0 else 'EFS']);a.scatter([0,1],[r.C_reduced,r.C_full],color=C['OS' if n==0 else 'EFS'],s=20);a.text(1.06,r.C_full,f"{'OS' if n==0 else 'EFS'}\nΔC {r.delta_C:.4f}",va='center')
 a.set_xticks([0,1],['Clinical + MRD','+ ZEB1']);a.set_xlim(-.2,1.7);a.set_ylabel('Apparent Harrell C')
 a=ax(f,'e',[.25,.08,.25,.20],'dotplot conditional tests',['joint_maturity','bmp','ancestry']);d=clinical_block_table();labels=['Joint maturity','+ BMP17','+ BMP17 + ancestry'];
 for n,end in enumerate(['OS','EFS']):
  t=d[d.endpoint==end];a.scatter(t.P,np.arange(3)+(n-.5)*.15,color=C[end],s=20,label=end)
 a.set_yticks(range(3),labels);a.set_xlabel('Joint ZEB1-block LRT P');a.axvline(.05,color='black',ls='--',lw=.6);a.legend(frameon=False,loc='lower right');a.set_xlim(-.01,.3)
 a=ax(f,'f',[.72,.08,.22,.20],'forest endpoint effect',['stjude_mrd']);t=data('stjude_mrd');effect(a,[r.outcome.replace('_',' ')+'\n'+('Base' if r.model=='base' else '+ proxy') for _,r in t.iterrows()],t.ZEB1_beta_log10_MRD_per_log2_FPKM_plus1.tolist(),t.HC3_t_CI95_low.tolist(),t.HC3_t_CI95_high.tolist(),'ZEB1 slope · log10 MRD');stamp(a,'Independent St Jude\nDifferent treatment-response endpoint');save(f,'Fig6')
def SFig1():
 f=start('SFig1',200)
 for col,k in enumerate(['hallmark_raw','hallmark_residual']):
  d=normalized_pathway(k);d=d.sort_values('pathway');a=ax(f,string.ascii_lowercase[col],[.28+col*.49,.15,.19,.76],'gsea pathway',[k],'HERO',parent='Fig1');a.scatter(d.NES,range(len(d)),c=np.where(d.q<.05,C['ZEB1'],'#B9B9B9'),s=15);a.set_yticks(range(len(d)),[short_pathway(x) for x in d.pathway]);a.set_xlabel('NES (+ = ZEB1-low)' if col==0 else 'NES (+ = ZEB1-high)');a.axvline(0,color='black',lw=.5);stamp(a,'Marginal' if col==0 else 'Subtype residual')
 f.legend([plt.Line2D([0],[0],marker='o',ls='',color=c) for c in [C['ZEB1'],'#B9B9B9']],['q < 0.05','q ≥ 0.05'],frameon=False,ncol=2,loc='lower center',bbox_to_anchor=(.5,.04))
 save(f,'SFig1')
def SFig2():
 f=start('SFig2',155);d=data('normal195').set_index('facs_stage').reindex(['DN1','DN2','DN3','ISP','DP_CD3neg','DP_CD3pos','CD4SP','CD8SP']);genes=['ZEB1','ZEB2','BCL11B','LYL1','CD34','CD1A','RAG1','RAG2','TCF7','GATA3','CD3D','CD8B'];m=d[['g_'+g for g in genes]].T;m.index=genes;m=m.div(m.max(axis=1),axis=0);a=ax(f,'a',[.12,.51,.43,.42],'gene heatmap',['normal195'],'HERO',parent='Fig2');matrix(f,a,m,cmap='YlOrBr',vmin=0,vmax=1,key='Gene-relative mean')
 d=data('normal206').set_index('inferred_stage');m=d[['g_'+g for g in genes]].T;m.index=genes;m=m.div(m.max(axis=1),axis=0);a=ax(f,'b',[.71,.51,.22,.42],'gene heatmap',['normal206'],'L',parent='Fig2');matrix(f,a,m,cmap='YlOrBr',vmin=0,vmax=1,key='Gene-relative mean');a.set_xticklabels(['DN early','DP imm.','DP mat.'],rotation=45,ha='right')
 d=data('normal_state_scores');p=d.pseudobulk_id.str.split('|',expand=True);d['donor']=p[0];d['stage']=p[1]
 for n,k in enumerate(['cortical_vs_early','SP_vs_cortical']):
  a=ax(f,string.ascii_lowercase[n+2],[.12+n*.47,.09,.31,.25],'paired donor',['normal_state_scores'],'M',parent='Fig2');order=['Early_DN','DP','CD4SP','CD8SP'];order=[x for x in order if x in d.stage.unique()];
  if not order:order=list(d.stage.unique())
  for donor,t in d.groupby('donor'):v=t.groupby('stage')[k].mean().reindex(order);a.plot(range(len(order)),v,'o-',alpha=.65,label=donor)
  a.set_xticks(range(len(order)),order,rotation=35,ha='right');a.set_ylabel(k.replace('_',' ')+' score');
  if n==1:
   handles,labels=a.get_legend_handles_labels();f.legend(handles,labels,frameon=False,ncol=5,loc='center',bbox_to_anchor=(.53,.39));save(f,'SFig2')
def SFig3():
 f=start('SFig3',150);d=data('sc_state')
 for n,g in enumerate(['BCL11B','MKI67','SPINK2','CD1A']):
  a=ax(f,string.ascii_lowercase[n],[.12+(n%2)*.49,.58-(n//2)*.45,.30,.32],'paired patient',['sc_state'],'HERO' if n==0 else 'M',parent='Fig3');paired(a,d,g)
 save(f,'SFig3')
def SFig4():
 f=start('SFig4',160);d=data('mutation_assoc').sort_values('delta');a=ax(f,'a',[.14,.12,.38,.79],'lesion dotplot',['mutation_assoc'],'HERO',parent='Fig4');a.scatter(d.delta,range(len(d)),c=np.where(d.mw_q<.05,C['ZEB1'],'#999999'),s=20);a.set_yticks(range(len(d)),d.gene);a.axvline(0,color='black',lw=.5);a.set_xlabel('Altered − unaltered median ZEB1')
 a=ax(f,'b',[.69,.12,.22,.79],'bar composition',['mutation_assoc'],'L',parent='Fig4');a.barh(range(len(d)),d.n_mut,color=C['early']);a.set_yticks(range(len(d)),d.gene);a.set_xlabel('Altered patients');f.legend([plt.Line2D([0],[0],marker='o',ls='',color=c) for c in [C['ZEB1'],'#999999']],['MW BH q < 0.05','MW BH q ≥ 0.05'],frameon=False,ncol=2,loc='lower center',bbox_to_anchor=(.5,.01));save(f,'SFig4')
def SFig5():
 f=start('SFig5',165);d=data('aiims_native')
 for n,k in enumerate(['cortical_vs_early_raw','SP_vs_cortical_raw']):
  a=ax(f,string.ascii_lowercase[n],[.15+n*.48,.58,.31,.34],'box distribution',['aiims_native'],'HERO' if n==0 else 'L',parent='Fig5');boxes(a,d,'source_disease_state',k,colors=[C['early'],C['cortical'],C['SP']]);a.set_xticklabels([x.replace('_TALL','') for x in d.source_disease_state.unique()],rotation=40,ha='right');a.set_ylabel(k.replace('_raw','').replace('_',' ')+' score')
 t=data('aiims_conditional');a=ax(f,'c',[.23,.10,.26,.29],'forest conditional effect',['aiims_conditional'],'M',parent='Fig5');effect(a,['SP / cortical','Cortical / early'],t.conditional_rank_rho.tolist(),t.stratified_patient_bootstrap_CI95_low.tolist(),t.stratified_patient_bootstrap_CI95_high.tolist(),'Within-phenotype rank rho');stamp(a,'AIIMS35 · BH2 q = 0.023 / 0.760')
 d=data('nopho_genes');a=ax(f,'d',[.68,.10,.26,.29],'scatter',['nopho_genes','nopho_gene_assoc'],'M',parent='Fig5');scatter(a,d,'ZEB1','BCL11B',xlabel='NOPHO ZEB1 (author scale)',ylabel='BCL11B (author scale)');stamp(a,'n = 108 · rho = 0.27\nBH q = 0.0078');save(f,'SFig5')
def SFig6():
 f=start('SFig6',170)
 for n,end in enumerate(['OS','EFS']):
  d=js('landmark')['results'][end];t=pd.DataFrame(d['HR_by_elapsed_time']);a=ax(f,string.ascii_lowercase[n],[.13+n*.49,.58,.30,.33],'time effect',['landmark'],'HERO' if n==0 else 'L',parent='Fig6');x=t.elapsed_days/365.25;a.plot(x,t.HR,'o-',color=C[end]);a.fill_between(x,t.lower95,t.upper95,color=C[end],alpha=.2);a.axhline(1,color='black',lw=.6,ls='--');a.set_xlabel('Years after day29 landmark');a.set_ylabel(end+' pointwise HR');stamp(a,f"Joint BH2 q = {d['ZEB1_block_LRT']['BH_q_two_endpoints']:.3g}\nTime-interaction P = {d['ZEB1_time_increment_LRT']['p']:.3g}")
 d=data('nested');a=ax(f,'c',[.26,.12,.25,.27],'forest incremental effect',['nested'],'M',parent='Fig6');s=d[d.boot_dC_lo.notna()];labels={'OS_clinical_mrd_plus_ZEB1':'OS · add ZEB1','EFS_clinical_mrd_plus_ZEB1':'EFS · add ZEB1','OS_mrd_increment_noZEB1':'OS · add MRD','OS_clinical_plus_MRD':'OS · add MRD','EFS_clinical_plus_MRD':'EFS · add MRD','NonETP_OS_mrd_plus_ZEB1':'Non-ETP OS · add ZEB1'};effect(a,[labels.get(x,x.replace('clinical','').replace('_',' ')) for x in s.pair],s.delta_C.tolist(),s.boot_dC_lo.tolist(),s.boot_dC_hi.tolist(),'Bootstrap Δ Harrell C');
 d=clinical_block_table();a=ax(f,'d',[.77,.12,.18,.27],'conditional dotplot',['joint_maturity','bmp','ancestry'],'M',parent='Fig6');a.scatter(d.q,range(len(d)),c=[C[x] for x in d.endpoint],s=20);a.set_yticks(range(len(d)),[r.endpoint+' '+r.context for _,r in d.iterrows()]);a.axvline(.05,ls='--',color='black',lw=.6);a.set_xlabel('Recorded family-specific BH q');save(f,'SFig6')
def SFig7():
 f=start('SFig7',180);d=data('cd4_hallmark');m=d.pivot(index='pathway',columns='condition',values='NES');m.index=[short_pathway(x) for x in m.index];a=ax(f,'a',[.30,.12,.28,.79],'perturbation pathway heatmap',['cd4_hallmark'],'HERO',parent='Fig6');matrix(f,a,m,vmin=-3,vmax=3,key='ZEB1 perturbation NES')
 d=data('depmap');a=ax(f,'b',[.75,.58,.20,.30],'origin scatter',['depmap'],'L',parent='Fig6');scatter(a,d,'ZEB1_log2_TPM_plus1','ZEB1_Chronos_effect',xlabel='Cell-line ZEB1 log2(TPM+1)',ylabel='Native ZEB1 Chronos');stamp(a,'T-ALL origins · descriptive')
 a=ax(f,'c',[.75,.11,.20,.30],'perturbation pathway dotplot',['cd4_hallmark'],'M',parent='Fig6');s=d=None;h=data('cd4_hallmark');s=h[h.pathway.isin(['HALLMARK_E2F_TARGETS','HALLMARK_MYC_TARGETS_V1','HALLMARK_G2M_CHECKPOINT'])];
 for n,(cond,t) in enumerate(s.groupby('condition')):a.scatter(t.NES,np.arange(len(t))+(n-1)*.15,label=cond,s=18)
 a.set_yticks(range(3),['E2F','G2M','MYC V1']);a.set_xlabel('NES · normal CD4 context');handles,labels=a.get_legend_handles_labels();f.legend(handles,labels,frameon=False,loc='lower right',bbox_to_anchor=(.99,.008),ncol=3);save(f,'SFig7')
def SFig8():
 f=start('SFig8',170);d=data('drug_slots');drugs=list(d.drug.unique())
 for n,end in enumerate(['native_normalized','fractional_midrank']):
  t=d[d.endpoint==end];m=t.pivot(index='drug',columns='model',values='beta').reindex(drugs);a=ax(f,string.ascii_lowercase[n],[.22+n*.47,.53,.20,.39],'patient drug heatmap',['drug_slots'],'HERO' if n==0 else 'L',parent='Fig6');matrix(f,a,m,vmin=-.3,vmax=.3,key='ZEB1 slope');a.set_xticklabels(['Base','+ fixed proxy'],rotation=40,ha='right');stamp(a,'Normalized LC50' if n==0 else 'Fractional midrank LC50')
 a=ax(f,'c',[.18,.16,.30,.19],'drug availability bar',['drug_slots'],'M',parent='Fig6');s=d.drop_duplicates('drug');a.bar(range(len(s)),s.n,color=C['early']);a.set_xticks(range(len(s)),s.drug,rotation=75,ha='right');a.set_ylabel('Complete patients');a.axhline(30,color='black',lw=.5,ls='--');f.legend([plt.Rectangle((0,0),1,1,color='#DDDDDD')],['UNESTIMATED (40 / 72 slots)'],frameon=False,loc='lower right',bbox_to_anchor=(.55,.39))
 d=data('main_mrd');a=ax(f,'d',[.77,.10,.18,.25],'forest response effect',['main_mrd'],'M',parent='Fig6');s=d[d.endpoint.str.contains('continuous',case=False)]
 if len(s)!=2:s=d.iloc[:2]
 effect(a,['Fixed proxy','+ BMP17 / ancestry'],s.target_beta.tolist(),s.CI95_low.tolist(),s.CI95_high.tolist(),'Main-cohort MRD slope');stamp(a,'Same AALL0434 patients');save(f,'SFig8')
def SFig9():
 f=start('SFig9',165);d=data('old_pair');a=ax(f,'a',[.11,.56,.30,.36],'paired patient',['old_pair'],'HERO',parent='Fig1');
 for _,r in d.iterrows():a.plot([0,1],[r.ZEB1_paired_delta,r.ZEB2_paired_delta],color='#888888',alpha=.5)
 a.scatter(np.zeros(len(d)),d.ZEB1_paired_delta,color=C['ZEB1'],s=15);a.scatter(np.ones(len(d)),d.ZEB2_paired_delta,color=C['ZEB2'],s=15);a.set_xticks([0,1],['ZEB1 programme','ZEB2 programme']);a.set_ylabel('Fixed programme paired contrast');stamp(a,'18 matched pairs · accepted original result')
 d=data('old_definitions');a=ax(f,'b',[.65,.56,.29,.36],'forest definition effect',['old_definitions'],'M',parent='Fig1');effect(a,d.definition.tolist(),d.estimate.tolist(),d.ci_low.tolist(),d.ci_high.tolist(),'Bilateral programme asymmetry')
 d=data('old_samecell');d=d[(d.method=='log_primary')&(d.status=='ESTIMATED')];a=ax(f,'c',[.11,.10,.30,.30],'patient distribution',['old_samecell'],'L',parent='Fig3');boxes(a,d,'cohort','r');a.set_ylabel('Fixed output within-cell partial r');a.axhline(0,color='black',lw=.5)
 d=data('old_fitness');a=ax(f,'d',[.65,.10,.29,.30],'origin distribution',['old_fitness'],'M',parent='Fig6');boxes(a,d,'analysis_lineage','relative_ZEB1_minus_ZEB2');a.set_ylabel('Native ZEB1 − ZEB2 fitness');stamp(a,'Exploratory distinct origins\nNot independent-release replication');save(f,'SFig9')

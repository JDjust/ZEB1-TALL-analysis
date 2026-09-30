from plot_common import *
def main():
    f=fig(218);atlas=read(V/'Yayon_existing_umap_ZEB1_ZEB2.tsv.gz');assert len(atlas)==391462
    def grp(s):
        return 'Early/DN' if 'ETP' in s or 'DN' in s else 'DP' if s.startswith('T_DP') else 'SP' if s in ['T_CD4','T_CD8'] else 'Other'
    group=atlas.author_cell_type.map(grp);pal={'Early/DN':Z2,'DP':Z1,'SP':GREY,'Other':'#D9D9D9'}
    for j,(letter,title,gene) in enumerate([('A','Normal thymocyte states',None),('B','ZEB1','ZEB1_count'),('C','ZEB2','ZEB2_count')]):
        a=ax(f,box(f,13+j*56,8,49,32),letter,title)
        if gene is None:
            for name in ['Other','SP','DP','Early/DN']:
                q=atlas[group==name];a.scatter(q.UMAP1,q.UMAP2,s=.04,color=pal[name],alpha=.4,rasterized=True,lw=0)
            for name in ['T_ETP','T_DN(Q)','T_DP(Q)','T_CD4','T_CD8']:
                q=atlas[atlas.author_cell_type==name];a.annotate(name.removeprefix('T_'),(q.UMAP1.median(),q.UMAP2.median()),xytext={'T_CD4':(-8,8),'T_CD8':(10,-7),'T_ETP':(12,-8)}.get(name,(0,0)),textcoords='offset points',ha='center',fontsize=8,bbox={'facecolor':'white','edgecolor':'none','alpha':.7,'pad':.1},arrowprops={'arrowstyle':'-','lw':.4,'color':INK} if name=='T_ETP' else None)
        else:
            a.scatter(atlas.UMAP1,atlas.UMAP2,s=.025,color='#D9D9D9',rasterized=True,lw=0)
            q=atlas[atlas[gene]>0];cmap=LinearSegmentedColormap.from_list(gene,[NA,Z1 if j==1 else Z2]);im=a.scatter(q.UMAP1,q.UMAP2,c=np.log1p(q[gene]),vmin=0,vmax=np.log1p(5),cmap=cmap,s=.15,rasterized=True,lw=0)
            cb=scale(f,im,box(f,75+(j-1)*56,43,34,1.2),('ZEB1' if j==1 else 'ZEB2')+' log1p count',[0,.9,np.log1p(5)]);cb.ax.set_xticklabels(['0','0.9','1.8'])
        a.set(xlim=(-3.75,13.9),ylim=(-3,9.35),xticks=[],yticks=[]);a.spines[:].set_visible(False)
    m=ax(f,box(f,25,63,144,26),'D','Six-gene state context')
    dot=read(V/'F1D_marker_dot_summary.tsv');st=['T_ETP','T_DN(early)','T_DN(P)','T_DN(Q)','T_DP(Q)-early','T_DP(Q)','T_DP(P)','T_CD4','T_CD8'];genes=['ZEB1','ZEB2','CD34','LYL1','CD1A','BCL11B']
    mean=dot.pivot(index='gene',columns='stage',values='mean_log1p').loc[genes,st];frac=dot.pivot(index='gene',columns='stage',values='fraction_detected').loc[genes,st];z=mean.sub(mean.mean(axis=1),axis=0).div(mean.std(axis=1),axis=0).clip(-2,2)
    for i,g in enumerate(genes):
        for j,s in enumerate(st):im=m.scatter(j,i,s=15+90*frac.loc[g,s],c=[z.loc[g,s]],vmin=-2,vmax=2,cmap=CMAP,lw=.25,edgecolor='white')
    m.set(ylim=(5.6,-1.8),xlim=(-.5,8.5));m.set_yticks(range(6),genes);m.set_xticks(range(9),['ETP','DN early','DN(P)','DN(Q)','DP(Q) early','DP(Q)','DP(P)','CD4','CD8'],rotation=25,ha='right');m.tick_params(length=0);m.spines[:].set_visible(False)
    for start,end,label,co in [(0,3,'Early / DN',Z2),(4,6,'DP',Z1),(7,8,'SP',GREY)]:
        m.add_patch(Rectangle((start-.45,-1.65),end-start+.9,.65,color=co,alpha=.25,lw=0));m.text((start+end)/2,-1.32,label,ha='center',va='center')
    for v,x in zip([.25,.5,.75],[.48,.59,.70]):
        m.scatter(x,1.21,s=15+90*v,transform=m.transAxes,facecolor='white',edgecolor=INK,lw=.5,clip_on=False);m.text(x+.025,1.21,f'{int(v*100)}%',transform=m.transAxes,va='center')
    f.text(*anchor(f,29,57),'Detected cells (%)',fontsize=8)
    scale(f,im,box(f,146,55,23,1.2),'Within-gene mean z',[-2,0,2])
    tracks=read(SD/'F1/F1D_stage_tracks.tsv');tracks=tracks[tracks.track=='Balance']
    for j,name in enumerate(['GSE142522','GSE195812','GSE206710','Park/HTA']):
        a=ax(f,box(f,24+(j%2)*83,108+(j//2)*38,64,22),'E' if j==0 else '',name)
        t=tracks[tracks.dataset==name].sort_values('order');a.plot(t.order,t.value,color=INK,lw=.8);a.scatter(t.order,t.value,s=16,color=INK);a.axhline(0,color=LIGHT,lw=.5)
        labels=[str(v).replace('early cortical','early\ncortical').replace('late cortical','late\ncortical').replace('DP_CD3neg','DP CD3−').replace('DP_CD3pos','DP CD3+').replace('_',' ') for v in t.stage]
        a.set_xticks(t.order,labels,rotation=38,ha='right');a.set_ylabel('Balance\n'+name);grid(a)
    d=ax(f,box(f,24,184,55,24),'F','Paired RNA donors')
    donor=read(D/'a4_yayon/a4a_donor_level_balance.tsv');assert len(donor)==5
    for r in donor.itertuples():d.plot([0,1],[r.early,r.cortical_DP],color=GREY,lw=.85);d.scatter([0,1],[r.early,r.cortical_DP],color=[Z2,Z1],s=20)
    d.set_xticks([0,1],['Early','DP']);d.set_xlim(-.2,1.2);d.set_ylim(-3.2,1.9);d.set_ylabel('ZEB balance');d.text(.02,.98,'n=5\nExact P=.0625',transform=d.transAxes,va='top');grid(d)
    a=ax(f,box(f,115,184,56,24),'G','State balance and tissue context')
    q=read(D/'a4_yayon/a4ab_stage_balance_and_cma.tsv').dropna(subset=['cma_weighted_median']);assert len(q)==19
    for r in q.itertuples():a.scatter(r.cma_weighted_median,r.balance_median,s=np.clip(r.n_cells/250,13,65),color=pal[grp(r.stage)],alpha=.8)
    offsets={'T_ETP':(6,0),'T_DN(early)':(7,-14),'T_DP(P)':(6,0),'T_DP(Q)':(-1,12),'T_CD4':(-35,6),'T_CD8':(-35,-15)}
    for name,off in offsets.items():
        r=q[q.stage==name].iloc[0];a.annotate(name.removeprefix('T_'),(r.cma_weighted_median,r.balance_median),xytext=off,textcoords='offset points',fontsize=8)
    a.set_xlabel('CMA position (cortex → medulla)');a.set_ylabel('RNA stage median\nZEB balance');a.set_xlim(-.8,.85);zeros(a)

    save(f,1)
if __name__=='__main__':main()

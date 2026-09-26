suppressPackageStartupMessages({library(data.table);library(ggplot2);library(patchwork)})
argv<-grep('^--file=',commandArgs(FALSE),value=TRUE)
script<-normalizePath(sub('^--file=','',argv),winslash='/')
root<-normalizePath(file.path(dirname(script),'../../..'),winslash='/')
src<-file.path(root,'total/rebuild_2026/data/source_data_rebuilt/revision4_targeted')
out<-file.path(root,'submission/figures/supplementary')
dir.create(out,recursive=TRUE,showWarnings=FALSE)
th<-theme_classic(base_family='Arial',base_size=8)+theme(
  plot.title=element_text(face='bold',size=9,margin=margin(b=5)),
  axis.title=element_text(size=8),axis.text=element_text(size=7),
  legend.title=element_text(size=7),legend.text=element_text(size=7),
  panel.grid.major.y=element_line(color='#e5eaed',linewidth=.25),
  plot.margin=margin(5,5,5,5))
save_plate<-function(p,stem,w=183,h=130){
  ggsave(file.path(out,paste0(stem,'.pdf')),p,width=w/25.4,height=h/25.4,device=cairo_pdf,bg='white')
  ggsave(file.path(out,paste0(stem,'.png')),p,width=w/25.4,height=h/25.4,dpi=360,bg='white')
  ggsave(file.path(out,paste0(stem,'.tiff')),p,width=w/25.4,height=h/25.4,dpi=360,
         compression='lzw',bg='white')
}
d<-fread(file.path(src,'spline_df_2_3_4.tsv'))
d[,subtype:=factor(subtype,levels=c('BCL11B','ETP-like','TLX3'))]
cols<-c('BCL11B'='#b85732','ETP-like'='#66849a','TLX3'='#28857f')
p1<-ggplot(d,aes(spline_df,median_residual,color=subtype,group=subtype))+
  geom_hline(yintercept=0,color='#bec8cc',linewidth=.35)+
  geom_line(linewidth=.6)+geom_point(size=2.3)+
  scale_x_continuous(breaks=2:4)+scale_color_manual(values=cols,name=NULL)+
  labs(title='A  Residual polarity across spline df',x='Natural spline degrees of freedom',
       y='Median within-cohort residual')+th+
  theme(legend.position='bottom')
p2<-unique(d[,.(spline_df,subtype_incremental_r2,matched_BCL11B_minus_ETP_mean)])
p2a<-ggplot(p2,aes(factor(spline_df),subtype_incremental_r2*100))+
  geom_col(width=.53,fill='#4d778b')+
  geom_text(aes(label=sprintf('%.2f',subtype_incremental_r2*100)),vjust=-.4,size=2.6)+
  coord_cartesian(ylim=c(0,16))+labs(title='B  Added subtype variance',
  x='Spline df',y=expression(Delta*R^2~'(percentage points)'))+th
n<-fread(file.path(src,'normal_proxy_donor_stage.tsv'))
n[,stage:=factor(stage,levels=c('DN_early','DP_immature','DP_mature'))]
p3<-ggplot(n,aes(stage,dev,group=donor,color=donor))+
  geom_line(linewidth=.65)+geom_point(size=2)+
  scale_color_manual(values=c(D1='#476f89',D2='#b67336',D3='#218b76'),name='Donor')+
  labs(title='C  Donor-level normal-thymus proxy',x=NULL,
       y='Three-gene expression score')+th+
  theme(axis.text.x=element_text(angle=25,hjust=1),legend.position='bottom')
save_plate((p1|p2a|p3)+plot_layout(widths=c(1.2,.8,1.2)),'FigureS10',183,100)

g<-fread(file.path(src,'matched_BCL11B_ETP_gene_results.tsv'))
g[,significant:=adj.P.Val<.05 & abs(logFC)>.5]
g[,direction:=fifelse(significant & logFC>0,'Higher in BCL11B',
                      fifelse(significant & logFC<0,'Lower in BCL11B','Other'))]
g[,neglogp:=-log10(pmax(P.Value,1e-300))]
top<-g[significant==TRUE][order(P.Value)][1:min(10,.N)]
pv<-ggplot(g,aes(logFC,neglogp,color=direction))+
  geom_point(size=.55,alpha=.5)+
  ggrepel::geom_text_repel(data=top,aes(label=symbol),size=2.2,
    max.overlaps=Inf,show.legend=FALSE,box.padding=.25,point.padding=.1,
    min.segment.length=0,segment.size=.2,color='#253f4e')+
  geom_vline(xintercept=c(-.5,.5),linetype=2,color='#afbfc6',linewidth=.3)+
  geom_hline(yintercept=-log10(.05),linetype=2,color='#afbfc6',linewidth=.3)+
  scale_color_manual(values=c('Higher in BCL11B'='#b85e34',
    'Lower in BCL11B'='#3b7190','Other'='#b4bec4'),name=NULL)+
  labs(title='A  Pair-adjusted expression contrast',x=expression('BCL11B minus ETP-like log'[2]*' fold change'),
       y=expression(-log[10]*'(nominal P)'))+th+theme(legend.position='bottom')
h<-fread(file.path(src,'matched_BCL11B_ETP_hallmark_camera.tsv'))
h<-h[order(FDR,PValue)][1:10]
h[,label:=gsub('HALLMARK_','',pathway,fixed=TRUE)]
h[,label:=gsub('_',' ',label,fixed=TRUE)]
h[,label:=factor(label,levels=rev(label))]
ph<-ggplot(h,aes(-log10(FDR),label,size=NGenes,color=Direction))+
  geom_segment(data=h,aes(x=0,xend=-log10(FDR),y=label,yend=label),
    inherit.aes=FALSE,linewidth=.45,color='#bdc8cc')+
  geom_point()+
  geom_vline(xintercept=-log10(.05),linetype=2,color='#83969e',linewidth=.4)+
  scale_color_manual(values=c(Up='#b85e34',Down='#3b7190'))+
  scale_size(range=c(2,4),name='Genes')+
  labs(title='B  Hallmark competitive gene-set tests',x=expression(-log[10]*'(BH FDR)'),y=NULL)+th+
  theme(legend.position='right')
save_plate((pv|ph)+plot_layout(widths=c(1,1.25)),'FigureS11',183,112)

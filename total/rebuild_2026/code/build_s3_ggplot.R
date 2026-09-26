# Supplementary Figure S3: frozen non-circular developmental coordinate.
suppressPackageStartupMessages({library(ggplot2);library(patchwork);library(ragg);library(readr)})
argv<-grep('^--file=',commandArgs(FALSE),value=TRUE)
script<-normalizePath(sub('^--file=','',argv),winslash='/')
root<-normalizePath(file.path(dirname(script),'../../..'),winslash='/')
v<-file.path(root,'total/data/validation')
out<-file.path(root,'total/rebuild_2026/edition_20260925/figures/supplementary')
sdir<-file.path(root,'total/rebuild_2026/data/source_data_rebuilt/S3')
dir.create(out,recursive=TRUE,showWarnings=FALSE);dir.create(sdir,recursive=TRUE,showWarnings=FALSE)
rd<-function(p) as.data.frame(read_tsv(file.path(v,p),show_col_types=FALSE,
  progress=FALSE,locale=locale(encoding='UTF-8')))
wt<-function(x,n) write.table(x,file.path(sdir,n),sep='\t',quote=FALSE,row.names=FALSE,na='NA')
freeze<-rd('zeb_developmental_residual/gene_freeze_decision.tsv')
direction<-rd('zeb_developmental_residual/gene_direction_by_dataset.tsv')
units<-rd('zeb_developmental_residual/normal_stage_units_scores.tsv')
pat<-rd('zeb_developmental_residual/polonen_patient_residual.tsv')
curve<-unique(pat[,c('dev','expected_balance')])
curve<-curve[order(curve$dev),]
nr<-range(units$dev)
curve$in_range<-curve$dev>=nr[1]&curve$dev<=nr[2]
stopifnot(nrow(freeze)==12,nrow(units)==20,nrow(pat)==1309)
ink<-'#152b39';teal<-'#168b82';blue<-'#246b96';red<-'#bc3446';purple<-'#776298';pale<-'#e5ecef'
theme_pub<-theme_classic(base_family='Arial',base_size=12.2)+theme(
 plot.title=element_text(colour=ink,face='bold',size=13.6,margin=margin(b=5)),
 plot.subtitle=element_text(colour=ink,size=10.3),
 axis.title=element_text(colour=ink,size=11),axis.text=element_text(colour=ink,size=10),
 legend.text=element_text(colour=ink,size=9.5),plot.margin=margin(8,9,8,9),
 panel.grid.major.y=element_line(colour=pale,linewidth=.35))
freeze$gene<-factor(freeze$gene,levels=rev(freeze$gene))
freeze$selection<-ifelse(freeze$frozen,'Frozen three-gene coordinate','Measured candidate')
pA<-ggplot(freeze,aes(n_measured,gene))+
 geom_vline(xintercept=3,linetype=2,colour='#8da3ad',linewidth=.55)+
 geom_segment(aes(x=0,xend=n_measured,yend=gene),colour=pale,linewidth=1.4)+
 geom_point(aes(colour=selection),size=3.6)+
 geom_text(aes(label=paste0(n_agree,'/',n_measured)),hjust=-.4,colour=ink,size=3.15)+
 scale_colour_manual(values=c('Frozen three-gene coordinate'=red,
                               'Measured candidate'=blue),name=NULL)+
 scale_x_continuous(limits=c(0,4.3),breaks=0:4)+
 labs(title='A  Marker eligibility across four atlases',
      subtitle='Label: datasets with expected direction / datasets measured',
      x='Atlases measuring marker',y=NULL)+theme_pub+
 theme(legend.position='bottom',panel.grid.major.y=element_blank())

direction$gene<-factor(direction$gene,levels=levels(freeze$gene))
direction$dataset<-factor(direction$dataset,
 levels=c('GSE142522','GSE195812','GSE206710','Park/HTA'))
direction$status<-ifelse(is.na(direction$mid_minus_early),'Unavailable',
 ifelse(direction$agrees,'Expected direction','Opposite direction'))
direction$symbol<-ifelse(is.na(direction$mid_minus_early),'',
 ifelse(direction$mid_minus_early>0,'+','-'))
pB<-ggplot(direction,aes(dataset,gene,fill=status))+
 geom_tile(colour='white',linewidth=.7)+
 geom_text(aes(label=symbol),colour=ink,size=3.5,fontface='bold')+
 scale_fill_manual(values=c('Expected direction'='#8dbab6',
                            'Opposite direction'='#ca7680',
                            'Unavailable'='#eef2f2'),name=NULL)+
 labs(title='B  Observed middle-minus-early directions',
      subtitle='CD34/LYL1 expected negative; CD1A expected positive',
      x=NULL,y=NULL)+theme_pub+
 theme(axis.text.x=element_text(angle=35,hjust=1,size=9),
       panel.grid=element_blank(),axis.line=element_blank(),axis.ticks=element_blank(),
       legend.position='bottom')

units$dataset<-factor(units$dataset,levels=c('GSE195812','GSE206710'))
pC<-ggplot(units,aes(dev,balance,colour=dataset))+
 geom_line(data=curve,aes(x=dev,y=expected_balance),inherit.aes=FALSE,
           colour=ink,linewidth=.9)+
 geom_point(size=3,alpha=.9)+
 scale_colour_manual(values=c('GSE195812'=teal,'GSE206710'=purple),name=NULL)+
 labs(title='C  Frozen normal-stage reference',
      subtitle='20 stage units; score = z(CD1A) - mean[z(CD34), z(LYL1)]; R2 = 0.728',
      x='Frozen developmental coordinate',y='Normal-stage balance')+theme_pub+
 theme(legend.position='bottom')

# D is a domain-of-reference audit, not a new score or model.
sub<-rd('zeb_developmental_residual/subtype_residual_summary.tsv')
sub<-sub[order(sub$dev_median),]
clean<-function(x){x<-as.character(x);x[grepl('^LMO2.*-like$',x)]<-'LMO2 gamma-delta-like';
 x[grepl('^TAL1.*-like$',x)&!grepl('DP',x)]<-'TAL1 alpha-beta-like';x}
pat$label<-clean(pat$subtype);sub$label<-clean(sub$subtype)
pat$label<-factor(pat$label,levels=rev(sub$label))
domain<-range(units$dev)
pD<-ggplot(pat,aes(dev,label))+
 annotate('rect',xmin=domain[1],xmax=domain[2],ymin=-Inf,ymax=Inf,
          fill='#e5f0ee',alpha=.72)+
 geom_vline(xintercept=domain,colour='#80aca7',linetype=2,linewidth=.5)+
 geom_point(position=position_jitter(height=.18,width=0,seed=19),
            colour='#8b9fa8',alpha=.38,size=.7)+
 stat_summary(fun=median,geom='point',colour=red,size=2.4)+
 labs(title='D  Patient coordinate against normal support',
      subtitle='1,309 diagnostic patients; shaded: normal-stage coordinate range',
      x='Frozen developmental coordinate',y=NULL)+theme_pub+
 theme(axis.text.y=element_text(size=8.2),panel.grid.major.y=element_blank())

fig<-(pA|pB)/(pC|pD)+plot_layout(heights=c(.95,1.15))
wt(freeze,'S3A_marker_freeze_decision.tsv');wt(direction,'S3B_marker_directions.tsv')
wt(units,'S3C_normal_stage_units.tsv');wt(curve,'S3C_frozen_curve.tsv')
wt(pat[,c('sample_id','subtype','dev','label')],'S3D_patient_coordinate.tsv')
base<-file.path(out,'S3')
ggsave(paste0(base,'.pdf'),fig,width=13.6,height=10.8,units='in',device=cairo_pdf,bg='white')
ggsave(paste0(base,'.svg'),fig,width=13.6,height=10.8,units='in',device=svg,bg='white')
ggsave(paste0(base,'.png'),fig,width=13.6,height=10.8,units='in',device=agg_png,dpi=300,bg='white')
ggsave(paste0(base,'.tiff'),fig,width=13.6,height=10.8,units='in',device=agg_tiff,dpi=300,bg='white')
writeLines(capture.output(sessionInfo()),paste0(base,'_R_sessionInfo.txt'))
manifest<-data.frame(figure='S3',panel=LETTERS[1:4],
 source=c('zeb_developmental_residual/gene_freeze_decision.tsv',
          'zeb_developmental_residual/gene_direction_by_dataset.tsv',
          'zeb_developmental_residual/normal_stage_units_scores.tsv + F3A-B_frozen_expected_curve.tsv',
          'zeb_developmental_residual/polonen_patient_residual.tsv'),
 displayed_data=c('S3A_marker_freeze_decision.tsv','S3B_marker_directions.tsv',
                  'S3C_normal_stage_units.tsv + S3C_frozen_curve.tsv',
                  'S3D_patient_coordinate.tsv'),
 display_transform=c('Count frozen marker availability and directional agreement',
                     'Frozen middle-minus-early sign and missingness',
                     'Reuse frozen normal reference curve; no refit',
                     'Raw patient scores, median markers and frozen normal domain'),
 script='total/rebuild_2026/code/build_s3_ggplot.R')
wt(manifest,'S3_panel_manifest.tsv');cat(base,'\n')

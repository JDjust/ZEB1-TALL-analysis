# Supplementary Figure S9: full GSE146901 expression series, including peripheral controls.
suppressPackageStartupMessages({library(ggplot2);library(patchwork);library(ragg);library(readr)})
argv<-grep('^--file=',commandArgs(FALSE),value=TRUE)
script<-normalizePath(sub('^--file=','',argv),winslash='/')
root<-normalizePath(file.path(dirname(script),'../../..'),winslash='/')
v<-file.path(root,'total/data/validation');out<-file.path(root,'total/rebuild_2026/edition_20260925/figures/supplementary')
sdir<-file.path(root,'total/rebuild_2026/data/source_data_rebuilt/S9')
dir.create(out,recursive=TRUE,showWarnings=FALSE);dir.create(sdir,recursive=TRUE,showWarnings=FALSE)
rd<-function(p) as.data.frame(read_tsv(file.path(v,p),show_col_types=FALSE,
 progress=FALSE,locale=locale(encoding='UTF-8')))
wt<-function(x,n) write.table(x,file.path(sdir,n),sep='\t',quote=FALSE,row.names=FALSE,na='NA')
d<-rd('gse146901/zeb_expression_by_sample.tsv')
test<-rd('gse146901/zeb_expression_etp_vs_nonetp.tsv')
stopifnot(nrow(d)==22,nrow(test)==4,
 sum(d$group=='ETP')==8,sum(d$group=='non-ETP')==10,sum(d$group=='normal T')==4)
ink<-'#152b39';red<-'#bc3446';blue<-'#246b96';teal<-'#168b82';gold<-'#c28830';pale<-'#e5ecee'
pal<-c('ETP'=teal,'non-ETP'=gold,'Peripheral T'=blue)
d$group_label<-ifelse(d$group=='normal T','Peripheral T',d$group)
d$group_label<-factor(d$group_label,levels=c('ETP','non-ETP','Peripheral T'))
d<-d[order(d$group_label,d$balance),]
d$sample<-factor(d$sample,levels=rev(d$sample))
theme_pub<-theme_classic(base_family='Arial',base_size=12.3)+theme(
 plot.title=element_text(colour=ink,face='bold',size=13.6,margin=margin(b=5)),
 plot.subtitle=element_text(colour=ink,size=10.2),axis.title=element_text(colour=ink,size=11),
 axis.text=element_text(colour=ink,size=10),legend.text=element_text(colour=ink,size=9.5),
 plot.margin=margin(8,9,8,9),panel.grid.major.y=element_line(colour=pale,linewidth=.35))
metrics<-c('log2_ZEB1','log2_ZEB2','balance','log2_ratio')
lab<-c('log2_ZEB1'='ZEB1','log2_ZEB2'='ZEB2','balance'='Balance','log2_ratio'='ZEB1/ZEB2 ratio')
heat<-do.call(rbind,lapply(metrics,function(m)data.frame(sample=d$sample,group=d$group_label,
 metric=m,value=d[[m]])))
heat$z<-ave(heat$value,heat$metric,FUN=function(x) as.numeric(scale(x)))
heat$metric<-factor(heat$metric,levels=metrics,labels=lab)
pA<-ggplot(heat,aes(metric,sample,fill=z))+
 geom_tile(colour='white',linewidth=.7)+
 geom_point(aes(x=.55,y=sample,colour=group),shape=15,size=3.3,inherit.aes=FALSE,
            data=unique(heat[,c('sample','group')]))+
 scale_colour_manual(values=pal,name='Group')+
 scale_fill_gradient2(low='#5482a2',mid='#f7f8f5',high='#bf6170',
                      midpoint=0,limits=c(-2.5,2.5),oob=scales::squish,name='Metric z')+
 labs(title='A  All 22 samples, four ZEB metrics',
      subtitle='Each metric z-scored across samples for display',
      x=NULL,y='Sample ID')+theme_pub+
 theme(axis.line=element_blank(),axis.ticks=element_blank(),
       panel.grid=element_blank(),legend.position='bottom')

long<-do.call(rbind,lapply(metrics,function(m)data.frame(sample=d$sample,group=d$group_label,
 metric=lab[[m]],value=d[[m]])))
long$metric<-factor(long$metric,levels=unname(lab))
pB<-ggplot(long,aes(group,value,colour=group))+
 geom_boxplot(width=.45,fill='white',colour='#647d88',outlier.shape=NA)+
 geom_point(position=position_jitter(width=.11,seed=11),size=2.2,alpha=.87)+
 facet_wrap(~metric,scales='free_y',ncol=2)+
 scale_colour_manual(values=pal,guide='none')+
 labs(title='B  Group distributions retain individual samples',
      subtitle='Four healthy controls are peripheral T cells, not normal thymocytes',
      x=NULL,y='Frozen expression metric')+theme_pub+
 theme(axis.text.x=element_text(angle=24,hjust=1,size=9.3),
       strip.text=element_text(colour=ink,face='bold',size=10.7))

test$metric_label<-factor(lab[test$metric],levels=rev(unname(lab)))
test$delta<-test$median_nonETP-test$median_ETP
test$p_label<-ifelse(test$mannwhitney_p<.05,
 sprintf('P = %.3f',test$mannwhitney_p),sprintf('P = %.2f',test$mannwhitney_p))
pC<-ggplot(test,aes(delta,metric_label))+
 geom_vline(xintercept=0,colour='#9cb0b9',linewidth=.5)+
 geom_segment(aes(x=0,xend=delta,yend=metric_label),colour=pale,linewidth=1.7)+
 geom_point(aes(colour=mannwhitney_p<.05),size=3.4)+
 geom_text(aes(label=p_label),hjust=ifelse(test$delta>0,-.16,1.16),
           colour=ink,size=3.3,family='Arial')+
 scale_colour_manual(values=c('TRUE'=red,'FALSE'=blue),guide='none')+
 scale_x_continuous(limits=c(-1.3,2.8))+
 labs(title='C  ETP versus non-ETP tests',
      subtitle='Median non-ETP minus ETP; Mann-Whitney P',
      x='Median difference in each original metric',y=NULL)+theme_pub+
 theme(panel.grid.major.y=element_blank())

fig<-wrap_plots(pA,pB/pC,ncol=2,widths=c(.87,1.13))
wt(d,'S9A_all_22_samples.tsv');wt(heat,'S9A_display_matrix.tsv')
wt(long,'S9B_patient_metrics_long.tsv');wt(test,'S9C_frozen_group_tests.tsv')
base<-file.path(out,'S9')
ggsave(paste0(base,'.pdf'),fig,width=14,height=10.8,units='in',device=cairo_pdf,bg='white')
ggsave(paste0(base,'.svg'),fig,width=14,height=10.8,units='in',device=svg,bg='white')
ggsave(paste0(base,'.png'),fig,width=14,height=10.8,units='in',device=agg_png,dpi=300,bg='white')
ggsave(paste0(base,'.tiff'),fig,width=14,height=10.8,units='in',device=agg_tiff,dpi=300,bg='white')
writeLines(capture.output(sessionInfo()),paste0(base,'_R_sessionInfo.txt'))
manifest<-data.frame(figure='S9',panel=LETTERS[1:3],
 source=c('gse146901/zeb_expression_by_sample.tsv','gse146901/zeb_expression_by_sample.tsv',
          'gse146901/zeb_expression_etp_vs_nonetp.tsv'),
 displayed_data=c('S9A_all_22_samples.tsv + S9A_display_matrix.tsv',
                  'S9B_patient_metrics_long.tsv','S9C_frozen_group_tests.tsv'),
 display_transform=c('All samples; per-metric display z; group sorted only',
                     'Raw frozen metrics with patient-level points',
                     'Frozen median contrasts and Mann-Whitney P; no new test'),
 script='total/rebuild_2026/code/build_s9_ggplot.R')
wt(manifest,'S9_panel_manifest.tsv');cat(base,'\n')

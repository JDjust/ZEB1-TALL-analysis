# Supplementary Figure S5: residual robustness and normal-domain extrapolation audit.
suppressPackageStartupMessages({library(ggplot2);library(patchwork);library(ragg);library(readr)})
argv<-grep('^--file=',commandArgs(FALSE),value=TRUE)
script<-normalizePath(sub('^--file=','',argv),winslash='/')
root<-normalizePath(file.path(dirname(script),'../../..'),winslash='/')
v<-file.path(root,'total/data/validation');out<-file.path(root,'total/rebuild_2026/edition_20260925/figures/supplementary')
sdir<-file.path(root,'total/rebuild_2026/data/source_data_rebuilt/S5')
dir.create(out,recursive=TRUE,showWarnings=FALSE);dir.create(sdir,recursive=TRUE,showWarnings=FALSE)
rd<-function(p) as.data.frame(read_tsv(file.path(v,p),show_col_types=FALSE,
 progress=FALSE,locale=locale(encoding='UTF-8')))
wt<-function(x,n) write.table(x,file.path(sdir,n),sep='\t',quote=FALSE,row.names=FALSE,na='NA')
pat<-rd('zeb_developmental_residual/polonen_patient_residual.tsv')
sub<-rd('zeb_developmental_residual/subtype_residual_summary.tsv')
units<-rd('zeb_developmental_residual/normal_stage_units_scores.tsv')
stopifnot(nrow(pat)==1309,nrow(sub)==17,sum(pat$subtype=='BCL11B')==18)
ink<-'#152b39';blue<-'#246b96';red<-'#bc3446';teal<-'#168b82';purple<-'#776298';pale<-'#e7eeef'
theme_pub<-theme_classic(base_family='Arial',base_size=12)+theme(
 plot.title=element_text(colour=ink,face='bold',size=13.3,margin=margin(b=5)),
 plot.subtitle=element_text(colour=ink,size=10),axis.title=element_text(colour=ink,size=11),
 axis.text=element_text(colour=ink,size=9.8),legend.text=element_text(colour=ink,size=9.5),
 plot.margin=margin(8,9,8,9),panel.grid.major.y=element_line(colour=pale,linewidth=.35))
domain<-range(units$dev);b<-pat[pat$subtype=='BCL11B',]
b$domain<-ifelse(b$dev_outside_normal_range,'Outside normal range','Within normal range')
b$sample_id<-factor(b$sample_id,levels=b$sample_id[order(b$dev)])
pA<-ggplot(b,aes(dev,residual_normal,colour=domain))+
 annotate('rect',xmin=domain[1],xmax=domain[2],ymin=-Inf,ymax=Inf,
          fill='#e5f0ee',alpha=.7)+
 geom_vline(xintercept=domain,colour='#78aaa4',linetype=2,linewidth=.5)+
 geom_point(size=3.5)+
 scale_colour_manual(values=c('Outside normal range'=red,'Within normal range'=teal),name=NULL)+
 labs(title='A  BCL11B patients and reference domain',
      subtitle='14/18 developmental coordinates outside normal-stage range',
      x='Frozen developmental coordinate',y='Normal-reference residual')+theme_pub+
 theme(legend.position='bottom')

long<-rbind(data.frame(sample_id=b$sample_id,estimate='Normal reference',residual=b$residual_normal),
            data.frame(sample_id=b$sample_id,estimate='Within-cohort spline',residual=b$residual_within))
long$estimate<-factor(long$estimate,levels=c('Normal reference','Within-cohort spline'))
pB<-ggplot(long,aes(estimate,residual,group=sample_id))+
 geom_hline(yintercept=0,colour='#aabcc3',linewidth=.5)+
 geom_line(colour='#9aaeb6',alpha=.7,linewidth=.6)+
 geom_point(aes(colour=estimate),size=2.5)+
 stat_summary(aes(group=estimate),fun=median,geom='point',shape=23,
              size=4,fill=ink,colour='white',stroke=.6)+
 scale_colour_manual(values=c('Normal reference'=red,'Within-cohort spline'=blue),guide='none')+
 labs(title='B  Same 18 cases, two frozen residuals',
      subtitle='Median -2.64 versus -1.76; within-cohort fit reduces extrapolation',
      x=NULL,y='Per-patient balance residual')+theme_pub+
 theme(axis.text.x=element_text(size=9.1))

clean<-function(x){x<-as.character(x);x[grepl('^LMO2.*-like$',x)]<-'LMO2 gamma-delta-like';
 x[grepl('^TAL1.*-like$',x)&!grepl('DP',x)]<-'TAL1 alpha-beta-like';x}
sub$label<-clean(sub$subtype)
sub$highlight<-ifelse(sub$label %in% c('BCL11B','ETP-like','SPI1','TLX3'),sub$label,'Other subtypes')
pC<-ggplot(sub,aes(residual_normal_median,residual_within_median))+
 geom_abline(slope=1,intercept=0,colour='#aabcc3',linetype=2,linewidth=.6)+
 geom_point(aes(colour=highlight,size=n),alpha=.88)+
 ggrepel::geom_text_repel(data=sub[sub$highlight!='Other subtypes',],aes(label=label),
   colour=ink,size=3.2,seed=8,max.overlaps=Inf,show.legend=FALSE)+
 scale_colour_manual(values=c('BCL11B'=purple,'ETP-like'=teal,'SPI1'=red,
                              'TLX3'='#c28830','Other subtypes'='#8ba1ab'),guide='none')+
 scale_size_continuous(range=c(2.1,5.3),name='Patients')+
 labs(title='C  Subtype medians across residual definitions',
      subtitle='The two frozen approaches preserve the extreme negative BCL11B state',
      x='Normal-reference residual median',y='Within-cohort residual median')+theme_pub+
 theme(legend.position='bottom')

sub<-sub[order(sub$n_outside_normal_dev_range/sub$n,decreasing=TRUE),]
sub$label<-factor(sub$label,levels=rev(sub$label))
sub$outside_fraction<-sub$n_outside_normal_dev_range/sub$n
pD<-ggplot(sub,aes(outside_fraction,label))+
 geom_segment(aes(x=0,xend=outside_fraction,yend=label),colour=pale,linewidth=1.7)+
 geom_point(aes(colour=label=='BCL11B'),size=3.3)+
 geom_text(aes(label=paste0(n_outside_normal_dev_range,'/',n)),hjust=-.25,
           colour=ink,size=3.1)+
 scale_colour_manual(values=c('TRUE'=purple,'FALSE'=blue),guide='none')+
 scale_x_continuous(labels=scales::percent_format(),limits=c(0,1.18),breaks=c(0,.25,.5,.75,1))+
 labs(title='E  Extrapolation frequency by subtype',
      subtitle='Numerator/denominator printed for every molecular subtype',
      x='Patients outside observed normal-stage coordinate range',y=NULL)+theme_pub+
 theme(axis.text.y=element_text(size=8.7),panel.grid.major.y=element_blank())

sub$label<-as.character(sub$label)
pE<-ggplot(sub,aes(n,residual_normal_median))+
 geom_hline(yintercept=0,colour='#aabcc3',linewidth=.5)+
 geom_point(aes(colour=highlight),size=3.3)+
 ggrepel::geom_text_repel(data=sub[sub$highlight!='Other subtypes',],aes(label=label),
  colour=ink,size=3.2,seed=8,max.overlaps=Inf,show.legend=FALSE)+
 scale_colour_manual(values=c('BCL11B'=purple,'ETP-like'=teal,'SPI1'=red,
                              'TLX3'='#c28830','Other subtypes'='#8ba1ab'),guide='none')+
 scale_x_log10(breaks=c(5,10,20,50,100,300))+
 labs(title='D  Effect size is shown with subtype size',
      subtitle='Small groups remain visible; no smoothed density is used',
      x='Subtype patients (log scale)',y='Normal-reference residual median')+theme_pub

fig<-(pA|pB)/(pC|pE)/pD+plot_layout(heights=c(.9,.9,1.15))
wt(b,'S5A-B_BCL11B_patients.tsv');wt(sub,'S5C-E_subtype_robustness.tsv')
wt(long,'S5B_paired_residuals_long.tsv')
base<-file.path(out,'S5')
ggsave(paste0(base,'.pdf'),fig,width=14.2,height=14.1,units='in',device=cairo_pdf,bg='white')
ggsave(paste0(base,'.svg'),fig,width=14.2,height=14.1,units='in',device=svg,bg='white')
ggsave(paste0(base,'.png'),fig,width=14.2,height=14.1,units='in',device=agg_png,dpi=300,bg='white')
ggsave(paste0(base,'.tiff'),fig,width=14.2,height=14.1,units='in',device=agg_tiff,dpi=300,bg='white')
writeLines(capture.output(sessionInfo()),paste0(base,'_R_sessionInfo.txt'))
manifest<-data.frame(figure='S5',panel=LETTERS[1:5],
 source=c('zeb_developmental_residual/polonen_patient_residual.tsv + normal_stage_units_scores.tsv',
          'zeb_developmental_residual/polonen_patient_residual.tsv',
          'zeb_developmental_residual/subtype_residual_summary.tsv',
          'zeb_developmental_residual/subtype_residual_summary.tsv',
          'zeb_developmental_residual/subtype_residual_summary.tsv'),
 displayed_data=c('S5A-B_BCL11B_patients.tsv','S5B_paired_residuals_long.tsv',
                  'S5C-E_subtype_robustness.tsv','S5C-E_subtype_robustness.tsv',
                  'S5C-E_subtype_robustness.tsv'),
 display_transform=c('All 18 BCL11B patients; frozen normal-domain flag',
                     'Two existing residual definitions paired by patient',
                     'Frozen subtype median residuals','Frozen subtype sizes and normal residual medians',
                     'Frozen outside-domain counts and denominators'),
 script='total/rebuild_2026/code/build_s5_ggplot.R')
wt(manifest,'S5_panel_manifest.tsv');cat(base,'\n')

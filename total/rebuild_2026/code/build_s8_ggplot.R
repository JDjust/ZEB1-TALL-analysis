# Supplementary Figure S8: Lim state contrast and cohort composition.
suppressPackageStartupMessages({library(ggplot2);library(patchwork);library(ragg);library(readr)})
argv<-grep('^--file=',commandArgs(FALSE),value=TRUE)
script<-normalizePath(sub('^--file=','',argv),winslash='/')
root<-normalizePath(file.path(dirname(script),'../../..'),winslash='/')
v<-file.path(root,'total/data/validation');out<-file.path(root,'total/rebuild_2026/edition_20260925/figures/supplementary')
sdir<-file.path(root,'total/rebuild_2026/data/source_data_rebuilt/S8')
dir.create(out,recursive=TRUE,showWarnings=FALSE);dir.create(sdir,recursive=TRUE,showWarnings=FALSE)
rd<-function(p) as.data.frame(read_tsv(file.path(v,p),show_col_types=FALSE,
 progress=FALSE,locale=locale(encoding='UTF-8')))
wt<-function(x,n) write.table(x,file.path(sdir,n),sep='\t',quote=FALSE,row.names=FALSE,na='NA')
pairs<-rd('gateA_lim2025/q3_state_day0_pairs.tsv')
day0<-rd('gateA_lim2025/eligible_patient_timepoint.tsv')
tests<-rd('gateA_lim2025/q3_state.tsv')
stopifnot(nrow(pairs)==41,sum(day0$timepoint=='Day0')==54)
ink<-'#152b39';teal<-'#168b82';purple<-'#776298';blue<-'#246b96';red<-'#bc3446';gold<-'#c28830';pale<-'#e5edef'
theme_pub<-theme_classic(base_family='Arial',base_size=12.2)+theme(
 plot.title=element_text(colour=ink,face='bold',size=13.5,margin=margin(b=5)),
 plot.subtitle=element_text(colour=ink,size=10.2),axis.title=element_text(colour=ink,size=11),
 axis.text=element_text(colour=ink,size=10),legend.text=element_text(colour=ink,size=9.6),
 plot.margin=margin(8,9,8,9),panel.grid.major.y=element_line(colour=pale,linewidth=.35))
pairs$induction<-ifelse(pairs$induction_pos=='IF','IF','Responsive')
pairs$etp<-ifelse(pairs$etp_status_pos=='ETP','ETP','non-ETP')
pairs$rank<-seq_len(nrow(pairs))[order(order(pairs$d_balance))]
pA<-ggplot(pairs,aes(rank,d_balance,colour=induction))+
 geom_hline(yintercept=0,colour='#98acb5',linewidth=.5)+
 geom_segment(aes(y=0,yend=d_balance,xend=rank),linewidth=.65,alpha=.72)+
 geom_point(size=2.7)+
 scale_colour_manual(values=c('IF'=purple,'Responsive'=teal),name=NULL)+
 labs(title='A  Patient-level ZBTB16 state deviations',
      subtitle='41 Day0 pairs; ZBTB16-positive minus negative balance',
      x='Patients ordered by paired difference',y='Paired balance difference')+theme_pub+
 theme(legend.position='bottom',panel.grid.major.x=element_blank())

long<-rbind(data.frame(patient_id=pairs$patient_id_pos,induction=pairs$induction,
                       metric='ZEB1',delta=pairs$d_ZEB1_logcpm),
            data.frame(patient_id=pairs$patient_id_pos,induction=pairs$induction,
                       metric='ZEB2',delta=pairs$d_ZEB2_logcpm),
            data.frame(patient_id=pairs$patient_id_pos,induction=pairs$induction,
                       metric='LMO2',delta=pairs$d_LMO2_logcpm))
long$metric<-factor(long$metric,levels=c('ZEB1','ZEB2','LMO2'))
pB<-ggplot(long,aes(metric,delta,fill=metric))+
 geom_hline(yintercept=0,colour='#98acb5',linewidth=.5)+
 geom_violin(width=.76,alpha=.26,colour=NA,scale='width')+
 geom_point(aes(colour=metric),position=position_jitter(width=.14,seed=8),
            size=1.6,alpha=.62)+
 stat_summary(fun=median,geom='point',shape=23,fill=ink,colour='white',
              stroke=.55,size=3.8)+
 scale_fill_manual(values=c('ZEB1'=red,'ZEB2'=blue,'LMO2'=gold),guide='none')+
 scale_colour_manual(values=c('ZEB1'=red,'ZEB2'=blue,'LMO2'=gold),guide='none')+
 labs(title='B  Constituent gene changes in paired states',
      subtitle='Raw patient differences; diamonds mark medians',
      x=NULL,y='ZBTB16-positive minus negative log2 CPM')+theme_pub

pairs$positive_fraction<-pairs$n_cells_pos/(pairs$n_cells_pos+pairs$n_cells_neg)
pC<-ggplot(pairs,aes(positive_fraction,d_balance,colour=etp,shape=induction))+
 geom_hline(yintercept=0,colour='#98acb5',linewidth=.5)+
 geom_point(size=2.7,alpha=.83)+
 scale_colour_manual(values=c('ETP'=teal,'non-ETP'=gold),name=NULL)+
 scale_shape_manual(values=c('IF'=16,'Responsive'=17),name=NULL)+
 scale_x_continuous(labels=scales::percent_format())+
 labs(title='C  State abundance and paired effect',
      subtitle='Each point is one patient; no fitted association is claimed',
      x='Fraction of malignant cells ZBTB16-positive',y='Paired balance difference')+theme_pub+
 theme(legend.position='bottom')

d0<-day0[day0$timepoint=='Day0',]
comp<-as.data.frame(table(induction=d0$induction,subtype=d0$subtype_l1),stringsAsFactors=FALSE)
comp<-comp[comp$Freq>0,];names(comp)[3]<-'n'
comp$induction<-ifelse(comp$induction=='IF','IF','Responsive')
comp$subtype<-factor(comp$subtype,levels=c('HOXA','LMO','TAL','TLX','other'))
pD<-ggplot(comp,aes(induction,n,fill=subtype))+
 geom_col(width=.62,colour='white',linewidth=.6)+
 geom_text(aes(label=ifelse(n>=2,n,'')),position=position_stack(vjust=.5),
           size=3.2,colour=ink)+
 scale_fill_manual(values=c('HOXA'='#a9849d','LMO'='#758e9d','TAL'=red,
                            'TLX'=blue,'other'='#c2a375'),name='Genomic class')+
 labs(title='D  Baseline genomic-class composition',
      subtitle='54 Day0 patients: 19 IF and 35 responsive',
      x=NULL,y='Patients')+theme_pub+
 theme(legend.position='right')

fig<-(pA|pB)/(pC|pD)
wt(pairs,'S8A-C_Lim_patient_state_pairs.tsv');wt(long,'S8B_gene_differences.tsv')
wt(comp,'S8D_Day0_class_composition.tsv');wt(tests,'S8_frozen_state_tests.tsv')
base<-file.path(out,'S8')
ggsave(paste0(base,'.pdf'),fig,width=13.8,height=10.3,units='in',device=cairo_pdf,bg='white')
ggsave(paste0(base,'.svg'),fig,width=13.8,height=10.3,units='in',device=svg,bg='white')
ggsave(paste0(base,'.png'),fig,width=13.8,height=10.3,units='in',device=agg_png,dpi=300,bg='white')
ggsave(paste0(base,'.tiff'),fig,width=13.8,height=10.3,units='in',device=agg_tiff,dpi=300,bg='white')
writeLines(capture.output(sessionInfo()),paste0(base,'_R_sessionInfo.txt'))
manifest<-data.frame(figure='S8',panel=LETTERS[1:4],
 source=c(rep('gateA_lim2025/q3_state_day0_pairs.tsv',3),
          'gateA_lim2025/eligible_patient_timepoint.tsv'),
 displayed_data=c('S8A-C_Lim_patient_state_pairs.tsv','S8B_gene_differences.tsv',
                  'S8A-C_Lim_patient_state_pairs.tsv','S8D_Day0_class_composition.tsv'),
 display_transform=c('Patient-level paired balance ranked for display',
                     'Frozen within-patient gene differences, violin density and raw points',
                     'Observed cell-state fraction versus paired difference; no model',
                     'Counts of frozen broad genomic class by induction response'),
 script='total/rebuild_2026/code/build_s8_ggplot.R')
wt(manifest,'S8_panel_manifest.tsv');cat(base,'\n')

# Supplementary Figure S7. Lim patient-level QC and Day28 sensitivity.
suppressPackageStartupMessages({library(ggplot2);library(patchwork);library(ragg);library(readr)})
argv<-grep('^--file=',commandArgs(FALSE),value=TRUE)
script<-normalizePath(sub('^--file=','',argv),winslash='/')
root<-normalizePath(file.path(dirname(script),'../../..'),winslash='/')
v<-file.path(root,'total/data/validation');out<-file.path(root,'total/rebuild_2026/edition_20260925/figures/supplementary')
sdir<-file.path(root,'total/rebuild_2026/data/source_data_rebuilt/S7')
dir.create(out,recursive=TRUE,showWarnings=FALSE);dir.create(sdir,recursive=TRUE,showWarnings=FALSE)
rd<-function(p) as.data.frame(read_tsv(file.path(v,p),show_col_types=FALSE,
 progress=FALSE,locale=locale(encoding='UTF-8')))
wt<-function(x,n) write.table(x,file.path(sdir,n),sep='\t',quote=FALSE,row.names=FALSE,na='NA')
inventory<-read.delim(file.path(v,'gateA_lim2025/inventory.tsv'),header=FALSE,
                     sep='\t',stringsAsFactors=FALSE,col.names=c('item','value'))
day0<-rd('gateA_lim2025/eligible_patient_timepoint.tsv')
state<-rd('gateA_lim2025/q3_state_day0_pairs.tsv')
q2<-rd('gateA_lim2025/q2_pairs.tsv')
q2stat<-rd('gateA_lim2025/q2_longitudinal.tsv')
stopifnot(nrow(state)==41,nrow(q2)==12,sum(day0$timepoint=='Day0')==54)
ink<-'#152b39';blue<-'#246b96';teal<-'#168b82';red<-'#bc3446';purple<-'#776298';pale<-'#e7edef'
theme_pub<-theme_classic(base_family='Arial',base_size=12.2)+theme(
 plot.title=element_text(colour=ink,face='bold',size=13.5,margin=margin(b=5)),
 plot.subtitle=element_text(colour=ink,size=10.2),
 axis.title=element_text(colour=ink,size=11),axis.text=element_text(colour=ink,size=10),
 legend.text=element_text(colour=ink,size=9.7),plot.margin=margin(8,9,8,9),
 panel.grid.major.y=element_line(colour=pale,linewidth=.35))

getinv<-function(x) as.numeric(inventory$value[inventory$item==x])
flow<-data.frame(stage=c('Joined cells','Malignant cells','Patients in metadata',
                         'Malignant biological samples','Eligible Day0 samples',
                         'Eligible Day28 samples','Day0 paired-state samples'),
 n=c(getinv('n_cells_joined'),getinv('n_malignant'),getinv('n_patients_meta'),
     getinv('n_biological_samples_malignant'),getinv('n_eligible_day0'),
     getinv('n_eligible_day28'),nrow(state)),
 unit=c('cells','cells','patients','samples','samples','samples','paired patients'))
flow$stage<-factor(flow$stage,levels=rev(flow$stage))
flow$group<-ifelse(flow$unit=='cells','Cells','Patients / samples')
pA<-ggplot(flow,aes(n,stage,colour=group))+
 geom_segment(aes(x=1,xend=n,yend=stage),colour=pale,linewidth=1.6)+
 geom_point(size=3.9)+
 geom_text(aes(label=paste0(scales::comma(n),' ',unit)),hjust=-.15,
           colour=ink,size=3.15,family='Arial')+
 scale_x_log10(limits=c(1,1e7),breaks=c(1,100,10000,1e6),labels=scales::label_comma())+
 scale_colour_manual(values=c('Cells'=red,'Patients / samples'=teal),guide='none')+
 labs(title='A  Lim source and eligible analysis units',
      subtitle='Cell counts and patient/sample counts are distinct units',
      x='Count (log scale)',y=NULL)+theme_pub+
 theme(panel.grid.major.y=element_blank())

d0<-day0[day0$timepoint=='Day0',]
d0$induction<-factor(d0$induction,levels=c('IF','responsive'),
                     labels=c('IF','Responsive'))
pB<-ggplot(d0,aes(induction,n_cells,colour=etp_status))+
 geom_hline(yintercept=getinv('min_blasts'),colour='#9eb2ba',linetype=2)+
 geom_boxplot(outlier.shape=NA,width=.4,colour='#526a76',fill='white')+
 geom_point(position=position_jitter(width=.14,seed=7),size=2.3,alpha=.85)+
 scale_y_log10(labels=scales::label_comma())+
 scale_colour_manual(values=c('ETP'=teal,'nonETP'='#c38d37'),name='ETP status')+
 labs(title='B  Eligible Day0 malignant pseudobulks',
      subtitle='54 patient samples; minimum 50 malignant cells',
      x=NULL,y='Malignant cells per patient')+theme_pub+
 theme(legend.position='bottom')

pC<-ggplot(state,aes(n_cells_neg,n_cells_pos,colour=etp_status_pos))+
 geom_vline(xintercept=getinv('min_state_cells'),colour='#9eb2ba',linetype=2)+
 geom_hline(yintercept=getinv('min_state_cells'),colour='#9eb2ba',linetype=2)+
 geom_point(size=2.7,alpha=.85)+
 scale_x_log10(labels=scales::label_comma())+scale_y_log10(labels=scales::label_comma())+
 scale_colour_manual(values=c('ETP'=teal,'nonETP'='#c38d37'),guide='none')+
 labs(title='C  Both ZBTB16 states pass the cell threshold',
      subtitle='41 paired Day0 malignant pseudobulks; dotted limits = 20 cells',
      x='ZBTB16-negative malignant cells',y='ZBTB16-positive malignant cells')+theme_pub

q2long<-rbind(data.frame(patient_id=q2$patient_id,time='Day0',balance=q2$day0_balance),
              data.frame(patient_id=q2$patient_id,time='Day28',balance=q2$day28_balance))
q2long$time<-factor(q2long$time,levels=c('Day0','Day28'))
pD<-ggplot(q2long,aes(time,balance,group=patient_id))+
 geom_hline(yintercept=0,colour='#9db2bb',linewidth=.45)+
 geom_line(colour='#a1b4bc',linewidth=.7,alpha=.85)+
 geom_point(aes(colour=time),size=2.6)+
 stat_summary(aes(group=time),fun=median,geom='point',shape=23,
              fill=ink,colour='white',stroke=.5,size=4)+
 scale_colour_manual(values=c('Day0'=teal,'Day28'=purple),guide='none')+
 labs(title='D  Day28 change is not a stable main result',x=NULL,
      subtitle='12 paired patients (9 IF, 3 responsive); frozen Wilcoxon P = 0.569',
      y='Frozen ZEB balance')+theme_pub

fig<-(pA|pB)/(pC|pD)
wt(inventory,'S7A_Lim_inventory.tsv');wt(flow,'S7A_eligible_units.tsv')
wt(d0,'S7B_Day0_pseudobulks.tsv');wt(state,'S7C_paired_state_cell_counts.tsv')
wt(q2,'S7D_Day0_Day28_pairs.tsv');wt(q2stat,'S7D_frozen_longitudinal_tests.tsv')
base<-file.path(out,'S7')
ggsave(paste0(base,'.pdf'),fig,width=13.8,height=10.1,units='in',device=cairo_pdf,bg='white')
ggsave(paste0(base,'.svg'),fig,width=13.8,height=10.1,units='in',device=svg,bg='white')
ggsave(paste0(base,'.png'),fig,width=13.8,height=10.1,units='in',device=agg_png,dpi=300,bg='white')
ggsave(paste0(base,'.tiff'),fig,width=13.8,height=10.1,units='in',device=agg_tiff,dpi=300,bg='white')
writeLines(capture.output(sessionInfo()),paste0(base,'_R_sessionInfo.txt'))
manifest<-data.frame(figure='S7',panel=LETTERS[1:4],
 source=c('gateA_lim2025/inventory.tsv',
          'gateA_lim2025/eligible_patient_timepoint.tsv + inventory.tsv',
          'gateA_lim2025/q3_state_day0_pairs.tsv + inventory.tsv',
          'gateA_lim2025/q2_pairs.tsv + q2_longitudinal.tsv'),
 displayed_data=c('S7A_eligible_units.tsv','S7B_Day0_pseudobulks.tsv',
                  'S7C_paired_state_cell_counts.tsv','S7D_Day0_Day28_pairs.tsv'),
 display_transform=c('Frozen inventory counts, biological units identified',
                     'Raw malignant cells per eligible Day0 pseudobulk',
                     'Raw positive/negative malignant cell counts per patient',
                     'Twelve frozen patient pairs displayed without a new test'),
 script='total/rebuild_2026/code/build_s7_ggplot.R')
wt(manifest,'S7_panel_manifest.tsv');cat(base,'\n')

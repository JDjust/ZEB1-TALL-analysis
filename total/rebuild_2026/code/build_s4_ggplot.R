# Supplementary Figure S4: full 17-subtype distributions for individual ZEB-family genes.
suppressPackageStartupMessages({library(ggplot2);library(patchwork);library(ragg);library(readr)})
argv<-grep('^--file=',commandArgs(FALSE),value=TRUE)
script<-normalizePath(sub('^--file=','',argv),winslash='/')
root<-normalizePath(file.path(dirname(script),'../../..'),winslash='/')
v<-file.path(root,'total/data/validation')
out<-file.path(root,'total/rebuild_2026/edition_20260925/figures/supplementary')
sdir<-file.path(root,'total/rebuild_2026/data/source_data_rebuilt/S4')
dir.create(out,recursive=TRUE,showWarnings=FALSE);dir.create(sdir,recursive=TRUE,showWarnings=FALSE)
rd<-function(p) as.data.frame(read_tsv(file.path(v,p),show_col_types=FALSE,
 progress=FALSE,locale=locale(encoding='UTF-8')))
wt<-function(x,n) write.table(x,file.path(sdir,n),sep='\t',quote=FALSE,row.names=FALSE,na='NA')
pat<-rd('polonen_round1/patient_level_frozen_balance.tsv')
res<-rd('zeb_developmental_residual/subtype_residual_summary.tsv')
stopifnot(nrow(pat)==1309,nrow(res)==17)
clean<-function(x){x<-as.character(x);x[grepl('^LMO2.*-like$',x)]<-'LMO2 gamma-delta-like';
 x[grepl('^TAL1.*-like$',x)&!grepl('DP',x)]<-'TAL1 alpha-beta-like';x}
res<-res[order(res$balance_median),]
res$label<-paste0(clean(res$subtype),'  (n=',res$n,')')
pat$label<-res$label[match(pat$subtype,res$subtype)]
stopifnot(!anyNA(pat$label))
pat$label<-factor(pat$label,levels=rev(res$label))
ink<-'#152b39';blue<-'#246b96';red<-'#bc3446';gold<-'#c28830';teal<-'#168b82'
theme_pub<-theme_classic(base_family='Arial',base_size=12.2)+theme(
 plot.title=element_text(colour=ink,face='bold',size=13.5,margin=margin(b=5)),
 plot.subtitle=element_text(colour=ink,size=10.1),
 axis.title=element_text(colour=ink,size=11),axis.text=element_text(colour=ink,size=9.7),
 plot.margin=margin(7,8,8,9),
 panel.grid.major.y=element_line(colour='#e9eef0',linewidth=.3))
draw<-function(field,col,panel,title,xlabel){
 d<-pat[,c('sample_id','subtype','label',field)]
 names(d)[4]<-'value'
 p<-ggplot(d,aes(value,label))+
   geom_point(position=position_jitter(height=.22,width=0,seed=42),
              colour=col,alpha=.30,size=.72)+
   geom_boxplot(width=.48,fill='white',colour='#314b59',linewidth=.48,
                outlier.shape=NA,alpha=.74)+
   stat_summary(fun=median,geom='point',shape=21,fill=col,colour='white',
                stroke=.5,size=2.7)+
   labs(title=paste0(panel,'  ',title),subtitle='One point per diagnostic patient; median and IQR overlaid',
        x=xlabel,y=NULL)+theme_pub+
   theme(axis.text.y=element_text(size=8.7),panel.grid.major.x=element_blank())
 list(plot=p,data=d)
}
a<-draw('ZEB1',red,'A','ZEB1 by molecular subtype','TMM log2 CPM: ZEB1')
b<-draw('ZEB2',blue,'B','ZEB2 by molecular subtype','TMM log2 CPM: ZEB2')
c<-draw('LMO2',gold,'C','LMO2 by molecular subtype','TMM log2 CPM: LMO2')
d<-draw('balance',teal,'D','ZEB1 - ZEB2 balance','z(ZEB1) - z(ZEB2), within cohort')
fig<-(a$plot|b$plot)/(c$plot|d$plot)
wt(a$data,'S4A_ZEB1_patients.tsv');wt(b$data,'S4B_ZEB2_patients.tsv')
wt(c$data,'S4C_LMO2_patients.tsv');wt(d$data,'S4D_balance_patients.tsv')
wt(res[,c('subtype','n','balance_median','label')],'S4_subtype_order_and_n.tsv')
base<-file.path(out,'S4')
ggsave(paste0(base,'.pdf'),fig,width=14.4,height=12.6,units='in',device=cairo_pdf,bg='white')
ggsave(paste0(base,'.svg'),fig,width=14.4,height=12.6,units='in',device=svg,bg='white')
ggsave(paste0(base,'.png'),fig,width=14.4,height=12.6,units='in',device=agg_png,dpi=300,bg='white')
ggsave(paste0(base,'.tiff'),fig,width=14.4,height=12.6,units='in',device=agg_tiff,dpi=300,bg='white')
writeLines(capture.output(sessionInfo()),paste0(base,'_R_sessionInfo.txt'))
manifest<-data.frame(figure='S4',panel=LETTERS[1:4],
 source=rep('polonen_round1/patient_level_frozen_balance.tsv + zeb_developmental_residual/subtype_residual_summary.tsv',4),
 displayed_data=c('S4A_ZEB1_patients.tsv','S4B_ZEB2_patients.tsv',
                  'S4C_LMO2_patients.tsv','S4D_balance_patients.tsv'),
 display_transform=rep('Patient-level raw values; box IQR and median; frozen subtype order by balance',4),
 script='total/rebuild_2026/code/build_s4_ggplot.R')
wt(manifest,'S4_panel_manifest.tsv');cat(base,'\n')

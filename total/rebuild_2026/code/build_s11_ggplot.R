# Supplementary Figure S11: every GSE162280 BCL11B-rearranged case and source counts.
suppressPackageStartupMessages({library(ggplot2);library(patchwork);library(ragg);library(readr)})
argv<-grep('^--file=',commandArgs(FALSE),value=TRUE)
script<-normalizePath(sub('^--file=','',argv),winslash='/')
root<-normalizePath(file.path(dirname(script),'../../..'),winslash='/')
v<-file.path(root,'total/data/validation');out<-file.path(root,'total/rebuild_2026/edition_20260925/figures/supplementary')
sdir<-file.path(root,'total/rebuild_2026/data/source_data_rebuilt/S11')
dir.create(out,recursive=TRUE,showWarnings=FALSE);dir.create(sdir,recursive=TRUE,showWarnings=FALSE)
rd<-function(p) as.data.frame(read_tsv(file.path(v,p),show_col_types=FALSE,
 progress=FALSE,locale=locale(encoding='UTF-8')))
wt<-function(x,n) write.table(x,file.path(sdir,n),sep='\t',quote=FALSE,row.names=FALSE,na='NA')
d<-rd('gse162280/case_manifest_and_expression.tsv')
summary<-rd('gse162280/partner_group_summary.tsv')
stopifnot(nrow(d)==12,sum(d$zeb2_partner=='yes')==5,sum(d$zeb2_partner=='no')==7,
 all(d$log2_ZEB1_over_ZEB2<0))
ink<-'#152b39';red<-'#bc3446';blue<-'#246b96';teal<-'#168b82';purple<-'#776298';gold<-'#c28830';pale<-'#e7edef'
d$lesion<-ifelse(d$zeb2_partner=='yes','ZEB2 fusion',
 ifelse(grepl('ARID1B',d$partner),'ARID1B enhancer','CDK6 enhancer'))
d$lesion<-factor(d$lesion,levels=c('ZEB2 fusion','ARID1B enhancer','CDK6 enhancer'))
d$case<-factor(d$case,levels=rev(d$case))
pal<-c('ZEB2 fusion'=purple,'ARID1B enhancer'=teal,'CDK6 enhancer'=gold)
theme_pub<-theme_classic(base_family='Arial',base_size=12.2)+theme(
 plot.title=element_text(colour=ink,face='bold',size=13.5,margin=margin(b=5)),
 plot.subtitle=element_text(colour=ink,size=10.2),axis.title=element_text(colour=ink,size=11),
 axis.text=element_text(colour=ink,size=9.9),legend.text=element_text(colour=ink,size=9.5),
 plot.margin=margin(8,9,8,9),panel.grid.major.y=element_line(colour=pale,linewidth=.35))

expr<-rbind(data.frame(case=d$case,lesion=d$lesion,gene='ZEB1',cpm=d$ZEB1_cpm),
            data.frame(case=d$case,lesion=d$lesion,gene='ZEB2',cpm=d$ZEB2_cpm))
expr$gene<-factor(expr$gene,levels=c('ZEB1','ZEB2'))
pA<-ggplot(d,aes(y=case))+
 geom_segment(aes(x=ZEB1_cpm,xend=ZEB2_cpm,yend=case),colour='#b3c4ca',linewidth=.9)+
 geom_point(data=expr,aes(x=cpm,colour=gene),size=3.1)+
 scale_x_log10(labels=scales::label_number())+
 scale_colour_manual(values=c('ZEB1'=red,'ZEB2'=blue),name=NULL)+
 labs(title='A  Case-level ZEB1 and ZEB2 RNA abundance',
      subtitle='Counts per million; connected values share one HTSeq library',
      x='Gene-level CPM (log scale)',y=NULL)+theme_pub+
 theme(legend.position='bottom',panel.grid.major.y=element_blank())

pB<-ggplot(d,aes(log2_ZEB1_over_ZEB2,case,colour=lesion))+
 geom_vline(xintercept=0,colour='#8fa5ae',linewidth=.6)+
 geom_segment(aes(x=0,xend=log2_ZEB1_over_ZEB2,yend=case),
              colour=pale,linewidth=1.5)+
 geom_point(size=3.5)+
 scale_colour_manual(values=pal,name=NULL)+
 labs(title='B  All 12 cases are ZEB2-dominant',
      subtitle='Five ZEB2 fusions and seven other partners; zero is equality',
      x='log2(ZEB1 CPM / ZEB2 CPM)',y=NULL)+theme_pub+
 theme(legend.position='bottom',panel.grid.major.y=element_blank())

group<-aggregate(log2_ZEB1_over_ZEB2~lesion,d,function(x)c(n=length(x),min=min(x),median=median(x),max=max(x)))
group<-data.frame(lesion=group$lesion,group$log2_ZEB1_over_ZEB2,check.names=FALSE)
pC<-ggplot(d,aes(lesion,log2_ZEB1_over_ZEB2,colour=lesion))+
 geom_hline(yintercept=0,colour='#8fa5ae',linewidth=.6)+
 geom_point(position=position_jitter(width=.09,seed=9),size=3.4)+
 stat_summary(fun=median,geom='point',shape=23,fill=ink,colour='white',
              stroke=.55,size=4)+
 scale_colour_manual(values=pal,guide='none')+
 labs(title='C  Convergence spans three rearrangement classes',
      subtitle='Non-ZEB2-partner groups comprise five ARID1B and two CDK6 cases',
      x=NULL,y='log2(ZEB1 CPM / ZEB2 CPM)')+theme_pub+
 theme(axis.text.x=element_text(angle=25,hjust=1,size=9.4))

pD<-ggplot(d,aes(library_counts/1e6,log2_ZEB1_over_ZEB2,
                 colour=lesion,shape=grepl('ETP-ALL',phenotype)))+
 geom_hline(yintercept=0,colour='#8fa5ae',linewidth=.6)+
 geom_point(size=3.7)+
 ggrepel::geom_text_repel(aes(label=case),colour=ink,size=3.0,
                          seed=12,max.overlaps=Inf,show.legend=FALSE)+
 scale_colour_manual(values=pal,guide='none')+
 scale_shape_manual(values=c('FALSE'=16,'TRUE'=17),labels=c('Other / unknown','ETP-ALL'),
                    name='Phenotype')+
 labs(title='D  Ratio is shown against library depth',
      subtitle='AML, MPAL and ETP-ALL are all represented; one phenotype is unlabeled',
      x='Library assigned counts (million)',y='log2(ZEB1 CPM / ZEB2 CPM)')+theme_pub+
 theme(legend.position='bottom')

fig<-(pA|pB)/(pC|pD)
wt(d,'S11_full_case_manifest_and_counts.tsv');wt(expr,'S11A_case_CPM_long.tsv')
wt(group,'S11C_partner_group_range.tsv');wt(summary,'S11_frozen_partner_summary.tsv')
base<-file.path(out,'S11')
ggsave(paste0(base,'.pdf'),fig,width=13.8,height=10.2,units='in',device=cairo_pdf,bg='white')
ggsave(paste0(base,'.svg'),fig,width=13.8,height=10.2,units='in',device=svg,bg='white')
ggsave(paste0(base,'.png'),fig,width=13.8,height=10.2,units='in',device=agg_png,dpi=300,bg='white')
ggsave(paste0(base,'.tiff'),fig,width=13.8,height=10.2,units='in',device=agg_tiff,dpi=300,bg='white')
writeLines(capture.output(sessionInfo()),paste0(base,'_R_sessionInfo.txt'))
manifest<-data.frame(figure='S11',panel=LETTERS[1:4],
 source=rep('gse162280/case_manifest_and_expression.tsv',4),
 displayed_data=c('S11A_case_CPM_long.tsv','S11_full_case_manifest_and_counts.tsv',
                  'S11C_partner_group_range.tsv','S11_full_case_manifest_and_counts.tsv'),
 display_transform=c('ZEB1/ZEB2 gene-level CPM paired by case',
                     'Frozen case-level ratio ordered by case',
                     'Frozen ratio by lesion group; group range and median only',
                     'Frozen library counts and phenotype annotation'),
 script='total/rebuild_2026/code/build_s11_ggplot.R')
wt(manifest,'S11_panel_manifest.tsv');cat(base,'\n')

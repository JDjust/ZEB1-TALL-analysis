# Supplementary Figure S12: scoped historical cohort contrasts only.
suppressPackageStartupMessages({library(ggplot2);library(patchwork);library(ragg);library(readr)})
argv<-grep('^--file=',commandArgs(FALSE),value=TRUE)
script<-normalizePath(sub('^--file=','',argv),winslash='/')
root<-normalizePath(file.path(dirname(script),'../../..'),winslash='/')
v<-file.path(root,'total/data/validation');out<-file.path(root,'total/rebuild_2026/edition_20260925/figures/supplementary')
sdir<-file.path(root,'total/rebuild_2026/data/source_data_rebuilt/S12')
dir.create(out,recursive=TRUE,showWarnings=FALSE);dir.create(sdir,recursive=TRUE,showWarnings=FALSE)
rd<-function(p) as.data.frame(read_tsv(file.path(v,p),show_col_types=FALSE,
 progress=FALSE,locale=locale(encoding='UTF-8')))
wt<-function(x,n) write.table(x,file.path(sdir,n),sep='\t',quote=FALSE,row.names=FALSE,na='NA')
co<-rd('zeb_family_state/cohort_effects.tsv');meta<-rd('zeb_family_state/random_effects_meta.tsv')
loo<-rd('zeb_family_state/leave_one_cohort_out.tsv')
keys<-data.frame(gene=c('ZEB1','LMO2','ZEB2'),subtype=c('TLX','TLX','TAL'),
 label=c('TLX : ZEB1 high','TLX : LMO2 low','TAL : ZEB2 low'))
pick<-function(d){x<-merge(d,keys,by=c('gene','subtype'));x$label<-factor(x$label,levels=keys$label);x}
co<-pick(co);meta<-pick(meta);loo<-pick(loo)
co$cohort<-factor(co$cohort,levels=rev(c('GSE110636','GSE26713','GSE62156','TARGET-ALL-P2')))
meta$method<-factor(meta$method,levels=c('DL_normal','REML_modified_HK'))
meta$ypos<-as.numeric(meta$label)+ifelse(meta$method=='DL_normal',-.14,.14)
ink<-'#152b39';red<-'#bc3446';blue<-'#246b96';teal<-'#168b82';pale<-'#e5ecee'
cols<-c('TLX : ZEB1 high'=red,'TLX : LMO2 low'=teal,'TAL : ZEB2 low'=blue)
theme_pub<-theme_classic(base_family='Arial',base_size=12.2)+theme(
 plot.title=element_text(colour=ink,face='bold',size=13.5,margin=margin(b=5)),
 plot.subtitle=element_text(colour=ink,size=10.1),axis.title=element_text(colour=ink,size=11),
 axis.text=element_text(colour=ink,size=9.8),legend.text=element_text(colour=ink,size=9.5),
 strip.text=element_text(colour=ink,face='bold',size=10.2),plot.margin=margin(8,9,8,9),
 panel.grid.major.y=element_line(colour=pale,linewidth=.35))
pA<-ggplot(co,aes(hedges_g,cohort,colour=label))+
 geom_vline(xintercept=0,colour='#9dafb7',linewidth=.5)+
 geom_segment(aes(x=ci95_lo,xend=ci95_hi,yend=cohort),linewidth=.8)+
 geom_point(size=3.1)+
 facet_wrap(~label,ncol=1)+
 scale_colour_manual(values=cols,guide='none')+
 labs(title='A  Four-cohort directional context',
      subtitle='Subtype versus rest Hedges g, 95% CI',
      x='Cohort effect size',y=NULL)+theme_pub+
 theme(panel.grid.major.y=element_blank())

pB<-ggplot(meta,aes(hedges_g,ypos,colour=method))+
 geom_vline(xintercept=0,colour='#9dafb7',linewidth=.5)+
 geom_segment(aes(x=ci95_lo,xend=ci95_hi,yend=ypos),linewidth=.9)+
 geom_point(size=3)+
 scale_y_continuous(breaks=seq_along(levels(meta$label)),
                    labels=levels(meta$label))+
 scale_colour_manual(values=c('DL_normal'=red,'REML_modified_HK'=blue),
                     labels=c('DerSimonian-Laird','REML / modified H-K'),name=NULL)+
 labs(title='B  Historical cohort meta-analysis',
      subtitle='TLX ZEB1 modified H-K interval crosses zero',
      x='Pooled Hedges g (95% interval)',y=NULL)+theme_pub+
 theme(legend.position='bottom',axis.text.y=element_text(size=9.2))

lo<-loo[loo$method=='DL_normal',]
lo$omitted_cohort<-factor(lo$omitted_cohort,levels=c('GSE110636','GSE26713',
                                                     'GSE62156','TARGET-ALL-P2'))
pC<-ggplot(lo,aes(hedges_g,omitted_cohort,colour=label))+
 geom_vline(xintercept=0,colour='#9dafb7',linewidth=.5)+
 geom_segment(aes(x=ci95_lo,xend=ci95_hi,yend=omitted_cohort),linewidth=.8)+
 geom_point(size=2.9)+
 facet_wrap(~label,ncol=1)+
 scale_colour_manual(values=cols,guide='none')+
 labs(title='C  Leave-one-out stability',
      subtitle='Each point omits the named cohort (DL)',
      x='Three-cohort pooled Hedges g',y='Omitted cohort')+theme_pub+
 theme(panel.grid.major.y=element_blank(),axis.text.y=element_text(size=9.1))

fig<-wrap_plots(pA,pB,pC,ncol=3,widths=c(1.08,1,1.08))
wt(co,'S12A_cohort_effects.tsv');wt(meta,'S12B_meta_method_sensitivity.tsv')
wt(lo,'S12C_leave_one_out.tsv')
base<-file.path(out,'S12')
ggsave(paste0(base,'.pdf'),fig,width=15.0,height=8.3,units='in',device=cairo_pdf,bg='white')
ggsave(paste0(base,'.svg'),fig,width=15.0,height=8.3,units='in',device=svg,bg='white')
ggsave(paste0(base,'.png'),fig,width=15.0,height=8.3,units='in',device=agg_png,dpi=300,bg='white')
ggsave(paste0(base,'.tiff'),fig,width=15.0,height=8.3,units='in',device=agg_tiff,dpi=300,bg='white')
writeLines(capture.output(sessionInfo()),paste0(base,'_R_sessionInfo.txt'))
manifest<-data.frame(figure='S12',panel=LETTERS[1:3],
 source=c('zeb_family_state/cohort_effects.tsv','zeb_family_state/random_effects_meta.tsv',
          'zeb_family_state/leave_one_cohort_out.tsv'),
 displayed_data=c('S12A_cohort_effects.tsv','S12B_meta_method_sensitivity.tsv',
                  'S12C_leave_one_out.tsv'),
 display_transform=c('Only frozen TLX ZEB1, TLX LMO2, TAL ZEB2 contrasts',
                     'Frozen DL and REML/modified HK models; no refit',
                     'Frozen DL leave-one-out summaries'),
 script='total/rebuild_2026/code/build_s12_ggplot.R')
wt(manifest,'S12_panel_manifest.tsv');cat(base,'\n')

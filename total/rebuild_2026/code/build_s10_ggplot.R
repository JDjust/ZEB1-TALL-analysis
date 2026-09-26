# Supplementary Figure S10: sample-level HiChIP and negative orthogonal assays.
suppressPackageStartupMessages({library(ggplot2);library(patchwork);library(ragg);library(readr)})
argv<-grep('^--file=',commandArgs(FALSE),value=TRUE)
script<-normalizePath(sub('^--file=','',argv),winslash='/')
root<-normalizePath(file.path(dirname(script),'../../..'),winslash='/')
v<-file.path(root,'total/data/validation');out<-file.path(root,'total/rebuild_2026/edition_20260925/figures/supplementary')
sdir<-file.path(root,'total/rebuild_2026/data/source_data_rebuilt/S10')
dir.create(out,recursive=TRUE,showWarnings=FALSE);dir.create(sdir,recursive=TRUE,showWarnings=FALSE)
rd<-function(p) as.data.frame(read_tsv(file.path(v,p),show_col_types=FALSE,
 progress=FALSE,locale=locale(encoding='UTF-8')))
wt<-function(x,n) write.table(x,file.path(sdir,n),sep='\t',quote=FALSE,row.names=FALSE,na='NA')
h<-rd('zeb_chromatin_reciprocal/hichip_zeb1_zeb2_contacts.tsv')
a<-rd('zeb_chromatin_reciprocal/scatac_zeb1_zeb2_peaks.tsv')
l<-rd('gse146901/zeb_merged_loop_rates.tsv')
stopifnot(nrow(h)==10,nrow(a)==7,nrow(l)==3)
ink<-'#152b39';red<-'#bc3446';blue<-'#246b96';teal<-'#168b82';gold<-'#c28830';purple<-'#776298';pale<-'#e5ecee'
group_pal<-c('non-ETP T-ALL'=blue,'ETP'=teal,'CD34'=gold,'thymus'=purple)
group_label<-function(x) ifelse(x=='non-ETP T-ALL','non-ETP',x)
theme_pub<-theme_classic(base_family='Arial',base_size=12)+theme(
 plot.title=element_text(colour=ink,face='bold',size=13.2,margin=margin(b=5)),
 plot.subtitle=element_text(colour=ink,size=10),axis.title=element_text(colour=ink,size=11),
 axis.text=element_text(colour=ink,size=9.8),legend.text=element_text(colour=ink,size=9.5),
 plot.margin=margin(8,9,8,9),panel.grid.major.y=element_line(colour=pale,linewidth=.35))
h$group<-factor(h$group,levels=c('non-ETP T-ALL','ETP','CD34','thymus'))
h$million_contacts<-h$n_contacts/1e6
pA<-ggplot(h,aes(million_contacts,reorder(sample,million_contacts),colour=group))+
 geom_segment(aes(x=0,xend=million_contacts,yend=reorder(sample,million_contacts)),
              colour=pale,linewidth=1.4)+
 geom_point(size=3.4)+
 scale_colour_manual(values=group_pal,name=NULL)+
 labs(title='A  HiChIP library depth and group',
      subtitle='10 usable contact libraries; one additional thymus archive was damaged',
      x='Filtered contacts (million)',y='Sample')+theme_pub+
 theme(legend.position='bottom',panel.grid.major.y=element_blank())

pB<-ggplot(h,aes(ZEB1_per_M,ZEB2_per_M,colour=group))+
 geom_abline(slope=1,intercept=0,colour='#aabcc3',linetype=2,linewidth=.5)+
 geom_point(size=3.4)+
 ggrepel::geom_text_repel(aes(label=sample),colour=ink,size=2.9,
                          seed=9,max.overlaps=Inf,show.legend=FALSE)+
 scale_colour_manual(values=group_pal,name=NULL)+
 labs(title='B  Both loci, normalized for library depth',
      subtitle='Contacts per million filtered contacts; each point is one sample',
      x='ZEB1 gene-body contacts / million',y='ZEB2 gene-body contacts / million')+theme_pub+
 theme(legend.position='bottom')

h$ZEB1_prom_per_M<-h$ZEB1_prom/h$n_contacts*1e6
h$ZEB2_prom_per_M<-h$ZEB2_prom/h$n_contacts*1e6
pC<-ggplot(h,aes(ZEB1_prom_per_M,ZEB2_prom_per_M,colour=group))+
 geom_abline(slope=1,intercept=0,colour='#aabcc3',linetype=2,linewidth=.5)+
 geom_point(size=3.4)+
 scale_colour_manual(values=group_pal,name=NULL)+
 labs(title='C  Promoter-window contacts are a separate layer',
      subtitle='Display normalization only; locus proximity does not establish promoter action',
      x='ZEB1 promoter-window contacts / million',
      y='ZEB2 promoter-window contacts / million')+theme_pub+
 theme(legend.position='bottom')

a$group<-factor(a$group,levels=c('non-ETP T-ALL','ETP'))
pD<-ggplot(a,aes(ZEB1_per_10k_peaks,ZEB2_per_10k_peaks,colour=group,shape=source))+
 geom_abline(slope=1,intercept=0,colour='#aabcc3',linetype=2,linewidth=.5)+
 geom_point(size=3.6)+
 ggrepel::geom_text_repel(aes(label=sample),colour=ink,size=2.9,
                          seed=9,max.overlaps=Inf,show.legend=FALSE)+
 scale_colour_manual(values=group_pal,name=NULL)+
 labs(title='D  scATAC peak sets do not independently replicate',
      subtitle='Five non-ETP and two ETP peak sets; mixed primary/PDX sources',
      x='ZEB1 peaks / 10,000 peaks',y='ZEB2 peaks / 10,000 peaks')+theme_pub+
 theme(legend.position='bottom')

features<-c('ZEB1_body_per_1000','ZEB2_body_per_1000',
            'ZEB1_locus_1Mb_per_1000','ZEB2_locus_1Mb_per_1000')
ll<-do.call(rbind,lapply(features,function(f)data.frame(group=l$group,
                         feature=f,rate=l[[f]],n_loops=l$n_loops)))
ll$feature<-factor(ll$feature,levels=features,
 labels=c('ZEB1 body','ZEB2 body','ZEB1 1-Mb locus','ZEB2 1-Mb locus'))
ll$group<-factor(ll$group,levels=c('ETP','non-ETP','normal T'))
pE<-ggplot(ll,aes(group,rate,colour=group))+
 geom_hline(yintercept=0,colour='#afc0c6',linewidth=.4)+
 geom_point(size=3.2)+
 facet_wrap(~feature,scales='free_y',nrow=1)+
 scale_colour_manual(values=c('ETP'=teal,'non-ETP'=blue,'normal T'=gold),guide='none')+
 labs(title='E  Pooled GSE146901 loops: no reciprocal 3D replication',
      subtitle='Rates per 1,000 loops; pooled group calls are not patient-level replicates',
      x=NULL,y='Loop calls / 1,000')+theme_pub+
 theme(strip.text=element_text(colour=ink,face='bold',size=10),
       axis.text.x=element_text(angle=35,hjust=1,size=8.9))

fig<-(pA|pB)/(pC|pD)/pE+plot_layout(heights=c(.9,.9,.8))
wt(h,'S10A-C_HiChIP_sample_contacts.tsv');wt(a,'S10D_scATAC_peak_sets.tsv')
wt(l,'S10E_pooled_loop_original.tsv');wt(ll,'S10E_pooled_loop_long.tsv')
base<-file.path(out,'S10')
ggsave(paste0(base,'.pdf'),fig,width=14.2,height=13.1,units='in',device=cairo_pdf,bg='white')
ggsave(paste0(base,'.svg'),fig,width=14.2,height=13.1,units='in',device=svg,bg='white')
ggsave(paste0(base,'.png'),fig,width=14.2,height=13.1,units='in',device=agg_png,dpi=300,bg='white')
ggsave(paste0(base,'.tiff'),fig,width=14.2,height=13.1,units='in',device=agg_tiff,dpi=300,bg='white')
writeLines(capture.output(sessionInfo()),paste0(base,'_R_sessionInfo.txt'))
manifest<-data.frame(figure='S10',panel=LETTERS[1:5],
 source=c(rep('zeb_chromatin_reciprocal/hichip_zeb1_zeb2_contacts.tsv',3),
          'zeb_chromatin_reciprocal/scatac_zeb1_zeb2_peaks.tsv',
          'gse146901/zeb_merged_loop_rates.tsv'),
 displayed_data=c(rep('S10A-C_HiChIP_sample_contacts.tsv',3),
                  'S10D_scATAC_peak_sets.tsv','S10E_pooled_loop_long.tsv'),
 display_transform=c('Frozen sample depth and group','Frozen gene-body contacts per million',
                     'Promoter-window counts divided by frozen library contacts',
                     'Frozen scATAC peak rates','Frozen pooled loop rates reshaped; no new calls'),
 script='total/rebuild_2026/code/build_s10_ggplot.R')
wt(manifest,'S10_panel_manifest.tsv');cat(base,'\n')

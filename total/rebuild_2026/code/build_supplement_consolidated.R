# Consolidates the frozen, individually reviewed R supplement panels into seven
# high-density figures. Original scripts remain the source of every panel.
suppressPackageStartupMessages({library(ggplot2);library(patchwork);library(ragg)})
argv<-grep('^--file=',commandArgs(FALSE),value=TRUE)
script<-normalizePath(sub('^--file=','',argv),winslash='/')
code<-dirname(script)
root<-normalizePath(file.path(code,'../../..'),winslash='/')
out<-file.path(root,'total/rebuild_2026/edition_20260925/figures/supplementary')
data_root<-file.path(root,'total/rebuild_2026/data/source_data_rebuilt')
dir.create(out,recursive=TRUE,showWarnings=FALSE)
source_panel_script<-function(n){
 e<-new.env(parent=.GlobalEnv)
 # The prototype scripts calculate their panel objects and source tables.
 # Suppress only their old numbered exports; final exports are below.
 e$ggsave<-function(...) invisible(NULL)
 sys.source(file.path(code,paste0('build_s',n,'_ggplot.R')),envir=e)
 e
}
memo<-new.env(parent=emptyenv())
getenv<-function(n){k<-as.character(n)
 if(!exists(k,envir=memo,inherits=FALSE)) assign(k,source_panel_script(n),envir=memo)
 get(k,envir=memo,inherits=FALSE)}
getpanel<-function(n,key,letter){
 e<-getenv(n)
 p<-if(grepl('\\$',key)){
   parts<-strsplit(key,'\\$')[[1]]
   e[[parts[1]]][[parts[2]]]
 } else e[[paste0('p',key)]]
 if(!inherits(p,'ggplot')) stop('Missing ggplot panel: ',n,'/',key)
 ttl<-p$labels$title
 if(is.null(ttl)) ttl<-paste0('Source S',n,' panel ',key)
 p<-p+labs(title=paste0(letter,'  ',sub('^[A-Z]  +','',ttl)))
 if(n==5 && key=='C') p<-p+scale_x_continuous(expand=expansion(mult=c(.20,.08)))+
   coord_cartesian(clip='off')
 p
}
spec<-list(
 S1=list(ref=c('2:A','2:B','3:A','2:C','2:D','2:E','3:B','3:D'),
         size=c(19,18),layout='normal_score'),
 S2=list(ref=c('4:a$plot','4:b$plot','4:c$plot','4:d$plot',
               '5:A','5:B','5:C','5:D'),size=c(19,19),layout='grid4'),
 S3=list(ref=c('7:A','7:B','7:C','7:D','8:A','8:B','8:C','8:D'),
         size=c(19,19),layout='grid4'),
 S4=list(ref=c('9:A','9:B','9:C'),size=c(18,11),layout='bulk'),
 S5=list(ref=c('10:A','10:B','10:C','10:D','10:E'),
         size=c(18,16),layout='chromatin'),
 S6=list(ref=c('11:A','11:B','11:C','11:D'),size=c(17,14),layout='grid2'),
 S7=list(ref=c('12:A','12:B','12:C'),size=c(18,7),layout='horizontal'))
for(fig_name in names(spec)){
 sp<-spec[[fig_name]]
 parts<-strsplit(sp$ref,':',fixed=TRUE)
 panels<-lapply(seq_along(parts),function(i){
   getpanel(as.integer(parts[[i]][1]),parts[[i]][2],LETTERS[i])})
 fig<-switch(sp$layout,
  normal_score=panels[[1]]/
   (panels[[2]]|panels[[3]])/
   (panels[[4]]|panels[[5]]|panels[[6]])/
   (panels[[7]]|panels[[8]])+
   plot_layout(heights=c(.80,1.12,.82,1.12)),
  grid4=(panels[[1]]|panels[[2]])/
   (panels[[3]]|panels[[4]])/
   (panels[[5]]|panels[[6]])/
   (panels[[7]]|panels[[8]]),
  bulk=panels[[1]]|(panels[[2]]/panels[[3]])+
   plot_layout(widths=c(.91,1.09)),
  chromatin=(panels[[1]]|panels[[2]])/
   (panels[[3]]|panels[[4]])/panels[[5]]+
   plot_layout(heights=c(.9,.9,.8)),
  grid2=(panels[[1]]|panels[[2]])/(panels[[3]]|panels[[4]]),
  horizontal=wrap_plots(panels,ncol=3,widths=c(1.08,1,1.08)))
 # All final panels use dark axis text even if an old prototype overrode its theme.
 fig<-fig & theme(axis.text=element_text(colour='#152b39'),
                  axis.title=element_text(colour='#152b39'))
 base<-file.path(out,fig_name)
 ggsave(paste0(base,'.pdf'),fig,width=sp$size[1],height=sp$size[2],
        units='in',device=cairo_pdf,bg='white',limitsize=FALSE)
 ggsave(paste0(base,'.svg'),fig,width=sp$size[1],height=sp$size[2],
        units='in',device=svg,bg='white',limitsize=FALSE)
 ggsave(paste0(base,'.png'),fig,width=sp$size[1],height=sp$size[2],
        units='in',device=agg_png,dpi=300,bg='white',limitsize=FALSE)
 ggsave(paste0(base,'.tiff'),fig,width=sp$size[1],height=sp$size[2],
        units='in',device=agg_tiff,dpi=300,bg='white',limitsize=FALSE)
 final_sd<-file.path(data_root,paste0('Supplement_',fig_name))
 dir.create(final_sd,recursive=TRUE,showWarnings=FALSE)
 entries<-vector('list',length(parts))
 for(i in seq_along(parts)){
   n<-as.integer(parts[[i]][1]);key<-parts[[i]][2]
   old_panel<-if(grepl('\\$',key)) toupper(substr(key,1,1)) else key
   legacy<-file.path(data_root,paste0('S',n))
   pm<-read.delim(file.path(legacy,paste0('S',n,'_panel_manifest.tsv')),
                  check.names=FALSE,stringsAsFactors=FALSE)
   row<-pm[pm$panel==old_panel,,drop=FALSE]
   if(nrow(row)!=1) stop('Ambiguous source manifest: ',n,'/',old_panel)
   old_names<-strsplit(row$displayed_data,' + ',fixed=TRUE)[[1]]
   final_files<-paste0('source_S',n,'_',basename(old_names))
   for(j in seq_along(old_names)){
     old_file<-file.path(legacy,old_names[j])
     if(!file.copy(old_file,file.path(final_sd,final_files[j]),overwrite=TRUE))
        stop('Missing source file: ',old_file)
   }
   entries[[i]]<-data.frame(figure=fig_name,panel=LETTERS[i],
     source=row$source,displayed_data=paste(final_files,collapse=' + '),
     display_transform=row$display_transform,
     panel_script=row$script,
     assembly_script='total/rebuild_2026/code/build_supplement_consolidated.R')
 }
 manifest<-do.call(rbind,entries)
 write.table(manifest,file.path(final_sd,paste0(fig_name,'_panel_manifest.tsv')),
             sep='\t',quote=FALSE,row.names=FALSE)
 write.table(manifest,paste0(base,'_panel_manifest.tsv'),
             sep='\t',quote=FALSE,row.names=FALSE)
 writeLines(capture.output(sessionInfo()),paste0(base,'_R_sessionInfo.txt'))
 cat('Built ',fig_name,' with ',length(panels),' panels\n',sep='')
}

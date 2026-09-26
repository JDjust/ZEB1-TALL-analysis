# Editorial relayout of frozen validation panels into three coherent main figures.
# No scientific model is recomputed; the panel-producing R scripts are sourced.
suppressPackageStartupMessages({library(ggplot2);library(patchwork);library(ragg)})
argv<-grep('^--file=',commandArgs(FALSE),value=TRUE)
script<-normalizePath(sub('^--file=','',argv),winslash='/')
code<-dirname(script)
root<-normalizePath(file.path(code,'../../..'),winslash='/')
out<-file.path(root,'total/rebuild_2026/edition_20260925/figures/main')
data_root<-file.path(root,'total/rebuild_2026/data/source_data_rebuilt')
dir.create(out,recursive=TRUE,showWarnings=FALSE)
sources<-c(F5='build_f5_ggplot.R',F6='build_f6_ggplot.R',
           S8='build_s8_ggplot.R',S9='build_s9_ggplot.R',
           S12='build_s12_ggplot.R')
envs<-list()
for(k in names(sources)){
 e<-new.env(parent=.GlobalEnv)
 e$ggsave<-function(...) invisible(NULL)
 sys.source(file.path(code,sources[[k]]),envir=e)
 envs[[k]]<-e
}
legacy_manifests<-lapply(names(sources),function(k){
 path<-if(startsWith(k,'F')) file.path(out,paste0(k,'_panel_manifest.tsv'))
       else file.path(data_root,k,paste0(k,'_panel_manifest.tsv'))
 read.delim(path,check.names=FALSE,stringsAsFactors=FALSE)
})
names(legacy_manifests)<-names(sources)
getpanel<-function(ref,letter){
 part<-strsplit(ref,':',fixed=TRUE)[[1]]
 p<-envs[[part[1]]][[paste0('p',part[2])]]
 if(!inherits(p,'ggplot')) stop('Missing panel ',ref)
 ttl<-p$labels$title
 p+labs(title=paste0(letter,'  ',sub('^[A-Z]  +','',ttl)))+
   theme(axis.text=element_text(colour='#152b39'),
         axis.title=element_text(colour='#152b39'))
}
spec<-list(
 F5=list(ref=c('F5:A','F5:B','F5:C','S9:A','S9:C','S12:B'),
         size=c(17,12.4),layout='bulk'),
 F6=list(ref=c('F5:D','S8:A','S8:B','S8:C',
               'F5:E','F5:F','F5:G','S8:D'),
         size=c(17,16),layout='lim'),
 F7=list(ref=paste0('F6:',LETTERS[1:8]),
         size=c(14.4,11.4),layout='chromatin'))
for(fig_name in names(spec)){
 sp<-spec[[fig_name]]
 panels<-lapply(seq_along(sp$ref),function(i)getpanel(sp$ref[i],LETTERS[i]))
 fig<-switch(sp$layout,
  bulk=(panels[[1]]|panels[[2]]|panels[[3]])/
       panels[[4]]/
       (panels[[5]]|panels[[6]])+
       plot_layout(heights=c(.9,1.10,.88)),
  lim=(panels[[1]]|panels[[2]])/
      (panels[[3]]|panels[[4]])/
      (panels[[5]]|panels[[6]])/
      (panels[[7]]|panels[[8]])+
      plot_layout(heights=c(1,1,1,.9)),
  chromatin=(panels[[1]]|panels[[2]]|panels[[3]])/
       (panels[[4]]|panels[[5]]|panels[[6]])/
       (panels[[7]]|panels[[8]])+
       plot_layout(heights=c(.95,1,.95)))
 base<-file.path(out,fig_name)
 ggsave(paste0(base,'.pdf'),fig,width=sp$size[1],height=sp$size[2],
        units='in',device=cairo_pdf,bg='white',limitsize=FALSE)
 ggsave(paste0(base,'.svg'),fig,width=sp$size[1],height=sp$size[2],
        units='in',device=svg,bg='white',limitsize=FALSE)
 ggsave(paste0(base,'.png'),fig,width=sp$size[1],height=sp$size[2],
        units='in',device=agg_png,dpi=300,bg='white',limitsize=FALSE)
 ggsave(paste0(base,'.tiff'),fig,width=sp$size[1],height=sp$size[2],
        units='in',device=agg_tiff,dpi=300,bg='white',limitsize=FALSE)
 final_sd<-file.path(data_root,paste0('Main_',fig_name))
 dir.create(final_sd,recursive=TRUE,showWarnings=FALSE)
 manifest<-vector('list',length(sp$ref))
 for(i in seq_along(sp$ref)){
  part<-strsplit(sp$ref[i],':',fixed=TRUE)[[1]]
  k<-part[1];panel<-part[2]
  legacy<-file.path(data_root,k)
  m<-legacy_manifests[[k]]
  row<-m[m$panel==panel,,drop=FALSE]
  if(nrow(row)!=1) stop('Manifest row missing: ',sp$ref[i])
  old_names<-strsplit(row$displayed_data,' + ',fixed=TRUE)[[1]]
  new_names<-paste0('source_',k,'_',basename(old_names))
  for(j in seq_along(old_names)){
   src<-file.path(legacy,old_names[j]);dst<-file.path(final_sd,new_names[j])
   if(!file.copy(src,dst,overwrite=TRUE)) stop('Missing source data: ',src)
  }
  manifest[[i]]<-data.frame(figure=fig_name,panel=LETTERS[i],
    original_panel=sp$ref[i],source=row$source,
    displayed_data=paste(new_names,collapse=' + '),
    display_transform=row$display_transform,
    panel_script=row$script,
    assembly_script='total/rebuild_2026/code/build_main_validation_relayout.R')
 }
 manifest<-do.call(rbind,manifest)
 write.table(manifest,file.path(final_sd,paste0(fig_name,'_panel_manifest.tsv')),
             sep='\t',quote=FALSE,row.names=FALSE)
 write.table(manifest,paste0(base,'_panel_manifest.tsv'),
             sep='\t',quote=FALSE,row.names=FALSE)
 writeLines(capture.output(sessionInfo()),paste0(base,'_R_sessionInfo.txt'))
 cat('Built ',fig_name,' with ',length(panels),' frozen panels\n',sep='')
}

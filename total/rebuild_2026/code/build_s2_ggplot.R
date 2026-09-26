# Supplementary Figure S2: complete normal-thymus stage profiles and donor contrast.
suppressPackageStartupMessages({library(ggplot2);library(patchwork);library(ragg);library(readr)})
argv <- grep('^--file=', commandArgs(FALSE), value=TRUE)
script <- normalizePath(sub('^--file=', '', argv), winslash='/')
root <- normalizePath(file.path(dirname(script),'../../..'),winslash='/')
valid <- file.path(root,'total/data/validation')
out <- file.path(root,'total/rebuild_2026/edition_20260925/figures/supplementary')
sdir <- file.path(root,'total/rebuild_2026/data/source_data_rebuilt/S2')
dir.create(out,recursive=TRUE,showWarnings=FALSE);dir.create(sdir,recursive=TRUE,showWarnings=FALSE)
rd <- function(p) as.data.frame(read_tsv(file.path(valid,p),show_col_types=FALSE,
  progress=FALSE,locale=locale(encoding='UTF-8')))
wt <- function(x,n) write.table(x,file.path(sdir,n),sep='\t',quote=FALSE,row.names=FALSE,na='NA')
normal <- rd('zeb_developmental_axis/normal_thymus_balance_trajectory.tsv')
park <- as.data.frame(read_tsv(file.path(root,'total/data/SourceData_ParkHTA_donor_celltype.tsv'),
  show_col_types=FALSE,progress=FALSE,locale=locale(encoding='UTF-8')))
stopifnot(nrow(normal)==21)
ink <- '#152b39'; blue <- '#246b96'; red <- '#bc3446'; teal <- '#168b82';
gold <- '#c28830'; purple <- '#776298'; pale <- '#e7edf0'
ds <- c('GSE142522','GSE195812','GSE206710','Park/HTA')
normal$dataset <- factor(normal$dataset,levels=ds)
normal$stage_key <- paste(normal$dataset,sprintf('%02d',normal$order),normal$stage,sep='|')
normal$stage_key <- factor(normal$stage_key,levels=unique(normal$stage_key[order(normal$dataset,normal$order)]))
theme_pub <- theme_classic(base_family='Arial',base_size=12.1)+
 theme(plot.title=element_text(colour=ink,face='bold',size=13.3,margin=margin(b=5)),
       plot.subtitle=element_text(colour=ink,size=10.2),
       axis.title=element_text(colour=ink,size=11),
       axis.text=element_text(colour=ink,size=10),
       legend.text=element_text(colour=ink,size=9.7),
       strip.text=element_text(colour=ink,face='bold',size=10.2),
       plot.margin=margin(8,9,8,9),
       panel.grid.major.y=element_line(colour=pale,linewidth=.35))

# A: all observed stage-level balance values; deliberate discontinuities remain visible.
normal$colour <- ifelse(normal$balance>=0,'Cortical side','Early/later side')
pA <- ggplot(normal,aes(order,balance))+
 geom_hline(yintercept=0,colour='#aebec6',linewidth=.5)+
 geom_segment(aes(xend=order,y=0,yend=balance,colour=colour),linewidth=1.1)+
 geom_point(aes(colour=colour),size=2.8)+
 facet_wrap(~dataset,scales='free_x',nrow=1)+
 scale_colour_manual(values=c('Cortical side'=red,'Early/later side'=blue),guide='none')+
 scale_x_continuous(breaks=1:8)+
 labs(title='A  Every observed stage in four normal thymus datasets',
      subtitle='Frozen balance; stage order is defined within each dataset',
      x='Measured stage order',y='ZEB1 - ZEB2 balance')+theme_pub+
 theme(panel.grid.major.x=element_blank())

# B: gene-resolved stage matrix, standardized only within dataset and gene for display.
genes <- c('ZEB1','ZEB2','LMO2','LYL1')
g <- do.call(rbind,lapply(genes,function(gene) data.frame(dataset=normal$dataset,
  stage=normal$stage,order=normal$order,gene=gene,value=normal[[gene]])))
g$z <- ave(g$value,interaction(g$dataset,g$gene,drop=TRUE),FUN=function(x) if(all(is.na(x)))
  rep(NA_real_,length(x)) else as.numeric(scale(x)))
g$stage <- factor(paste(g$dataset,sprintf('%02d',g$order),g$stage,sep='|'),
                 levels=rev(levels(normal$stage_key)))
g$stage_label <- as.character(g$stage)
g$stage_label <- sub('^[^|]*\\|[0-9]+\\|','',g$stage_label)
g$gene <- factor(g$gene,levels=genes)
pB <- ggplot(g,aes(gene,stage,fill=z))+
 geom_tile(colour='white',linewidth=.7)+
 geom_text(aes(label=ifelse(is.na(z),'NA','')),colour=ink,size=2.8)+
 facet_grid(dataset~.,scales='free_y',space='free_y')+
 scale_y_discrete(labels=function(x) sub('^[^|]*\\|[0-9]+\\|','',x))+
 scale_fill_gradient2(low='#4a7fa0',mid='#f7f8f4',high='#bd5868',midpoint=0,
                      limits=c(-2,2),oob=scales::squish,na.value='#edf1f2',name='Within-atlas z')+
 labs(title='B  Full gene-by-stage map',
      subtitle='Each gene is standardized within its own atlas; NA was unmeasured',
      x=NULL,y=NULL)+theme_pub+
 theme(axis.line=element_blank(),axis.ticks=element_blank(),
       panel.grid=element_blank(),axis.text.y=element_text(size=8.7),
       strip.text.y=element_text(angle=0,size=9),legend.position='bottom',
       legend.key.width=grid::unit(13,'pt'))

# C: the GSE142522 CD8 single-positive ZEB2 rebound is explicit.
g142 <- normal[normal$dataset=='GSE142522',]
g142$stage <- factor(g142$stage,levels=g142$stage)
pC <- ggplot(g142,aes(stage,ZEB2))+
 geom_segment(aes(xend=stage,y=0,yend=ZEB2),linewidth=1.2,colour=pale)+
 geom_point(aes(colour=stage=='CD8 SP'),size=3.7)+
 geom_text(data=g142[g142$stage=='CD8 SP',],aes(label=sprintf('%.2f',ZEB2)),
           nudge_y=.55,colour=ink,size=3.3,family='Arial')+
 scale_colour_manual(values=c('FALSE'=blue,'TRUE'=red),guide='none')+
 scale_y_continuous(limits=c(0,11),expand=expansion(mult=c(0,.02)))+
 labs(title='C  CD8 SP ZEB2 recovery',subtitle='GSE142522 sorted bulk stages',
      x=NULL,y='Frozen ZEB2 expression')+theme_pub+
 theme(axis.text.x=element_text(angle=45,hjust=1,size=8.6),panel.grid.major.x=element_blank())

# D: the GSE195812 ISP local balance dip; no monotonic interpolation.
g195 <- normal[normal$dataset=='GSE195812',]
g195$stage <- factor(g195$stage,levels=g195$stage)
pD <- ggplot(g195,aes(stage,balance))+
 geom_hline(yintercept=0,colour='#aebec6',linewidth=.5)+
 geom_segment(aes(xend=stage,y=0,yend=balance),linewidth=1.25,colour=pale)+
 geom_point(aes(colour=stage=='ISP'),size=3.8)+
 geom_text(data=g195[g195$stage=='ISP',],aes(label=sprintf('%.2f',balance)),
           nudge_y=-.35,colour=ink,size=3.3,family='Arial')+
 scale_colour_manual(values=c('FALSE'=teal,'TRUE'=red),guide='none')+
 labs(title='D  ISP is a local reversal',subtitle='GSE195812 pooled FACS stage libraries',
      x=NULL,y='Frozen balance')+theme_pub+
 theme(axis.text.x=element_text(angle=45,hjust=1,size=8.6),panel.grid.major.x=element_blank())

# E: each line is a real eligible Park donor; no cell-level pseudoreplication.
don <- park[park$celltype %in% c('double negative thymocyte',
  'CD4-positive, alpha-beta T cell') & park$n_cells>=20,]
don$stage <- ifelse(don$celltype=='double negative thymocyte','DN','CD4 SP')
wide <- reshape(don[,c('donor','stage','ZEB1')],idvar='donor',timevar='stage',direction='wide')
wide <- wide[complete.cases(wide),]
don <- don[don$donor %in% wide$donor,]
don$stage <- factor(don$stage,levels=c('DN','CD4 SP'))
stopifnot(nrow(wide)==15)
pE <- ggplot(don,aes(stage,ZEB1,group=donor))+
 geom_line(colour='#94a9b2',linewidth=.65)+
 geom_point(aes(colour=stage),size=2.5)+
 scale_colour_manual(values=c('DN'=red,'CD4 SP'=blue),guide='none')+
 labs(title='E  Paired donor-level ZEB1 contrast',
      subtitle='Park/HTA: 15 donors with >=20 cells at both stages',
      x=NULL,y='Donor mean ZEB1')+theme_pub

fig <- pA / (pB | (pC / pD / pE)) + plot_layout(heights=c(.7,1.8),widths=c(1,1))
wt(normal,'S2A_complete_stage_balance.tsv');wt(g,'S2B_gene_by_stage_display_z.tsv')
wt(g142,'S2C_GSE142522_CD8_recovery.tsv');wt(g195,'S2D_GSE195812_ISP_dip.tsv')
wt(don[,c('donor','stage','ZEB1','n_cells')],'S2E_Park_paired_donors.tsv')
base <- file.path(out,'S2')
ggsave(paste0(base,'.pdf'),fig,width=13.6,height=12.0,units='in',device=cairo_pdf,bg='white')
ggsave(paste0(base,'.svg'),fig,width=13.6,height=12.0,units='in',device=svg,bg='white')
ggsave(paste0(base,'.png'),fig,width=13.6,height=12.0,units='in',device=agg_png,dpi=300,bg='white')
ggsave(paste0(base,'.tiff'),fig,width=13.6,height=12.0,units='in',device=agg_tiff,dpi=300,bg='white')
writeLines(capture.output(sessionInfo()),paste0(base,'_R_sessionInfo.txt'))
manifest <- data.frame(figure='S2',panel=LETTERS[1:5],
 source=c(rep('zeb_developmental_axis/normal_thymus_balance_trajectory.tsv',4),
          'total/data/SourceData_ParkHTA_donor_celltype.tsv'),
 displayed_data=c('S2A_complete_stage_balance.tsv','S2B_gene_by_stage_display_z.tsv',
                  'S2C_GSE142522_CD8_recovery.tsv','S2D_GSE195812_ISP_dip.tsv',
                  'S2E_Park_paired_donors.tsv'),
 display_transform=c('All 21 frozen stage entries','Within-atlas per-gene z for display only',
                     'Frozen six stage ZEB2 values','Frozen eight stage balance values',
                     'DN and CD4 SP matched donors with >=20 cells per stage'),
 script='total/rebuild_2026/code/build_s2_ggplot.R')
wt(manifest,'S2_panel_manifest.tsv');cat(base,'\n')

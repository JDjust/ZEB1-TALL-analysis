# Supplementary Figure S1: actual source units, cohort composition and evidence scope.
suppressPackageStartupMessages({library(ggplot2); library(patchwork); library(ragg); library(readr)})
argv <- grep('^--file=', commandArgs(FALSE), value=TRUE)
script <- normalizePath(sub('^--file=', '', argv), winslash='/')
root <- normalizePath(file.path(dirname(script), '../../..'), winslash='/')
valid <- file.path(root, 'total/data/validation')
out <- file.path(root, 'total/rebuild_2026/edition_20260925/figures/supplementary')
sdir <- file.path(root, 'total/rebuild_2026/data/source_data_rebuilt/S1')
dir.create(out, recursive=TRUE, showWarnings=FALSE)
dir.create(sdir, recursive=TRUE, showWarnings=FALSE)
rd <- function(p) as.data.frame(read_tsv(file.path(valid,p), show_col_types=FALSE,
                                           progress=FALSE, locale=locale(encoding='UTF-8')))
wt <- function(x,n) write.table(x,file.path(sdir,n),sep='\t',quote=FALSE,row.names=FALSE,na='NA')
normal <- rd('zeb_developmental_axis/normal_thymus_balance_trajectory.tsv')
audit <- rd('polonen_round1/sample_audit.tsv')
bulk <- rd('gse146901/zeb_expression_by_sample.tsv')
bcl <- rd('gse162280/case_manifest_and_expression.tsv')
lim <- rd('gateA_lim2025/q3_state_day0_pairs.tsv')
hichip <- rd('zeb_chromatin_reciprocal/hichip_zeb1_zeb2_contacts.tsv')
atac <- rd('zeb_chromatin_reciprocal/scatac_zeb1_zeb2_peaks.tsv')
stopifnot(nrow(normal)==21,nrow(bulk)==22,nrow(bcl)==12,nrow(lim)==41,
          nrow(hichip)==10,nrow(atac)==7)
ink <- '#152b39'; teal <- '#168b82'; blue <- '#246b96'; red <- '#bc3446';
purple <- '#776298'; gold <- '#c28830'; pale <- '#e6ecee'
theme_pub <- theme_classic(base_family='Arial',base_size=12.5)+
 theme(plot.title=element_text(colour=ink,face='bold',size=14,margin=margin(b=5)),
       plot.subtitle=element_text(colour=ink,size=10.5),
       plot.caption=element_text(colour=ink,size=9.5),
       axis.title=element_text(colour=ink,size=11.5),
       axis.text=element_text(colour=ink,size=10.5),
       legend.text=element_text(colour=ink,size=10),
       legend.title=element_text(colour=ink,size=10),
       plot.margin=margin(8,9,8,9),
       panel.grid.major.y=element_line(colour=pale,linewidth=.35))

# A. Stage units are not independent donors; place sample unit alongside count.
normal_n <- aggregate(stage ~ dataset,normal,function(x) length(unique(x)))
names(normal_n)[2] <- 'n_stage_units'
unit_map <- c('GSE142522'='sorted bulk stage','GSE195812'='pooled FACS stage',
              'Park/HTA'='donor-stage aggregate','GSE206710'='donor-stage aggregate')
normal_n$unit <- unname(unit_map[normal_n$dataset])
normal_n$dataset <- factor(normal_n$dataset,levels=rev(normal_n$dataset[order(normal_n$n_stage_units)]))
pA <- ggplot(normal_n,aes(n_stage_units,dataset))+
 geom_segment(aes(x=0,xend=n_stage_units,yend=dataset),colour=pale,linewidth=2)+
 geom_point(size=4.2,colour=blue)+
 geom_text(aes(label=paste0(n_stage_units,'  ',unit)),hjust=-.08,size=3.45,
           colour=ink,family='Arial')+
 scale_x_continuous(limits=c(0,20),breaks=c(0,5,10,15),expand=c(0,0))+
 labs(title='A  Normal atlases: stage labels and sampling units',
      subtitle='Stage units are plotted as provided by each frozen atlas',
      x='Distinct stage labels',y=NULL)+theme_pub+
 theme(panel.grid.major.y=element_blank())

# B. Each glyph has an explicit biological/sample unit; no false cross-study pooling.
cohorts <- data.frame(dataset=c('Polonen diagnostic','Lim Day 0 paired states',
  'GSE146901 bulk RNA','GSE162280 BCL11B-R','HiChIP','scATAC'),
  n=c(1309,nrow(lim),nrow(bulk),nrow(bcl),nrow(hichip),nrow(atac)),
  unit=c('patients','paired patients','RNA samples','BCL11B-R cases',
         'contact libraries','peak sets'),
  scope=c('Discovery','Cell-state localization','Expression replication',
          'Within-lesion convergence','Orthogonal support','Negative sensitivity'))
cohorts$dataset <- factor(cohorts$dataset,levels=rev(cohorts$dataset))
cohorts$label <- paste0('n=',cohorts$n,'  ',cohorts$unit)
pB <- ggplot(cohorts,aes(n,dataset,colour=scope))+
 geom_segment(aes(x=1,xend=n,yend=dataset),colour=pale,linewidth=1.7)+
 geom_point(size=4.1)+
 geom_text(aes(label=label),hjust=-.08,colour=ink,size=3.45,family='Arial')+
 scale_x_log10(limits=c(1,6000),breaks=c(1,10,100,1000),labels=scales::label_comma())+
 scale_colour_manual(values=c('Discovery'=blue,'Cell-state localization'=teal,
    'Expression replication'=red,'Within-lesion convergence'=purple,
    'Orthogonal support'=gold,'Negative sensitivity'='#81919a'),guide='none')+
 labs(title='B  Discovery and validation units',
      subtitle='Counts are descriptive; axes are logarithmic',
      x='Number of independent displayed units',y=NULL)+theme_pub+
 theme(panel.grid.major.y=element_blank())

# C. Actual subdivisions of four validation resources.
comp <- rbind(
  data.frame(source='GSE146901',group=as.character(bulk$group)),
  data.frame(source='GSE162280',group=ifelse(bcl$zeb2_partner=='yes','ZEB2 fusion','Other partner')),
  data.frame(source='HiChIP',group=as.character(hichip$group)),
  data.frame(source='scATAC',group=as.character(atac$group)))
comp <- as.data.frame(table(comp$source,comp$group),stringsAsFactors=FALSE)
names(comp) <- c('source','group','n');comp <- comp[comp$n>0,]
comp$source <- factor(comp$source,levels=c('GSE146901','GSE162280','HiChIP','scATAC'))
comp$colour <- ifelse(grepl('ETP',comp$group,ignore.case=TRUE)&
  !grepl('non',comp$group,ignore.case=TRUE),teal,
  ifelse(grepl('fusion',comp$group,ignore.case=TRUE),purple,
  ifelse(grepl('CD34|thymus|normal',comp$group,ignore.case=TRUE),gold,blue)))
pC <- ggplot(comp,aes(n,reorder(group,n),fill=colour))+
 geom_col(width=.68)+geom_text(aes(label=n),hjust=-.15,size=3.2,colour=ink,family='Arial')+
 facet_wrap(~source,scales='free_y',ncol=2)+
 scale_fill_identity()+scale_x_continuous(expand=expansion(mult=c(0,.18)))+
 labs(title='C  Validation series: source composition',
      subtitle='Distinct case/library groups; no pooled sample is counted as a patient',
      x='Samples or cases',y=NULL)+theme_pub+
 theme(strip.text=element_text(colour=ink,face='bold',size=10.5),
       axis.text.y=element_text(size=9.4),panel.grid.major.y=element_blank())

# D. The evidence matrix records direct support versus an explicit sensitivity/limit.
claims <- c('Normal stage pattern','17-subtype landscape','Residual reconfiguration',
            'ETP expression polarity','BCL11B convergence','Chromatin organization')
sources <- c('Normal thymus','Polonen','Lim','GSE146901','GSE162280','HiChIP','scATAC')
ev <- expand.grid(source=sources,claim=claims,stringsAsFactors=FALSE)
ev$role <- 'Not tested'
set_role <- function(claim,src,role) ev$role[ev$claim==claim & ev$source %in% src] <<- role
set_role('Normal stage pattern','Normal thymus','Direct')
set_role('17-subtype landscape','Polonen','Direct')
set_role('Residual reconfiguration',c('Normal thymus','Polonen'),'Direct')
set_role('ETP expression polarity',c('Polonen','GSE146901'),'Direct')
set_role('ETP expression polarity','Lim','Context')
set_role('BCL11B convergence',c('Polonen','GSE162280'),'Direct')
set_role('Chromatin organization','HiChIP','Direct')
set_role('Chromatin organization',c('GSE146901','scATAC'),'Limit')
ev$source <- factor(ev$source,levels=sources)
ev$claim <- factor(ev$claim,levels=rev(claims))
ev$glyph <- c('Not tested'='', 'Context'='C','Direct'='D','Limit'='L')[ev$role]
pD <- ggplot(ev,aes(source,claim,fill=role))+
 geom_tile(colour='white',linewidth=.9)+
 geom_text(aes(label=glyph),colour=ink,size=3.5,fontface='bold',family='Arial')+
 scale_fill_manual(values=c('Direct'='#7bb4b0','Context'='#b9cbd3',
                            'Limit'='#e9c9ad','Not tested'='#f4f6f6'),
                   breaks=c('Direct','Context','Limit'),name='Evidence role')+
 labs(title='D  Claim-by-source evidence scope',
      subtitle='D direct; C context; L measured limitation',x=NULL,y=NULL)+theme_pub+
 theme(axis.text.x=element_text(angle=35,hjust=1,size=9.2),
       axis.line=element_blank(),axis.ticks=element_blank(),
       panel.grid=element_blank(),legend.position='bottom')

fig <- (pA | pB) / (pC | pD) + plot_layout(heights=c(.9,1.15))
wt(normal_n,'S1A_normal_stage_units.tsv');wt(cohorts,'S1B_source_units.tsv')
wt(comp,'S1C_validation_composition.tsv');wt(ev,'S1D_evidence_scope.tsv')
base <- file.path(out,'S1')
ggsave(paste0(base,'.pdf'),fig,width=13.6,height=9.8,units='in',device=cairo_pdf,bg='white')
ggsave(paste0(base,'.svg'),fig,width=13.6,height=9.8,units='in',device=svg,bg='white')
ggsave(paste0(base,'.png'),fig,width=13.6,height=9.8,units='in',device=agg_png,dpi=300,bg='white')
ggsave(paste0(base,'.tiff'),fig,width=13.6,height=9.8,units='in',device=agg_tiff,dpi=300,bg='white')
writeLines(capture.output(sessionInfo()),paste0(base,'_R_sessionInfo.txt'))
manifest <- data.frame(figure='S1',panel=LETTERS[1:4],
 source=c('zeb_developmental_axis/normal_thymus_balance_trajectory.tsv',
          'polonen_round1/sample_audit.tsv + gateA_lim2025/q3_state_day0_pairs.tsv + gse146901/zeb_expression_by_sample.tsv + gse162280/case_manifest_and_expression.tsv + zeb_chromatin_reciprocal/hichip_zeb1_zeb2_contacts.tsv + zeb_chromatin_reciprocal/scatac_zeb1_zeb2_peaks.tsv',
          'gse146901/zeb_expression_by_sample.tsv + gse162280/case_manifest_and_expression.tsv + zeb_chromatin_reciprocal/hichip_zeb1_zeb2_contacts.tsv + zeb_chromatin_reciprocal/scatac_zeb1_zeb2_peaks.tsv',
          'Frozen study design and negative-validation interpretation'),
 displayed_data=c('S1A_normal_stage_units.tsv','S1B_source_units.tsv',
                  'S1C_validation_composition.tsv','S1D_evidence_scope.tsv'),
 display_transform=c('Count distinct frozen stage labels per atlas',
                     'Display frozen source units on log axis',
                     'Count observed rows by frozen group',
                     'Editorial evidence-role mapping, not statistical analysis'),
 script='total/rebuild_2026/code/build_s1_ggplot.R')
wt(manifest,'S1_panel_manifest.tsv')
cat(base,'\n')

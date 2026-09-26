# Supplementary Figure S6. Forest grammar adapted from code1/38森林图/bioR38.forest.R.
suppressPackageStartupMessages({library(ggplot2);library(patchwork);library(ragg);library(readr)})
argv<-grep('^--file=',commandArgs(FALSE),value=TRUE)
script<-normalizePath(sub('^--file=','',argv),winslash='/')
root<-normalizePath(file.path(dirname(script),'../../..'),winslash='/')
v<-file.path(root,'total/data/validation');out<-file.path(root,'total/rebuild_2026/edition_20260925/figures/supplementary')
sdir<-file.path(root,'total/rebuild_2026/data/source_data_rebuilt/S6')
dir.create(out,recursive=TRUE,showWarnings=FALSE);dir.create(sdir,recursive=TRUE,showWarnings=FALSE)
rd<-function(p) as.data.frame(read_tsv(file.path(v,p),show_col_types=FALSE,
 progress=FALSE,locale=locale(encoding='UTF-8')))
wt<-function(x,n) write.table(x,file.path(sdir,n),sep='\t',quote=FALSE,row.names=FALSE,na='NA')
IF<-rd('polonen_round1/model_induction_failure.tsv')
M<-rd('polonen_round1/model_poor_morphology.tsv')
MRD<-rd('polonen_round1/model_mrd.tsv')
S<-rd('polonen_round1/model_survival.tsv')
evt<-rd('polonen_round1/events_by_subtype.tsv')
clean<-function(x){x<-as.character(x);x[grepl('^LMO2.*-like$',x)]<-'LMO2 gamma-delta-like';
 x[grepl('^TAL1.*-like$',x)&!grepl('DP',x)]<-'TAL1 alpha-beta-like';x}
evt$subtype<-clean(evt$subtype)
mk<-function(x,endpoint,rows,measure='OR'){
 d<-x[rows,];d$endpoint<-endpoint;d$estimate<-if(measure=='HR') d$hr else d$or
 d$measure<-measure;d$adjustment<-c('Unadjusted','Subtype-adjusted');
 d[,c('endpoint','adjustment','n','events','estimate','ci_low','ci_high','p','measure')]
}
clin<-rbind(mk(IF,'Induction failure',1:2),mk(M,'M2/M3 morphology',1:2),
 mk(MRD,'MRD >=0.1%',1:2),mk(MRD,'MRD >=0.01%',3:4),
 mk(S,'Event-free survival',1:2,'HR'),mk(S,'Overall survival',3:4,'HR'))
clin$endpoint<-factor(clin$endpoint,levels=rev(c('Induction failure','M2/M3 morphology',
 'MRD >=0.1%','MRD >=0.01%','Event-free survival','Overall survival')))
clin$adjustment<-factor(clin$adjustment,levels=c('Unadjusted','Subtype-adjusted'))
stopifnot(nrow(clin)==12,all(is.finite(clin$estimate)),nrow(evt)==17)
ink<-'#152b39';red<-'#bc3446';teal<-'#168b82';blue<-'#246b96';pale<-'#e5ecef'
theme_pub<-theme_classic(base_family='Arial',base_size=12.2)+theme(
 plot.title=element_text(colour=ink,face='bold',size=13.6,margin=margin(b=5)),
 plot.subtitle=element_text(colour=ink,size=10.2),
 axis.title=element_text(colour=ink,size=11),axis.text=element_text(colour=ink,size=10),
 legend.text=element_text(colour=ink,size=9.5),plot.margin=margin(8,9,8,9),
 panel.grid.major.y=element_line(colour=pale,linewidth=.35))

# The code1 forest interval grammar (reference line, CI segment and point) is
# retained; the hard-coded toy HR values are replaced by six frozen endpoints.
pA<-ggplot(clin,aes(estimate,endpoint,colour=adjustment))+
 geom_vline(xintercept=1,colour='#8ea4ae',linetype=2,linewidth=.55)+
 geom_segment(aes(x=ci_low,xend=ci_high,yend=endpoint),
              position=position_dodge(width=.52),linewidth=.8)+
 geom_point(position=position_dodge(width=.52),size=2.8)+
 scale_x_log10(limits=c(.46,1.32),breaks=c(.5,.75,1,1.25),
               labels=scales::number_format(accuracy=.01))+
 scale_colour_manual(values=c('Unadjusted'=red,'Subtype-adjusted'=teal),name=NULL)+
 labs(title='A  All frozen clinical endpoints',
      subtitle='OR for binary endpoints; HR for survival; 95% confidence intervals',
      x='Effect per unit higher ZEB balance (ratio scale)',y=NULL)+theme_pub+
 theme(legend.position='bottom')

clin$adj_key<-ifelse(clin$adjustment=='Unadjusted','raw','adjusted')
wide<-reshape(clin[,c('endpoint','adj_key','estimate','p','measure','n','events')],
 idvar='endpoint',timevar='adj_key',direction='wide')
wide$endpoint<-factor(wide$endpoint,levels=levels(clin$endpoint))
pB<-ggplot(wide,aes(y=endpoint))+
 geom_vline(xintercept=0,colour='#8ea4ae',linetype=2,linewidth=.55)+
 geom_segment(aes(x=log(estimate.raw),xend=log(estimate.adjusted),
                  yend=endpoint),linewidth=1.25,colour='#8aa5ad',
              arrow=grid::arrow(length=grid::unit(3,'pt'),type='closed'))+
 geom_point(aes(x=log(estimate.raw)),colour=red,size=3)+
 geom_point(aes(x=log(estimate.adjusted)),colour=teal,size=3)+
 scale_x_continuous(breaks=log(c(.5,.75,1,1.25)),labels=c('.50','.75','1.00','1.25'))+
 labs(title='B  Attenuation after subtype adjustment',
      subtitle='Arrow points from unadjusted to subtype-adjusted effect',
      x='Effect ratio (displayed on log scale)',y=NULL)+theme_pub

# Frozen subtype events are displayed as observed rates, not new outcome models.
evt$if_rate<-evt$induction_failure/evt$n
evt$mrd_rate<-evt$mrd_ge_0.1/evt$mrd_known
long<-rbind(data.frame(subtype=evt$subtype,endpoint='Induction failure',rate=evt$if_rate,
                       events=evt$induction_failure,denominator=evt$n,median_balance=evt$median_balance),
            data.frame(subtype=evt$subtype,endpoint='MRD >=0.1%',rate=evt$mrd_rate,
                       events=evt$mrd_ge_0.1,denominator=evt$mrd_known,median_balance=evt$median_balance))
order<-evt$subtype[order(evt$median_balance)]
long$subtype<-factor(long$subtype,levels=rev(order))
long$endpoint<-factor(long$endpoint,levels=c('Induction failure','MRD >=0.1%'))
long$label<-paste0(long$events,'/',long$denominator)
pC<-ggplot(long,aes(endpoint,subtype,fill=rate))+
 geom_tile(colour='white',linewidth=.8)+
 geom_text(aes(label=label),colour=ink,size=3.05,family='Arial')+
 scale_fill_gradient(low='#f6f8f6',high='#cb7783',limits=c(0,.85),
                     oob=scales::squish,labels=scales::percent_format(),name='Observed rate')+
 labs(title='C  Endpoint composition across molecular subtypes',
      subtitle='Events/known denominators printed; subtypes ordered by median balance',
      x=NULL,y=NULL)+theme_pub+
 theme(panel.grid=element_blank(),axis.line=element_blank(),axis.ticks=element_blank(),
       axis.text.y=element_text(size=9.1),legend.position='right')

fig<-(pA|pB)/pC+plot_layout(heights=c(.93,1.12))
wt(clin,'S6A_all_endpoint_models.tsv');wt(wide,'S6B_attenuation_pairs.tsv')
wt(long,'S6C_events_by_subtype.tsv')
base<-file.path(out,'S6')
ggsave(paste0(base,'.pdf'),fig,width=13.8,height=11.0,units='in',device=cairo_pdf,bg='white')
ggsave(paste0(base,'.svg'),fig,width=13.8,height=11.0,units='in',device=svg,bg='white')
ggsave(paste0(base,'.png'),fig,width=13.8,height=11.0,units='in',device=agg_png,dpi=300,bg='white')
ggsave(paste0(base,'.tiff'),fig,width=13.8,height=11.0,units='in',device=agg_tiff,dpi=300,bg='white')
writeLines(capture.output(sessionInfo()),paste0(base,'_R_sessionInfo.txt'))
manifest<-data.frame(figure='S6',panel=LETTERS[1:3],
 source=c('polonen_round1/model_induction_failure.tsv + model_poor_morphology.tsv + model_mrd.tsv + model_survival.tsv',
          'Same frozen model tables as S6A','polonen_round1/events_by_subtype.tsv'),
 displayed_data=c('S6A_all_endpoint_models.tsv','S6B_attenuation_pairs.tsv',
                  'S6C_events_by_subtype.tsv'),
 display_transform=c('code1 forest grammar with frozen OR/HR and intervals',
                     'Paired log-effect positions; no refitting',
                     'Frozen subtype events divided by frozen known denominators'),
 script='total/rebuild_2026/code/build_s6_ggplot.R')
wt(manifest,'S6_panel_manifest.tsv');cat(base,'\n')

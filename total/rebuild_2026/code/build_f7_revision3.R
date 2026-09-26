# Figure 7: post-review robustness and evidence-boundary synthesis.
suppressPackageStartupMessages({library(ggplot2);library(patchwork);library(data.table)})
argv <- grep("^--file=", commandArgs(FALSE), value=TRUE)
script <- normalizePath(sub("^--file=", "", argv), winslash="/")
root <- normalizePath(file.path(dirname(script), "../../.."), winslash="/")
src <- file.path(root,"total/rebuild_2026/data/source_data_rebuilt/revision3_robustness")
rd <- function(x) fread(file.path(src,x),encoding="UTF-8")
cf <- rd("crossfit_subtype_summary.tsv")
alt <- rd("alternative_coordinates_subtype_summary.tsv")
mod <- rd("model_comparison.tsv")
pair <- rd("bcl11b_etp_matched_pairs.tsv")
ms <- rd("bcl11b_etp_matched_summary.tsv")
con <- rd("bcl11b_pairwise_bootstrap.tsv")
ink <- "#1c3645"; grey <- "#a6b7be"; purple <- "#76519c"
teal <- "#16877f"; orange <- "#c57732"; blue <- "#467899"
base <- theme_classic(base_family="Arial",base_size=9) +
  theme(plot.title=element_text(face="bold",size=10.5,colour=ink,margin=margin(b=4)),
        plot.subtitle=element_text(size=8.4,colour="#526873"),
        axis.title=element_text(size=8.7,colour=ink),
        axis.text=element_text(size=8.2,colour=ink),
        panel.grid.major.y=element_line(colour="#e7edef",linewidth=.25),
        plot.margin=margin(4,5,4,5),legend.position="bottom",
        legend.title=element_blank(),legend.text=element_text(size=8))

cf[,label:=ifelse(grepl("^LMO2",subtype),"LMO2 gamma-delta-like",
             ifelse(grepl("^TAL1",subtype)&grepl("like$",subtype)&subtype!="TAL1 DP-like",
                    "TAL1 alpha-beta-like",subtype))]
cf[,label:=factor(label,levels=rev(label[order(median_crossfit)]))]
cf[,focus:=ifelse(subtype=="BCL11B","BCL11B",
             ifelse(subtype=="TLX3","TLX3","Other subtypes"))]
pA <- ggplot(cf,aes(y=label))+
  geom_vline(xintercept=0,colour="#8fa3ab",linewidth=.35)+
  geom_segment(aes(x=median_primary,xend=median_crossfit,yend=label),
               colour="#b7c6cb",linewidth=.65)+
  geom_point(aes(x=median_primary),colour=grey,size=1.75)+
  geom_point(aes(x=median_crossfit,colour=focus),size=2.15)+
  scale_colour_manual(values=c(BCL11B=purple,TLX3=orange,
                                "Other subtypes"=ink),breaks=c("BCL11B","TLX3","Other subtypes"))+
  labs(title="A  Leave-one-subtype-out residuals",
       subtitle="Each subtype predicted from the other 16; grey = original cohort fit",
       x="Median observed - expected balance",y=NULL)+base+
  theme(axis.text.y=element_text(size=7.5),legend.position="none")

sel <- alt[subtype %in% c("BCL11B","ETP-like","TLX3")]
sel[,definition:=factor(definition,levels=rev(c("primary","minus_CD34","minus_LYL1",
                                              "minus_CD1A","expanded")),
  labels=rev(c("3-gene primary","without CD34","without LYL1","without CD1A",
               "expanded marker set")))]
sel[,subtype:=factor(subtype,levels=c("BCL11B","ETP-like","TLX3"))]
pB <- ggplot(sel,aes(median_residual,definition,colour=subtype))+
  geom_vline(xintercept=0,colour="#a3b3b9",linewidth=.35)+
  geom_point(size=2.45,position=position_dodge(width=.45))+
  scale_colour_manual(values=c(BCL11B=purple,"ETP-like"=teal,TLX3=orange))+
  labs(title="B  Marker-definition sensitivity",
       subtitle="BCL11B and TLX3 rank 1st or 2nd",
       x="Median residual",y=NULL)+base+
  theme(axis.text.y=element_text(size=7.7))

mm <- melt(mod[,.(model,adjusted_r2,cv_r2_mean)],id.vars="model",
           variable.name="metric",value.name="r2")
mm[,model:=factor(model,levels=rev(c("immunophenotype","development","subtype","combined")),
  labels=rev(c("Immunophenotype","Development","17 subtypes","Combined")))]
mm[,metric:=factor(metric,levels=c("adjusted_r2","cv_r2_mean"),
  labels=c("Adjusted R-squared","5-fold CV R-squared"))]
pC <- ggplot(mm,aes(r2,model,colour=metric))+
  geom_point(position=position_dodge(width=.35),size=2.6)+
  scale_colour_manual(values=c("Adjusted R-squared"=blue,"5-fold CV R-squared"=orange))+
  scale_x_continuous(limits=c(0,.43),breaks=c(0,.1,.2,.3,.4))+
  labs(title="C  Fair model comparison",
       subtitle="Repeated stratified five-fold CV",
       x="Explained variance",y=NULL)+base+
  theme(axis.text.y=element_text(size=7.8))

long <- rbindlist(list(pair[,.(id=.I,group="BCL11B",residual=BCL11B_residual)],
                       pair[,.(id=.I,group="Matched ETP-like",residual=ETP_residual)]))
long[,group:=factor(group,levels=c("BCL11B","Matched ETP-like"))]
pD <- ggplot(long,aes(group,residual,group=id))+
  geom_hline(yintercept=0,colour="#98abb2",linewidth=.35)+
  geom_line(colour="#b8c7cb",alpha=.75,linewidth=.45)+
  geom_point(aes(colour=group),size=1.65,alpha=.9)+
  scale_colour_manual(values=c(BCL11B=purple,"Matched ETP-like"=teal),guide="none")+
  labs(title="D  Development-matched ETP contrast",
       subtitle=sprintf("18 pairs; mean residual difference %.2f (95%% CI %.2f to %.2f)",
                        ms$mean_residual_difference,ms$ci_low,ms$ci_high),
       x=NULL,y="Primary residual")+base+
  theme(axis.text.x=element_text(size=8))

con[,comparison:=gsub("BCL11B vs ","",comparison,fixed=TRUE)]
con[,comparison:=ifelse(grepl("^LMO2",comparison),"LMO2 gamma-delta-like",comparison)]
con[,comparison:=factor(comparison,levels=rev(c("ETP-like","SPI1","LMO2 gamma-delta-like")))]
pE <- ggplot(con,aes(y=comparison,x=median_difference))+
  geom_vline(xintercept=0,colour="#8fa3ab",linewidth=.4)+
  geom_segment(aes(x=ci_low,xend=ci_high,yend=comparison),linewidth=.8,colour=ink)+
  geom_point(size=2.8,colour=purple)+
  scale_x_continuous(limits=c(-2.65,1.6),breaks=c(-2,-1,0,1))+
  labs(title="E  Fixed pairwise contrasts",
       subtitle="Bootstrap 95% CI",
       x="Difference in median primary residual",y=NULL)+base+
  theme(axis.text.y=element_text(size=7.7))

rankdat <- copy(alt)
rankdat[,definition:=factor(definition,
  levels=rev(c("primary","minus_CD34","minus_LYL1","minus_CD1A","expanded")),
  labels=rev(c("3-gene primary","without CD34","without LYL1",
               "without CD1A","expanded marker set")))]
rankdat[,focus:=fifelse(subtype=="BCL11B","BCL11B",
                 fifelse(subtype=="ETP-like","ETP-like",
                 fifelse(subtype=="TLX3","TLX3","Other subtypes")))]
pF <- ggplot(rankdat,aes(rank_negative,definition))+
  geom_point(data=rankdat[focus=="Other subtypes"],
             colour="#c8d3d7",size=1.45)+
  geom_point(data=rankdat[focus!="Other subtypes"],
             aes(colour=focus),size=2.8)+
  scale_colour_manual(values=c(BCL11B=purple,"ETP-like"=teal,TLX3=orange),
                      breaks=c("BCL11B","ETP-like","TLX3"))+
  scale_x_continuous(breaks=c(1,5,9,13,17),limits=c(.5,17.5))+
  labs(title="F  Subtype ranks across five marker definitions",
       subtitle="Each dot is one of 17 subtypes; ranks use median residual",
       x="Rank (1 = most ZEB2-skewed; 17 = most ZEB1-skewed)",y=NULL)+base+
  theme(axis.text.y=element_text(size=7.6),
        panel.grid.major.x=element_line(colour="#eef1f2",linewidth=.22),
        legend.position="top",legend.justification="right",
        legend.key.width=grid::unit(5,"mm"),legend.margin=margin(0,0,0,0))

F7 <- pA / (pB|pC) / (pD|pE) / pF +
  plot_layout(heights=c(2.7,1.75,1.65,1.25))
out <- file.path(root,"total/rebuild_2026/edition_20260925/figures/main")
dir.create(out,recursive=TRUE,showWarnings=FALSE)
ggsave(file.path(out,"F7.pdf"),F7,width=7.2,height=8.5,device=cairo_pdf,bg="white")
ggsave(file.path(out,"F7.png"),F7,width=7.2,height=8.5,dpi=250,bg="white")

library(data.table)
library(ggplot2)
library(cowplot)

## load in iCHIP GUESSFM datasets

DISEASES<-c('T1D','ATD','RA','CEL')

gfm<-fread("~/DATA/JAVIERRE_ICHIP/out/GUESSFM_FINAL/gene_prioritisation_final_fix_june_0.01.csv")
gfm<-subset(gfm, disease %in% DISEASES)

## load in PMI datasets
data.dir<-'~/DATA/JAVIERRE_ICHIP/out/tnact_hierarchical_geneScore/'
suffix<-'_IC.pmi_full.tab'

pmi<-do.call("rbind",lapply(DISEASES,function(d){
  t<-fread(file.path(data.dir,paste0(d,suffix)))
}))

pmi$disease<-sub("\\_IC","",pmi$disease)

## next we merge the scores

gfm<-gfm[,.(ensg,overall_ppi,disease)]
gfm$k<-with(gfm,paste(ensg,disease,sep=":"))
gfm<-gfm[,.(overall_ppi,k)]
setkey(gfm,k)
pmi<-pmi[,.(disease,ensg,name,biotype,overall_gene_score)]
pmi$k<-with(pmi,paste(ensg,disease,sep=":"))
setkey(pmi,k)

m<-pmi[gfm]

## why are genes missing from the PMI dataset ?

m<-subset(m,!is.na(disease) & biotype=='protein_coding')

pdf("~/gitr/fyr/pics/gfm_vs_pmi_gs.pdf")
ggplot(m,aes(x=overall_gene_score,y=overall_ppi)) + geom_point() + facet_wrap(~disease) + theme_bw() + geom_abline(slope=1,intercept=0,alpha=0.5) +
  xlab("PMI Gene Score") + ylab("GUESSFM Gene Score") + geom_smooth(method = "glm", se = FALSE)
dev.off()

## statistics of those that are called in PMI vs GUESSFM

thresh<-0.5

library(xtable)
c<-rbindlist(lapply(split(m,m$disease),function(d){
 both<-nrow(subset(d,overall_gene_score>0.5 & overall_ppi>0.5))
 gfm<-nrow(subset(d,overall_gene_score<0.5 & overall_ppi>0.5))
 pmi<-nrow(subset(d,overall_gene_score>0.5 & overall_ppi<0.5))
 data.table(disease=unique(d$disease),gfm=gfm,pmi=pmi,both=both,all=nrow(d))
}))

xtable(c)
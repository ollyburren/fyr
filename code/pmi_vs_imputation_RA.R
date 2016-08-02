# Program to compare PMI to classical imputation using RA Okada et al dataset

library(data.table)

GWAS.DIR<-'/Users/oliver/DATA/PMI_COMPARISON/OKADA/'
CHROMOSOME<-'1'
## Eyre dataset
PMI.FILE.RA<-file.path(GWAS.DIR,'combined.bed')
res<-fread(PMI.FILE.RA,header=FALSE)
setnames(res,c('chr','start','stop','info','p.pmi','chra','starta','startb','name','p.imp'))
res<-res[,.(chr,start,stop,info,p.pmi,name,p.imp)]
res.imp<-res[-grep("NA",res$info),]
res.imp<-subset(res.imp,p.imp>0)
## what is the correlation coefficient
with(res.imp,cor(p.pmi,p.imp))

library(ggplot2)
library(cowplot)

pdf("/Users/oliver/gitr/fyr/pics/pmi_vs_imp.pdf")
ggplot(res.imp,aes(y=-log10(p.pmi),x=-log10(p.imp))) + 
  geom_hex(bins=100) + theme_bw() + xlab("-log10(p.imputed)") + geom_abline(alpha=0.5)
dev.off()
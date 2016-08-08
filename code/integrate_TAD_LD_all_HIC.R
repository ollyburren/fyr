library(data.table)
library(GenomicRanges)
library(parallel)

processTAD<-function(tad.f,ld.gr){
  tad<-fread(tad.f)
  setnames(ld,c('chr','start','end','id'))
  tads.gr<-with(tad,GRanges(seqnames=Rle(chr),ranges=IRanges(start=start,end=end)))
  tads.gr$tid<-1:length(tads.gr)
  m<-as.data.table(mergeByOverlaps(ld.gr,tads.gr))
  ## cut up TADs based on LD - slow code but appears to work change mc.cores for when not in RStudio
  tad.ld<-mclapply(split(m,m$tid),function(x){
    print(paste("processing",unique(x$tid)))
    n<-nrow(x)
    if(n>1){
      x[1]$tads.gr.end<-x[1]$ld.gr.end
      x[n]$tads.gr.start<-x[n]$ld.gr.start
      ##rest just match the ld blocks
      if(n>2){
        s<-2:(n-1)
        x[s]$tads.gr.end<-x[s]$ld.gr.end
        x[s]$tads.gr.start<-x[s]$ld.gr.start
      }
    }
    return(x[,.(tads.gr.seqnames,tads.gr.start,tads.gr.end,tid)])
  },mc.cores=1)
  tad.ld<-rbindlist(tad.ld)
  setnames(tad.ld,c('chr','start','end','tid'))
  tad.ld.gr<-with(tad.ld,GRanges(seqnames=Rle(chr),ranges=IRanges(start=start,end=end),tid=tid))
  tad.ld.gr$tid.ld<-paste(tad.ld.gr$tid,1:length(tad.ld.gr),sep='.')
  tad.ld.gr
}

annotateTAD<-function(tad.f,h3.gr){
  tad<-fread(tad.f)
  tads.gr<-with(tad,GRanges(seqnames=Rle(chr),ranges=IRanges(start=start,end=end)))
  tads.gr$tid<-1:length(tads.gr)
  h3.tad<-as.data.table(mergeByOverlaps(h3.gr,tads.gr))
  h3.tad<-h3.tad[,.(info,tid,id)]
  setkey(h3.tad,tid)
  inf<-as.data.table(do.call('rbind',strsplit(h3.tad$info,':')))
  setnames(inf,c('tname','ensg','tbiotype','tstrand'))
  tad.genes<-cbind(h3.tad,inf)
  tad.genes$info<-NULL
  tad.genes
}

DATA.DIR<-'/Users/oliver/DATA/JAVIERRE_GWAS/'
tad.dir<-file.path(DATA.DIR,'support/TAD')
tad.files<-list.files(path=tad.dir,pattern = "*.bed",full.names = TRUE)

tad.ld.file<-file.path(DATA.DIR,'support/TAD/ld_tad.RData')
if(!file.exists(tad.ld.file)){
  ld.blocks<-file.path(DATA.DIR,'support','0.1cM_regions.b37.bed')
  ld<-fread(ld.blocks)
  setnames(ld,c('chr','start','end','id'))
  ld.gr<-with(ld,GRanges(seqnames=Rle(chr),ranges=IRanges(start=start,end=end),id=id))
  tissue<-lapply(tad.files,function(x) processTAD(x,ld.gr))
  names(tissue)<-sub('_mean_merged_75overlap.bed','',basename(tad.files))
  save(tissue,file=tad.ld.file)
}else{
  load(tad.ld.file)
}

### next load in coding SNPs
csnps.file<-file.path(DATA.DIR,'support','cSNPs_w_ENSG.e75.bed')
cs<-fread(csnps.file)
cs$V1<-NULL
setnames(cs,c('chr','start','end'))
cs.gr<-with(cs,GRanges(seqnames=Rle(chr),ranges=IRanges(start=start+1,end=end)))
cs.gr<-reduce(cs.gr)

### next load pmi data

pmi.file<-file.path(DATA.DIR,'out/pmi','RA_OKADA_IMB.pmi')
pmi<-fread(pmi.file)
pmi.gr<-with(pmi,GRanges(seqnames=Rle(chr),ranges=IRanges(start=start+1,end=end),ppi=ppi,rs=rsid))

## remove coding SNPs
ol<-as.matrix(findOverlaps(pmi.gr,cs.gr))
pmif.gr<-pmi.gr[-ol[,1],]

## remove MHC

mhc.gr<-GRanges(seqnames=Rle('6'),ranges=IRanges(start=25e6,end=35e6))
ol<-as.matrix(findOverlaps(pmif.gr,mhc.gr))
pmif.gr<-pmif.gr[-ol[,1],]

## load in annotations
h3.file<-file.path(DATA.DIR,'support','HindIII_baits_e75.bed')
h3<-fread(h3.file)
setnames(h3,c('chr','start','end','id','info'))
## there are approx 480 baits which don't appear to overlap a gene in e75 remove these
h3.p<-subset(h3,info!='')
## only use protein_coding
h3.p<-h3[grep("protein_coding",h3$info),]
h3.gr<-with(h3.p,GRanges(seqnames=Rle(chr),ranges=IRanges(start=start,end=end),id=id,info=info))
an.tissue<-lapply(tad.files,function(x) annotateTAD(x,h3.gr))
names(an.tissue)<-sub('_mean_merged_75overlap.bed','',basename(tad.files))
## next compute TAD scores for each tissue TAD.

tad.tissue<-lapply(names(tissue),function(t){
  print(paste("processing",t))
  tad.ld.gr<-tissue[[t]]
  mer<-as.data.table(mergeByOverlaps(pmif.gr,tad.ld.gr))
  mer.sm<-mer[,list(sum=sum(pmif.gr.ppi)),by=tid.ld]
  tdf<-as.data.table(do.call("rbind",strsplit(mer.sm$tid.ld,'\\.')))
  setnames(tdf,c('tid','tidld'))
  mer.sm<-cbind(mer.sm,tdf)
  cogs<-mer.sm[,list(cogs.score=prod(1-sum)),by=tid]
  cogs$cogs.score<-pmin(1-cogs$cogs.score,1)
  cogs$tid<-as.numeric(cogs$tid)
  setkey(cogs,tid)
  an<-an.tissue[[t]]
  res<-an[cogs]
  res$tissue<-t
  res[res[,.I[which.max(cogs.score)],by=ensg]$V1]
})

tad.tissue<-rbindlist(tad.tissue)

## there are some TAD's that have no baits in. We remove these as cannot assign gene
tads.no.baits<-subset(tad.tissue,is.na(ensg))
tad.tissue<-subset(tad.tissue,!is.na(ensg))

library(reshape2)

tm<-melt(as.data.frame(tad.tissue),id.vars=c('tname','ensg','tissue'),measure.vars = 'cogs.score')
tad.scores<-dcast(tm,ensg+tname~tissue+variable,fill=0)


indx <- max.col(tad.scores[,3:10], ties.method='first') + 2
tad.scores$max.score<-tad.scores[cbind(1:nrow(tad.scores), indx)]
tad.scores$max.score.tissue<-sub("_cogs.score","",names(tad.scores)[indx])

## we use COGS to create a dataset that is comparable but uses PCHiC information
cogs<-fread("/Users/oliver/gitr//cd4chic/gwas_paper/DATA//out/geneScore_hic//RA_OKADA_IMB_NO_CODING_SNPS.pmi.tab")
cogs<-subset(cogs,biotype=='protein_coding')

## coverage is quite similar 16478 for cogs and 16091 for tad.scores
tad.scores<-as.data.table(tad.scores)
ts<-tad.scores[,.(ensg,tname,max.score)]
setkey(ts,ensg)
cogs.s<-cogs[,.(ensg,name,all_gene_score)]
setkey(cogs.s,ensg)

### merge

comp<-ts[cogs.s]
setnames(comp,c('ensg','tname','tad','name','cogs'))
comp$tad<-as.numeric(comp$tad)
comp.f<-subset(comp,!is.na(tad) & !is.na(cogs))

library(ggplot2)
library(cowplot)

ggplot(comp.f,aes(x=tad,y=cogs)) + geom_point(colour='grey') + geom_abline(slope=1,intercept=0,colour='red',linetype='dashed') + 
  theme_bw(base_size = 30) + xlab('TAD Scores') + ylab('COGS Scores') + geom_vline(xintercept = 0.5 , colour = ' black') + geom_hline(yintercept=0.5 , colour = 'black') +
  annotate("text", x = 0.25, y = 0.25, label = "14759", size =10) + annotate("text", x = 0.25, y = 0.75, label = "26", size =10) +
  annotate("text", x = 0.75, y = 0.25, label = "466", size =10) + annotate("text", x = 0.75, y = 0.75, label = "103", size =10)
  

## Chris asked what number of genes in each of 16 quadrants

comp.f$tad.cat<-cut(comp.f$tad,breaks=seq(0,1,by=0.25),include.lowest = TRUE)
comp.f$cogs.cat<-cut(comp.f$cogs,breaks=seq(0,1,by=0.25),include.lowest = TRUE)
with(comp.f,table(tad.cat,cogs.cat))

## MS wants a Venn diagram
library(VennDiagram)
library(xtable)
comp.f$tad.cat<-cut(comp.f$tad,breaks=seq(0,1,by=0.5),include.lowest = TRUE)
comp.f$cogs.cat<-cut(comp.f$cogs,breaks=seq(0,1,by=0.5),include.lowest = TRUE)
tc<-with(comp.f,table(tad.cat,cogs.cat))
xtable(tc)


grid.newpage()
draw.pairwise.venn(26+103,466+103,103, category = c("COGS", "TAD"), lty = rep("blank",2), 
                   fill = c("light blue", "pink"), alpha = rep(0.5, 2), cat.pos = c(0,0), cat.dist = rep(0.025, 2),cex = rep(3, 3),cat.cex = rep(3,2))

### What is the distribution of distances from TAD boundaries for COGS prioritised genes only compared to those 
## prioritised by both ?

## For a bait what is the distance to the nearest tad boundary over all tissues ?

dist.tad<-lapply(tad.files,function(f){
  tad<-fread(f)
  tads.gr<-with(tad,GRanges(seqnames=Rle(chr),ranges=IRanges(start=start,end=end)))
  tads.gr$tid<-1:length(tads.gr)
  h3.tad<-as.data.table(mergeByOverlaps(h3.gr,tads.gr))
  mid.bait<- with(h3.tad,((h3.gr.end - h3.gr.start)/2) + h3.gr.start)
  ## compute distance to start and distance to end
  dists<-with(h3.tad,cbind(abs(tads.gr.start-mid.bait),abs(tads.gr.end-mid.bait)))
  h3.tad$dist.close.tad<-round(apply(dists,1,min))
  ret<-h3.tad[,.(id,info,tid,dist.close.tad)]
  ret<-ret[ret[,.I[which.min(dist.close.tad)],by=info]$V1]
  inf<-as.data.table(do.call('rbind',strsplit(ret$info,':')))
  setnames(inf,c('tname','ensg','tbiotype','tstrand'))
  ret<-cbind(ret,inf)
  ret$info<-NULL
  ret[,.(ensg,dist.close.tad)]
})
names(dist.tad)<-sub('_mean_merged_75overlap.bed','',basename(tad.files))

tm<-split(tad.scores$ensg,tad.scores$max.score.tissue)
ds<-lapply(names(tm),function(n){
  p<-tm[[n]]
  k<-dist.tad[[n]]
  setkey(k,ensg)
  k[p]
})
ds<-rbindlist(ds)
setkey(ds,ensg)
setkey(tad.scores,ensg)
tad.scores<-tad.scores[ds]
ts<-tad.scores[,.(ensg,max.score,dist.close.tad)]
ts<-ts[comp.f]
ts$ol.cat<-with(ts,paste0(tad.cat,cogs.cat))
ts[ts$ol.cat=='(0.5,1](0.5,1]',]$ol.cat='Both'
ts[ts$ol.cat=='[0,0.5][0,0.5]',]$ol.cat='NEITHER'
ts[ts$ol.cat=='(0.5,1][0,0.5]',]$ol.cat='TAD'
ts[ts$ol.cat=='[0,0.5](0.5,1]',]$ol.cat='COGS'
ts.f<-subset(ts,ol.cat  %in% c('Both','COGS','TAD'))
ggplot(ts.f,aes(factor(ol.cat,levels=c('COGS','TAD','Both')),dist.close.tad/1e6)) + 
  geom_boxplot() + xlab('Method') + ylab('Distance to closest TAD boundary (Mb)') +
  theme_bw(base_size = 30)

## violin plots for MS -  bin

# comp.f$tad.cat<-cut(comp.f$tad,breaks=seq(0,1,by=0.2),include.lowest = TRUE)
# ggplot(comp.f,aes(tad.cat,cogs)) + geom_violin() + theme_bw() +
#   theme(axis.text.x=element_text(angle = -90, hjust = 0)) +
#   xlab('TAD score bin') + ylab('COGs Score')



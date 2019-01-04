#!/usr/bin/env Rscript
# This file is for prioritizing T1D SNPs

## Command Arg Parameters ##
args = commandArgs(trailingOnly=T)
hmsg = 'Rscript roadmap.r [data/roadmap_dist.tsv]
  - Argument [data/roadmap_dist.tsv] is a mendatory file path including SNPs distance from RoadMap enhancers.'
if(length(args)<1|length(args)>1) stop(hmsg)
rd.path = unlist(args[1])

# System parameter
library(pbapply)
source('src/pdtime.r')
t0 = Sys.time()

#########################
## Function start here ##
#########################
# 1. RoadMap data compile
rd = unlist(read.table(rd.path,sep='\n'))
n=length(rd); gwas=NULL; enh=NULL
tmp=pblapply(rd,function(x) {
	rdc = as.character(x)
	rdc.v = unlist(strsplit(rdc,'\t'))
	if(length(rdc.v)==4) gwas = rbind(gwas,rdc.v)
	else if(length(rdc.v)==6) {
		enh_1 = rdc.v[2:6]
		pos   = paste0(enh_1[1],':',enh_1[2],'-',enh_1[3])
		enh_  = c(pos,enh_1[4:5])
		enh   = rbind(enh,enh_)
	}
})

#pb = winProgressBar(title="Loop progress",
#     label="Ready to read table..",min=0,max=n,width=500)
#for(i in 1:n) {
#	rdc = as.character(rd[i])
#	if(i%%2==1) gwas = rbind(gwas,unlist(strsplit(rdc,'\t')))
#	else {
#		enh_1 = unlist(strsplit(rdc,'\t'))[2:6]
#		pos   = paste0(enh_1[1],':',enh_1[2],'-',enh_1[3])
#		enh_ = c(pos,enh_1[4:5])
#		enh   = rbind(enh,enh_)
#	}
#    ## Progress time ##
#    setWinProgressBar(pb,i,label=paste0(round(i/n*100,0),
#              " % (",i,"/",n,") done for ",pdtime(t0,2)))
#    ###################
#}
#close(pb)

gwas = as.data.frame(gwas); print(dim(gwas))
colnames(gwas) = c('chr','start','end','rsid')
enh = as.data.frame(enh)  ; print(dim((enh)))
colnames(enh)  = c('enh_pos','num','dis')
dis = as.numeric(as.character(enh$dis))
rd.df = cbind(gwas,enh[1:2],dis)

# 2. Filtering SNPs by distance of closest enhancers
rd.df_ = subset(rd.df,dis==0)
cat(paste0('Enhancer occupied by SNPs = ',length(unique(rd.df_$rsid)),'\n'))
cat(paste0('SNPs in RoadMap enhancers = ',length(unique(rd.df_$enh_pos)),'\n'))
f.name1 = 'data/roadmap_dist_df2.tsv'
write.table(rd.df_,f.name1,row.names=F,quote=F,sep='\t')
cat(paste0('File write: ',f.name1,'\n'))

snp.enh = unique(rd.df_[,1:4])
n1 = length(unique(snp.enh$rsid))
f.name2 = paste0('data/snp_enh_',n1,'_2.bed')
write.table(snp.enh,f.name2,row.names=F,col.names=F,quote=F,sep='\t')
cat(paste0('File write: ',f.name2,'\n'))
cat(paste0('>> ',pdtime(t0,1),'\n'))
##################
## Function end ##
##################
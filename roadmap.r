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
# 1. Compiling SNP data distant from RoadMap enhancers
rd = unlist(read.table(rd.path,sep='\n'))
cat(paste0('>> Row number = ',length(rd)/2,'\n'))
rd.odd  = as.character(rd[c(TRUE,FALSE)])
rd.even = as.character(rd[c(FALSE,TRUE)])

gwas = as.data.frame(do.call(rbind,strsplit(rd.odd, '\t',fixed=T)))
enh_ = as.data.frame(do.call(rbind,strsplit(rd.even,'\t',fixed=T)))[2:6]
pos  = unlist(apply(enh_,1,function(x) paste0(x[1],':',x[2],'-',x[3])))
enh  = cbind(pos,enh_[,4:5])

colnames(gwas) = c('chr','start','end','rsid')
colnames(enh)  = c('enh_pos','num','dis')
dis = as.numeric(as.character(enh$dis))
rd.df = cbind(gwas,enh[1:2],dis)

# 2. Filtering SNPs by distance of closest enhancers
rd.df_ = unique(subset(rd.df,dis==0))
cat(paste0('\nEnhancer occupied by SNPs = ',length(unique(rd.df_$rsid)),'\n'))
cat(paste0('SNPs in RoadMap enhancers = ',length(unique(rd.df_$enh_pos)),'\n'))
f.name1 = 'data/roadmap_dist_df.tsv'
write.table(rd.df_,f.name1,row.names=F,quote=F,sep='\t')
cat(paste0('\nFile write: ',f.name1,'\n'))

snp.enh = unique(rd.df_[,1:4])
n1 = length(unique(snp.enh$rsid))
f.name2 = paste0('data/snp_enh_',n1,'.bed')
write.table(snp.enh,f.name2,row.names=F,col.names=F,quote=F,sep='\t')
cat(paste0('File write: ',f.name2,'\n'))
cat(paste0('\n>> ',pdtime(t0,1),'\n'))
##################
## Function end ##
##################
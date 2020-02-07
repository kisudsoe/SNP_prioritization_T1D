#!/usr/bin/env Rscript
# This file is for identifying random distributions from GTEx data

## Command Arg Parameters ##
# Usage: Rscript gtex_random.r
args = commandArgs(trailingOnly=T)
hmsg = 'Rscript gtex_random.r
  - No argument is needed for this function.'
if(length(args)>0) stop(hmsg)

# System parameter
suppressMessages(library(plyr))
suppressMessages(library(data.table))
source('src/pdtime.r')
t0 = Sys.time()

#########################
## Function start here ##
#########################
# 1. Loading GTEx significant file
cat("\n(1/2) Loading GTEx significant file\n")
gtex = readRDS('db/gtex_signif_5e-8.tsv.rds')
cat(paste0(' - GTEx table, rows= ',dim(gtex)[1],' cols= ',dim(gtex)[2],'\n'))

# 2. eQTL SNP filter iteration
cat("\n(2/2) eQTL SNP filter iteration\n")
dir_out = 'random/gtex'
f_name = 'random/gtex_dist.tsv'
dir.create(file.path(dir_out))

rsid = c(1:10000)
rd_name = paste0('gtex_rsid',rsid,'.bed')
rd_path = paste0('random/seeds/rsid',rsid,'.bed')
colnames(gtex)[9] = 'rsid'
rd.li = list(); n = length(rd_path)
rd.li = lapply(c(1:n),function(i) {
    rd = as.data.frame(fread(rd_path[i])) # data.table
    colnames(rd) = c('chr','start','end','rsid')
    #cat(paste0('\t(',i,'/',n,') ',rd_path[i],'\n'))
    gtex_ = subset(gtex,rsid%in%rd$rsid)
    write.table(gtex_,paste0(dir_out,'/gtex_rsid',rsid[i],'.tsv'),
                row.names=F,quote=F,sep='\t')
    
    eqtl_snp_n  = length(unique(gtex_$rsid))
    asso_gene_n = length(unique(gtex_$gene_id))
    rd.row = c(rd_name[i],eqtl_snp_n,asso_gene_n)
    if(i%%1000==0) cat(paste0('  - ',round((i/n)*100,1),'% Working process\n'))
    return(t(rd.row))
})
rd_df = ldply(rd.li,data.frame) # plyr
colnames(rd_df) = c('random','SNP_n','Gene_n')
write.table(rd_df,f_name,row.names=F,quote=F,sep='\t')
cat(paste0('\nFile write: ',f_name,'\n'))
cat(paste0('>> ',pdtime(t0,1),'\n'))
##################
## Function end ##
##################
#!/usr/bin/env Rscript
# This file is for identifying random distributions from GTEx data

## Command Arg Parameters ##
# Usage: Rscript lncrnasnp_random.r
args = commandArgs(trailingOnly=T)
hmsg = 'Rscript lncrnasnp_random.r_random.r
  - No argument is needed for this function.'
if(length(args)>0) stop(hmsg)

# System parameters
suppressMessages(library(parallel)) # for multicore process
suppressMessages(library(plyr))
suppressMessages(library(data.table))
source('src/pdtime.r')
t0 = Sys.time()
dir_out = 'random/lncrna'
dir.create(file.path(dir_out))
f_name = 'random/lncrna_dist.tsv'

#########################
## Function start here ##
#########################
# 1. Read DB files..
cat("\n(1/2) Read DB files..\n")
snplnc_path = 'db/lncRNASNP2_snplist.txt.rds'
snplnc = readRDS(snplnc_path)
cat(paste0('  - ',snplnc_path,'; ',pdtime(t0,2),'\n'))

# 2. Random SNP iteration
rsid = c(1:10000)
rd_name = paste0('lncrna_rsid',rsid,'.bed')
rd_path = paste0('random/seeds/rsid',rsid,'.bed')

cat("\n(2/2) Random SNP iteration\n")
n=length(rd_path)
rd_li=mclapply(c(1:n),function(i) { # parallel cpu process
    # Read SNPs
    snp_df = read.delim(rd_path[i],header=F)
    colnames(snp_df) = c('chr','start','end','dbsnp')
    #cat(paste0('\t(',i,'/',n,') ',rd_path[i],'\n'))
    
    # Merge and annotate SNPs within lncRNAs
    #snp_lnc = merge(snp_df,snplnc,by='dbsnp') # merge process takes long time
    snp_lnc = subset(snplnc,dbsnp%in%snp_df$dbsnp)
    f_name1 = paste0(dir_out,'/',rd_name[i])
    write.table(snp_lnc,f_name1,row.names=F,quote=F,sep='\t')
    
    # Count numbers for return
    df_n = data.frame(random  =rd_name[i],
                      SNP_n   =c(length(unique(snp_lnc$dbsnp))),
                      lncRNA_n=c(length(unique(snp_lnc$lncRNA))))
    if(i%%1000==0) cat(paste0('  - ',round((i/n)*100,1),'% ',
                              pdtime(t0,2),'\n'))
    return(df_n)
})
rd_df = ldply(rd_li,data.frame) # plyr
colnames(rd_df) = c('random','SNP_n','lncRNA_n')
write.table(rd_df,f_name,row.names=F,quote=F,sep='\t')
cat(paste0('\nFile write: ',f_name,'\n'))
cat(paste0(pdtime(t0,1),'\n'))
##################
## Function end ##
##################
#!/usr/bin/env Rscript
# This file is for filtering LDlink data

## Command Arg Parameters ##
# T1D_gwas.bat: Rscript T1D.r [SNP_file_path] [LDlink_data_folder_path]
args = commandArgs(trailingOnly=T)
snp.path = unlist(args[1])
ld.path = unlist(args[2])

# System parameters
date = Sys.Date()
t0 = Sys.time()


## Function start here ##
snpids = unlist(read.table(snp.path))
cat(paste0('Input SNP list number = ',length(snpids),'\n'))

ldlink = unlist(lapply(snpids, function(x) { paste0(ld.path,'/',x,'.tsv') }))
snptb  = data.frame(snpids=snpids, ldlink=ldlink)

library(plyr)
source('src/pdtime.r')
ldlink.li = NULL # List variable for ldlink data
n = nrow(snptb); t1 = Sys.time()
pb = winProgressBar(title="Loop progress",
     label="Ready to read files..",min=0,max=n,width=500)
for(i in 1:n) {
    t = try(read.table(as.character(snptb[i,2]),header=T))
    if("try-error" %in% class(t)) ldlink[[i]] = NULL # get error from empty file
    else { # No errors
        tb1 = read.table(as.character(snptb[i,2]),header=T)
        tb2 = data.frame(SNPid=rep(snptb[i,1],nrow(tb1)),tb1)
        ldlink.li[[i]] = tb2
    }
    ## Progress time ##
    setWinProgressBar(pb,i,label=paste0(round(i/n*100,0),
                      " % (",i,"/",n,") done for ",
                      pdtime(t1,2)))
    ###################
}
close(pb); #print(length(ldlink.li))
ldlink.df = ldply(ldlink.li, data.frame)
#print(dim(ldlink.df))
ldlink_ = unique(subset(ldlink.df,R2>0.6 & Dprime==1)$RS_Number

snp.t1 = unique(snpids); cat(paste0('SNP Tier1 = ',length(snp.t1),'\n'))
snp.t2 = setdiff(ldlink_, snp.t1)); cat(paste0('SNP Tier2 = ',length(snp.t2),'\n'))
snp.seed = union(snp.t1,snp.t2); cat(paste0('SNP seed  = ',length(snp.seed),'\n'))
pdtime(t0,1)

## Function end ##
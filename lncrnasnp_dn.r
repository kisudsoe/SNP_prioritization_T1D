#!/usr/bin/env Rscript
# This file is for downloading lncRNASNP2 data

## Command Arg Parameters ##
# roadmap.bat: Rscript lncrnasnp_dn.r
args = commandArgs(trailingOnly=T)
hmsg = 'Rscript lncrnasnp_dn.r has no need any argument.'
if(length(args)>0) stop(hmsg)

# System parameter
source('src/saveasrds.r')
source('src/pdtime.r')
t0 = Sys.time()

#########################
## Function start here ##
#########################
input = list(
    #snplist=c(url='http://bioinfo.life.hust.edu.cn/static/lncRNASNP2/downloads/snps_mod.txt',
    #        path='db/lncRNASNP2_snplist.txt'),
    lncrna =c(url='http://bioinfo.life.hust.edu.cn/static/lncRNASNP2/downloads/lncrnas.txt',
            path='db/lncrnas.txt'),
    disease= c(url='http://bioinfo.life.hust.edu.cn/static/lncRNASNP2/downloads/lncRNA_associated_disease_experiment.txt',
            path='db/lncrna-diseases_experiment.txt'))

# 1. Download lncRNASNP2 data
k1=lapply(input,function(l) {
    tb     = try(download.file(l[1],destfile=l[2]))
    if("try-error"%in% class(tb)) stop('The file url seems to be changed. Please check the lncRNASNP2 homepasge.')
})
cat(paste0('>> ',pdtime(t0,2),'\n'))

# 2. Generate rds files
k2=lapply(input,function(l) {
    saveasrds(l[2])
})
cat(paste0('>> ',pdtime(t0,1),'\n'))
##################
## Function end ##
##################
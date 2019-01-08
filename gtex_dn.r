#!/usr/bin/env Rscript
# This file is for filtering GTEx data

## Command Arg Parameters ##
# Usage: Rscript gtex_dn.r
args = commandArgs(trailingOnly=T)
hmsg = 'Rscript gtex.r
  - No argument is needed.'
if(length(args) > 0) stop(hmsg)

# System parameter
library(tools)
source('src/pdtime.r')
t0 = Sys.time()

#########################
## Function start here ##
#########################
url = 'https://storage.googleapis.com/gtex_analysis_v7/single_tissue_eqtl_data/GTEx_Analysis_v7_eQTL.tar.gz'
dir = 'db/'
name= paste0(dir,basename(url))
tb = try(download.file(url,destfile=name))
if('try-error' %in% class(tb)) stop('File address seems to be changed.')
f.li = untar(name,exdir=dir)

#f.li = list.files(paste0(dir,file_path_sans_ext(basename(url))))
f.li = untar(name,list=TRUE)
f.name = paste0(dir,'gtex_files.txt')
write.table(f.li,f.name,row.names=F,col.names=F,quote=F,sep='\t')
cat(paste0('File write: ',f.name,'\n'))
cat(paste0('>> ',pdtime(t0,1),'\n'))
##################
## Function end ##
##################
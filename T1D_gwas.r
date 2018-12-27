#!/usr/bin/env Rscript
# This file is for filtering GWAS Catalog data

## Command Arg Parameters ##
# T1D_gwas.bat: Rscript T1D.r [GWAS_file_path] [p-value_criteria]
args = commandArgs(trailingOnly=T)
path = unlist(args[1])
pval = unlist(args[2])

# System parameters
t0 = Sys.time()


## Function start here ##
source('src/pdtime.r')
cat(paste0('Process initiate at ',t0,'\n'))
if(is.na(pval)) pval=5*10^-8
if(!is.na(path)) {
    cat(paste0(' >> Loading src/gwas_filt.r function <<\n'))
    source('src/gwas_filt.r')
    gwas_filt(path,pval)
    #cat(paste0(pdtime(t0,2),'\n'))
} else { print("path parameter doesn't exist.") }

## Function end ##
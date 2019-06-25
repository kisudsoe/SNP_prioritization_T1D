#!/usr/bin/env Rscript
# This file is for filtering GWAS Catalog data

## Command Arg Parameters ##
# T1D_gwas.bat: Rscript T1D_gwas.r [GWAS_file_path] [p-value_criteria]
args = commandArgs(trailingOnly=T)
msg = 'Rscript T1D_gwas.r [GWAS_file_path] [p-value_criteria]
  - Argument [GWAS_file_path] is mendatory.
  - Argument [p-value_criteria] is optional.
  - Default p-value is set as 5e-08.\n'
if(length(args)<1|length(args)>2) stop(msg)
path = unlist(args[1])
pval = as.numeric(unlist(args[2]))

# System parameters
source('src/pdtime.r')
t0 = Sys.time()

#########################
## Function start here ##
#########################
cat(paste0('Process initiate at ',t0,'\n'))
if(is.na(pval)) pval=5e-08
if(!is.na(path)) {
    cat(paste0(' >> Loading src/gwas_filt.r function <<\n\n'))
    cat(paste0('  - Input file: ',path,'\n'))
    
    snp = read.delim(path)
    cat(paste0('  - Table;    rows= ',dim(snp)[1],', cols= ',dim(snp)[2],'\n'))
    cat(paste0('  - Criteria      = ',pval,'\n'))

    snp_   = subset(snp,`P.VALUE`<pval)
    cat(paste0('  - Subset;   rows= ',dim(snp_)[1],', cols= ',dim(snp_)[2],'\n'))
    coord  = paste0('chr',snp_$`CHR_ID`,':',snp_$`CHR_POS`); head(coord)
    snp_df = data.frame(rsid=snp_$`SNPS`,coord=coord,pval=snp_$`P.VALUE`,cytoband=snp_$`REGION`)
    snp_df = unique(snp_df)
    snp_n  = length(unique(snp_df$rsid))
    cat(paste0('  - IDs (p < 5e-8)= ',snp_n,'\n\n'))
    
    f_name = paste0('data/gwas_',pval,'_',snp_n,'.tsv')
    write.table(snp_df,f_name,sep='\t',quote=F,row.names=F) #,col.names=F
    cat(paste0(' >> Write tsv file: ',f_name,' <<\n'))
    cat(paste0(pdtime(t0,1),'\n'))
} else { print("path parameter doesn't exist.") }
##################
## Function end ##
##################
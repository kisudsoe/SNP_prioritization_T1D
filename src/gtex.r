#!/usr/bin/env Rscript
# This file is for filtering RoadMap data

## Command Arg Parameters ##
# Usage: Rscript gtex.r data/seedSNP_1817.bed db/gtex_signif_3e-4.tsv.gz
# Usage: Rscript gtex.r data/seedSNP_1817.bed db/gtex_signif_5e-8.tsv.gz
args = commandArgs(trailingOnly=T)
hmsg = 'Rscript gtex_filt.r [SNP_file_path] [GTEx_download_target_dir]
  - Argument [SNP_file_path] is mandatory.
  - Argument [GTEx_download_target_dir] is mandatory.'
if(length(args)<2|length(args)>2) stop(hmsg)
snp.path = args[1]
gte.path = args[2]

# System parameter
library(tools)
library(plyr)
source('src/pdtime.r')
t0 = Sys.time()

#########################
## Function start here ##
#########################
snp = read.delim(snp.path,header=F)
colnames(snp) = c('chr','start','end','ann')
snpids = gsub('(.*)_.*','\\1',snp[,4])
len.snpids = prettyNum(length(snpids),big.mark=',',preserve.width='none')
cat(paste0('Input SNPs number = ',len.snpids,'\n'))

# 1. Loading GTEx significant file
cat("\n1. Loading GTEx significant file\n")
if(file_ext(gte.path)=='gz') { gte = read.delim(gzfile(gte.path),header=T)
} else if(file_ext(gte.path)=='rds') gte = readRDS(gte.path)
colnames(gte)[9] = 'rsid'
dim.gte = prettyNum(dim(gte),big.mark=',',preserve.width='none')
cat(paste0('  ',basename(gte.path),': rows= ',dim.gte[1],' cols= ',dim.gte[2],'\n'))
cat(paste0('  ',pdtime(t0,2),'\n'))

# 2. eQTL SNP filteration
cat("\n2. eQTL SNP filteration\n")
gte_ = subset(gte,rsid%in%snpids)
dim.gte_ = prettyNum(dim(gte_),big.mark=',',preserve.width='none')
cat(paste0('  Overlapped table, rows= ',dim.gte_[1],' cols= ',dim.gte_[2],'\n'))
len.gte.rsid = prettyNum(length(unique(gte_$rsid)),big.mark=',',preserve.width='none')
cat(paste0('  eQTL SNPs = ',len.gte.rsid,'\n'))
len.gte.gene = prettyNum(length(unique(gte_$gene_id)),big.mark=',',preserve.width='none')
cat(paste0('  Associated genes = ',len.gte.gene,'\n'))
f.name1 = paste0('data/gtex_5e-08_',length(unique(gte_$rsid)),'.tsv')
write.table(gte_,f.name1,row.names=F,col.names=T,quote=T,sep='\t')
cat(paste0('File write: ',f.name1,'\n'))

# 3. eQTL SNP BED file generation
cat("\n3. eQTL SNP BED file generation\n")
snp.df = data.frame(snp,rsid=snpids)
snp_ = subset(snp.df,rsid%in%gte$rsid)[1:4]
dim.snp_ = prettyNum(dim(snp_),big.mark=',',preserve.width='none')
cat(paste0('  GTEx SNP BED, rows= ',dim.snp_[1],' cols= ',dim.snp_[2],'\n'))
len.snp = prettyNum(length(unique(snp_$ann)),big.mark=',',preserve.width='none')
cat(paste0('  eQTL SNPs = ',len.snp,'\n'))
f.name2 = paste0('data/snp_',length(unique(snp_$ann)),'_gtex.bed')
write.table(snp_,f.name2,row.names=F,col.names=F,quote=T,sep='\t')
cat(paste0('File write: ',f.name2,'\n'))
cat(paste0('\n',pdtime(t0,1),'\n'))
##################
## Function end ##
##################
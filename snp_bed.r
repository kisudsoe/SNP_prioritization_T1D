#!/usr/bin/env Rscript
# This file is for annotating SNPs using Ensembl biomaRt

## Command Arg Parameters ##
# T1D.bat: Rscript snp_bed.r [seedSNP_annotation_file_path]
args = commandArgs(trailingOnly=T)
if(length(args)<1) stop("[seedSNP_annotation_file_path] argument must be supplied.")
if(length(args)>2) stop("Too many arguments was detected. Argument [seedSNP_annotation_file_path] is mendatory.")
path = unlist(args[1])
if(length(args)>1) ann.path  = unlist(args[2])

# system parameters
library(biomaRt)
source('src/pdtime.r') # pdtime(time,1/2); 1= Job done, 2= Job process
source('src/mktsv.r')
t0 = Sys.time(); cat(paste0('>> Process initiate at ',t0,'\n\n'))


#########################
## Function start here ##
#########################
# Ensembl biomart (grch37)
#snps = read.delim(path); colnames(snps)[1]='rsid'
#cat(paste0('Input contents, rows= ',dim(snps)[1],' cols= ',dim(snps)[2],'\n'))
#hg19_snp = useMart(biomart="ENSEMBL_MART_SNP", host="grch37.ensembl.org",
#                   dataset="hsapiens_snp", path="/biomart/martservice")
#snp_attr = c("refsnp_id","chr_name","chrom_start") #,"chrom_end"
#snps_ = getBM(attributes = snp_attr,
#              filters    = "snp_filter",
#              values     = snps,
#              mart       = hg19_snp)
#snps_[,2] = paste0('chr',snps_[,2])
#snps_[,3] = snps_[,3]
#colnames(snps_) = c('rsid','chr','pos')
#snp_1 = merge(snps,snps_,by='rsid',all.x=TRUE)

# Optional annotation file such as `seedSNP_1817_ldlink`
ann = read.delim(ann.path)
snp_2 = unique(ann[2:ncol(ann)])
colnames(snp_2)[1] = "rsid"
snp_2a= strsplit(as.character(snp_2$coord),":")
snp_2b= as.data.frame(do.call(rbind,snp_2a))
colnames(snp_2b) = c("chr","pos")
snp_2c= data.frame(rsid =snp_2$rsid,
                   chr  =snp_2b$chr,
                   pos  =snp_2b$pos)
cat(paste0('Annotations,    rows= ',dim(snp_2)[1],' cols= ',dim(snp_2)[2],'\n'))

snp_bed = unique(rbind(snp_1,snp_2c))
cat(paste0('Merge table,    rows= ',dim(snp_bed)[1],' cols= ',dim(snp_bed)[2],'\n'))
f.name = paste0(path,'.bed')
mktsv(file=snp_bed, path=f.name)

cat(paste0('\n>> ',pdtime(t0,1),'\n'))
##################
## Function end ##
##################
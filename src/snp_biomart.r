#!/usr/bin/env Rscript
# This file is for annotating SNPs using Ensembl biomaRt

## Command Arg Parameters ##
# T1D.bat: Rscript snp_biomart.r
args = commandArgs(trailingOnly=T)
if(length(args)>0) stop("Rscript snp_biomart.r, no argument is needed.")

# system parameters
library(biomaRt)
source('src/pdtime.r') # pdtime(time,1/2); 1= Job done, 2= Job process
t0 = Sys.time(); cat(paste0('>> Process initiate at ',t0,'\n\n'))


#########################
## Function start here ##
#########################
# Ensembl biomart (grch37)
snps = read.delim(path); colnames(snps)[1]='rsid'
cat(paste0('Input contents, rows= ',dim(snps)[1],' cols= ',dim(snps)[2],'\n'))
hg19_snp = useMart(biomart="ENSEMBL_MART_SNP", host="grch37.ensembl.org",
                   dataset="hsapiens_snp", path="/biomart/martservice")
snp_attr = c("refsnp_id","chr_name","chrom_start") #,"chrom_end"
snps_ = getBM(attributes = snp_attr,
              filters    = "snp_filter",
              values     = snps,
              mart       = hg19_snp)
snps_[,2] = paste0('chr',snps_[,2])
snps_[,3] = snps_[,3]
colnames(snps_) = c('rsid','chr','pos')
snp_bed = merge(snps,snps_,by='rsid',all.x=TRUE)
snp_bed_= data.frame(
	chr  =snp_bed$chr
	start=as.numeric(as.character(snp_bed$pos))-1
	end  =snp_bed$pos
	rsid =snp_bed$rsid
)

cat(paste0('Table, rows= ',dim(snp_bed_)[1],' cols= ',dim(snp_bed_)[2],'\n'))
f.name = paste0(path,'.bed')
write.table(snp_bed_,f.name,row.names=F,col.names=F,quote=F,sep='\t')
cat(paste0('\n>> ',pdtime(t0,1),'\n'))
##################
## Function end ##
##################
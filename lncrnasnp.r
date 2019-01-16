#!/usr/bin/env Rscript
# This file is for filtering RoadMap data

## Command Arg Parameters ##
# CMD command: Rscript lncrnasnp.r data/seedSNP_1817.bed db/lncRNASNP2_snplist.txt.gz db/lncrnas.txt.gz db/lncrna-diseases_experiment.txt.gz
args = commandArgs(trailingOnly=T)
hmsg = 'Rscript lncrnasnp.r [SNP_BED_file_path] [lncRNAsnp2_SNP_list_file_path] [lncRNAsnp2_lncRNA_list_file_path] [lncRNAsnp2_diseases_list_file_path]
  - [SNP_BED_file_path] is a mendatory argument for your SNP list with BED format.
  - [lncRNAsnp2_SNP_list_file_path] is a mendatory argument for SNP-lncRNA pair list.
  - [lncRNAsnp2_lncRNA_list_file_path] is a mendatory argument for lncRNA annotation.
  - [lncRNAsnp2_diseases_list_file_path] is an optional argument for experimental validate lncRNA-associated diseases.'
if(length(args)<3|length(args)>4) stop(hmsg)
snp_path    = args[1]
snplnc_path = args[2]
ann_path    = args[3]
if(length(args)>3) dis_path = args[4]
if(length(args)>3) print(paste0('>> ',dis_path))

# System parameter
source('src/data_table.r') # data_table(df)
source('src/pdtime.r')
t0 = Sys.time()

#########################
## Function start here ##
#########################
# 1. Read files..
cat("\n(1/3) Read files..\n")
snp = read.delim(snp_path,header=F)
colnames(snp) = c('chr','start','end','snp_id')
rsid = gsub('(.*)_.*','\\1',snp[,4])
snp.df = cbind(snp,rsid)
cat(paste0('  - ',snp_path,';\t',pdtime(t0,2),'\n'))
snplnc = read.delim(gzfile(snplnc_path),header=T)
cat(paste0('  - ',snplnc_path,';\t',pdtime(t0,2),'\n'))
ann = read.delim(gzfile(ann_path),header=T)
cat(paste0('  - ',ann_path,';\t',pdtime(t0,2),'\n'))
df = data.frame(path=c(snp_path,   snplnc_path,   ann_path),
				nrow=c(dim(snp)[1],dim(snplnc)[1],dim(ann)[1]),
				ncol=c(dim(snp)[2],dim(snplnc)[2],dim(ann)[2]))

if(length(args)>3) {
	dis = read.delim(gzfile(dis_path),header=T)
	cat(paste0('  - ',dis_path,';\t',pdtime(t0,2),'\n'))
	df = rbind(df,c(dis_path,dim(dis)[1],dim(dis)[2]))
}
knitr::kable(df)

# 2. Overlapping lncRNA to my SNP list and binding annotation..
cat("\n(2/3) Overlapping lncRNA to my SNP list and binding annotation..\n")



# 3. 
cat("\n(3/3) ..\n")

cat(paste0(pdtime(t0,1),'\n'))
##################
## Function end ##
##################
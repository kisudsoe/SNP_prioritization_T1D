#!/usr/bin/env Rscript
# This file is for filtering RoadMap data

## Command Arg Parameters ##
# Usage: Rscript gtex_filt.r 3e-04
args = commandArgs(trailingOnly=T)
hmsg = 'Rscript gtex_filt.r [p-value_criteria]
  - An argument [p-value_criteria] is needed no argument.
  - The default p-value criteria is p < 3e-04.'
if(length(args)<1|length(args)>1) stop(hmsg)
dir  = 'db/GTEx_Analysis_v7_eQTL'
pval = args[1]

# System parameter
library(tools)
library(plyr)
source('src/pdtime.r')
t0 = Sys.time()

#########################
## Function start here ##
#########################
# 1. Load BED files
cat(">> Loading GTEx BED files\n")
path   = unlist(read.delim('db/gtex_files.txt',header=F))
name   = file_path_sans_ext(file_path_sans_ext(path))
#surfix = file_ext(name)
surfix = sub('.*/.*.v7.','',name)
f.df   = data.frame(path=paste0('db/',path),surfix)
f.sig  = subset(f.df,surfix=='signif_variant_gene_pairs')$path
f.sig  = as.character(unlist(f.sig))

gte.li=list(); n=length(f.sig)
pb = winProgressBar(title="Loop progress",
     label="Ready to read files..",min=0,max=n,width=500)
cat('File reading...\n')
for(i in 1:3) {
	cat(paste0('(',i,'/',n,') ',f.sig[i],'\n'))
	if(file.exists(f.sig[i])) {
		gte.li[[i]] = read.delim(gzfile(f.sig[i]),header=T)
	} else {
		cat(paste0('No such file: ',f.sig[i]))
	}
	## Progress time ##
    setWinProgressBar(pb,i,label=paste0(round(i/n*100,0),
              " % (",i,"/",n,") done for ",pdtime(t0,2)))
    ###################
}
close(pb)
gte.df = ldply(gte.li,data.frame) # plyr
cat(paste0(' - GTEx table, rows= ',dim(gte.df)[1],' cols= ',dim(gte.df)[2],'\n'))
cat(paste0(' - BED file read complete. ',pdtime(t0,2),'\n'))

# 2. Loading annotation file
cat("\n>> Loading GTEx SNP annotation file\n")
path2 = 'db/GTEx_Analysis_2016-01-15_v7_WholeGenomeSeq_635Ind_PASS_AB02_GQ20_HETX_MISS15_PLINKQC.lookup_table.txt.gz'
ann = read.delim(gzfile(path2),header=T)
cat(paste0(' - GTEx annotation, rows= ',dim(ann)[1],' cols= ',dim(ann)[2],'\n'))
gte.ann = merge(gte.df[,c(1:2,7,10:12)],ann[,c(1:3,7)],by='variant_id')

# 3. Filtering GTEx data by Normalized p-value
cat("\n>> Filtering GTEx data by SNP-gene expression association p-value\n")
gte.sig = subset(gte.ann,pval_nomial<pval)
f.name = paste0('db/gtex_signif.tsv')
write.table(gte.sig,f.name,row.names=F,col.names=T,quote=F,sep='\t')
cat(paste0(' - GTEx significant, rows= ',dim(gte.sig)[1],' cols= ',dim(gte.sig)[2],'\n'))
cat(paste0('\n>> ',pdtime(t0,1),'\n'))
##################
## Function end ##
##################
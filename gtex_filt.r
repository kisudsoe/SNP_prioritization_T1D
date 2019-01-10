#!/usr/bin/env Rscript
# This file is for filtering RoadMap data

## Command Arg Parameters ##
# CMD Usage: Rscript gtex_filt.r 0.0003
# CMD Usage: Rscript gtex_filt.r 0.00000005
args = commandArgs(trailingOnly=T)
hmsg = 'Rscript gtex_filt.r [p-value_criteria]
  - An argument [p-value_criteria] is needed no argument.
  - The default p-value criteria is p < 0.0003.'
if(length(args)<1|length(args)>1) stop(hmsg)
dir  = 'db/GTEx_Analysis_v7_eQTL'
pval = args[1]
pval = as.numeric(as.character(pval))
cat(paste0('p-value threshold = ',pval,'\n'))

# System parameter
library(tools)
library(plyr)
source('src/pdtime.r')
t0 = Sys.time()

#########################
## Function start here ##
#########################
# 1. Load BED files
cat("\n(1/3) Loading GTEx BED files\n")
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
cat('  - File reading...\n')
for(i in 1:n) {
	tissue = sub('.*/(.*).v7..*','\\1',f.sig[i])
	cat(paste0('  (',i,'/',n,') ',tissue,'\n'))
	if(file.exists(f.sig[i])) {
		#gte.li[[i]] = read.delim(gzfile(f.sig[i]),header=T)
		tb1 = read.delim(gzfile(f.sig[i]),header=T)
		tb2 = data.frame(tb1,tissue=rep(tissue,nrow(tb1)))
		gte.li[[i]] = tb2
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
print(summary(data.frame(gte.df$pval_nominal)))
cat(paste0(' - GTEx table, rows= ',dim(gte.df)[1],' cols= ',dim(gte.df)[2],'\n'))
cat(paste0(' - BED file read complete. ',pdtime(t0,2),'\n'))

# 2. Filtering GTEx data by Normalized p-value
cat("\n(2/3) Filtering by nominal p-value\n")
gte.sig = subset(gte.df,pval_nominal<pval)
print(summary(data.frame(gte.sig$pval_nominal)))
cat(paste0(' - GTEx significant, rows= ',dim(gte.sig)[1],' cols= ',dim(gte.sig)[2],'\n'))
cat(paste0(' - ',pdtime(t0,2),'\n'))

# 3. Loading and mering annotation file
cat("\n(3/3) Loading annotation file\n")
path2 = 'db/GTEx_Analysis_2016-01-15_v7_WholeGenomeSeq_635Ind_PASS_AB02_GQ20_HETX_MISS15_PLINKQC.lookup_table.txt.gz'
ann = read.delim(gzfile(path2),header=T)
cat(paste0(' - Annotation file read complete. ',pdtime(t0,2),'\n'))
cat(paste0(' - Annotation file, rows= ',dim(ann)[1],' cols= ',dim(ann)[2],'\n'))
gte.ann = merge(gte.sig[,c(1:2,7:9,13)],ann[,c(1:3,7)],by='variant_id')
cat(paste0(' - GTEx annotation, rows= ',dim(gte.ann)[1],' cols= ',dim(gte.ann)[2],'\n'))

f.name = paste0('db/gtex_signif.tsv')
write.table(gte.ann,f.name,row.names=F,col.names=T,quote=F,sep='\t')
cat(paste0('\n>> File write: ',f.name,'\n'))
cat(paste0('>> ',pdtime(t0,1),'\n'))
##################
## Function end ##
##################
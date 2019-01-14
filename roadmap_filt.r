#!/usr/bin/env Rscript
# This file is for filtering RoadMap data

## Command Arg Parameters ##
# roadmap.bat: Rscript roadmap_filt.r
args = commandArgs(trailingOnly=T)
hmsg = 'Rscript roadmap_filt.r is needed no argument.'
if(length(args)>0) stop(hmsg)
dir = 'db/roadmap/'

# System parameter
library(rtracklayer)
library(plyr)
source('src/pdtime.r')
t0 = Sys.time()

#########################
## Function start here ##
#########################
# 1. Load BED files
cid = as.character(formatC(c(1:129),width=3,flag='0'))
path = paste0(dir,'E',cid,'_25_imputed12marks_hg38lift_dense.bed.gz')
road.li=list(); n=length(path)
pb = winProgressBar(title='Loading files',
	 label='Ready to read BED files..',min=0,max=n,width=500)
for(i in 1:n) {
	tb1 = try(import(path[i],format='bed')) # rtracklayer
	if('try-error' %in% class(tb1)) road.li[[i]] = NULL # get error from no file
	else road.li[[i]] = tb1 # No errors
	## Progress time ##
    setWinProgressBar(pb,i,label=paste0(round(i/n*100,0),
              " % (",i,"/",n,") done for ",pdtime(t0,2)))
    ###################
}
close(pb)
cat(paste0('\n>> BED file read, ',pdtime(t0,2),'\n'))
road.df = ldply(road.li, data.frame)

# 2. Filtering by enhancers and save as BED file
cat(">> Filtering RoadMap data by enhancers\n")
road.enh= subset(road.df,name %in% c("13_EnhA1","14_EnhA2","15_EnhAF","16_EnhW1","17_EnhW2","18_EnhAc"))
f.name = 'db/roadmap_enh.bed'
write.table(road.enh[,c(1:3,6)],f.name,row.names=F,col.names=F,quote=F,sep='\t')
cat(paste0('File write: ',f.name,'\n'))
cat(paste0('>> ',pdtime(t0,1),'\n'))
##################
## Function end ##
##################
#!/usr/bin/env Rscript
# This file is for downloading RoadMap data

## Command Arg Parameters ##
# roadmap.bat: Rscript roadmap_dn.r
args = commandArgs(trailingOnly=T)
hmsg = 'Rscript roadmap_dn.r is needed no argument.'
if(length(args)>0) stop(hmsg)

# mkdir db/roadmap
dir = 'db/roadmap/'
if(file.exists(dir)) { cat(paste0('Directory exists: ',dir,'\n'))
} else {
	dir.create(file.path(dir))
	cat(paste0('Directory generated: ',dir,'\n'))
}

# System parameter
source('src/pdtime.r')
t0 = Sys.time()

#########################
## Function start here ##
#########################
cid = as.character(formatC(c(1:129),width=3,flag='0'))
urls = paste0('https://egg2.wustl.edu/roadmap/data/byFileType/chromhmmSegmentations/ChmmModels/imputed12marks/jointModel/final/E',cid,'_25_imputed12marks_hg38lift_dense.bed.gz')
#head(urls)
t=lapply(urls,function(url) {
	f.name = paste0(dir,basename(url))
	tb = try(download.file(url,destfile=f.name)) # debug no file error
	if("try-error" %in% class(tb)) tb=NULL
})
cat(paste0('>> ',pdtime(t0,1),'\n'))
##################
## Function end ##
##################
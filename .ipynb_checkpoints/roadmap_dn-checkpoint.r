#!/usr/bin/env Rscript
# This file is for downloading RoadMap data.

## Command Arg Parameters ##
# roadmap.bat: Rscript roadmap_dn.r
args = commandArgs(trailingOnly=T)
hmsg = 'Rscript roadmap_dn.r is needed no argument.'
if(length(args)>0) stop(hmsg)

# System parameter
source('src/pdtime.r')
source('src/saveasrds.r')
t0 = Sys.time()

# mkdir db/roadmap
dir = 'db/roadmap/'
if(file.exists(dir)) { cat(paste0('Directory exists: ',dir,'\n'))
} else {
	dir.create(file.path(dir))
	cat(paste0('Directory generated: ',dir,'\n'))
}

#########################
## Function start here ##
#########################
cid = as.character(formatC(c(1:129),width=3,flag='0'))
urls = paste0('https://egg2.wustl.edu/roadmap/data/byFileType/chromhmmSegmentations/ChmmModels/imputed12marks/jointModel/final/E',cid,'_25_imputed12marks_hg38lift_dense.bed.gz')
#head(urls)
t=lapply(urls,function(url) {
	f_name = paste0(dir,basename(url))
	tb = try(download.file(url,destfile=f_name)) # debug no file error
	if("try-error" %in% class(tb)) tb=NULL
    cat(paste0('\n>> ',f_name,'\n'))
    try(R.utils::gunzip(f_name))
    f_name2 = tools::file_path_sans_ext(f_name)
    cat(paste0('>> ',f_name2,'\n'))
    try(bedasrds(f_name2))
})
cat(paste0('\n>> ',pdtime(t0,1),'\n'))
##################
## Function end ##
##################
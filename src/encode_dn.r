#!/usr/bin/env Rscript
# This file is for downloading RoadMap data.

## Command Arg Parameters ##
# roadmap.bat: Rscript encode_dn.r
args = commandArgs(trailingOnly=T)
hmsg = 'Rscript encode_dn.r is needed no argument.'
if(length(args)>0) stop(hmsg)

# System parameter
source('src/pdtime.r')
source('src/saveasrds.r')
t0 = Sys.time()

# mkdir db/encode
dir = 'db/encode/'
if(file.exists(dir)) { cat(paste0('Directory exists: ',dir,'\n'))
} else {
	dir.create(file.path(dir))
	cat(paste0('Directory generated: ',dir,'\n'))
}

#########################
## Function start here ##
#########################
url = paste0('http://hgdownload.cse.ucsc.edu/goldenpath/hg19/encodeDCC/wgEncodeRegTfbsClustered/wgEncodeRegTfbsClusteredV3.bed.gz')
f_name = paste0(dir,basename(url))
tb = try(download.file(url,destfile=f_name)) # debug no file error
if("try-error" %in% class(tb)) tb=NULL
cat(paste0('\n>> ',f_name,'\n'))
cat(paste0('\n>> ',pdtime(t0,1),'\n'))
##################
## Function end ##
##################
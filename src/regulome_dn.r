#!/usr/bin/env Rscript
# This file is for downloading RegulomeDB data

## Command Arg Parameters ##
# roadmap.bat: Rscript regulome_dn.r
args = commandArgs(trailingOnly=T)
hmsg = 'Rscript regulome_dn.r has no need any argument.'
if(length(args)>0) stop(hmsg)

# System parameter
source('src/saveasrds.r')
source('src/pdtime.r')
t0  = Sys.time()

# mkdir
dir = 'db_gwas/regulome/'
url_base = 'http://www.regulomedb.org/downloads/RegulomeDB.'
if(file.exists(dir)) {
    cat(paste0('Directory exists: ',dir,'\n'))
} else {
    dir.create(file.path(dir))
    cat(paste0('Directory created: ',dir,'\n'))
}

#########################
## Function start here ##
#########################
input = list(
    #total= c(url=paste0(url_base,'dbSNP141.txt.gz'),
    #    path=paste0(dir,'dbSNP141.txt.gz')),
    cat1 = c(url=paste0(url_base,'dbSNP132.Category1.txt.gz'),
        path=paste0(dir,'dbSNP132.Category1.txt.gz')),
    cat2 = c(url=paste0(url_base,'dbSNP132.Category2.txt.gz'),
        path=paste0(dir,'dbSNP132.Category2.txt.gz'))
)

cat(paste0('Download RegulomeDB data'))
for(i in 1:length(input)) {
    tb = try(download.file(input[[i]][1],destfile=input[[i]][2]))
    if('try-error' %in% class(tb)) stop('The file url seems to be changed. Please check the lncRNASNP2 homepasge: http://www.regulomedb.org/downloads')

    # Save as RDS file
    saveasrds(input[[i]][2])
}
cat(paste0('>> ',pdtime(t0,1),'\n'))
##################
## Function end ##
##################
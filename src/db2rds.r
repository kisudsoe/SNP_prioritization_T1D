#!/usr/bin/env Rscript
# This file is for converting DB data to rds files.

## Command Arg Parameters ##
# roadmap.bat: Rscript src/db2rds.r
args = commandArgs(trailingOnly=T)
hmsg = 'Rscript src.db2rds.r is no need any argument.'
if(length(args)>0) stop(hmsg)

# System parameter
source('src/saveasrds.r')
source('src/pdtime.r')
t0 = Sys.time()

#########################
## Function start here ##
#########################
db = list(
    regulomedb=c('db/RegulomeDB.dbSNP132.Category1.txt',
                 'db/RegulomeDB.dbSNP132.Category2.txt'),
    encode    =c('db/wgEncodeRegTfbsClusteredV3.bed')
    gtex      =c('db/gtex_signif_5e-8.tsv'),
    lncRNA    =c('db/lncRNASNP2_snplist.txt',
                 'db/lncrnas.txt',
                 'db/lncrna-diseases_experiment.txt')
)

k1=lapply(db,function(paths) {
    k2=lapply(paths,function(f_path) {
        cat(paste0('\nGenerate RDS: ',f_path,'\n'))
        saveasrds(f_path)
    })
})
cat(paste0('>> ',pdtime(t0,1),'\n'))
##################
## Function end ##
##################
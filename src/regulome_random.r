#!/usr/bin/env Rscript
# This file is for converting Regulome data as RDS format file.

## Command Arg Parameters ##
# CMD command: Rscript regulome_random.r db/RegulomeDB.dbSNP132.Category1.txt.rds db/RegulomeDB.dbSNP132.Category2.txt.rds
args = commandArgs(trailingOnly=T)
hmsg = 'Rscript regulome.r [regulome_file_path_1] [regulome_file_path_2] [regulome_file_path_3]
  - At least one and up to three [regulome_file_path] arguments are needed to be supplied.'
if(length(args)<1) stop(hmsg)
f_path = NULL
for(i in 1:length(args)) {
    f_path = c(f_path,unlist(args[i]))
}

# System parameter
library(tools)
library(plyr)
source('src/pdtime.r')
source('src/saveasrds.r')
t0 = Sys.time()

#########################
## Function start here ##
#########################
# 1. Read RegulomeDB files
cat("\n(1/2) Read RegulomeDB files\n")
n = length(f_path); cat(paste0('Input regulome files = ',n,'\n'))
reg.li = list()
for(i in 1:n) {
	cat(paste0('\nRead input file ',i,': ',f_path[i],'\n'))
	if(file_ext(f_path[i])=='gz') reg.li[[i]] = try(read.delim(gzfile(f_path[i]),header=F))
    else if(file_ext(f_path[i])=='rds') reg.li[[i]] = try(readRDS(f_path[i]))
	else reg.li[[i]] = try(read.delim(f_path[i],header=F))
	colnames(reg.li[[i]]) = c('chr','id','rsid','description','level')
	cat(paste0('Table row = ',dim(reg.li[[i]])[1],', col = ',dim(reg.li[[i]])[2],'\n'))
	cat(paste0(pdtime(t0,2),'\n'))
}
reg.df = ldply(reg.li) # plyr
reg_1f_only = subset(reg.df,level%in%c("1f"))
reg_2b = subset(reg.df,level%in%c("1a","1b","1c","1d","1e","1f","2a","2b"))

# 2. Process random SNP files..
cat("\n(2/2) Process random SNP files..\n")
ids = c(1:10000)
snp_name  = paste0('rsid',ids,'.bed')
snp_path  = paste0('random/seeds/',snp_name)
snp_fname = paste0('regulome_rsid',ids,'.tsv')

n = length(snp_path); #dist.li = list()
pb = winProgressBar(title='Loading files',
	 label='Ready to read RDS files..',min=0,max=n,width=500)
dist.li=lapply(ids,function(i) {
    snp = read.delim(snp_path[i],header=F)
    colnames(snp) = c('chr','start','end','rsid')
    snpids = snp$rsid
    
    snpids.1f_only = subset(reg_1f_only,rsid%in%snpids)
    snpids.2b = subset(reg_2b,rsid%in%snpids)
    #dist.li[[i]] = data.frame(random=snp_fname[i],
    out = data.frame(random=snp_fname[i],
                              regulome_2b_n=nrow(snpids.2b),
                              func_Motifs=nrow(snpids.2b)-nrow(snpids.1f_only))
    
    f_name1 = paste0('random/regulome/regulome_rsid',i,'.tsv')
    write.table(snpids.2b,f_name1,row.names=F,col.names=T,quote=T,sep='\t')
    ## Progress time ##
    setWinProgressBar(pb,i,label=paste0(round(i/n*100,0),
              " % (",i,"/",n,") done for ",pdtime(t0,2)))
    ###################
    return(out)
})
dist.df = ldply(dist.li,data.frame) # plyr
close(pb)

f_name2 = paste0('random/regulome_dist.tsv')
write.table(dist.df,f_name2,row.names=F,col.names=T,quote=T,sep='\t')
cat(paste0('>> File write: ',f_name2,'\n'))
cat(paste0(pdtime(t0,1),'\n'))
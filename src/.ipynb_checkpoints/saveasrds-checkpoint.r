# Read file and save as rds file
library(data.table)
library(tools)
saveasrds = function(f_path) {
    df = fread(f_path,sep='\t',header=T,stringsAsFactors=F)
    f_name = paste0(f_path,'.rds')
    saveRDS(df,file=f_name)
    cat(paste0('>> File write: ',f_name,'\n'))
}

library(rtracklayer)
bedasrds = function(f_path) {
    df = import(f_path,format='bed')
    f_name = paste0(f_path,'.rds')
    saveRDS(df,file=f_name)
    cat(paste0('>> File write: ',f_name,'\n'))
}
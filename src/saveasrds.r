# Read file and save as rds file
suppressMessages(library(data.table))
suppressMessages(library(tools))
saveasrds = function(f_paths) {
    n = length(f_paths); df.li = list()
    for(i in 1:n) {
        if(file_ext(f_paths[i])=='gz') {
            try(R.utils::gunzip(f_paths[i]))
            f_path = file_path_sans_ext(f_paths[i])
            df.li  = fread(f_path,sep='\t',header=T,stringsAsFactors=F)
        } else {
            df.li  = fread(f_paths[i],sep='\t',header=T,stringsAsFactors=F)
        }
    }
    df = ldply(df.li) # plyr
    f_name = paste0(f_path,'.rds')
    saveRDS(df,file=f_name)
    cat(paste0('>> File write: ',f_name,'\n'))
}

suppressMessages(library(rtracklayer))
bedasrds = function(f_path) {
    df = as.data.frame(import(f_path,format='bed'))
    f_name = paste0(f_path,'.rds')
    saveRDS(df,file=f_name)
    cat(paste0('>> File write: ',f_name,'\n'))
}
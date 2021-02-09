#!/usr/bin/env Rscript
# This file is for identifying overlap numbers of SNPs between databases.

## Command Arg Parameters ##
# Usage 1: Rscript src/overlap_random.r random/roadmap random/encode random/overlap1
# Usage 2: Rscript src/overlap_random.r random/regulome random/overlap1 random/overlap2
args = commandArgs(trailingOnly=T)
hmsg = 'Rscript src/overlap_random.r [folder_path_1] [folder_path_2] [output_path]
  - [folder_path_1, _2] are mendatory arguments for the target folder path.
  - [output_path] is a mendatory argument for the file output path.'
if(length(args)<3|length(args)>3) stop(hmsg)

# System parameter
suppressMessages(library(data.table))
suppressMessages(library(plyr))
source('src/pdtime.r')
t0 = Sys.time()
paths = c(args[1],args[2])
dir_out = args[3]
dir.create(file.path(dir_out))

#########################
## Function start here ##
#########################
paths_li=list(); n=length(args)
for(i in 1:2) paths_li[[i]] = list.files(paths[i])

m = length(paths_li[[1]])
rd_li=lapply(c(1:m), function(j) {
    # 1. Read files
    df=list(); snp_ids=list()
    for(i in 1:2) {
        path = paste0(paths[i],'/',paths_li[[i]][j])
        if_ro = grep('roadmap',paths_li[[i]][j])
        if_en = grep('encode',paths_li[[i]][j])
        if_ov1= grep('overlap',paths_li[[i]][j])
        if_re = grep('regulome',paths_li[[i]][j])
        
        if(length(if_en)>0|length(if_ro)>0|length(if_ov1)>0) {
            df_= as.data.frame(fread(path),header=F)
            df[[i]]= df_[which(df_[,9]==0),]
            snp_ids[[i]] = df[[i]][,4]
        } else if(length(if_re)>0) {
            df_= as.data.frame(fread(path),header=T)
            snp_ids[[i]] = df_[,3]
        } else stop('There are unknown files..')
    }
    
    # 2. Overlap
    snp_inter = unique(intersect(snp_ids[[1]],snp_ids[[2]]))
    snp_inter_df = data.frame(rsid=snp_inter)
    
    f_name = paste0('overlap_rsid',j,'.tsv')
    f_path = paste0(dir_out,'/',f_name)
    
    # 3. Write result
    if_ro = grep('roadmap',paths_li[[2]][j])
    if_en = grep('encode',paths_li[[2]][j])
    if_ov1= grep('overlap',paths_li[[2]][j])
    if_re = grep('regulome',paths_li[[2]][j])
    if(length(if_en)>0|length(if_ro)>0|length(if_ov1)>0) {
        colnames(df[[2]]) = c('chr','start','end','rsid',
                              'enh_chr','enh_start','enh_end',
                              'num','dist')
        df_merge = merge(snp_inter_df,df[[2]],by='rsid')
        enh_pos_ = unique(unlist(apply(df_merge,1,function(x) paste0(x[5],':',x[6],'-',x[7]))))
        #print(paste0(paths_li[[2]][j],', snp_n = ',length(snp_inter),', enh_n = ',length(enh_pos_)))
        if(nrow(df_merge)>0) write.table(df_merge,f_path,sep='\t',col.names=F,row.names=F,quote=F)
        rd_row = c(f_name,length(snp_inter),length(enh_pos_))
    } else if(length(if_re)>0) {
        stop('Need to code here..')
    } else stop('Unknown file format..')
    
    if(j%%1000==0) cat(paste0('  - ',round((j/m)*100,1),'% Working process\n'))
    return(t(rd_row))
})
rd_df = ldply(rd_li,data.frame) # plyr
colnames(rd_df) = c('random','SNP_n','Enh_n')

f_path = paste0(dir_out,'_dist.tsv')
write.table(rd_df,f_path,row.names=F,quote=F,sep='\t')
cat(paste0('\nFile write: ',f_path,'\n'))
cat(paste0('\n>> ',pdtime(t0,1),'\n'))
##################
## Function end ##
##################
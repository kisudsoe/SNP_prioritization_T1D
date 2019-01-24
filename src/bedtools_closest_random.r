#!/usr/bin/env Rscript
# This file is for identifying T1D SNPs occupied in RoadMap enhancers

## Command Arg Parameters ##
# Rscript src/bedtools_closest_random.r roadmap
# Rscript src/bedtools_closest_random.r encode
args = commandArgs(trailingOnly=T)
hmsg = 'Rscript src/bedtools_closest_random.r [roadmap/encode]
  - Argument [roadmap/encode] is mendatory.
  - You should choose one of the roadmap/encode.'
if(length(args)<1|length(args)>1) stop(hmsg)

# System parameter
suppressMessages(library(data.table))
suppressMessages(library(plyr))
source('src/pdtime.r')
t0 = Sys.time()
rsid = c(1:10000)
if(args[1]=='roadmap') {
    dir = 'random/roadmap/'
    path = paste0(dir,'roadmap_rsid',rsid,'.tsv')
    f_path = paste0('random/roadmap_dist.tsv')
} else if(args[1]=='encode') {
    dir = 'random/encode/'
    path = paste0(dir,'encode_rsid',rsid,'.tsv')
    f_path = paste0('random/encode_dist.tsv')
} else stop(hmsg)

#########################
## Function start here ##
#########################
rd.li=list(); n=length(path)
rd.li = lapply(c(1:n),function(i) {
    # 1. Compiling SNP data distant from RoadMap enhancers
    rd = as.data.frame(fread(path[i])) # data.table
    #cat(paste0(' (',i,'/',n,')',path[i],'\n'))
    pos  = unlist(apply(rd,1,function(x) paste0(x[5],':',x[6],'-',x[7])))
    rd_  = cbind(rd[1:4],pos,rd[,8:9])
    colnames(rd_) = c('chr','start','end','rsid','enh_pos','num','dist')
    dis = as.numeric(as.character(rd_$dist))
    rd.df = cbind(rd_[1:6],dis); head(rd.df)
    
    # 2. Filtering SNPs by distance of closest enhancers
    rd.df_ = unique(subset(rd.df,dis==0))
    rd.df_snp = length(unique(rd.df_$rsid))
    rd.df_enh = length(unique(rd.df_$enh_pos))
    
    rd.row = c(basename(path[i]),rd.df_snp,rd.df_enh)
    if(i%%1000==0) cat(paste0('  - ',round((i/n)*100,1),'% Working process\n'))
    return(t(rd.row))
})
rd.df_ = ldply(rd.li,data.frame) # plyr
colnames(rd.df_) = c('random','SNP_n','Enh_n')
                  
write.table(rd.df_,f_path,row.names=F,quote=F,sep='\t')
cat(paste0('\nFile write: ',f_path,'\n'))
cat(paste0('\n>> ',pdtime(t0,1),'\n'))
##################
## Function end ##
##################
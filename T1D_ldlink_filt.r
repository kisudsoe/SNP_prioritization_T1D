#!/usr/bin/env Rscript
# This file is for filtering LDlink data

## Command Arg Parameters ##
# T1D_gwas.bat: Rscript T1D_ldlink_filt.r [SNP_file_path] [LDlink_download_target_dir]
args = commandArgs(trailingOnly=T)
<<<<<<< HEAD:T1D_ldlink_filt.r
if(length(args)<2) stop("Two arguments [SNP_file_path] [LDlink_download_target_dir] are mendatory.")
if(length(args)>2) stop("Too many arguments. Two arguments [SNP_file_path] [LDlink_download_target_dir] are needed.")
snp.path = unlist(args[1])
ld.path = unlist(args[2])

# System parameter
library(plyr)
source('src/pdtime.r')
=======
snp.path = unlist(args[1])
ld.path = unlist(args[2])

# System parameters
date = Sys.Date()
>>>>>>> parent of 0c366a1... write: db/roadmap - core functions:T1D_ldlink.r
t0 = Sys.time()


## Function start here ##
<<<<<<< HEAD:T1D_ldlink_filt.r
#########################
snpdf = unique(read.delim(snp.path)[1:2])
snpids = snpdf$rsid
cat(paste0('Input SNP list number = ',length(snpids),'\n\n'))
=======
snpids = unlist(read.table(snp.path))
cat(paste0('Input SNP list number = ',length(snpids),'\n'))
>>>>>>> parent of 0c366a1... write: db/roadmap - core functions:T1D_ldlink.r

ldlink = unlist(lapply(snpids, function(x) { paste0(ld.path,'/',x,'.tsv') }))
snptb  = data.frame(snpids=snpids, ldlink=ldlink)

ldlink.li = NULL # List variable for ldlink data
n = nrow(snptb); t1 = Sys.time()
pb = winProgressBar(title="Loop progress",
     label="Ready to read files..",min=0,max=n,width=500)
for(i in 1:n) {
    t = try(read.table(as.character(snptb[i,2]),header=T))
    if("try-error" %in% class(t)) ldlink[[i]] = NULL # get error from empty file
    else { # No errors
        tb1 = read.table(as.character(snptb[i,2]),header=T)
        tb2 = data.frame(SNPid=rep(snptb[i,1],nrow(tb1)),tb1)
        ldlink.li[[i]] = tb2
    }
    ## Progress time ##
    setWinProgressBar(pb,i,label=paste0(round(i/n*100,0),
                      " % (",i,"/",n,") done for ",
                      pdtime(t1,2)))
    ###################
}
<<<<<<< HEAD:T1D_ldlink_filt.r
close(pb)
ldlink.df = ldply(ldlink.li, data.frame) # plyr
ldlink_1 = subset(ldlink.df,R2>0.6 & Dprime==1)
ldlink_2 = unique(data.frame(gwasSNPs=ldlink_1$`SNPid`,
                             ldSNPs   =ldlink_1$`RS_Number`,
                             coord    =ldlink_1$`Coord`))
ldlink_ = ldlink_2[!ldlink_2$`ldSNPs` %in% c("."),] # Exclude no rsid elements
print(paste0("ldlink_ nrow = ",nrow(ldlink_)))
=======
close(pb); #print(length(ldlink.li))
ldlink.df = ldply(ldlink.li, data.frame)
#print(dim(ldlink.df))
ldlink_ = unique(subset(ldlink.df,R2>0.6 & Dprime==1)$RS_Number
>>>>>>> parent of 0c366a1... write: db/roadmap - core functions:T1D_ldlink.r

snp.t1 = unique(snpids); cat(paste0('SNP Tier1 = ',length(snp.t1),'\n'))
snp.t2 = setdiff(ldlink_, snp.t1)); cat(paste0('SNP Tier2 = ',length(snp.t2),'\n'))
snp.seed = union(snp.t1,snp.t2); cat(paste0('SNP seed  = ',length(snp.seed),'\n'))
pdtime(t0,1)

<<<<<<< HEAD:T1D_ldlink_filt.r
snpids_ = data.frame(gwasSNPs=snp.t1,ldSNPs=snpdf$rsid,coord=snpdf$coord) # coord data needed later.
snp.t2 = unique(setdiff(ldlink_$ldSNPs, snpids_$gwasSNPs))
n1 = length(snp.t2); cat(paste0('SNP Tier2 = ',n1,'\n'))

snp.seed = unique(rbind(snpids_,ldlink_))
n2 = length(unique(snp.seed$ldSNPs))
cat(paste0('SNP seed  = ',n2,'\n\n'))

cat(">> TSV files generation <<\n")
f.name1 = paste0('data/seedSNP_',n2,'_ldlink.tsv')
write.table(snp.seed,f.name1,row.names=F,quote=F,sep='\t')
cat(paste0('File write: ',f.name1,'\n'))

f.name2 = paste0('data/seedSNP_',n2,'.tsv')
write.table(unique(snp.seed$ldSNPs),f.name2,row.names=F,col.names=F,quote=F,sep='\t')
cat(paste0('File write: ',f.name2,'\n\n'))
cat(paste0('>> ',pdtime(t0,1),'\n'))

##################
## Function end ##
##################
=======
## Function end ##
>>>>>>> parent of 0c366a1... write: db/roadmap - core functions:T1D_ldlink.r

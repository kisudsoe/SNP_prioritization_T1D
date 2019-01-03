#!/usr/bin/env Rscript
# This file is for filtering LDlink data

## Command Arg Parameters ##
# gwas.bat: Rscript ldlink_filt.r [SNP_file_path] [LDlink_download_target_dir]
args = commandArgs(trailingOnly=T)
hmsg = 'Rscript T1D_ldlink_filt.r [SNP_file_path] [LDlink_download_target_dir]
  - Arguments [SNP_file_path] and [LDlink_download_target_dir] are mendatory.'
if(length(args)<2|length(args)>2) stop(hmsg)
snp.path = unlist(args[1])
ld.path = unlist(args[2])

# System parameter
library(plyr)
source('src/pdtime.r')
t0 = Sys.time()

#########################
## Function start here ##
#########################
snpdf = unique(read.delim(snp.path)[1:2]) # rsid, coord
snpids = snpdf$rsid
cat(paste0('Input SNP list number = ',length(snpids),'\n\n'))

ldlink = unlist(lapply(snpids, function(x) { paste0(ld.path,'/',x,'.tsv') }))
snptb  = data.frame(snpids=snpids, ldlink=ldlink)

n=nrow(snptb); ldlink.li=NULL # List variable for ldlink data
pb = winProgressBar(title="Loop progress",
     label="Ready to read files..",min=0,max=n,width=500)
for(i in 1:n) {
    tb1 = try(read.table(as.character(snptb[i,2]),header=T))
    if("try-error" %in% class(tb1)) ldlink.li[[i]] = NULL # get error from empty file
    else { # No errors then,
        tb2 = data.frame(SNPid=rep(snptb[i,1],nrow(tb1)),tb1)
        ldlink.li[[i]] = tb2
    }
    ## Progress time ##
    setWinProgressBar(pb,i,label=paste0(round(i/n*100,0),
              " % (",i,"/",n,") done for ",pdtime(t0,2)))
    ###################
}
close(pb)
ldlink.df = ldply(ldlink.li, data.frame) # plyr
ldlink_1 = subset(ldlink.df,R2>0.6 & Dprime==1)
ldlink_2 = unique(data.frame(gwasSNPs=ldlink_1$`SNPid`,
                             ldSNPs  =ldlink_1$`RS_Number`))
                             #coord   =ldlink_1$`Coord`))
ldlink_ = ldlink_2[!ldlink_2$`ldSNPs` %in% c("."),] # Exclude no rsid elements
ex = nrow(ldlink_2[ldlink_2$`ldSNPs` %in% c("."),])
cat(paste0('::Expluded no rsid elements = ',ex,'\n'))

# 1/3. Numbers of SNPs
cat(paste0('\n>> 1/3. Numbers of SNPs\n'))
snp.t1 = unique(snpids); cat(paste0('SNP Tier1 = ',length(snp.t1),'\n'))
snpids_ = data.frame(gwasSNPs=snp.t1,ldSNPs=snpdf$rsid)#,coord=snpdf$coord)
snp.t2 = unique(setdiff(ldlink_$ldSNPs, snpids_$gwasSNPs))
n1 = length(snp.t2); cat(paste0('SNP Tier2 = ',n1,'\n'))

snp.seed = unique(rbind(snpids_,ldlink_))
snp.coord= data.frame(ldSNPs=ldlink_1$`RS_Number`,coord=ldlink_1$`Coord`)
snp.seed_= merge(snp.seed,snp.coord,by='ldSNPs',all.x=T)[,c(2,1,3)]
n2 = length(unique(snp.seed_$ldSNPs))
cat(paste0('SNP seed  = ',n2,'\n'))

# 2/3. Generation of a result TSV file
cat("\n>> 2/3. Generation of a result TSV file\n")
f.name1 = paste0('data/seedSNP_',n2,'_ldlink.tsv')
write.table(snp.seed_,f.name1,row.names=F,quote=F,sep='\t')
cat(paste0('File write: ',f.name1,'\n'))

# 3/3. Generation of a result BED file
cat("\n>> 3/3. Generation of a result BED file\n")
coord = strsplit(as.character(snp.seed_$coord),':')
coord.df = as.data.frame(do.call(rbind,coord))
colnames(coord.df) = c('chr','pos')
start = as.numeric(as.character(coord.df$pos))-1
snp_bed = data.frame(chr  =coord.df$chr,
	                 start=start,
	                 end  =coord.df$pos,
	                 rsid =snp.seed_$ldSNPs)
unique(snp_bed)
cat(paste0('Table, rows= ',dim(snp_bed)[1],' cols= ',dim(snp_bed)[2],'\n'))
f.name2 = paste0('data/seedSNP_',n2,'.bed')
write.table(snp_bed,f.name2,row.names=F,col.names=F,quote=F,sep='\t')
cat(paste0('File write: ',f.name2,'\n'))
cat(paste0('>> ',pdtime(t0,1),'\n'))
##################
## Function end ##
##################

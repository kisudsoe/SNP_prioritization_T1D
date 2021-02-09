#!/usr/bin/env Rscript
# This file is for filtering LDlink data

## Command Arg Parameters ##
# ldlink.bat; Rscript ldlink_filt.r data/gwas_5e-08_129.tsv db/ldlink
args = commandArgs(trailingOnly=T)
hmsg = 'Rscript ldlink_filt.r [SNP_file_path] [LDlink_download_target_dir] [LDlink_filter_option]
  - Argument [SNP_file_path] is a mandatory for indicating root GWAS data file path.
  - Argument [LDlink_download_target_dir] is a mandatory for indicating LDlink data folder path.
  - Argument [LDlink_filter_option] is a mandatory. Choose one of the following option numbers.
    1) "r2>0.6 or Dprime=1"
    2) "r2>0.6"					<-This is usual choice to get LD associated SNPs.
    3) "Dprime=1"
    4) "r2>0.6 and Dprime=1"	<-This is the most stringent criteria.'
if(length(args)<3|length(args)>3) stop(hmsg)
snp_path    = unlist(args[1])
ld_path     = unlist(args[2])
filt_option = unlist(args[3])

# System parameter
library(plyr)
source('src/pdtime.r')
t0 = Sys.time()

#########################
## Function start here ##
#########################
snpdf = unique(read.delim(snp_path)[,1:2]) # rsid, coord
snpids = snpdf$rsid
cat(paste0('Input SNP list number = ',length(snpids),'\n\n'))

ldlink = unlist(lapply(snpids, function(x) { paste0(ld_path,'/',x,'.tsv') }))
snptb  = data.frame(snpids=snpids, ldlink=ldlink)

n=nrow(snptb); ldlink.li=NULL # List variable for ldlink data
pb = winProgressBar(title="Loop progress",
     label="Ready to read files..",min=0,max=n,width=500)
for(i in 1:n) {
    tb1 = try(read.table(as.character(snptb[i,2]),header=T))
    if("try-error" %in% class(tb1)) ldlink.li[[i]] = NULL # get error from empty file
    else { # If no errors,
        tb2 = data.frame(SNPid=rep(snptb[i,1],nrow(tb1)),tb1)
        ldlink.li[[i]] = tb2
    }
    ## Progress time ##
    setWinProgressBar(pb,i,label=paste0(round(i/n*100,0),
              " % (",i,"/",n,") done for ",pdtime(t0,2)))
    ###################
}
close(pb)
ldlink.df= ldply(ldlink.li, data.frame) # plyr
if(filt_option==1) {
	cat('Filtering option, r2 > 0.6 or Dprime = 1 was chosen.\n')
	ldlink_1 = subset(ldlink.df,R2>0.6 | Dprime==1) # r2 > 0.6 or D' = 1
} else if(filt_option==2) {
	cat('Filtering option, r2 > 0.6 was chosen.\n')
	ldlink_1 = subset(ldlink.df,R2>0.6) # r2 > 0.6
} else if(filt_option==3) {
	cat('Filtering option, Dprime = 1 was chosen.\n')
	ldlink_1 = subset(ldlink.df,Dprime==1) # D' = 1
} else if(filt_option==4) {
	cat('Filtering option, r2 > 0.6 and Dprime = 1 was chosen.\n')
	ldlink_1 = subset(ldlink.df,R2>0.6 & Dprime==1) # r2 > 0.6 and D' = 1
}
ldlink_2 = unique(data.frame(
	gwasSNPs=ldlink_1$`SNPid`,
    ldSNPs  =ldlink_1$`RS_Number`
    #coord   =ldlink_1$`Coord`
))
ldlink_ = ldlink_2[!ldlink_2$`ldSNPs` %in% c("."),] # Exclude no rsid elements
ex = nrow(ldlink_2[ldlink_2$`ldSNPs` %in% c("."),])
cat(paste0('::Excluded no rsid elements = ',ex,'\n'))

# 1/3. Numbers of SNPs
cat(paste0('\n1/3. Numbers of SNPs\n'))
snp.t1 = unique(snpids); cat(paste0('  SNP Tier1 = ',length(snp.t1),'\n'))
snpids_ = data.frame(gwasSNPs=snp.t1,ldSNPs=snpdf$rsid)#,coord=snpdf$coord)
snp.t2 = unique(setdiff(ldlink_$ldSNPs, snpids_$gwasSNPs))
n1 = length(snp.t2); cat(paste0('  SNP Tier2 = ',n1,'\n'))

snp.seed  = unique(rbind(snpids_,ldlink_))
snp.coord = data.frame(ldSNPs=ldlink_1$`RS_Number`,coord=ldlink_1$`Coord`)
snp.seed_ = merge(snp.seed,snp.coord,by='ldSNPs',all.x=T)[,c(2,1,3)]
snp.seed_ = unique(snp.seed_)
n2 = length(unique(snp.seed_$ldSNPs))
cat(paste0('  SNP seed  = ',n2,'\n'))

# 2/3. Generation of a result TSV file
cat("\n2/3. Generation of a result TSV file\n")
f.name1 = paste0('data/seedSNP_',n2,'_ldlink.tsv')
write.table(snp.seed_,f.name1,row.names=F,quote=F,sep='\t')
cat(paste0('  File write: ',f.name1,'\n'))

# 3/3. Generation of a result BED file
cat("\n3/3. Generation of a result BED file\n")
coord = strsplit(as.character(snp.seed_$coord),':')
coord.df = as.data.frame(do.call(rbind,coord))
colnames(coord.df) = c('chr','pos')
start = as.numeric(as.character(coord.df$pos))-1
snp_bed = data.frame(chr  =coord.df$chr,
	                 start=start,
	                 end  =coord.df$pos,
	                 rsid =snp.seed_$ldSNPs)
snp_bed = unique(snp_bed)
cat(paste0('  Table, rows= ',dim(snp_bed)[1],' cols= ',dim(snp_bed)[2],'\n'))
f.name2 = paste0('data/seedSNP_',n2,'.bed')
write.table(snp_bed,f.name2,row.names=F,col.names=F,quote=F,sep='\t')
cat(paste0('  File write: ',f.name2,'\n'))
cat(paste0(pdtime(t0,1),'\n'))
##################
## Function end ##
##################
#!/usr/bin/env Rscript
# This file is for selecting random SNP sets from dbSNP database.

## Command Arg Parameters ##
# CMD usage: Rscript src/rsid_random.r 1817 10000
args = commandArgs(trailingOnly=T)
hmsg = 'Rscript src/rsid_random.r [SNP_number] [SNP_set_number]
  - [SNP_number] is a mendatory numeric argument for a number of SNPs.
  - [SNP_set_number] is a mendatory numeric argument for iteration
  number of random SNP selection. Default is 1000.'
if(length(args)<2|length(args)>2) stop(hmsg)
n_snp = as.numeric(as.character(args[1]))
n_snp_set = as.numeric(as.character(args[2]))

# System parameter
suppressMessages(library(plyr))
suppressMessages(library(data.table))
source('src/pdtime.r')
source('src/saveasrds.r')
t0 = Sys.time()

# mkdir 'db/dbSNP' and 'random' folders
dir = 'db/dbSNP/'
dir.create(file.path(dir))
dir2 = 'random/'
dir.create(file.path(dir2))

#########################
## Function start here ##
#########################
# 1. Downloading dbSNP data
cat('\n(1/3) Downloading dbSNP data..\n')
chr = c(1:22,'X','Y','MT')
urls = paste0('ftp://ftp.ncbi.nih.gov/snp/organisms/human_9606_b151_GRCh37p13/BED/bed_chr_',chr,'.bed.gz')
t=lapply(urls,function(url) {
	#f_name1 = paste0(dir,basename(url))
	#cat(paste0('  - ',f_name1))
	#tb = try(download.file(url,destfile=f_name1,quiet=T))
	#if('try-error' %in% class(tb)) stop('dbSNP151 data address seems to be changed.')
    #cat(paste0(' > unzipping..'))
    #try(R.utils::gunzip(f_name1))
    #cat(paste0(' > converting into RDS..\n'))
    #f_name1_= tools::file_path_sans_ext(f_name1)
    #try(bedasrds(f_name1_))
})
cat(paste0('\n >> ',pdtime(t0,2),'\n'))

# 2. dbSNP BED file reading..
cat('\n(2/3) Reading dbSNP data RDS files..\n')
f_name2 = paste0(dir,'bed_chr_',chr,'.bed.rds')
n=length(f_name2)
snp.li = lapply(c(1:n),function(i) {
    f = f_name2[i]
    if(file.exists(f)) {
        tb = try(as.data.frame(readRDS(f)))
        if('try-error' %in% class(tb)) stop('There is no dbSNP RDS files before this process.') # check no file error
        cat(paste0('\t(',i,'/',n,') ',
               basename(f),', rows= ',dim(tb)[1],' cols= ',dim(tb)[2],', ',
               pdtime(t0,2),'\n'))
        return(tb)
    } else cat(paste0('None such file: ',f))
})
snp.df = ldply(snp.li,data.frame)[,c(1:3,6)] # plyr
cat(paste0('\n >> SNP table, rows= ',dim(snp.df)[1],' cols= ',dim(snp.df)[2],'\n'))
cat(paste0(' >> ',pdtime(t0,2),'\n'))

# 3. Generating random-sets of SNPs
cat('\n(3/3) Generating random-sets of SNPs\n')
cat(paste0('  - SNP_number = ',n_snp,'\n  - SNP_set_number = ',n_snp_set,'\n'))
t=lapply(c(1:n_snp_set),function(i) {
	rd.row = sample(nrow(snp.df),n_snp)
	rd.df  = snp.df[rd.row,]
    #rd.df  = rd.df[sort(rownames(rd.df)),]
    m=nrow(rd.df)
    rd.li=lapply(c(1:m),function(j) {
        row   = t(unname(rd.df[j,]))
        row_2 = as.numeric(row[2])
        row_3 = as.numeric(row[3])
        if(row_3-row_2<0) row_ = row[c(1,3,2,4)]
        else row_ = row
        return(t(row_))
    })
    rd.df_ = ldply(rd.li,data.frame)
    f_name = paste0(dir2,'rsid',i,'.bed')
	write.table(rd.df_,f_name,row.names=F,col.names=F,quote=F,sep='\t')
    #cat(paste0(' >> File write: ',f_name,'\n'))
})
cat(paste0('\n >> ',pdtime(t0,1),'\n'))
##################
## Function end ##
##################
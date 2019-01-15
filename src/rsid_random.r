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
library(rtracklayer)
library(plyr)
source('src/pdtime.r')
t0 = Sys.time()

# mkdir 'db/dbSNP' and 'random' folders
dir = 'db/dbSNP/'
if(file.exists(dir)) { cat(paste0('Directory exists: ',dir,'\n'))
} else {
	dir.create(file.path(dir))
	cat(paste0('Directory generated: ',dir,'\n'))
}
dir2 = 'random/'
if(file.exists(dir2)) { cat(paste0('Directory exists: ',dir2,'\n'))
} else {
	dir.create(file.path(dir2))
	cat(paste0('Directory generated: ',dir2,'\n'))
}

#########################
## Function start here ##
#########################
# 1. Downloading dbSNP data
cat('\n(1/3) Downloading dbSNP data\n')
chr = c(1:22,'X','Y','MT')
urls = paste0('ftp://ftp.ncbi.nih.gov/snp/organisms/human_9606_b151_GRCh37p13/BED/bed_chr_',chr,'.bed.gz')
t=lapply(urls,function(url) {
	f.name = paste0(dir,basename(url))
	#cat(paste0('\t',f.name,'\n'))
	#tb = try(download.file(url,destfile=f.name,quiet=T))
	#if('try-error' %in% class(tb)) stop('dbSNP151 data address seems to be changed.')
})
cat(paste0('\n>> ',pdtime(t0,2),'\n'))

# 2. Loading dbSNP data
cat('\n(2/3) Loading dbSNP data\n')
f.name = paste0(dir,'bed_chr_',chr,'.bed')
snp.li=list(); n=length(f.name)
pb = winProgressBar(title=paste0("Loop progress start at ",Sys.time()),
     label="Ready to read files..",min=0,max=n,width=500)
cat('  - dbSNP BED file reading...\n')
for(i in 1:n) {
	if(file.exists(f.name[i])) {
		tb = try(as.data.frame(import(f.name[i],format='bed'))) # rtracklayer
		if('try-error' %in% class(tb)) stop('You should decompress gz files before this process.') # get error from no file
		snp.li[[i]] = tb
	} else { cat(paste0('None such file: ',f.name[i])) }
	cat(paste0('  (',i,'/',n,')',
		basename(f.name[i]),', rows= ',dim(tb)[1],' cols= ',dim(tb)[2],', ',
				 pdtime(t0,2),'\n'))
	## Progress time ##
    setWinProgressBar(pb,i,label=paste0(round(i/n*100,0),
              " % (",i,"/",n,") done for ",pdtime(t0,2)))
    ###################
}
close(pb)
snp.df = ldply(snp.li,data.frame)[,c(1:3,6)] # plyr
cat(paste0('\n >> SNP table, rows= ',dim(snp.df)[1],' cols= ',dim(snp.df)[2],'\n'))
cat(paste0(' >> ',pdtime(t0,2),'\n'))

# 3. Generating random-sets of SNPs
cat('\n(3/3) Generating random-sets of SNPs\n')
cat(paste0('  - SNP_number = ',n_snp,
		 '\n  - SNP_set_number = ',n_snp_set,'\n'))
num_snp_set = c(1:n_snp_set)
t=lapply(num_snp_set,function(i) {
	rd.row = sample(nrow(snp.df),n_snp)
	rd.df  = snp.df[rd.row,]
	write.table(rd.df,paste0(dir2,'rsid',i,'.bed'),
				row.names=F,col.names=F,quote=F,sep='\t')
})
cat(paste0('\n>> ',pdtime(t0,1),'\n'))
##################
## Function end ##
##################
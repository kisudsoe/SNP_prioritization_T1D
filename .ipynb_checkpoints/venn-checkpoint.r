#!/usr/bin/env Rscript
# This file is for venn analysis of SNP lists

## Command Arg Parameters ##
# CMD command1: Rscript venn.r data/seedSNP_1817.bed data/snp_484_roadmap_dist.bed data/snp_364_encode_dist.bed data/snp_94_regulome2b.bed
# CMD command2: Rscript venn.r data/seedSNP_1817.bed data/snp_140_roadmap_encode.bed data/snp_26_core.bed data/snp_745_gtex.bed
args = commandArgs(trailingOnly=T)
hmsg = 'Rscript regulome.r [SNP_BED_file_1] [SNP_BED_file_2] [SNP_BED_file_3] [...]
  - Arguments [SNP_BED_file_1], [SNP_BED_file_2], and so on are needed.'
n = length(args)
if(n < 1) stop(hmsg)
f.path = NULL
for(i in 1:n) {
	f.path = c(f.path,unlist(args[i]))
}

# System parameter
library(tools)
library(limma)
source('src/pdtime.r')
t0 = Sys.time()

#########################
## Function start here ##
#########################
snp.li=list(); snpids.li=list(); names=NULL
for(i in 1:n) {
	tb = read.delim(f.path[i],header=F)
	colnames(tb) = c('chr','start','end','ann')
	snp.li[[i]] = tb
	snpids.li[[i]] = tb$ann
	names = c(names,file_path_sans_ext(basename(f.path[i])))
	print(file_path_sans_ext(basename(f.path[i])))
}

unionlist = Reduce(union,snpids.li)
# Collapse the list to one dataFrame list
unionPr=NULL; subtitle=NULL
for(i in 1:length(snpids.li)) {
	unionPr = cbind(unionPr,snpids.li[[i]][match(unionlist,snpids.li[[i]])])
	subtitle= c(subtitle,paste0(names[i],'\n',length(snpids.li[[i]])))
}
rownames(unionPr) = unionlist

# Generate binary table to match with ID list
union = (unionPr != '')			# Transcform values to TRUE, if ID exists.
union[is.na(union)] = FALSE		# Transform NA to FALSE value
union = as.data.frame(union)	# Make 'union' to data.frame from
colnames(union) = subtitle			# Names attach to venn diagram

# Draw Venn Diagram
f.name1a= paste0('venn_',names[1],'_',names[2],'.png')
f.name1 = paste0('fig/',f.name1a)
title1 = paste0('Venn analysis of ',nrow(union),' SNPs')
png(f.name1,width=10,height=10,units='in',res=100)
vennDiagram(union[2:4], main=title1, circle.col=rainbow(length(subtitle))) # limma
dev.off()
cat(paste0('\nFigure draw: ',f.name1,'\n'))

# Write files
colnames(union) = names
union.df = cbind(union[2:4],ann=rownames(union))
union.out= merge(snp.li[[1]],union.df,by="ann")
venn.li  = list(list=union.out, vennCounts=vennCounts(union[2:4]))
f.name2  = paste0('data/venn.tsv')
write.table(venn.li[[1]],f.name2,row.names=F,col.names=T,quote=F,sep='\t')
cat(paste0('File write: ',f.name2,'\n'))

f.name3 = paste0('data/vennCounts.tsv')
write.table(venn.li[[2]],f.name3,row.names=F,col.names=T,quote=F,sep='\t')
cat(paste0('File write: ',f.name3,'\n'))

core.df = union.out[(union.out[,5]==TRUE&
					 union.out[,6]==TRUE&
					 union.out[,7]==TRUE),c(2:4,1)]
f.name4 = paste0('data/snp_',nrow(core.df),'_core.bed')
write.table(core.df,f.name4,row.names=F,col.names=F,quote=F,sep='\t')
cat(paste0('File write: ',f.name4,'\n'))
cat(paste0('\n>> ',pdtime(t0,1),'\n'))
##################
## Function end ##
##################
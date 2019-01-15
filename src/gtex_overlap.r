#!/usr/bin/env Rscript
# CMD: Rscript src/gtex_overlap.r
path1 = 'data/snp_140_roadmap_encode.bed'
path2 = 'data/gtex_5e-08_745.tsv'
path3 = 'data/gtex_nearest_df.tsv'

cat(paste0('(1/3) Read files..\n'))
snp_140   = read.delim(path1,head=F)
snpids    = gsub('(.*)_.*','\\1',snp_140[,4])
snp_140_  = data.frame(snp_140,rsid=snpids)
cat(paste0(' - ',path1,', rows= ',dim(snp_140)[1],' cols= ',dim(snp_140)[2],'\n'))
gtex      = read.delim(path2)
cat(paste0(' - ',path2,', rows= ',dim(gtex)[1],' cols= ',dim(gtex)[2],'\n'))
gtex_rsid = length(levels(gtex$rsid)); gtex_gid = length(levels(gtex$gene_id))
cat(paste0(' - ',path2,', SNPs= ',gtex_rsid,' Genes= ',gtex_gid,'\n'))
gtex_nearest = read.delim(path3)
작성중..

cat(paste0('\n(2/3) Overlap these two files..\n'))
gtex_ = subset(gtex,rsid%in%snpids)
cat(paste0(' - TFBS overlap, rows= ',dim(gtex_)[1],' cols= ',dim(gtex_)[2],'\n'))
gtex_rsid_ = length(unique(gtex_$rsid)); gtex_gid_ = length(unique(gtex_$gene_id))
cat(paste0(' - TFBS overlap, SNPs= ',gtex_rsid_,' Genes= ',gtex_gid_,'\n'))
gtex_bl = subset(gtex_,tissue=='Whole_Blood')

cat(paste0('\n - Whole_Blood, rows= ',dim(gtex_bl)[1],' cols= ',dim(gtex_bl)[2],'\n'))
gtex_bl_rsid_ = length(unique(gtex_bl$rsid)); gtex_bl_gid_ = length(unique(gtex_bl$gene_id))
cat(paste0(' - Whole_Blood, SNPs= ',gtex_bl_rsid_,' Genes= ',gtex_bl_gid_,'\n'))

snp_74 = subset(snp_140_,rsid%in%gtex_$rsid)[1:4]
f.name = paste0('data/snp_',gtex_rsid_,'_gtex_enh.bed')
write.table(snp_74,f.name,row.names=F,col.names=F,quote=T,sep='\t')
cat(paste0('\n>> File write: ',f.name,'\n'))

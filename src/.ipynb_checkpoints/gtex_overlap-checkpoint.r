#!/usr/bin/env Rscript
# CMD: Rscript src/gtex_overlap.r
path1 = 'data/snp_140_roadmap_encode.bed'
path2 = 'data/gtex_5e-08_745.tsv'
path3 = 'data/gtex_nearest_df.tsv'

# 1. Read files..
cat(paste0('(1/2) Read files..\n'))
snp_140   = read.delim(path1,head=F)
colnames(snp_140) = c('chr','start','end','ann')
cat(paste0('  - ',path1,', rows= ',dim(snp_140)[1],' cols= ',dim(snp_140)[2],'\n'))
snpids    = gsub('(.*)_.*','\\1',snp_140[,4])
snp_140_  = cbind(snp_140,snpids)

gtex      = read.delim(path2)
cat(paste0('  - ',path2,', rows= ',dim(gtex)[1],' cols= ',dim(gtex)[2],'\n'))
ENSGids   = unlist(lapply(as.character(gtex$gene_id),function(x) {
	unlist(strsplit(x,'.',fixed=T))[1]
}))
gtex_ENSGid= cbind(gtex,ENSGids)
gtex_rsid = length(unique(gtex$rsid)); gtex_gid = length(unique(ENSGids))

nearest   = read.delim(path3)
cat(paste0('  - ',path3,', rows= ',dim(nearest)[1],' cols= ',dim(nearest)[2],'\n'))
near_ENSGids = intersect(ENSGids,unique(nearest$ENSGid))
near_gid = length(unique(near_ENSGids))
cat(paste0(' >> SNPs= ',gtex_rsid,' Genes= ',gtex_gid,' (Nearest= ',near_gid,')\n'))

# 2. Overlap these two files..
cat(paste0('\n(2/2) Overlap these two files..\n'))
gtex_ = subset(gtex_ENSGid,rsid%in%snpids)
cat(paste0('  - TFBS overlap, rows= ',dim(gtex_)[1],' cols= ',dim(gtex_)[2],'\n'))
gtex_rsid_ = length(unique(gtex_$rsid)); gtex_gid_ = length(unique(gtex_$ENSGid))
gtex_near_gid = length(unique(intersect(gtex_$ENSGids,near_ENSGids)))
cat(paste0(' >> SNPs= ',gtex_rsid_,' Genes= ',gtex_gid_,
		   ' (Nearest= ',gtex_near_gid,')\n'))

gtex_bl = subset(gtex_,tissue=='Whole_Blood')
cat(paste0('\n  - Whole_Blood, rows= ',dim(gtex_bl)[1],' cols= ',dim(gtex_bl)[2],'\n'))
gtex_bl_rsid_ = length(unique(gtex_bl$rsid)); gtex_bl_gid_ = length(unique(gtex_bl$ENSGid))
gtex_near_bl_gid  = length(unique(intersect(gtex_bl$ENSGid,near_ENSGids)))
cat(paste0(' >> SNPs= ',gtex_bl_rsid_,' Genes= ',gtex_bl_gid_,
		   ' (Nearest= ',gtex_near_bl_gid,')\n'))

snp_74 = subset(snp_140_,snpids%in%unique(gtex$rsid))[1:4]
f.name = paste0('data/snp_',nrow(snp_74),'_gtex_enh.bed')
write.table(snp_74,f.name,row.names=F,col.names=F,quote=T,sep='\t')
cat(paste0('\n>> File write: ',f.name,'\n'))

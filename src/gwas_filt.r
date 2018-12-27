gwas_filt = function(path,pval=5*10^-8) {
    cat(paste0('  - Criteria= ',pval,'\n'))
    cat(paste0('  - Input file: ',path,'\n'))
    snp = read.table(path,header=T,sep='\t')
    cat(paste0('  - Table rows= ',dim(snp)[1],', cols= ',dim(snp)[2],'\n'))
    snp_ = subset(snp,`P.VALUE`<pval)
    snp_id = unique(snp_$`SNPS`)
    cat(paste0('  - Unique IDs (p<5E-8)= ',length(snp_id),'; ','\n\n'))
    
    write.table(snp_id,paste0('SNP_',pval,'_',length(snp_id),'.tsv'),
                quote=F,row.names=F,col.names=F)
}
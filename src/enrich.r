suppressMessages(library(argparser))

## Parsing Arguments ##
p = arg_parser("Function for enrichment analysis (Fisher/Permutation)")

### Examples
p = add_argument(p, '--example', help="See examples to run this tool.", flag=T)
example_msg = '
Rscript src/enrich.r --roadmap --dir db/roadmap_dist --out db
Rscript src/enrich.r --permu --gwassnp --enhsnp --groups --permn --out
'

### Shared Arguments
p = add_argument(p, '--out',     help="[Path] Target directory for output files.")

### Arguments for roadmap_merge function
p = add_argument(p, '--roadmap', help="Filter and merge Roadmap dist restuls", flag=T)
p = add_argument(p, '--dir',     help="[Path] Driectory including roadmap dist results.")

### Arguments for perm_test function
p = add_argument(p, '--permu',   help="Run fisher test", flag=T)
p = add_argument(p, '--gwassnp', help="[Path] GWAS list file. Columns: <SNPS> <MAPPED_TRAIT> <P.VALUE>")
p = add_argument(p, '--enhsnp',  help="[Path] BED file including enhancer residing SNP-GWAS/PheWAS traits for fisher and permtest functions.")
p = add_argument(p, '--permn',   help="[Number] Set permutation number. Default=500")
p = add_argument(p, '--groups',  help="[Path] TXT file including age-related trait list.")

argv = parse_args(p)


## Functions ##
suppressMessages(library(dplyr))
suppressMessages(library(tidyr))

perm_test = function(
    f_gwas_snp = NULL,
    f_enh_snp  = NULL,
    out        = 'perm_test',
    perm_n     = 5000,
    f_groups   = NULL
) {
    paste0('\n** Run perm_test function in enrich.r **\n\n') %>% cat
    source('src/pdtime.r')
    suppressMessages(library(regioneR))
    
    # Read files
    gwas_snp = read.delim(f_gwas_snp,header=F,stringsAsFactors=F)
    colnames(gwas_snp) = c('Chr','Start','End','snp-phenotype')
    paste0('* gwas snp = ') %>% cat; dim(gwas_snp) %>% print

    enh_snp = read.delim(f_enh_snp,stringsAsFactors=F)
    paste0('* enh_snp = ') %>% cat; dim(enh_snp) %>% print
    
    # 1-1. Extract Rsids from gwas_snp
    t1=Sys.time()
    paste0('* Extract Rsids:\n') %>% cat
    gwas_rsid = lapply(gwas_snp$`snp-phenotype`,function(x) {
        row = strsplit(x,"\\_")[[1]]
        pheno = paste0(row[2],"_",row[3])
        return(data.frame(Rsid=row[1],Phenotype=pheno))
    })
    gwas_rsid_df1 = data.table::rbindlist(gwas_rsid)
    gwas_rsid_df = data.frame(
        gwas_snp[,1:3],
        gwas_rsid_df1
    )
    paste0('  gwas_snp = ') %>% cat
    dim(gwas_rsid_df) %>% print
    paste0(pdtime(t1,2),'\n') %>% cat

    # 1-2. Extract Rsids from enh_snp
    t1=Sys.time()
    enh_snp_li = apply(enh_snp[,c(4,8)],1,function(x) {
        row1 = strsplit(x[1],'\\_')[[1]]
        row2 = strsplit(x[2],'\\.')[[1]]
        pheno = paste0(row1[2],'_',row1[3])
        return(data.frame(Rsid=row1[1],Phenotype=pheno,Group=row2[1],Enh_peak=row2[2]))
    })
    enh_snp_df1 = data.table::rbindlist(enh_snp_li)
    enh_snp_df = data.frame(
        enh_snp[,1:3],
        enh_snp_df1
    )
    enh_snp_df = data.frame(lapply(enh_snp_df,as.character),stringsAsFactors=F)
    paste0('  enh_snp = ') %>% cat
    dim(enh_snp_df) %>% print
    paste0(pdtime(t1,2),'\n') %>% cat
    
    # 2. Extract rsids by groups and phenotypes as BED format
    t1=Sys.time()
    groups = enh_snp_df$Group %>% unique
    n = length(groups)
    
    # Read age-related trait list
    age = read.delim(f_groups,stringsAsFactors=F) %>% unlist

    ## Input all gwas rsid set
    all_gwas = toGRanges(gwas_rsid_df[,1:4],format="BED")

    permPval_df = NULL
    for(i in 37:n) { # by group
        t2=Sys.time()
        paste0('  ',i,'/',n,' ',groups[i],'... ') %>% cat
        ## Enhancer Rsids by group
        enh_snp_df_sub = subset(enh_snp_df,Group==groups[i])
        phenos1 = enh_snp_df_sub$Phenotype %>% unique
        phenos  = intersect(phenos1, age)

        m = length(phenos)
        paste0(m,'.. ') %>% cat
        pval = lapply(c(1:m),function(j) { # by phenotype
            ## Input enhancer rsid set and gwas phenotype rsid set
            enh_snp_df_sub_sub = subset(enh_snp_df_sub,Phenotype==phenos[j])
            enhancer = toGRanges(enh_snp_df_sub_sub[,1:4],format="BED")
            gwas_rsid_df_sub = subset(gwas_rsid_df,Phenotype==phenos[j])
            gwas_pheno = toGRanges(gwas_rsid_df_sub[,1:4],format="BED")

            ## Perform a permutation test
            perm_n = as.numeric(perm_n)
            pt = permTest(
                A                  = enhancer,
                ntime              = perm_n,
                randomize.function = resampleRegions,
                universe           = all_gwas,
                evaluate.function  = numOverlaps,
                B                  = gwas_pheno,
                verbose=F )
            pt$numOverlaps$pval
        }) %>% unlist
        paste0(length(pval),'.. ') %>% cat

        ## Collecting the results as df
        permPval_df1 = data.frame(
            Group     = rep(groups[i],m),
            Phenotype = phenos,
            Pvalue    = pval
        )
        permPval_df = rbind(permPval_df,permPval_df1)
        paste0(pdtime(t2,2),'\n') %>% cat
    }
    permPval_df$fdr = p.adjust(permPval_df$Pvalue,method='BH') # Add FDR
    permPval_df$bonf = p.adjust(permPval_df$Pvalue,method='bonferroni') # Add Bonferroni
    paste0('  permPval dim = ') %>% cat
    dim(permPval_df) %>% print
    paste0(pdtime(t1,2),'\n') %>% cat

    # 3. Save as TSV file
    f_name = paste0(out,'.tsv')
    write.table(permPval_df,f_name,row.names=F,sep='\t',quote=F)
    paste0('Write file: ',f_name,'\n') %>% cat
    return(permPval_df)
}


read_roadmap_file = function(path) {
    roadmap_dist = read.delim(path, header=F)
    roadmap_dist_sub = roadmap_dist[,c(1:4,8,14)]
    colnames(roadmap_dist_sub) = c('Chr','Start','End','Rsid','State','Dist')
    roadmap_out = subset(roadmap_dist_sub, Dist==0) # Filter by dist=0
    return(roadmap_out)
}

roadmap_merge = function(
    f_dir = NULL,
    out   = 'db'
) {
    paste0('\n** Run roadmap_merge function in enrich.r **\n\n') %>% cat
    file_paths = list.files(f_dir, full.names=T)
    n = length(file_paths)

    # Read files
    err_count = 0
    paste0('Read roadmap dist files... ',n,'.. ') %>% cat
    road_li = lapply(c(1:n),function(i) {
        cell_type = strsplit(file_paths[i] %>% basename,'_')[[1]][1]
        roadmap   = read_roadmap_file(file_paths[i]) # Read dist file and filter by dist=0
        roadmap1  = data.frame(
            roadmap, Cell=rep(cell_type,nrow(roadmap))
        ) %>% unique
        return(roadmap1)
    })
    road_df = data.table::rbindlist(road_li)
    dim(road_df) %>% print

    # Write result as a file
    f_name = paste0(out,'/roadmap_merge.tsv')
    write.table(road_df,f_name,sep='\t',row.names=F,quote=F)
    paste0('Write TSV file: ',f_name,'\n') %>% cat
}

## Functions End ##

## Run Function ##
source('src/pdtime.r')
t0 = Sys.time()

if(argv$example) {
    cat(example_msg)
} else if(argv$roadmap) {
    roadmap_merge(
        f_dir = argv$dir,
        out   = argv$out
    )
} else if(argv$permu) {
    perm_test(
        f_gwas_snp = argv$gwassnp,
        f_enh_snp  = argv$enhsnp,
        out        = argv$out,
        perm_n     = argv$permn,
        f_groups   = argv$groups
    )
}
paste0('\n',pdtime(t0,1),'\n') %>% cat
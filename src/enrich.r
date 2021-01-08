suppressMessages(library(argparser))

## Parsing Arguments ##
p = arg_parser("Function for enrichment analysis (Fisher/Permutation)")

### Examples
p = add_argument(p, '--example', help="See examples to run this tool.", flag=T)
example_msg = '
Rscript src/enrich.r --permu --verbose \
    --gwassnp data/seedSNP_1817_bm.bed \
    --chrstatus db/roadmap_bed/E001_25_imputed12marks_dense.bed \
    --dbsource roadmap_bed \
    --permn 5000 \
    --out enrich

Rscript src/enrich.r --permu \
    --gwassnp data/seedSNP_1817_bm.bed \
    --chrstatus db/roadmap_bed \
    --dbsource roadmap_bed \
    --permn 5000 \
    --out enrich
'

### Shared Arguments
p = add_argument(p, '--out',       help="[Path] Target directory for output files.")
p = add_argument(p, '--verbose',   help="[Flag] Show detailed process.", flag=T)

### Arguments for perm_test function
p = add_argument(p, '--permu',     help="[Flag] Run fisher test", flag=T)
p = add_argument(p, '--gwassnp',   help="[Path] GWAS list file. Columns: <SNPS> <MAPPED_TRAIT> <P.VALUE>")
p = add_argument(p, '--chrstatus', help="[Path] BED file including chromosome status for permtest functions.")
p = add_argument(p, '--dbsource',  help="[roadmap_bed/encode] Type in one of these options for source of your chromosome status.")
p = add_argument(p, '--permn',     help="[Number] Set permutation number. Default=500", default=500, type="numeric")

### Arguments for draw_heatmap function
p = add_argument(p, '--heatmap',     help="[Flag] Draw a heatmap from permu results.", flag=T)
p = add_argument(p, '--zscoretb',   help="[Path] Z-score table file.")

argv = parse_args(p)

## Load Common Libraries ##

suppressMessages(library(dplyr))
suppressMessages(library(tidyr))

## Functions ##
draw_heatmap = function(
    f_zscore = NULL,
    out      = 'enrich'
) {}


perm_test = function(
    f_gwas_snp = NULL,
    f_status   = NULL,
    db_src     = NULL,
    out        = 'enrich',
    perm_n     = 5000,
    verbose    = NULL
) {
    paste0('\n** Run perm_test function in enrich.r **\n\n') %>% cat
    suppressMessages(library(regioneR))
    ifelse(!dir.exists(out), dir.create(out), "")

    # Read file
    paste0('* Gwas snp = ') %>% cat
    gwas_snp = read.delim(f_gwas_snp, header=F)
    colnames(gwas_snp) = c('Chr','Start','End','Rsid')
    gwas_snp_bed = toGRanges(gwas_snp[,1:4], format="BED")
    dim(gwas_snp) %>% print

    # Run perm_test for each status files
    f_status_paths = list.files(f_status, full.names=T, include.dirs=T)
    n = length(f_status_paths)
    if(n==0) { f_status_paths = f_status; n=1 }
    paste0('* ',n,' files were found from ',f_status,'.\n\n') %>% cat
    cell_types=NULL; zscore_li=list(); pval_li=list()
    source('src/pdtime.r')
    for(i in 1:n) {
        t1=Sys.time()
        file_name = basename(f_status_paths[i])
        cell_type = strsplit(file_name,'_')[[1]][1] # Roadmap file-specific
        cell_types = c(cell_types,cell_type)

        ## Read file
        paste0(i,' Load ',cell_type,': ') %>% cat
        status = read_status_file(f_status_paths[i],db_src)
        dim(status) %>% print

        ## Run perm_test
        pt_df = perm_test_calc(gwas_snp_bed,status,perm_n,verbose)
        zscore_li[[i]] = data.frame(Status=pt_df$Status, Zscore=pt_df$Zscore)
        pval_li[[i]]   = data.frame(Status=pt_df$Status, Pval  =pt_df$Pval)
        paste0(pdtime(t1,2),'\n') %>% cat
    }
    pm_merge = function(x,y) {
        merge(x=x,y=y,by='Status',all=T)
    }
    zscore_df = Reduce(pm_merge, zscore_li)
    pval_df   = Reduce(pm_merge, pval_li)
    colnames(zscore_df) = c('Status',cell_types)
    colnames(pval_df)   = c('Status',cell_types)
    
    # Write result as files
    gwas_base = basename(f_gwas_snp)
    gwas_f_name = tools::file_path_sans_ext(gwas_base)
    
    f_name1 = paste0(out,'/',db_src,'-',gwas_f_name,'-zscore.tsv')
    write.table(zscore_df,f_name1,sep='\t',row.names=F,quote=F)
    paste0('\n* Write file: ',f_name1,'\n') %>% cat
    
    f_name2 = paste0(out,'/',db_src,'-',gwas_f_name,'-pval.tsv')
    write.table(pval_df,f_name2,sep='\t',row.names=F,quote=F)
    paste0('* Write file: ',f_name2,'\n') %>% cat
}

read_status_file = function(path,db_src) {
    paste0(db_src,' = ') %>% cat
    if(db_src=='roadmap_bed') {
        status_raw = read.delim(path,header=F,skip=1)
        status_raw = status_raw[,1:4]
        colnames(status_raw) = c('Chr','Start','End','Ann')
    } else if(db_src=='encode') {
        paste0('encode = ') %>% cat
    } else {
        paste0('\n\n[Error] Unknown option is flagged. Please choose one of [roadmap_bed/encode].\n') %>% cat
        stop()
    }
    return(status_raw)
}

perm_test_calc = function(
    gwas_snp_bed = NULL,
    status   = NULL,
    perm_n   = 5000,
    verbose  = FALSE
) {
    #paste0('\n** Run perm_test function in enrich.r **\n\n') %>% cat
    #suppressMessages(library(regioneR))

    # Calculate permTest to get z-scores and p-values by chromosome status
    paste0('  Run permTest: ') %>% cat
    status_ann = status$Ann %>% unique %>% sort
    paste0(length(status_ann),' annotations, [') %>% cat
    status_all_bed = toGRanges(status[,1:4],format="BED")
    pt_li = lapply(status_ann,function(x) {
        '.' %>% cat
        # Convert data.frame to GRanges format
        status_sub = subset(status,Ann==x)
        status_sub_bed = toGRanges(status_sub,format="BED")
        
        pt = permTest(
            A                  = status_sub_bed,
            B                  = gwas_snp_bed,
            ntime              = perm_n,
            randomize.function = resampleRegions,
            universe           = status_all_bed,
            evaluate.function  = numOverlaps,
            force.parallel     = T,
            verbose            = verbose
        )
        #pt$numOverlaps$pval
        pt_df = data.frame(Status=x,Zscore=pt$numOverlaps$zscore,Pval=pt$numOverlaps$pval)
        return(pt_df)
    })
    paste0('] done. ') %>% cat
    pt_df = data.table::rbindlist(pt_li)
    return(pt_df)
}

## Functions End ##

## Run Function ##
source('src/pdtime.r')
t0 = Sys.time()

if(argv$example) {
    cat(example_msg)
} else if(argv$permu) {
    perm_test(
        f_gwas_snp = argv$gwassnp,
        f_status   = argv$chrstatus,
        db_src     = argv$dbsource,
        perm_n     = argv$permn,
        out        = argv$out,
        verbose    = argv$verbose
    )
} else if(argv$heatmap) {
    draw_heatmap()
}
paste0('\n',pdtime(t0,1),'\n') %>% cat
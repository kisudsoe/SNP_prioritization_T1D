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

Rscript src/enrich.r --heatmap \
    --pmdata enrich/roadmap_bed-snp_484_roadmap_dist-permn_100-zscore.tsv \
    --meta db/roadmap_meta.tsv \
    --out enrich \
    --range -3,3 \
    --annot BLOOD,PANCREAS,THYMUS \
    --fileext png

Rscript src/enrich.r --splittfbs \
    --tfbs db/wgEncodeRegTfbsClusteredWithCellsV3.bed \
    --out db/tfbs_cell
'

### Shared Arguments
p = add_argument(p,'--out',
    help="[Path] Target directory for output files.")
p = add_argument(p,'--verbose',flag=T,
    help="[Flag] Show detailed process.")

### Arguments for perm_test function
p = add_argument(p,'--permu',flag=T,
    help="[Flag] Run fisher test")
p = add_argument(p,'--gwassnp',
    help="[Path] GWAS list file. Columns: <SNPS> <MAPPED_TRAIT> <P.VALUE>")
p = add_argument(p,'--chrstatus',
    help="[Path] BED file including chromosome status for permtest functions.")
p = add_argument(p,'--dbsource', 
    help="[roadmap_bed/encode_bed/gtex_tsv] Type in one of these options for source of your chromosome status.")
p = add_argument(p,'--permn',default=1000, type="numeric",
    help="[Number] Set permutation number. Default=1000")

### Arguments for draw_heatmap function
p = add_argument(p,'--heatmap',flag=T,
    help="[Flag] Draw a heatmap from permu results.")
p = add_argument(p,'--pmdata',
    help="[Path] Z-score table file.")
p = add_argument(p,'--meta',
    help="[Path] Add roadmap meta-info file for heatmap annotation.")
p = add_argument(p,'--range',default="-4,4",
    help="[-4,4] Set coloring Z-score range to display.")
p = add_argument(p,'--annot',
    help="[BLOOD,PANCREAS,THYMUS ...] Choose ANATOMYs of roadmap meta-info to annotate heatmap.")
p = add_argument(p,'--fileext',
    help="[png/sgv] Choose output figure format.")

### Arguments for split_tfbs function
p = add_argument(p,'--splittfbs',flag=T,
    help="[Flag] Split ENCODE TFBS BED file by cell types.")
p = add_argument(p,'--tfbs',
    help="[Path] ENCODE TFBS BED file.")


argv = parse_args(p)

## Load Common Libraries ##

suppressMessages(library(dplyr))
suppressMessages(library(tidyr))

## Functions ##
split_tfbs = function(
    f_tfbs = NULL,
    out    = NULL
) {
    paste0('\n** Run split_tfbs function in enrich.r **\n\n') %>% cat
    ifelse(!dir.exists(out), dir.create(out), "")

    # Read file
    paste0('* ENCODE TFBS table = ') %>% cat
    tfbs = read.delim(f_tfbs, stringsAsFactors=F)
    colnames(tfbs) = c('Chr','Start','End','TF','ids','Cells')
    dim(tfbs) %>% print

    # Get unique cell types
    cell_types = strsplit(tfbs$Cells,',') %>% unlist %>% unique %>% sort
    n = length(cell_types)
    paste0('* ',n,' unique cell types are found.\n\n') %>% cat

    # Expand data by unique cell types
    source('src/pdtime.r')
    o=lapply(c(1:n), function(i) {
        t1=Sys.time()

        # Subset ENCODE data by cell type
        paste0(i,' ',cell_types[i],':\t') %>% cat
        tfbs_sub = subset(tfbs, grepl(cell_types[i],tfbs$Cells))
        m = nrow(tfbs_sub)
        progress = m%/%10
        paste0(m,' regions, filtering [') %>% cat

        # Split file by unique cell types
        tfbs_filt_li = lapply(c(1:m), function(j) {
            if(j%%progress==0) { '.' %>% cat }

            tfbs_sub_j = tfbs_sub[j,1:4]
            tfbs_cell = tfbs_sub[j,]$Cells
            cells = strsplit(tfbs_cell,',')[[1]]
            which_cells = which(cells==cell_types[i])
            if(length(which_cells)>0) {
                tfbs_sub_df = data.frame(tfbs_sub_j, Cell=cells[which_cells]) # debug
            } else tfbs_sub_df = NULL
            return(tfbs_sub_df)
        })
        paste0('] row = ') %>% cat
        tfbs_filt_df = data.table::rbindlist(tfbs_filt_li)
        tf_len = tfbs_filt_df$TF %>% unique %>% length
        paste0(nrow(tfbs_filt_df),', TFs = ',tf_len,'. ') %>% cat

        # Save file as BED format
        f_name = paste0(out,'/',cell_types[i],'.bed')
        write.table(tfbs_filt_df[,1:4],f_name,sep='\t',row.names=F,col.names=F,quote=F)
        paste0('Save: ',f_name,'; ',pdtime(t1,2),'\n') %>% cat
        return(NULL)
    })
}


draw_heatmap = function(
    f_pmdata = NULL,
    f_meta   = NULL,
    range    = NULL,
    out      = 'enrich',
    annot    = NULL,
    fileext  = 'png'
) {
    paste0('\n** Run draw_heatmap function in enrich.r **\n\n') %>% cat
    suppressMessages(library(ComplexHeatmap))
    suppressMessages(library(circlize))
    ifelse(!dir.exists(out), dir.create(out), "")

    # Read file
    paste0('* Permutation result table = ') %>% cat
    pmdata = read.delim(f_pmdata, stringsAsFactors=F)
    dim(pmdata) %>% print

    # Prepare heatmap table
    pmdata_mat = pmdata[,-1]
    rownames(pmdata_mat) = pmdata$Status

    # [Optional] Read meta-info. file
    if(length(f_meta)>0) {
        paste0('* [Optional] Add meta-info. table = ') %>% cat
        meta = read.delim(f_meta, stringsAsFactors=F)
        dim(meta) %>% print

        # Match meta-info. with original colnames
        meta_db = tools::file_path_sans_ext(f_meta %>% basename)
        if(meta_db=='roadmap_meta') {
            pmdata_mat = pmdata_mat[,meta$EID]
            pmdata_col = paste0(meta$EDACC_NAME," (",meta$EID,")")
        } else if(meta_db=='encode_meta') {
            #pmdata_mat = pmdata_mat[,meta$Cell_Name] # Error: undefined columns selected
            pmdata_col = paste0(meta$Cell_Name," (",meta$Symbol,")")
        }
        ha1 = HeatmapAnnotation(Name=anno_text(pmdata_col))
        show_column_names = FALSE

        # [Optional] Add column annotation to heatmap by ANATOMY
        if(length(f_meta)>0 & length(annot)>0) {
            if(meta_db=='roadmap_meta') { anatomy = meta$ANATOMY
            } else if(meta_db=='encode_meta') { anatomy = meta$Tissue }
            annots = strsplit(annot,',')[[1]]
            `%notin%` = Negate(`%in%`) # define %notin% operator
            anatomy[meta$ANATOMY %notin% annots] = 'Other' # A bug to change as NA
            anatomy_uq = anatomy %>% unique %>% sort
            ana_num = length(anatomy_uq)

            if(ana_num<5) { ann_cols = terrain.colors(length(anatomy_uq))
            } else ann_cols = rainbow(length(anatomy_uq))

            ## Set column annotation color
            names(ann_cols) = anatomy_uq
            ha2 = HeatmapAnnotation(
                Anatomy=anatomy,
                col=list(Anatomy=ann_cols),
                gp=gpar(col="black", lwd=.5)
            )
        } else { ha2 = NULL }
    } else {
        paste0('* [Notice] Meta-info table is not found.') %>% cat
        show_column_names = TRUE
        ha1 = NULL
        ha2 = NULL
    }
    
    # Configuration for heatmap
    z_range = strsplit(range,',')[[1]] %>% as.numeric
    pmdata_mat = as.matrix(pmdata_mat)
    my_col = colorRamp2(c(z_range[1],0,z_range[2]),c("#2E86C1","#FEF9E7","#C0392B"))
    if(meta_db=='roadmap_meta') {
        cluster_rows = TRUE
        cluster_columns = TRUE
        wh = c(22,12)
    } else if(meta_db=='encode_meta') {
        cluster_rows = FALSE
        cluster_columns = FALSE
        wh = c(22,25)
    }

    # Draw heatmap
    file_base = basename(f_pmdata)
    file_name = tools::file_path_sans_ext(file_base)
    f_name = paste0(out,'/',file_name,'.png')

    if(fileext=='png') { png(f_name,width=wh[1],height=wh[2],units='in',res=300)
    } else if(fileext=='svg') {}
    p=Heatmap(
        pmdata_mat,
        name = "Enrich.\nz-score",
        col  = my_col,
        rect_gp = gpar(col = "black", lwd = .5),
        column_dend_height = unit(1, "in"),
        row_dend_width = unit(1, "in"),
        column_names_max_height = unit(7,"in"),
        show_column_names = show_column_names,
        bottom_annotation = ha1,
        top_annotation = ha2,
        cluster_rows = cluster_rows,
        cluster_columns = cluster_columns
    )
    print(p)
    dev.off()
    paste0('\nSave as ',f_name,'\n') %>% cat
}


perm_test = function(
    f_gwas_snp = NULL,
    f_status   = NULL,
    db_src     = NULL,
    out        = 'enrich',
    perm_n     = 1000,
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
    cell_types=NULL; zscore_li=list(); pval_li=list(); overlap_li=list()
    source('src/pdtime.r')
    for(i in 1:n) {
        t1=Sys.time()
        file_name = basename(f_status_paths[i])
        if(db_src=='roadmap_bed') {
            cell_type = strsplit(file_name,'_')[[1]][1]
        } else if(db_src=='encode_bed') {
            cell_type = tools::file_path_sans_ext(file_name)
        }
        cell_types = c(cell_types,cell_type)

        ## Read file
        paste0(i,' ',cell_type,': ') %>% cat
        status = read_status_file(f_status_paths[i],db_src)
        paste0(nrow(status),'; ') %>% cat

        ## Run perm_test
        pt_df = perm_test_calc(gwas_snp_bed,status,db_src,perm_n,verbose)
        zscore_li[[i]]  = data.frame(Status=pt_df$Status, Zscore =pt_df$Zscore)
        pval_li[[i]]    = data.frame(Status=pt_df$Status, Pval   =pt_df$Pval)
        overlap_li[[i]] = data.frame(Status=pt_df$Status, Overlap=pt_df$Overlap)
        paste0(pdtime(t1,2),'\n') %>% cat
    }
    pm_merge = function(x,y) { merge(x=x,y=y,by='Status',all=T) }
    zscore_df  = Reduce(pm_merge, zscore_li)
    pval_df    = Reduce(pm_merge, pval_li)
    overlap_df = Reduce(pm_merge, overlap_li)
    colnames(zscore_df)  = c('Status',cell_types)
    colnames(pval_df)    = c('Status',cell_types)
    colnames(overlap_df) = c('Status',cell_types)
    
    # Write result as files
    gwas_base = basename(f_gwas_snp)
    gwas_f_name = tools::file_path_sans_ext(gwas_base)
    
    f_name1 = paste0(out,'/',db_src,'-',gwas_f_name,'-permn_',perm_n,'-zscore.tsv')
    write.table(zscore_df,f_name1,sep='\t',row.names=F,quote=F)
    paste0('\n* Write file: ',f_name1,'\n') %>% cat
    
    f_name2 = paste0(out,'/',db_src,'-',gwas_f_name,'-permn_',perm_n,'-pval.tsv')
    write.table(pval_df,f_name2,sep='\t',row.names=F,quote=F)
    paste0('* Write file: ',f_name2,'\n') %>% cat

    f_name3 = paste0(out,'/',db_src,'-',gwas_f_name,'-permn_',perm_n,'-overlap.tsv')
    write.table(overlap_df,f_name3,sep='\t',row.names=F,quote=F)
    paste0('* Write file: ',f_name3,'\n') %>% cat
}

read_status_file = function(path,db_src) {
    paste0(db_src,' = ') %>% cat
    if(db_src=='roadmap_bed') {
        status_raw = read.delim(path,header=F,skip=1)
        status_raw = status_raw[,1:4]
        colnames(status_raw) = c('Chr','Start','End','Ann')
    } else if(db_src=='encode_bed') {
        status_raw = read.delim(path,header=F)
        status_raw = status_raw[,1:4]
        colnames(status_raw) = c('Chr','Start','End','Ann')
    } else if(db_src=='gtex_tsv') {
        #
    } else {
        paste0('\n\n[Error] Unknown option is flagged. Please choose one of [roadmap_bed/encode].\n') %>% cat
        stop()
    }
    return(status_raw)
}

perm_test_calc = function(
    gwas_snp_bed = NULL,
    status   = NULL,
    db_src   = NULL,
    perm_n   = 5000,
    verbose  = FALSE
) {
    #paste0('\n** Run perm_test function in enrich.r **\n\n') %>% cat
    #suppressMessages(library(regioneR))

    # Calculate permTest to get z-scores and p-values by chromosome status
    paste0('permTest - ') %>% cat
    status_ann = status$Ann %>% unique %>% sort
    n = length(status_ann)
    paste0(n,' annots [') %>% cat

    # Set background status as universe
    status_all_bed = toGRanges(status[,1:4],format="BED")
    # if(db_src=='roadmap_bed') {
    #     status_all_bed = toGRanges(status[,1:4],format="BED")
    # } else if (db_src=='encode_bed') {
    #     encode_tfbs_v3 = readRDS('db/wgEncodeRegTfbsClusteredV3.bed.rds')
    #     status_all_bed = toGRanges(encode_tfbs_v3[,1:4],format="BED")
    # }
    
    progress = n%/%5
    pt_li = lapply(c(1:n),function(i) {
        if(n>=5 & i%%progress==0) { '.' %>% cat
        } else '.' %>% cat
        # Convert data.frame to GRanges format
        status_sub = subset(status,Ann==status_ann[i])
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
        pt_df = data.frame(
            Status  = status_ann[i],
            Zscore  = pt$numOverlaps$zscore,
            Pval    = pt$numOverlaps$pval,
            Overlap = pt$numOverlaps$observed
        )
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
    draw_heatmap(
        f_pmdata = argv$pmdata,
        f_meta   = argv$meta,
        range    = argv$range,
        out      = argv$out,
        annot    = argv$annot,
        fileext  = argv$fileext
    )
} else if(argv$splittfbs) {
    split_tfbs(
        f_tfbs = argv$tfbs,
        out    = argv$out
    )
}
paste0('\n',pdtime(t0,1),'\n') %>% cat
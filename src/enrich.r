suppressMessages(library(argparser))

## Parsing Arguments ##
p = arg_parser("Function for enrichment analysis (Fisher/Permutation)")

### Examples
p = add_argument(p, '--example', help="See examples to run this tool.", flag=T)
example_msg = '
Rscript src/enrich.r --splittfbs \
    --tfbs db/wgEncodeRegTfbsClusteredWithCellsV3.bed \
    --out db/tfbs_cell

Rscript src/enrich.r --splitgtex \
    --gtex db/gtex_signif_5e-8.tsv.rds \
    --out db/gtex_tsv

Rscript src/enrich.r --permu \
    --gwassnp data/seedSNP_1817_bm.bed \
    --chrstatus db/roadmap_bed \
    --dbsource roadmap_bed \
    --permn 1000 \
    --out enrich

Rscript src/enrich.r --permu \
    --gwassnp data/gwas_5e-08_129_hg19.bed \
    --chrstatus db/encode_bed \
    --dbsource encode_bed \
    --permn 100 \
    --out enrich

Rscript src/enrich.r --gtexperm \
    --gtex_base data/gtex_5e-08_745.tsv \
    --snp_filt db/gtex_tsv \
    --gtex_median_tpm db/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_median_tpm.gct
    --permn 100 \
    --out enrich

Rscript src/enrich.r --heatmap \
    --pmdata enrich/roadmap_bed-snp_484_roadmap_dist-permn_100-zscore.tsv \
    --meta db/roadmap_meta.tsv \
    --out enrich \
    --range -3,3 \
    --annot BLOOD,PANCREAS,THYMUS \
    --fileext png
'

### Shared Arguments
p = add_argument(p,'--gwas_snp',
    help="[Path] GWAS list BED file. Columns: <Chr> <Start> <End> <Rsid>")
p = add_argument(p,'--out',
    help="[Path] Target directory for output files.")
p = add_argument(p,'--verbose',flag=T,
    help="[Flag] Show detailed process.")
p = add_argument(p,'--perm_n',default=1000, type="numeric",
    help="[Number] Set permutation number. Default=1000")

### Arguments for perm_test function
p = add_argument(p,'--permu',flag=T,
    help="[Flag] Run permutation test for Roadmap/ENCODE data.")
p = add_argument(p,'--chr_status',
    help="[Path] BED file including chromosome status for permtest functions.")
p = add_argument(p,'--db_source', 
    help="[roadmap_bed/encode_bed] Type in one of these options for source of your chromosome status.")

### Arguments for draw_heatmap function
p = add_argument(p,'--heatmap',flag=T,
    help="[Flag] Draw a heatmap from permu results.")
p = add_argument(p,'--pm_data',
    help="[Path] Z-score table file.")
p = add_argument(p,'--meta',
    help="[Path] Add roadmap meta-info file for heatmap annotation.")
p = add_argument(p,'--range',default="-4,4",
    help="[-4,4] Set coloring Z-score range to display.")
p = add_argument(p,'--annot',
    help="[BLOOD,PANCREAS,THYMUS ...] Choose ANATOMYs of roadmap meta-info to annotate heatmap.")
p = add_argument(p,'--file_ext',
    help="[png/sgv] Choose output figure format.")

### Arguments for split_tfbs function
p = add_argument(p,'--split_tfbs',flag=T,
    help="[Flag] Split ENCODE TFBS BED file by cell types.")
p = add_argument(p,'--tfbs',
    help="[Path] ENCODE TFBS BED file.")

### Arguments for split_gtex function
p = add_argument(p,'--split_gtex',flag=T,
    help="[Flag] Split GTEx eQTL file by tissue types.")
p = add_argument(p,'--gtex',
    help="[Path] GTEx eQTL RDS file.")

### Arguments for gtex_perm_test function
p = add_argument(p,'--gtex_perm',flag=T,
    help="[Flag] Run permutation test for GTEx eQTL data.")
p = add_argument(p,'--gtex_base',
    help="[Path] GTEx overlapped TSV file.")
p = add_argument(p,'--gtex_median_tpm',
    help="[Path] GTEx gene median tpm GCT file.")


argv = parse_args(p)

## Load Common Libraries ##

suppressMessages(library(dplyr))
suppressMessages(library(tidyr))

## Functions ##
split_gtex = function(
    f_gtex = NULL,
    out    = NULL
) {
    paste0('\n** Run split_gtex function in enrich.r **\n\n') %>% cat
    ifelse(!dir.exists(out), dir.create(out), "")

    # Read file
    paste0('* GTEx table = ') %>% cat
    gtex = readRDS(f_gtex)
    dim(gtex) %>% print

    # Get unique tissue types
    tissue_types = gtex$tissue %>% unique %>% sort
    n = length(tissue_types)
    paste0('* ',n,' unique tissue types are found.\n\n') %>% cat

    # Split file by unique tissue types
    source('src/pdtime.r')
    o=lapply(c(1:n), function(i) {
        t1 = Sys.time()
        paste0(i,' ',tissue_types[i],':\t') %>% cat
        gtex_sub = subset(gtex,tissue==tissue_types[i])
        m = nrow(gtex_sub)
        paste0(m,' pairs, SNPs = ') %>% cat
        snps_len = gtex_sub$variant_id %>% unique %>% length
        gene_len = gtex_sub$gene_id %>% unique %>% length
        paste0(snps_len,', genes = ',gene_len,'. ') %>% cat

        # Save file as BED format
        f_name = paste0(out,'/',tissue_types[i],'.tsv')
        write.table(gtex_sub[,c(1:5,9)],f_name,sep='\t',row.names=F,quote=F)
        paste0('Save: ',f_name,'; ',pdtime(t1,2),'\n') %>% cat
        return(NULL)
    })
}


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


gtex_perm_test = function(
    f_gwas_snp  = NULL,
    f_gtex_base = NULL,
    f_gtex_median_tpm = NULL,
    out         = 'enrich',
    perm_n      = 1000
) {
    paste0('\n** Run gtex_perm_test function in enrich.r **\n\n') %>% cat
    suppressMessages(library(fgsea))
    suppressMessages(library(data.table))
    suppressMessages(library(ggplot2))
    suppressMessages(library(tidyr))
    ifelse(!dir.exists(out), dir.create(out), "")

    # Read files
    paste0('* Input SNPs = ') %>% cat
    snps = read.delim(f_gwas_snp, header=F, stringsAsFactors=F)
    colnames(snps) = c('Chr','Start','End','Rsid')
    dim(snps) %>% print

    paste0('* GTEx eQTL pairs = ') %>% cat
    gtex_pairs = readRDS(f_gtex_base)
    colnames(gtex_pairs) = c('Variant_id','Gene_id','Pval','Slope','Slope_se','Tissue','Chr','Pos','Rsid')
    gtex_pairs$Ensgid = sapply(gtex_pairs$Gene_id,function(x) strsplit(x,'\\.')[[1]][1])
    dim(gtex_pairs) %>% print

    paste0('* GTEx gene median tpm GCT = ') %>% cat
    gct = readRDS(f_gtex_median_tpm)
    nrow(gct) %>% cat
    #gct_negate = data.frame(gct[,1:2],-gct[,3:ncol(gct)])
    #' -> negate ' %>% cat
    gct_gather = gct %>% gather("Tissue","TPM",-Ensgid,-Description)
    ' -> transform ' %>% cat
    eqtl_genes = gtex_pairs$Gene_id %>% unique
    gct_filt_gene = subset(gct_gather, Ensgid %in% eqtl_genes)$Ensgid %>% unique
    paste0('-> Filtering= ',length(gct_filt_gene),'\n') %>% cat

    # Filter SNP eQTL genes and their tissues
    paste0('* Filtering eQTLs by SNP: ') %>% cat
    gtex_pairs_sub = subset(gtex_pairs,Rsid %in% snps$Rsid)
    gene_eqtl = gtex_pairs_sub$Gene_id %>% unique %>% sort
    paste0('genes = ',length(gene_eqtl)) %>% cat
    eqtl_tissue = gtex_pairs_sub$Tissue %>% unique %>% sort
    n = length(eqtl_tissue)
    paste0(', tissues = ',n) %>% cat
    
    # Preparing gene sets by tissue eQTL genes
    ' [' %>% cat
    geneset_li = lapply(c(1:n),function(i) {
        if(i%%5==0) '.' %>% cat
        subset(gtex_pairs_sub,Tissue==eqtl_tissue[i])$Ensgid %>% unique
    })
    names(geneset_li) = eqtl_tissue
    '] done.\n\n' %>% cat

    # Run by each tissue
    fgseaRes_li = list()
    for(i in 18:n) {
        paste0(i,' ',eqtl_tissue[i],': ') %>% cat
        gtex_pairs$Tissue %>% table %>% print #<- debugging..

        # Prepare data
        eqtlgene_tissue_all = subset(gtex_pairs,Tissue==eqtl_tissue[i])$Ensgid %>% unique
        gene_n = length(geneset_li[[i]])
        paste0('genes = ',gene_n,' / ',length(eqtlgene_tissue_all),', ') %>% cat
        
        ranked_ts = subset(gct_gather, Tissue==eqtl_tissue[i])
        ranked_ts_gene = subset(ranked_ts, Ensgid %in% eqtlgene_tissue_all)
        ranked = ranked_ts_gene$TPM
        names(ranked) = ranked_ts_gene$Ensgid
        ranked = ranked[order(ranked_ts_gene$TPM)] # sort by rank
        paste0('exp. = ',length(ranked),'\n') %>% cat

        # Run fgsea (fgsea)
        inter = intersect(geneset_li[[i]],names(ranked))
        result = fgsea(
            pathways = geneset_li[i],
            stats    = ranked,
            nperm    = perm_n
        )
        fgseaRes_li[[i]] = data.frame(
            result[,1:7],
            queried_gene=gene_n,
            total=length(ranked)
        )
    }
    fgseaRes = data.table::rbindlist(fgseaRes_li) %>% as.data.frame
    # Save result as a file
    f_base = basename(f_gwas_snp)
    file_base = tools::file_path_sans_ext(f_base)
    f_name1 = paste0(out,'/gtex-',file_base,'-permn_',perm_n,'.tsv')
    write.table(fgseaRes[order(fgseaRes$pval), ],f_name1,sep='\t',row.names=F)
    paste0('\nWrite file: ',f_name1,'\n') %>% cat


    # Prepare tissues for plot
    topTissueUp_df = fgseaRes[fgseaRes$ES>0,]
    if(nrow(topTissueUp_df)>0) { 
        topTissueUp = topTissueUp_df[order(topTissueUp_df$pval),]$pathway
    } else topTissueUp = NULL
    topTissueDn_df = fgseaRes[fgseaRes$ES<0,]
    if(nrow(topTissueDn_df)>0) { 
        topTissueDn = topTissueDn_df[order(topTissueDn_df$pval),]$pathway
    } else topTissueDn = NULL
    topTissues  = c(topTissueUp,rev(topTissueDn))

    # Draw plot
    f_name2 = paste0(out,'/gtex-',file_base,'-permn_',perm_n,'.png')
    png(f_name2,width=10,height=20,units='in',res=200)
    plotGseaTable(geneset_li[topTissues],ranked,fgseaRes, gseaParam=0.5)
    dev.off()
    paste0('Draw figure: ',f_name2,'\n') %>% cat

    # paste0('\n* SNP overlapped eQTL gene-tissue pairs = ') %>% cat
    
    # gsa_df = GSAsummaryTable(gsaRes)
    # zscore_df  = data.frame(Name=c('Stat (dist.dir)'),
    #     gsa_df %>% select('Stat (dist.dir)') %>% t)
    # pval_df    = data.frame(Name=c('p (dist.dir.up)','p (dist.dir.dn)'),
    #     gsa_df %>% select('p (dist.dir.up)','p (dist.dir.dn)') %>% t)
    # overlap_df = data.frame(Name=c('Genes (tot)','Genes (up)','Genes (down)'),
    #     gsa_df %>% select('Genes (tot)','Genes (up)','Genes (down)') %>% t)
    # colnames(zscore_df)  = c('Name',gsa_df$Name)
    # colnames(pval_df)    = c('Name',gsa_df$Name)
    # colnames(overlap_df) = c('Name',gsa_df$Name)
    
    # # Write TSV file
    # file_name = tools::file_path_sans_ext(f_gtex_base %>% basename)

    # f_name1 = paste0(out,'/',file_name,'-permn_',perm_n,'-zscore.tsv')
    # write.table(zscore_df,f_name1,sep='\t',row.names=F,quote=F)
    # paste0('\n* Write file: ',f_name1,'\n') %>% cat
    
    # f_name2 = paste0(out,'/',file_name,'-permn_',perm_n,'-pval.tsv')
    # write.table(pval_df,f_name2,sep='\t',row.names=F,quote=F)
    # paste0('* Write file: ',f_name2,'\n') %>% cat

    # f_name3 = paste0(out,'/',file_name,'-permn_',perm_n,'-overlap.tsv')
    # write.table(overlap_df,f_name3,sep='\t',row.names=F,quote=F)
    # paste0('* Write file: ',f_name3,'\n') %>% cat
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
    if(n==1) { f_status_paths = f_status; n=1 }
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
    pm_merge = function(x,y) { merge(x=x, y=y, by='Status', all=T) }
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
    perm_n   = 1000,
    verbose  = FALSE
) {
    # Calculate permTest to get z-scores and p-values by chromosome status
    paste0('permTest - ') %>% cat
    status_ann = status$Ann %>% unique %>% sort
    n = length(status_ann)
    paste0(n,' annots [') %>% cat

    # Set background status as universe
    status_all_bed = toGRanges(status[,1:4],format="BED")
    
    check_n = 5; progress = n%/%check_n
    pt_li = lapply(c(1:n),function(i) {
        if(n>=check_n & i%%progress==0) { '.' %>% cat } 
        else if(n<check_n) { '.' %>% cat }
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

## SQLite code sniffet ##
library(RSQLite)

path1 = 'db/gtex_analysis_v8_rnaseq_gene_median_tpm_ensgid.gct.rds'
path2 = 'db/gtex_signif_5e-8.tsv.rds'
f_db  = 'db/gtex_signif_5e-8.db'
db_name = "gtex"

gct   = readRDS(path1); dim(gct)
pair  = readRDS(path2); dim(pair)
colnames(pair) = c('Variant_id','Gene_id','Pval','Slope','Slope_se','Tissue','Chr','Pos','Rsid')
pair$Ensgid = sapply(pair$Gene_id,function(x) strsplit(x,'\\.')[[1]][1])

conn = dbConnect(RSQLite::SQLite(),f_db)
dbWriteTable(conn, db_name, pair)

my_query = paste0('SELECT * FROM ',db_name,' LIMIT 10')
dbGetQuery(conn, my_query) %>% print
## Sniffet End ##

## Run Function ##
source('src/pdtime.r')
t0 = Sys.time()

if(argv$example) {
    cat(example_msg)
} else if(argv$permu) {
    perm_test(
        f_gwas_snp = argv$gwas_snp,
        f_status   = argv$chr_status,
        db_src     = argv$db_source,
        perm_n     = argv$perm_n,
        out        = argv$out,
        verbose    = argv$verbose
    )
} else if(argv$heatmap) {
    draw_heatmap(
        f_pmdata = argv$pm_data,
        f_meta   = argv$meta,
        range    = argv$range,
        out      = argv$out,
        annot    = argv$annot,
        fileext  = argv$file_ext
    )
} else if(argv$split_tfbs) {
    split_tfbs(
        f_tfbs = argv$tfbs,
        out    = argv$out
    )
} else if(argv$split_gtex) {
    split_gtex(
        f_gtex = argv$gtex,
        out    = argv$out
    )
} else if(argv$gtex_perm) {
    gtex_perm_test(
        f_gwas_snp  = argv$gwas_snp,
        f_gtex_base = argv$gtex_base,
        f_gtex_median_tpm = argv$gtex_median_tpm,
        out         = argv$out,
        perm_n      = argv$perm_n
    )
}

paste0('\n',pdtime(t0,1),'\n') %>% cat
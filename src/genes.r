## Change log ##
# Written by Seungsoo Kim, PhD
# 

## Parsing Arguments ##
suppressMessages(library(argparser))
p = arg_parser("Function for enrichment analysis (Fisher/Permutation)")

### Command examples
p = add_argument(p, '--example',flag=T,
    help="See command examples by functions.")
example_msg = '
Rscript src/genes.r --disgenet \
    --f_gene data/disgenet/TFs_157.tsv \
    --out data/disgenet

Rscript src/genes.r --cprofiler \
    --
'

### Shared Arguments
p = add_argument(p,'--disgenet',flag=T,
    help="[Function] Retrieve gene associated diseases from DisGeNET.\nSee details at https://www.disgenet.org/disgenet2r\nor https://www.disgenet.org/static/disgenet2r/disgenet2r.html")
p = add_argument(p,'--f_gene',
    help="[Path] Gene list TSV file. Column: <Gene> <...>")
p = add_argument(p,'--out',
    help="[Path] ")

p = add_argument(p,'--cprofiler',flag=T,
    help="[Function] Calculate gene set enrichment to ontologies.\nSee details at https://www.bioconductor.org/packages/release/bioc/html/clusterProfiler.html\nor https://www.bioconductor.org/packages/release/bioc/vignettes/clusterProfiler/inst/doc/clusterProfiler.html")

argv = parse_args(p)

## Load Common Libraries ##
suppressMessages(library(dplyr))

## Functions ##
disgenet = function(
    f_gene = NULL,
    out    = NULL
) {
    # Load library
    suppressMessages(library(disgenet2r))
    ifelse(!dir.exists(out), dir.create(out), "")

    # Configurations
    gda_plot_wh = c(8.5,8) # inch
    enrich_plot_wh = c(10,10) # inch

    f_base = tools::file_path_sans_ext(f_gene %>% basename)
    f_name_gda = paste0(out,'/',f_base,"_gda.tsv")
    f_name_gda_plot = paste0(tools::file_path_sans_ext(f_name_gda),'.png')
    f_name_enrich = paste0(out,'/',f_base,"_enrich.tsv")
    f_name_enrich_plot = paste0(tools::file_path_sans_ext(f_name_enrich),'.png')

    # Read file
    paste0('\n* Genes = ') %>% cat
    gene = read.delim(f_gene,stringsAsFactors=F)$Gene %>% unique
    length(gene) %>% print

    # Retrieving gene-disease associations (GDAs)
    paste0('* Retrieve from DisGeNet = ') %>% cat
    gdaRes = gene2disease(gene=gene,vocabulary="HGNC",database="CURATED",verbose=T)
    gdaRes_df = extract(gdaRes)
    dim(gdaRes_df) %>% print

    write.table(gdaRes_df,f_name_gda,row.names=F,quote=F,sep='\t')
    paste0('  Write file: ',f_name_gda,'\n') %>% cat

    # Draw result as plot
    png(f_name_gda_plot,width=gda_plot_wh[1],height=gda_plot_wh[2],units='in',res=150)
    plot(gdaRes,class="Heatmap",prop=10) # DiseaseClass, Heatmap, Network
    dev.off()
    paste0('  Draw plot: ',f_name_gda_plot,'\n') %>% cat

    # Performing a disease enrichment
    paste0('* Enrichment = ') %>% cat
    enrichRes = disease_enrichment(entities=gene,vocabulary="HGNC",database="CURATED",universe="HUMAN_CODING")
    enrichRes_df = enrichRes@qresult[,c("Description","FDR","Ratio","BgRatio")]
    dim(enrichRes_df) %>% print

    write.table(enrichRes_df,f_name_enrich,row.names=F,quote=F,sep='\t')
    paste0('  Write file: ',f_name_enrich,'\n') %>% cat

    # Draw enrich result as plot
    png(f_name_enrich_plot,width=enrich_plot_wh[1],height=enrich_plot_wh[2],units='in',res=150)
    plot(enrichRes,class="Enrichment",count=3,cutoff=0.05)
    dev.off()
    paste0('  Draw plot: ',f_name_enrich_plot,'\n') %>% cat
}

cprofiler = function(
    f_gene = NULL,
    out    = NULL
) {
    # Load library
    suppressMessages(library(clusterProfiler))
    ifelse(!dir.exists(out), dir.create(out), "")

    # Configuration
    f_base = tools::file_path_sans_ext(f_gene %>% basename)
    f_name = paste0(out,'/',f_base,'.tsv')

    # Read file
    paste0('\n* Genes = ') %>% cat
    gene = read.delim(f_gene,stringsAsFactors=F)$Gene %>% unique
    length(gene) %>% print

    # Universal enrichment analysis
    paste0('1. WikiPathways analysis: ') %>% cat
    
}

## Functions End ##

## Run Function ##
source('src/pdtime.r')
t0 = Sys.time()

if(argv$example) {
    cat(example_msg)
} else if(argv$disgenet) {
    disgenet(
        f_gene = argv$f_gene,
        out    = argv$out
    )
}

paste0('\n',pdtime(t0,1),'\n') %>% cat
## SQLite code snippet ##
library(dplyr)
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
## snippet End ##
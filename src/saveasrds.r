# Read file and save as rds file
library(data.table)
saveasrds = function(f_path) {
    df = fread(f_path,sep='\t',header=T,stringsAsFactors=F)
    f_name = paste0(f_path,'.rds')
    save(df,file=f_name)
}
@echo off
: Filtering GWAS Catalog data
:: Rscript T1D_gwas.r [GWAS_file_path] [p-value_criteria]
Rscript T1D_gwas.r db/GWAS_EFO0001359.tsv 5e-08

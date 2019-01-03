@echo off

: 1. Filtering GWAS Catalog data
:: Rscript T1D_gwas.r [GWAS_file_path] [p-value_criteria]
<<<<<<< HEAD
Rscript T1D_gwas.r db/GWAS_EFO0001359.tsv 5e-08
=======
Rscript T1D_gwas.r db/GWAS_EFO0001359.tsv

: 2. Download LDlink data
:: python T1D_ldlink.py [SNP_file_path.txt]
python T1D_ldlink.py SNP_5e-08_129.tsv

: 3. Filtering LDlink data
:: Rscript T1D_ldlink.r [SNP_file_path.txt] [LDlink_data_folder_path]
Rscript T1D_ldlink.r SNP_5e-08_129.tsv db/ldlink
>>>>>>> parent of 0c366a1... write: db/roadmap - core functions

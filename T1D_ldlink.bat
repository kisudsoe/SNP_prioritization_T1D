@echo off
<<<<<<< HEAD

: pre2-1. Download LDlink data
:: python T1D_ldlink.py [SNP_file_path.txt]
python T1D_ldlink_dn.py data/gwas_5e-08_129.tsv

: pre2-2. Filtering LDlink data
:: Rscript T1D_ldlink.r [SNP_file_path.txt] [LDlink_download_target_dir]
Rscript T1D_ldlink_filt.r data/gwas_5e-08_129.tsv db/ldlink
=======
:: Rscript T1D_ldlink.r [SNP_file_path.txt] [LDlink_data_folder_path]

Rscript T1D_ldlink.r SNP_5e-08_129.tsv db/ldlink
>>>>>>> parent of 0c366a1... write: db/roadmap - core functions

@echo off
: 1. Download LDlink data
:: python T1D_ldlink.py [SNP_file_path]
python ldlink_dn.py data/gwas_5e-08_129.tsv

: 2. Filtering LDlink data
:: Rscript T1D_ldlink.r [SNP_file_path] [LDlink_download_target_dir]
Rscript ldlink_filt.r data/gwas_5e-08_129.tsv db/ldlink

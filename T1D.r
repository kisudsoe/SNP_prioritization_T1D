#!/usr/bin/env Rscript
# This file is for prioritizing T1D SNPs

## Command Arg Parameters ##
# T1D.bat: Rscript T1D.r [data/rodamap_inter.bed] []
args = commandArgs(trailingOnly=T)
hmsg = 'Rscript T1D_ldlink_filt.r [SNP_file_path] [LDlink_download_target_dir]
  - Arguments [SNP_file_path] and [LDlink_download_target_dir] are mendatory.'
if(length(args)<2|length(args)>2) stop(hmsg)
snp.path = unlist(args[1])
ld.path = unlist(args[2])

# System parameter
source('src/pdtime.r')
t0 = Sys.time()

#########################
## Function start here ##
#########################


##################
## Function end ##
##################
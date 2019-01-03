@echo off
: 1. Download RoadMap data
:: No argument in roadmap_dn.r
:: Target download foler is `db/roadmap`. Total file size is 1.62 GB.
Rscript roadmap_dn.r

: 2. Compile RoadMap data
:: No argument in roadmap_filt.r
:: Result BED file is generated at `db/roadmap_enh.bed`. The BED file size is 139 MB.
Rscript roadmap_filt.r

: 3. Prioritizing SNPs
:: Rscript roadmap.r [data/roadmap_dist.tsv]
Rscript roadmap.r data/roadmap_dist.tsv
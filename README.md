# PriSNP: Prioritization of SNPs

This is a project for prioritization of SNPs associated with Type 1 diabetes.



## 1. SNP data download and filter

### 1-1. GWAS Catalog data download and filter

To download **GWAS Catalog data**, you can [search certain disease](https://www.ebi.ac.uk/gwas/). In this study, we downloaded [SNP-sets for type 1 diabetes](https://www.ebi.ac.uk/gwas/efotraits/EFO_0001359). Then you can run R code file for filtering the GWAS Catalog data as below command line:

```cmd
Rscript T1D_gwas.r [GWAS_file_path] [p-value_criteria]
```



### 1-2. LDlink data download and compile

To download **LDlink data**, you can run `T1D_ldlink.py` as below `CMD` command line:

- To run the code, you need list of SNP RS IDs of dbSNP database as txt file

```CMD
python ldlink.py [SNP_file_path.txt]
```

To Filter the LDlink data, you can run `T1D_ldlink.r` as below `CMD` command line:

```CMD
Rscript ldlink.r [SNP_file_path.txt] [LDlink_data_folder_path]
```

To use bedtools later, you have to prepare SNP list as [bed format](https://genome.ucsc.edu/FAQ/FAQformat.html). If you have simple dbSNP rsid list, you can run `src/snp_biomart.r` for generate bed file. But you should check `NA` values and fill it manually.

```CMD
Rscript src/snp_bed.r
```



## 2. RoadMap data download and filter

The [RoadMAP project](https://egg2.wustl.edu/roadmap/web_portal/imputed.html) provides epigenome annotations such as [12-mark/127-reference epigenome/25-state Imputation Based Chromatin State Model](https://egg2.wustl.edu/roadmap/data/byFileType/chromhmmSegmentations/ChmmModels/imputed12marks/jointModel/final/). We downloaded the 127 files by their [cell types](https://github.com/mdozmorov/genomerunner_web/wiki/Roadmap-cell-types) (e.g., `E001_25_imputed12marks_hg38lift_dense.bed.gz` and etc) using R code (`T1D_roadmap.r`). And then we filtered the data by [annotation code](https://egg2.wustl.edu/roadmap/web_portal/imputed.html) (see db/[roadmap] ) including 13_EnhA1, 14_EnhA2, 15_EnhAF, 16_EnhW1, 17_EnhW2, 18_EnhAc.

To download RoadMap data, you need to install `AnnotationHub` and `rtracklayer` in `BiocManaer` as below R code:

```R
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("AnnotationHub", version = "3.8")
BiocManager::install("rtracklayer", version = "3.8")
```

To download **RoadMap data**, you can run `roadmap_dn.r` as below `CMD` command line:

- The RoadMap BED files will download at `db/roadmap` folder
- This process takes about ~10 min that depends on your download speed.

```CMD
Rscript roadmap_dn.r
```

To filter the RoadMap data by **Enhancers**, you can run `roadmap_filt.r` as below `CMD` command line:

- The result file would be saved as `data/roadmap_enh.bed`
- This process takes >12 min that depends on your hard drive read speed.

```CMD
Rscript roadmap_filt.r
```

To reduce file size and faster process, merge RoadMap enhancer information using a `BASH` tool `bedtools`. Here is the `BASH` pipeline for `bedtools sort` and `bedtools merge`. Then, to identify T1D SNPs occupied in RoadMap enhancers, you can use `BASH` tool `bedtools intersect` as below code:

- Compressed file size of `roadmap_enh.bed.gz` is >139 MB.
- Compressed file size of `roadmap_enh_merer.bed.gz` is about 3.7 MB.

```BASH
bedtools sort -i db/roadmap_enh.bed | bedtools merge -i stdin -c 1 -o count > db/roadmap_enh_merge.bed
bedtools closest -d -a data/seedSNP_1817.bed -b db/roadmap_enh_merge.bed > data/roadmap_dist.bed
```



```BASH

```





## 3. ENCODE ChIP-seq data download and filter


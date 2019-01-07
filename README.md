# Disease SNP prioritization

This is an protocol for prioritization of SNPs associated certain phenotype/disease. Here is a study for prioritization of SNPs associated with Type 1 diabetes. You can follow the below analysis steps.



## 1. Seed SNPs preparation for type 1 diabetes (T1D)

### gwas.r

To download **GWAS Catalog data**, you can [search certain disease](https://www.ebi.ac.uk/gwas/). In this study, we downloaded [SNP-sets for type 1 diabetes](https://www.ebi.ac.uk/gwas/efotraits/EFO_0001359). Then you can run R code file for filtering the GWAS Catalog data as below command line:

- Instead of `[ ]`, you have to put the arguments `file path` or `value` by the options.

```cmd
Rscript gwas.r [GWAS_file_path] [p-value_criteria]
```

### ldlink.py/ ldlink.r

To download **LDlink data**, you can run `T1D_ldlink.py` as below `CMD` command line:

- To run the code, you need list of SNP RS IDs of dbSNP database as txt file

```CMD
python ldlink.py [SNP_file_path.txt]
```

To Filter the LDlink data, you can run `T1D_ldlink.r` as below `CMD` command line:

```CMD
Rscript ldlink.r [SNP_file_path.txt] [LDlink_data_folder_path]
```

### Q1. How to make private SNP list BED file?

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

### roadmap_dn.r

To download **RoadMap data**, you can run `roadmap_dn.r` as below `CMD` command line:

- The RoadMap BED files will download at `db/roadmap` folder
- This process takes about ~10 min that depends on your download speed.

```CMD
Rscript roadmap_dn.r
```

### roadmap_filt.r

To filter the RoadMap data by **Enhancers**, you can run `roadmap_filt.r` as below `CMD` command line:

- The result file would be saved as `data/roadmap_enh.bed`
- This process takes >12 min that depends on your hard drive read speed.

```CMD
Rscript roadmap_filt.r
```

### bedtools merge/ bedtools closest

To avoid multiple count of enhancers as well as to reduce file size and to achieve faster process, merge RoadMap enhancer information using a `BASH` tool `bedtools`. Here is the `BASH` pipeline for `bedtools sort` and `bedtools merge`. Then, to identify T1D SNPs occupied in RoadMap enhancers, you can use `BASH` tool `bedtools intersect` as below code:

- Compressed file size of `roadmap_enh.bed.gz` is >139 MB.
- Compressed file size of `roadmap_enh_merer.bed.gz` is about 3.7 MB.

```BASH
bedtools sort -i db/roadmap_enh.bed | bedtools merge -i stdin -c 1 -o count > db/roadmap_enh_merge.bed
bedtools closest -d -a data/seedSNP_1817.bed -b db/roadmap_enh_merge.bed > data/roadmap_dist.tsv
```

### src/bedtools_closest.r

To prioritize RoadMap enhancer occupied SNPs, you can run `roadmap.r` as below `CMD` command line:

- `data/roadmap_dist_df.tsv` file is obtained that is for enhancer annotated file .
- `data/snp_484_roadmap_dist.bed` file is obtained that is for `BED` format file for USCS browser.

```BASH
::Usage: Rscript src/bedtools_closest.r [bedtools_closest_result_file_path]
Rscript src/bedtools_closest.r data/roadmap_dist.tsv
```

### Q2. How about just use not merged roadmap_enh.bed file?

Instead of merge file, when you use original `db/roadmap_enh.bed` file, you can find a lot of duplicated enhancers regions.

```BASH
bedtools sort -i db/roadmap_enh.bed | bedtools closest -d -a data/seedSNP_1817.bed -b stdin > data/roadmap_dist2.tsv
```



## 3. ENCODE ChIP-seq data download and filter

The **ENCODE ChIP-seq** for regulatory transcription factor binding site (Reg-TFBS) cluster data can downloaded <u>wgEncodeRegTfbsClusteredV3</u> data from [UCSC FTP](http://hgdownload.cse.ucsc.edu/goldenpath/hg19/encodeDCC/wgEncodeRegTfbsClustered/) (68 MB) or [bioconductor `data("wgEncodeTfbsV3")`](https://www.bioconductor.org/packages/devel/bioc/vignettes/ChIPpeakAnno/inst/doc/ChIPpeakAnno.html). Here, we assume having downloaded UCSC FTP file `wgEncodeRegTfbsClusteredV3.bed.gz` (81 MB).

### bedtools merge/ bedtools closest

To identify TFBS occupied SNPs, you can use `bedtools merge` and `bedtools closest` as following code:

- Merging the ENCODE TFBS data give you benefits such as avoiding multiple count of enhancers as well as reducing file size and achieving faster process

```BASH
bedtools merge -i db/wgEncodeRegTfbsClusteredV3.bed.gz -c 1 -o count > db/encode_tfbs_merge.bed
bedtools closest -d -a data/seedSNP_1817.bed -b db/encode_tfbs_merge.bed > data/encode_dist.tsv
```

### src/bedtools_closest.r

To prioritize RoadMap enhancer occupied SNPs, you can run `roadmap.r` as below `CMD` command line:

- `data/roadmap_dist_df.tsv` file is obtained that is for enhancer annotated file .
- `data/snp_enh_484.bed` file is obtained that is for `BED` format file for USCS browser.

```CMD
::Usage: Rscript src/bedtools_closest.r [bedtools_closest_result_file_path]
Rscript src/bedtools_closest.r data/encode_dist.tsv
```



## 4. Regulome DB data download and filter

The [Regulome DB](http://www.regulomedb.org/downloads) provides [category scores for SNPs by evidences](http://www.regulomedb.org/help) (see `Regulome score.txt`), including eQTL, TF binding, matched TF motif, matched DNase Footprint, and DNase peak. In this study, we stringently filtered and used high-score (`≥ 2b`) SNPs for our study.

- `RegulomeDB.dbSNP132.Category1.txt.gz` (2 MB)
- `RegulomeDB.dbSNP132.Category2.txt.gz` (39.3 MB)
- Or you can download total dataset: `RegulomDB.` (2.8 GB)

```CMD
Rscript regulome.r data/seedSNP_1817_ldlink.tsv db/RegulomeDB.dbSNP132.Category1.txt.gz db/RegulomeDB.dbSNP132.Category2.txt.gz
```

The result files are save at `data/` folder:

- `data/regulome_#.tsv`
- `data/snp_#_regulome2b.bed`



## 5. Venn analysis of functional motif SNPs

Summary for SNPs with RoadMap annotation, ENCODE ChIP-seq, and RegulomeDB.

```CMD

```


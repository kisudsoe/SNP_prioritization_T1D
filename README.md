# Disease SNP prioritization

This is an protocol for prioritization of SNPs associated certain phenotype/disease. Here is a study for prioritization of SNPs associated with Type 1 diabetes. You can follow the below analysis steps.



## 1. Seed SNPs preparation for type 1 diabetes (T1D)

### > gwas.r

To download **GWAS Catalog data** (MacArthur et al, 2017, Nucleic Acids Research, pmid 27899670), you can [search certain disease](https://www.ebi.ac.uk/gwas/). In this study, we downloaded [SNP-sets for type 1 diabetes](https://www.ebi.ac.uk/gwas/efotraits/EFO_0001359). Then you can run R code file for filtering the GWAS Catalog data as below `CMD` command line:

- Bellow functions are run under the windows or linux command console environment.
- Instead of `[ ]`, you have to put the arguments `file path` or `value` by the options.
- Usage: `Rscript gwas.r [GWAS_file_path] [p-value_criteria]`

```cmd
Rscript src/gwas.r ^
	db/GWAS_EFO0001359.tsv ^
	5e-08
```

### > ldlink_dn.py and ldink_filt.r

To download **LDlink data** (version 3.3.0 12/24/2018) (Machiela et al, 2015, Bioinformatics, pmid 26139635), you can run `T1D_ldlink.py` as below `CMD` command line:

- To run the code, you need list of SNP RS IDs of dbSNP database as txt file
- Usage: `python ldlink.py [SNP_file_path.txt]`

```CMD
python src/ldlink_dn.py ^
	data/gwas_5e-08_129.tsv
```

> ...
> 129/129 = rs11580078
>   status_code = 200
>   line number = 950
>   file saved = db/ldlink/rs11580078.tsv
>
> Download process completed.
> Job time= 00:42:49

To filter the LDlink data by statistical criteria, both r<sup>2</sup> >0.6 and D'=1, you can run `T1D_ldlink.r` as below `CMD` command line:

- Usage: `Rscript ldlink.r [SNP_file_path.txt] [LDlink_data_folder_path] [LDlink_filter_option]`
- The `LDlink_filter_option` is a mandatory. Choose one of the following option numbers.
  1. `r2>0.6 or Dprime=1`
  2. `r2>0.6`
  3. `Dprime=1`
  4. `r2>0.6 and Dprime=1`

```CMD
Rscript ldlink_filt.r ^
	data/gwas_5e-08_129.tsv db/ldlink 4
```

> Input SNP list number = 129
>
> Error in read.table(as.character(snptb[i, 2]), header = T) :
>     more columns than column names
> In addition: Warning message:
> In read.table(as.character(snptb[i, 2]), header = T) :
>     incomplete final line found by readTableHeader on 'db/ldlink/rs75793288.tsv'
> NULL
> Filtering option, r2 > 0.6 and Dprime = 1 was chosen.
> ::Excluded no rsid elements = 10
>
> 1/3. Numbers of SNPs
> SNP Tier1 = 129
> SNP Tier2 = 1688
> SNP seed  = 1817
>
> 2/3. Generation of a result TSV file
> File write: data/seedSNP_1817_ldlink.tsv
>
> 3/3. Generation of a result BED file
> Table, rows= 1817 cols= 4
> File write: data/seedSNP_1817.bed
> Job done for 6.1 sec

* `data/seedSNP_5245_ldlink.tsv` -> **Supplementary Table 1** and **Supplementary Table 2**

### Q1. Generation of private SNP list (rsids) to BED file format?

To use bedtools later, you have to prepare SNP list as [bed format](https://genome.ucsc.edu/FAQ/FAQformat.html). If you have simple dbSNP rsid list, you can run `src/biomart_snp.r` for generate bed file. But you should check `NA` values and fill it manually.

* Usage: `Rscript src/biomart_snp.r [rsid_list_file_path]`

```CMD
Rscript src/biomart_snp.r ^
	data/seedSNP_5245_biomart.txt
```

> Input contents, rows= 5244 cols= 1
> Table, rows= 10280 cols= 4
>
> Job done for 12.7 sec

**[IMPORTANT]** Before you move to next step, we make sure that your `seedSNP_#.bed` file has no NA values.



## 2. Roadmap data download and filter

The [RoadMAP project](https://egg2.wustl.edu/roadmap/web_portal/imputed.html) provides epigenome annotations such as [12-mark/127-reference epigenome/25-state Imputation Based Chromatin State Model](https://egg2.wustl.edu/roadmap/data/byFileType/chromhmmSegmentations/ChmmModels/imputed12marks/jointModel/final/) by using ChromHMM algorithm (Ernst and Kellis, 2012, Nature Methods, pmid 22373907). We downloaded the 127 files by their [cell types](https://github.com/mdozmorov/genomerunner_web/wiki/Roadmap-cell-types) (e.g., `E001_25_imputed12marks_dense.bed.gz` and etc) using R code (`T1D_roadmap.r`). And then we filtered the data by [annotation code](https://egg2.wustl.edu/roadmap/web_portal/imputed.html) (see db/[roadmap] ) including 13_EnhA1, 14_EnhA2, 15_EnhAF, 16_EnhW1, 17_EnhW2, 18_EnhAc.

To download RoadMap data, you need to install `AnnotationHub` and `rtracklayer` in `BiocManaer` as below R code:

```R
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("AnnotationHub", version = "3.8")
BiocManager::install("rtracklayer", version = "3.8")
```

### > roadmap_dn.r

To download **RoadMap data**, you can run `roadmap_dn.r` as below `CMD` command line:

- The 127 cell type-specific RoadMap BED files will download at `db/roadmap` folder
- This process takes about ~30 min that depends on your download speed.

```CMD
Rscript src/roadmap_dn.r
```

### > roadmap_filt.r

To filter the RoadMap data by **Enhancers**, you can run `roadmap_filt.r` as below `CMD` command line:

- The result file would be saved as `data/roadmap_enh.bed`
- This process takes ~3 min that depends on your computer processor speed.

```CMD
Rscript src/roadmap_filt.r
```

### Q2. If you have memory problem..

When running the `roadmap_filt.r` function, it stop with not enough memory error, You can use `roadmap_filt_dtr.r` function for limited memory usage (~3.8 GB).

```CMD
Rscript src/roadmap_filt_dtr.r
```

### $ bedtools merge/ bedtools closest

To avoid multiple count of enhancers as well as to reduce file size and to achieve faster process, merge RoadMap enhancer information using a `BASH` tool `bedtools`. Here is the `BASH` pipeline for `bedtools sort` and `bedtools merge`. Then, to identify T1D SNPs occupied in RoadMap enhancers, you can use `BASH` tool `bedtools intersect` as below code:

- Compressed file size of `roadmap_enh.bed.gz` is >139 MB.
- Compressed file size of `roadmap_enh_merge.bed.gz` is about 3.7 MB.
- Removing NA values, `data/seedSNP_1817_bm.bed` file is updated version from the `data/seedSNP_1817.bed` file.

```SHELL
bedtools sort -i db/roadmap_enh.bed | bedtools merge -i stdin -c 1 -o count > db/roadmap_enh_merge.bed
bedtools sort -i data/seedSNP_1817_bm.bed | bedtools closest -d -a stdin -b db/roadmap_enh_merge.bed > data/roadmap_dist.tsv
```

### > src/bedtools_closest.r

To prioritize RoadMap enhancer occupied SNPs, you can run `src/bedtools_closestroadmap.r` as below `CMD` command line:

- `data/roadmap_dist_df.tsv` file is obtained that is for enhancer annotated file .
- `data/snp_484_roadmap_dist.bed` file is obtained that is for `BED` format file for USCS browser.
- Usage: `Rscript src/bedtools_closest.r [bedtools_closest_result_file_path] [double_line_result]`

```CMD
Rscript src/bedtools_closest.r ^
	data/roadmap_dist.tsv ^
	False
```

> Row number = 1817
>   Enhancer occupied by SNPs = 188
>   SNPs in RoadMap enhancers = 484
>
> File write: data/roadmap_dist_df.tsv
> File write: data/snp_484_roadmap_dist.bed
>
> Job done for 0.1 sec

* `data/roadmap_dist_df.tsv` -> **Supplementary Table 2**

### Q3. How about just use not merged roadmap_enh.bed file?

Instead of merge file, when you use original `db/roadmap_enh.bed` file, you can find a lot of duplicated enhancers regions.

```SHELL
bedtools sort -i db/roadmap_enh.bed | bedtools closest -d -a data/seedSNP_1817.bed -b stdin > data/roadmap_dist2.tsv
```



## 3. ENCODE ChIP-seq data download and filter

The **ENCODE ChIP-seq** for regulatory transcription factor binding site (Reg-TFBS) cluster data can downloaded <u>wgEncodeRegTfbsClusteredV3</u> data from [UCSC FTP](http://hgdownload.cse.ucsc.edu/goldenpath/hg19/encodeDCC/wgEncodeRegTfbsClustered/) (68 MB) or [bioconductor `data("wgEncodeTfbsV3")`](https://www.bioconductor.org/packages/devel/bioc/vignettes/ChIPpeakAnno/inst/doc/ChIPpeakAnno.html). Here, we assume having downloaded UCSC FTP file `wgEncodeRegTfbsClusteredV3.bed.gz` (81 MB).

```CMD
Rsciprt src/encode_dn.r
```

> Directory generated: db/encode/
> trying URL 'http://hgdownload.cse.ucsc.edu/goldenpath/hg19/encodeDCC/wgEncodeRegTfbsClustered/wgEncodeRegTfbsClusteredV3.bed.gz'
> Content type 'application/x-gzip' length 84986946 bytes (81.0 MB)
>
> downloaded 81.0 MB
>
> db/encode/wgEncodeRegTfbsClusteredV3.bed.gz
>
> Job done for 1.8 min

### $ bedtools merge | bedtools closest

To identify TFBS occupied SNPs, you can use `bedtools merge` and `bedtools closest` as following code:

- Merging the ENCODE TFBS data give you benefits such as avoiding multiple count of enhancers as well as reducing file size and achieving faster process

```SHELL
bedtools merge -i db/wgEncodeRegTfbsClusteredV3.bed.gz -c 1 -o count > db/encode_tfbs_merge.bed
bedtools sort -i data/seedSNP_1817_bm.bed | bedtools closest -d -a stdin -b db/encode_tfbs_merge.bed > data/encode_dist.tsv
```

### > src/bedtools_closest.r

To prioritize ENCODE Reg-TFBS occupied SNPs, you can run `src/bedtools_closestroadmap.r` as below `CMD` command line:

- `data/roadmap_dist_df.tsv` file is obtained that is for enhancer annotated file .
- `data/snp_enh_484.bed` file is obtained that is for `BED` format file for USCS browser.
- Usage: `Rscript src/bedtools_closest_roadmap.r [bedtools_closest_result_file_path] [double_line_result]`

```CMD
Rscript src/bedtools_closest.r ^
	data/encode_dist.tsv ^
	False
```

> Row number = 1817
>   Enhancer occupied by SNPs = 232
>   SNPs in RoadMap enhancers = 364
>
> File write: data/encode_dist_df.tsv
> File write: data/snp_364_encode_dist.bed
>
> Job done for 0.1 sec

* `data/encode_dist_df.tsv` -> **Supplementary Table 2**



## 4. Regulome DB data download and filter

The [**RegulomeDB**](http://www.regulomedb.org/index) provides [category scores for SNPs by evidences](http://www.regulomedb.org/help) (see `Regulome score.txt`), including eQTL, TF binding, matched TF motif, matched DNase Footprint, and DNase peak. In this study, we stringently filtered and used high-score (`≥ 2b`) SNPs for our study. Before you start, you can download the files from [Regulome DB download page](http://www.regulomedb.org/downloads). To make this process faster, you can convert the downloaded files to RDS format.

- `RegulomeDB.dbSNP132.Category1.txt.gz` (2 MB)
- `RegulomeDB.dbSNP132.Category2.txt.gz` (39.3 MB)
- Or you can download total dataset: `RegulomeDB.dbSNP141.txt.gz` (2.8 GB)

```CMD
Rscript src/regulome_dn.r
```

> In dir.create(file.path(dir)) : 'db\regulome' already exists
> trying URL 'http://www.regulomedb.org/downloads/RegulomeDB.dbSNP132.Category1.txt.gz'
> Content type 'application/gzip' length 2096454 bytes (2.0 MB)
>
> downloaded 2.0 MB
>
> trying URL 'http://www.regulomedb.org/downloads/RegulomeDB.dbSNP132.Category2.txt.gz'
> Content type 'application/gzip' length 41253483 bytes (39.3 MB)
>
> downloaded 39.3 MB
>
> Job process: 52.2 sec
> File write: db/regulome/RegulomeDB.dbSNP132.Category1.txt.gz.rds
> File write: db/regulome/RegulomeDB.dbSNP132.Category2.txt.gz.rds
> Job done for 1.2 min

Here we converted the download files to RDS format files to achieve fast loading speed. Use the RegulomeDB RDS files, you can filter and analyze the dataset by using `regulome.r` as following command line:

```CMD
Rscript src/regulome.r ^
	data/seedSNP_1817_bm.bed ^
	db/RegulomeDB.dbSNP132.Category1.txt.rds ^
	db/RegulomeDB.dbSNP132.Category2.txt.rds
```

> Input SNPs number = 1817
> Input regulome files = 2
>
> Read input file 1: db/RegulomeDB.dbSNP132.Category1.txt.rds
> Table row = 39432, col = 5
> Job process: 0.2 sec
>
> Read input file 2: db/RegulomeDB.dbSNP132.Category2.txt.rds
> Table row = 407796, col = 5
> Job process: 5 sec
>
> Regulome score >=2b, SNPs = 430528
> Functional motifs (2b-1f_only) = 395823
>   Regulome >=2b SNPs = 94
>   SNPs with functional motifs (1f_only-2b) = 45
>
> File write: data/regulome_94.tsv
> File write: data/snp_94_regulome2b.bed
> Job done for 5.3 sec

The result files are save at `data/` folder:

- `data/regulome_94.tsv` -> **Supplementary Table 2**
- `data/snp_94_regulome2b.bed`



### Venn analysis to identify core SNPs

Summary for SNPs with RoadMap annotation, ENCODE ChIP-seq, and RegulomeDB. This R code for Venn analysis uses **Bioconductor** `limma` R package. The installation of the `limma` package as below:

```R
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("limma", version = "3.8")
```

To prioritize the SNPs, you can run `venn.r` as below `CMD` command line with these files:

- `data/snp_484_roadmap_dist.bed`
- `data/snp_1253_encode_dist.bed`
- `data/snp_301_regulome2b.bed`

### > venn.r

```CMD
Rscript src/venn.r ^
	data/snp_484_roadmap_dist.bed ^
	data/snp_364_encode_dist.bed ^
	data/snp_94_regulome2b.bed
```

> package 'eulerr' was built under R version 3.6.2
>   snp_484_roadmap_dist
>   snp_364_encode_dist
>   snp_94_regulome2b
> 
>File write: data/venn.tsv
> 
> ** Euler fitting... done.
> 
>Figure draw: fig/euler_snp_484_roadmap_dist_snp_364_encode_dist.png
> File write: data/snp_26_core.bed

The result files are generated as below:

- `venn_tfbs.tsv`: binary SNP overlap table
- `snp_79_core.bed`

The result figure is generated as below:

![Venn analysis of 1817 SNPs](./fig/euler_snp_484_roadmap_dist_snp_364_encode_dist.png)

## 5. GTEx eQTL data download and filter

The [Genotype-Tissue Expression (GTEx)](https://gtexportal.org/home/) project is a public resource to study tissue-specific gene expression and their regulation by SNPs. GTEx version 7 includes 11,688 samples, 53 tissues and 714 donors. You can download [GTEx eQTL data](https://gtexportal.org/home/datasets) `GTEx_Analysis_v7_eQTL.tar.gz` (915 MB) and filter by statistical criteria `p < 3e-04`. The `GTEx_Analysis_v7_eQTL.tar.gz` compressed file includes:

- 48 files with `db/GTEx_Analysis_v7_eQTL/*.egenes.txt` extensions
- 48 files with `db/GTEx_Analysis_v7_eQTL/*.signif_variant_gene_pairs.txt` extensions

And we need SNP annotations to achieve Rsid for GTEx ids.

- `db/GTEx_Analysis_2016-01-15_v7_WholeGenomeSeq_635Ind_PASS_AB02_GQ20_HETX_MISS15_PLINKQC.lookup_table.txt.gz` (440 MB)
- [Nominal p-values](https://gtexportal.org/home/documentationPage) from GTEx data were generated for each variant-gene pair by testing the alternative hypothesis that the slope of a linear regression model between genotype and expression deviates from 0.

### > gtex_dn.r [and] gtex_filt.r

```CMD
Rscript src/gtex_dn.r
```

> (1/2) Download eQTL data
> trying URL 'https://storage.googleapis.com/gtex_analysis_v7/single_tissue_eqtl_data/GTEx_Analysis_v7_eQTL.tar.gz'
> Content type 'application/x-tar' length 959746583 bytes (915.3 MB)
>
> downloaded 915.3 MB
>
> File write: db/gtex_files.txt
>
> (2/2) Download SNPid annotation file
> trying URL 'https://storage.googleapis.com/gtex_analysis_v7/reference/GTEx_Analysis_2016-01-15_v7_WholeGenomeSeq_635Ind_PASS_AB02_GQ20_HETX_MISS15_PLINKQC.lookup_table.txt.gz'
> Content type 'application/gzip' length 461926948 bytes (440.5 MB)
>
> downloaded 440.5 MB
>
> Job done for 1.2 min

`> 50 min` The downloaded files were converted to RDS files.

```CMD
Rscript src/gtex_rds.r
```

>  - Loading GTEx BED files
>  - File reading...
>     (1/48) Adipose_Subcutaneous
>     ...
>   (48/48) Whole_Blood
> NULL
>  gte.df.pval_nominal
>  Min.   :0.000e+00
>  1st Qu.:0.000e+00
>  Median :1.111e-07
>  Mean   :7.398e-06
>  3rd Qu.:5.219e-06
>  Max.   :8.005e-04
>  - GTEx table, rows= 36781356 cols= 13
>  - BED file read complete. Job process: 9.3 min
>  - Loading annotation file - Annotation file read complete. Job process: 35.5 min
>  - Annotation file, rows= 40738696 cols= 7
>  - GTEx annotation, rows= 36781356 cols= 9
>
> Job process: 50.7 min
>
>  - Saving as RDS file..
>
> Job process: 52 min
>
> Job done for 52 min

To filter the GTEx data by p <5e-08, I executed following code:

```CMD
Rscript src/gtex_filt.r ^
	5e-08
```

> p-value threshold = 5e-08
>
> (1/3) Loading GTEx RDS file
>  - GTEx table, rows= 36781356 cols= 9
>  - BED file read complete. Job process: 53.3 sec
>
> (2/3) Filtering by nominal p-value
>  gte.sig.pval_nominal
>  Min.   :0.000e+00
>  1st Qu.:0.000e+00
>  Median :3.120e-12
>  Mean   :3.897e-09
>  3rd Qu.:1.429e-09
>  Max.   :5.000e-08
>  - GTEx significant, rows= 17113536 cols= 9
>  - Job process: 1.1 min
>
> File write: db/gtex_signif.tsv
> Job done for 3.8 min

The result file size are huge and the process takes long time (~50 min)

- `gtex_signif_5e-8.tsv.rds` (322 MB)
- This file was compressed by `zip` as three separated <100 MB files.
  1. `gtex_signif_5e-8.tsv.zip`
  2. `gtex_signif_5e-8.tsv.z01`
  3. `gtex_signif_5e-8.tsv.gz.z02`

To identify T1D SNPs 

```CMD
Rscript src/gtex.r ^
	data/seedSNP_1817_bm.bed ^
	db/gtex_signif_5e-8.tsv.rds
```

> Input SNPs number = 1,817
>
> (1/3) Loading GTEx significant file
>   gtex_signif_5e-8.tsv.rds: rows= 17,113,536 cols= 9
>   Job process: 23.9 sec
>
> (2/3) eQTL SNP filteration
>   Overlapped table, rows= 29,785 cols= 9
>   eQTL SNPs = 745
>   Associated genes = 159
> File write: data/gtex_5e-08_745.tsv
>
> (3/3) eQTL SNP BED file generation
>   GTEx SNP BED, rows= 745 cols= 4
>   eQTL SNPs = 745
> File write: data/snp_745_gtex.bed
>
> Job done for 26.7 sec

The result files of criteria 5e-08 are here:

- `gtex_5e-08_745.tsv` -> **Supplementary Table 3**
- `snp_745_gtex.bed`

### Roadmap & ENCODE

Preparing a file: `data_gwas/snp_140_roadmap_encode.bed`

```Bash
bedtools sort -i data/snp_484_roadmap_dist.bed | bedtools closest -d -a stdin -b db/encode_tfbs_merge.bed > data/roadmap_encode.tsv
```

```CMD
Rscript src/bedtools_closest.r ^
	data/roadmap_encode_dist.tsv ^
	False
```

### Venn analysis and overlap SNPs

To prioritize the eQTL SNPs among the 26 high-probability causal enhancer SNPs, you can run `venn.r` as below `CMD` command line with these files:

- `data/snp_140_roadmap_encode.bed` - Enhancer occupied SNP list **<- manually generated**
- `data/snp_26_core_tfbs.bed` - High-probability causal enhancer SNP list
- `data/snp_745_gtex.bed` - eQTL SNP list

### > venn.r

```CMD
Rscript src/venn.r ^
	data/snp_140_roadmap_encode.bed ^
	data/snp_26_core_tfbs.bed ^
	data/snp_745_gtex.bed
```

> package 'eulerr' was built under R version 3.6.2
>   snp_140_roadmap_encode
>   snp_26_core_tfbs
>   snp_745_gtex
>
> File write: data/venn.tsv
>
> ** Euler fitting... done.
>
> Figure draw: fig/euler_snp_140_roadmap_encode_snp_26_core_tfbs.png
> File write: data/snp_15_core.bed
> Job done for 0.2 sec

The result files are generated as below:

- `venn.tsv` > `venn_gtex.tsv`: binary SNP overlap table
- `snp_15_core.bed` > `snp_15_core_gtex.bed`: SNP `BED` format file

The result figure is generated as below:

![](./fig/euler_snp_140_roadmap_encode_snp_26_core_tfbs.png)



## 5-1. Nearest gene approach

### downloading Ensembl gene location data

To identify nearest genes from the eQTL SNPs, firstly you need to download gene location data from Ensembl database biomart (version=Grch37). 

```CMD
Rscript src/biomart_gene.r
```

>  - Ensembl table, rows= 63677 cols= 5
>  - File write: db/ensembl_gene_ann.tsv
>
>  - Filter result, rows= 57736 cols= 5
>  - File write: db/ensembl_gene.bed
>
> Job done for 17.1 sec

### $ bedtools closest

To identify nearest genes from the eQTL SNPs, you can use `bedtools merge` and `bedtools closest` as following `BASH` codes:

```SHELL
bedtools sort -i db/ensembl_gene.bed | bedtools closest -d -a data/seedSNP_5245_bm.bed -b stdin > data/seedSNP_nearest.tsv
```

```SHELL
bedtools sort -i db/ensembl_gene.bed | bedtools closest -d -a data/snp_2676_gtex.bed -b stdin > data/gtex_nearest.tsv
```

### > src/bedtools_closest.r

To prioritize RoadMap enhancer occupied SNPs, you can run `src/bedtools_closestroadmap.r` as below `CMD` command line:

- Usage: `Rscript src/bedtools_closest_gtex.r [bedtools_closest_result_file_path]`

```CMD
Rscript src/bedtools_closest_gtex.r data/seedSNP_nearest.tsv
```
> Row number = 2114
> Input SNPs = 1817
> Nearest genes = 175
>
> File write: data/seedSNP_nearest_df.tsv
>
> Job done for 0.1 sec


```CMD
Rscript src/bedtools_closest_gtex.r data/gtex_nearest.tsv
```

> Row number = 2909
> Input SNPs = 2676
> Nearest genes = 308
>
> File write: data/gtex_nearest_df.tsv
>
> Job done for 0.2 sec

* `data/gtex_nearest_df.tsv` -> **Supplementary Table 4**

### > src/gtex_overlap.r

To identify the eQTL SNPs occupied on TFBS binding enhancers, you can run `src/gtex_overlap.r` as below `CMD` command line:

```CMD
Rscript src/gtex_overlap.r ^
	data/snp_140_roadmap_encode.bed ^
	data/gtex_5e-08_745.tsv ^
	data/gtex_nearest_df.tsv
```

> 1. Read files..
>     data/snp_140_roadmap_encode.bed, rows= 140 cols= 4
>     data/gtex_5e-08_745.tsv, rows= 29785 cols= 9
>     data/gtex_nearest_df.tsv, rows= 788 cols= 7
>
>   SNPs= 745 Genes= 159 (Nearest= 44)
>
> 2. Overlap the two files..
>     TFBS overlap, rows= 5301 cols= 10
>    SNPs= 74 Genes= 94 (Nearest= 25)
> 
>File write: data/snp_74_gtex_enh.bed
> Job done for 0.3 sec



## 6. lncRNASNP2 data download and filter

Human SNPs located in long non-coding RNAs (lncRNAs) are archived in [**lncRNASNP2 database**](http://bioinfo.life.hust.edu.cn/lncRNASNP#!/). You can download these data at the [download page](http://bioinfo.life.hust.edu.cn/lncRNASNP#!/download):

- `lncRNASNP2_snplist.txt.gz` - **SNP list** includes the list of human SNPs in lncRNASNP database.
- `lncrnas.txt.gz` - **lncRNA list** includes the list of human lncRNAs in lncRNASNP database.
- `lncrna-diseases_experiment.txt.gz` - **Experimental validated lncRNA-associated diseases** includes all experiment validated lncRNA-associated diseases.
- `Rscript lncrnasnp.r [SNP_BED_file_path] [lncRNAsnp2_SNP_list_file_path] [lncRNAsnp2_lncRNA_list_file_path] [lncRNAsnp2_diseases_list_file_path]`

```CMD
Rscript lncrnasnp_dn.r
```

> 1: package 'data.table' was built under R version 3.5.2
>
> 2: package 'GenomeInfoDb' was built under R version 3.5.2
>
> trying URL 'http://bioinfo.life.hust.edu.cn/static/lncRNASNP2/downloads/snps_mod.txt'
>
> Content type 'text/plain; charset=GBK' length 477785336 bytes (455.7 MB)
>
> downloaded 455.7 MB
>
> trying URL 'http://bioinfo.life.hust.edu.cn/static/lncRNASNP2/downloads/lncrnas.txt'
>
> Content type 'text/plain; charset=GBK' length 7005411 bytes (6.7 MB)
>
> downloaded 6.7 MB
>
> trying URL 'http://bioinfo.life.hust.edu.cn/static/lncRNASNP2/downloads/lncRNA_associated_disease_experiment.txt'
>
> Content type 'text/plain; charset=GBK' length 31542 bytes (30 KB)
>
> downloaded 30 KB
>
> Job process: 1.4 min
>
> File write: db/lncRNASNP2_snplist.txt.rds
>
> File write: db/lncrnas.txt.rds
>
> File write: db/lncrna-diseases_experiment.txt.rds
>
> Job done for 2.2 min

To identify lncRNA overlapped longevity SNPs:

```CMD
Rscript src/lncrnasnp.r ^
	data/seedSNP_1817_bm.bed ^
	db/lncRNASNP2_snplist.txt.rds ^
	db/lncrnas.txt.rds db/lncrna-diseases_experiment.txt.rds
```

> 1. Read files..
>     data/seedSNP_1817_bm.bed...                   Job process: 0 sec
>     db/lncRNASNP2_snplist.txt.rds...              Job process: 13.3 sec
>     db/lncrnas.txt.rds...                         Job process: 13.5 sec
>     db/lncrna-diseases_experiment.txt.rds...      Job process: 13.5 sec
> 
> 
> | path                                  |     nrow | ncol |
> | :------------------------------------ | -------: | ---: |
> | data/seedSNP_1817_bm.bed              |     1817 |    4 |
> | db/lncRNASNP2_snplist.txt.rds         | 10205295 |    2 |
> | db/lncrnas.txt.rds                    |   141271 |    4 |
> | db/lncrna-diseases_experiment.txt.rds |      753 |    3 |
>  Job process: 13.6 sec
> 2. Overlapping lncRNA to my SNP list and binding annotation..
> 
> | lncRNA | SNPs |
> | -----: | ---: |
> |     42 |   78 |
> 
>  File write: data/snp_78_lncrnasnp.bed
> 
>3. Annotating SNPs in lncRNAs
>     File write: data/lncrnasnp_78.tsv
>
> Job done for 17.1 sec

- `data/snp_78_lncrnasnp.bed` - 78 SNPs `BED` file
- `data/lncrnasnp_78.tsv` -> **Supplementary Table 5**

### > venn.r

```CMD
Rscript src/venn.r ^
	data/snp_26_core_tfbs.bed ^
	data/snp_74_gtex_enh.bed ^
	data/snp_78_lncrnasnp.bed
```

> package 'eulerr' was built under R version 3.6.2
>   snp_26_core_tfbs
>   snp_74_gtex_enh
>   snp_78_lncrnasnp
>
> File write: data/venn.tsv
>
> ** Euler fitting... done.
>
> Figure draw: fig/euler_snp_26_core_tfbs_snp_74_gtex_enh.png
> File write: data/snp_2_core.bed
> Job done for 0.2 sec

- `venn.tsv` -> `venn_lncrnasnp.tsv`: binary SNP overlap table
- `vennCounts.tsv` -> `vennCounts_lncrnasnp.tsv`: overlapped SNP numbers
- `snp_2_core.bed` -> `snp_2_core_lncrnasnp.bed`: SNP `BED` format file

![](fig/euler_snp_26_core_tfbs_snp_74_gtex_enh.png)

### > Venn.r; Summary Table

```cmd
Rscript src/venn.r ^
	data/seedSNP_1817_bm.bed ^
	data/snp_745_gtex.bed ^
	data/snp_364_encode_dist.bed ^
	data/snp_484_roadmap_dist.bed ^
	data/snp_94_regulome2b.bed ^
	data/snp_78_lncrnasnp.bed
```

> package 'eulerr' was built under R version 3.6.2
>   seedSNP_1817_bm
>   snp_745_gtex
>   snp_364_encode_dist
>   snp_484_roadmap_dist
>   snp_94_regulome2b
>   snp_78_lncrnasnp
>
> [Message] Can't plot Venn diagram for more than 5 sets.
>
> File write: data/venn.tsv
>
> [Message] If you need snp_#_core.bed file, please input three groups.
> Job done for 0.1 sec

- `venn.tsv` -> `summary.tsv` : Summary file for this analysis.



# *Enrichment analysis



## Roadmap ChromHMM data

### Run permutation test

To calculate enrichment, run below command function at `bash`.

* T1D GWAS SNP file: `data/seedSNP_1817_bm.bed`
* Roadmap 127 BED files at `db/roadmap_bed` directory

```bash
Rscript src/enrich.r --permu \
    --gwassnp data/seedSNP_1817_bm.bed \
    --chrstatus db/roadmap_bed \
    --dbsource roadmap_bed \
    --permn 1000 \
    --out enrich
```

> ** Run perm_test function in enrich.r **
>
> * Gwas snp = [1] 1817    4
> * 127 files were found from db/roadmap_bed.
>
> 1 Load E001: roadmap_bed = [1] 933206      4
>   Run permTest: 25 annotations, [.........................] done. Job process: 8.3 min
> 2 Load E002: roadmap_bed = [1] 837982      4
>   Run permTest: 25 annotations, [.........................] done. Job process: 8.3 min
> ...
> 126 Load E128: roadmap_bed = [1] 835642      4
>   Run permTest: 25 annotations, [.........................] done. Job process: 7.7 min
> 127 Load E129: roadmap_bed = [1] 912567      4
>   Run permTest: 25 annotations, [.........................] done. Job process: 9.6 min
>
> * Write file: enrich/roadmap_bed-seedSNP_1817_bm-zscore.tsv
> * Write file: enrich/roadmap_bed-seedSNP_1817_bm-pval.tsv
> There were 50 or more warnings (use warnings() to see the first 50)
>
> Job done: 2021-01-08 15:59:44 for 17.1 hr



Run this function for 129 seed gwas SNPs.

```bash
Rscript src/enrich.r --permu \
    --gwassnp data/gwas_5e-08_129_hg19.bed \
    --chrstatus db/roadmap_bed \
    --dbsource roadmap_bed \
    --permn 1000 \
    --out enrich
```

> ** Run perm_test function in enrich.r **
>
> * Gwas snp = [1] 129   4
> * 127 files were found from db/roadmap_bed.
>
> 1 Load E001: roadmap_bed = 933206;  permTest for25 annotations, [.........................] done. Job process: 6.9 min
> 2 Load E002: roadmap_bed = 837982;  permTest for25 annotations, [.........................] done. Job process: 6.6 min
> ...
> 126 Load E128: roadmap_bed = 835642;  permTest for25 annotations, [.........................] done. Job process: 7.2 min
> 127 Load E129: roadmap_bed = 912567;  permTest for25 annotations, [.........................] done. Job process: 6.8 min
>
> * Write file: enrich/roadmap_bed-gwas_5e-08_129_hg19-permn_1000-zscore.tsv
> * Write file: enrich/roadmap_bed-gwas_5e-08_129_hg19-pval.tsv
> There were 50 or more warnings (use warnings() to see the first 50)
>
> Job done: 2021-01-09 13:11:18 for 14.1 hr



Run this function for 484 roadmap enhancer SNPs.

```bash
Rscript src/enrich.r --permu \
    --gwassnp data/snp_484_roadmap_dist.bed \
    --chrstatus db/roadmap_bed \
    --dbsource roadmap_bed \
    --permn 1000 \
    --out enrich
```

> ** Run perm_test function in enrich.r **
>
> * Gwas snp = [1] 484   4
> * 127 files were found from db/roadmap_bed.
>
> 1 E001: roadmap_bed = 933206; permTest - 25 annots = [............] done. Job process: 6.9 min
> 2 E002: roadmap_bed = 837982; permTest - 25 annots = [............] done. Job process: 6.3 min
> ...
> 126 E128: roadmap_bed = 835642; permTest - 25 annots = [............] done. Job process: 6.4 min
> 127 E129: roadmap_bed = 912567; permTest - 25 annots = [............] done. Job process: 6.8 min
>
> * Write file: enrich/roadmap_bed-snp_484_roadmap_dist-permn_1000-zscore.tsv
> * Write file: enrich/roadmap_bed-snp_484_roadmap_dist-permn_1000-pval.tsv
> * Write file: enrich/roadmap_bed-snp_484_roadmap_dist-permn_1000-overlap.tsv
> There were 50 or more warnings (use warnings() to see the first 50)
>
> Job done: 2021-01-10 17:23:15 for 13.7 hr



Run this function for 364 encode tfbs SNPs.

```bash
Rscript src/enrich.r --permu \
    --gwassnp data/snp_364_encode_dist.bed \
    --chrstatus db/roadmap_bed \
    --dbsource roadmap_bed \
    --permn 1000 \
    --out enrich
```

> ** Run perm_test function in enrich.r **
>
> * Gwas snp = [1] 364   4
> * 127 files were found from db/roadmap_bed.
>
> 1 E001: roadmap_bed = 933206; permTest - 25 annots [.........................] done. Job process: 1.9 min
> 2 E002: roadmap_bed = 837982; permTest - 25 annots [.........................] done. Job process: 1.7 min
> ...
> 126 E128: roadmap_bed = 835642; permTest - 25 annots [.........................] done. Job process: 1.4 min
> 127 E129: roadmap_bed = 912567; permTest - 25 annots [.........................] done. Job process: 1.4 min
>
> * Write file: enrich/roadmap_bed-snp_364_encode_dist-permn_100-zscore.tsv
> * Write file: enrich/roadmap_bed-snp_364_encode_dist-permn_100-pval.tsv
> * Write file: enrich/roadmap_bed-snp_364_encode_dist-permn_100-overlap.tsv



Run this function for 745 gtex eQTL SNPs.

```bash
Rscript src/enrich.r --permu \
    --gwassnp data/snp_745_gtex.bed \
    --chrstatus db/roadmap_bed \
    --dbsource roadmap_bed \
    --permn 1000 \
    --out enrich
```

> ** Run perm_test function in enrich.r **
>
> * Gwas snp = [1] 745   4
> * 127 files were found from db/roadmap_bed.
>
> 1 E001: roadmap_bed = 933206; permTest - 25 annots [.....] done. Job process: 1.5 min
> 2 E002: roadmap_bed = 837982; permTest - 25 annots [.....] done. Job process: 1.1 min
> ...
> 126 E128: roadmap_bed = 835642; permTest - 25 annots [.....] done. Job process: 1.5 min
> 127 E129: roadmap_bed = 912567; permTest - 25 annots [.....] done. Job process: 1.5 min
>
> * Write file: enrich/roadmap_bed-snp_745_gtex-permn_100-zscore.tsv
> * Write file: enrich/roadmap_bed-snp_745_gtex-permn_100-pval.tsv
> * Write file: enrich/roadmap_bed-snp_745_gtex-permn_100-overlap.tsv
>   There were 50 or more warnings (use warnings() to see the first 50)
>
> Job done: 2021-01-11 08:56:39 for 2.8 hr



### Draw heatmap

To draw heatmap by using the z-scores calculated from the permutation test, run below command function at `bash`.

```bash
Rscript src/enrich.r --heatmap \
    --pmdata enrich/roadmap_bed-seedSNP_1817_bm-permn1000-zscore.tsv \
    --meta db/roadmap_meta.tsv \
    --out enrich \
    --range -3,3 \
    --annot BLOOD,PANCREAS,THYMUS \
    --fileext png
```

> ** Run draw_heatmap function in enrich.r **
>
> * Permutation result table = [1]  25 128
> * Prepare table... done
>
> Save as enrich/roadmap_bed-seedSNP_1817_bm-zscore.png
>
> Job done: 2021-01-08 21:05:43 for 1.6 sec

```bash
Rscript src/enrich.r --heatmap \
    --pmdata enrich/roadmap_bed-gwas_5e-08_129_hg19-permn_1000-zscore.tsv \
    --meta db/roadmap_meta.tsv \
    --out enrich \
    --range -4,4 \
    --annot BLOOD,PANCREAS,THYMUS \
    --fileext png
```

> ** Run draw_heatmap function in enrich.r **
>
> * Permutation result table = [1]  25 128
> * Prepare table... done
>
> Save as enrich/roadmap_bed-gwas_5e-08_129_hg19-zscore.png
>
> Job done: 2021-01-08 21:06:40 for 1.7 sec



Draw heatmap for 484 roadmap enhancer SNPs:

```bash
Rscript src/enrich.r --heatmap \
    --pmdata enrich/roadmap_bed-snp_484_roadmap_dist-permn_100-zscore.tsv \
    --meta db/roadmap_meta.tsv \
    --out enrich \
    --range -3,3 \
    --annot BLOOD,PANCREAS,THYMUS \
    --fileext png
```

> ** Run draw_heatmap function in enrich.r **
>
> * Permutation result table = [1]  25 128
> * [Optional] Add meta-info. table = [1] 127   9
>
> Save as enrich/roadmap_bed-snp_484_roadmap_dist-permn_100-zscore.png
>
> Job done: 2021-01-09 21:37:26 for 3.4 sec



Draw heatmap for 364 tfbs residing SNPs:

```
Rscript src/enrich.r --heatmap \
    --pmdata enrich/roadmap_bed-snp_364_encode_dist-permn_100-zscore.tsv \
    --meta db/roadmap_meta.tsv \
    --out enrich \
    --range -4,4 \
    --annot BLOOD,PANCREAS,THYMUS \
    --fileext png
```

> ** Run draw_heatmap function in enrich.r **
>
> * Permutation result table = [1]  25 128
> * [Optional] Add meta-info. table = [1] 127   9
>
> Save as enrich/roadmap_bed-snp_364_encode_dist-permn_100-zscore.png



Draw heatmap for 745 GTEx eQTL SNPs:

```
Rscript src/enrich.r --heatmap \
    --pmdata enrich/roadmap_bed-snp_745_gtex-permn_100-zscore.tsv \
    --meta db/roadmap_meta.tsv \
    --out enrich \
    --range -4,4 \
    --annot BLOOD,PANCREAS,THYMUS \
    --fileext png
```

> ** Run draw_heatmap function in enrich.r **
>
> * Permutation result table = [1]  25 128
> * [Optional] Add meta-info. table = [1] 127   9
>
> Save as enrich/roadmap_bed-snp_745_gtex-permn_100-zscore.png
>
> Job done: 2021-01-11 21:51:22 for 4.2 sec



## ENCODE TFBS data

Data browse and download: http://genome.ucsc.edu/ENCODE/downloads.html

* [wgEncodeRegTfbsClusteredWithCellsV3.bed](http://hgdownload.soe.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeRegTfbsClustered/wgEncodeRegTfbsClusteredWithCellsV3.bed.gz): TFBS clusters together with input cell sources (BED 5+1 format: standard 5 fields of BED followed by comma-separated list of cell types)
* `encode_meta.tsv`: Cell line or tissue used as the source of experimental material [[ref](https://genome.ucsc.edu/ENCODE/cellTypes.html)].

### Split BED file by cell types

To split ENCODE TFBS file by cell types, run below code:

```bash
Rscript src/enrich.r --splittfbs \
    --tfbs db/wgEncodeRegTfbsClusteredWithCellsV3.bed \
    --out db/encode_bed
```

> ** Run split_tfbs function in enrich.r **
>
> * ENCODE TFBS table = [1] 4380443       6
> * 91 unique cell types are found.
>
> 1 A549: 380569 regions, filtering [..........] row = 380569, TFs = 24. Save: db/tfbs_cell/A549.bed; Job process: 2.7 min
> 2 AG04449:      49378 regions, filtering [..........] row = 49378, TFs = 1. Save: db/tfbs_cell/AG04449.bed; Job process: 23.2 sec
> 3 AG04450:      46148 regions, filtering [..........] row = 46148, TFs = 1. Save: db/tfbs_cell/AG04450.bed; Job process: 22.9 sec
> 4 AG09309:      43614 regions, filtering [..........] row = 43614, TFs = 1. Save: db/tfbs_cell/AG09309.bed; Job process: 21.2 sec
> 5 AG09319:      49235 regions, filtering [..........] row = 49235, TFs = 1. Save: db/tfbs_cell/AG09319.bed; Job process: 23.5 sec
> 6 AG10803:      47269 regions, filtering [..........] row = 47269, TFs = 1. Save: db/tfbs_cell/AG10803.bed; Job process: 23.2 sec
> 7 AoAF: 45395 regions, filtering [..........] row = 45395, TFs = 1. Save: db/tfbs_cell/AoAF.bed; Job process: 21.8 sec
> 8 BE2_C:        54249 regions, filtering [..........] row = 54249, TFs = 1. Save: db/tfbs_cell/BE2_C.bed; Job process: 25 sec
> 9 BJ:   44331 regions, filtering [..........] row = 44331, TFs = 1. Save: db/tfbs_cell/BJ.bed; Job process: 19.7 sec
> 10 Caco-2:      46436 regions, filtering [..........] row = 46436, TFs = 1. Save: db/tfbs_cell/Caco-2.bed; Job process: 21.2 sec
> 11 Dnd41:       51005 regions, filtering [..........] row = 51005, TFs = 2. Save: db/tfbs_cell/Dnd41.bed; Job process: 23.6 sec
> 12 ECC-1:       70393 regions, filtering [..........] row = 70393, TFs = 5. Save: db/tfbs_cell/ECC-1.bed; Job process: 31.4 sec
> 13 Fibrobl:     45188 regions, filtering [..........] row = 45188, TFs = 1. Save: db/tfbs_cell/Fibrobl.bed; Job process: 21 sec
> 14 GM06990:     45215 regions, filtering [..........] row = 45215, TFs = 1. Save: db/tfbs_cell/GM06990.bed; Job process: 22.3 sec
> 15 GM08714:     441 regions, filtering [..........] row = 441, TFs = 1. Save: db/tfbs_cell/GM08714.bed; Job process: 1.3 sec
> 16 GM10847:     17131 regions, filtering [..........] row = 17131, TFs = 2. Save: db/tfbs_cell/GM10847.bed; Job process: 8.5 sec
> 17 GM12801:     2882 regions, filtering [..........] row = 2882, TFs = 1. Save: db/tfbs_cell/GM12801.bed; Job process: 2.5 sec
> 18 GM12864:     46474 regions, filtering [..........] row = 46474, TFs = 1. Save: db/tfbs_cell/GM12864.bed; Job process: 25.8 sec
> 19 GM12865:     43737 regions, filtering [..........] row = 43737, TFs = 1. Save: db/tfbs_cell/GM12865.bed; Job process: 21.9 sec
> 20 GM12872:     46803 regions, filtering [..........] row = 46803, TFs = 1. Save: db/tfbs_cell/GM12872.bed; Job process: 23.5 sec
> 21 GM12873:     50614 regions, filtering [..........] row = 50614, TFs = 1. Save: db/tfbs_cell/GM12873.bed; Job process: 25.7 sec
> 22 GM12874:     37269 regions, filtering [..........] row = 37269, TFs = 1. Save: db/tfbs_cell/GM12874.bed; Job process: 19.5 sec
> 23 GM12875:     38921 regions, filtering [..........] row = 38921, TFs = 1. Save: db/tfbs_cell/GM12875.bed; Job process: 20.9 sec
> 24 GM12878:     1061985 regions, filtering [..........] row = 1061985, TFs = 76. Save: db/tfbs_cell/GM12878.bed; Job process: 8 min25 GM12891:     179010 regions, filtering [..........] row = 179010, TFs = 8. Save: db/tfbs_cell/GM12891.bed; Job process: 1.3 min
> 26 GM12892:     109701 regions, filtering [..........] row = 109701, TFs = 6. Save: db/tfbs_cell/GM12892.bed; Job process: 45.6 sec27 GM15510:     25455 regions, filtering [..........] row = 25455, TFs = 2. Save: db/tfbs_cell/GM15510.bed; Job process: 11.5 sec
> 28 GM18505:     23947 regions, filtering [..........] row = 23947, TFs = 2. Save: db/tfbs_cell/GM18505.bed; Job process: 10.9 sec
> 29 GM18526:     15322 regions, filtering [..........] row = 15322, TFs = 2. Save: db/tfbs_cell/GM18526.bed; Job process: 7.4 sec
> 30 GM18951:     28090 regions, filtering [..........] row = 28090, TFs = 2. Save: db/tfbs_cell/GM18951.bed; Job process: 12.9 sec
> 31 GM19099:     23281 regions, filtering [..........] row = 23281, TFs = 2. Save: db/tfbs_cell/GM19099.bed; Job process: 10.9 sec
> 32 GM19193:     20116 regions, filtering [..........] row = 20116, TFs = 2. Save: db/tfbs_cell/GM19193.bed; Job process: 9.5 sec
> 33 GM19238:     48654 regions, filtering [..........] row = 48654, TFs = 1. Save: db/tfbs_cell/GM19238.bed; Job process: 22.2 sec
> 34 GM19239:     40161 regions, filtering [..........] row = 40161, TFs = 1. Save: db/tfbs_cell/GM19239.bed; Job process: 18.8 sec
> 35 GM19240:     44690 regions, filtering [..........] row = 44690, TFs = 1. Save: db/tfbs_cell/GM19240.bed; Job process: 21.1 sec
> 36 Gliobla:     75041 regions, filtering [..........] row = 75041, TFs = 2. Save: db/tfbs_cell/Gliobla.bed; Job process: 35 sec
> 37 H1-hESC:     579539 regions, filtering [..........] row = 579539, TFs = 50. Save: db/tfbs_cell/H1-hESC.bed; Job process: 4.2 min38 HA-sp:       46196 regions, filtering [..........] row = 46196, TFs = 1. Save: db/tfbs_cell/HA-sp.bed; Job process: 22 sec
> 39 HAc: 45162 regions, filtering [..........] row = 45162, TFs = 1. Save: db/tfbs_cell/HAc.bed; Job process: 22.1 sec
> 40 HBMEC:       57865 regions, filtering [..........] row = 57865, TFs = 1. Save: db/tfbs_cell/HBMEC.bed; Job process: 27.1 sec
> 41 HCFaa:       41307 regions, filtering [..........] row = 41307, TFs = 1. Save: db/tfbs_cell/HCFaa.bed; Job process: 19.6 sec
> 42 HCM: 49849 regions, filtering [..........] row = 49849, TFs = 1. Save: db/tfbs_cell/HCM.bed; Job process: 23.3 sec
> 43 HCPEpiC:     60523 regions, filtering [..........] row = 60523, TFs = 1. Save: db/tfbs_cell/HCPEpiC.bed; Job process: 28.1 sec
> 44 HCT-116:     108478 regions, filtering [..........] row = 108478, TFs = 5. Save: db/tfbs_cell/HCT-116.bed; Job process: 47.9 sec
> 45 HEEpiC:      45911 regions, filtering [..........] row = 45911, TFs = 1. Save: db/tfbs_cell/HEEpiC.bed; Job process: 21.2 sec
> 46 HEK293:      110724 regions, filtering [..........] row = 83501, TFs = 5. Save: db/tfbs_cell/HEK293.bed; Job process: 41.9 sec 
> 47 HEK293-T-REx:        27223 regions, filtering [..........] row = 27223, TFs = 1. Save: db/tfbs_cell/HEK293-T-REx.bed; Job process:s: 12.6 sec
> 48 HFF: 44893 regions, filtering [..........] row = 34712, TFs = 1. Save: db/tfbs_cell/HFF.bed; Job process: 18.6 sec
> 49 HFF-Myc:     43263 regions, filtering [..........] row = 43263, TFs = 1. Save: db/tfbs_cell/HFF-Myc.bed; Job process: 21.6 sec
> 50 HL-60:       16646 regions, filtering [..........] row = 16646, TFs = 1. Save: db/tfbs_cell/HL-60.bed; Job process: 9.6 sec
> 51 HMEC:        59572 regions, filtering [..........] row = 59572, TFs = 2. Save: db/tfbs_cell/HMEC.bed; Job process: 28.9 sec       
> 52 HMF: 53800 regions, filtering [..........] row = 53800, TFs = 1. Save: db/tfbs_cell/HMF.bed; Job process: 28.7 sec
> 53 HPAF:        56233 regions, filtering [..........] row = 56233, TFs = 1. Save: db/tfbs_cell/HPAF.bed; Job process: 31 sec        
> 54 HPF: 45785 regions, filtering [..........] row = 45785, TFs = 1. Save: db/tfbs_cell/HPF.bed; Job process: 22.8 sec                 
> 55 HRE: 42047 regions, filtering [..........] row = 42047, TFs = 1. Save: db/tfbs_cell/HRE.bed; Job process: 20.5 sec                 
> 56 HRPEpiC:     52715 regions, filtering [..........] row = 52715, TFs = 1. Save: db/tfbs_cell/HRPEpiC.bed; Job process: 24.3 sec 
> 57 HSMM:        57052 regions, filtering [..........] row = 51799, TFs = 2. Save: db/tfbs_cell/HSMM.bed; Job process: 28.6 sec 
> 58 HSMMtube:    49162 regions, filtering [..........] row = 49162, TFs = 2. Save: db/tfbs_cell/HSMMtube.bed; Job process: 25.3 sec 
> 59 HUVEC:       199307 regions, filtering [..........] row = 199307, TFs = 8. Save: db/tfbs_cell/HUVEC.bed; Job process: 1.8 min 
> 60 HVMF:        46055 regions, filtering [..........] row = 46055, TFs = 1. Save: db/tfbs_cell/HVMF.bed; Job process: 29.6 sec        
> 61 HeLa-S3:     706062 regions, filtering [..........] row = 706062, TFs = 55. Save: db/tfbs_cell/HeLa-S3.bed; Job process: 6.1 min
> 62 HepG2:       977490 regions, filtering [..........] row = 977490, TFs = 59. Save: db/tfbs_cell/HepG2.bed; Job process: 7.7 min 
> 63 IMR90:       207461 regions, filtering [..........] row = 207461, TFs = 5. Save: db/tfbs_cell/IMR90.bed; Job process: 1.6 min     
> 64 K562:        1338658 regions, filtering [..........] row = 1338658, TFs = 100. Save: db/tfbs_cell/K562.bed; Job process: 11.9 min
> 65 MCF-7:       199597 regions, filtering [..........] row = 199597, TFs = 7. Save: db/tfbs_cell/MCF-7.bed; Job process: 2.1 min   
> 66 MCF10A-Er-Src:       240871 regions, filtering [..........] row = 240871, TFs = 5. Save: db/tfbs_cell/MCF10A-Er-Src.bed; Job process: 2.4 min
> 67 NB4: 108920 regions, filtering [..........] row = 108920, TFs = 4. Save: db/tfbs_cell/NB4.bed; Job process: 1.2 min
> 68 NH-A:        41714 regions, filtering [..........] row = 41714, TFs = 2. Save: db/tfbs_cell/NH-A.bed; Job process: 29.6 sec     
> 69 NHDF-Ad:     51628 regions, filtering [..........] row = 51628, TFs = 2. Save: db/tfbs_cell/NHDF-Ad.bed; Job process: 35.5 sec
> 70 NHDF-neo:    45555 regions, filtering [..........] row = 45555, TFs = 1. Save: db/tfbs_cell/NHDF-neo.bed; Job process: 31.7 sec 
> 71 NHEK:        73479 regions, filtering [..........] row = 73479, TFs = 3. Save: db/tfbs_cell/NHEK.bed; Job process: 49.3 sec     
> 72 NHLF:        46063 regions, filtering [..........] row = 46063, TFs = 2. Save: db/tfbs_cell/NHLF.bed; Job process: 31.4 sec     
> 73 NT2-D1:      7534 regions, filtering [..........] row = 7534, TFs = 3. Save: db/tfbs_cell/NT2-D1.bed; Job process: 5 sec        
> 74 Osteobl:     52928 regions, filtering [..........] row = 52928, TFs = 1. Save: db/tfbs_cell/Osteobl.bed; Job process: 32.6 sec  
> 75 PANC-1:      32124 regions, filtering [..........] row = 32124, TFs = 4. Save: db/tfbs_cell/PANC-1.bed; Job process: 19.5 sec   
> 76 PBDE:        30254 regions, filtering [..........] row = 30175, TFs = 2. Save: db/tfbs_cell/PBDE.bed; Job process: 17.2 sec     
> 77 PBDEFetal:   2125 regions, filtering [..........] row = 2125, TFs = 1. Save: db/tfbs_cell/PBDEFetal.bed; Job process: 2 sec     
> 78 PFSK-1:      40286 regions, filtering [..........] row = 40286, TFs = 4. Save: db/tfbs_cell/PFSK-1.bed; Job process: 23 sec     
> 79 ProgFib:     55609 regions, filtering [..........] row = 55609, TFs = 2. Save: db/tfbs_cell/ProgFib.bed; Job process: 35.8 sec  
> 80 RPTEC:       58235 regions, filtering [..........] row = 58235, TFs = 1. Save: db/tfbs_cell/RPTEC.bed; Job process: 40.8 sec    
> 81 Raji:        13738 regions, filtering [..........] row = 13738, TFs = 1. Save: db/tfbs_cell/Raji.bed; Job process: 8.8 sec      
> 82 SAEC:        42144 regions, filtering [..........] row = 42144, TFs = 1. Save: db/tfbs_cell/SAEC.bed; Job process: 33.4 sec
> 83 SH-SY5Y:     49309 regions, filtering [..........] row = 49309, TFs = 2. Save: db/tfbs_cell/SH-SY5Y.bed; Job process: 34.9 sec
> 84 SK-N-MC:     36916 regions, filtering [..........] row = 36916, TFs = 2. Save: db/tfbs_cell/SK-N-MC.bed; Job process: 25.4 sec
> 85 SK-N-SH:     269921 regions, filtering [..........] row = 77541, TFs = 4. Save: db/tfbs_cell/SK-N-SH.bed; Job process: 1.8 min
> 86 SK-N-SH_RA:  192380 regions, filtering [..........] row = 192380, TFs = 5. Save: db/tfbs_cell/SK-N-SH_RA.bed; Job process: 2 min
> 87 T-47D:       132022 regions, filtering [..........] row = 132022, TFs = 5. Save: db/tfbs_cell/T-47D.bed; Job process: 1.3 min
> 88 U2OS:        32015 regions, filtering [..........] row = 32015, TFs = 2. Save: db/tfbs_cell/U2OS.bed; Job process: 19.6 sec
> 89 U87: 30009 regions, filtering [..........] row = 30009, TFs = 2. Save: db/tfbs_cell/U87.bed; Job process: 18.1 sec
> 90 WERI-Rb-1:   49892 regions, filtering [..........] row = 49892, TFs = 1. Save: db/tfbs_cell/WERI-Rb-1.bed; Job process: 30.4 sec
> 91 WI-38:       31216 regions, filtering [..........] row = 31216, TFs = 1. Save: db/tfbs_cell/WI-38.bed; Job process: 20.9 sec

### Run permutation test

Run below code for T1D GWAS 129 SNP data:

```bash
Rscript src/enrich.r --permu \
    --gwassnp data/gwas_5e-08_129_hg19.bed \
    --chrstatus db/encode_bed \
    --dbsource encode_bed \
    --permn 100 \
    --out enrich
```

> ** Run perm_test function in enrich.r **
>
> * Gwas snp = [1] 129   4
> * 91 files were found from db/encode_bed.
>
> 1 A549: encode_bed = 380569; permTest - 24 annots [......] done. Job process: 1.8 min
> 2 AG04449: encode_bed = 49378; permTest - 1 annots [Error in if (i%%progress == 0) { : missing value where TRUE/FALSE needed       
> Calls: perm_test -> perm_test_calc -> lapply -> lapply -> FUN
> Execution halted
> root@f1e31bda3122:/git# Rscript src/enrich.r --permu     --gwassnp data/gwas_5e-08_129_hg19.bed     --chrstatus db/encode_bed     --dbsource encode_bed     --permn 100     --out enrich
>
> ** Run perm_test function in enrich.r **
>
> * Gwas snp = [1] 129   4
> * 91 files were found from db/encode_bed.
>
> 1 A549: encode_bed = 380569; permTest - 24 annots [......] done. Job process: 42.7 sec
> 2 AG04449: encode_bed = 49378; permTest - 1 annots [] done. Job process: 2.4 sec
> ...
> 90 WERI-Rb-1: encode_bed = 49892; permTest - 1 annots [] done. Job process: 4.7 sec
> 91 WI-38: encode_bed = 31216; permTest - 1 annots [] done. Job process: 3.7 sec
>
> * Write file: enrich/encode_bed-gwas_5e-08_129_hg19-permn_100-zscore.tsv
> * Write file: enrich/encode_bed-gwas_5e-08_129_hg19-permn_100-pval.tsv
> * Write file: enrich/encode_bed-gwas_5e-08_129_hg19-permn_100-overlap.tsv



Run below code for T1D candidate1817 SNP data:

```bash
Rscript src/enrich.r --permu \
    --gwassnp data/seedSNP_1817_bm.bed \
    --chrstatus db/encode_bed \
    --dbsource encode_bed \
    --permn 100 \
    --out enrich
```

> * Gwas snp = [1] 1817    4
> * 91 files were found from db/encode_bed.
>
> 1 A549: encode_bed = 380569; permTest - 24 annots [........................] done. Job process: 1.2 min
> 2 AG04449: encode_bed = 49378; permTest - 1 annots [.] done. Job process: 3.9 sec
> 3 AG04450: encode_bed = 46148; permTest - 1 annots [.] done. Job process: 4.8 sec
> 4 AG09309: encode_bed = 43614; permTest - 1 annots [.] done. Job process: 3.6 sec
> 5 AG09319: encode_bed = 49235; permTest - 1 annots [.] done. Job process: 4.8 sec
> 6 AG10803: encode_bed = 47269; permTest - 1 annots [.] done. Job process: 4.9 sec
> 7 AoAF: encode_bed = 45395; permTest - 1 annots [.] done. Job process: 4.3 sec
> 8 BE2_C: encode_bed = 54249; permTest - 1 annots [.] done. Job process: 4.5 sec
> 9 BJ: encode_bed = 44331; permTest - 1 annots [.] done. Job process: 4.6 sec
> 10 Caco-2: encode_bed = 46436; permTest - 1 annots [.] done. Job process: 6.5 sec
> 11 Dnd41: encode_bed = 51005; permTest - 2 annots [..] done. Job process: 7.5 sec
> 12 ECC-1: encode_bed = 70393; permTest - 5 annots [.....] done. Job process: 17.7 sec
> 13 Fibrobl: encode_bed = 45188; permTest - 1 annots [.] done. Job process: 4.8 sec
> 14 GM06990: encode_bed = 45215; permTest - 1 annots [.] done. Job process: 4.8 sec
> 15 GM08714: encode_bed = 441; permTest - 1 annots [.] done. Job process: 4.1 sec
> 16 GM10847: encode_bed = 17131; permTest - 2 annots [..] done. Job process: 6.4 sec
> 17 GM12801: encode_bed = 2882; permTest - 1 annots [.] done. Job process: 2.3 sec
> 18 GM12864: encode_bed = 46474; permTest - 1 annots [.] done. Job process: 3.6 sec
> 19 GM12865: encode_bed = 43737; permTest - 1 annots [.] done. Job process: 3.1 sec
> 20 GM12872: encode_bed = 46803; permTest - 1 annots [.] done. Job process: 3.4 sec
> 21 GM12873: encode_bed = 50614; permTest - 1 annots [.] done. Job process: 3.8 sec
> 22 GM12874: encode_bed = 37269; permTest - 1 annots [.] done. Job process: 3.2 sec
> 23 GM12875: encode_bed = 38921; permTest - 1 annots [.] done. Job process: 3.4 sec
> 24 GM12878: encode_bed = 1061985; permTest - 76 annots [............................................................................] done. Job process: 4.1 min
> 25 GM12891: encode_bed = 179010; permTest - 8 annots [........] done. Job process: 27.7 sec
> 26 GM12892: encode_bed = 109701; permTest - 6 annots [......] done. Job process: 19 sec
> 27 GM15510: encode_bed = 25455; permTest - 2 annots [..] done. Job process: 5.8 sec
> 28 GM18505: encode_bed = 23947; permTest - 2 annots [..] done. Job process: 5.9 sec
> 29 GM18526: encode_bed = 15322; permTest - 2 annots [..] done. Job process: 7 sec
> 30 GM18951: encode_bed = 28090; permTest - 2 annots [..] done. Job process: 6 sec
> 31 GM19099: encode_bed = 23281; permTest - 2 annots [..] done. Job process: 5.8 sec
> 32 GM19193: encode_bed = 20116; permTest - 2 annots [..] done. Job process: 5.3 sec
> 33 GM19238: encode_bed = 48654; permTest - 1 annots [.] done. Job process: 5 sec
> 34 GM19239: encode_bed = 40161; permTest - 1 annots [.] done. Job process: 4.5 sec
> 35 GM19240: encode_bed = 44690; permTest - 1 annots [.] done. Job process: 3.6 sec
> 36 Gliobla: encode_bed = 75041; permTest - 2 annots [..] done. Job process: 8.1 sec
> 37 H1-hESC: encode_bed = 579539; permTest - 50 annots [..................................................] done. Job process: 2.6 min
> 38 HA-sp: encode_bed = 46196; permTest - 1 annots [.] done. Job process: 4.6 sec
> 39 HAc: encode_bed = 45162; permTest - 1 annots [.] done. Job process: 3.9 sec
> 40 HBMEC: encode_bed = 57865; permTest - 1 annots [.] done. Job process: 4.9 sec
> 41 HCFaa: encode_bed = 41307; permTest - 1 annots [.] done. Job process: 4.1 sec
> 42 HCM: encode_bed = 49849; permTest - 1 annots [.] done. Job process: 3.9 sec
> 43 HCPEpiC: encode_bed = 60523; permTest - 1 annots [.] done. Job process: 4.3 sec
> 44 HCT-116: encode_bed = 108478; permTest - 5 annots [.....] done. Job process: 16.7 sec
> 45 HEEpiC: encode_bed = 45911; permTest - 1 annots [.] done. Job process: 4.7 sec
> 46 HEK293-T-REx: encode_bed = 27223; permTest - 1 annots [.] done. Job process: 4.3 sec
> 47 HEK293: encode_bed = 83501; permTest - 5 annots [.....] done. Job process: 17 sec
> 48 HFF-Myc: encode_bed = 43263; permTest - 1 annots [.] done. Job process: 4.3 sec
> 49 HFF: encode_bed = 34712; permTest - 1 annots [.] done. Job process: 3.7 sec
> 50 HL-60: encode_bed = 16646; permTest - 1 annots [.] done. Job process: 4.8 sec
> 51 HMEC: encode_bed = 59572; permTest - 2 annots [..] done. Job process: 7.5 sec
> 52 HMF: encode_bed = 53800; permTest - 1 annots [.] done. Job process: 4.6 sec
> 53 HPAF: encode_bed = 56233; permTest - 1 annots [.] done. Job process: 5.8 sec
> 54 HPF: encode_bed = 45785; permTest - 1 annots [.] done. Job process: 4.1 sec
> 55 HRE: encode_bed = 42047; permTest - 1 annots [.] done. Job process: 4.2 sec
> 56 HRPEpiC: encode_bed = 52715; permTest - 1 annots [.] done. Job process: 5 sec
> 57 HSMM: encode_bed = 51799; permTest - 2 annots [..] done. Job process: 7.2 sec
> 58 HSMMtube: encode_bed = 49162; permTest - 2 annots [..] done. Job process: 8 sec
> 59 HUVEC: encode_bed = 199307; permTest - 8 annots [........] done. Job process: 31.9 sec
> 60 HVMF: encode_bed = 46055; permTest - 1 annots [.] done. Job process: 3.8 sec
> 61 HeLa-S3: encode_bed = 706062; permTest - 55 annots [.......................................................] done. Job process: 
> 3 min
> 62 HepG2: encode_bed = 977490; permTest - 59 annots [...........................................................] done. Job process: 4.4 min
> 63 IMR90: encode_bed = 207461; permTest - 5 annots [.....] done. Job process: 21.2 sec
> 64 K562: encode_bed = 1338658; permTest - 100 annots [....................................................................................................] done. Job process: 6.1 min
> 65 MCF-7: encode_bed = 199597; permTest - 7 annots [.......] done. Job process: 32.2 sec
> 66 MCF10A-Er-Src: encode_bed = 240871; permTest - 5 annots [.....] done. Job process: 24.6 sec
> 67 NB4: encode_bed = 108920; permTest - 4 annots [....] done. Job process: 17.2 sec
> 68 NH-A: encode_bed = 41714; permTest - 2 annots [..] done. Job process: 7.8 sec
> 69 NHDF-Ad: encode_bed = 51628; permTest - 2 annots [..] done. Job process: 8.4 sec
> 70 NHDF-neo: encode_bed = 45555; permTest - 1 annots [.] done. Job process: 4.3 sec
> 71 NHEK: encode_bed = 73479; permTest - 3 annots [...] done. Job process: 12.2 sec
> 72 NHLF: encode_bed = 46063; permTest - 2 annots [..] done. Job process: 7.8 sec
> 73 NT2-D1: encode_bed = 7534; permTest - 3 annots [...] done. Job process: 8.6 sec
> 74 Osteobl: encode_bed = 52928; permTest - 1 annots [.] done. Job process: 4.3 sec
> 75 PANC-1: encode_bed = 32124; permTest - 4 annots [....] done. Job process: 12.1 sec
> 76 PBDE: encode_bed = 30175; permTest - 2 annots [..] done. Job process: 8.1 sec
> 77 PBDEFetal: encode_bed = 2125; permTest - 1 annots [.] done. Job process: 3.7 sec
> 78 PFSK-1: encode_bed = 40286; permTest - 4 annots [....] done. Job process: 14.9 sec
> 79 ProgFib: encode_bed = 55609; permTest - 2 annots [..] done. Job process: 9.7 sec
> 80 RPTEC: encode_bed = 58235; permTest - 1 annots [.] done. Job process: 5.1 sec
> 81 Raji: encode_bed = 13738; permTest - 1 annots [.] done. Job process: 4.1 sec
> 82 SAEC: encode_bed = 42144; permTest - 1 annots [.] done. Job process: 5.1 sec
> 83 SH-SY5Y: encode_bed = 49309; permTest - 2 annots [..] done. Job process: 7.8 sec
> 84 SK-N-MC: encode_bed = 36916; permTest - 2 annots [..] done. Job process: 10.3 sec
> 85 SK-N-SH: encode_bed = 77541; permTest - 4 annots [....] done. Job process: 15 sec
> 86 SK-N-SH_RA: encode_bed = 192380; permTest - 5 annots [.....] done. Job process: 28.6 sec
> 87 T-47D: encode_bed = 132022; permTest - 5 annots [.....] done. Job process: 23 sec
> 88 U2OS: encode_bed = 32015; permTest - 2 annots [..] done. Job process: 8.4 sec
> 89 U87: encode_bed = 30009; permTest - 2 annots [..] done. Job process: 6.8 sec
> 90 WERI-Rb-1: encode_bed = 49892; permTest - 1 annots [.] done. Job process: 4.5 sec
> 91 WI-38: encode_bed = 31216; permTest - 1 annots [.] done. Job process: 4.8 sec
>
> * Write file: enrich/encode_bed-seedSNP_1817_bm-permn_100-zscore.tsv
> * Write file: enrich/encode_bed-seedSNP_1817_bm-permn_100-pval.tsv
> * Write file: enrich/encode_bed-seedSNP_1817_bm-permn_100-overlap.tsv



### Draw heatmap

Run below code for T1D GWAS 129 SNP data:

```bash
Rscript src/enrich.r --heatmap \
    --pmdata enrich/encode_bed-gwas_5e-08_129_hg19-permn_100-zscore.tsv \
    --meta db/encode_meta.tsv \
    --out enrich \
    --range -3,3 \
    --annot blood,pancreas \
    --fileext png
```

> ** Run draw_heatmap function in enrich.r **
>
> * Permutation result table = [1] 161  92
> * [Optional] Add meta-info. table = [1] 91  7
>
> Save as enrich/encode_bed-gwas_5e-08_129_hg19-permn_100-zscore.png

Run below code for T1D candidate1817 SNP data:

```bash
Rscript src/enrich.r --heatmap \
    --pmdata enrich/encode_bed-seedSNP_1817_bm-permn_100-zscore.tsv \
    --meta db/encode_meta.tsv \
    --out enrich \
    --range -3,3 \
    --annot blood,pancreas \
    --fileext png
```

> ** Run draw_heatmap function in enrich.r **
>
> * Permutation result table = [1] 161  92
> * [Optional] Add meta-info. table = [1] 91  7
>
> Save as enrich/encode_bed-seedSNP_1817_bm-permn_100-zscore.png



## GTEx eQTL data

To split ENCODE TFBS file by cell types, run below code:

```bash
Rscript src/enrich.r --splitgtex \
    --gtex db/gtex_signif_5e-8.tsv.rds \
    --out db/gtex_tsv
```

> ** Run split_gtex function in enrich.r **
>
> * GTEx table = [1] 17113536        9
> * 48 unique tissue types are found.
>
> 1 Adipose_Subcutaneous: 732572 pairs, SNPs = 479346, genes = 6458. Save: db/gtex_tsv/Adipose_Subcutaneous.tsv; Job process: 14.3 sec
> 2 Adipose_Visceral_Omentum:     488060 pairs, SNPs = 343596, genes = 4640. Save: db/gtex_tsv/Adipose_Visceral_Omentum.tsv; Job process: 10.3 sec
> 3 Adrenal_Gland:        251774 pairs, SNPs = 178087, genes = 3013. Save: db/gtex_tsv/Adrenal_Gland.tsv; Job process: 5.4 sec
> 4 Artery_Aorta: 497489 pairs, SNPs = 347636, genes = 4889. Save: db/gtex_tsv/Artery_Aorta.tsv; Job process: 9.1 sec
> 5 Artery_Coronary:      170735 pairs, SNPs = 121719, genes = 2055. Save: db/gtex_tsv/Artery_Coronary.tsv; Job process: 4.4 sec
> 6 Artery_Tibial:        745555 pairs, SNPs = 503781, genes = 6548. Save: db/gtex_tsv/Artery_Tibial.tsv; Job process: 13.2 sec
> 7 Brain_Amygdala:       75907 pairs, SNPs = 47669, genes = 892. Save: db/gtex_tsv/Brain_Amygdala.tsv; Job process: 2.4 sec
> 8 Brain_Anterior_cingulate_cortex_BA24: 130324 pairs, SNPs = 88664, genes = 1568. Save: db/gtex_tsv/Brain_Anterior_cingulate_cortex_BA24.tsv; Job process: 3.2 sec
> 9 Brain_Caudate_basal_ganglia:  193909 pairs, SNPs = 136644, genes = 2305. Save: db/gtex_tsv/Brain_Caudate_basal_ganglia.tsv; Job process: 4.6 sec
> 10 Brain_Cerebellar_Hemisphere: 255466 pairs, SNPs = 161056, genes = 2902. Save: db/gtex_tsv/Brain_Cerebellar_Hemisphere.tsv; Job process: 5.9 sec
> 11 Brain_Cerebellum:    364347 pairs, SNPs = 237434, genes = 4156. Save: db/gtex_tsv/Brain_Cerebellum.tsv; Job process: 7.5 sec
> 12 Brain_Cortex:        210939 pairs, SNPs = 148013, genes = 2581. Save: db/gtex_tsv/Brain_Cortex.tsv; Job process: 5.1 sec
> 13 Brain_Frontal_Cortex_BA9:    151424 pairs, SNPs = 103230, genes = 1964. Save: db/gtex_tsv/Brain_Frontal_Cortex_BA9.tsv; Job process: 4.3 sec
> 14 Brain_Hippocampus:   105435 pairs, SNPs = 71524, genes = 1250. Save: db/gtex_tsv/Brain_Hippocampus.tsv; Job process: 2.6 sec
> 15 Brain_Hypothalamus:  110627 pairs, SNPs = 73904, genes = 1321. Save: db/gtex_tsv/Brain_Hypothalamus.tsv; Job process: 3.6 sec
> 16 Brain_Nucleus_accumbens_basal_ganglia:       167615 pairs, SNPs = 116190, genes = 2016. Save: db/gtex_tsv/Brain_Nucleus_accumbens_basal_ganglia.tsv; Job process: 4.6 sec
> 17 Brain_Putamen_basal_ganglia: 132821 pairs, SNPs = 90976, genes = 1625. Save: db/gtex_tsv/Brain_Putamen_basal_ganglia.tsv; Job process: 4.1 sec
> 18 Brain_Spinal_cord_cervical_c-1:      79575 pairs, SNPs = 49284, genes = 1033. Save: db/gtex_tsv/Brain_Spinal_cord_cervical_c-1.tsv; Job process: 2.5 sec
> 19 Brain_Substantia_nigra:      62312 pairs, SNPs = 38290, genes = 738. Save: db/gtex_tsv/Brain_Substantia_nigra.tsv; Job process: 
> 2.1 sec
> 20 Breast_Mammary_Tissue:       359737 pairs, SNPs = 248378, genes = 3630. Save: db/gtex_tsv/Breast_Mammary_Tissue.tsv; Job process: 7.3 sec
> 21 Cells_EBV-transformed_lymphocytes:   133183 pairs, SNPs = 89175, genes = 1690. Save: db/gtex_tsv/Cells_EBV-transformed_lymphocytes.tsv; Job process: 3.5 sec
> 22 Cells_Transformed_fibroblasts:       593793 pairs, SNPs = 416044, genes = 5771. Save: db/gtex_tsv/Cells_Transformed_fibroblasts.tsv; Job process: 11.3 sec
> 23 Colon_Sigmoid:       322044 pairs, SNPs = 221996, genes = 3431. Save: db/gtex_tsv/Colon_Sigmoid.tsv; Job process: 7.6 sec
> 24 Colon_Transverse:    380068 pairs, SNPs = 261772, genes = 3929. Save: db/gtex_tsv/Colon_Transverse.tsv; Job process: 8 sec
> 25 Esophagus_Gastroesophageal_Junction: 337678 pairs, SNPs = 235643, genes = 3562. Save: db/gtex_tsv/Esophagus_Gastroesophageal_Junction.tsv; Job process: 8.2 sec
> 26 Esophagus_Mucosa:    726202 pairs, SNPs = 468871, genes = 6378. Save: db/gtex_tsv/Esophagus_Mucosa.tsv; Job process: 14.7 sec
> 27 Esophagus_Muscularis:        678726 pairs, SNPs = 456865, genes = 6025. Save: db/gtex_tsv/Esophagus_Muscularis.tsv; Job process: 14.2 sec
> 28 Heart_Atrial_Appendage:      402932 pairs, SNPs = 287951, genes = 4053. Save: db/gtex_tsv/Heart_Atrial_Appendage.tsv; Job process: 8.7 sec
> 29 Heart_Left_Ventricle:        365179 pairs, SNPs = 260492, genes = 3713. Save: db/gtex_tsv/Heart_Left_Ventricle.tsv; Job process: 8.3 sec
> 30 Liver:       157448 pairs, SNPs = 111629, genes = 1877. Save: db/gtex_tsv/Liver.tsv; Job process: 4.3 sec
> 31 Lung:        694407 pairs, SNPs = 453798, genes = 6102. Save: db/gtex_tsv/Lung.tsv; Job process: 13.2 sec
> 32 Minor_Salivary_Gland:        66481 pairs, SNPs = 45632, genes = 941. Save: db/gtex_tsv/Minor_Salivary_Gland.tsv; Job process: 3 sec
> 33 Muscle_Skeletal:     724430 pairs, SNPs = 485201, genes = 5791. Save: db/gtex_tsv/Muscle_Skeletal.tsv; Job process: 13.6 sec
> 34 Nerve_Tibial:        890851 pairs, SNPs = 571923, genes = 7782. Save: db/gtex_tsv/Nerve_Tibial.tsv; Job process: 19 sec
> 35 Ovary:       134741 pairs, SNPs = 89668, genes = 1670. Save: db/gtex_tsv/Ovary.tsv; Job process: 4.4 sec
> 36 Pancreas:    345819 pairs, SNPs = 247140, genes = 3764. Save: db/gtex_tsv/Pancreas.tsv; Job process: 7.4 sec
> 37 Pituitary:   262721 pairs, SNPs = 176991, genes = 2846. Save: db/gtex_tsv/Pituitary.tsv; Job process: 5.8 sec
> 38 Prostate:    143429 pairs, SNPs = 94202, genes = 1647. Save: db/gtex_tsv/Prostate.tsv; Job process: 4 sec
> 39 Skin_Not_Sun_Exposed_Suprapubic:     607740 pairs, SNPs = 405351, genes = 5689. Save: db/gtex_tsv/Skin_Not_Sun_Exposed_Suprapubic.tsv; Job process: 12.8 sec
> 40 Skin_Sun_Exposed_Lower_leg:  821197 pairs, SNPs = 532886, genes = 7065. Save: db/gtex_tsv/Skin_Sun_Exposed_Lower_leg.tsv; Job process: 16 sec
> 41 Small_Intestine_Terminal_Ileum:      137388 pairs, SNPs = 95423, genes = 1734. Save: db/gtex_tsv/Small_Intestine_Terminal_Ileum.tsv; Job process: 4.6 sec
> 42 Spleen:      232394 pairs, SNPs = 165671, genes = 3025. Save: db/gtex_tsv/Spleen.tsv; Job process: 4.9 sec
> 43 Stomach:     305539 pairs, SNPs = 216267, genes = 3157. Save: db/gtex_tsv/Stomach.tsv; Job process: 6.7 sec
> 44 Testis:      677525 pairs, SNPs = 473444, genes = 7154. Save: db/gtex_tsv/Testis.tsv; Job process: 12.8 sec
> 45 Thyroid:     997164 pairs, SNPs = 614306, genes = 8063. Save: db/gtex_tsv/Thyroid.tsv; Job process: 17.7 sec
> 46 Uterus:      95351 pairs, SNPs = 58512, genes = 1100. Save: db/gtex_tsv/Uterus.tsv; Job process: 3.5 sec
> 47 Vagina:      91635 pairs, SNPs = 59912, genes = 1015. Save: db/gtex_tsv/Vagina.tsv; Job process: 3.1 sec
> 48 Whole_Blood: 500848 pairs, SNPs = 344421, genes = 4847. Save: db/gtex_tsv/Whole_Blood.tsv; Job process: 9.6 sec
>
> Job done: 2021-01-10 06:53:10 for 7 min



### Run permutation test

Run below code for T1D candidate 1817 SNPs:

```bash
f=(seedSNP_1817_bm.bed gwas_5e-08_129_hg19.bed snp_364_encode_dist.bed snp_484_roadmap_dist.bed snp_745_gtex.bed)

for i in {2..4}
do
Rscript src/enrich.r --gtex_perm \
	--gwas_snp data/${f[i]} \
    --gtex_base db/gtex_signif_5e-8.db \
    --gtex_median_tpm db/gtex_analysis_v8_rnaseq_gene_median_tpm_ensgid.gct.rds \
    --perm_n 1000 \
    --out enrich \
    > enrich/log_${f[i]}.txt
done
```

> ** Run gtex_perm_test function in enrich.r **
>
> * Input SNPs = [1] 1817    4
> * db/gtex_signif_5e-8.db, Ensgid = [1] 21404
> * GTEx gene median tpm GCT = 56156 -> transform -> Filt gene = 20691
> * Filtering eQTLs by SNP: genes = 159, tissues = 48 [.........] done.
>
> 1 Adipose_Subcutaneous: genes = 43 / 6458, exp. = 6231
> 2 Adipose_Visceral_Omentum: genes = 32 / 4640, exp. = 4450
> 3 Adrenal_Gland: genes = 22 / 3013, exp. = 2895
> 4 Artery_Aorta: genes = 36 / 4889, exp. = 4709
> 5 Artery_Coronary: genes = 14 / 2055, exp. = 1966
> 6 Artery_Tibial: genes = 46 / 6548, exp. = 6324
> 7 Brain_Amygdala: genes = 7 / 892, exp. = 855
> 8 Brain_Anterior_cingulate_cortex_BA24: genes = 9 / 1568, exp. = 1499
> 9 Brain_Caudate_basal_ganglia: genes = 16 / 2305, exp. = 2203
> 10 Brain_Cerebellar_Hemisphere: genes = 9 / 2902, exp. = 2780
> 11 Brain_Cerebellum: genes = 18 / 4156, exp. = 3991
> 12 Brain_Cortex: genes = 13 / 2581, exp. = 2480
> 13 Brain_Frontal_Cortex_BA9: genes = 7 / 1964, exp. = 1885
> 14 Brain_Hippocampus: genes = 8 / 1250, exp. = 1194
> 15 Brain_Hypothalamus: genes = 7 / 1321, exp. = 1251
> 16 Brain_Nucleus_accumbens_basal_ganglia: genes = 10 / 2016, exp. = 1940
> 17 Brain_Putamen_basal_ganglia: genes = 6 / 1625, exp. = 1560
> 18 Brain_Spinal_cord_cervical_c-1: genes = 4 / 1033, exp. = 0
> 19 Brain_Substantia_nigra: genes = 3 / 738, exp. = 713
> 20 Breast_Mammary_Tissue: genes = 22 / 3630, exp. = 3476
> 21 Cells_EBV-transformed_lymphocytes: genes = 14 / 1690, exp. = 0
> 22 Cells_Transformed_fibroblasts: genes = 28 / 5771, exp. = 0
> 23 Colon_Sigmoid: genes = 26 / 3431, exp. = 3285
> 24 Colon_Transverse: genes = 33 / 3929, exp. = 3755
> 25 Esophagus_Gastroesophageal_Junction: genes = 23 / 3562, exp. = 3427
> 26 Esophagus_Mucosa: genes = 47 / 6378, exp. = 6156
> 27 Esophagus_Muscularis: genes = 49 / 6025, exp. = 5811
> 28 Heart_Atrial_Appendage: genes = 29 / 4053, exp. = 3898
> 29 Heart_Left_Ventricle: genes = 26 / 3713, exp. = 3571
> 30 Liver: genes = 18 / 1877, exp. = 1788
> 31 Lung: genes = 51 / 6102, exp. = 5874
> 32 Minor_Salivary_Gland: genes = 6 / 941, exp. = 889
> 33 Muscle_Skeletal: genes = 43 / 5791, exp. = 5609
> 34 Nerve_Tibial: genes = 48 / 7782, exp. = 7495
> 35 Ovary: genes = 6 / 1670, exp. = 1590
> 36 Pancreas: genes = 30 / 3764, exp. = 3627
> 37 Pituitary: genes = 16 / 2846, exp. = 2718
> 38 Prostate: genes = 13 / 1647, exp. = 1555
> 39 Skin_Not_Sun_Exposed_Suprapubic: genes = 42 / 5689, exp. = 5483
> 40 Skin_Sun_Exposed_Lower_leg: genes = 54 / 7065, exp. = 6816
> 41 Small_Intestine_Terminal_Ileum: genes = 14 / 1734, exp. = 1647
> 42 Spleen: genes = 18 / 3025, exp. = 2893
> 43 Stomach: genes = 27 / 3157, exp. = 3028
> 44 Testis: genes = 40 / 7154, exp. = 6896
> 45 Thyroid: genes = 58 / 8063, exp. = 7757
> 46 Uterus: genes = 4 / 1100, exp. = 1039
> 47 Vagina: genes = 3 / 1015, exp. = 947
> 48 Whole_Blood: genes = 43 / 4847, exp. = 4680
>
> Write file: enrich/gtex-seedSNP_1817_bm-permn_1000.tsv
> Write file: enrich/gtex-seedSNP_1817_bm-permn_1000.rds

```bash
Rscript src/enrich.r --gtex_pm_plot \
	--fgsea_rds enrich/gtex-seedSNP_1817_bm-permn_1000.rds \
	--out enrich
```

> 



Run below code for T1D GWAS 129 SNP data:

```bash
Rscript src/enrich.r --gtex_perm \
	--gwas_snp data/gwas_5e-08_129_hg19.bed \
    --gtex_base db/gtex_signif_5e-8.db \
    --gtex_median_tpm db/gtex_analysis_v8_rnaseq_gene_median_tpm_ensgid.gct.rds \
    --perm_n 1000 \
    --out enrich
```

> 



# *Obsolete codes

## Prepare Roadmap enhancers

SNP residing on regulatory elements in Roadmap data: Run below command function at `bash`

Run a bash file `db/roadmap_dist.sh`:

```bash
sortBed -i db/roadmap/E001_25_imputed12marks_dense.bed | closestBed -d -a data/seedSNP_1817_bm.bed -b stdin > db/roadmap_dist/E001_snp_dist.tsv
...
sortBed -i db/roadmap/E129_25_imputed12marks_dense.bed | closestBed -d -a data/seedSNP_1817_bm.bed -b stdin > db/roadmap_dist/E129_snp_dist.tsv
```

```bash
bash roadmap_dist.sh
```

* Result files are stored at `db/roadmap_dist` folder.



Filtering by dist 0 and merge results from cell types

```bash
Rscript src/enrich.r --roadmap --dir db/roadmap_dist --out db
```

> ** Run roadmap_merge function in enrich.r **
>
> Read roadmap dist files... 127.. [1] 230759      7
> Write TSV file: dbroadmap_merge.tsv
>
> Job done: 2021-01-06 22:12:10 for 3.1 sec
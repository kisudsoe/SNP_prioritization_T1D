# PriSNP: Prioritization of SNPs

This is a project for prioritization of SNPs associated with Type 1 diabetes.



## Preparation: Data download

To download **GWAS Catalog data**, you can [search certain disease](https://www.ebi.ac.uk/gwas/). In this study, we downloaded [SNP-sets for type 1 diabetes](https://www.ebi.ac.uk/gwas/efotraits/EFO_0001359). Then you can run R code file for filtering the GWAS Catalog data as below command line:

```cmd
Rscript T1D_gwas.r [GWAS_file_path] [p-value_criteria]
```

<<<<<<< HEAD
To download **LDlink data**, you can run `T1D_ldlink.py` as below command line:

- To run the code, you need list of SNP RS IDs of dbSNP database as txt file

```CMD
python T1D_ldlink.py [SNP_file_path.txt]
```

To Filter the LDlink data, you can run `T1D_ldlink.r` as below command line:

```CMD
Rscript T1D_ldlink.r [SNP_file_path.txt] [LDlink_data_folder_path]
```

To use bedtools later, you have to prepare seedSNP file as [bed format](https://genome.ucsc.edu/FAQ/FAQformat.html). You can run `snp_bed.r` for generate bed file. But you should check `NA` values and fill it manually.

```CMD
Rscript snp_bed.r [seedSNP_file_path] [annotation_file_path;optional]
```



## Preparation 2. RoadMap data download and filter

The [RoadMAP project](https://egg2.wustl.edu/roadmap/web_portal/imputed.html) provides epigenome annotations such as [12-mark/127-reference epigenome/25-state Imputation Based Chromatin State Model](https://egg2.wustl.edu/roadmap/data/byFileType/chromhmmSegmentations/ChmmModels/imputed12marks/jointModel/final/). We downloaded the 127 files (e.g., E001_25_imputed12marks_hg38lift_dense.bed.gz and etc) using R code (`T1D_roadmap.r`). And then we filtered the data by [annotation code](https://egg2.wustl.edu/roadmap/web_portal/imputed.html) (see db/[roadmap] ) including 13_EnhA1, 14_EnhA2, 15_EnhAF, 16_EnhW1, 17_EnhW2, 18_EnhAc.

To download RoadMap data, you need to install `AnnotationHub` and `rtracklayer` in `BiocManaer` as below CMD code:


To download **LDlink data**, you can run Python code file as below command line:

- To run the code, you need list of SNP RS IDs of dbSNP database as txt file
>>>>>>> parent of 0c366a1... write: db/roadmap - core functions

```cmd
python T1D_ldlink.py [SNP_file_path.txt]
```



To Filter the LDlink data, you can run R code file as below command line:

```cmd
Rscript T1D_ldlink.r []
```




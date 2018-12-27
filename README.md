# PriSNP: Prioritization of SNPs

This is a project for prioritization of SNPs associated with Type 1 diabetes.



## Preparation: Data download

To download **GWAS Catalog data**, you can [search certain disease](https://www.ebi.ac.uk/gwas/). In this study, we downloaded [SNP-sets for type 1 diabetes](https://www.ebi.ac.uk/gwas/efotraits/EFO_0001359). Then you can run R code file for filtering the GWAS Catalog data as below command line:

```cmd
Rscript T1D_gwas.r [GWAS_file_path] [p-value_criteria]
```



To download **LDlink data**, you can run Python code file as below command line:

- To run the code, you need list of SNP RS IDs of dbSNP database as txt file

```cmd
python T1D_ldlink.py [SNP_file_path.txt]
```



To Filter the LDlink data, you can run R code file as below command line:

```cmd
Rscript T1D_ldlink.r []
```




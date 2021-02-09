# Monte Carlo permutation for random SNP sets

## 1. Generation of random SNP sets

To calculate SNP backgrounds for this pipeline, we generated 10,000 random control SNP sets from [dbSNP database (version 151, GRCh37p13)](ftp://ftp.ncbi.nih.gov/snp/organisms/human_9606_b151_GRCh37p13/BED/). To generate the random controls, you can run `src/rsid_random.r` as below `CMD` command line:

- This code takes very long time to complete (>hours)

```CMD
Rscript src/rsid_random.r 1817 10000
```

- Pre-generated random SNP sets are archived as zip files (e.g. random/seeds.zip).

To use the generated random SNP sets, you should sort the contents by using `bedtools sort` function as below `SHELL` command line:

- `bedtools sort -i random/rsid1.bed > random/seeds/rsid1_sort.bed`
- `bedtools sort -i random/rsid2.bed > random/seeds/rsid2_sort.bed`
- ...

```SHELL
$ bash random_startend.sh
```



## 2. Achieving random distribution

To archive the random distribution, you can use below R scripts and BASH commands:

- `$ roadmap_random.sh`
- `$ encode_random.sh`
- `src/bedtools_closest_random.r`
- `regulome_random.r`
- `gtex_random.r`
- `lncrnasnp_random.r`

### $ roadmap_random.sh

The `SHELL` script includes iterations of this code: `bedtools closest -d -a random/seeds/rsid1_sort.bed -b db/roadmap_enh_merge.bed > random/roadmap/roadmap_rsid1.tsv`

```SHELL
$ bash roadmap_random.sh
```

### > src/bedtools_closest_random.r

```CMD
> Rscript src/bedtools_closest_random.r roadmap
```

- Result file is `random/roadmap_dist.tsv`

### $ encode_random.sh

```SHELL
$ bash encode_random.sh
```

### > src/bedtools_closest_random.r

```CMD
> Rscript src/bedtools_closest_random.r encode
```

### > regulome_random.r

```CMD
Rscript regulome_random.r db/RegulomeDB.dbSNP132.Category1.txt.rds db/RegulomeDB.dbSNP132.Category2.txt.rds
```

### > gtex_random.r

```CMD
Rscript gtex_random.r
```

### > lncrnasnp_random.r

```CMD
Rscript lncrnasnp_random.r
```



## 3. Monte Carlo simulation

We used a Monte Carlo method to test for an enrichment for T1D SNPs occupied in enhancers. T1D SNPs were identified in 10,000 sets of 1,817 SNPs that were randomly selected from the [Single Nucleotide Polymorphism database](https://www.ncbi.nlm.nih.gov/projects/SNP/snp_summary.cgi) (dbSNP build 151, 10/06/2017).
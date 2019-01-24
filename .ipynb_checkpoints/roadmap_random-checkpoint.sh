#!/bin/sh
for i in $(seq 1 10000); do
    path='random/seeds/rsid'$i'.bed'
    name='random/roadmap/roadmap_rsid'$i'.tsv'
    #echo '    '$path' > '$name
    `bedtools closest -d -a $path -b db/roadmap_enh_merge.bed > $name`
done
#bedtools closest -d -a random/seeds/rsid1_sort.bed -b db/roadmap_enh_merge.bed > random/roadmap/roadmap_rsid1.tsv
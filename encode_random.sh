#!/bin/sh
for i in $(seq 1 10000); do
    path='random/seeds/rsid'$i'.bed'
    name='random/encode/encode_rsid'$i'.tsv'
    #echo '    '$path' > '$name
    `bedtools closest -d -a $path -b db/encode_tfbs_merge.bed > $name`
done
#bedtools closest -d -a data/seedSNP_1817.bed -b db/encode_tfbs_merge.bed > data/encode_dist.tsv
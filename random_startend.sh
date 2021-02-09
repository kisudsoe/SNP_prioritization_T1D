#!/bin/sh
#for path in `find random/* -maxdepth 1 -not -type d`
for i in $(seq 1 10000); do
    path='random/rsid'$i'.bed'
    #name=${path%%.*}'_.bed'
    name='random/seeds/rsid'$i'.bed'
    #echo '    '$path' > '$name
    `awk '$2>$3{print $1"\t"$3"\t"$2"\t"$4; next}{print $0}' $path | 
        bedtools sort -i stdin > $name`
done
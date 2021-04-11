#!/bin/bash

#set -x

for filename in communities/*.bed; do
    sort -k1,1 -k2,2n $filename > "communities/$(basename "$filename" .bed)_sort.bed"
    ./bedToBigBed "communities/$(basename "$filename" .bed)_sort.bed" hg19.chrom.sizes "communities/$(basename "$filename" .bed).bb"
    rm $filename
    rm "communities/$(basename "$filename" .bed)_sort.bed"
done


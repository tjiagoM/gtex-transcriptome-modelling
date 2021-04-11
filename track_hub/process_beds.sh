#!/bin/bash

#set -x

for filename in hg38/*.bed; do
    sort -k1,1 -k2,2n $filename > "hg38/$(basename "$filename" .bed)_sort.bed"
    ./bedToBigBed "hg38/$(basename "$filename" .bed)_sort.bed" hg38.chrom.sizes "hg38/$(basename "$filename" .bed).bb"
    rm $filename
    rm "hg38/$(basename "$filename" .bed)_sort.bed"
done


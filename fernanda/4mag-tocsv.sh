#!/bin/bash

file=$1
prefix=$2

cp $1 ${prefix}checkm.csv #make a copy of the checkm report to reformat it to readable csv file

sed -i 's/[ \t]/ /g' ${prefix}checkm.csv

sed -i 's/ {/, {/g' ${prefix}checkm.csv

sed -i 's/{//g' ${prefix}checkm.csv

sed -i 's/}//g' ${prefix}checkm.csv

cut -d, -f1,12,13 ${prefix}checkm.csv > small${prefix}check.csv #make a small version with only bin name, completeness and contamination

#rm ${prefix}checkm.csv #if you want to remove the big csv file remove the hashtag

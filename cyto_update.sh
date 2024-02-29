#!/bin/bash

#----
# Download
#----

wget "http://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/cytoBandIdeo.txt.gz" -O cytogenetic_bands_hg19.tsv.gz
wget "http://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/cytoBandIdeo.txt.gz" -O cytogenetic_bands_hg38.tsv.gz

truncate -s 0 cyto.csv
chrs="1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X Y M MT"

function extract_bands(){
raw_cyto_file=$1
chrs=$2
build=$3
zcat $1 | sed 's/chr//g' | cut -f1-4 | awk -v OFS="," -v chrs="${chrs}" -v build=${build} \
'BEGIN {split(chrs, chrArr, " "); for (i = 1; i in chrArr; ++i) {chrMap[chrArr[i]] = i}} # Assign each chromosome a number based on its position
$1 in chrMap {print chrMap[$1],$2,$3,$1$4,build}' | tr "\t" "," | sort -k1g,1 -k2g,2 -k3g,3 
}

for raw_cyto_file in cytogenetic_bands_hg19.tsv.gz cytogenetic_bands_hg38.tsv.gz
do
  build=`echo $raw_cyto_file | grep -q hg19 && echo hg19 || echo hg38`
  extract_bands $raw_cyto_file "${chrs}" $build >> cyto.csv
done

rm cytogenetic_bands_hg19.tsv.gz
rm cytogenetic_bands_hg38.tsv.gz


#----
# Add to database
#----

rm cytobands.db

sqlite3 -header -csv cytobands.db "CREATE TABLE cyto_pos (chr TINYINT, start INT, end INT, cyto TEXT, build TEXT)"
sqlite3 -header -csv cytobands.db ".import cyto.csv cyto_pos"

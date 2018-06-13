#!/bin/bash

#----
# Download
#----

echo "chr,start,end,cyto" > cyto.csv
wget "http://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/cytoBandIdeo.txt.gz" -O cytogenetic_bands.tsv.gz
zcat cytogenetic_bands.tsv.gz | sed 's/chr//g' | cut -f1-4 | awk -v OFS="," '$1 ~ /^[0-9]+$/ {print $1,$2,$3,$1$4}' | tr "\t" "," | sort -k1n,1 -k2n,2 -k3n,3 >> cyto.csv
rm cytogenetic_bands.tsv.gz

#----
# Add to database
#----
rm cytobands.db
sqlite3 -header -csv cytobands.db ".import cyto.csv cyto_pos"
#sqlite3 cytobands.db "PRAGMA foreign_keys=off; BEGIN TRANSACTION; ALTER TABLE cyto_pos MODIFY COLUMN chr TINYINT;COMMIT; PRAGMA foreign_keys=on;"

sqlite3 cytobands.db \
"PRAGMA foreign_keys=off;

BEGIN TRANSACTION;

ALTER TABLE cyto_pos RENAME TO _cyto_pos_old;

CREATE TABLE cyto_pos (
  chr TINYINT NULL,
  start INT NULL,
  end INT NULL,
  cyto TEXT
);

INSERT INTO cyto_pos (chr, start, end, cyto)
  SELECT chr, start, end, cyto
  FROM _cyto_pos_old;

COMMIT;
DROP TABLE _cyto_pos_old;
PRAGMA foreign_keys=on"